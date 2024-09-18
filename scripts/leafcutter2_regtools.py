#!/usr/bin/env python


__author__ = "Yang Li, Chao Dai, Quinn Hauck"
__status__ = "Development"
__version__ = "v2.0.0"

import argparse
import gzip
import os
import re
import shutil
import sys
import tempfile
from statistics import mean, median, stdev
from datetime import datetime

import BackwardSpliceJunctionClassifier as sjcb
import ForwardSpliceJunctionClassifier as sjcf
import pandas as pd
import pyfastx


def natural_sort(l: list):
    """Natural sort a list of string/tuple, similar to bash `sort -V`

    Parameters:
    -----------
    l : list
        l can be a list of string or numerics; or a list of varing length of tuples

    Returns:
    --------
    return : a sorted list
    """

    untuple = lambda tup: "".join([str(e) for e in tup])
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", untuple(key))]
    return sorted(l, key=alphanum_key)


def cluster_intervals(E: list):
    """Clusters intervals together

    Parameters:
    -----------
    E : list
        list of tuples, e.g. [(start, end)]

    Returns: tuple
    ---------
    Eclusters : list of list
        Each element is a list of tuples (introns) clustered into 1 cluster
    cluster: list
        A list of tuples of introns
    """
    E = natural_sort(E)
    current = E[0]
    Eclusters, cluster = [], []

    i = 0
    while i < len(E):

        if overlaps(E[i], current):
            cluster.append(E[i])
        else:
            Eclusters.append(cluster)
            cluster = [E[i]]
        current = (E[i][0], max([current[1], E[i][1]]))
        i += 1

    if len(cluster) > 0:

        Eclusters.append(cluster)

    return Eclusters, E


def overlaps(A: tuple, B: tuple):
    """Checks if A and B overlaps

    Parameters:
    -----------
    A : tuple
        start and end coordinates of 1 intron
    B : tuple
        start and end coordinates of another intron

    Returns:
    --------
    return : boolean
        Indicates whether genomic ranges of A and B overlap or not.
    """

    if A[1] < B[0] or B[1] < A[0]:
        return False
    else:
        return True


def pool_junc_reads(flist, options):
    """Pool junction reads

    Parameters
    ----------
    flist : str
        The file list
    options : argparse object
        Passed arguments

    Returns:
    --------
    return
        No returns. Use side effects.


    Side-effects:
    -------------
    write introns and counts by clusters. Output file is NOT versions sorted.
    format: [chrom]:[strand] [start]:[end]:[reads]
            e.g. chr17:+ 410646:413144:3 410646:413147:62
    """

    global chromLst

    outPrefix = options.outprefix
    rundir = options.rundir
    maxIntronLen = int(options.maxintronlen)
    checkchrom = options.checkchrom
    print(f"Max Intron Length: {maxIntronLen}")
    outFile = f"{rundir}/{outPrefix}_pooled"

    if not os.path.exists(rundir):
        os.mkdir(rundir)

    # store introns in `by_chrom`, a nested dictionary
    by_chrom = {}  # { k=(chrom, strand) : v={ k=(start, end) : v=reads } }

    for libl in flist:
        lib = libl.strip()
        if not os.path.isfile(lib):
            continue

        if options.verbose:
            sys.stderr.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} scanning {lib}...\n")

        if ".gz" in lib:
            F = gzip.open(lib)
        else:
            F = open(lib)

        for ln in F:

            if type(ln) == bytes:
                ln = ln.decode("utf-8")  # convert bytes to string

            lnsplit = ln.split()
            if len(lnsplit) < 6:
                sys.stderr.write(f"Error in {lib} \n")
                continue

            if len(lnsplit) == 12:  # 12 fields regtools junc file
                (
                    chrom,
                    A,
                    B,
                    dot,
                    counts,
                    strand,
                    rA,
                    rb,
                    rgb,
                    blockCount,
                    blockSize,
                    blockStarts,
                ) = lnsplit
                if int(blockCount) > 2:
                    print(ln, "ignored...")
                    continue
                Aoff, Boff = blockSize.split(",")[:2]
                A, B = int(A) + int(Aoff), int(B) - int(Boff)  # get intron

            elif len(lnsplit) == 6:
                # old leafcutter junctions
                chrom, A, B, dot, counts, strand = lnsplit
                A, B = int(A), int(B)

            if checkchrom and (chrom not in chromLst):
                continue

            A, B = int(A), int(B) + int(options.offset)

            if B - A > int(maxIntronLen):
                continue

            # sum up all the reads at the same junctions if junctions already exist
            try:
                by_chrom[(chrom, strand)][(A, B)] = (
                    int(counts) + by_chrom[(chrom, strand)][(A, B)]
                )
            except:
                try:
                    by_chrom[(chrom, strand)][(A, B)] = int(
                        counts
                    )  # when only 1 junction
                except:
                    by_chrom[(chrom, strand)] = {
                        (A, B): int(counts)
                    }  # when only 1 junction

    with open(outFile, "w") as fout:
        Ncluster = 0
        sys.stderr.write("Parsing...\n")

        for chrom in by_chrom:
            read_ks = [
                k for k, v in by_chrom[chrom].items() if v >= 3
            ]  # read-keys, require junction reads > 3
            read_ks.sort()  # sort read-keys: (start, end)

            sys.stderr.write(f"{chrom[0]}:{chrom[1]}..\n")

            if read_ks:
                clu = cluster_intervals(read_ks)[
                    0
                ]  # clusters of introns, [[(start, end),..],..]

                for cl in clu:
                    if len(cl) > 0:  # 1 if cluster has more than one intron
                        buf = f"{chrom[0]}:{chrom[1]} "  # chr:strand
                        for interval, count in [(x, by_chrom[chrom][x]) for x in cl]:
                            buf += (
                                f"{interval[0]}:{interval[1]}" + f":{count}" + " "
                            )  # start:end:reads
                        fout.write(buf + "\n")
                        Ncluster += 1

        sys.stderr.write(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Wrote {Ncluster} clusters..\n")


def refine_linked(clusters):
    """Re-cluster introns into clusters of linked introns

    Linked introns are introns that share either 5' or 3' splice site

    Parameters:
    -----------
    clusters : tuple
        format is [((start, end), reads)], eg.
        [((413430, 423479), 3), ((410646, 413144), 3), ((410646, 413147), 62), ((410646, 420671), 4)]

    Returns:
    --------
    return : list of list
        base element is a tuple, format: ((start, end), reads). e.g.
        [[((413430, 423479), 3)], [((410646, 413144), 3), ((410646, 413147), 62), ((410646, 420671), 4)]]
    """

    unassigned = [x for x in clusters[1:]]
    current = [clusters[0]]
    splicesites = set([current[0][0][0], current[0][0][1]])  # start & end of intron
    newClusters = []
    while len(unassigned) > 0:
        finished = False

        while not finished:
            finished = True
            torm = []
            for intron in unassigned:
                (start, end), count = intron
                if start in splicesites or end in splicesites:
                    current.append(intron)
                    splicesites.add(start)
                    splicesites.add(end)
                    finished = False
                    torm.append(intron)
            for intron in torm:
                unassigned.remove(intron)
        newClusters.append(current)
        current = []
        if len(unassigned) > 0:
            current = [unassigned[0]]
            splicesites = set([current[0][0][0], current[0][0][1]])
            if len(unassigned) == 1:  # assign the last intron
                newClusters.append(current)
            unassigned = unassigned[1:]
    return newClusters


def refine_cluster(clu: list, cutoff: float, readcutoff: int):
    """Filter introns based on cutoffs

    Parameters:
    -----------
    clu : list of tuples, a single cluster
        list of tuples, each tuple an intron of the cluster
    cutoff : float
        reads ratio cutoff, passed in from option --mincluratio
    readcutoff : int
        minimum reads cutoff, passed in from option --minreads

    Filters:
    --------
        1. compute ratio of reads for each intron in a cluster
        2. remove intron if:
            - ratio < ratio_cutoff
            - OR reads of the intron < readcutoff
        3. re-cluster with remaining introns

    Returns:
    --------
        return : list of list
        list of refined (filtered) clusters

    """

    remove = []
    dic = {}
    intervals = []

    reCLU = False  # re-cluster flag
    totN = 0

    for inter, count in clu:
        totN += count
    for inter, count in clu:
        if count / float(totN) >= cutoff and count >= readcutoff:
            intervals.append(inter)
            dic[inter] = count  # {(start, end): reads}
        else:
            reCLU = True  # any intron not passing filters will enforce reCLU

    if len(intervals) == 0:
        return []  # base case

    # Below makes sure that after trimming/filtering, the clusters are still good
    # afterwards - each clusters have linked introns that pass filters.

    Atmp, _ = cluster_intervals(intervals)
    A = []

    # A: a list of linked intron clusters
    if len(Atmp) > 0:
        for cl in Atmp:  # Atmp is a list of list
            if len(cl) == 1:  # 1 intron
                A.append(cl)
            for c in refine_linked([(x, 0) for x in cl]):  # >1 introns
                if len(c) > 0:
                    A.append([x[0] for x in c])

    if len(A) == 1:  # A has 1 cluster of introns
        rc = [(x, dic[x]) for x in A[0]]

        if len(rc) > 0:
            if reCLU:  # recompute because ratio changed after removal of some introns
                return refine_cluster(
                    [(x, dic[x]) for x in A[0]], cutoff, readcutoff
                )  # recursive
            else:
                return [[(x, dic[x]) for x in A[0]]]

    NCs = []  # As in N Clusters, here A has more than 1 clusters of introns
    for c in A:
        if len(c) > 1:  # c has more than 1 introns
            NC = refine_cluster(
                [(x, dic[x]) for x in c], cutoff, readcutoff
            )  # recursive
            NCs += NC

    return NCs


def refine_clusters(options):
    """Refine clusters.

    Refine clusters such that kept clusters that are written to file meets
    the following criteria:
        * introns a linked (share either 5' or 3' splice site)
        * minimum total cluster reads cutoff
        * minimum intron reads cutoff
        * minimum reads ratio per cluster cutoff

    However, if constitutive flag `const` is on, then non-linked introns
    are also written out, and they do not subject to cluster ratio and
    cluster reads cutoff filters.

    Parameters:
    -----------
        options : argparse object

    Returns:
    --------
        return : no returns. Use side-effects

    Side-effects:
    -------------
        write refined clusters to file - `*_refined`. Output file is NOT
        version sorted.

    """

    outPrefix = options.outprefix
    rundir = options.rundir
    minratio = float(options.mincluratio)
    minclureads = int(options.minclureads)
    minreads = int(options.minreads)

    inFile = f"{rundir}/{outPrefix}_pooled"
    outFile = f"{rundir}/{outPrefix}_refined"
    fout = open(outFile, "w")

    sys.stderr.write(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Refine clusters from {inFile}...\n")

    Ncl = 0
    for ln in open(inFile):  # pooled juncs
        if type(ln) == bytes:
            ln = ln.decode("utf-8")  # convert bytes to string

        clu = []  # each cluster: [((start, end), reads),..]
        totN = 0  # total cluster reads
        chrom = ln.split()[0]
        for ex in ln.split()[1:]:  # for an exon
            A, B, N = ex.split(":")
            clu.append(((int(A), int(B)), int(N)))
            totN += int(N)

        if totN < minclureads:
            continue

        if (
            options.const
        ):  # include constitutive introns. These are clusters that only have 1 intron, hence "constitutive"
            if len(clu) == 1:
                buf = f"{chrom} "
                for interval, count in clu:
                    buf += f"{interval[0]}:{interval[1]}" + f":{count}" + " "
                Ncl += 1
                fout.write(buf + "\n")  # e.g. 'chr:strand start:end:reads'

        for cl in refine_linked(clu):  # only linked intron clusters
            rc = refine_cluster(cl, minratio, minreads)
            if len(rc) > 0:
                for clu in rc:
                    buf = f"{chrom} "
                    for interval, count in clu:
                        buf += f"{interval[0]}:{interval[1]}" + f":{count}" + " "
                    Ncl += 1
                    fout.write(buf + "\n")
    sys.stderr.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Split into {Ncl} clusters...\n")
    fout.close()


def addlowusage(options):
    """Add low usage introns to refined clusters

    Parameters:
    -----------
    options : argparse object
        pass in command options


    Returns:
    --------
    return : null
        no returns. Write files in side-effects.

    Side-effects:
    ------------
        written files:
            - [out_prefix]_lowusage_introns : file stores low usage introns (by cluster).
              Output file is version sorted.
            - [out_prefix]_refined_noisy    : file stores all usage introns (by cluster),
              although each cluster must pass min cluster reads cutoff. Output file is
              version sorted.

    """

    global chromLst

    sys.stderr.write(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Add low usage introns...\n")

    outPrefix = options.outprefix
    rundir = options.rundir
    pooled = f"{rundir}/{outPrefix}_pooled"
    minclureads = int(options.minclureads)
    minreads = int(options.minreads)

    if options.cluster == None:
        refined_cluster = f"{rundir}/{outPrefix}_refined"
    else:
        refined_cluster = options.cluster

    outFile = (
        f"{rundir}/{outPrefix}_refined_noisy"  # out file that includes noisy introns
    )
    outFile_lowusageintrons = (
        f"{rundir}/{outPrefix}_lowusage_introns"  # out file for lowusage introns
    )

    fout = open(outFile, "w")
    fout_lowusage = open(outFile_lowusageintrons, "w")

    # get 5' sites, 5' sites, and clusters of introns from refined file, see data structure below
    exons5, exons3, cluExons = {}, {}, {}
    cluN = 0  # clusterID

    # construct 5' sites, 3' sites, and clusters dict from refined
    for ln in open(refined_cluster):
        if type(ln) == bytes:
            ln = ln.decode("utf-8")  # convert bytes to string

        chrom = ln.split()[0]
        cluN += 1
        cluExons[(chrom, cluN)] = []  # keys are (chrom, cluster_number)
        for exon in ln.split()[1:]:
            A, B, count = exon.split(":")  # start, end, reads
            if chrom not in exons5:
                exons5[chrom] = {}
                exons3[chrom] = {}
            exons5[chrom][int(A)] = (
                chrom,
                cluN,
            )  # 5' sites, { k=chrom, v={ k=start, v=(chrom, clusterID) } }
            exons3[chrom][int(B)] = (
                chrom,
                cluN,
            )  # 3' sites, { k=chrom, v={ k=end, v=(chrom, clusterID) } }
            cluExons[(chrom, cluN)].append(
                exon
            )  # introns, { k=(chrom, clusterID), v=['start:end:reads'] }

    # Below for loop adds back clusters (previously filtered out in refined_clusters)
    # in the pooled junc file. These previously removed introns are added to all
    # cluExons, as well as to lowusage_intron
    lowusage_intron = {}  # { k=(chrom, clusterID), v=['start:end:reads'...]}
    for ln in open(pooled):  # each cluster/line from pool_juncs file

        if type(ln) == bytes:
            ln = ln.decode("utf-8")  # convert bytes to string

        clu = []
        totN = 0
        chrom = ln.split()[0]

        if chrom in exons5:  # ensure chrom is in exons5 level-1 keys

            # Below for loop adds introns that were filtered out in refined, aka noisy_introns,
            # back to a total intron cluster dict and to a lowusage (noisy) intron cluster dict
            for exon in ln.split()[1:]:
                A, B, N = exon.split(":")  # start, end, reads

                # when 5' site in refined
                if int(A) in exons5[chrom]:
                    clu = exons5[chrom][
                        int(A)
                    ]  # set clu=(chrom, clusterID), key for cluExons
                    if exon not in cluExons[clu]:  # exon was filtered out by refined
                        cluExons[clu].append(exon)  # add it to cluExons
                        if clu not in lowusage_intron:
                            lowusage_intron[clu] = []
                        lowusage_intron[clu].append(exon)  # also add it to lowusage

                # else when 3' site in refined, perform same procedure
                elif int(B) in exons3[chrom]:  # when 3' site is in refined
                    clu = exons3[chrom][int(B)]
                    if exon not in cluExons[clu]:
                        cluExons[clu].append(exon)
                        if clu not in lowusage_intron:
                            lowusage_intron[clu] = []
                        lowusage_intron[clu].append(exon)

                # neither 5' nor 3' splice site in refined, only add cluster if intron meets minreads requirement
                # because of the minreads requirement, this intron is not noisy, thus do not add to lowusage_intron
                else:
                    if int(N) > minreads:
                        cluN += 1
                        cluExons[(chrom, cluN)] = [exon]

    # write low usage introns
    ks = natural_sort(
        lowusage_intron.keys()
    )  # e.g. { k=(chrom, clusterID), v=['start:end:reads'...]}
    for clu in ks:  # e.g. (chrom, clusterID)
        fout_lowusage.write(clu[0] + " " + " ".join(lowusage_intron[clu]) + "\n")
    fout_lowusage.close()

    # write all intron clusters
    cluLst = natural_sort(cluExons.keys())
    for clu in cluLst:
        if not options.const:  # if -C flag not set, do not write constitutive introns
            if len(cluExons[clu]) == 1:
                continue  # skip write out if only 1 intron in cluster, aka, constitutive

        # only write introns if minimum cluster reads criteria is met
        if sum([int(ex.split(":")[-1]) for ex in cluExons[clu]]) < minclureads:
            continue
        chrom = clu[0]
        buf = f"{chrom}"
        for ex in cluExons[clu]:
            buf += " " + ex
        fout.write(buf + "\n")
    fout.close()


def sort_junctions(libl, options):
    """Sort junctions by cluster

    For each intron cluster, sort introns. Write both numerator (intron reads) and
    denominator (cluster reads) into output file.

    Parameters:
    -----------
        libl : str
            A list of junction file paths
        options: argparse object
            Attributes store command line options

    Returns:
    --------
        return : no returns. Use site effect.

    Side-effects:
    -------------
        text file : '{rundir}/{outPrefix}_sortedLibs'
            store junfile names that are processed/sorted
        text file:  '{rundir}/{outPrefix}...sorted.gz'
            a series of sorted input junction files, sorted.
    """

    global chromLst

    outPrefix = options.outprefix
    rundir = options.rundir
    checkchrom = options.checkchrom

    if options.cluster == None:  # if not providing refined clusters externally
        refined_cluster = (
            f"{rundir}/{outPrefix}_refined_noisy"  # note refined noisy intron clusters
        )
        sys.stderr.write(f"\nUsing {refined_cluster} as refined cluster...\n")
    else:
        refined_cluster = options.cluster

    runName = f"{rundir}/{outPrefix}"

    # exons:  { k=chrom : v={ k=(start, end) : v=clusterID } }
    # cluExons:  { k=clusterID : v=[(chrom, start, end)] }
    exons, cluExons = {}, {}
    cluN = 0  # clusterID

    # fill in exons, cluExons dict from `*refined_noisy` intron cluster file
    for ln in open(
        refined_cluster
    ):  # e.g. ln = "chr10:+ 135203:179993:5 135302:179993:29"

        if type(ln) == bytes:
            ln = ln.decode("utf-8")  # convert bytes to string

        chrom = ln.split()[0]  # e.g. "chr10:+"
        cluN += 1
        for exon in ln.split()[1:]:  # e.g. "135203:179993:5 135302:179993:29"
            A, B, count = exon.split(":")

            if chrom not in exons:
                exons[chrom] = {}
            if (int(A), int(B)) not in exons[chrom]:
                exons[chrom][(int(A), int(B))] = [cluN]
            else:
                exons[chrom][(int(A), int(B))].append(cluN)
            if cluN not in cluExons:
                cluExons[cluN] = []
            cluExons[cluN].append((chrom, A, B))

    merges = {}  # stores junc file names as dict { k=filename : v=[filename] }
    for ll in libl:
        lib = ll.rstrip()  # 1 junc file path
        libN = lib.split("/")[
            -1
        ]  # get library name from junc file name eg. GTEX-1117F-0226-SM-5GZZ7.tsv.gz
        if not os.path.isfile(lib):
            continue
        if libN not in merges:
            merges[libN] = []  # why use list, `libN` should always be one element
        merges[libN].append(lib)

    fout_runlibs = open(
        os.path.join(rundir, outPrefix) + "_sortedlibs", "w"
    )  # intermediate file to store sorted junc file names to be written

    # operate on each libN (library), each libN can have more than 1+ junc files
    for libN in merges:
        by_chrom = {}  # to store junctions from original unsorted junc file

        # write sorted junc file names into intermediate file
        foutName = os.path.join(
            rundir, outPrefix + "_" + libN + ".junc.sorted.gz"
        )  # 'test/gtex_w_clu/gtex_GTEX-1IDJU-0006-SM-CMKFK.junc.sorted.gz'
        fout_runlibs.write(foutName + "\n")  # e.g. 'test/gtex_w_clu/gtex_sortedlibs'

        if options.verbose:
            sys.stderr.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Sorting {libN}..\n")
        if len(merges[libN]) > 1:
            if options.verbose:
                sys.stderr.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} merging {' '.join(merges[libN])}...\n")
        else:
            pass
        fout = gzip.open(
            foutName, "wt"
        )  # e.g. 'test/gtex_w_clu/gtex_GTEX-1IDJU-0006-SM-CMKFK.junc.sorted.gz'

        # -------- Process and write junction files --------

        # write header
        fout.write(f"chrom {libN}\n")  # 'chrom GTEX-111VG-0526-SM-5N9BW\n'

        # -------- Gather counts from all junc files of library --------
        # store in by_chrom: { ('chr1', '+') : { (100, 300) : 5, (500, 700): 10, ... } }
        for lib in merges[libN]:
            if ".gz" in lib:
                F = gzip.open(lib)
            else:
                F = open(lib)

            for ln in F:  # 1 line: e.g. "chr17\t81701131\t81701534\t.\t1\t+"

                if type(ln) == bytes:
                    ln = ln.decode("utf-8")  # convert bytes to string
                lnsplit = ln.split()

                if len(lnsplit) < 6:
                    sys.stderr.write(f"Error in {lib} \n")
                    continue

                if len(lnsplit) == 12:
                    (
                        chrom,
                        A,
                        B,
                        dot,
                        counts,
                        strand,
                        rA,
                        rb,
                        rgb,
                        blockCount,
                        blockSize,
                        blockStarts,
                    ) = lnsplit
                    if int(blockCount) > 2:
                        print(ln, "ignored...")
                        continue
                    Aoff, Boff = blockSize.split(",")[:2]
                    A, B = int(A) + int(Aoff), int(B) - int(Boff)

                elif len(lnsplit) == 6:
                    # old leafcutter junctions
                    chrom, A, B, dot, counts, strand = lnsplit
                    A, B = int(A), int(B)

                A, B = int(A), int(B) + int(options.offset)  # start, end + offset

                chrom = (chrom, strand)
                if chrom not in by_chrom:
                    by_chrom[chrom] = (
                        {}
                    )  # store introns from junc file, key: ('chr1', '+')

                intron = (A, B)
                if intron in by_chrom[chrom]:  # sum up reads by intron from junc files
                    by_chrom[chrom][intron] += int(counts)
                else:
                    by_chrom[chrom][intron] = int(counts)

        # ------- Take clusters from refined_noisy, assign reads -------
        # reads are from by_chrom (junc files)
        # For each intron cluster, write fraction for each intron (one intron per line).
        for clu in cluExons:  # cluExons: { k=cluID : v=[(chrom, start, end)...]}
            buf = []
            ks = cluExons[clu]  # eg: [('chr1:+', 827776, 829002), ..] introns of a clu
            ks.sort()  # no need to version sort within cluster

            # Step 1: sum cluster level reads from each intron
            # gather (sum) reads for each cluster in refined_noisy, read counts are from junc file (by_chrom)
            tot = 0  # sum of total read counts per cluster
            usages = []
            for exon in ks:
                chrom, start, end = exon
                chrom = tuple(
                    chrom.split(":")
                )  # convert 'chr3:+' to ('chr3', '+') as in by_chrom
                start, end = int(start), int(end)

                if chrom not in by_chrom:
                    pass
                elif (start, end) in by_chrom[chrom]:
                    tot += by_chrom[chrom][(start, end)]

            # Step 2: append intron usage fraction to stream buffer
            for exon in ks:
                chrom, start, end = exon
                start, end = int(start), int(end)
                chrom = tuple(chrom.split(":"))
                chromID, strand = chrom  # chromID eg: 'chr3'

                intron = (
                    chromID,
                    start,
                    end + 1,
                    strand,
                )  # converting to 1-based coordinates

                if chrom not in by_chrom:
                    # if refined exon chrom is not found in junc file, write 0/cluster_total
                    buf.append(f"{chromID}:{start}:{end}:clu_{clu}_{strand} 0/{tot}\n")
                elif (start, end) in by_chrom[chrom]:
                    # if refind exon is in junc file, write exon reads / cluster_total
                    buf.append(
                        f"{chromID}:{start}:{end}:clu_{clu}_{strand} {by_chrom[chrom][(start,end)]}/{tot}\n"
                    )
                else:
                    # if refined exon is not found in junc file, write 0/cluster_total
                    buf.append(f"{chromID}:{start}:{end}:clu_{clu}_{strand} 0/{tot}\n")

            fout.write("".join(buf))
        fout.close()
    fout_runlibs.close()


def merge_files(fnames, fout, options):
    """Merge a list of files into a gzip file

    Parameters:
    -----------
    fnames : list
        list of file path. This is a list of files within a single batch
    fout : a single opened file io
        specifically, a gzip.open('w') io object
    options :
        argparse object

    Returns:
    --------
        return : no returns. Use side-effects.

    Side-effects:
    -------------
        After merging all files from `fnames` list, write out into a gzipped
    file io, as opened in `fout`.
    """

    fopen = []
    for fname in fnames:  # each file to be merged
        if fname[-3:] == ".gz":
            fopen.append(gzip.open(fname))
        else:
            fopen.append(open(fname))

    finished = False
    N = 0
    while not finished:  # cycle through files in batch
        N += 1
        if N % 50000 == 0:
            sys.stderr.write(".")
        buf = []
        for f in fopen:  # each opened file
            ln = f.readline().decode().split()  # read 1 line
            if len(ln) == 0:  # end of line = finish
                finished = True
                break
            chrom = ln[0]  # e.g. "chrom" or "chr1:825552:829002:clu_1_+"
            data = ln[1:]  # e.g. "GTEX-1117F-0626-SM-5N9CS.leafcutter" or "0/0"
            if len(buf) == 0:
                buf.append(chrom)
            buf += data  # e.g. ['chrom', 'GTEX-111VG-0526-SM-5N9BW.leafcutter', 'GTEX-1117F-0626-SM-5N9CS.leafcutter'] for first lines, or ['chr1:825552:829002:clu_1_+', '0/0', '0/0'] for 2+ lines
            # each file the exact same chromosome coordinates, effectively we are collecting counts into columns 2 and after

        if len(buf) > 0:
            if buf[0] == "chrom":
                if options.verbose:
                    sys.stderr.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} merging {len(buf)-1} files")
            fout.write(" ".join(buf) + "\n")  # combining sample counts into columns
        else:
            break

    sys.stderr.write(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}  done.\n")
    for fin in fopen:
        fin.close()


def merge_junctions(options):
    """Merge junctions

    Merge a list of sorted junction files into a single merged junction file.
    Each input sorted junction files must have the same introns, i.e. first
    column of each row must be the same across all files to be merged.

    Parameters:
    -----------
    options : argparse object

    Returns:
    ---------
    return : null
        No returns. Use side effect.

    Side-effects:
    -------------
        Collect previously sorted junction files. Merge junction files in batches.
    And finally, all batches are merged into a single file. Reads fractions are in
    columns.
        row1  : col1=`chrom`, col2 and beyond are input file names merged
        row2+ : col1=`intron identifier`, reads fraction from each input file
    """

    outPrefix = options.outprefix
    rundir = options.rundir

    fnameout = os.path.join(f"{rundir}/{outPrefix}")
    flist = fnameout + "_sortedlibs"  # sorted juncs file list
    lsts = []  # = flist
    for ln in open(flist):
        if type(ln) == bytes:
            ln = ln.decode("utf-8")  # convert bytes to string
        lsts.append(ln.strip())

    if options.verbose:
        sys.stderr.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} \nMerging {len(lsts)} junction files...\n")

    # Change 300 if max open file is < 300
    # set up batch N per batch
    N = min([300, max([100, int(len(lsts) ** (0.5))])])

    # tmpfiles = []
    while len(lsts) > 1:  # initial list of sorted junc files

        # convert lsts (list of file paths) to clst (list of lists)
        # each sublist is a batch of upto 100 files.
        clst = []  # list of batches, each batch has up to 100 sorted junc files
        for i in range(
            0, int(len(lsts) / N) + 1
        ):  # merge in batches of max(100, len(lsts))
            lst = lsts[N * i : N * (i + 1)]
            if len(lst) > 0:
                clst.append(lst)
        lsts = (
            []
        )  # clear initial file list, now repurposed to store merged file names (temp)

        for lst in clst:  # run by batch
            if len(lst) == 0:
                continue
            tmpfile = tempfile.mktemp()
            os.mkdir(tmpfile)
            foutname = tmpfile + "/tmpmerge.gz"
            fout = gzip.open(
                foutname, "wt"
            )  # create a temp file for the batch of files to be merged

            merge_files(lst, fout, options)  # merge the batch into `fout`
            lsts.append(foutname)  # save the temp merged file name
            # tmpfiles.append(foutname) # this line is not needed.
            fout.close()

    if not options.const:
        shutil.move(lsts[0], fnameout + "_perind.counts.gz")
    else:
        shutil.move(lsts[0], fnameout + "_perind.constcounts.gz")


def get_numers(options):
    """Get numerators from merged count table

    Parameters:
    -----------
    options : argparse object

    Returns:
    --------
    return : null
        No returns. Use side-effect to write out file.

    Side-effects:
    -------------
        Take in count tables, extract numerators for each sample.

    """

    outPrefix = options.outprefix
    rundir = options.rundir

    if not options.const:
        fname = f"{rundir}/{outPrefix}_perind.counts.gz"
        fnameout = f"{rundir}/{outPrefix}_perind_numers.counts.gz"
    else:
        fname = f"{rundir}/{outPrefix}_perind.constcounts.gz"
        fnameout = f"{rundir}/{outPrefix}_perind_numers.constcounts.gz"

    input_file = gzip.open(fname, "r")
    fout = gzip.open(fnameout, "wt")
    first_line = True

    sys.stderr.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} \nExtracting numerators (read counts) from {fname}...")

    for l in input_file:
        if first_line:
            fout.write(
                " ".join(l.decode().strip().split(" ")[1:]) + "\n"
            )  # print the sample names
            first_line = False
        else:
            l = l.decode().strip()
            words = l.split(" ")
            fout.write(
                words[0] + " " + " ".join([g.split("/")[0] for g in words[1:]]) + "\n"
            )  # write intron and numerators

    input_file.close()
    fout.close()
    sys.stderr.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} done.\n")


# def loadIntronAnnotations(fn):
#     """parse and load intron annotations from a single file"""
#
#     dic_noise = (
#         {}
#     )  # noisy intron annotation: { k=(chr, start, end, strand) : v=classfication }
#     if fn[-3:] == ".gz":
#         for ln in gzip.open(fn):
#             chrom, s, e, strand, classification = ln.decode().split()
#             # dic_noise[(chrom,int(s)-1,int(e),strand)] = classification # intron annotation
#             dic_noise[(chrom, int(s), int(e), strand)] = (
#                 classification  # intron annotation
#             )
#     else:
#         for ln in open(fn):
#             if type(ln) == bytes:
#                 ln = ln.decode("utf-8")  # convert bytes to string
#             chrom, s, e, strand, classification = ln.split()
#             # dic_noise[(chrom,int(s)-1,int(e),strand)] = classification
#             dic_noise[(chrom, int(s), int(e), strand)] = classification
#     return dic_noise


def pick_sjc_label(df):
    '''Pick one label for each intron based on priroty'''
    sjc_priority = {'PR': 1, 'UP': 2, 'NE':3} # PR > UP > NE
    df['priority'] = df['SJClass'].map(sjc_priority)
    chosen = df.sort_values(by = 'priority', ascending = True).iloc[0]
    return chosen.drop('priority')


def merge_discordant_logics(sjc_file: str):
    '''some junctions have multiple classifications. Use conservative approach
    to merge them.
    '''

    classifier = {
        # each bit represents:
        # is_GTFAnnotatedCoding?  is_GTFAnnotated?  is_LF2AnnotatedCoding?  is_ClosetoUTR?
        '0000': 'UP', # UnProductive,
        '0001': 'NE', # NEither productive nor unproductive
        '0010': 'PR', # PRoductive
        '0011': 'PR', # PRoductive
        '0100': 'UP', # UnProductive
        '0101': 'NE', # Neither Productive nor UnProductive
        '0110': 'PR', # PRoductive
        '0111': 'PR', # PRodutive
        '1000': 'PR', # PRoductive
        '1001': 'PR', # PRoductive
        '1010': 'PR', # PRoductive
        '1011': 'PR', # PRoductive
        '1100': 'PR', # PRoductive
        '1101': 'PR', # PRoductive
        '1110': 'PR', # PRoductive
        '1111': 'PR'  # PRoductive
        }

    sjc = pd.read_csv(sjc_file, sep = "\t")
    classifer_3bits = {
        # each bit represents:
        # is_GTFAnnotated?  is_LF2AnnotatedCoding?  is_ClosetoUTR?
        '000': 'UP', # UnProductive,
        '001': 'NE', # NEither productive nor unproductive
        '010': 'PR', # PRoductive
        '011': 'PR', # PRoductive
        '100': 'UP', # UnProductive
        '101': 'NE', # Neither Productive nor UnProductive
        '110': 'PR', # PRoductive
        '111': 'PR', # PRodutive
    }


    # group dt; NOTE:ForwardSpliceJunctionClassifier has an extra Strand column, backward doesn't
    if 'Strand' in sjc.columns:
        sjc = sjc[['Gene_name', 'Intron_coord', 'Strand', 'GencodePC', 'Annot', 'Coding', 'UTR']]

        # convert Annotation, Coding, UTR status to SJ categories
        sjc['SJClass'] = sjc[['GencodePC', 'Annot', 'Coding', 'UTR'
                              ]].astype(int).astype(str).agg(''.join, axis=1).map(classifier)

        # if multiple classifications, take the one with highest priority
        sjc = sjc.groupby(['Intron_coord', 'Strand']).apply(pick_sjc_label).reset_index(drop=True)
        sjc = sjc[['Intron_coord', 'Strand', 'SJClass', 'Gene_name']]

        # convert df to dict
        sjc = sjc.set_index(['Intron_coord', 'Strand']).to_dict(orient='index')
    else:
        sjc = sjc[['Gene_name', 'Intron_coord', 'Annot', 'Coding', 'UTR']]
        
        # convert Annotation, Coding, UTR status to SJ categories
        sjc['SJClass'] = sjc[['Annot', 'Coding', 'UTR'
                              ]].astype(int).astype(str).agg(''.join, axis=1).map(classifer_3bits)
        
        # if multiple classifications, take the one with highest priority
        sjc = sjc.groupby(['Intron_coord', 'Strand']).apply(pick_sjc_label).reset_index(drop=True)
        sjc = sjc[['Intron_coord', 'Strand', 'SJClass', 'Gene_name']]

        # convert df to dict
        sjc = sjc.set_index('Intron_coord').to_dict(orient='index')

    sjc = {flatten_tuple(k): v for k, v in sjc.items()}
    sys.stderr.write(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Loaded junction annotation lookup.\n")

    # sjc is a dcitionary with:
    # - keys: intron coordinates, e.g. ('chr1', 1000, 2000, '+') or ('chr1', 1000, 2000) for backward
    # - values: a dictionary e.g. {'SJClass': 'UP', 'Gene_name': 'DNMBP'}
    return sjc

 
def boolean_to_bit(bool_vec):
    # Convert boolean vector to string of "1"s and "0"s
    bin_str = "".join(["1" if b else "0" for b in bool_vec])

    return bin_str


def flatten_tuple(t):
    # t: tuple like ('chr1:100-200', '+')' or str 'chr1:100-200'
    if isinstance(t, tuple):
        c, ab = t[0].split(":")
    elif isinstance(t, str):
        c, ab = t.split(":")
    a, b = ab.split("-")
    a, b = int(a), int(b)
    if isinstance(t, tuple):
        s = t[1]
        return((c, a, b, s)) # e.g. ('chr1', 100, 200, '+')
    else:
        return((c, a, b))


def annotate_noisy(options):
    """Annotate introns

    Produces 3 files. Details in side-effects.

    Parameters:
    -----------
        options : argparse object

    Returns:
    --------
        return : no return. Use side-effects to write output files.

    Side-effects:
    -------------
        noisy counts :  output result file.
            Count table that has noisy introns marked
            with `*`. Each intron cluster only output up to 1 `pseudo` noisy
            intron coordinates. This `pseudo` coordinates are constructed by
            taking the min(start) and max(end) of all functional intron
            coordinates, to create the longest range. Reads from all noisy
            introns of a cluster are added together.

        noisy numerators : output result file.
            Same as noisy counts, except here only reports the numerators.
            Also, used `*` to mark both start and end pseudo coordinates
            of noisy intron.

        noisy counts by intron : output file for diagnostics.
            Same count table as the output of sorted juncs count table,
            except that each intron is annotated with `F`, `N`, or `PN` to
            denote `putative_functional`, `noisy`, or `putative_noisy`.

    """

    outPrefix = options.outprefix
    rundir = options.rundir
    minreadstd = float(options.minreadstd)
    fnameout = f"{rundir}/{outPrefix}"
    sjc_f = f"{rundir}/{outPrefix}_junction_classifications.txt"  # Classified junction annotations

    sys.stderr.write(
        f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} \nAnnotating introns with custom-classified annotations {sjc_f}...\n"
    )

    if sjc_f != None:

        if options.verbose:
            sys.stderr.write(
                f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Loading {sjc_f} for (un/)productive splicing classification..\n"
            )
        sjc = merge_discordant_logics(sjc_f)

        if options.verbose:
            sys.stderr.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Loaded..\n")

    if not options.const:
        fname = fnameout + "_perind.counts.gz"  # no constitutive introns
    else:
        fname = fnameout + "_perind.constcounts.gz"

    noisydiag = fname.replace(
        ".gz", ".noise_by_intron.gz"
    )  # eg: run/out_perind.counts.noise_by_intron.gz
    numersdiag = fname.replace(".gz", ".noise_by_intron.gz").replace(
        "perind", "perind_numers"
    )  # eg: run/out_perind_numers.counts.noise.gz

    foutdiag = gzip.open(noisydiag, "wt")
    foutdiagnumers = gzip.open(numersdiag, "wt")

    F = gzip.open(fname)
    ln = F.readline().decode()
    foutdiag.write(ln)
    foutdiagnumers.write(ln)

    N_skipped_introns = 0

    for ln in F:
        if type(ln) == bytes:
            ln = ln.decode("utf-8")  # convert bytes to string
        ln = ln.split()
        intron = ln[0]
        chrom, s, e, clu = intron.split(":")  # chr, start, end, clu_1_+
        strand = clu.split("_")[-1]

        # check if Strand is in sjc keys:
        if len(list(sjc.keys())[0]) == 4:
            intronid = chrom, int(s), int(e), strand
        elif len(list(sjc.keys())[0]) == 3:
            intronid = chrom, int(s), int(e)

        usages = [
            int(x.split("/")[0]) / (float(x.split("/")[1]) + 0.1) for x in ln[1:]
        ]  # intron usage ratios
        reads = [int(x.split("/")[0]) for x in ln[1:]]  # numerators
        sdreads = stdev(reads)  # standard deviation of read counts across samples

        # remove intron if read count SD < 0.5 and usage ratios are all 0
        if sum(usages) == 0 or sdreads < minreadstd:
            N_skipped_introns += 1
            continue

        # annotate using custom classification
        if intronid in sjc:
            classification = sjc[intronid]["SJClass"]
        else:
            classification = "IN"  # IN: INtergenic

        # add class flag and write to *_perind.noise_by_intron.gz, eg: 'chr1:825552:829002:clu_1_+:F 1/14 0/25 1/33 1/14 1/33'
        foutdiag.write(intron + f":{classification}" + " " + " ".join(ln[1:]) + "\n")
        foutdiagnumers.write(
            intron
            + f":{classification}"
            + " "
            + " ".join([str(x) for x in reads])
            + "\n"
        )

    foutdiag.close()
    sys.stderr.write(
        f"Filtered out {N_skipped_introns} introns with SD < {minreadstd} or zero usage.\n"
    )
    sys.stderr.write(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Annotation done.\n")


def main(options, libl):

    if options.cluster == None:
        pool_junc_reads(libl, options)
        refine_clusters(options)
        addlowusage(options)

    sort_junctions(libl, options)
    merge_junctions(options)
    get_numers(options)

    if not options.const:
        perind_file = f"{options.rundir}/{options.outprefix}_perind.counts.gz"
    else:
        perind_file = f"{options.rundir}/{options.outprefix}_perind.constcounts.gz"

    if options.annot != None and options.genome != None:

        sys.stdout.write(f"Loading genome {options.genome} ...")
        fa = pyfastx.Fasta(options.genome)
        sys.stdout.write("done!\n")

        if options.backward:
            sys.stdout.write("Using backward splicing junction classifier...\n")
            sjcb.ClassifySpliceJunction(
                perind_file=perind_file,
                gtf_annot=options.annot,
                fa = fa,
                rundir=options.rundir,
                outprefix=options.outprefix,
                verbose=options.verbose,
            )
        else:
            sys.stdout.write("Using forward splicing junction classifier...\n")
            sjcf.ClassifySpliceJunction(
                perind_file=perind_file,
                gtf_annot=options.annot,
                fa = fa,
                rundir=options.rundir,
                outprefix=options.outprefix,
                verbose=options.verbose,
            )
        annotate_noisy(options)

    else:
        sys.stderr.write("Skipping annotation step...\n")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-j",
        "--juncfiles",
        dest="juncfiles",
        type=str,
        required=True,
        help="text file with all junction files to be processed",
    )

    parser.add_argument(
        "-o",
        "--outprefix",
        dest="outprefix",
        default="leafcutter",
        help="output prefix (default leafcutter)",
    )

    parser.add_argument(
        "-q",
        "--quiet",
        dest="verbose",
        default=True,
        action="store_false",
        help="don't print status messages to stdout",
    )

    parser.add_argument(
        "-r",
        "--rundir",
        dest="rundir",
        default="./",
        help="write to directory (default ./)",
    )

    parser.add_argument(
        "-l",
        "--maxintronlen",
        dest="maxintronlen",
        default=100000,
        help="maximum intron length in bp (default 100,000bp)",
    )

    parser.add_argument(
        "-m",
        "--minclureads",
        dest="minclureads",
        default=30,
        help="minimum reads in a cluster (default 30 reads)",
    )

    parser.add_argument(
        "-M",
        "--minreads",
        dest="minreads",
        default=5,
        help="minimum reads for a junction to be considered for \
                        clustering(default 5 reads)",
    )

    parser.add_argument(
        "-D",
        "--minreadstd",
        dest="minreadstd",
        default=0.5,
        help="minimum standard deviation of reads across samples for a \
                        junction to be included in output (default 0.5)",
    )

    parser.add_argument(
        "-p",
        "--mincluratio",
        dest="mincluratio",
        default=0.001,
        help="minimum fraction of reads in a cluster that support a junction \
              (default 0.001)",
    )

    parser.add_argument(
        "-c",
        "--cluster",
        dest="cluster",
        default=None,
        help="refined cluster file when clusters are already made",
    )

    parser.add_argument(
        "-k",
        "--checkchrom",
        dest="checkchrom",
        action="store_true",
        default=False,
        help="check that the chromosomes are well formated e.g. chr1, chr2, \
              ..., or 1, 2, ...",
    )

    parser.add_argument(
        "-C",
        "--includeconst",
        dest="const",
        action="store_true",
        default=False,
        help="also include constitutive introns",
    )

    parser.add_argument(
        "-A",
        "--annot",
        dest="annot",
        default=None,
        help="Gencode GTF annotation file, e.g. gencode.v37.annotation.gtf.gz",
    )

    parser.add_argument(
        "-G",
        "--genome",
        dest="genome",
        default=None,
        help="Genome fasta file, e.g. hg38.fa",
    )

    # parser.add_argument("-N", "--noise", dest="noiseclass", default = None,
    #     help="Use provided intron_class.txt.gz to help identify noisy junction. \
    #           Multiple annotation files delimited by `,`. ")

    parser.add_argument(
        "-f",
        "--offset",
        dest="offset",
        default=0,
        help="Offset sometimes useful for off by 1 annotations. (default 0)",
    )

    parser.add_argument(
        "-T",
        "--keeptemp",
        dest="keeptemp",
        action="store_true",
        default=False,
        help="keep temporary files. (default false)",
    )

    parser.add_argument(
        "-B",
        "--Backward",
        dest="backward",
        action="store_true",
        default=False,
        help="Use backword splicing junction classifier. (default false)"
    )

    options = parser.parse_args()

    if options.juncfiles == None:
        sys.stderr.write("Error: no junction file provided...\n")
        exit(0)

    # Get the junction file list
    libl = []
    for junc in open(options.juncfiles):
        junc = junc.strip()
        try:
            open(junc)
        except:
            sys.stderr.write(f"{junc} does not exist... check your junction files.\n")
            exit(0)
        libl.append(junc)

    chromLst = (
        [f"chr{x}" for x in range(1, 23)]
        + ["chrX", "chrY"]
        + [f"{x}" for x in range(1, 23)]
        + ["X", "Y"]
    )

    main(options, libl)

    if not options.keeptemp:
        sys.stderr.write("Remove generated temp files... \n")
        with open(os.path.join(options.rundir, options.outprefix) + "_sortedlibs") as f:
            for tmp in [ln.strip() for ln in f.readlines()]:
                os.remove(tmp)
        os.remove(os.path.join(options.rundir, options.outprefix) + "_sortedlibs")
        sys.stderr.write(f"\n{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Done.\n")

