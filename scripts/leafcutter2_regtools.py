#!/usr/bin/env python

'''Identify noisy splicing

This script takes in junction files produced by regtools, then construct intron 
clusters. Next, constructed intron clusters are processed to identify rarely
spliced introns. Such rarely spliced introns are based on certain filtering
cut-offs. The result is a text file, following the same format as the standard 
leafcutter tool. The first column of each row indicates genome coordinates of
introns and the labeling of N[: noisy] or F[: Functional].

Junction files are processed using regtools:

    * https://github.com/griffithlab/regtools
    * /home/yangili1/tools/regtools/build/regtools junctions extract -a 8 -i \
      50 -I 500000 bamfile.bam -o outfile.junc
    * Using regtools speeds up the junction extraction step by an order of 
      magnitude or more

Procedure sequence: 

    Steps 1-3 may be skipped if running with provided cluster file. (-c path/to/clusterfile)
    Otherwise, the script generates cluster file that includes noisy introns.

    1. pool_junc_reads:
            Pool introns from a list of junction files. 
            Output text file - file name: {out_prefix}_pooled.
            Each line stores a cluster of introns: 
                e.g. chr1:+ 779093:803918:3 793042:795469:4
            Each cluster include both linked (sharing splice sites) and unlinked 
            introns.
    
    2. refine_clusters:
            Take the pooled junction reads from the 1st procedure and refine it.
            Output text file - file name: {out_prefix}_refined.
            such that filters such as minimal junction reads, minimal cluster
            reads, or ratios are met. Importantly, it ensures introns within 
            a cluster must be linked - sharing either a 5' or a 3' splice site.
            Each line is a linked intron cluster, eg:
                chr1:- 5754:5878:5
                chr1:- 5810:6470:405 5890:6470:9 6173:6470:150

    3. addlowusage:
            Step 2 produces refined clusters that pass filters. To detect noisy
            splicing, we add lowly expressed introns back into these filtered
            intron clusters. Step 3 takes input from step 2, {output_prefix}_refined, 
            then add back noisy introns to each cluster. Note that each cluster
            would still pass cluster reads filter.
            Output 2 sorted text file:
            - {out_prefix}_lowusage_introns: noisy introns only (intermediate) 
            - {out_prefix}_refined_noisy: refined clusters including noisy introns (final output)

    4. sort_junctions
            Sort junction files using step 3 {out_prefix}_refined_noisy file, or if -c specified, use
            provided cluster file to sort input junctions following the same order in cluster file.
            Out a sorted junction file {out_prefix}_{libName}.junc.sorted.gz for each library in junc list.
    
    5. merge_junctions
            Merge sorted junc files, specified in {out_prefix}_sortedlibs, aka. 
            each {out_prefix}_{libName}.junc.sorted.gz. The result is a count table
            with rows are introns, columns are counts for each library.
            Write out: {out_prefix}_perind.counts.gz

    6. get_numers
            Write out the numerators from {out_prefix}_perind.counts.gz. 
            Output file: {out_prefix}_perind_numers.counts.gz
    
    7. annotate_noisy
            Annotate {out_prefix}_perind.counts.gz. Introns are marked with noisy vers functional. 
            `*` appended coordinates indicate noisy. 
            Output several files:
                - {out_prefix}_perind.counts.noise.gz: output functional introns (intact), and 
                  noisy introns. Note the start and end coordinates of noisy introns are recalibrated
                  to the min(starts) and max(ends) of all functional introns within cluster.
                - {out_prefix}_perind_numers.counts.noise.gz: same as above, except write numerators.
                - {out_prefix}_perind.counts.noise_by_intron.gz: same as the first output, except here
                  noisy introns' coordinates are kept as their original coordinates. 

NOTE: 
    * Minimum version requirement - python v3.6

'''

import sys
import tempfile
import os
import gzip
import shutil
from statistics import mean, median, stdev
import pickle
import argparse
import pandas as pd
from Bio.Seq import Seq
import pyfastx



__author__    = "Yang Li, Chao Dai"
__email__     = "chaodai@uchicago.edu"
__status__    = "Development"
__version__   =  "v1.0.0"


def natural_sort(l: list): 
    '''Natural sort a list of string/tuple, similar to bash `sort -V`
    
    Parameters:
    -----------
    l : list
        l can be a list of string or numerics; or a list of varing length of tuples
    
    Returns:
    --------
    return : a sorted list
    '''
    import re
    
    untuple = lambda tup: ''.join([str(e) for e in tup])
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', untuple(key))]
    return sorted(l, key=alphanum_key)


def cluster_intervals(E: list):
    '''Clusters intervals together
    
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
    '''
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
    '''Checks if A and B overlaps
    
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
    '''

    if A[1] < B[0] or B[1] < A[0]:
        return False
    else: return True


def pool_junc_reads(flist, options):
    '''Pool junction reads

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
    '''

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
    by_chrom = {} # { k=(chrom, strand) : v={ k=(start, end) : v=reads } }

    for libl in flist:
        lib = libl.strip()
        if not os.path.isfile(lib):
            continue

        if options.verbose:
            sys.stderr.write(f"scanning {lib}...\n")

        if ".gz" in lib:
            F = gzip.open(lib)
        else:
            F = open(lib)

        for ln in F:

            if type(ln) == bytes:
                ln = ln.decode('utf-8') # convert bytes to string
            
            lnsplit=ln.split()
            if len(lnsplit) < 6: 
                sys.stderr.write(f"Error in {lib} \n")
                continue

            if len(lnsplit) == 12: # 12 fields regtools junc file
                chrom, A, B, dot, counts, strand, rA,rb, rgb, blockCount, blockSize, blockStarts = lnsplit
                if int(blockCount) > 2:  
                    print(ln, "ignored...")
                    continue
                Aoff, Boff = blockSize.split(",")[:2]
                A, B = int(A)+int(Aoff), int(B)-int(Boff) # get intron

            elif len(lnsplit) == 6:
                # old leafcutter junctions
                chrom, A, B, dot, counts, strand = lnsplit
                A, B = int(A), int(B)
            
            if checkchrom and (chrom not in chromLst): 
                continue

            A, B = int(A), int(B)+int(options.offset)

            if B-A > int(maxIntronLen): continue
            
            # sum up all the reads at the same junctions if junctions already exist
            try: by_chrom[(chrom,strand)][(A,B)] = int(counts) + \
                    by_chrom[(chrom,strand)][(A,B)]                                                                  
            except:
                try: by_chrom[(chrom,strand)][(A,B)] = int(counts) # when only 1 junction
                except: by_chrom[(chrom,strand)] = {(A,B):int(counts)} # when only 1 junction


    with open(outFile, 'w') as fout:
        Ncluster = 0
        sys.stderr.write("Parsing...\n")
    
        for chrom in by_chrom:
            read_ks = [k for k,v in by_chrom[chrom].items() if v >= 3] # read-keys, require junction reads > 3
            read_ks.sort() # sort read-keys: (start, end)
        
            sys.stderr.write(f"{chrom[0]}:{chrom[1]}..\n")
            
            if read_ks:
                clu = cluster_intervals(read_ks)[0] # clusters of introns, [[(start, end),..],..]
                
                for cl in clu:
                    if len(cl) > 0: # 1 if cluster has more than one intron  
                        buf = f'{chrom[0]}:{chrom[1]} ' # chr:strand
                        for interval, count in [(x, by_chrom[chrom][x]) for x in cl]:
                            buf += f"{interval[0]}:{interval[1]}" + f":{count}" + " " # start:end:reads
                        fout.write(buf+'\n')
                        Ncluster += 1
                        
        sys.stderr.write(f"\nWrote {Ncluster} clusters..\n")


def refine_linked(clusters):
    '''Re-cluster introns into clusters of linked introns

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
    '''

    unassigned = [x for x in clusters[1:]]
    current = [clusters[0]]
    splicesites = set([current[0][0][0],current[0][0][1]]) # start & end of intron
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
            splicesites = set([current[0][0][0],current[0][0][1]])
            if len(unassigned) == 1: # assign the last intron
                newClusters.append(current)
            unassigned = unassigned[1:]
    return newClusters


def refine_cluster(clu: list, cutoff: float, readcutoff: int):
    '''Filter introns based on cutoffs

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

    '''
    
    remove = []
    dic = {}
    intervals = []

    reCLU = False # re-cluster flag
    totN = 0

    for inter, count in clu:
        totN += count
    for inter, count in clu:
        if (count/float(totN) >= cutoff and count >= readcutoff):
            intervals.append(inter)
            dic[inter] = count # {(start, end): reads}
        else:
            reCLU = True # any intron not passing filters will enforce reCLU

    
    if len(intervals) == 0: return [] # base case
    
    # Below makes sure that after trimming/filtering, the clusters are still good
    # afterwards - each clusters have linked introns that pass filters.

    Atmp, _ = cluster_intervals(intervals)
    A = []

    # A: a list of linked intron clusters
    if len(Atmp) > 0:
        for cl in Atmp: # Atmp is a list of list
            if len(cl) == 1: # 1 intron
                A.append(cl)
            for c in refine_linked([(x,0) for x in cl]): # >1 introns
                if len(c) > 0:
                    A.append([x[0] for x in c])
            

    if len(A) == 1: # A has 1 cluster of introns
        rc = [(x, dic[x]) for x in A[0]]
    
        if len(rc) > 0:
            if reCLU: # recompute because ratio changed after removal of some introns
                return refine_cluster([(x, dic[x]) for x in A[0]], cutoff, readcutoff) # recursive
            else:
                return [[(x, dic[x]) for x in A[0]]]
    
    NCs = [] # As in N Clusters, here A has more than 1 clusters of introns
    for c in A: 
        if len(c) > 1: # c has more than 1 introns
            NC = refine_cluster([(x, dic[x]) for x in c], cutoff, readcutoff) # recursive
            NCs += NC
    
    return NCs 



def refine_clusters(options):
    '''Refine clusters.
    
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

    '''

    outPrefix = options.outprefix
    rundir = options.rundir
    minratio = float(options.mincluratio)
    minclureads = int(options.minclureads)
    minreads = int(options.minreads)

    inFile = f"{rundir}/{outPrefix}_pooled"
    outFile = f"{rundir}/{outPrefix}_refined"
    fout = open(outFile,'w')

    sys.stderr.write(f"\nRefine clusters from {inFile}...\n")

    Ncl = 0
    for ln in open(inFile): # pooled juncs
        if type(ln) == bytes:
            ln = ln.decode('utf-8') # convert bytes to string
        
        clu = [] # each cluster: [((start, end), reads),..]
        totN = 0 # total cluster reads
        chrom = ln.split()[0]
        for ex in ln.split()[1:]: # for an exon
            A, B, N = ex.split(":")
            clu.append(((int(A),int(B)), int(N)))
            totN += int(N)
        
        if totN < minclureads: continue

        if options.const: # include constitutive introns. These are clusters that only have 1 intron, hence "constitutive"
            if len(clu) == 1:
                buf = f'{chrom} '
                for interval, count in clu:
                    buf += f"{interval[0]}:{interval[1]}" + f":{count}" + " "
                Ncl += 1
                fout.write(buf+'\n') # e.g. 'chr:strand start:end:reads'
        
        for cl in refine_linked(clu): # only linked intron clusters
            rc = refine_cluster(cl, minratio, minreads)
            if len(rc) > 0:
                for clu in rc:
                    buf = f'{chrom} '
                    for interval, count in clu:
                        buf += f"{interval[0]}:{interval[1]}" + f":{count}" + " "
                    Ncl += 1
                    fout.write(buf+'\n')
    sys.stderr.write(f"Split into {Ncl} clusters...\n")
    fout.close()



def addlowusage(options):
    '''Add low usage introns to refined clusters
    
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
              
    '''

    global chromLst

    sys.stderr.write("\nAdd low usage introns...\n")

    outPrefix = options.outprefix
    rundir = options.rundir
    pooled = f"{rundir}/{outPrefix}_pooled"
    minclureads = int(options.minclureads)
    minreads = int(options.minreads)

    if options.cluster == None:
        refined_cluster = f"{rundir}/{outPrefix}_refined"
    else:
        refined_cluster = options.cluster

    outFile = f"{rundir}/{outPrefix}_refined_noisy" # out file that includes noisy introns
    outFile_lowusageintrons = f"{rundir}/{outPrefix}_lowusage_introns" # out file for lowusage introns

    fout = open(outFile,'w')
    fout_lowusage = open(outFile_lowusageintrons,'w')
    

    # get 5' sites, 5' sites, and clusters of introns from refined file, see data structure below
    exons5,exons3, cluExons = {}, {}, {}
    cluN = 0 # clusterID

    # construct 5' sites, 3' sites, and clusters dict from refined
    for ln in open(refined_cluster):
        if type(ln) == bytes:
            ln = ln.decode('utf-8') # convert bytes to string

        chrom = ln.split()[0]
        cluN += 1
        cluExons[(chrom,cluN)] = [] # keys are (chrom, cluster_number)
        for exon in ln.split()[1:]:
            A, B, count = exon.split(":") # start, end, reads
            if chrom not in exons5:
                exons5[chrom] = {}
                exons3[chrom] = {}
            exons5[chrom][int(A)] = (chrom, cluN) # 5' sites, { k=chrom, v={ k=start, v=(chrom, clusterID) } }
            exons3[chrom][int(B)] = (chrom, cluN) # 3' sites, { k=chrom, v={ k=end, v=(chrom, clusterID) } }
            cluExons[(chrom, cluN)].append(exon) # introns, { k=(chrom, clusterID), v=['start:end:reads'] }

    
    # Below for loop adds back clusters (previously filtered out in refined_clusters)
    # in the pooled junc file. These previously removed introns are added to all 
    # cluExons, as well as to lowusage_intron
    lowusage_intron = {} # { k=(chrom, clusterID), v=['start:end:reads'...]}
    for ln in open(pooled): # each cluster/line from pool_juncs file

        if type(ln) == bytes:
            ln = ln.decode('utf-8') # convert bytes to string

        clu = []
        totN = 0
        chrom = ln.split()[0]

        if chrom in exons5: # ensure chrom is in exons5 level-1 keys

            # Below for loop adds introns that were filtered out in refined, aka noisy_introns, 
            # back to a total intron cluster dict and to a lowusage (noisy) intron cluster dict
            for exon in ln.split()[1:]:
                A, B, N = exon.split(":") # start, end, reads

                # when 5' site in refined
                if int(A) in exons5[chrom]:
                    clu = exons5[chrom][int(A)] # set clu=(chrom, clusterID), key for cluExons
                    if exon not in cluExons[clu]: # exon was filtered out by refined
                        cluExons[clu].append(exon) # add it to cluExons
                        if clu not in lowusage_intron:
                            lowusage_intron[clu] = []
                        lowusage_intron[clu].append(exon) # also add it to lowusage
                
                # else when 3' site in refined, perform same procedure
                elif int(B) in exons3[chrom]: # when 3' site is in refined
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
    ks = natural_sort(lowusage_intron.keys()) # e.g. { k=(chrom, clusterID), v=['start:end:reads'...]}
    for clu in ks: # e.g. (chrom, clusterID)
        fout_lowusage.write(clu[0] + " " + " ".join(lowusage_intron[clu])+'\n')
    fout_lowusage.close()

    # write all intron clusters
    cluLst = natural_sort(cluExons.keys())
    for clu in cluLst:
        if not options.const: # if -C flag not set, do not write constitutive introns
            if len(cluExons[clu]) == 1: continue # skip write out if only 1 intron in cluster, aka, constitutive

        # only write introns if minimum cluster reads criteria is met
        if sum([int(ex.split(":")[-1]) for ex in cluExons[clu]]) < minclureads:
            continue
        chrom = clu[0]
        buf = f'{chrom}'
        for ex in cluExons[clu]:
            buf += " " + ex
        fout.write(buf+'\n')
    fout.close()


def sort_junctions(libl, options):
    '''Sort junctions by cluster

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
    '''
    
    global chromLst


    outPrefix = options.outprefix
    rundir = options.rundir
    checkchrom = options.checkchrom

    if options.cluster == None: # if not providing refined clusters externally
        refined_cluster = f"{rundir}/{outPrefix}_refined_noisy" # note refined noisy intron clusters
        sys.stderr.write(f"\nUsing {refined_cluster} as refined cluster...\n")
    else:
        refined_cluster = options.cluster

    runName = f"{rundir}/{outPrefix}"

    # exons:  { k=chrom : v={ k=(start, end) : v=clusterID } }
    # cluExons:  { k=clusterID : v=[(chrom, start, end)] }
    exons, cluExons = {}, {} 
    cluN = 0 # clusterID
    
    # fill in exons, cluExons dict from `*refined_noisy` intron cluster file
    for ln in open(refined_cluster): # e.g. ln = "chr10:+ 135203:179993:5 135302:179993:29"

        if type(ln) == bytes:
            ln = ln.decode('utf-8') # convert bytes to string

        chrom = ln.split()[0] # e.g. "chr10:+"
        cluN += 1
        for exon in ln.split()[1:]: # e.g. "135203:179993:5 135302:179993:29"
            A, B, count = exon.split(":")
            
            if chrom not in exons:
                exons[chrom] = {}
            if (int(A),int(B)) not in exons[chrom]:
                exons[chrom][(int(A),int(B))] = [cluN]
            else:
                exons[chrom][(int(A),int(B))].append(cluN)
            if cluN not in cluExons:
                cluExons[cluN] = []
            cluExons[cluN].append((chrom, A, B))

    merges = {} # stores junc file names as dict { k=filename : v=[filename] }
    for ll in libl:
        lib = ll.rstrip() # 1 junc file path
        libN = lib.split('/')[-1].split('.')[0] # get library name from junc file name eg. GTEX-1117F-0226-SM-5GZZ7
        if not os.path.isfile(lib):
            continue
        if libN not in merges:
            merges[libN] = [] # why use list, `libN` should always be one element
        merges[libN].append(lib)

    fout_runlibs = open(os.path.join(rundir, outPrefix) + '_sortedlibs', 'w') # intermediate file to store sorted junc file names to be written

    # operate on each libN (library), each libN can have more than 1+ junc files
    for libN in merges: 
        by_chrom = {} # to store junctions from original unsorted junc file
        
        
        # write sorted junc file names into intermediate file
        foutName = os.path.join(rundir, outPrefix + '_' + libN + '.junc.sorted.gz') # 'test/gtex_w_clu/gtex_GTEX-1IDJU-0006-SM-CMKFK.junc.sorted.gz'
        fout_runlibs.write(foutName + '\n') # e.g. 'test/gtex_w_clu/gtex_sortedlibs'

        if options.verbose:   
            sys.stderr.write(f"Sorting {libN}..\n")
        if len(merges[libN]) > 1: 
            if options.verbose:   
                sys.stderr.write(f"merging {' '.join(merges[libN])}...\n")
        else: pass
        fout = gzip.open(foutName,'wt') # e.g. 'test/gtex_w_clu/gtex_GTEX-1IDJU-0006-SM-CMKFK.junc.sorted.gz'

        #-------- Process and write junction files --------
        
        # write header
        fout.write(f'chrom {libN}\n') # 'chrom GTEX-111VG-0526-SM-5N9BW\n'


        #-------- Gather counts from all junc files of library --------
        # store in by_chrom: { ('chr1', '+') : { (100, 300) : 5, (500, 700): 10, ... } }
        for lib in merges[libN]:
            if ".gz" in lib: 
                F = gzip.open(lib)
            else: 
                F = open(lib)
        
            for ln in F: # 1 line: e.g. "chr17\t81701131\t81701534\t.\t1\t+"

                if type(ln) == bytes:
                    ln = ln.decode('utf-8') # convert bytes to string
                lnsplit = ln.split()

                if len(lnsplit) < 6:
                    sys.stderr.write(f"Error in {lib} \n")
                    continue

                if len(lnsplit) == 12:
                    chrom, A, B, dot, counts, strand, rA,rb, rgb, blockCount, blockSize, blockStarts = lnsplit
                    if int(blockCount) > 2:
                        print(ln, "ignored...")
                        continue
                    Aoff, Boff = blockSize.split(",")[:2]
                    A, B = int(A)+int(Aoff), int(B)-int(Boff)

                elif len(lnsplit) == 6:
                    # old leafcutter junctions                                                                                                                       
                    chrom, A, B, dot, counts, strand = lnsplit
                    A, B = int(A), int(B)
                    
                A, B = int(A), int(B) + int(options.offset) # start, end + offset
            
                chrom = (chrom, strand)
                if chrom not in by_chrom: 
                    by_chrom[chrom] = {} # store introns from junc file, key: ('chr1', '+')
                
                intron = (A, B)
                if intron in by_chrom[chrom]: # sum up reads by intron from junc files
                    by_chrom[chrom][intron] += int(counts)
                else:
                    by_chrom[chrom][intron] = int(counts)
        
        
        #------- Take clusters from refined_noisy, assign reads -------
        # reads are from by_chrom (junc files)
        # For each intron cluster, write fraction for each intron (one intron per line).
        for clu in cluExons: # cluExons: { k=cluID : v=[(chrom, start, end)...]}
            buf = []
            ks = cluExons[clu] # eg: [('chr1:+', 827776, 829002), ..] introns of a clu
            ks.sort() # no need to version sort within cluster

            # Step 1: sum cluster level reads from each intron
            # gather (sum) reads for each cluster in refined_noisy, read counts are from junc file (by_chrom)
            tot = 0 # sum of total read counts per cluster
            usages = []
            for exon in ks:
                chrom, start, end = exon
                chrom = tuple(chrom.split(":")) # convert 'chr3:+' to ('chr3', '+') as in by_chrom
                start, end = int(start), int(end)

                if chrom not in by_chrom:
                    pass
                elif (start, end) in by_chrom[chrom]:
                    tot += by_chrom[chrom][(start,end)]

            # Step 2: append intron usage fraction to stream buffer
            for exon in ks:
                chrom, start, end = exon
                start, end = int(start), int(end)
                chrom = tuple(chrom.split(":"))
                chromID, strand = chrom # chromID eg: 'chr3'

                intron = chromID, start, end+1, strand # converting to 1-based coordinates

                if chrom not in by_chrom: 
                    # if refined exon chrom is not found in junc file, write 0/cluster_total
                    buf.append(f"{chromID}:{start}:{end}:clu_{clu}_{strand} 0/{tot}\n")
                elif (start,end) in by_chrom[chrom]: 
                    # if refind exon is in junc file, write exon reads / cluster_total
                    buf.append(f"{chromID}:{start}:{end}:clu_{clu}_{strand} {by_chrom[chrom][(start,end)]}/{tot}\n")
                else:
                    # if refined exon is not found in junc file, write 0/cluster_total
                    buf.append(f"{chromID}:{start}:{end}:clu_{clu}_{strand} 0/{tot}\n")
        
            fout.write("".join(buf))
        fout.close()
    fout_runlibs.close()


def merge_files(fnames, fout, options):
    '''Merge a list of files into a gzip file

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
    '''

    fopen = []
    for fname in fnames: # each file to be merged
        if fname[-3:] == ".gz":
            fopen.append(gzip.open(fname))
        else:
            fopen.append(open(fname))

    finished = False
    N = 0
    while not finished: # cycle through files in batch
        N += 1
        if N % 50000 == 0: 
            sys.stderr.write(".")
        buf = []
        for f in fopen: # each opened file
            ln = f.readline().decode().split() # read 1 line
            if len(ln) == 0: # end of line = finish
                finished = True
                break
            chrom = ln[0] # e.g. "chrom" or "chr1:825552:829002:clu_1_+"
            data = ln[1:] # e.g. "GTEX-1117F-0626-SM-5N9CS.leafcutter" or "0/0"
            if len(buf) == 0:
                buf.append(chrom)
            buf += data # e.g. ['chrom', 'GTEX-111VG-0526-SM-5N9BW.leafcutter', 'GTEX-1117F-0626-SM-5N9CS.leafcutter'] for first lines, or ['chr1:825552:829002:clu_1_+', '0/0', '0/0'] for 2+ lines
            # each file the exact same chromosome coordinates, effectively we are collecting counts into columns 2 and after
        
        if len(buf) > 0:
            if buf[0] == "chrom":
                if options.verbose:
                    sys.stderr.write(f"merging {len(buf)-1} files")
            fout.write(" ".join(buf)+'\n') # combining sample counts into columns
        else:
            break

    sys.stderr.write(" done.\n")
    for fin in fopen:
        fin.close()


def merge_junctions(options):    
    '''Merge junctions

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
    '''

    outPrefix = options.outprefix
    rundir = options.rundir
    
    fnameout = os.path.join(f'{rundir}/{outPrefix}')
    flist = fnameout + '_sortedlibs' # sorted juncs file list
    lsts = [] # = flist
    for ln in open(flist):
        if type(ln) == bytes:
            ln = ln.decode('utf-8') # convert bytes to string
        lsts.append(ln.strip())

    if options.verbose:
        sys.stderr.write(f"\nMerging {len(lsts)} junction files...\n")
    
    # Change 300 if max open file is < 300
    # set up batch N per batch
    N = min([300, max([100, int(len(lsts)**(0.5))])])

    # tmpfiles = []
    while len(lsts) > 1: # initial list of sorted junc files
        
        # convert lsts (list of file paths) to clst (list of lists)
        # each sublist is a batch of upto 100 files.
        clst = [] # list of batches, each batch has up to 100 sorted junc files
        for i in range(0, int(len(lsts)/N)+1): # merge in batches of max(100, len(lsts))
            lst = lsts[N*i:N*(i+1)]
            if len(lst) > 0:
                clst.append(lst)
        lsts = [] # clear initial file list, now repurposed to store merged file names (temp)
    
        for lst in clst: # run by batch
            if len(lst) == 0: 
                continue
            tmpfile = tempfile.mktemp()
            os.mkdir(tmpfile)
            foutname = tmpfile+"/tmpmerge.gz"
            fout = gzip.open(foutname,'wt') # create a temp file for the batch of files to be merged
            
            merge_files(lst, fout, options) # merge the batch into `fout`
            lsts.append(foutname) # save the temp merged file name
            # tmpfiles.append(foutname) # this line is not needed.
            fout.close()
            
    if not options.const:
        shutil.move(lsts[0], fnameout+"_perind.counts.gz") 
    else:
        shutil.move(lsts[0], fnameout+"_perind.constcounts.gz") 


def get_numers(options):
    '''Get numerators from merged count table

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

    '''

    outPrefix = options.outprefix
    rundir = options.rundir

    if not options.const:                                                                                                                                                                                                                                                                                                                           
        fname = f"{rundir}/{outPrefix}_perind.counts.gz"
        fnameout = f"{rundir}/{outPrefix}_perind_numers.counts.gz"
    else:
        fname = f"{rundir}/{outPrefix}_perind.constcounts.gz"
        fnameout = f"{rundir}/{outPrefix}_perind_numers.constcounts.gz"
    
    input_file=gzip.open(fname, 'r')
    fout = gzip.open(fnameout,'wt')
    first_line=True

    sys.stderr.write(f"\nExtracting numerators (read counts) from {fname}...")
    
    for l in input_file:
        if first_line:
            fout.write(" ".join(l.decode().strip().split(" ")[1:])+'\n') # print the sample names
            first_line=False
        else:
            l=l.decode().strip()
            words=l.split(" ")
            fout.write(words[0]+ " "+ " ".join( [ g.split("/")[0] for g in words[1:] ] ) +'\n') # write intron and numerators

    input_file.close()
    fout.close()
    sys.stderr.write(" done.\n")


def loadIntronAnnotations(fn):
    '''parse and load intron annotations from a single file'''

    dic_noise = {} # noisy intron annotation: { k=(chr, start, end, strand) : v=classfication }
    if fn[-3:] == '.gz':
        for ln in gzip.open(fn):
            chrom, s, e, strand, classification = ln.decode().split()
            # dic_noise[(chrom,int(s)-1,int(e),strand)] = classification # intron annotation
            dic_noise[(chrom,int(s),int(e),strand)] = classification # intron annotation
    else:
        for ln in open(fn):
            if type(ln) == bytes:
                ln = ln.decode('utf-8') # convert bytes to string
            chrom, s, e, strand, classification = ln.split()
            # dic_noise[(chrom,int(s)-1,int(e),strand)] = classification 
            dic_noise[(chrom,int(s),int(e),strand)] = classification 
    return dic_noise

def annotate_noisy(options):
    '''Annotate introns

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
    
    '''

    outPrefix = options.outprefix
    rundir = options.rundir
    minreadstd = float(options.minreadstd)
    fnameout = f"{rundir}/{outPrefix}"
    # noisy_annotations = options.noiseclass # intron annotation - noisy or funciton, etc. 
    sjc_f = f"{rundir}/{outPrefix}_junction_classifications.txt" # Classified junction annotations

    dic_class = {"putative_functional":"F",
                 "productive": "F",
                 "Productive": "F",
                 "putative_noisy":"PN", 
                 "noisy":"N"}
    
    sys.stderr.write(f"\nAnnotating introns with custom-classified annotations {sjc_f}...\n")
    
    if sjc_f != None:

        if options.verbose:
            sys.stderr.write(f"Loading {sjc_f} for (un/)productive splicing classification..\n")
        sjc = merge_discordant_logics(sjc_f) 
        # noisy_annotations = [x.strip() for x in noisy_annotations.split(',')]
        # dic_noise = [loadIntronAnnotations(f) for f in noisy_annotations]
        # dic_noise = {k:v for dic in dic_noise for k,v in dic.items()}

        if options.verbose:
            sys.stderr.write("Loaded..\n")


    if not options.const:
        fname =  fnameout+"_perind.counts.gz" # no constitutive introns
    else:
        fname = fnameout+"_perind.constcounts.gz"

    noisydiag = fname.replace(".gz",".noise_by_intron.gz") # eg: run/out_perind.counts.noise_by_intron.gz
    numersdiag = fname.replace(".gz",".noise_by_intron.gz").replace("perind",'perind_numers') # eg: run/out_perind_numers.counts.noise.gz

    foutdiag = gzip.open(noisydiag,'wt')
    foutdiagnumers = gzip.open(numersdiag, 'wt')


    F = gzip.open(fname)
    ln = F.readline().decode()
    foutdiag.write(ln)
    foutdiagnumers.write(ln)

    N_skipped_introns = 0

    for ln in F:
        if type(ln) == bytes:
            ln = ln.decode('utf-8') # convert bytes to string
        ln = ln.split()
        intron = ln[0]
        chrom, s, e, clu = intron.split(":") # chr, start, end, clu_1_+
        strand = clu.split("_")[-1]
    
        intronid = chrom, int(s), int(e)
        usages = [int(x.split("/")[0]) / (float(x.split("/")[1])+0.1) for x in ln[1:]] # intron usage ratios
        reads = [int(x.split("/")[0]) for x in ln[1:]] # numerators
        sdreads = stdev(reads) # standard deviation of read counts across samples

        # remove intron if read count SD < 0.5 and usage ratios are all 0
        if sum(usages) == 0 or sdreads < minreadstd:
            N_skipped_introns += 1
            continue

        # annotate using custom classification
        if intronid in sjc:
            classification = sjc[intronid]['SJClass']
        else:
            classification = "IN" # IN: INtergenic
            
        # if dic_class[classification] in ["N","PN"]:
        #     # noisy_clusters[clu] = '' 
        #     clusters[clu]["noise"].append(ln) # add intron to clusters under key eg: ['clu_1_+']['noise']
        # else:
        #     clusters[clu]["functional"].append(ln) # otherwise under key eg: ['clu_1_+]['functional']

        # add class flag and write to *_perind.noise_by_intron.gz, eg: 'chr1:825552:829002:clu_1_+:F 1/14 0/25 1/33 1/14 1/33'
        foutdiag.write(intron + f":{classification}" + ' ' + ' '.join(ln[1:]) + '\n')
        foutdiagnumers.write(intron + f":{classification}" + ' ' + ' '.join([str(x) for x in reads]) + '\n')


    foutdiag.close()
    # fout.close()
    sys.stderr.write(f"Filtered out {N_skipped_introns} introns with SD < {minreadstd} or zero usage.\n")
    sys.stderr.write(f"Annotation done.\n")





#-------------------------------------------------------------------
# functions for Splice Junction Classification using GTF and Refseq


def check_utrs(junc,utrs):
    '''
    checks if junction is close or within 100bp of UTRs
    '''
    for s1,s2 in list(utrs):
        if abs(junc[0]-s1) < 100 or abs(junc[1]-s2) < 100:
            return True
    return False

def solve_NMD(chrom, strand, junc, start_codons, stop_codons,gene_name, 
              verbose = False, exonLcutoff = 1000):
    '''
    Compute whether there is a possible combination that uses the junction without
    inducing a PTC. We start with all annotated stop codon and go backwards.
    '''

    global fa
    
    seed = []

    junc.sort()
    if strand == "+":
        junc.reverse()
        
    for c in stop_codons:
        if strand == "+":
            seed.append([c[1]])
        else:
            seed.append([c[0]])

    # seed starts with just stop codon and then a possible 3'ss-5'ss junction
    # without introducing a PTC [stop_codon,3'ss, 5'ss, 3'ss, ..., start_codon]
    
    seq_db = {}
    junc_pass = {}
    junc_fail = {}
    path_pass = []
    proteins = []
    
    dic_terminus = {}
    dic_paths = {}

    depth = 0
    while len(seed) > 0:
        new_seed = []
        final_check = []
        depth += 1
        if verbose:
            sys.stdout.write("Depth %s, Seed L = %s\n"%(depth, len(seed)))
        #print(start_codons, [s[-1] for s in seed][-10:], len(junc))
        framepos = {}
                    
        for s in seed:
            # first check that the seed paths are good        
            bool_ptc = False
            leftover = ''
            if len(s) > 0:                
                leftover = Seq("")
                allprot = Seq("")
                for i in range(0, len(s)-1, 2):
                    exon_coord = s[i:i+2]
                    exon_coord.sort()
                    exon_coord = tuple(exon_coord)
                    exlen = exon_coord[1]-exon_coord[0]

                    startpos = (len(leftover)+exlen+1)%3
                    if strand == '+':
                        seq = Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))+leftover 
                        prot = seq[startpos:].translate()
                        leftover = seq[:startpos]                                                                                                               
                        allprot = prot+allprot  
                    else:
                        seq = leftover+Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))
                        aseq = seq
                        if startpos > 0:
                            leftover = seq[-startpos:]
                        else:
                            leftover = Seq("")
                        seq = seq.reverse_complement()
                        prot = seq[startpos:].translate()
                        allprot = prot+allprot
                    bool_ptc = "*" in allprot[:-1]
                    
            if bool_ptc:
                #This transcript failed
                for i in range(1, len(s)-1, 2):                                                                                                                  
                    j_coord = s[i:i+2]                                                                                                                           
                    j_coord.sort()                                                                                                                             
                    j_coord = tuple(j_coord)                                                                                                                     
                    if j_coord not in junc_fail:                                                                                                                 
                        junc_fail[j_coord] = 0                                                                                                                   
                    junc_fail[j_coord] += 1  

                continue
        
            # passed
            if len(s) > 2:
                terminus = (s[-2],s[-1],leftover)
                
                if terminus in dic_terminus:
                    dic_terminus[terminus].append(tuple(s))
                    continue
                else:
                    dic_terminus[terminus] = [tuple(s)]
            
            last_pos = s[-1]
            
            for start in start_codons:                
                #print("start", start, abs(last_pos-start[0]))
                if strand == "+" and last_pos > start[0] and abs(last_pos-start[0]) < exonLcutoff:
                    final_check.append(s+[start[0]])
                elif strand == "-" and last_pos < start[1] and abs(last_pos-start[1]) < exonLcutoff:
                    final_check.append(s+[start[1]]) 
            for j0,j1 in junc:                
                if strand == "+" and last_pos > j1 and abs(last_pos-j1) < exonLcutoff:
                    new_seed.append(s+[j1,j0])

                #print("junction", (j0,j1), abs(last_pos-j0))
                if strand == "-" and last_pos < j0 and abs(last_pos-j0) < exonLcutoff: 
                    new_seed.append(s+[j0,j1])
                    
        # check that the possible final paths are good
        for s in final_check:
            leftover = Seq("")
            allprot = Seq("")
            for i in range(0, len(s)-1, 2):
                exon_coord = s[i:i+2]
                exon_coord.sort()
                exon_coord = tuple(exon_coord)
                exlen = exon_coord[1]-exon_coord[0]
                startpos = (len(leftover)+exlen+1)%3
                if strand == "+":
                    seq = Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))+leftover
                    leftover = seq[:startpos]  
                    prot = seq[startpos:].translate()
                    allprot = prot+allprot
                else:
                    seq = leftover+Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))
                    if startpos > 0:                                                                                                    
                        leftover = seq[-startpos:]                                    
                    else:
                        leftover = Seq("")
                    seq = seq.reverse_complement()                                                                                                           
                    prot = seq[startpos:].translate()                                                                                                        
                    allprot = prot+allprot                    
            bool_ptc = "*" in allprot[:-1]
        
        
            if not bool_ptc:
                # all pass
                proteins.append("\t".join([gene_name,chrom,strand, "-".join([str(x) for x in s]), str(allprot)])+'\n')
                #print("ALL PASS %s"%(s))
                path_pass.append(tuple(s))
                for i in range(1, len(s), 2):
                    j_coord = s[i:i+2]
                    j_coord.sort()
                    j_coord = tuple(j_coord)
                    if j_coord not in junc_pass:
                        junc_pass[j_coord] = 0
                    junc_pass[j_coord] += 1

        seed = new_seed
    
    
    while True:
        new_paths = []
        for terminus in dic_terminus:
            terminus_pass = False
            for path_subset in dic_terminus[terminus]:
                for path in path_pass:
                    if path[:len(path_subset)] == path_subset:
                        terminus_pass = True
                        break
            #print(terminus, terminus_pass)
            if terminus_pass:
                for path_subset in dic_terminus[terminus]:
                    if path_subset in path_pass: continue
                    new_paths.append(path_subset)
                    path_pass.append(path_subset)
                    for i in range(1, len(path_subset), 2):
                        j_coord = list(path_subset[i:i+2])
                        j_coord.sort()
                        j_coord = tuple(j_coord)
                        if j_coord not in junc_pass:
                            junc_pass[j_coord] = 0
                            if verbose:
                                sys.stdout.write("junction %s pass\n"%j_coord)
        if len(new_paths) == 0:
            break
            
    return junc_pass,junc_fail,proteins

def parse_gtf(gtf: str):
    '''Lower level function to parse GTF file
    - gtf: str : path to GTF annotation file
    - returns: dictionary with keys: 
        chrom, source, type, start, end, strand, frame, info, gene_name, transcript_type, transcript_name, gene_type
    '''
    fields = ["chrom", "source", "type", "start", "end", ".","strand","frame", "info"]
    open_gtf = lambda x: gzip.open(x) if ".gz" in x else open(x)
    for ln in open_gtf(gtf):
        ln = ln.decode('ascii') if ".gz" in gtf else ln
        dic = {}
        if ln[0] == "#": continue
        ln = ln.strip().split('\t')
        for i in range(len(fields)):
            dic[fields[i]] = ln[i]

        # add 4 additional fields, parsed from info field
        for ks in ['gene_name', "transcript_type","transcript_name", "gene_type"]:
            info_fields = [{x.split()[0]: x.split()[1].replace('"', '')} 
                          for x in dic['info'].split(';') if len(x.split()) > 1]
            info_fields = {k: v for d in info_fields for k, v in d.items()}
            try: 
                dic[ks] = info_fields[ks]
            except:
                dic[ks] = None # if line is a gene, then wont have transcript info
        yield dic
         

def parse_annotation(gtf_annot: str):
    '''
    Used `parse_gtf` to first parse a gtf file, then extract and return further
    information, including gene coordinates, intron info, and splice site to 
    gene-name dictionary.
    
    - gtf_anno: str : path to GTF annotation file
    - returns: 
        - genes_coords: gene coordinates, grouped by chromosome and strand 
        - introns_info: a dictionary with these keys: junctions (ie all introns),
                        start_codon, stop_codon, utrs, pcjunctions (ie only protein coding introns) 
        - ss2gene: a dictionary with splice site as keys, and gene name as values, eg. ('chr1', 11869): 'DDX11L1'
    '''
    genes_info = {}
    introns_info = {}
    genes_coords = {}
    ss2gene = {}

    for dic in parse_gtf(gtf_annot):
        chrom = dic['chrom']
        gname = dic['gene_name']
        tname = dic['transcript_name'], gname
        anntype  = dic['type']
        if dic['type'] == 'gene': 
            dic['transcript_type'] = "gene"
        if dic['transcript_type'] == "nonsense_mediated_decay" and anntype == "stop_codon": 
            continue 
        if dic['transcript_type'] != "protein_coding" and anntype == "UTR": 
            continue
        start, end = int(dic['start']) - 1, int(dic['end']) # convert to BED format
        strand = dic['strand']
    
        if (chrom, strand) not in genes_coords:
            genes_coords[(chrom, strand)] = []

        if tname not in genes_info: # tname is (transcript_name, gene_name)
            genes_info[tname] = {'exons':[],
                                 'start_codon':set(), 
                                 'stop_codon':set(), 
                                 'utrs':set(),
                                 'type':dic['transcript_type']
                                 }
        
        if anntype in ["start_codon", "stop_codon"]:
            genes_info[tname][anntype].add((start,end))
        elif anntype in ["gene"]:
            genes_coords[(chrom,strand)].append(((start,end), gname))
        elif anntype in ['exon']:
            genes_info[tname]['exons'].append((start,end))
            # Store gene info for splice sites
            ss2gene[(chrom, int(dic['start']) - 1)] = dic['gene_name'] # BED
            ss2gene[(chrom, int(dic['end']))] = dic['gene_name'] # BED

        elif anntype in ['UTR']:
            genes_info[tname]['utrs'].add((start,end))

    for gene in genes_info: # this is actually transcript level!!! (transcript_name, gene_name)
        gene_name = gene[1]
        exons = genes_info[gene]['exons']
        exons.sort()

        pc = False # Basically use transcript_type == protein_coding to set flag
        if genes_info[gene]['type'] == "protein_coding": # AGAIN here gene = (transcript_name, gene_name)
            pc = True

        junctions = set() # all junctions/introns
        pcjunctions = set() # only protein coding junctions/introns

        for i in range(len(exons)-1):
            intron = (exons[i][1], exons[i+1][0])
            junctions.add(intron)
            if pc:
                pcjunctions.add(intron)

        if gene_name not in introns_info: # here only gene_name, does not incl. transcript_name
            introns_info[gene_name] = {'junctions':set(),
                                       'start_codon':set(),
                                       'stop_codon':set(),
                                       'utrs':set(),
                                       'pcjunctions':set(),
                                 }

        introns_info[gene_name]['start_codon'] = introns_info[gene_name]['start_codon'].union(genes_info[gene]['start_codon'])
        introns_info[gene_name]['stop_codon'] = introns_info[gene_name]['stop_codon'].union(genes_info[gene]['stop_codon'])
        introns_info[gene_name]['junctions'] = introns_info[gene_name]['junctions'].union(junctions)
        introns_info[gene_name]['utrs'] = introns_info[gene_name]['utrs'].union(genes_info[gene]['utrs'])
        introns_info[gene_name]['pcjunctions'] = introns_info[gene_name]['pcjunctions'].union(pcjunctions)

    return genes_coords, introns_info, ss2gene


def get_feature(fname: str, feature: str = "exon"):
    ss2gene = {}
    if ".gz" in fname:
        F = gzip.open(fname)
    else:
        F = open(fname)
    for ln in F:
        if ln[0] == "#": continue
        ln = ln.split('\t')
        gID = ln[-1].split('gene_name "')[1].split('"')[0]

        if ln[2] != feature:
            continue
        ss2gene[(ln[0], int(ln[3]))] = gID
        ss2gene[(ln[0], int(ln[4]))] = gID
    return ss2gene


def get_overlap_stream(L1: list, L2: list, relax = 0):
    '''                                                                                                                                                                
    L1 and L2 are sorted lists of tuples with key, values                                                                                                              
    '''
    i, j = 0, 0
    while i < len(L1) and j < len(L2):
        if L1[i][0][1] < L2[j][0][0]:
            i += 1
            continue
        elif L2[j][0][1] < L1[i][0][0]:
            j += 1
            continue
        else:
            k = 0
            # hits overlapping, check all L2 that may overlap with intron                                                                                              
            while L2[j+k][0][0] <= L1[i][0][1]:
                if overlaps(L1[i][0], L2[j+k][0]):
                    yield L1[i], L2[j+k]
                k += 1
                if j+k == len(L2): break
            i += 1


def ClassifySpliceJunction(options):
    '''
        - perind_file: str : path to counts file, e.g. leafcutter_perind.counts.gz
        - gtf_annot: str : Annotation GTF file, for example gencode.v37.annotation.gtf.gz
        - rundir: str : run directory, default is current directory
    '''

    gtf_annot, rundir, outprefix = options.annot, options.rundir, options.outprefix
    verbose = False or options.verbose
    perind_file = f"{rundir}/{outprefix}_perind.counts.gz"

    # read leafcutter perind file and store junctions in dictionary: dic_junc
    # key = (chrom,strand), value = list of junctions [(start,end)]
    dic_junc = {}
    sys.stdout.write(f"Processing junction counts {perind_file}...")
    for ln in gzip.open(perind_file):
        junc_info = ln.decode('ascii').split()[0] # first column
        if junc_info == "chrom": continue # header
        
        chrom, start, end, clu_strand = junc_info.split(":")
        strand = clu_strand.split("_")[-1]
        if (chrom,strand) not in dic_junc: 
            dic_junc[(chrom,strand)] = []
        dic_junc[(chrom,strand)].append((int(start), int(end)))

    sys.stdout.write("done!\n")
    if verbose:
        sys.stdout.write("Processed: ")
        for chrstrand in dic_junc:
            sys.stdout.write(f"{len(dic_junc[chrstrand])} jxns on {chrstrand[0]} ({chrstrand[1]}).")

    
    # load or parse gtf annotations
    # g_coords: gene coordinates, grouped by chromosome and strand
    # g_info: a dictionary with (transcript_name, gene_name) as keys, and intron info as values
    try: 
        sys.stdout.write("Loading annotations...\n")
        parsed_gtf = f"{rundir}/{gtf_annot.split('/')[-1].split('.gtf')[0]}_SJC_annotations.pckle"
        with open(parsed_gtf, 'rb') as f:
            g_coords, g_info = pickle.load(f)
        sys.stdout.write("done!\n")
    except:
        sys.stdout.write("Parsing annotations for the first time...\n")
        g_coords, g_info, ss2gene = parse_annotation(gtf_annot)
        
        for chrom,strand in g_coords:
            to_remove_gcoords = set()
            to_remove_ginfo = set()
            for gene in g_coords[(chrom,strand)]:
                if gene[1] in g_info: # check if gene_name is in introns_info keys
                    if len(g_info[gene[1]]['stop_codon']) == 0: # if there are no stop codons
                        to_remove_gcoords.add(gene) # mark gene for remove
                        to_remove_ginfo.add(gene[1]) # mark gene_name for removal

            for g in to_remove_gcoords:
                g_coords[(chrom,strand)].remove(g)
            for g in to_remove_ginfo:
                g_info.pop(g)
        sys.stdout.write("Saving parsed annotations...\n")
        with open(parsed_gtf, 'wb') as f:
            pickle.dump((g_coords, g_info), f)

    gene_juncs = {}
    for chrom,strand in dic_junc:
        if (chrom,strand) not in g_coords: 
            sys.stderr.write(f"Could not find {chrom} ({strand}) in annotations...\n")
            continue
        juncs = [(x,x) for x in dic_junc[(chrom,strand)]]
        juncs.sort()

        coords = g_coords[(chrom,strand)]
        coords.sort()
        
        # save junctions that overlapping a gene in gene_juncs dictionary: (gene_name, chrom, strand) : [junctions]
        for junc, geneinfo in get_overlap_stream(juncs,coords): 
            info = (geneinfo[1], chrom,strand)
            if info not in gene_juncs:
                gene_juncs[info] = []
            gene_juncs[info].append(junc[0])

    fout = open(f"{rundir}/{outprefix}_junction_classifications.txt",'w')
    fout.write("\t".join(["Gene_name","Intron_coord","Annot","Coding", "UTR"])+'\n')
    
    for gene_name, chrom, strand in gene_juncs:
        sys.stdout.write(f"Processing {gene_name} ({chrom}:{strand})\n")
        
        query_juncs = gene_juncs[(gene_name,chrom,strand)] # from LeafCutter perind file
        if gene_name not in g_info: continue
        junctions = g_info[gene_name]['junctions'] # from annotation
        
        # classify all junctions in gene
        junctions = list(junctions.union(query_juncs))

        start_codons = g_info[gene_name]['start_codon'] 
        stop_codons = g_info[gene_name]['stop_codon']

        if verbose:
            sys.stdout.write(f"LeafCutter junctions ({len(query_juncs)}) All junctions ({len(junctions)}) Start codons ({len(start_codons)}) Stop codons ({len(stop_codons)}) \n")

        junc_pass, junc_fail, proteins = solve_NMD(chrom,strand,junctions, 
                                                   start_codons, stop_codons, 
                                                   gene_name)
        
        for j in junctions:
            bool_pass = j in junc_pass or j in g_info[gene_name]['pcjunctions']
            bool_fail = j in junc_fail
            utr = False
            if not bool_pass:
                # Check that it's not in UTR                
                utr = check_utrs(j,g_info[gene_name]['utrs'])

            if bool_fail or bool_pass:
                tested = True
            else:
                tested = False
            annotated = j in g_info[gene_name]['junctions']
            #if not bool_pass and annotated:
            #print("%s %s %s junction: %s tested: %s utr: %s coding: %s annotated: %s "%(chrom, strand, gene_name, j, tested,utr, bool_pass, annotated))
            
            fout.write('\t'.join([gene_name, f'{chrom}:{j[0]}-{j[1]}',
                                  str(annotated), str(bool_pass), str(utr)])+'\n')
        




def boolean_to_bit(bool_vec: list):
    '''Convert boolean vector to string of "1"s and "0"s
    - bool_vec: list : boolean vector
    '''
    # Convert boolean vector to string of "1"s and "0"s
    bin_str = ''.join(['1' if b else '0' for b in bool_vec])
    
    # Convert this binary string into an integer
    # bit_num = int(bin_str, 2)
    
    return bin_str
           
            
def merge_discordant_logics(sjc_file: str):
    '''some junctions have multiple classifications. Use conservative approach
    to merge them.
    '''
    sjc = pd.read_csv(sjc_file, sep = "\t")

    classifier = {
        # each bit represents [ is annotated, is coding, is UTR ]
        '000': 'UP', # UnProductive,
        '001': 'NE', # NEither, not productive, but not considered unprod. due to close to UTR
        '010': 'PR', # PRoductive
        '011': 'PR', # PRoductive
        '100': 'UP', # UnProductive
        '101': 'PR', # PRoductive
        '110': 'PR', # PRoductive
        '111': 'PR' # Produtive
        }
    
    # group dt
    sjc = sjc.groupby('Intron_coord').agg('max').reset_index()
    # convert Annotation, Coding, UTR status to binary strings then to SJ categories
    sjc['SJClass'] = sjc.apply(lambda x: boolean_to_bit(x[2:5]), axis=1).map(classifier)
    
    # convert df to dict
    sjc = sjc.set_index('Intron_coord').to_dict(orient='index')
    sjc = {tuple([k.split(':')[0]]) + tuple(k.split(':')[1].split('-')): v for k, v in sjc.items()}
    sjc = {(k[0], int(k[1]), int(k[2])): v for k, v in sjc.items()}


    return sjc
    # sjc is a dcitionary with:
    # - keys: intron coordinates, e.g. ('chr1', 1000, 2000)
    # - values: a dictionary e.g. {'Gene_name': 'DNMBP', 'Annot': False, 'Coding': False, 'UTR': False, 'SJClass': 'UP'})

    


#-------------------------------------------------------------------

def main(options, libl):
    
    if options.cluster == None:
        pool_junc_reads(libl, options)
        refine_clusters(options)
        addlowusage(options)
    
    
    sort_junctions(libl, options)
    merge_junctions(options)
    get_numers(options)

    if options.annot != None and options.genome != None:
        global fa
        sys.stdout.write(f"Loading genome {options.genome} ...\n")
        fa = pyfastx.Fasta(options.genome)
        ClassifySpliceJunction(options)
        annotate_noisy(options)
    
    else: 
        sys.stderr.write("Skipping annotation step...\n")



#-------------------------------------------------------------------

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-j", "--juncfiles", dest="juncfiles",
        type=str, required=True,
        help="text file with all junction files to be processed")

    parser.add_argument("-o", "--outprefix", dest="outprefix",
        default = 'leafcutter',
        help="output prefix (default leafcutter)")

    parser.add_argument("-q", "--quiet", dest="verbose", default=True,
        action="store_false", help="don't print status messages to stdout")

    parser.add_argument("-r", "--rundir", dest="rundir", default='./',
                  help="write to directory (default ./)")
    
    parser.add_argument("-l", "--maxintronlen", dest="maxintronlen",
        default = 100000, 
        help="maximum intron length in bp (default 100,000bp)")

    parser.add_argument("-m", "--minclureads", dest="minclureads", default = 30,
        help="minimum reads in a cluster (default 30 reads)")

    parser.add_argument("-M", "--minreads", dest="minreads", default = 5,
                  help="minimum reads for a junction to be considered for \
                        clustering(default 5 reads)")

    parser.add_argument("-D", "--minreadstd", dest="minreadstd", default = 0.5,
                  help="minimum standard deviation of reads across samples for a \
                        junction to be included in output (default 0.5)")

    parser.add_argument("-p", "--mincluratio", dest="mincluratio", 
        default = 0.001,
        help="minimum fraction of reads in a cluster that support a junction \
              (default 0.001)")

    parser.add_argument("-c", "--cluster", dest="cluster", default = None,
        help="refined cluster file when clusters are already made")

    parser.add_argument("-k", "--checkchrom", dest="checkchrom",
        action="store_true",default = False,
        help="check that the chromosomes are well formated e.g. chr1, chr2, \
              ..., or 1, 2, ...")
    
    parser.add_argument("-C", "--includeconst", dest="const", \
        action="store_true", default = False, 
        help="also include constitutive introns")
    
    parser.add_argument("-A", "--annot", dest="annot", default = None,
        help="Gencode GTF annotation file, e.g. gencode.v37.annotation.gtf.gz")
    
    parser.add_argument("-G", "--genome", dest="genome", default = None,
        help="Genome fasta file, e.g. hg38.fa")

    # parser.add_argument("-N", "--noise", dest="noiseclass", default = None,
    #     help="Use provided intron_class.txt.gz to help identify noisy junction. \
    #           Multiple annotation files delimited by `,`. ")

    parser.add_argument("-f", "--offset", dest="offset", default = 0,
        help="Offset sometimes useful for off by 1 annotations. (default 0)")
    
    parser.add_argument("-T", "--keeptemp", dest="keeptemp", \
        action="store_true", default = False, 
        help="keep temporary files. (default false)")

    options = parser.parse_args()

    if options.juncfiles == None:
        sys.stderr.write("Error: no junction file provided...\n")
        exit(0)
    # if options.noiseclass == None:
    #     sys.stderr.write("Error: no intron class annotation provided...\nRequired for current implementation...\n")
    #     exit(0)
    
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

    chromLst = [f"chr{x}" for x in range(1,23)]+['chrX','chrY'] + \
        [f"{x}" for x in range(1,23)]+['X','Y']

    main(options, libl)

    if not options.keeptemp:
        sys.stderr.write('Remove generated temp files... \n')
        with open(os.path.join(options.rundir, options.outprefix) + '_sortedlibs') as f:
            for tmp in [ln.strip() for ln in f.readlines()]:
                os.remove(tmp)
        os.remove(os.path.join(options.rundir, options.outprefix) + '_sortedlibs')
        sys.stderr.write('Done.\n')




