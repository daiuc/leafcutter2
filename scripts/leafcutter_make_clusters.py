#!/usr/bin/env python

'''Cluster introns(junctions)

This script take junction files and output intron clusters. By default, only 
clusters that have linked introns (sharing 5' or 3' splice site) are written
out. If --includeconst is on, then both linked intron clusters and constitutive
clusters are written out. Regardless, all introns within a cluster with linked
introns are written out, including both linked introns and unlinked introns.

By default, no read filters are applied.


Junction files are processed using regtools:

    * https://github.com/griffithlab/regtools
    * /home/yangili1/tools/regtools/build/regtools junctions extract -a 8 -i \
      50 -I 500000 bamfile.bam -o outfile.junc
    * Using regtools speeds up the junction extraction step by an order of 
      magnitude or more

Procedure sequence: 

    1. pool_junc_reads:
        Pool introns from a list of junction files. Output file is a text file,
        where each line stores introns that belong to a single cluster, eg:
            chr1:- 5754:5878:5 5810:6470:405 5890:6470:9 6173:6470:150
        First column is `chrom:strand`; subsequent columns are `start:end:reads ...`
        Any junctions that overlap are considered within an intron cluster.
    
    2. refine_clusters:
        Take the pooled junction reads from the 1st procedure and refine it,
        such that filters such as minimal junction reads, minimal cluster
        reads, or ratios are met. Importantly, it ensures introns within 
        a cluster must be linked - sharing either a 5' or a 3' splice site.
        Output file is a text file, where each line stores linked intron cluster:
            chr1:- 5754:5878:5
            chr1:- 5810:6470:405 5890:6470:9 6173:6470:150

NOTE: 
    * Minimum version requirement - python v3.6

'''

import sys
import tempfile
import os
import gzip
import shutil


__author__    = "Yang Li, Chao Dai"
__email__     = "chaodai@uchicago.edu"
__status__    = "Development"
__version__   =  "v0.0.1"


def natural_sort(l): 
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


def cluster_intervals(E):
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


def overlaps(A,B):
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
            # read_ks = [k for k,v in by_chrom[chrom].items()] # read-keys, do not enforce min junction reads
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
    lastIntron = False
    while len(unassigned) > 0:
        finished = False
    
        while not finished:
            torm = [] # list fo introns to remove from unassigned
            finished = True
            for intron in unassigned:
                (start, end), count = intron
                if start in splicesites or end in splicesites: # if overlapping 5' or 3'
                    current.append(intron) # store intron 
                    splicesites.add(start) # store 5' and 3' to splicesites
                    splicesites.add(end)
                    finished = False
                    torm.append(intron) # mark to remove after check linkage
            
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


def refine_cluster(clu, cutoff, readcutoff):
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
        * each cluster's introns are linked (share either 5' or 3' splice site)
        * minimum total cluster reads cutoff
        * minimum intron reads cutoff
        * minimum reads ratio per cluster cutoff
    
    However, if constitutive flag `const` is on, then non-linked singleton introns
    are also written out, and they do not subject to cluster ratio and cluster 
    reads cutoff filters.

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
    if minclureads == 0: 
        sys.stderr.write(f"\n minimum cluster reads set to 0, skip any filtering on reads or ratios...\n")

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

        #[[((413430, 423479), 3)], [((410646, 413144), 3), ((410646, 413147), 62), ((410646, 420671), 4)]] 
        for cl in refine_linked(clu): # only linked intron clusters
            if minclureads != 0:
                rc = refine_cluster(cl, minratio, minreads)
            else: rc = [cl]
            
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

    refined_cluster = f"{rundir}/{outPrefix}_refined"

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

#-------------------------------------------

def main(options, libl):
    
    pool_junc_reads(libl, options)
    refine_clusters(options)
    addlowusage(options)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="Cluster introns into clusters of linked introns. Recommend \
            using default parameters, which do not enforce any reads/ratio filters."
    )
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
                        clustering. (default 5 reads)")

    parser.add_argument("-p", "--mincluratio", dest="mincluratio", 
        default = 0.001,
        help="minimum fraction of reads in a cluster that support a junction \
              (default 0.001)")

    parser.add_argument("-k", "--checkchrom", dest="checkchrom",
        action="store_false", default = True,
        help="check that the chromosomes are well formated e.g. chr1, chr2, \
              ..., or 1, 2, ... (default: true)")
    
    parser.add_argument("-C", "--includeconst", dest="const", \
        action="store_true", default = False, 
        help="also include constitutive introns. (default: false)")

    parser.add_argument("-f", "--offset", dest="offset", default = 0,
        help="Offset sometimes useful for off by 1 annotations. (default 0)")

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

    chromLst = [f"chr{x}" for x in range(1,23)]+['chrX','chrY'] + \
        [f"{x}" for x in range(1,23)]+['X','Y']

    main(options, libl)
