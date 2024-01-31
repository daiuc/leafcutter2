# SpliceJunctionClassified V0.1 (Updated Jan 2024)
# Written by Yang Li Nov-2023


def main(perind_file,gtf_annot = "gencode.v37.annotation.gtf", verbose = False):

    # Load junctions from leafcutter output


    dic_junc = {}
    sys.stdout.write("Processing junctions...")
    for ln in gzip.open(perind_file):
        junc_info = ln.decode('ascii').split()[0]
        if junc_info == "chrom": continue
        
        chrom, start, end, clu_strand = junc_info.split(":")
        strand = clu_strand.split("_")[-1]
        if (chrom,strand) not in dic_junc: 
            dic_junc[(chrom,strand)] = []
        dic_junc[(chrom,strand)].append((int(start), int(end)))

    sys.stdout.write("done!\n")
    if verbose:
        sys.stdout.write("Processed: ")
        for chrstrand in dic_junc:
            sys.stdout.write("%s jxns on %s (%s). "%(len(dic_junc[chrstrand]), chrstrand[0], chrstrand[1]))

    

    try: 
        sys.stdout.write("Loading annotations...")
        g_coords, g_info = pickle.load(open(gtf_annot.replace(".gtf","")+"_SJC_annotations.pckle",'rb'))
        sys.stdout.write("done!\n")
    except:
        sys.stdout.write("Parsing annotations for the first time...\n")
        g_coords, g_info, ss2gene = parse_annotation(gtf_annot)
        
        for chrom,strand in g_coords:
            to_remove_gcoords = set()
            to_remove_ginfo = set()
            for gene in g_coords[(chrom,strand)]:
                if gene[1] in g_info:
                    if len(g_info[gene[1]]['stop_codon']) == 0:
                        to_remove_gcoords.add(gene)
                        to_remove_ginfo.add(gene[1])

            for g in to_remove_gcoords:
                g_coords[(chrom,strand)].remove(g)
            for g in to_remove_ginfo:
                g_info.pop(g)

        pickle.dump((g_coords,g_info), open(gtf_annot.replace(".gtf","")+"_SJC_annotations.pckle",'wb'))

    gene_juncs = {}
    for chrom,strand in dic_junc:
        if (chrom,strand) not in g_coords: 
            sys.stderr.write("Could not find %s (%s) in annotations...\n"%(chrom,strand))
            continue
        juncs = [(x,x) for x in dic_junc[(chrom,strand)]]
        juncs.sort()

        coords = g_coords[(chrom,strand)]
        coords.sort()
        
        for junc, geneinfo in get_overlap_stream(juncs,coords): 
            info = (geneinfo[1], chrom,strand)
            if info not in gene_juncs:
                gene_juncs[info] = []
            gene_juncs[info].append(junc[0])

    fout = open(options.outprefix+"_junction_classifications.txt",'w')
    #fout_prot = open("FLP_annotations.txt",'w') 
    fout.write("\t".join(["Gene_name","Intron_coord","Annot","Coding", "UTR"])+'\n')
    
    for gene_name, chrom, strand in gene_juncs:
        sys.stdout.write("Processing %s (%s:%s)\n"%(gene_name,chrom,strand))
        
        querry_juncs = gene_juncs[(gene_name,chrom,strand)]
        if gene_name not in g_info: continue
        junctions = g_info[gene_name]['junctions']
        
        # classify all junctions in gene
        junctions = list(junctions.union(querry_juncs))

        start_codons = g_info[gene_name]['start_codon'] 
        stop_codons = g_info[gene_name]['stop_codon']

        if verbose:
            sys.stdout.write("LeafCutter junctions (%s) All junctions (%s) Start codons (%s) Stop codons (%s) \n"%(len(querry_juncs),  len(junctions), len(start_codons), len(stop_codons)))

        junc_pass, junc_fail, proteins = solve_NMD(chrom,strand,junctions, start_codons, stop_codons, gene_name)        
        
        #for prot in proteins:
        #    fout_prot.write(prot)

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
            
            fout.write('\t'.join([gene_name, '%s:%s-%s'%(chrom,j[0],j[1]), str(annotated), str(bool_pass), str(utr)])+'\n')
        
            

def check_utrs(junc,utrs):
    '''
    checks if junction is close or within UTRs
    '''
    for s1,s2 in list(utrs):
        if abs(junc[0]-s1) < 100 or abs(junc[1]-s2) < 100:
            return True
    return False

def solve_NMD(chrom, strand, junc, start_codons, stop_codons,gene_name, verbose = False, exonLcutoff = 1000):
    '''
    Compute whether there is a possible combination that uses the junction without
    inducing a PTC. We start with all annotated stop codon and go backwards.
    '''

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

                #print("Fail",s)
                #print(allprot)
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
                    
                    #print(dic_juncbookmark)

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

def parse_gtf(gtf):
    fields = ["chrom", "source", "type", "start", "end", ".","strand","frame", "info"]
    for ln in open(gtf):
        dic = {}
        if ln[0] == "#": continue
        ln = ln.strip().split('\t')
        for i in range(len(fields)):
            dic[fields[i]] = ln[i]

        for ks in ['gene_name', "transcript_type","transcript_name", "gene_type"]:
            dic[ks] =  dic['info'].split(ks)[-1].split('"')[1].split('"')[0]
        yield dic
         

def parse_annotation(gtf_annot):
    ''' 
    Parse annotation by gene, and record junctions and Start/Stop codons.
    returns genes_coords with genes coordinates by chromosome and strand
    as well as info on start/stop locations
    '''
    genes_info = {}
    introns_info = {}
    genes_coords = {}
    ss2gene = {}

    for dic in parse_gtf(gtf_annot):
        chrom = dic['chrom']
        gname = dic['gene_name']
        tname = (dic['transcript_name'],gname)
        anntype  = dic['type']
        if dic['type'] == 'gene': dic['transcript_type'] = "gene"
        if dic['transcript_type'] == "nonsense_mediated_decay" and anntype == "stop_codon": continue 
        if dic['transcript_type'] != "protein_coding" and anntype == "UTR": continue
        start, end = int(dic['start']), int(dic['end'])
        strand = dic['strand']
    
        if (chrom,strand) not in genes_coords:
            genes_coords[(chrom,strand)] = []

        if tname not in genes_info:
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
            ss2gene[(chrom, int(dic['start']))] = dic['gene_name']
            ss2gene[(chrom, int(dic['end']))] = dic['gene_name']

        elif anntype in ['UTR']:
            genes_info[tname]['utrs'].add((start,end))

    for gene in genes_info:
        gene_name = gene[1]
        exons = genes_info[gene]['exons']
        exons.sort()

        pc = False
        if genes_info[gene]['type'] == "protein_coding":
            pc = True

        junctions = set()
        pcjunctions = set()

        for i in range(len(exons)-1):
            intron = (exons[i][1], exons[i+1][0])
            junctions.add(intron)
            if pc:
                pcjunctions.add(intron)

        if gene_name not in introns_info:
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

    return (genes_coords,introns_info,ss2gene)



def getmean(lst):
    return sum(lst)/float(len(lst))

def getmedian(lst):
    lst.sort()
    n = len(lst)
    if n < 1:
        return None
    if n % 2 == 1:
        return lst[n//2]
    else:
        return sum(lst[n//2-1:n//2+1])/2.0

def get_feature(fname, feature = "exon"):
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


def get_overlap_stream(L1, L2, relax = 0):
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

def overlaps(A,B):
    '''                                                                                                                                                                
    Checks if A and B overlaps                                                                                                                                         
    '''
    if A[1] < B[0] or B[1] < A[0]:
        return False
    else: return True

if __name__ == "__main__":


    from optparse import OptionParser
    parser = OptionParser()

    import pickle
    import sys
    import pyfastx
    import gzip
    import pickle
    import numpy as np
    from Bio.Seq import Seq
    global fa

    parser = OptionParser()

    parser.add_option("-j", "--juncfile", dest="juncfile",
                  help="LeafCutter perind file with all junction files to be processed")

    parser.add_option("-o", "--outprefix", dest="outprefix", default = 'Leaf2',
                      help="output prefix (default: Leaf2)")

    parser.add_option("-a", "--annotation", dest="annot", default = 'gencode.v37.annotation.gtf',
                  help="Annotation file (default: gencode.v37.annotation.gtf)")
    
    (options, args) = parser.parse_args()
    
    if options.juncfile is None:
        sys.stderr.write("Error: no LeafCutter junction file provided...\npython SpliceJunctionClassifier.py -j Leafcutter_perind.counts.gz\n")
        exit(0)

    genome_fasta = 'genome.fa'
    sys.stdout.write("Loading genome...")
    fa = pyfastx.Fasta(genome_fasta)
    sys.stdout.write("done!\n")
    main(options.juncfile, options.annot)

