# SpliceJunctionClassified V0.1 (Updated Jan 2024)
# Written by Yang Li Nov-2023
# NOTE: this version works with 1-based coordinates. GTF junctions are kept as its
# original 1-based GTF format. BUT junctions from previous leafcutter perind.counts.gz
# files are converted from original 0-based BED format to non-standard format like so: 
# (start, end) -> (start, end+1).

import argparse
import gzip
import pickle
import sys

import pyfastx
from Bio.Seq import Seq
from tables.undoredo import attr_to_shadow


def check_utrs(junc, utrs):
    """
    checks if junction is close or within 100bp of UTRs
    """
    for s1, s2 in list(utrs):
        if abs(junc[0] - s1) < 100 or abs(junc[1] - s2) < 100:
            return True
    return False


def solve_NMD(
    chrom,
    strand,
    junc,
    start_codons,
    stop_codons,
    gene_name,
    fa,
    verbose=False,
    exonLcutoff=1000,
    annBase1=True, # 1-based annotation: start-codons, stop_codons
    juncBase1=False, # is junction 1-based? Default leafcutter2 is 0-based (like BED)
):
    """
    Compute whether there is a possible combination that uses the junction without
    inducing a PTC. We start with all annotated stop codon and go backwards.

    """

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

    junc_pass = {}
    junc_fail = {}
    path_pass = []
    proteins = []

    dic_terminus = {}

    depth = 0
    while len(seed) > 0:
        new_seed = []
        final_check = []
        depth += 1
        if verbose:
            sys.stdout.write("Depth %s, Seed L = %s\n" % (depth, len(seed)))
        # print(start_codons, [s[-1] for s in seed][-10:], len(junc))
        framepos = {}

        for s in seed:
            # first check that the seed paths are good
            bool_ptc = False
            leftover = ""
            if len(s) > 0:
                leftover = Seq("")
                allprot = Seq("")
                for i in range(0, len(s) - 1, 2):
                    exon_coord = s[i : i + 2]
                    exon_coord.sort()
                    exon_coord = tuple(exon_coord)
                    exlen = exon_coord[1] - exon_coord[0]

                    startpos = (len(leftover) + exlen + 1) % 3
                    if strand == "+":
                        seq = (
                            Seq(fa.fetch(chrom, (exon_coord[0], exon_coord[1])))
                            + leftover
                        )
                        prot = seq[startpos:].translate()
                        leftover = seq[:startpos]
                        allprot = prot + allprot
                    else:
                        seq = leftover + Seq(
                            fa.fetch(chrom, (exon_coord[0], exon_coord[1]))
                        )
                        aseq = seq
                        if startpos > 0:
                            leftover = seq[-startpos:]
                        else:
                            leftover = Seq("")
                        seq = seq.reverse_complement()
                        prot = seq[startpos:].translate()
                        allprot = prot + allprot
                    bool_ptc = "*" in allprot[:-1]

            if bool_ptc:
                # This transcript failed
                for i in range(1, len(s) - 1, 2):
                    j_coord = s[i : i + 2]
                    j_coord.sort()
                    j_coord = tuple(j_coord)
                    if j_coord not in junc_fail:
                        junc_fail[j_coord] = 0
                    junc_fail[j_coord] += 1

                continue

            # passed
            if len(s) > 2:
                terminus = (s[-2], s[-1], leftover)

                if terminus in dic_terminus:
                    dic_terminus[terminus].append(tuple(s))
                    continue
                else:
                    dic_terminus[terminus] = [tuple(s)]

            last_pos = s[-1]

            for start in start_codons:
                # print("start", start, abs(last_pos-start[0]))
                if (
                    strand == "+"
                    and last_pos > start[0]
                    and abs(last_pos - start[0]) < exonLcutoff
                ):
                    final_check.append(s + [start[0]])
                elif (
                    strand == "-"
                    and last_pos < start[1]
                    and abs(last_pos - start[1]) < exonLcutoff
                ):
                    final_check.append(s + [start[1]])
            for j0, j1 in junc:
                if strand == "+" and last_pos > j1 and abs(last_pos - j1) < exonLcutoff:
                    new_seed.append(s + [j1, j0])

                # print("junction", (j0,j1), abs(last_pos-j0))
                if strand == "-" and last_pos < j0 and abs(last_pos - j0) < exonLcutoff:
                    new_seed.append(s + [j0, j1])

        # check that the possible final paths are good
        for s in final_check:
            leftover = Seq("")
            allprot = Seq("")
            for i in range(0, len(s) - 1, 2):
                exon_coord = s[i : i + 2]
                exon_coord.sort()
                exon_coord = tuple(exon_coord)
                exlen = exon_coord[1] - exon_coord[0]
                startpos = (len(leftover) + exlen + 1) % 3
                if strand == "+":
                    seq = (
                        Seq(fa.fetch(chrom, (exon_coord[0], exon_coord[1]))) + leftover
                    )
                    leftover = seq[:startpos]
                    prot = seq[startpos:].translate()
                    allprot = prot + allprot
                else:
                    seq = leftover + Seq(
                        fa.fetch(chrom, (exon_coord[0], exon_coord[1]))
                    )
                    if startpos > 0:
                        leftover = seq[-startpos:]
                    else:
                        leftover = Seq("")
                    seq = seq.reverse_complement()
                    prot = seq[startpos:].translate()
                    allprot = prot + allprot
            bool_ptc = "*" in allprot[:-1]

            if not bool_ptc:
                # all pass
                proteins.append(
                    "\t".join(
                        [
                            gene_name,
                            chrom,
                            strand,
                            "-".join([str(x) for x in s]),
                            str(allprot),
                        ]
                    )
                    + "\n"
                )
                # print("ALL PASS %s"%(s))
                path_pass.append(tuple(s))
                for i in range(1, len(s), 2):
                    j_coord = s[i : i + 2]
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
                    if path[: len(path_subset)] == path_subset:
                        terminus_pass = True
                        break
            # print(terminus, terminus_pass)
            if terminus_pass:
                for path_subset in dic_terminus[terminus]:
                    if path_subset in path_pass:
                        continue
                    new_paths.append(path_subset)
                    path_pass.append(path_subset)
                    for i in range(1, len(path_subset), 2):
                        j_coord = list(path_subset[i : i + 2])
                        j_coord.sort()
                        j_coord = tuple(j_coord)
                        if j_coord not in junc_pass:
                            junc_pass[j_coord] = 0
                            if verbose:
                                sys.stdout.write("junction %s pass\n" % j_coord)
        if len(new_paths) == 0:
            break

    return junc_pass, junc_fail, proteins


def parse_gtf(gtf: str):
    """Lower level function to parse GTF file
    - gtf: str : path to GTF annotation file
    - returns: dictionary with keys:
        chrom, source, type, start, end, strand, frame, info, gene_name, transcript_type, transcript_name, gene_type
    """
    fields = ["chrom", "source", "type", "start", "end", ".", "strand", "frame", "info"]
    open_gtf = lambda x: gzip.open(x) if ".gz" in x else open(x)
    for ln in open_gtf(gtf):
        ln = ln.decode("ascii") if ".gz" in gtf else ln
        dic = {}
        if ln[0] == "#":
            continue
        ln = ln.strip().split("\t")
        for i in range(len(fields)):
            dic[fields[i]] = ln[i]

        # add 4 additional fields, parsed from info field
        for ks in ["gene_name", "transcript_type", "transcript_name", "gene_type"]:
            info_fields = [
                {x.split()[0]: x.split()[1].replace('"', "")}
                for x in dic["info"].split(";")
                if len(x.split()) > 1
            ]
            info_fields = {k: v for d in info_fields for k, v in d.items()}
            try:
                dic[ks] = info_fields[ks]
            except:
                dic[ks] = None  # if line is a gene, then wont have transcript info
        yield dic


def parse_annotation(gtf_annot: str):
    """
    Used `parse_gtf` to first parse a gtf file, then extract and return further
    information, including gene coordinates, intron info, and splice site to
    gene-name dictionary.

    - gtf_anno: str : path to GTF annotation file
    - returns:
        - genes_coords: gene coordinates, grouped by chromosome and strand
        - introns_info: a dictionary with these keys: junctions (ie all introns),
                        start_codon, stop_codon, utrs, pcjunctions (ie only protein coding introns)
        - ss2gene: a dictionary with splice site as keys, and gene name as values, eg. ('chr1', 11869): 'DDX11L1'
    """
    genes_info = {}
    introns_info = {}
    genes_coords = {}
    ss2gene = {}

    for dic in parse_gtf(gtf_annot):
        chrom = dic["chrom"]
        gname = dic["gene_name"]
        tname = dic["transcript_name"], gname
        anntype = dic["type"]
        if dic["type"] == "gene":
            dic["transcript_type"] = "gene"
        if (
            dic["transcript_type"] == "nonsense_mediated_decay"
            and anntype == "stop_codon"
        ):
            continue
        if dic["transcript_type"] != "protein_coding" and anntype == "UTR":
            continue
        start, end = int(dic["start"]), int(dic["end"])  # keep GTF 1-based
        strand = dic["strand"]

        if (chrom, strand) not in genes_coords:
            genes_coords[(chrom, strand)] = []

        if tname not in genes_info:  # tname is (transcript_name, gene_name)
            genes_info[tname] = {
                "exons": [],
                "start_codon": set(),
                "stop_codon": set(),
                "utrs": set(),
                "type": dic["transcript_type"],
            }

        if anntype in ["start_codon", "stop_codon"]:
            genes_info[tname][anntype].add((start, end))
        elif anntype in ["gene"]:
            genes_coords[(chrom, strand)].append(((start, end), gname))
        elif anntype in ["exon"]:
            genes_info[tname]["exons"].append((start, end))
            # Store gene info for splice sites
            ss2gene[(chrom, int(dic["start"]))] = dic["gene_name"]  # GTF
            ss2gene[(chrom, int(dic["end"]))] = dic["gene_name"]  # GTF

        elif anntype in ["UTR"]:
            genes_info[tname]["utrs"].add((start, end))

    for (
        gene
    ) in (
        genes_info
    ):  # this is actually transcript level!!! (transcript_name, gene_name)
        gene_name = gene[1]
        exons = genes_info[gene]["exons"]
        exons.sort()

        pc = False  # Basically use transcript_type == protein_coding to set flag
        if (
            genes_info[gene]["type"] == "protein_coding"
        ):  # AGAIN here gene = (transcript_name, gene_name)
            pc = True

        junctions = set()  # all junctions/introns
        pcjunctions = set()  # only protein coding junctions/introns

        for i in range(len(exons) - 1):
            intron = (exons[i][1], exons[i + 1][0])
            junctions.add(intron)
            if pc:
                pcjunctions.add(intron)

        if (
            gene_name not in introns_info
        ):  # here only gene_name, does not incl. transcript_name
            introns_info[gene_name] = {
                "junctions": set(),
                "start_codon": set(),
                "stop_codon": set(),
                "utrs": set(),
                "pcjunctions": set(),
            }

        introns_info[gene_name]["start_codon"] = introns_info[gene_name][
            "start_codon"
        ].union(genes_info[gene]["start_codon"])
        introns_info[gene_name]["stop_codon"] = introns_info[gene_name][
            "stop_codon"
        ].union(genes_info[gene]["stop_codon"])
        introns_info[gene_name]["junctions"] = introns_info[gene_name][
            "junctions"
        ].union(junctions)
        introns_info[gene_name]["utrs"] = introns_info[gene_name]["utrs"].union(
            genes_info[gene]["utrs"]
        )
        introns_info[gene_name]["pcjunctions"] = introns_info[gene_name][
            "pcjunctions"
        ].union(pcjunctions)

    return genes_coords, introns_info, ss2gene


def get_feature(fname, feature="exon"):
    ss2gene = {}
    if ".gz" in fname:
        F = gzip.open(fname)
    else:
        F = open(fname)
    for ln in F:
        if ln[0] == "#":
            continue
        ln = ln.split("\t")
        gID = ln[-1].split('gene_name "')[1].split('"')[0]

        if ln[2] != feature:
            continue
        ss2gene[(ln[0], int(ln[3]))] = gID
        ss2gene[(ln[0], int(ln[4]))] = gID
    return ss2gene


def get_overlap_stream(L1, L2, relax=0):
    """
    L1 and L2 are sorted lists of tuples with key, values
    """
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
            while L2[j + k][0][0] <= L1[i][0][1]:
                if overlaps(L1[i][0], L2[j + k][0]):
                    yield L1[i], L2[j + k]
                k += 1
                if j + k == len(L2):
                    break
            i += 1


def overlaps(A: tuple, B: tuple):
    """
    Checks if A and B overlaps
    A: tuple : (start, end)
    B: tuple : (start, end)
    """
    if A[1] < B[0] or B[1] < A[0]:
        return False
    else:
        return True


def ClassifySpliceJunction(
    perind_file: str, 
    gtf_annot: str, 
    fa: pyfastx.Fasta,
    rundir: str = ".",
    outprefix: str = "Leaf2",
    verbose: bool = False):
    """
    perind_file: str : LeafCutter perind counts file, e.g. leafcutter_perind.counts.gz
    gtf_annot: str : Annotation GTF file, for example gencode.v37.annotation.gtf.gz
    fa: pyfastx.Fasta : pyfastx.Fasta object for genome fasta file
    rundir: str : run directory, default is current directory
    outprefix: str : output prefix (default: Leaf2)
    verbose: bool : verbose mode (default: False)

    NOTE (Aug 2024): This script expects all coordinates to be 1-based. That means GTF annotations and processed intron junctions should be 1-based going into this the solve_NMDfunctionfunction.
    """

    # read leafcutter perind file and store junctions in dictionary: dic_junc
    # key = (chrom,strand), value = list of junctions [(start,end)]
    dic_junc = {}
    sys.stdout.write(f"Processing junction counts {perind_file}...")
    for ln in gzip.open(perind_file):
        junc_info = ln.decode("ascii").split()[0]  # first column
        if junc_info == "chrom":
            continue  # header

        chrom, start, end, clu_strand = junc_info.split(":")
        strand = clu_strand.split("_")[-1]
        if (chrom, strand) not in dic_junc:
            dic_junc[(chrom, strand)] = []
        dic_junc[(chrom, strand)].append((int(start), int(end)))

    #NOTE: convert start, end coordinates in dic_junc to (leafcutter1) 1-based from 0 based (like BED)
    for chrom, strand in dic_junc:
        dic_junc[(chrom, strand)] = [(x[0] , x[1] + 1) for x in dic_junc[(chrom, strand)]]

    sys.stdout.write("done!\n")
    if verbose:
        sys.stdout.write("Processed: ")
        for chrstrand in dic_junc:
            sys.stdout.write(
                f"{len(dic_junc[chrstrand])} jxns on {chrstrand[0]} ({chrstrand[1]})."
            )

    # load or parse gtf annotations
    # g_coords: gene coordinates, grouped by chromosome and strand
    # g_info: a dictionary with (transcript_name, gene_name) as keys, and intron info as values

    try:
        sys.stdout.write("Loading annotations...\n")
        parsed_gtf = f"{rundir}/{gtf_annot.split('/')[-1].split('.gtf')[0]}_SJC_annotations.pckle"
        with open(parsed_gtf, "rb") as f:
            g_coords, g_info = pickle.load(f)
        sys.stdout.write("done!\n")
    except:
        sys.stdout.write("Parsing annotations for the first time...\n")
        g_coords, g_info, ss2gene = parse_annotation(gtf_annot)

        for chrom, strand in g_coords:
            to_remove_gcoords = set()
            to_remove_ginfo = set()
            for gene in g_coords[(chrom, strand)]:
                if gene[1] in g_info:  # check if gene_name is in introns_info keys
                    if (
                        len(g_info[gene[1]]["stop_codon"]) == 0
                    ):  # if there are no stop codons
                        to_remove_gcoords.add(gene)  # mark gene for remove
                        to_remove_ginfo.add(gene[1])  # mark gene_name for removal

            for g in to_remove_gcoords:
                g_coords[(chrom, strand)].remove(g)
            for g in to_remove_ginfo:
                g_info.pop(g)
        sys.stdout.write("Saving parsed annotations...\n")
        with open(parsed_gtf, "wb") as f:
            pickle.dump((g_coords, g_info), f)

    gene_juncs = {} # store junctions that overlap with a gene
    for chrom, strand in dic_junc:
        if (chrom, strand) not in g_coords:
            sys.stderr.write(f"Could not find {chrom} ({strand}) in annotations...\n")
            continue
        juncs = [(x, x) for x in dic_junc[(chrom, strand)]]
        juncs.sort()

        coords = g_coords[(chrom, strand)]
        coords.sort()

        # save junctions that overlapping a gene in gene_juncs dictionary: (gene_name, chrom, strand) : [junctions]
        for junc, geneinfo in get_overlap_stream(juncs, coords):
            info = (geneinfo[1], chrom, strand)
            if info not in gene_juncs:
                gene_juncs[info] = []
            gene_juncs[info].append(junc[0])

    sys.stdout.write("Classifying junctions...\n")
    fout = open(f"{rundir}/{outprefix}_junction_classifications.txt", "w")
    fout.write(
        "\t".join(["Gene_name", "Intron_coord", "Annot", "Coding", "UTR"]) + "\n"
    )

    for gene_name, chrom, strand in gene_juncs:
        sys.stdout.write(f"Processing {gene_name} ({chrom}:{strand})\n")

        query_juncs = gene_juncs[
            (gene_name, chrom, strand)
        ]  # junctions from perind files that overlap with gene, 1 based here
        if gene_name not in g_info:
            continue
        junctions = g_info[gene_name]["junctions"]  # juncs from gtf

        # classify a union set of juncs from gtf and from perind file (overlapping with gene)
        junctions = list(junctions.union(query_juncs))

        start_codons = g_info[gene_name]["start_codon"]
        stop_codons = g_info[gene_name]["stop_codon"]

        if verbose:
            sys.stdout.write(
                f"LeafCutter junctions ({len(query_juncs)}) All junctions ({len(junctions)}) Start codons ({len(start_codons)}) Stop codons ({len(stop_codons)}) \n"
            )

        junc_pass, junc_fail, proteins = solve_NMD(
            chrom, strand, junctions, start_codons, stop_codons, gene_name, fa
        )

        for j in junctions:
            bool_pass = j in junc_pass or j in g_info[gene_name]["pcjunctions"]
            bool_fail = j in junc_fail
            utr = False
            if not bool_pass:
                # Check that it's not in UTR
                utr = check_utrs(j, g_info[gene_name]["utrs"])

            if bool_fail or bool_pass:
                tested = True
            else:
                tested = False
            annotated = j in g_info[gene_name]["junctions"]
            # if not bool_pass and annotated:
            # print("%s %s %s junction: %s tested: %s utr: %s coding: %s annotated: %s "%(chrom, strand, gene_name, j, tested,utr, bool_pass, annotated))

            #NOTE: revert to bed format (0-based) for output
            fout.write(
                "\t".join(
                    [
                        gene_name,
                        f"{chrom}:{j[0]}-{j[1]-1}", # revert back to 0-based
                        str(annotated),
                        str(bool_pass),
                        str(utr),
                    ]
                )
                + "\n"
            )
    fout.close()
    sys.stdout.write(f"done!Classified junctions written to {rundir}/{outprefix}_junction_classifications.txt\nNote coordinates are 0-based.")



def main(options):

    if options.countfile is None:
        sys.stderr.write(
            "Error: no LeafCutter junction file provided...\npython SpliceJunctionClassifier.py -j Leafcutter_perind.counts.gz\n"
        )
        exit(0)

    if options.genome is None:
        sys.stderr.write("Error: no genome fasta file selected...\npython SpliceJunctionClassifier.py -G genome.fa\n")
        exit(0)

    if options.annot is None:
        sys.stderr.write("Error: no annotation file with gene start and stop codon...\npython SpliceJunctionClassifier.py -A gencode.gtf/ensembl.gtf\n")
        exit(0)
 
    sys.stdout.write(f"Loading genome {options.genome} ...")
    fa = pyfastx.Fasta(options.genome)
    sys.stdout.write("done!\n")

    ClassifySpliceJunction(
        perind_file=options.countfile,
        gtf_annot=options.annot,
        fa=fa,
        rundir=options.rundir,
        outprefix=options.outprefix,
        verbose=options.verbose,
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="SpliceJunctionClassifier")

    parser.add_argument(
        "-c", "--countfile", dest="countfile", help="LeafCutter perind counts file, e.g. leafcutter_perind.counts.gz",
    )

    parser.add_argument(
        "-o", "--outprefix", dest="outprefix", default="Leaf2", help="output prefix (default: Leaf2)",
    )

    parser.add_argument(
        "-r", "--rundir", dest="rundir", default=".", help="run directory (default: .)"
    )

    parser.add_argument(
        "-A", "--annotation", dest="annot", help="Annotation GTF file",
    )

    parser.add_argument(
        "-G", "--genome", dest="genome", help="Reference genome fasta file."
    )

    parser.add_argument(
        "-v", "--verbose", dest="verbose", action="store_true", default=False, help="verbose mode",
    )

    options = parser.parse_args()

    main(options)
