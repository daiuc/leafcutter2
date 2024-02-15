Identify noisy splicing

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

'