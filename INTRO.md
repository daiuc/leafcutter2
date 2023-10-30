## NOTE

- Make sure intron junction annotation files are BED formatted (0 based left close right open). This is different from the standard leafcutter.
- leafcutter2 outputs BED formatted coordinates.
- key differences from leafcutter:
    - leafcutter2 use the same set of filters to construct intron clusters.
    - However, when counting junction reads towards predefined or on-demand-run intron clusters, no read filter is applied, essentially all junction reads are counted towards introns.
    - The filtering options (--MINCLUREADS, --MINREADS, --MINCLURATIO) are only used if you are not providing a pre-defined set of intron clusters. They are not used against counting junction reads towards introns.

## Introduciton

This script takes junction files produced by regtools as input and constructs intron clusters from them. Then, it processes the intron clusters to identify rarely spliced introns based on certain filtering cut-offs. The output is a text file that follows the same format as the standard leafcutter tool. The first column of each row shows the genome coordinates of introns and labels them as N[: noisy] or F[: functional/productive].


Main input files:
    - junction files, processed using `regtools extract junctions` from bam files. Details can be found [here](https://regtools.readthedocs.io/en/latest/commands/junctions-extract/).


Main output files:

- `{out_prefix}_perind.counts.noise.gz`: output functional introns (intact), and noisy introns. Note the start and end coordinates of noisy introns are merged to the min(starts) and max(ends) of all functional introns within cluster.

- `{out_prefix}_perind_numers.counts.noise.gz`: same as above, except write numerators.

- `{out_prefix}_perind.counts.noise_by_intron.gz`: same as the first output, except here noisy introns' coordinates are kept as their original coordinates. This is useful for diagnosis purposes. 

- Since 2023-10-15 an update was made to the script such that any outputs in the `*counts.noise.gz` and `*noise_by_intron.gz` only contains introns that have passed a standard deviation threshold (default 0.5). In other words, each intron must have std >= 0.5 across samples in order to qualify for an output.


## strategy

![splicing demo](leafcutter-demo.png)

## Functions explained

### Pool junction reads `pool_junc_reads`

This step is to pool junction reads from multiple junction files. The output is a text file that stores introns and their junction reads. Each line stores a cluster of introns,
which can be linked or unlinked. The format is like this:

```
chr1:+ 779093:803918:3 793042:795469:4
```

### Refine clusters `refine_clusters`

This step is to refine clusters from the pooled junction reads step. The output is a text file that stores refined clusters. Refined means that each intron within a cluster must
be linked (sharing a splice site). And both intron clusters and introns subject to some
minimum read count thresholds. The format is like this:
```
chr1:- 5754:5878:5
chr1:- 5810:6470:405 5890:6470:9 6173:6470:150
```

### Add back introns with low usage `addlowusage`

Because we are interested in studying unproductive introns, which often has very
low read counts. We want to add these lowly used introns back to the refined intron 
list. The output `{out_prefix}_refined_noisy` is a list of introns, organized within
clusters that we will use in subsequent steps for annotating introns / junctions that 
we want to analyze.


### Sort junctions `sort_junctions`

Note the purpose of the previous steps is to create a list of introns to be used as
a source of annotation. We could skip the previous steps and provide a list of introns
with appropriate format for use from this step onward. The junctions here are the
input junctions we will be analyzing for our analyses. We first sort our input 
junctions using the intron list we created in the previous steps or provided by the
user (with `-c` flag). The output is a sorted junction file for each library.

### Merge junctions `merge_junctions`

Each junction file represents one sample. Here we combine all samples into one file.
Rows are introns, columns are counts for each library/sample.

### Annotate junctions `annotate_junctions`

This step annotate our input junctions with the intron list (either from first 3 steps
 or provided by user). The goal is to annotate each junction with a cluster ID, and 
 denote whether this intron is productive or unproductive. At the moment, we use a 
 curated GENCODE annotation to determine whether an intron is productive or unprouctive.
 Each intron is subject a minimum standard deviation threshold (default 0.5) across
 samples. This step output two files `{out_prefix}_perind.counts.noise.gz` and `{out_prefix}_perind.counts.noise_by_intron.gz`. The `*counts.noise.gz` file contains
 productive introns, which uses its original coordinates, and unproductive introns,
 which uses the min(starts) and max(ends) of all productive introns within the cluster.
 So the unproductive counts here represent an aggretation of all unproductive introns 
 within the cluster. The `*counts.noise_by_intron.gz` file contains the same information
as the `*counts.noise.gz` file, except the unproductive introns are not aggregated. And
their coordinates are kept as their original coordinates.

