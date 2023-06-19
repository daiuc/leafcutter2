# leafcutter2 - (*alpha*)

## Prerequisites

- Minimum python version - `python v3.6`
- `regtools` - not required to run the script per se, but it is recommended to use regtools extracted junction files. For instance, we used this command to extract our junction files from bams `regtools junctions extract -a 8 -i 50 -I 500000 bamfile.bam -o outfile.junc` . See detailed regtools documentations [here](https://regtools.readthedocs.io/en/latest/commands/junctions-extract/).

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
    - junction files, processed using `regtools extract junctions`. 


Main output files:

- `{out_prefix}_perind.counts.noise.gz`: output functional introns (intact), and  noisy introns. Note the start and end coordinates of noisy introns are merged to the min(starts) and max(ends) of all functional introns within cluster.

- `{out_prefix}_perind_numers.counts.noise.gz`: same as above, except write numerators.

- `{out_prefix}_perind.counts.noise_by_intron.gz`: same as the first output, except here noisy introns' coordinates are kept as their original coordinates. This is useful for diagnosis purposes. 




## Run script

Recommended parameters for running the script:

```
python scripts/leafcutter2_regtools.py \
    -j junction_file_list_as_text_file \
    -N intron_junction_annotation_file \
    -o leafcutter2 \
    -r dir_to_house_leafcutter2_output
```

This mode first generate intron clusters based on the junction files. Then it counts junction reads towards each classified
introns.
-    `junction_file_list_as_text_file` should be a text file listing path to each junction file, one path per line
-    `intron_junction_annotation_file` is an annotation file for introns. A gencode based annotation file is provided at 
`data/gencode_v43_plus_v37_productive.intron_by_transcript_BEDlike.txt.gz`. Users can provide their own annotation file, 
but the format need to follow exactly, with last column being 'functional' or 'productive', and coordinates follow BED 
format. Users may also provide multiple annotation files, with file paths delimited by ",".
-    `-r` specify the directory of output, while `-o` specify the prefix of output file names (not including directory name)


In some cases, perhaps you want to first generate intron clusters, then run multiple rounds of noisy counts towards the 
same set of intron clusters, then you can first generate intron clusters. The script `scripts/leafcutter_make_clusters.py`
is a helpful script to make intron clusters separately. For instance, you may first run the clustering script to generate
intron clusters (filename: leafcutter_refined_noisy): 

```
python scripts/leafcutter_make_clusters.py \
    -j junction_file_list_as_text_file \
    -o leafcutter2 \
    -r dir_to_house_leafcutter2_output
```

Then run the leafcutter2 script to generate leafcutter2 outputs:

```
python scripts/leafcutter2_regtools.py \
    -j junction_file_list_as_text_file \
    -N intron_junction_annotation_file \
    -o leafcutter2 \
    -r dir_to_house_leafcutter2_output \
    -c leafcutter2_refined_noisy
```


### Parameters

```
python scripts/leafcutter2_regtools.py -h

usage: leafcutter2_regtools.py [-h] -j JUNCFILES 
                                                  [-o OUTPREFIX] [-q] [-r RUNDIR]
                                                  [-l MAXINTRONLEN] [-m MINCLUREADS]
                                                  [-M MINREADS] [-p MINCLURATIO]
                                                  [-c CLUSTER] [-k] [-C]
                                                  [-N NOISECLASS] [-f OFFSET] [-T]

optional arguments:
  -h, --help            show this help message and exit
  -j JUNCFILES, --juncfiles JUNCFILES
                        text file with all junction files to be processed
  -o OUTPREFIX, --outprefix OUTPREFIX
                        output prefix (default leafcutter)
  -q, --quiet           don't print status messages to stdout
  -r RUNDIR, --rundir RUNDIR
                        write to directory (default ./)
  -l MAXINTRONLEN, --maxintronlen MAXINTRONLEN
                        maximum intron length in bp (default 100,000bp)
  -m MINCLUREADS, --minclureads MINCLUREADS
                        minimum reads in a cluster (default 30 reads)
  -M MINREADS, --minreads MINREADS
                        minimum reads for a junction to be considered for
                        clustering(default 5 reads)
  -p MINCLURATIO, --mincluratio MINCLURATIO
                        minimum fraction of reads in a cluster that support a
                        junction (default 0.001)
  -c CLUSTER, --cluster CLUSTER
                        refined cluster file when clusters are already made
  -k, --checkchrom      check that the chromosomes are well formated e.g. chr1,
                        chr2, ..., or 1, 2, ...
  -C, --includeconst    also include constitutive introns
  -N NOISECLASS, --noise NOISECLASS
                        Use provided intron_class.txt.gz to help identify noisy
                        junction. Multiple annotation files delimited by `,`.
  -f OFFSET, --offset OFFSET
                        Offset sometimes useful for off by 1 annotations.
                        (default 0)
  -T, --keeptemp        keep temporary files. (default false)

```

## Examples

The `example` directory includes a snakemake noisy splicing QTL calling pipeline.

To run `snakemake`, change directory to `examples`. Run: 

```
snakemake -c1 results/noisy/leafcutter_perind.counts.noise.gz

```


