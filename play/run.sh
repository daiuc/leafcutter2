#!/bin/zsh

# run leafcutter

RUN1=run1

# python leafcutter2_regtools.py \
#   -j <(ls -1 NA*.junc) \
#   -o $RUN1 \
#   -r $RUN1 \
#   -k \
#   -A ./gencode.v43.pri.protein.gtf.gz \
#   -G chr1.fa.gz

RUN2=run2

# python leafcutter2_regtools_v2.py \
#   -j <(ls -1 NA*.junc) \
#   -o $RUN2 \
#   -r $RUN2 \
#   -k \
#   -A ./gencode.v43.pri.protein.gtf.gz \
#   -G chr1.fa.gz

RUN2=run3
# leafcutter2 now fixed with Quinn's script
python leafcutter2_regtools_v2.py \
  -j <(ls -1 NA*.junc) \
  -o $RUN2 \
  -r $RUN2 \
  -k \
  -A ./gencode.v43.pri.protein.gtf.gz \
  -G chr1.fa.gz

# run spliceclass v2 on bed formatte data
# python ./SpliceJunctionClassifier_v2.py -c ./qn/toyBED_perind.counts.gz -o v2onBed -r qn -A ./gencode.v43.pri.protein.gtf.gz -G ./chr1.fa.gz

# run forward on bed + 1 on end data
# python ./ForwardSpliceJunctionClassifier.py -c ./qn/toy_perind.counts.gz -o FsjcOnLf1 -r qn -A ./gencode.v43.pri.protein.gtf.gz -G ./chr1.fa.gz

# run forward v2 on bed format
# python ./ForwardSpliceJunctionClassifier_BED.py -c ./qn/toyBED_perind.counts.gz -o Fsjc2OnBed -r qn -A ./gencode.v43.pri.protein.gtf.gz -G ./chr1.fa.gz
