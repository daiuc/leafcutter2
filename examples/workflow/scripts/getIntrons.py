#!/usr/bin/env python

'''
Get intron bed file from provided transcript and exon bed files
'''

import os
import sys

__author__    = "Chao Dai"
__email__     = "chaodai@uchicago.edu"
__status__    = "Development"
__version__   =  "v0.0.1"


import pandas as pd
import pybedtools as pb


def getIntrons(transcripts, exons):
    '''Get introns from BED like dataframes of transcripts and exons
    
       transcripts : df, transcripts dataframe with the 6th column being
                     'geneID|transcriptID'
       exons       : df, same except for exons
    '''
    transcripts.columns = ['chrom', 'start', 'end', 'id', 'score', 'strand']
    exons.columns = ['chrom', 'start', 'end', 'id', 'score', 'strand']
    tei = {}
    for gid, tr in transcripts.groupby('id'):
        exon = exons[exons['id'] == gid]
        a = pb.BedTool.from_dataframe(tr)
        b = pb.BedTool.from_dataframe(exon)
        c = a.subtract(b, s=True)
        tei[gid] = {'transcript': a, 'exon': b, 'intron': c}
    return tei


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-T", "--TRANSCRIPT", dest="transcript", type=str,
        required=True, help="Transcript bed file")
    
    parser.add_argument("-E", "--EXON", dest="exon", type=str,
        required=True, help="exon bed file")
    
    parser.add_argument("-O", "--OUT", dest="out", type=str,
        required=True, help="output intron bed file")

    options = parser.parse_args()

    colnames = ['chrom', 'start', 'end', 'id', 'score', 'strand']

    sys.stdout.write(f'Load transcript annotation BED file {options.transcript}...\n')
    trpts = pd.read_csv(options.transcript, sep='\t', names=colnames)
    
    sys.stdout.write(f'Load exon annotation BED file {options.exon}...\n')
    exons = pd.read_csv(options.exon, sep='\t', names=colnames)

    sys.stdout.write(f'Subtract exons from transcripts to get introns...\n')
    tei = getIntrons(trpts, exons) # transcript, exons, introns

    introns = [x['intron'].to_dataframe() for x in tei.values()]
    introns = pd.concat(introns, axis=0).reset_index(drop=True)

    sys.stdout.write(f'Write introns to {options.out}...\n')
    introns.to_csv(options.out, sep='\t', header=False, index=False)
    sys.stdout.write('Done.\n')

