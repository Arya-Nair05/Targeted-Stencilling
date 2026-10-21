#!/bin/bash

## using the peaks called by FIBERTOOLS: converting to stanadard BED format
 gunzip  WT_WGS_FDR_peaks.bed | cut -f1-3 | sed '1d' > WT_WGS_peaks.bed
 gunzip  KO_WGS_FDR_peaks.bed | cut -f1-3 | sed '1d' > KO_WGS_peaks.bed

## Merging peaks using bedtools
 cat WT_WGS_peaks.bed KO_WGS_peaks.bed | sort -k1,1 -k2,2n > combined_peaks.bed
 bedtools merge -i combined_peaks.bed > KO_WT_merged_peaks.bed

## adding header line to the bed file
 awk 'BEGIN {print "chrom\tstart\tstop"} {print $0}' KO_WT_merged_peaks.bed > trial && mv trial KO_WT_merged_peaks.bed

## using SMGSTATS tool available at: https://github.com/exetergrdd/SMGStats.jl/tree/main 

## Finding differentially accesible regions using SMGSTATS diff function:
 WT_WGS=0307_H1_60min.fire.cram
 KO_WGS=051225_znf808.fire.cram
 KO_AS=101225_znf808_as_hac.haplotagged.bam

 smgstats diff --mods -p KO_WT_merged_peaks.bed \
  -a ${WT_WGS} \
  -b ${KO_WGS} ${KO_AS} \
  -o znf808_diff_access_CORRECTED.tsv.gz

## Haplotype specific DA regions 
 # ZNF808 KO AS:
 smgstats diff --mods --hap -p KO_WT_merged_peaks.bed \
  -a ${KO_AS} \
  -o znf808_AS_hap_DA.tsv.gz

 # ZNF808 KO WGS:
  smgstats diff -mods --hap -p KO_WT_merged_peaks.bed \
  -a ${KO_WGS} \
  -o znf808_WGS_hap_DA.tsv.gz

   
