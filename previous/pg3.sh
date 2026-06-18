#!/bin/bash


join -1 1 -2 1 -a 1 -o 1.2,1.3,1.4,1.5,2.2 <(awk -v OFS='\t' '{print gensub(/_[0-9]+$/,"","g",$1),$0}' pg3_genes_with_sp100_and_sp095.txt) /scratch/schudoma/pg3_speci_info.txt.gene_with_speci | tr " " "\t" > pg3_genes_with_sp100_and_sp095.txt.with_speci


# cut -f 1,7- /g/scb2/bork/ckim/HGT_MGE_project_v3/1.combine_SPIRE/SP100/sp100_emapper_result.tsv | sed "s/ \+/_/g" | awk -v OFS='\t' '{ $1=gensub(/^SP100_0*([0-9]+)/, "\\1", "g", $1); print $0 }' | sed "s/'//g" | sort -k1,1 > sp100_emapper_for_joining.txt


sort -k5,5 -k1,1 -T tmp pg3_genes_with_sp100_and_sp095.txt.with_speci > pg3_genes_with_sp100_and_sp095.txt.with_speci.sorted_by_speci

mkdir -p speci
awk -v OFS='\t' 'BEGIN { s=""; f=""; } NF>4 { if (s!=$5) { if (s!="") { close(f); }; s=$5; f=sprintf("speci/specI_v4_%05d.txt", s); } print $0 >> f; } END { close(f); } ' pg3_genes_with_sp100_and_sp095.txt.with_speci.sorted_by_speci
