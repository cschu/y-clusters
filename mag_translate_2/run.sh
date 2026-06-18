#!/bin/bash

# prepare MAG translation directory
# mkdir -p mag_translate_2 && cd mag_translate_2 && mkdir -p tmp

# get Y-contig id mapping
if [[ ! -f contig2ycontig.txt ]]; then
	awk -v OFS='\t' '{ printf("%s:%s\t%s\n", $1, $2, $3); }' /g/scb2/bork/ckim/HGT_MGE_project/HGT_MGE_project_re/data/SPIRE_contig_mapping.tsv | sort -T tmp -k1,1 > contig2ycontig.txt
fi


# sort SPIRE gene database dump by contig id
if [[ ! -f spire_genes.txt.sorted_by_contig ]]; then
	sort -k2,2 -T tmp /g/bork6/schudoma/projects/mge/spire_genes.txt > spire_genes.txt.sorted_by_contig
fi

# sort SPIRE bin database dump by bin id
if [[ ! -f spire_bins.txt.sorted_by_bin_id ]]; then
	sort -k1,1 /g/bork6/schudoma/projects/mge/spire_bins.txt > spire_bins.txt.sorted_by_bin_id
fi

# sort SPIRE contig database dump by SPIRE bin id
# bin_id, contig_id, kmer_size, ordinal
# 100008	1653229156	141	100713
if [[ ! -f spire_contigs.txt.sorted_by_bin ]]; then
	sort -T tmp -k1,1 /g/bork6/schudoma/projects/mge/spire_contigs.txt.2 | awk -v OFS='\t' '{ printf("%s\t%s\tk%s_%s\n", $1, $2, $3, $4); }' > spire_contigs.txt.sorted_by_bin
fi

# add SPIRE bin data to contig data
if [[ ! -f spire_contigs.txt.with_bins ]]; then
	join -1 1 -2 1 -a 1 spire_contigs.txt.sorted_by_bin spire_bins.txt.sorted_by_bin_id | tr " " "\t" > spire_contigs.txt.with_bins
fi

# generate joint-keys for adding in Y-MAG contig ids
if [[ ! -f spire_contigs.txt.with_bins.sorted_by_sample_contig ]]; then
	awk -v OFS='\t' '{ printf("%s:%s\t%s\n", gensub(/.psa_megahit.psb_metabat2.[0-9]{5}$/, "", "g", $4), $3, $0); }' spire_contigs.txt.with_bins | sort -T tmp -k1,1 > spire_contigs.txt.with_bins.sorted_by_sample_contig
fi

# add Y-MAG contig ids
if [[ ! -f spire_contigs.txt.with_bins.sorted_by_sample_contig.with_ycontig ]]; then
	join -1 1 -2 1 spire_contigs.txt.with_bins.sorted_by_sample_contig contig2ycontig.txt | tr " " "\t" | cut -f 2- > spire_contigs.txt.with_bins.sorted_by_sample_contig.with_ycontig
fi

# sort contigs by SPIRE bin name
if [[ ! -f spire_contigs.txt.with_bins.sorted_by_sample_contig.with_ycontig.sorted_by_bin_name ]]; then
	sort -T tmp -k4,4 spire_contigs.txt.with_bins.sorted_by_sample_contig.with_ycontig > spire_contigs.txt.with_bins.sorted_by_sample_contig.with_ycontig.sorted_by_bin_name
fi

# add specI to contigs
if [[ ! -f spire_contigs.txt.with_bins.sorted_by_sample_contig.with_ycontig.sorted_by_bin_name.with_speci ]]; then
	join -1 4 -2 1 spire_contigs.txt.with_bins.sorted_by_sample_contig.with_ycontig.sorted_by_bin_name ../../spire_speci_info.txt.gene_with_speci | tr " " "\t" > spire_contigs.txt.with_bins.sorted_by_sample_contig.with_ycontig.sorted_by_bin_name.with_speci
fi

# sort contigs by contig
if [[ ! -f spire_contigs.txt.with_bins.sorted_by_sample_contig.with_ycontig.sorted_by_bin_name.with_speci.sorted_by_contig ]]; then
	sort -T tmp -k3,3 spire_contigs.txt.with_bins.sorted_by_sample_contig.with_ycontig.sorted_by_bin_name.with_speci > spire_contigs.txt.with_bins.sorted_by_sample_contig.with_ycontig.sorted_by_bin_name.with_speci.sorted_by_contig
fi

# join gene data with annotated contigs
#join -1 2 -2 3 -a 1 spire_genes.txt.sorted_by_contig spire_contigs.txt.with_bins.sorted_by_sample_contig.with_ycontig.sorted_by_bin_name.with_speci.sorted_by_contig | tr " " "\t" > spire_genes_and_contigs.txt
if [[ ! -f spire_genes_and_contigs.txt ]]; then
	join -1 2 -2 3 -o 2.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.4,2.5,2.6,2.7 spire_genes.txt.sorted_by_contig spire_contigs.txt.with_bins.sorted_by_sample_contig.with_ycontig.sorted_by_bin_name.with_speci.sorted_by_contig | tr " " "\t" > spire_genes_and_contigs.txt
fi

# reducing size and ordering by ymag in preparation of addition of sp100/sp095 cluster 
if [[ ! -f spire_genes_and_contigs.txt.sorted_by_ymag_id ]]; then
	awk -v OFS='\t' '{ printf("%s\t%s_%s\t%s_%s\t%s\n", $1, $9, $3, $11, $3, $12); }' spire_genes_and_contigs.txt > spire_genes_and_contigs.txt.x && sort -T tmp -k3,3 spire_genes_and_contigs.txt.x > spire_genes_and_contigs.txt.sorted_by_ymag_id && rm spire_genes_and_contigs.txt.x
fi

# add Y-cluster information
if [[ ! -f spire_genes_and_contigs.txt.sorted_by_ymag_id.with_yclusters ]]; then
	join -1 3 -2 1 spire_genes_and_contigs.txt.sorted_by_ymag_id ../spire_genes_with_sp100_and_sp095.txt | tr " " "\t" > spire_genes_and_contigs.txt.sorted_by_ymag_id.with_yclusters
fi

# sort by speci
if [[ ! -f spire_genes_and_contigs.txt.sorted_by_ymag_id.with_yclusters.sorted_by_speci ]]; then
	sort -T tmp -k4,4 spire_genes_and_contigs.txt.sorted_by_ymag_id.with_yclusters > spire_genes_and_contigs.txt.sorted_by_ymag_id.with_yclusters.sorted_by_speci

	# break into speci clusters
	mkdir -p speci/ && awk -v OFS='\t' 'BEGIN { s=""; f=""; } { if (s!=$4) { if (s!="") { close(f); }; s=$4; f=sprintf("speci/specI_v4_%05d.txt", s); } print $0 >> f; } END { close(f); } ' spire_genes_and_contigs.txt.sorted_by_ymag_id.with_yclusters.sorted_by_speci

	# sort by sp095 cluster
	ls speci/specI_v4_?????.txt | xargs -I{} -P8 sh -c 'echo {}; sort -k6,6 -k1,1 {} > {}.by_sp095'
fi
