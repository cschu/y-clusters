#!/bin/bash


## PREPARE SP100

# shrink sp100 cluster ids, remove cluster size column (2)
awk -v OFS='\t' '{ print gensub(/^SP100_0*([0-9]+)/, "\\1", "g", $1),$3; }' /g/scb2/bork/ckim/HGT_MGE_project_v3/1.combine_SPIRE/SP100/SP100_members.tsv > SP100_members.short.tsv
chmod a-w SP100_members.short.tsv

# remove PG3 (=isolate) genes and (now) empty PG3 singleton clusters 
# GCA_000006665_04263
sed "s/GCA_[0-9]\+_[0-9]\+\+;\?//g" SP100_members.short.tsv | awk -v OFS='\t' 'NF>1 { print $0 }' > SP100_members.short.no_isolates_no_isolate_singletons.tsv
chmod a-w SP100_members.short.no_isolates_no_isolate_singletons.tsv

# separate MAG singleton and non-singleton clusters (not sure if necessary, but seems reasonable to assume to facilitate processing)
awk -F ';' 'NF==1 { print $0 >> "SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv"; next; } { print $0 >> "SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv" } END { close("SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv"); close("SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv") }' SP100_members.short.no_isolates_no_isolate_singletons.tsv

# sort sp100 singletons
sort -k1,1 SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv > SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted

# sort sp100 multi clusters
sort -k1,1 SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv > SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv.sorted

## PREPARE SP095

# shrink sp095/sp100 cluster ids, remove size columns (2,3)
awk -v OFS='\t' '{ print gensub(/^SP095_0*([0-9]+)/, "\\1", "g", $1),gensub(/SP100_0*([0-9]+)/, "\\1", "g", $4); }' /g/scb2/bork/ckim/HGT_MGE_project_v3/2.linclust/SP095/SP095_members.tsv > SP095_members.short.tsv

# linearize sp095 -> sp100\tsp095
awk -v OFS='\t' '{ split($2,a,";"); for (i in a) { print a[i],$1; }; }' SP095_members.short.tsv > SP095_members.short.tsv.linear

# sort sp095
sort -k1,1 SP095_members.short.tsv.linear > SP095_members.short.tsv.linear.sorted


## JOIN

# join sp095 -> sp100 singletons 
join -1 1 -2 1 -o 1.1,2.2,1.2 SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted SP095_members.short.tsv.linear.sorted | tr " " "\t" > SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095

# join sp095 -> sp100 multi clusters
join -1 1 -2 1 -o 1.1,2.2,1.2 SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv.sorted SP095_members.short.tsv.linear.sorted | tr " " "\t" > SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv.sorted.with_sp095

# add gene cluster type and unravel non-singletons
awk -F '\t' -v OFS='\t' '{n=split($3,a,";"); ctype="N"; for (i in a) { if (a[i] != "") print a[i],$1,$2,ctype;}}' SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv.sorted.with_sp095 > SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv.sorted.with_sp095.by_gene

# sort non-singletons
sort -k1,1 -k2,2 -k3,3 SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv.sorted.with_sp095.by_gene > SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv.sorted.with_sp095.by_gene.sorted

# get all SP100 clusters that contain isolate genes
cut -f 1 SP100_members.short.no_mags_no_mag_singletons.isolate_clusters.tsv.sorted > SP100_clusters_with_isolates.mc.tsv
cut -f 1 SP100_members.short.no_mags_no_mag_singletons.isolate_singletons.tsv.sorted > SP100_clusters_with_isolates.sc.tsv
sort -k1,1 -m SP100_clusters_with_isolates.mc.tsv SP100_clusters_with_isolates.sc.tsv > SP100_clusters_with_isolates.tsv
# add "S"ingleton type column
awk -v OFS='\t' '{print$0,"S"}' SP100_clusters_with_isolates.tsv > SP100_clusters_with_isolates.tsv.with_ctype

# add ctype column to singletons 
join -1 1 -2 1 -a 2 -o 2.1,2.2,2.3,1.2 SP100_clusters_with_isolates.tsv.with_ctype SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095 | tr " " "\t" | awk -v OFS='\t' 'NF==3 { $4="U" } { print $0 }' > SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095.with_ctype

# unravel singletons
awk -v OFS='\t' '{ split($3,a,";"); for (i in a) { print a[i],$1,$2,$4; }; }' SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095.with_ctype > SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095.with_ctype.by_gene

# sort singletons
sort -k1,1 SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095.with_ctype.by_gene > SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095.with_ctype.by_gene.sorted

# merge
sort -k1,1 -m SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv.sorted.with_sp095.by_gene.sorted SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095.with_ctype.by_gene.sorted > spire_genes_with_sp100_and_sp095.txt 



# GET MAG DATA

# shrink metadata table
# potentially requires editing if study names have changed
cut -f 1-4,10 /g/scb2/bork/ckim/HGT_MGE_project_v3/data/SPIRE_genome_annotations.tsv > spire_mags.txt
cp spire_mags.txt spire_mags.txt.original  # then modify spire_mags.txt

# create translation table
cut -f 1,2 spire_mags.txt | sort -k2,2 > spire_mags.txt.y2bin

# extract all gene lists
mkdir mag_translate && cd mag_translate
awk -v OFS='\t' '{ a[$4"/"$3"_per_bin_genecalls.tar.gz"]=1 } END {for (i in a) print i}' spire_mags.txt | xargs -I{} -P8 sh -c $'sample=$(basename {} _per_bin_genecalls.tar.gz); echo $sample; tar xzf {} \'*.fna\'; grep -H "^>" $sample*.fna | sed "s/.extracted.fna//;s/:>/\t/" | cut -f 1 -d " " > $sample.genes.txt; rm $sample*.fna'

cd -

# PREP FUNCTIONAL

# shrink SP100 cluster ids, remove irrelevant emapper annotation columns, and add func ann hash
cut -f 1,7- /g/scb2/bork/ckim/HGT_MGE_project_v3/1.combine_SPIRE/SP100/sp100_emapper_result.tsv | sed "s/ \+/_/g" | awk -v OFS='\t' '{ $1=gensub(/^SP100_0*([0-9]+)/, "\\1", "g", $1); print $0 }' | sed "s/'//g" | xargs -I{} -P8 sh -c $'cid=$(printf "{}\n" | cut -f 1); ann=$(printf "{}\n" | cut -f 2- | tr ";" "," | tr "\t" ";"); h=$(printf "%s\n" "$ann" | sha256sum - | cut -f 1 -d " "); ann=$(printf "%s\n" "$ann" | tr ";" "\t"); printf "%s\t%s\t%s\n" "$cid" "$h" "$ann"' > sp100_emapper_short.tsv





# remove empty genes -- no longer necessary
# awk -v OFS='\t' 'NF>2' pg3_genes_with_sp100_and_sp095.txt > pg3_genes_with_sp100_and_sp095.txt.valid && mv pg3_genes_with_sp100_and_sp095.txt.valid pg3_genes_with_sp100_and_sp095.txt

#SP095_000000000 1876    502313  SP100_0392236849;SP100_1049731408;SP100_1161872006

# 0       GCA_000005845_03294;GCA_000006665_04263;GCA_000006925_03662
