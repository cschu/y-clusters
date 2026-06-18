#!/bin/bash

set -e -o pipefail

SP100_MEMBERS=/g/scb2/bork/ckim/HGT_MGE_project/HGT_MGE_project_re/1.combine_SPIRE/SP100/SP100_members.tsv
SP095_MEMBERS=/g/scb2/bork/ckim/HGT_MGE_project/HGT_MGE_project_re/2.linclust/SP095_highcov/SP095_members.tsv

## PREPARE SP100

if [[ ! -s SP100_members.short.tsv ]]; then
	# echo "shrink sp100 cluster ids, remove cluster size column (2)"
	awk -v OFS='\t' '{ print gensub(/^SP100_0*([0-9]+)/, "\\1", "g", $1),$3; }' ${SP100_MEMBERS} > SP100_members.short.tsv
	chmod a-w SP100_members.short.tsv
fi

if [[ ! -f CLUSTER_SORT_STAGE_DONE ]]; then
	# echo "remove SPIRE (=MAG) genes and (now) empty SPIRE singleton clusters"
	sed "s/GD_MAG_[0-9]\+_[0-9]\+_[0-9]\+;\?//g" SP100_members.short.tsv | awk -v OFS='\t' 'NF>1 { print $0 }' | awk -F ';' 'NF==1 { print $0 >> "SP100_members.short.no_mags_no_mag_singletons.isolate_singletons.tsv"; next; } NF==2 && $2="" { print $1 >> "SP100_members.short.no_mags_no_mag_singletons.isolate_singletons.tsv"; next; } { print $0 >> "SP100_members.short.no_mags_no_mag_singletons.isolate_clusters.tsv" } END { close("SP100_members.short.no_mags_no_mag_singletons.isolate_singletons.tsv"); close("SP100_members.short.no_mags_no_mag_singletons.isolate_clusters.tsv") }'

	# echo sort sp100 singletons / clusters
	sort -k1,1 SP100_members.short.no_mags_no_mag_singletons.isolate_clusters.tsv > SP100_members.short.no_mags_no_mag_singletons.isolate_clusters.tsv.sorted &
	sort -k1,1 SP100_members.short.no_mags_no_mag_singletons.isolate_singletons.tsv > SP100_members.short.no_mags_no_mag_singletons.isolate_singletons.tsv.sorted &


	sed "s/GCA_[0-9]\+_[0-9]\+\+;\?//g" SP100_members.short.tsv | awk -v OFS='\t' 'NF>1 { print $0 }' | awk -F ';' 'NF==1 { print $0 >> "SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv"; next; } { print $0 >> "SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv" } END { close("SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv"); close("SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv") }'

	sort -k1,1 SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv > SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv.sorted &
	sort -k1,1 SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv > SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted &

	chmod a-w SP100_members.short.no_mags_no_mag_singletons.isolate_clusters.tsv.sorted SP100_members.short.no_mags_no_mag_singletons.isolate_singletons.tsv.sorted SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv.sorted SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted

	touch CLUSTER_SORT_STAGE_DONE

fi


## PREPARE SP095

if [[ ! -s SP095_members.short.tsv.linear.sorted ]]; then
	# echo "shrink sp095/sp100 cluster ids, remove size columns (2,3)"
	awk -v OFS='\t' '{ print gensub(/^SP095_0*([0-9]+)/, "\\1", "g", $1),gensub(/SP100_0*([0-9]+)/, "\\1", "g", $4); }' ${SP095_MEMBERS} > SP095_members.short.tsv
	chmod a-w SP095_members.short.tsv

	# echo "linearize sp095 -> sp100\\tsp095"
	awk -v OFS='\t' '{ split($2,a,";"); for (i in a) { print a[i],$1; }; }' SP095_members.short.tsv > SP095_members.short.tsv.linear

	# echo "sort sp095"
	sort -k1,1 SP095_members.short.tsv.linear > SP095_members.short.tsv.linear.sorted
	chmod a-w SP095_members.short.tsv.linear.sorted
fi


if [[ ! -f SP095_JOIN_DONE ]]; then
	# SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted
	# SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv.sorted
	# SP100_members.short.no_mags_no_mag_singletons.isolate_singletons.tsv.sorted
	# SP100_members.short.no_mags_no_mag_singletons.isolate_clusters.tsv.sorted
	echo "join sp095 -> sp100"
	for f in SP100_members.short.*.tsv.sorted; do
		join -1 1 -2 1 -o 1.1,2.2,1.2 $f SP095_members.short.tsv.linear.sorted | tr " " "\t" > $f.with_sp095 &
	done
	wait

	chmod a-w *.with_sp095
	touch SP095_JOIN_DONE
fi


if [[ ! -s SP100_members.short.no_mags_no_mag_singletons.isolate_clusters.tsv.sorted.with_sp095.by_gene ]]; then
	# echo "add gene cluster type and unravel non-singletons"
	awk -F '\t' -v OFS='\t' '{n=split($3,a,";"); if (n==2 && a[2]=="") ctype="S"; else ctype="N"; for (i in a) { if (a[i] != "") print a[i],$1,$2,ctype;}}' SP100_members.short.no_mags_no_mag_singletons.isolate_clusters.tsv.sorted.with_sp095 > SP100_members.short.no_mags_no_mag_singletons.isolate_clusters.tsv.sorted.with_sp095.by_gene
fi

if [[ ! -s SP100_members.short.no_mags_no_mag_singletons.isolate_singletons.tsv.sorted.with_sp095.by_gene ]]; then 
	# echo "add gene cluster type and unravel singletons"
	awk -v OFS='\t' '{ split($3,a,";"); for (i in a) { print a[i],$1,$2,"U"; }; }' SP100_members.short.no_mags_no_mag_singletons.isolate_singletons.tsv.sorted.with_sp095 > SP100_members.short.no_mags_no_mag_singletons.isolate_singletons.tsv.sorted.with_sp095.by_gene
fi

if [[ ! -s pg3_genes_with_sp100_and_sp095.txt ]]; then 
	# echo "sort isolate data"
	sort -k1,1 SP100_members.short.no_mags_no_mag_singletons.isolate_singletons.tsv.sorted.with_sp095.by_gene > SP100_members.short.no_mags_no_mag_singletons.isolate_singletons.tsv.sorted.with_sp095.by_gene.sorted &
	sort -k1,1 -k2,2 -k3,3 SP100_members.short.no_mags_no_mag_singletons.isolate_clusters.tsv.sorted.with_sp095.by_gene > SP100_members.short.no_mags_no_mag_singletons.isolate_clusters.tsv.sorted.with_sp095.by_gene.sorted &
	wait

	# echo "merge isolate data"
	sort -k1,1 -m SP100_members.short.no_mags_no_mag_singletons.isolate_clusters.tsv.sorted.with_sp095.by_gene.sorted SP100_members.short.no_mags_no_mag_singletons.isolate_singletons.tsv.sorted.with_sp095.by_gene.sorted > pg3_genes_with_sp100_and_sp095.txt
	chmod a-w pg3_genes_with_sp100_and_sp095.txt
fi

# echo add gene cluster type and unravel non-singletons
awk -F '\t' -v OFS='\t' '{n=split($3,a,";"); ctype="N"; for (i in a) { if (a[i] != "") print a[i],$1,$2,ctype;}}' SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv.sorted.with_sp095 > SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv.sorted.with_sp095.by_gene


#echo "sort non-singletons"
sort -k1,1 -k2,2 -k3,3 SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv.sorted.with_sp095.by_gene > SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv.sorted.with_sp095.by_gene.sorted

#echo "get all SP100 clusters that contain isolate genes"
cut -f 1 SP100_members.short.no_mags_no_mag_singletons.isolate_clusters.tsv.sorted > SP100_clusters_with_isolates.mc.tsv
cut -f 1 SP100_members.short.no_mags_no_mag_singletons.isolate_singletons.tsv.sorted > SP100_clusters_with_isolates.sc.tsv
sort -k1,1 -m SP100_clusters_with_isolates.mc.tsv SP100_clusters_with_isolates.sc.tsv > SP100_clusters_with_isolates.tsv
#echo "add "S"ingleton type column"
awk -v OFS='\t' '{print$0,"S"}' SP100_clusters_with_isolates.tsv > SP100_clusters_with_isolates.tsv.with_ctype

#echo "add ctype column to singletons"
join -1 1 -2 1 -a 2 -o 2.1,2.2,2.3,1.2 SP100_clusters_with_isolates.tsv.with_ctype SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095 | tr " " "\t" | awk -v OFS='\t' 'NF==3 { $4="U" } { print $0 }' > SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095.with_ctype

#echo "unravel singletons"
awk -v OFS='\t' '{ split($3,a,";"); for (i in a) { print a[i],$1,$2,$4; }; }' SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095.with_ctype > SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095.with_ctype.by_gene

#echo "sort singletons"
sort -k1,1 SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095.with_ctype.by_gene > SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095.with_ctype.by_gene.sorted

#echo "merge"
sort -k1,1 -m SP100_members.short.no_isolates_no_isolate_singletons.tsv.mag_clusters.tsv.sorted.with_sp095.by_gene.sorted SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095.with_ctype.by_gene.sorted > spire_genes_with_sp100_and_sp095.txt

cut -f 3- -d _ spire_genes_with_sp100_and_sp095.txt | cut -f 1,3 > spire_genes_with_sp100_and_sp095.txt.gene_sp095only



















