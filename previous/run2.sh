#!/bin/bash

SP095_ANNOTATED=SP095_members.short.tsv.linear.sorted.with_isolate_clusters.with_mag_clusters.with_hybrid.sorted


# get all SP095 hybrid clusters ...
cut -f 2- ${SP095_ANNOTATED} | uniq | cut -f 1 | uniq -c | sed "s/^ \+//" | awk -v OFS='\t' '$1>1 { print $2,"sp095_hybrid"; }' > sp095_hybrid_clusters.txt

# ... and join them to the SP095 data
join -1 2 -2 1 -a 1 -o 1.1,1.2,1.3,2.2 ${SP095_ANNOTATED} sp095_hybrid_clusters.txt | tr " " "\t" | awk -v OFS='\t' 'NF==3 { $4=$3 } NF==4 { $4=gensub(/sp095_/, "", "g", $4) } { print $0 }' > ${SP095_ANNOTATED}.with_sp095_type

# extract SP095 cluster types
cut -f 2,4 ${SP095_ANNOTATED}.with_sp095_type | uniq > sp095_clustertypes.txt


# count isolate/mag members per SP100 cluster
awk -v OFS='\t' '{split($2,a,";"); b["GC"]=0; b["GD"]=0; for (i in a) b[substr(a[i],1,2)]++; print $1,b["GC"],b["GD"];}' SP100_members.short.tsv > SP100_members.short.tsv.sizes


# add cluster type
awk -v OFS='\t' '$2==0 { if ($3==1) { $4="MU"; print $0; } else { $4="MM"; print $0; }; next; } $3==0 { if ($2==1) { $4="IU"; print $0; } else { $4="IM"; print $0; }; next; } { $4="HM"; print $0; }' SP100_members.short.tsv.sizes > SP100_members.short.tsv.sizes.with_ctype

# sort by sp100
sort -k1,1 SP100_members.short.tsv.sizes.with_ctype > SP100_members.short.tsv.sizes.with_ctype.sorted

# add sp100 size/ctypes to sp095
join -1 1 -2 1 SP095_members.short.tsv.linear.sorted SP100_members.short.tsv.sizes.with_ctype.sorted | tr " " "\t" > SP095_members.short.tsv.linear.sorted.with_sp100_sizes_ctype

# sort by sp095 id and sp100 ctype to prepare sp095 evaluation
sort -k2,2 -k5,5 SP095_members.short.tsv.linear.sorted.with_sp100_sizes_ctype > SP095_members.short.tsv.linear.sorted.with_sp100_sizes_ctype.sorted_for_sp095_eval

# annotate sp095 clusters
awk -v OFS='\t' 'BEGIN { cid=-1; ni=0; nm=0; a["MU"]=0; a["MM"]=0; a["IU"]=0; a["IM"]=0; a["HM"]=0; } { if (cid != $2) { if (cid != -1) { if (ni==0) { ctype="mag"; } else if (nm==0) { ctype="isolate"; } else { ctype="hybrid"; }; print cid,a["MU"]+a["MM"]+a["IU"]+a["IM"]+a["HM"],ni+nm,ni,nm,a["MU"],a["MM"],a["IU"],a["IM"],a["HM"],ctype; }; cid=$2; ni=0; nm=0; a["MU"]=0; a["MM"]=0; a["IU"]=0; a["IM"]=0; a["HM"]=0; } ni=(ni + $3); nm=(nm + $4); a[$5]++; } END { if (ni==0) { ctype="mag"; } else if (nm==0) { ctype="isolate"; } else { ctype="hybrid"; }; print cid,a["MU"]+a["MM"]+a["IU"]+a["IM"]+a["HM"],ni+nm,ni,nm,a["MU"],a["MM"],a["IU"],a["IM"],a["HM"],ctype; }' SP095_members.short.tsv.linear.sorted.with_sp100_sizes_ctype.sorted_for_sp095_eval > sp095_sizes_and_ctypes.txt
