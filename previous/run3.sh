#!/bin/bash

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
~
