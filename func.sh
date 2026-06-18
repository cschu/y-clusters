#!/bin/bash

cut -f 1,7- /g/scb2/bork/ckim/HGT_MGE_project_v3/1.combine_SPIRE/SP100/sp100_emapper_result.tsv | sed "s/ \+/_/g" | awk -v OFS='\t' '{ $1=gensub(/^SP100_0*([0-9]+)/, "\\1", "g", $1); print $0 }' | sed "s/'//g" | xargs -I{} -P8 sh -c $'cid=$(printf "{}\n" | cut -f 1); ann=$(printf "{}\n" | cut -f 2- | tr ";" "," | tr "\t" ";"); h=$(printf "%s\n" "$ann" | sha256sum - | cut -f 1 -d " "); ann=$(printf "%s\n" "$ann" | tr ";" "\t"); printf "%s\t%s\t%s\n" "$cid" "$h" "$ann"' > sp100_emapper_short.tsv
