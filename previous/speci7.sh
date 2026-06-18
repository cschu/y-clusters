#!/bin/bash



awk -v OFS='\t' '{print $1,gensub(/(GCA_[0-9]+)_[0-9]+/, "\\1", "g", $1),$3 }' speci/specI_v4_00007.txt > speci_7.txt
awk -v OFS='\t' '{ printf("%s:%s\t%s\t%s\n", $1, $3, $2, $6) }' mag_translate/speci/specI_v4_00007.txt >> speci_7.txt
