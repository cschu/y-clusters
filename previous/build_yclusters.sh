#!/bin/bash

set -e -o pipefail


SPECI=$1
PG3_SPECI_DATA=/scratch/schudoma/pg3_sp95/previous/speci/${SPECI}.txt
SPIRE_SPECI_DATA=/scratch/schudoma/pg3_sp95/previous/mag_translate/speci/${SPECI}.txt

mkdir -p yclusters/${SPECI} yclusters/tmp

TARGET=yclusters/${SPECI}/${SPECI}.y095.txt

if [[ ! -f ${TARGET}.tmp ]]; then
	if [[ -f ${PG3_SPECI_DATA} ]]; then 
		awk -v OFS='\t' '{print $1,gensub(/(GCA_[0-9]+)_[0-9]+/, "\\1", "g", $1),$3 }' ${PG3_SPECI_DATA} > ${TARGET}.tmp
	fi
	if [[ -f ${SPIRE_SPECI_DATA} ]]; then
		awk -v OFS='\t' '{ printf("%s:%s\t%s\t%s\n", $1, $3, $2, $6) }' ${SPIRE_SPECI_DATA} >> ${TARGET}.tmp
	fi
fi

N_GENOMES=$(cut -f 2 ${TARGET}.tmp | uniq | sort -u | wc -l)

if [[ ! -f ${TARGET}.tmp.by_sp095 ]]; then 
	sort -T yclusters/tmp/ -k3,3 ${TARGET}.tmp > ${TARGET}.tmp.by_sp095
fi

if [[ ! -f yclusters/${SPECI}/${SPECI}.y095.txt.csizes ]]; then
	cut -f 3 ${TARGET}.tmp.by_sp095 | uniq -c | sed "s/^ \+//" | awk -v OFS='\t' '{ print $2,$1 }' > yclusters/${SPECI}/${SPECI}.y095.txt.csizes.tmp
	sort -T yclusters/tmp -k1,1 yclusters/${SPECI}/${SPECI}.y095.txt.csizes.tmp > yclusters/${SPECI}/${SPECI}.y095.txt.csizes
fi

if [[ ! -f ${TARGET}.tmp.by_sp095.with_csize ]]; then
	join -1 3 -2 1 -o 1.1,1.2,1.3,2.2 ${TARGET}.tmp.by_sp095 yclusters/${SPECI}/${SPECI}.y095.txt.csizes | tr " " "\t" > ${TARGET}.tmp.by_sp095.with_csize
fi


if [[ ! -f ${TARGET} ]]; then
	awk -v OFS='\t' -v n=${N_GENOMES} '!/^GCA/ { split($1,gene,":"); printf("%s:%s\tSP095_%09i\t%i\t???\t%i\t%.5f\n", $2, gene[2], $3, $4, n, $4/(n+0)); }' ${TARGET}.tmp.by_sp095.with_csize > ${TARGET}
	cut -f 2- ${TARGET} | uniq -c | sed "s/^ \+//" | awk -v OFS='\t' '{ print $2,$3,$1,$5,$6 }' | sort -k2,2gr -k1,1 > yclusters/${SPECI}/${SPECI}.cluster_info.txt
fi





# SAMN07501112.psa_megahit.psb_metabat2.00002:k141_3428_38	SP095_000099977	122	???	15327	6.52293


# GCA_000163455_00001     GCA_000163455   32736772
# GCA_000163455_00002     GCA_000163455   10961509
# GCA_000163455_00003     GCA_000163455   521863305
# 
# GD_MAG_1080091_9_24:k141_14085_24       SAMEA14100538.psa_megahit.psb_metabat2.00014    4893900
# GD_MAG_1080091_9_25:k141_14085_25       SAMEA14100538.psa_megahit.psb_metabat2.00014    479682024
# GD_MAG_1080091_9_2:k141_14085_2 SAMEA14100538.psa_megahit.psb_metabat2.00014    31708
# GD_MAG_1080091_9_3:k141_14085_3 SAMEA14100538.psa_megahit.psb_metabat2.00014    46023
# GD_MAG_1080091_9_4:k141_14085_4 SAMEA14100538.psa_megahit.psb_metabat2.00014    59182
# 
# 
# SAMN07982914.psa_megahit.psb_metabat2.00003.k141_7352_6	SP095_000000000	14374	78	15413	0.93259
# SAMN07982946.psa_megahit.psb_metabat2.00002.k141_1326_26	SP095_000000000	14374	78	15413	0.93259
# SAMEA4385485.psa_megahit.psb_metabat2.00014.k141_58882_22	SP095_000000000	14374	78	15413	0.93259
# SAMEA5617852.psa_megahit.psb_metabat2.00003.k141_7403_26	SP095_000000000	14374	78	15413	0.93259
