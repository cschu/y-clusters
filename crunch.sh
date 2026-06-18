#!/bin/bash

NCPUS=${SLURM_CPUS_PER_TASK}

find . -maxdepth 1 -mindepth 1 -type l -name '*.genes.txt' | xargs -I{} -P${NCPUS} sh -c $'if [[ ! -f {}.LOCK ]] && [[ -f {} ]]; then touch {}.LOCK; sample=$(basename {} .genes.txt); echo $sample; awk -v OFS=\'\t\' \'{ split($2,a,"_"); print $0,a[1]"_"a[2],a[3]; }\' {} | sort -k3,3 -k4,4n > {}.tmp; awk -v OFS=\'\t\' -v sample=$sample \'$1==sample\' ../spire_contigs.txt > {}.contigs.tmp; if [[ -s {}.tmp ]] && [[ -s {}.contigs.tmp ]] && [[ ! -s DONE/{}.with_contigs.txt ]]; then join -1 3 -2 2 -o 1.1,1.2,2.3,2.4,2.5,2.6 {}.tmp {}.contigs.tmp | tr " " "\t" > DONE/{}.with_contigs.txt; rm -fv {}.tmp {}.contigs.tmp; rm -v {}; rm {}.LOCK; fi; fi'

exit
find . -maxdepth 1 -mindepth 1 -name '*.genes.txt' | xargs -I{} -P${NCPUS} sh -c \
	$'set -e -o pipefail; \
	if [[ ! -f {}.LOCK ]]; then \
		touch {}.LOCK; sample=$(basename {} .genes.txt); \
		echo $sample; \
		awk -v OFS=\'\t\' \'{ split($2,a,"_"); print $0,a[1]"_"a[2],a[3]; }\' {} | sort -k3,3 -k4,4n > {}.tmp; \
		awk -v OFS=\'\t\' -v sample=$sample \'$1==sample\' ../spire_contigs.txt > {}.contigs.tmp; \
		join -1 3 -2 2 -o 1.1,1.2,2.3,2.4,2.5,2.6 {}.tmp {}.contigs.tmp | tr " " "\t" > DONE/{}.with_contigs.txt; \
		rm -fv {}.tmp {}.contigs.tmp; \
		mv -v {} DONE/; \
		rm {}.LOCK; \
	fi'
