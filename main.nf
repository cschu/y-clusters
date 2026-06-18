params.sp100_members = "/g/scb2/bork/ckim/HGT_MGE_project_v3/1.combine_SPIRE/SP100/SP100_members.tsv"
params.sp095_members = "/g/scb2/bork/ckim/HGT_MGE_project_v3/2.linclust/SP095/SP095_members.tsv"


process preprocess_sp100 {
	// this removes the count column 2 and the SP100_ prefixes including any padding zeroes
	// -> file size reduction -> quicker file ops

	input:
	path(sp100_members)

	output:
	path("SP100_members.short.tsv"), emit: sp100
	tuple path("SP100_members.short.tsv.isolates"), val("isolates"), emit: isolates
	tuple path("SP100_members.short.tsv.mags"), val("mags"), emit: mags

	script:
	"""
	set -e -o pipefail
	awk -v OFS='\\t' '{ print gensub(/^SP100_0*([0-9]+)/, "\\\\1", "g", \$1),\$3; }' ${sp100_members} > SP100_members.short.tsv

	sed "s/GD_MAG_[0-9]\\+_[0-9]\\+_[0-9]\\+;\\?//g" SP100_members.short.tsv | awk -v OFS='\\t' 'NF>1' > SP100_members.short.tsv.isolates
	sed "s/GCA_[0-9]\\+_[0-9]\\+\\+;\\?//g" SP100_members.short.tsv | awk -v OFS='\\t' 'NF>1' > SP100_members.short.tsv.mags
	"""

}

process split_by_clustersize {
	input:
	path(sp100_members)
	val(genome_type)

	output:
	path("*.singletons.tsv"), emit: singletons
	path("*.non_singletons.tsv"), emit: non_singletons

	script:
	"""
	set -e -o pipefail
	mkdir -p tmp/

	awk -F ';' \
	 'NF==1 { print \$0 >> "singletons.txt"; next; } \
	  NF==2 && \$2="" { print \$1 >> "singletons.txt"; next; } \
	  { print \$0 >> "non_singletons.txt"; } \
	  END { close("non_singletons.txt"); close("singletons.txt"); }

	sort -T tmp/ -k1,1 singletons.txt > SP100_members.${genome_type}.singletons.tsv
	sort -T tmp/ -k1,1 non_singletons.txt > SP100_members.${genome_type}.non_singletons.tsv

	rm -fv singletons.txt non_singletons.txt
	"""

}







process preprocess_sp095 {
	// this removes the count columns 2,3 and the SP100_/SP095_ prefixes including any padding zeroes,
	// then linearises the table -> sp095: sp100_a, sp100_b -> sp095, sp100_a; sp095, sp100_b
	// and sorts by sp095

	input:
	path(sp095_members)

	output:
	path("SP095_members.short.linear.by_sp095.tsv"), emit: sp095

	script:
	"""
	set -e -o pipefail
	mkdir -p tmp/

	awk -v OFS='\\t' '{ print gensub(/^SP095_0*([0-9]+)/, "\\\\1", "g", \$1),gensub(/SP100_0*([0-9]+)/, "\\\\1", "g", \$4); }' ${sp095_members} > SP095_members.short.tsv
	awk -v OFS='\\t' '{ split(\$2,a,";"); for (i in a) { print a[i],\$1; }; }' SP095_members.short.tsv > SP095_members.short.tsv.linear
	sort -T tmp/ -k1,1 SP095_members.short.tsv.linear > SP095_members.short.linear.by_sp095.tsv

	rm -fv SP095_members.short.tsv SP095_members.short.tsv.linear
	"""
	



}





workflow {

	sp100_ch = Channel.fromPath(params.sp100_members)
	preprocess_sp100(sp100_ch)

	// remove_mag_genes(preprocess_sp100.out.sp100)
	split_by_clustersize(preprocess_sp100.out.isolates.mix(preprocess_sp100.out.mags))

	sp095_ch = Channel.fromPath(params.sp095_members)
	preprocess_sp095(sp095_ch)

}



// process remove_mag_genes {
// 	path(sp100_members)

// 	output:
// 	path("SP100_members.isolate.singletons"), emit: sp100_isolate_singletons
// 	path("SP100_members.isolate.non_singletons"), emit: sp100_isolate_non_singletons

// 	script:
// 	"""
// 	set -e -o pipefail
// 	mkdir -p tmp/

// 	sed "s/GD_MAG_[0-9]\\+_[0-9]\\+_[0-9]\\+;\\?//g" ${sp100_members} | awk -v OFS='\\t' 'NF>1' > SP100_members.short.tsv.no_mags
	
// 	awk -F ';' \
// 		'NF==1 { print \$0 >> "SP100_members.isolate.singletons.1"; next; } \
// 		 NF==2 && $2="" { print \$1 >> "SP100_members.isolate.singletons.1"; next; } \
// 		 { print \$0 >> "SP100_members.isolate.non_singletons.1"; } \
// 		 END { close("SP100_members.isolate.non_singletons.1"); close("SP100_members.isolate.singletons.1"); } \
// 		'

// 	sort -T tmp/ -k1,1 SP100_members.isolate.singletons.1 > SP100_members.isolate.singletons
// 	sort -T tmp/ -k1,1 SP100_members.isolate.non_singletons.1 > SP100_members.isolate.non_singletons

// 	rm -fv SP100_members.short.tsv.no_mags SP100_members.isolate.singletons.1 SP100_members.isolate.non_singletons.1
// 	"""

// }