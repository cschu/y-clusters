params.sp100_members = "/g/scb2/bork/ckim/HGT_MGE_project_v3/1.combine_SPIRE/SP100/SP100_members.tsv"
params.sp095_members = "/g/scb2/bork/ckim/HGT_MGE_project_v3/2.linclust/SP095/SP095_members.tsv"


process preprocess_sp100 {
	// this removes the count column 2 and the SP100_ prefixes including any padding zeroes
	// -> file size reduction -> quicker file ops

	input:
	path(sp100_members)

	output:
	path("SP100_members.short.tsv"), emit: sp100

	script:
	"""
	set -e -o pipefail
	awk -v OFS='\\t' '{ print gensub(/^SP100_0*([0-9]+)/, "\\\\1", "g", \$1),\$3; }' ${sp100_members} > SP100_members.short.tsv
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

	sp095_ch = Channel.fromPath(params.sp095_members)
	preprocess_sp095(sp095_ch)

}