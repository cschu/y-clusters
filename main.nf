params.sp100_members = "/g/scb2/bork/ckim/HGT_MGE_project_v3/1.combine_SPIRE/SP100/SP100_members.tsv"
params.sp095_members = "/g/scb2/bork/ckim/HGT_MGE_project_v3/2.linclust/SP095/SP095_members.tsv"
params.contig_map = "/g/scb2/bork/ckim/HGT_MGE_project_v3/data/SPIRE_contig_mapping.tsv"

params.spire_bins = "/g/bork6/schudoma/projects/mge/y_spire_bins.txt.bin_data"
params.spire_contigs = "/g/bork6/schudoma/projects/mge/spire_contigs.txt.y"
params.spire_genes = "/g/bork6/schudoma/projects/mge/spire_genes.txt.y10k"

params.spire_speci_info = "/g/bork6/tmp/schudoma/promge2_recovery/pg3_sp95/spire_speci_info.txt"



process create_contig2ycontig_map {
	input:
	path(contig_data)

	output:
	path("contig2ycontig.txt"), emit: ycontig_map

	script:
	"""
	set -e -o pipefail
	mkdir -p tmp/

	awk -v OFS='\\t' '{ printf("%s:%s\\t%s\\n", \$1, \$2, \$3); }' ${contig_data} | sort -T tmp/ -k1,1 > contig2ycontig.txt
	"""
}

process prepare_spire_contigs {
	// sort SPIRE contig database dump by SPIRE bin id
	// # bin_id, contig_id, kmer_size, ordinal
	// # 100008	1653229156	141	100713
	// add SPIRE bin data to contig data
	// generate joint-keys for adding in Y-MAG contig ids
	// add Y-MAG contig ids
	// sort contigs by SPIRE bin name
	// add specI to contigs
	// sort contigs by contig
	
	input:
	path(contigs)
	path(bins)
	path(ycontig_map)
	path(speci_info)

	output:
	path("spire_contigs.txt"), emit: contigs

	script:
	"""
	set -e -o pipefail
	mkdir -p tmp/

	awk -v OFS='\\t' '{ printf("%s\\t%s\\tk%s_%s\\n", \$1, \$2, \$3, \$4); }' ${contigs} > contigs.small
	sort -T tmp/ -k1,1 contigs.small > contigs.txt

	sort -T tmp/ -k1,1 ${bins} > bins.txt

	join -1 1 -2 1 -a 1 contigs.txt bins.txt | awk -v OFS='\\t' '{ printf("%s:%s\\t%s\\n", gensub(/.psa_megahit.psb_metabat2.[0-9]{5}\$/, "", "g", \$4), \$3, \$0); }' | sort -T tmp -k1,1 > contigs_with_bins.txt
	
	join -1 1 -2 1 contigs_with_bins.txt ${ycontig_map} | tr " " "\\t" | cut -f 2- | sort -T tmp/ -k4,4 > contigs_with_bins.with_ycontig.txt

	join -1 4 -2 1 contigs_with_bins.with_ycontig.txt ${speci_info} | tr " " "\\t" | sort -T tmp/ -k3,3 > spire_contigs.txt

	rm -fv contigs.small contigs.txt bins.txt contigs_with_bins.txt contigs_with_bins.with_ycontig.txt
	"""
}

process prepare_spire_genes {
	// sort SPIRE gene database dump by contig id
	input:
	path(genes)

	output:
	path("spire_genes.txt"), emit: genes

	script:
	"""
	set -e -o pipefail
	mkdir -p tmp/

	sort -T tmp/ -k2,2 ${genes} > spire_genes.txt
	"""
}

// process prepare_spire_bins {
// 	// sort SPIRE bin database dump by bin id
// 	input:
// 	path(bins)

// 	output:
// 	path("spire_bins.txt"), emit: bins

// 	script:
// 	"""
// 	set -e -o pipefail
// 	mkdir -p tmp/

// 	sort -T tmp/ -k1,1 ${bins} > spire_bins.txt
// 	"""
// }





process preprocess_sp100 {
	// this removes the count column 2 and the SP100_ prefixes including any padding zeroes
	// -> file size reduction -> quicker file ops
	// then, it splits sp100 into isolate and mag sets

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
	// splits sp100 <genometype> into isolate/non-isolate cluster sets

	input:
	tuple path(sp100_members), val(genome_type)

	output:
	tuple path("*.singletons.tsv"), val(genome_type), val("singleton"), emit: singletons
	tuple path("*.non_singletons.tsv"), val(genome_type), val("non_singleton"), emit: non_singletons

	script:
	"""
	set -e -o pipefail
	mkdir -p tmp/

	awk -F ';' \
	 'NF==1 { print \$0 >> "singletons.txt"; next; } \
	  NF==2 && \$2="" { print \$1 >> "singletons.txt"; next; } \
	  { print \$0 >> "non_singletons.txt"; } \
	  END { close("non_singletons.txt"); close("singletons.txt"); } \
	 ' ${sp100_members}

	touch singletons.txt non_singletons.txt

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
	path("SP095_members.short.linear.by_sp100.tsv"), emit: sp095

	script:
	"""
	set -e -o pipefail
	mkdir -p tmp/

	awk -v OFS='\\t' '{ print gensub(/^SP095_0*([0-9]+)/, "\\\\1", "g", \$1),gensub(/SP100_0*([0-9]+)/, "\\\\1", "g", \$4); }' ${sp095_members} > SP095_members.short.tsv
	awk -v OFS='\\t' '{ split(\$2,a,";"); for (i in a) { print a[i],\$1; }; }' SP095_members.short.tsv > SP095_members.short.tsv.linear
	sort -T tmp/ -k1,1 SP095_members.short.tsv.linear > SP095_members.short.linear.by_sp100.tsv

	rm -fv SP095_members.short.tsv SP095_members.short.tsv.linear
	"""
}

process add_sp095_clusters {
	input:
	tuple path(sp095_members), path(sp100_members), val(genome_type), val(cluster_type)	

	output:
	tuple path("SP100_members.${genome_type}.${cluster_type}.with_sp095.tsv"), val(genome_type), val(cluster_type), emit: sp100

	script:
	def csubtype_cmd = ""
	def sort_cmd = ""
	if (genome_type == "isolates") {
		if (cluster_type == "non_singleton") {
			csubtype_cmd += """awk -F '\\t' -v OFS='\\t' '{n=split(\$3,a,";"); if (n==2 && a[2]=="") ctype="S"; else ctype="N"; for (i in a) { if (a[i] != "") print a[i],\$1,\$2,ctype;}}'"""
			sort_cmd += """sort -k1,1 -k2,2 -k3,3 -T tmp/"""
		} else {
			csubtype_cmd += """awk -v OFS='\\t' '{ split(\$3,a,";"); for (i in a) { print a[i],\$1,\$2,"U"; }; }'"""
			sort_cmd += """sort -k1,1"""
		}
	} else {
		if (cluster_type == "non_singleton") {
			csubtype_cmd += """awk -F '\\t' -v OFS='\\t' '{n=split(\$3,a,";"); ctype="N"; for (i in a) { if (a[i] != "") print a[i],\$1,\$2,ctype;}}'"""
			sort_cmd += """sort -k1,1 -k2,2 -k3,3 -T tmp/"""
		} else {
			csubtype_cmd += """awk -v OFS='\\t' '{ split(\$3,a,";"); for (i in a) { print a[i],\$1,\$2,"S"; }; }'"""
			sort_cmd += """sort -k1,1"""
		}
	}

	"""
	set -e -o pipefail
	mkdir -p tmp/

	join -1 1 -2 1 -o 1.1,2.2,1.2 ${sp100_members} ${sp095_members} | tr " " "\t" > with_sp095.tsv

	${csubtype_cmd} with_sp095.tsv > with_sp095.tsv.1

	${sort_cmd} with_sp095.tsv.1 > SP100_members.${genome_type}.${cluster_type}.with_sp095.tsv

	rm -fv with_sp095.tsv*
	"""
}

process merge_isolate_clustertypes {
	input:
	path(files)

	output:
	tuple path("SP100_isolate_clusters.tsv"), val("isolates"), emit: sp100_isolates
	tuple path("SP100_isolate_clusters.tsv.by_sp100"), val("isolates"), emit: sp100_isolates_bysp100
	
	script:
	"""
	set -e -o pipefail
	mkdir -p tmp/

	sort -T tmp/ -m -k1,1 ${files} > SP100_isolate_clusters.tsv
	cut -f 2,4 SP100_isolate_clusters.tsv | sort -T tmp/ -k1,1 -u > SP100_isolate_clusters.tsv.by_sp100	
	"""
}

process correct_mag_singleton_clustertype {
	input:
	tuple path(mag_singletons), path(isolate_clusters_bysp100)

	output:
	path("SP100_members.mags.singletons.corrected_ctype.tsv"), emit: sp100_mag_singletons

	script:
	"""
	set -e -o pipefail
	mkdir -p tmp/

	sort -T tmp/ -k2,2 ${mag_singletons} > mag_singletons.by_sp100

	join -1 2 -2 1 -o 1.1,1.2,1.3,1.4,2.2 -a 1 mag_singletons.by_sp100 ${isolate_clusters_bysp100} | awk -v OFS='\\t' 'NF==4 { \$4="U" } { print \$0 }' | cut -f 1-4 > mag_singletons.by_gene

	sort -T tmp/ -k1,1 mag_singletons.by_gene > SP100_members.mags.singletons.corrected_ctype.tsv
	
	rm -fv mag_singletons.by_sp100 mag_singletons.by_gene
	"""
	// join -1 1 -2 1 -a 2 -o 2.1,2.2,2.3,1.2 SP100_clusters_with_isolates.tsv.with_ctype SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095 | tr " " "\t" | 
	// awk -v OFS='\t' 'NF==3 { $4="U" } { print $0 }' > SP100_members.short.no_isolates_no_isolate_singletons.mag_singletons.tsv.sorted.with_sp095.with_ctype
	// isolates: GCA_000003215_00074	341723245	234754752	U
	// mags: GD_MAG_0000000_10_12	1179883448	108643454	S
	
}

process merge_mag_clustertypes {
	input:
	path(files)

	output:
	tuple path("SP100_mag_clusters.tsv"), val("mags"), emit: sp100_mags
	// tuple path("SP100_isolate_clusters.tsv.by_sp100"), val("isolates"), emit: sp100_isolates_bysp100
	
	script:
	"""
	set -e -o pipefail
	mkdir -p tmp/

	sort -T tmp/ -m -k1,1 ${files} > SP100_mag_clusters.tsv
	"""
}



workflow {

	sp100_ch = Channel.fromPath(params.sp100_members)
	preprocess_sp100(sp100_ch)

	split_by_clustersize(preprocess_sp100.out.isolates.mix(preprocess_sp100.out.mags))

	sp095_ch = Channel.fromPath(params.sp095_members)
	preprocess_sp095(sp095_ch)

	contig_data_ch = Channel.fromPath(params.contig_map)
	create_contig2ycontig_map(contig_data_ch)

	// prepare_spire_bins(Channel.fromPath(params.spire_bins))
	prepare_spire_contigs(
		Channel.fromPath(params.spire_contigs),
		Channel.fromPath(params.spire_bins),
		create_contig2ycontig_map.out.ycontig_map,
		Channel.fromPath(params.spire_speci_info)
	)

	prepare_spire_genes(Channel.fromPath(params.spire_genes))

	add_sp095_clusters(
		preprocess_sp095.out.sp095.combine(split_by_clustersize.out.non_singletons.mix(split_by_clustersize.out.singletons)),		
	)

	merge_isolate_clustertypes(
		add_sp095_clusters.out.sp100
			.filter { it[1] == "isolates" }
			.map { it -> it[0] }
			.collect()
	)

	mag_singletons_ch = add_sp095_clusters.out.sp100
			.filter { it[1] == "mags" && it[2] == "singleton" }
			.map { it -> it[0] }
			.combine(merge_isolate_clustertypes.out.sp100_isolates_bysp100.map { it -> it[0] })

	mag_singletons_ch.dump(pretty: true, tag: "mag_singletons_ch")

	correct_mag_singleton_clustertype(
		mag_singletons_ch

		// tuple path(mag_singletons), path(isolate_clusters_bysp100)
	)

	merge_mag_clustertypes(
		add_sp095_clusters.out.sp100
			.filter { it[1] == "mags" && it[2] == "non_singleton" }
			.map { it -> it[0] }
			.mix(correct_mag_singleton_clustertype.out.sp100_mag_singletons)
			.collect()
	)

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