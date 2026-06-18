select bins.bin_name,CONCAT('k',contigs.kmer_size,'_',contigs.contig_ordinal,'_',genes.gene_ordinal) as gene_kid, genes.*, contigs.* from contigs join versions.spirev1_bins on versions.spirev1_bins.bin_id=contigs.bin_id join bins on bins.id=versions.spirev1_bins.bin_id join genes on contigs.id=genes.contig_id





select studies.id,studies.study_name,samples.sample_name,bins.bin_name from versions.spirev1_bins join bins on bins.id=versions.spirev1_bins.bin_id join samples on bins.sample_id=samples.id join studies on studies.id=samples.study_id where versions.spirev1_bins.included_in_spire





select bins.bin_name,CONCAT('k',contigs.kmer_size,'_',contigs.contig_ordinal,'_',genes.gene_ordinal) as gene_kid,genes.id from studies join samples on samples.study_id=studies.id join bins on bins.sample_id=samples.id join versions.spirev1_bins on versions.spirev1_bins.bin_id=bins.id join contigs on contigs.bin_id=bins.id join genes on genes.contig_id=contigs.id where studies.id = 2 and versions.spirev1_bins.included_in_spire
