import "nmdc-metag_full.wdl" as metag

workflow validate_finish {
   String input_bundle="https://portal.nersc.gov/project/m3408/test_data/validate_imports.zip"   
   String container="microbiomedata/workflowmeta:1.0.6"
   String proj="test"
   String informed_by="test"
   String resource="NERSC"
   String url_root="http"
   String git_url="test"
   String outdir="/tmp/testout"
  call stage {
	input: imports=input_bundle, 
	   container="microbiomedata/workflowmeta:1.0.6"
	}

  call metag.finish as fin {
	input: container="microbiomedata/workflowmeta:1.0.6",
           proj=proj,
           start=stage.start,
           resource=resource,
           url_root=url_root,
           git_url=git_url,
           informed_by=informed_by,
           read=stage.read,
           filtered=stage.filtered,
           filtered_stats=stage.filtered_stats,
           filtered_stats2=stage.filtered_stats2,
           fasta=stage.fasta,
           scaffold=stage.scaffold,
           agp=stage.agp,
           bam=stage.bam,
           samgz=stage.samgz,
           covstats=stage.covstats,
           asmstats=stage.asmstats,
           proteins_faa=stage.proteins_faa,
           functional_gff=stage.functional_gff,
           structural_gff=stage.structural_gff,
           ko_tsv=stage.ko_tsv,
           ec_tsv=stage.ec_tsv,
           cog_gff=stage.cog_gff,
           pfam_gff=stage.pfam_gff,
           tigrfam_gff=stage.tigrfam_gff,
           smart_gff=stage.smart_gff,
           supfam_gff=stage.supfam_gff,
           cath_funfam_gff=stage.cath_funfam_gff,
           crt_gff=stage.crt_gff,
           genemark_gff=stage.genemark_gff,
           prodigal_gff=stage.prodigal_gff,
           trna_gff=stage.trna_gff,
           misc_bind_misc_feature_regulatory_gff=stage.misc_bind_misc_feature_regulatory_gff,
           rrna_gff=stage.rrna_gff,
           ncrna_tmrna_gff=stage.ncrna_tmrna_gff,
           crt_crisprs=stage.crt_crisprs,
           product_names_tsv=stage.product_names_tsv,
           gene_phylogeny_tsv=stage.gene_phylogeny_tsv,
           ko_ec_gff=stage.ko_ec_gff,
           stats_tsv=stage.stats_tsv,
           stats_json=stage.stats_json,
           gottcha2_report_tsv=stage.gottcha2_report_tsv,
           gottcha2_full_tsv=stage.gottcha2_full_tsv,
           gottcha2_krona_html=stage.gottcha2_krona_html,
           centrifuge_classification_tsv=stage.centrifuge_classification_tsv,
           centrifuge_report_tsv=stage.centrifuge_report_tsv,
           centrifuge_krona_html=stage.centrifuge_krona_html,
           kraken2_classification_tsv=stage.kraken2_classification_tsv,
           kraken2_report_tsv=stage.kraken2_report_tsv,
           kraken2_krona_html=stage.kraken2_krona_html,
           short=stage.short,
           lowdepth=stage.lowdepth,
           unbinned=stage.unbinned,
           checkm=stage.checkm,
           hqmq_bin_fasta_files=stage.hqmq_bin_fasta_files,
           bin_fasta_files=stage.bin_fasta_files,
           mags_stats_json=stage.mags_stats_json,
           outdir=outdir
  }
}

task stage {

    String imports
    String container

    command {
         wget ${imports}
	 unzip validate_imports.zip
         date --iso-8601=seconds > start.txt

    }

    output {
         String start = read_string("start.txt")
   	 File read="nmdc_filtered.fastq.gz"
   	 File filtered="nmdc_filtered.fastq.gz"
	 File? filtered_stats="nmdc_filterStats.txt"
	 File? filtered_stats2="nmdc_filterStats2.txt"
	 File fasta="nmdc_contigs.fna"
	 File scaffold="nmdc_contigs.fna"
	 File? agp="nmdc_assembly.agp"
	 File bam="nmdc_pairedMapped_sorted.bam"
	 File? samgz="nmdc_pairedMapped.sam.gz"
	 File? covstats="nmdc_covstats.txt"
	 File asmstats="nmdc_asm_stats.json"
	 File proteins_faa="nmdc_proteins.faa"
	 File structural_gff="nmdc_structural_annotation.gff"
	 File functional_gff="nmdc_functional_annotation.gff"
	 File ko_tsv="nmdc_ko.tsv"
	 File ec_tsv="nmdc_ec.tsv"
	 File cog_gff="nmdc_cog.gff"
	 File pfam_gff="nmdc_pfam.gff"
	 File tigrfam_gff="nmdc_tigrfam.gff"
	 File smart_gff="nmdc_smart.gff"
	 File supfam_gff="nmdc_supfam.gff"
	 File cath_funfam_gff="nmdc_cath_funfam.gff"
	 File crt_gff="nmdc_crt.gff"
         File genemark_gff="nmdc_genemark.gff"
         File prodigal_gff="nmdc_prodigal.gff"
	 File trna_gff="nmdc_trna.gff"
	 File misc_bind_misc_feature_regulatory_gff="nmdc_rfam_misc_bind_misc_feature_regulatory.gff"
	 File rrna_gff="nmdc_rfam_rrna.gff"
	 File ncrna_tmrna_gff="nmdc_rfam_ncrna_tmrna.gff"
	 File ko_ec_gff="nmdc_ko_ec.gff"
	 File? stats_tsv="nmdc_stats.tsv"
	 File stats_json="nmdc_stats.json"
	 File gene_phylogeny_tsv="nmdc_gene_phylogeny.tsv"
	 File product_names_tsv="nmdc_product_names.tsv"
	 File crt_crisprs="nmdc_crt.crisprs"
	 File? short="nmdc_bins.tooShort.fa"
	 File? lowdepth="nmdc_bins.lowDepth.fa"
	 File? unbinned="nmdc_bins.unbinned.fa"
	 File? checkm="nmdc_checkm_qa.out"
	 Array[File] hqmq_bin_fasta_files=["nmdc_bin.fa"]
	 Array[File] bin_fasta_files=["nmdc_bin.fa"]
	 File? mags_stats_json="nmdc_mags_stats.json"
	 File? gottcha2_report_tsv="nmdc_kraken2_report.tsv"
	 File? gottcha2_full_tsv="nmdc_gottcha2_report_full.tsv"
	 File? gottcha2_krona_html="nmdc_gottcha2_krona.html"
	 File? centrifuge_classification_tsv="nmdc_centrifuge_classification.tsv"
	 File? centrifuge_report_tsv="nmdc_centrifuge_report.tsv"
	 File? centrifuge_krona_html="nmdc_centrifuge_krona.html"
	 File? kraken2_classification_tsv="nmdc_kraken2_classification.tsv"
	 File? kraken2_report_tsv="nmdc_kraken2_report.tsv"
	 File? kraken2_krona_html="nmdc_kraken2_krona.html"
	}

   runtime {
     memory: "10 GiB"
     cpu:  4
     maxRetries: 1
     docker: container
   }
}
