import "rqcfilter.wdl" as rqc
import "jgi_assembly.wdl" as assembly
import "annotation_full.wdl" as awf
import "ReadbasedAnalysis.wdl" as rba
import "mbin_nmdc.wdl" as mags


workflow nmdc_metag {
  String  container="bfoster1/img-omics:0.1.9"
  String  proj
  String  input_file
  String  outdir
  String  database="/refdata/img/"
  String  resource
  String  informed_by
  String?  git_url="https://github.com/microbiomedata/mg_annotation/releases/tag/0.1"
  String?  url_root="https://data.microbiomedata.org/data/"
  String  url_base="${url_root}${proj}/annotation/"

  call stage {
    input: container=container,
           input_file=input_file,
           proj=proj
  }
  # Estimate RQC runtime at an hour per compress GB
  call rqc.jgi_rqcfilter as qc {
    input: input_files=[stage.read],
           outdir="${outdir}/qa/",
           threads=16,
           memory="60G"
  }
  call assembly.jgi_metaASM as asm {
    input: input_file=qc.filtered,
           rename_contig_prefix=proj,
           outdir="${outdir}/assembly/"
  }

  call awf.annotation {
    input: imgap_project_id=stage.pref,
           imgap_input_fasta=asm.contig,
           database_location=database
  }
  call split_interleaved_fastq {
    input:
      reads=qc.filtered[0],
      container="microbiomedata/bbtools:38.90"
  }

  call rba.ReadbasedAnalysis {
    input: enabled_tools = { "gottcha2": true, "kraken2": true, "centrifuge": true},
           db = { "gottcha2": "/refdata/gottcha2/RefSeq-r90.cg.BacteriaArchaeaViruses.species.fna",
                  "kraken2": "/refdata/kraken2",
                  "centrifuge": "/refdata/centrifuge/p_compressed"
                },
           reads = split_interleaved_fastq.outFastq,
           paired = true,
           prefix = stage.pref,
           outdir = "${outdir}/ReadbasedAnalysis",
           cpu = 4
  }

  call mags.nmdc_mags {
    input:
      proj_name=proj,
      contig_file=asm.contig,
      sam_file=asm.bam,
      gff_file=annotation.functional_gff,
      gtdbtk_database="/refdata/GTDBTK_DB",
      outdir="${outdir}/MAGs",
      container="microbiomedata/nmdc_mbin:0.1.2"
  }


  call finish {
    input: container="microbiomedata/workflowmeta:1.0.0",
           start=stage.start,
           resource=resource,
           url_base=url_base,
           git_url=git_url,
           informed_by=informed_by,
           read = stage.read,
           filtered = qc.filtered[0],
           filtered_stats = qc.stats[0],
           fasta=asm.contig,
           scaffold=asm.scaffold,
           agp=asm.agp,
           bam=asm.bam,
           covstats=asm.covstats,
           proteins_faa=annotation.proteins_faa,
           functional_gff=annotation.functional_gff,
           structural_gff=annotation.structural_gff,
           ko_tsv=annotation.ko_tsv,
           ec_tsv=annotation.ec_tsv,
           cog_gff=annotation.cog_gff,
           pfam_gff=annotation.pfam_gff,
           tigrfam_gff=annotation.tigrfam_gff,
           smart_gff=annotation.smart_gff,
           supfam_gff=annotation.supfam_gff,
           cath_funfam_gff=annotation.cath_funfam_gff,
           ko_ec_gff=annotation.ko_ec_gff,
           stats_tsv=annotation.stats_tsv,
           stats_json=annotation.stats_json,
           gottcha2_report_tsv = ReadbasedAnalysis.gottcha2_report_tsv,
           gottcha2_full_tsv = ReadbasedAnalysis.gottcha2_full_tsv,
           gottcha2_krona_html = ReadbasedAnalysis.gottcha2_krona_html,
           centrifuge_classification_tsv = ReadbasedAnalysis.centrifuge_classification_tsv,
           centrifuge_report_tsv = ReadbasedAnalysis.centrifuge_report_tsv,
           centrifuge_krona_html = ReadbasedAnalysis.centrifuge_krona_html,
           kraken2_classification_tsv = ReadbasedAnalysis.kraken2_classification_tsv,
           kraken2_report_tsv = ReadbasedAnalysis.kraken2_report_tsv,
           kraken2_krona_html = ReadbasedAnalysis.kraken2_krona_html,
           short=nmdc_mags.short,
           lowdepth=nmdc_mags.low,
           unbinned=nmdc_mags.unbinned,
           checkm=nmdc_mags.checkm,
           outdir=outdir
  }

  meta {
    author: "Shane Canon"
    email: "scanon@lbl.gov"
    version: "1.0.0"
  }
}

task stage {
   String container
   String proj
   String prefix=sub(proj, ":", "_")
   String target="${prefix}.fastq.gz"
   String input_file

   command <<<
       set -e
       if [ $( echo ${input_file}|egrep -c "https*:") -gt 0 ] ; then
           wget ${input_file} -O ${target}
       else
           ln ${input_file} ${target} || cp ${input_file} ${target}
       fi
       # Capture the start time
       date --iso-8601=seconds > start.txt

   >>>

   output{
      File read = "${target}"
      String start = read_string("start.txt")
      String pref = "${prefix}"
   }
   runtime {
     memory: "1 GiB"
     cpu:  2
     maxRetries: 1
     docker: container
   }
}

task split_interleaved_fastq{
    File reads
    String container
    String? memory = "4G"
    String output1 = "input.left.fastq.gz"
    String output2 = "input.right.fastq.gz"

    runtime {
        docker: container
        mem: "4 GiB"
        cpu:  1
    }
    command {
         reformat.sh -Xmx${default="10G" memory} in=${reads} out1=${output1} out2=${output2}
    }

    output {
            Array[File] outFastq = [output1, output2]
    }
}



task finish {
   String container
   String start
   String informed_by
   String resource
   String url_base
   String git_url
   File read
   File filtered
   File filtered_stats
   File fasta
   File scaffold
   File agp
   File bam
   File covstats
   File proteins_faa
   File structural_gff
   File functional_gff
   File ko_tsv
   File ec_tsv
   File cog_gff
   File pfam_gff
   File tigrfam_gff
   File smart_gff
   File supfam_gff
   File cath_funfam_gff
   File ko_ec_gff
   File stats_tsv
   File stats_json
# Future
#    File gene_phylogeny_tsv
#    File proteins_cog_domtblout
#    File proteins_pfam_domtblout
#    File proteins_tigrfam_domtblout
#    File proteins_smart_domtblout
#    File proteins_supfam_domtblout
#    File proteins_cath_funfam_domtblout
#    File product_names_tsv
#    File crt_crisprs
   File short
   File lowdepth
   File unbinned
   File checkm
   File? gottcha2_report_tsv
   File? gottcha2_full_tsv
   File? gottcha2_krona_html
   File? centrifuge_classification_tsv
   File? centrifuge_report_tsv
   File? centrifuge_krona_html
   File? kraken2_classification_tsv
   File? kraken2_report_tsv
   File? kraken2_krona_html
   String outdir
   String qadir="${outdir}/qa/"
   String assemdir="${outdir}/assembly/"
   String annodir="${outdir}/annotation/"
   String magsdir="${outdir}/MAGs/"
   String rbadir="${outdir}/ReadbasedAnalysis/"

   command{
       set -e
       mkdir -p ${annodir}
       end=`date --iso-8601=seconds`

       # Generate QA objects
       /scripts/rqcstats.py ${filtered_stats} > stats.json
       /scripts/generate_objects.py --type "qa" --id ${informed_by} \
             --start ${start} --end $end \
             --resource '${resource}' --url ${url_base} --giturl ${git_url} \
             --extra stats.json \
             --inputs ${read} \
             --outputs \
             ${filtered} 'Filtered Reads' \
             ${filtered_stats} 'Filtered Stats'
       cp activity.json data_objects.json ${qadir}/

       # Generate assembly objects
       /scripts/generate_objects.py --type "assembly" --id ${informed_by} \
             --start ${start} --end $end \
             --resource '${resource}' --url ${url_base} --giturl ${git_url} \
             --inputs ${filtered} \
             --outputs \
             ${covstats} 'Metagenome Contig Coverage Stats' \
             ${fasta} 'Assembled contigs fasta' \
             ${scaffold} 'Assembled scaffold fasta' \
             ${agp} 'Assembled AGP file' \
             ${bam} 'Metagenome Alignment BAM file'
       cp activity.json data_objects.json ${assemdir}/

       # Generate annotation objects
       nmdc gff2json ${functional_gff} -of features.json -oa annotations.json -ai ${informed_by}

       /scripts/generate_objects.py --type "annotation" --id ${informed_by} \
             --start ${start} --end $end \
             --resource '${resource}' --url ${url_base} --giturl ${git_url} \
             --inputs ${fasta} \
             --outputs \
             ${proteins_faa} 'Protein FAA' \
             ${structural_gff} 'Structural annotation GFF file' \
             ${functional_gff} 'Functional annotation GFF file' \
             ${ko_tsv} 'KO TSV file' \
             ${ec_tsv} 'EC TSV file' \
             ${cog_gff} 'COG GFF file' \
             ${pfam_gff} 'PFAM GFF file' \
             ${tigrfam_gff} 'TigrFam GFF file' \
             ${smart_gff} 'SMART GFF file' \
             ${supfam_gff} 'SuperFam GFF file' \
             ${cath_funfam_gff} 'Cath FunFam GFF file' \
             ${ko_ec_gff} 'KO_EC GFF file'

       cp ${proteins_faa} ${structural_gff} ${functional_gff} \
          ${ko_tsv} ${ec_tsv} ${cog_gff} ${pfam_gff} ${tigrfam_gff} \
          ${smart_gff} ${supfam_gff} ${cath_funfam_gff} ${ko_ec_gff} \
          ${stats_tsv} ${stats_json} \
          ${annodir}/
       cp features.json annotations.json activity.json data_objects.json ${annodir}/


       /scripts/generate_objects.py --type "MAGs" --id ${informed_by} \
             --start ${start} --end $end \
             --resource '${resource}' --url ${url_base} --giturl ${git_url} \
             --inputs ${fasta} ${bam} ${functional_gff} \
             --outputs \
             ${short} "tooShort (< 3kb) filtered contigs fasta file by metaBat2" \
             ${lowdepth} "lowDepth (mean cov <1 )  filtered contigs fasta file by metabat2" \
             ${unbinned} "unbinned fasta file from metabat2" \
             ${checkm} "metabat2 bin checkm quality assessment result"
#             "gtdbtk.bac120.summary.tsv" "gtdbtk bacterial assignment result summary table" \
#             "gtdbtk.ar122.summary.tsv" "gtdbtk archaea assignment result summary table" \
       cp activity.json data_objects.json ${magsdir}/

       /scripts/generate_objects.py --type "ReadbasedAnalysis" --id ${informed_by} \
             --start ${start} --end $end \
             --resource '${resource}' --url ${url_base} --giturl ${git_url} \
             --inputs ${filtered} \
             --outputs \
             ${gottcha2_report_tsv} "Gottcha2 TSV report" \
             ${gottcha2_full_tsv} "Gottcha2 full TSV report" \
             ${gottcha2_krona_html} "Gottcha2 Krona HTML report" \
             ${centrifuge_classification_tsv} "Centrifuge classification TSV report" \
             ${centrifuge_report_tsv} "Centrifuge TSV report" \
             ${centrifuge_krona_html} "Centrifuge Krona HTML report" \
             ${kraken2_classification_tsv} "Kraken classification TSV report" \
             ${kraken2_report_tsv} "Kraken2 TSV report" \
             ${kraken2_krona_html} "Kraken2 Krona HTML report"
       cp activity.json data_objects.json ${rbadir}/
   }

   runtime {
     memory: "10 GiB"
     cpu:  4
     maxRetries: 1
     docker: container
   }
}


