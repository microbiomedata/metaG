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
           input_file=input_file
  }
  # Estimate RQC runtime at an hour per compress GB
  call rqc.jgi_rqcfilter as qc {
    input: input_files=[stage.read],
           threads=16,
           memory="60G"
  }
  call assembly.jgi_metaASM as asm {
    input: input_file=qc.filtered
  }

  call awf.annotation {
    input: imgap_project_id="nmdc_",
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
           prefix = "nmdc_",
           cpu = 4
  }

  call mags.nmdc_mags {
    input:
      proj_name="nmdc",
      contig_file=asm.contig,
      sam_file=asm.bam,
      gff_file=annotation.functional_gff,
      gtdbtk_database="/refdata/GTDBTK_DB",
      container="microbiomedata/nmdc_mbin:0.1.2"
  }


  call finish {
    input: container="microbiomedata/workflowmeta:1.0.4",
           proj=proj,
           start=stage.start,
           resource=resource,
           url_root=url_root,
           git_url=git_url,
           informed_by=informed_by,
           read = stage.read,
           filtered = qc.filtered[0],
           filtered_stats = qc.stats[0],
           filtered_stats2 = qc.stats2[0],
           fasta=asm.contig,
           scaffold=asm.scaffold,
           agp=asm.agp,
           bam=asm.bam,
           samgz=asm.samgz,
           covstats=asm.covstats,
           asmstats=asm.asmstats,
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
           hqmq_bin_fasta_files=nmdc_mags.hqmq_bin_fasta_files,
           bin_fasta_files=nmdc_mags.bin_fasta_files,
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
   String target="raw.fastq.gz"
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
   String proj
   String prefix=sub(proj, ":", "_")
   String start
   String informed_by
   String resource
   String url_root
   String git_url
   File read
   File filtered
   File filtered_stats
   File filtered_stats2
   File fasta
   File scaffold
   File agp
   File bam
   File samgz
   File covstats
   File asmstats
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
   Array[File] hqmq_bin_fasta_files
   Array[File] bin_fasta_files
   Int n_hqmq=length(hqmq_bin_fasta_files)
   Int n_bin=length(bin_fasta_files)
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
   String orig_prefix="scaffold"
   String sed="s/${orig_prefix}_/${proj}_/g"

   command{
       set -e
       end=`date --iso-8601=seconds`
       # cleanup any previous runs
       rm -rf ${qadir}
       rm -rf ${assemdir}
       rm -rf ${annodir}
       rm -rf ${magsdir}
       rm -rf ${rbadir}

       # Generate QA objects
       mkdir -p ${qadir}
       cp ${filtered} ${qadir}/${prefix}_filtered.fastq.gz
       cp ${filtered_stats} ${qadir}/${prefix}_filterStats.txt
       cp ${filtered_stats2} ${qadir}/${prefix}_filterStats2.txt
       /scripts/rqcstats.py ${filtered_stats} > stats.json

       /scripts/generate_objects.py --type "nmdc:ReadQCAnalysisActivity" --id ${informed_by} \
             --name "Read QC Activity for ${proj}" --part ${proj} \
             --start ${start} --end $end \
             --resource '${resource}' --url ${url_root}${proj}/qa/ --giturl ${git_url} \
             --extra stats.json \
             --inputs ${read} \
             --outputs \
             ${qadir}/${prefix}_filtered.fastq.gz 'Filtered Reads' \
             ${qadir}/${prefix}_filterStats.txt 'Filtered Stats'
       cp activity.json data_objects.json ${qadir}/

       # Generate assembly objects
       mkdir -p ${assemdir}
       cat ${fasta} | sed ${sed} > ${assemdir}/${prefix}_contigs.fna
       cat ${scaffold} | sed ${sed} > ${assemdir}/${prefix}_scaffolds.fna
       cat ${covstats} | sed ${sed} > ${assemdir}/${prefix}_covstats.txt
       cp ${asmstats} ${assemdir}/${prefix}_stats.json
       cat ${agp} | sed ${sed} > ${assemdir}/${prefix}_assembly.agp
       # TODO: Fix up IDs
       ## Bam file     
       samtools view -h ${bam} | sed ${sed} | \
          samtools view -hb -o ${assemdir}/${prefix}_pairedMapped_sorted.bam
       ## Sam.gz file
       samtools view -h ${samgz} | sed ${sed} | \
          gzip -c - > ${assemdir}/${prefix}_pairedMapped.sam.gz

       /scripts/generate_objects.py --type "nmdc:MetagenomeAssembly" --id ${informed_by} \
             --name "Assembly Activity for ${proj}" --part ${proj} \
             --start ${start} --end $end \
             --resource '${resource}' --url ${url_root}${proj}/assembly/ --giturl ${git_url} \
             --extra ${asmstats} \
             --inputs ${qadir}/${prefix}_filtered.fastq.gz \
             --outputs \
             ${assemdir}/${prefix}_contigs.fna 'Assembled contigs fasta' \
             ${assemdir}/${prefix}_scaffolds.fna 'Assembled scaffold fasta' \
             ${assemdir}/${prefix}_covstats.txt 'Metagenome Contig Coverage Stats' \
             ${assemdir}/${prefix}_assembly.agp 'Assembled AGP file' \
             ${assemdir}/${prefix}_pairedMapped_sorted.bam 'Metagenome Alignment BAM file'
       cp activity.json data_objects.json ${assemdir}/

       # Generate annotation objects
       mkdir -p ${annodir}
       cat ${proteins_faa} | sed ${sed} > ${annodir}/${prefix}_proteins.faa
       cat ${structural_gff} | sed ${sed} > ${annodir}/${prefix}_structural_annotation.gff
       cat ${functional_gff} | sed ${sed} > ${annodir}/${prefix}_functional_annotation.gff
       cat ${ko_tsv} | sed ${sed} > ${annodir}/${prefix}_ko.tsv
       cat ${ec_tsv} | sed ${sed} > ${annodir}/${prefix}_ec.tsv
       cat ${cog_gff} | sed ${sed} > ${annodir}/${prefix}_cog.gff
       cat ${pfam_gff} | sed ${sed} > ${annodir}/${prefix}_pfam.gff
       cat ${tigrfam_gff} | sed ${sed} > ${annodir}/${prefix}_tigrfam.gff
       cat ${smart_gff} | sed ${sed} > ${annodir}/${prefix}_smart.gff
       cat ${supfam_gff} | sed ${sed} > ${annodir}/${prefix}_supfam.gff
       cat ${cath_funfam_gff} | sed ${sed} > ${annodir}/${prefix}_cath_funfam.gff
       cat ${ko_ec_gff} | sed ${sed} > ${annodir}/${prefix}_ko_ec.gff
       cat ${stats_tsv} | sed ${sed} > ${annodir}/${prefix}_stats.tsv
       cat ${stats_json} | sed ${sed} > ${annodir}/${prefix}_stats.json
       nmdc gff2json ${annodir}/${prefix}_functional_annotation.gff -of features.json -oa annotations.json -ai ${informed_by}

       /scripts/generate_objects.py --type "nmdc:MetagenomeAnnotationActivity" --id ${informed_by} \
             --name "Annotation Activity for ${proj}" --part ${proj} \
             --start ${start} --end $end \
             --resource '${resource}' --url ${url_root}${proj}/annotation/ --giturl ${git_url} \
             --inputs ${assemdir}/${prefix}_contigs.fna \
             --outputs \
             ${annodir}/${prefix}_proteins.faa 'Protein FAA' \
             ${annodir}/${prefix}_structural_annotation.gff 'Structural annotation GFF file' \
             ${annodir}/${prefix}_functional_annotation.gff 'Functional annotation GFF file' \
             ${annodir}/${prefix}_ko.tsv 'KO TSV file' \
             ${annodir}/${prefix}_ec.tsv 'EC TSV file' \
             ${annodir}/${prefix}_cog.gff 'COG GFF file' \
             ${annodir}/${prefix}_pfam.gff 'PFAM GFF file' \
             ${annodir}/${prefix}_tigrfam.gff 'TigrFam GFF file' \
             ${annodir}/${prefix}_smart.gff 'SMART GFF file' \
             ${annodir}/${prefix}_supfam.gff 'SuperFam GFF file' \
             ${annodir}/${prefix}_cath_funfam.gff 'Cath FunFam GFF file' \
             ${annodir}/${prefix}_ko_ec.gff 'KO_EC GFF file'
       cp features.json annotations.json activity.json data_objects.json ${annodir}/

       # MAGS
       mkdir -p ${magsdir}
       cat ${lowdepth} | sed ${sed} > ${magsdir}/${prefix}_bins.lowDepth.fa
       cat ${short} | sed ${sed} > ${magsdir}/${prefix}_bins.tooShort.fa
       cat ${unbinned} | sed ${sed} > ${magsdir}/${prefix}_bins.unbinned.fa
       cp ${checkm} ${magsdir}/${prefix}_checkm_qa.out
       #TODO: Check naming
       mkdir -p hqmq
       if [ ${n_hqmq} -gt 0 ] ; then
           (cd hqmq && cp ${sep=" " hqmq_bin_fasta_files} .)
           (cd hqmq && sed -i ${sed} *.fa && zip ../${prefix}_hqmq_bin.zip *.fa)
       else
           (cd hqmq && touch no_hqmq_mags.txt)
           (cd hqmq && zip ../${prefix}_hqmq_bin.zip *.txt)
       fi
       cp ${prefix}_hqmq_bin.zip ${magsdir}

       mkdir meta
       (cd meta && cp ${sep=" " bin_fasta_files} .)
       (cd meta && sed -i ${sed} *.fa && zip ../${prefix}_metabat_bin.zip *.fa)
       cp ${prefix}_metabat_bin.zip ${magsdir}

       /scripts/generate_objects.py --type "nmdc:MAGsAnalysisActivity" --id ${informed_by} \
             --name "MAGs Analysis Activity for ${proj}" --part ${proj} \
             --start ${start} --end $end \
             --resource '${resource}' --url ${url_root}${proj}/MAGs/ --giturl ${git_url} \
             --inputs ${assemdir}/${prefix}_contigs.fna \
                      ${assemdir}/${prefix}_pairedMapped_sorted.bam \
                      ${annodir}/${prefix}_functional_annotation.gff \
             --outputs \
             ${magsdir}/${prefix}_bins.lowDepth.fa "lowDepth (mean cov <1 )  filtered contigs fasta file by metabat2" \
             ${magsdir}/${prefix}_bins.tooShort.fa "tooShort (< 3kb) filtered contigs fasta file by metaBat2" \
             ${magsdir}/${prefix}_bins.unbinned.fa "unbinned fasta file from metabat2" \
             ${magsdir}/${prefix}_checkm_qa.out "metabat2 bin checkm quality assessment result" \
             ${magsdir}/${prefix}_hqmq_bin.zip "high-quality and medium-quality bins" \
             ${magsdir}${prefix}_metabat_bin.zip "metabat2 bins"
#             "gtdbtk.bac120.summary.tsv" "gtdbtk bacterial assignment result summary table" \
#             "gtdbtk.ar122.summary.tsv" "gtdbtk archaea assignment result summary table" \
       cp activity.json data_objects.json ${magsdir}/

       #Readbased Analysis
       mkdir -p ${rbadir}
       cp ${gottcha2_report_tsv} ${rbadir}/${prefix}_gottcha2_report.tsv
       cp ${gottcha2_full_tsv} ${rbadir}/${prefix}_gottcha2_report_full.tsv
       cp ${gottcha2_krona_html} ${rbadir}/${prefix}_gottcha2_krona.html

       cp ${centrifuge_classification_tsv} ${rbadir}/${prefix}_centrifuge_classification.tsv
       cp ${centrifuge_report_tsv} ${rbadir}/${prefix}_centrifuge_report.tsv
       cp ${centrifuge_krona_html} ${rbadir}/${prefix}_centrifuge_krona.html
          
       cp ${kraken2_classification_tsv} ${rbadir}/${prefix}_kraken2_classification.tsv
       cp ${kraken2_report_tsv} ${rbadir}/${prefix}_kraken2_report.tsv
       cp ${kraken2_krona_html} ${rbadir}/${prefix}_kraken2_krona.html
 
       /scripts/generate_objects.py --type "nmdc:ReadbasedAnalysisActivity" --id ${informed_by} \
             --name "Readbased Analysis Activity for ${proj}" --part ${proj} \
             --start ${start} --end $end \
             --resource '${resource}' --url ${url_root}${proj}/ReadbasedAnalysis/ --giturl ${git_url} \
             --inputs ${qadir}/${prefix}_filtered.fastq.gz \
             --outputs \
             ${rbadir}/${prefix}_gottcha2_report.tsv "Gottcha2 TSV report" \
             ${rbadir}/${prefix}_gottcha2_report_full.tsv "Gottcha2 full TSV report" \
             ${rbadir}/${prefix}_gottcha2_krona.html "Gottcha2 Krona HTML report" \
             ${rbadir}/${prefix}_centrifuge_classification.tsv "Centrifuge classification TSV report" \
             ${rbadir}/${prefix}_centrifuge_report.tsv "Centrifuge TSV report" \
             ${rbadir}/${prefix}_centrifuge_krona.html "Centrifuge Krona HTML report" \
             ${rbadir}/${prefix}_kraken2_classification.tsv "Kraken classification TSV report" \
             ${rbadir}/${prefix}_kraken2_report.tsv "Kraken2 TSV report" \
             ${rbadir}/${prefix}_kraken2_krona.html "Kraken2 Krona HTML report"
       cp activity.json data_objects.json ${rbadir}/

       # Top-level Container
       /scripts/generate_objects.py --type "nmdc:MetagenomeAnalysisActivity" --id ${informed_by} \
             --name "Metagenome Analysis Activity for ${proj}" \
             --activityid=${proj} \
             --start ${start} --end $end \
             --resource '${resource}' --url ${url_root}${proj}/ --giturl ${git_url}
       cp activity.json ${outdir}/
   }

   runtime {
     memory: "10 GiB"
     cpu:  4
     maxRetries: 1
     docker: container
   }
}


