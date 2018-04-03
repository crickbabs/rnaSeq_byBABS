#!/bin/sh

# output directories
html_dir=html
dot_dir=dot
txt_dir=txt
if [[ ! -d $html_dir ]]; then mkdir -v $html_dir; fi
if [[ ! -d $dot_dir ]]; then mkdir -v $dot_dir; fi
if [[ ! -d $txt_dir ]]; then mkdir -v $txt_dir; fi

# fastq screen configuration file
FSCREEN=conf/fastq_screen.conf

# multiqc configuration file (not mandatory)
MULTIQC_FILE=yml/multiqc_config.yml

nextflow run rnaSeq-byBABS.nf \
	--design "csv/single_end.csv" \
	--output_directory "nextflow_single_end" \
	--binomial_nomenclature "homo sapiens" \
	--strandedness "none" \
	--read_length "77" \
	--single_end \
	--modules "yml/modules.yml" \
	--genome_type "ensembl" \
	--genome_version "GRCh38" \
	--genome_release "86" \
	--genome_sequence_extension "primary_assembly.fa" \
	--fastq_screen_conf $FSCREEN \
	--publish_directory_mode "symlink" \
	--publish_directory_overwrite \
	--r_script "r/analysis_example.r" \
	--multiqc_config $MULTIQC_FILE \
	-with-timeline $html_dir/timeline_paired_end.html \
	-with-dag $dot_dir/dag_paired_end.dot \
	-with-trace $txt_dir/trace_paired_end.txt \
	-resume

