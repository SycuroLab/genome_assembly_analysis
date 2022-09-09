# ***************************************
# * Snakefile for genome_assembly_analysis pipeline *
# ***************************************

# **** Variables ****

configfile: "config.yaml"

# **** Imports ****

import pandas as pd
import os

os.environ["GTDBTK_DATA_PATH"] = "/bulk/IMCshared_bulk/shared/dbs/gtdbtk-1.5.0/db"

SAMPLES = pd.read_csv(config["list_files"], header = None)
SAMPLES = SAMPLES[0].tolist()

# **** Rules ****

rule all:
    input:
	#expand(config["input_dir"]+"/assembly_analysis/{sample}.fa", sample=SAMPLES), 
        expand(config["output_dir"]+"/assembly_analysis/{sample}/quast/transposed_report.tsv", sample=SAMPLES),
        expand(config["output_dir"]+"/assembly_analysis/{sample}/prokka/{sample}.fna",sample=SAMPLES),
        expand(config["output_dir"]+"/assembly_analysis/{sample}/prokka/{sample}.gff",sample=SAMPLES),
        expand(config["output_dir"]+"/assembly_analysis/{sample}/extracted_sequences/cpn60_metadata.csv",sample=SAMPLES),
        expand(config["output_dir"]+"/assembly_analysis/{sample}/metaerg/data/all.gff",sample=SAMPLES),
        expand(config["output_dir"]+"/assembly_analysis/{sample}/checkm/checkm.tsv",sample=SAMPLES),
#        expand(config["output_dir"]+"/assembly_analysis/{sample}/gtdbtk/gtdbtk.bac120.summary.tsv",sample=SAMPLES)


rule quast:
    input:
        assembly_file = os.path.join(config["input_dir"], "{sample}.fa")
    output:
        quast_transposed_report_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","quast","transposed_report.tsv")
    params:
        quast_dir = os.path.join(config["output_dir"],"assembly_analysis","{sample}","quast"),
        threads = config["quast_threads"]
    conda: "utils/envs/quast_env.yaml"
    shell:
       "quast.py --output-dir {params.quast_dir} --threads {params.threads} {input.assembly_file}"

rule prokka:
    input:
        assembly_file = os.path.join(config["input_dir"],"{sample}.fa")
    output:
        prokka_fna_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","prokka","{sample}.fna"),
        prokka_gff_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","prokka","{sample}.gff")
    params:
        prokka_dir = os.path.join(config["output_dir"],"assembly_analysis","{sample}","prokka"),
        threads = config["prokka_threads"],
	prefix = "{sample}"
    conda: "utils/envs/prokka_env.yaml"
    shell:
       "prokka --metagenome --outdir {params.prokka_dir} --prefix {params.prefix} {input.assembly_file} --cpus {params.threads} --rfam 1 --force"

rule extract_marker_sequences:
    input:
        prokka_fna_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","prokka","{sample}.fna"),
        prokka_gff_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","prokka","{sample}.gff")
    output:
        extracted_marker_seqs_csv_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","extracted_sequences","cpn60_metadata.csv"),    
    params:
        extracted_sequences_dir = os.path.join(config["output_dir"],"assembly_analysis","{sample}","extracted_sequences"),
    conda: "utils/envs/biopython_env.yaml"
    shell:
        "python utils/scripts/extract_marker_sequences.py --fasta_infile {input.prokka_fna_file} --gff_infile {input.prokka_gff_file} --output_dir {params.extracted_sequences_dir}"

rule metaerg:
    input:
        assembly_file = os.path.join(config["input_dir"],"{sample}.fa")
    output:
#        metaerg_fna_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","assembly","metaerg","{sample}_metagenome.fna"),
        metaerg_gff_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","metaerg","data","all.gff")
    params:
        metaerg_dir = os.path.join(config["output_dir"],"assembly_analysis","{sample}","metaerg"),
        metaerg_database_path = config["metaerg_database_path"],
        locustag = "{sample}",
        threads = config["metaerg_threads"]
    shell:
       "singularity run -H $HOME -B {params.metaerg_database_path}:/NGStools/metaerg/db -B /work:/work -B /bulk:/bulk /global/software/singularity/images/software/metaerg2.sif /NGStools/metaerg/bin/metaerg.pl --mincontiglen 200 --gcode 11 --gtype meta --minorflen 180 --cpus {params.threads} --evalue 1e-05 --identity 20 --coverage 70 --locustag {params.locustag} --force --outdir {params.metaerg_dir} {input.assembly_file}"

rule checkm:
    input:
        assembly_file = os.path.join(config["input_dir"],"{sample}.fa")
    output:
        checkm_table_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","checkm","checkm.tsv")
    params:
        checkm_database = config["checkm_database_path"],
        checkm_dir = os.path.join(config["output_dir"],"assembly_analysis","{sample}","checkm"),
	threads = config["checkm_threads"]
    conda: "utils/envs/checkm_env.yaml"
    shell:
       "checkm data setRoot {params.checkm_database}; "
       "mkdir -p {params.checkm_dir}; "
       "filename=$(basename {input.assembly_file}); "
       "cp {input.assembly_file} {params.checkm_dir}/$filename; "
       "checkm lineage_wf -t {params.threads} -x fa --tab_table --file {output.checkm_table_file} {params.checkm_dir} {params.checkm_dir}; "

#rule gtdbtk:
#    input:
#        assembly_file = os.path.join(config["input_dir"],"{sample}.fa")
#    output:
#        gtdbtk_file = os.path.join(config["output_dir"],"assembly_analysis","{sample}","gtdbtk","gtdbtk.bac120.summary.tsv")
#    params:
#       gtdbtk_data_path = config["gtdbtk_database_path"],
#       gtdbtk_dir = os.path.join(config["output_dir"],"assembly_analysis","{sample}","gtdbtk"),
#       threads = config["gtdbtk_threads"]
#    conda: "utils/envs/gtdbtk_env.yaml"
#    shell:
#       "GTDBTK_DATA_PATH=\"{params.gtdbtk_data_path}\"; "       
#       "filename=$(basename {input.assembly_file}); "
#       "cp {input.assembly_file} {params.gtdbtk_dir}/$filename; "
#       "gtdbtk classify_wf --genome_dir {params.gtdbtk_dir} --extension \"fa\" --cpus {params.threads} --out_dir {params.gtdbtk_dir}; "

