#!/usr/bin/python
import os
import sys
import re
import csv
import argparse

parser = argparse.ArgumentParser()

# Sample command
# python merge_assembly_analysis_data.py --input_dir /work/sycuro_lab/kevin/Pasolli_Vaginal_MAGs_analysis_2022-04-26/assembly_analysis --output_dir /work/sycuro_lab/kevin/Pasolli_Vaginal_MAGs_analysis_2022-04-26


# /work/sycuro_lab/kevin/Pasolli_Vaginal_MAGs_analysis_2022-04-26/assembly_analysis
input_dir = None
output_dir = None

parser.add_argument('--input_dir', action='store', dest='input_dir',
                    help='input directory as input. (i.e. $HOME)')
parser.add_argument('--output_dir', action='store', dest='output_dir',
                    help='output directory as input. (i.e. $HOME)')
parser.add_argument('--version', action='version', version='%(prog)s 1.0')

results = parser.parse_args()

input_dir = results.input_dir
output_dir = results.output_dir

if(input_dir == None):
    print('\n')
    print('error: please use the --input_dir option to specify the input directory as input')
    print('input_dir =' + ' ' + str(input_dir))
    print('\n')
    parser.print_help()
    sys.exit(1)
if(output_dir == None):
    print('\n')
    print('error: please use the --output_dir option to specify the output directory as input')
    print('output_dir =' + ' ' + str(output_dir))
    print('\n')
    parser.print_help()
    sys.exit(1)

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def parse_quast_data(quast_infile):

    quast_data = {}
    i = 0
    with open(quast_infile, "r") as quast_input_file:
        csv_reader = csv.reader(quast_input_file, delimiter='\t')
        for row in csv_reader:

            #print(row)
            #sys.exit()
            if(i != 0):
                #print(row)
                #sys.exit()
                (assembly_name, num_contigs_le_0, num_contigs_le_1000, num_contigs_le_5000, num_contigs_le_10000, num_contigs_le_25000, num_contigs_le_50000, total_length_le_0, total_length_le_1000, total_length_le_5000, total_length_le_10000, total_length_le_25000, total_length_le_50000, num_contigs, largest_contig, total_length, percent_gc, n50, n75, l50, l75, num_ns_per_100_kbp) = row
                quast_data = [assembly_name, num_contigs_le_0, num_contigs_le_1000, num_contigs_le_5000, num_contigs_le_10000, num_contigs_le_25000, num_contigs_le_50000, total_length_le_0, total_length_le_1000, total_length_le_5000, total_length_le_10000, total_length_le_25000, total_length_le_50000, num_contigs, largest_contig, total_length, percent_gc, n50, n75, l50, l75, num_ns_per_100_kbp]
            i += 1
    return quast_data

def parse_checkm_data(checkm_infile):

    checkm_data = {}
    i = 0
    with open(checkm_infile, "r") as checkm_input_file:
        csv_reader = csv.reader(checkm_input_file, delimiter='\t')
        for row in csv_reader:

            #print(row)
            #sys.exit()
            if(i != 0):
                #print(row)
                #sys.exit()
                (bin_id,marker_lineage,num_genomes,num_markers,num_marker_sets,zero,one,two,three,four,five_plus,completeness,contamination,strain_heterogeneity) = row
                checkm_data = [bin_id,marker_lineage,num_genomes,num_markers,num_marker_sets,zero,one,two,three,four,five_plus,completeness,contamination,strain_heterogeneity]
            i += 1
    return checkm_data

def parse_gtdbtk_data(gtdbtk_infile):

    gtdbtk_data = {}
    i = 0
    with open(gtdbtk_infile, "r") as gtdbtk_input_file:
        csv_reader = csv.reader(gtdbtk_input_file, delimiter='\t')
        for row in csv_reader:

            #print(row)
            #sys.exit()
            if(i != 0):
                #print(row)
                #sys.exit()
                (user_genome,classification,fastani_reference,fastani_reference_radius,fastani_taxonomy,fastani_ani,fastani_af,closest_placement_reference,closest_placement_radius,closest_placement_taxonomy,closest_placement_ani,closest_placement_af,pplacer_taxonomy,classification_method,note,other_related_references,msa_percent,translation_table,red_value,warnings) = row
                gtdbtk_data = [user_genome,classification,fastani_reference,fastani_reference_radius,fastani_taxonomy,fastani_ani,fastani_af,closest_placement_reference,closest_placement_radius,closest_placement_taxonomy,closest_placement_ani,closest_placement_af,pplacer_taxonomy,classification_method,note,other_related_references,msa_percent,translation_table,red_value,warnings]
            i += 1
    return gtdbtk_data

dirs = [ f.path for f in os.scandir(input_dir) if f.is_dir() ]

quast_tsv_outfile = os.path.join(output_dir, "all_quast_transposed_report.tsv")
quast_tsv_output_file = open(quast_tsv_outfile, 'w+')
quast_tsv_writer = csv.writer(quast_tsv_output_file, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
quast_tsv_writer.writerow(["Assembly","# contigs (>= 0 bp)","# contigs (>= 1000 bp)","# contigs (>= 5000 bp)","# contigs (>= 10000 bp)","# contigs (>= 25000 bp)","# contigs (>= 50000 bp)","Total length (>= 0 bp)","Total length (>= 1000 bp)","Total length (>= 5000 bp)","Total length (>= 10000 bp)","Total length (>= 25000 bp)","Total length (>= 50000 bp)","# contigs","Largest contig","Total length","GC (%)","N50","N75","L50","L75","# N's per 100 kbp"])

checkm_tsv_outfile = os.path.join(output_dir, "all_checkm.tsv")
checkm_tsv_output_file = open(checkm_tsv_outfile, 'w+')
checkm_tsv_writer = csv.writer(checkm_tsv_output_file, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
checkm_tsv_writer.writerow(["Bin Id","Marker lineage","# genomes","# markers","# marker sets","0","1","2","3","4","5+","Completeness","Contamination","Strain heterogeneity"])

gtdbtk_tsv_outfile = os.path.join(output_dir, "all_gtdbtk.bac120.summary.tsv")
gtdbtk_tsv_output_file = open(gtdbtk_tsv_outfile, 'w+')
gtdbtk_tsv_writer = csv.writer(gtdbtk_tsv_output_file, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
gtdbtk_tsv_writer.writerow(["user_genome","classification","fastani_reference","fastani_reference_radius","fastani_taxonomy","fastani_ani","fastani_af","closest_placement_reference","closest_placement_radius","closest_placement_taxonomy","closest_placement_ani","closest_placement_af","pplacer_taxonomy","classification_method","note","other_related_references(genome_id,species_name,radius,ANI,AF)","msa_percent","translation_table","red_value","warnings"])

merged_assembly_analysis_tsv_outfile = os.path.join(output_dir, "all_merged_assembly_analysis_metadata.tsv")
merged_assembly_analysis_tsv_output_file = open(merged_assembly_analysis_tsv_outfile, 'w+')
merged_assembly_analysis_tsv_writer = csv.writer(merged_assembly_analysis_tsv_output_file, delimiter='\t', quotechar='', quoting=csv.QUOTE_NONE)
merged_assembly_analysis_tsv_writer.writerow(["Genome ID", "GTDB-tk Classification","Total Number of Contigs","Largest Contig","Total Length","GC (%)","N50","Completeness","Contamination","Strain heterogeneity","GTDB-tk Lineage","FastANI Reference"," FastANI Reference Radius","FastANI Taxonomy","FastANI ANI","FastANI AF","Closest Placement Reference","Closest Placement Radius","Closest Placement Taxonomy","Closest Placement ANI","Closest Placement AF","pplacer Taxonomy","Classification Method","Note","Other Related References ( Genome ID, Species Name, Radius, ANI, AF )","MSA Percent","Translation Table","Red Value","Warnings"])

for assembly_dir in dirs:
    genome_name = os.path.basename(assembly_dir)
    print(genome_name)
    print(os.path.join(output_dir,genome_name))
    print(assembly_dir)
    #print(assembly_dir)
    #sys.exit()
    quast_infile = os.path.join(assembly_dir,"quast/transposed_report.tsv")
    print(quast_infile)
    quast_data = parse_quast_data(quast_infile)
    quast_tsv_writer.writerow(quast_data)

    checkm_infile = os.path.join(assembly_dir,"checkm/checkm.tsv")
    print(checkm_infile)
    checkm_data = parse_checkm_data(checkm_infile)
    checkm_tsv_writer.writerow(checkm_data)

    gtdbtk_infile = os.path.join(assembly_dir,"gtdbtk/gtdbtk.bac120.summary.tsv")
    print(gtdbtk_infile)
    gtdbtk_data = parse_gtdbtk_data(gtdbtk_infile)
    gtdbtk_tsv_writer.writerow(gtdbtk_data)
    
    (assembly_name, num_contigs_le_0, num_contigs_le_1000, num_contigs_le_5000, num_contigs_le_10000, num_contigs_le_25000, num_contigs_le_50000, total_length_le_0, total_length_le_1000, total_length_le_5000, total_length_le_10000, total_length_le_25000, total_length_le_50000, num_contigs, largest_contig, total_length, percent_gc, n50, n75, l50, l75, num_ns_per_100_kbp) = quast_data
    (bin_id,marker_lineage,num_genomes,num_markers,num_marker_sets,zero,one,two,three,four,five_plus,completeness,contamination,strain_heterogeneity) = checkm_data
    (user_genome,classification,fastani_reference,fastani_reference_radius,fastani_taxonomy,fastani_ani,fastani_af,closest_placement_reference,closest_placement_radius,closest_placement_taxonomy,closest_placement_ani,closest_placement_af,pplacer_taxonomy,classification_method,note,other_related_references,msa_percent,translation_table,red_value,warnings) = gtdbtk_data

    print(classification)
    #d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus jensenii
    
    taxonomy = classification.split(";")[-1]
    print(taxonomy)

    #sys.exit()

    merged_assembly_analysis_tsv_writer.writerow([genome_name,taxonomy,num_contigs,largest_contig,total_length,percent_gc,n50,completeness,contamination,strain_heterogeneity,classification,fastani_reference,fastani_reference_radius,fastani_taxonomy,fastani_ani,fastani_af,closest_placement_reference,closest_placement_radius,closest_placement_taxonomy,closest_placement_ani,closest_placement_af,pplacer_taxonomy,classification_method,note,other_related_references,msa_percent,translation_table,red_value,warnings])
