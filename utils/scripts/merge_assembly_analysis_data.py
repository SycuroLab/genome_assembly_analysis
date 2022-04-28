#!/usr/bin/python
import os
import sys
import re
import csv
import argparse

parser = argparse.ArgumentParser()

# Sample command
#


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



