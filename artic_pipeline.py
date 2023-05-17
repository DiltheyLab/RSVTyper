### This program executes the artic pipeline with all necessary variables.
import os
variable_list = []
with open("variables_for_artic.txt", "r") as fin:
    for line in fin:
        line = line.rstrip()
        variable = line.split()[-1]
        variable_list.append(variable)


path_to_reads = variable_list[0]
if path_to_reads[-1] == "/":
    path_to_reads = path_to_reads[:-1]
barcode_list = path_to_reads.split("/")
for element in barcode_list:
    if "barcode" in element:
        barcode = element
medaka_model = variable_list[1]
path_to_primer_scheme = variable_list[2]
version = variable_list[3]
sample = variable_list[4]
output_dir = variable_list[5]

os.system(f"mv variables_for_artic.txt {output_dir}")
os.chdir(output_dir)

# Read filtering
os.system(f"artic guppyplex --skip-quality-check --min-length 350 --max-length 900 --directory {path_to_reads} --prefix demultiplexed")

# Medaka
os.system(f"artic minion --medaka --medaka-model {medaka_model} --normalise 200 --threads 4 --scheme-directory {path_to_primer_scheme[:-1]} --read-file demultiplexed_{barcode}.fastq RSV-2023/{version} {sample}")
