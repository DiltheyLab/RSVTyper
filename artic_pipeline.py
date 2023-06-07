### This program executes the artic pipeline with all necessary variables.
import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", help = "Path to basecalled, demultiplexed fastq-files. It should end with the barcode directory (e.g. barcode15/)", required = True)
parser.add_argument("-s", "--sample", help = "Name of the sample", required = True)
parser.add_argument("-o", "--outputDir", help = "Output directory", required = True)
parser.add_argument("-m", "--medakaModel", help = "Medaka model that should be used for the artic pipeline (depends on basecaller used)", required = True)
parser.add_argument("-a", "--schemeDir", help = "Path to primer scheme if the location of it was changed", default = path_to_primer_scheme)
parser.add_argument("-v", "--version", required = True)

args = parser.parse_args()

path_to_reads = args.input
if path_to_reads[-1] == "/":
    path_to_reads = path_to_reads[:-1]

output_dir = args.outputDir
medaka_model = args.medakaModel
path_to_primer_scheme = args.schemeDir
barcode_list = path_to_reads.split("/")
for element in barcode_list:
    if "barcode" in element:
        barcode = element
version = args.version
sample = args.sample

# Read filtering
os.system(f"artic guppyplex --skip-quality-check --min-length 350 --max-length 900 --directory {path_to_reads} --prefix demultiplexed")

# Medaka
os.system(f"artic minion --medaka --medaka-model {medaka_model} --normalise 200 --threads 4 --scheme-directory {path_to_primer_scheme[:-1]} --read-file demultiplexed_{barcode}.fastq RSV-2023/{version} {sample}")
