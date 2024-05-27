#!/usr/bin/env python
"""This script runs the artic pipeline with Nanopore reads that have already been basecalled and demultiplexed. Variant calling is done by the experimental medaka pipeline.
It is specificially made for RSV samples as it determines the subtype and chooses the reference that matches the given reads the most."""

import os
import sys
import argparse
import re

def detect_subtype():

    reference_subtype = "all_consensus_references.fasta"
    coverage_summary = f"{sample}_coverage_summary.txt"
    alignment = f"{sample}_vs_{reference_subtype[:-6]}"
    minimap_command = f"minimap2 -ax map-ont {path_to_reference}/{reference_subtype} {path_to_reads}/*fastq* > {path_to_subtype_detection}/{alignment}.sam"
    os.system(minimap_command)

    convert_to_bam = f"samtools view -1 -S -b -F 256 {path_to_subtype_detection}/{alignment}.sam > {path_to_subtype_detection}/{alignment}.bam"
    sort = f"samtools sort {path_to_subtype_detection}/{alignment}.bam -o {path_to_subtype_detection}/{alignment}.sorted.bam"
    index = f"samtools index {path_to_subtype_detection}/{alignment}.sorted.bam"
    command_chain = f"{convert_to_bam} && {sort} && {index}"
    os.system(command_chain)

    # Determining the number of mapped reads for every reference in the region 3500-5500.
    # This region has been selected as it seems to have diverged between the subtypes.
    with open(f"{path_to_subtype_detection}/{coverage_summary}", "w") as fout:
        for i in range(1, 11):
            consensus_seq = "cluster_A_" + str(i) + "_consensus"
            samtools_coverage = f"samtools coverage -r {consensus_seq}:3500-5500 {path_to_subtype_detection}/{alignment}.sorted.bam"
            creating_coverage_summary = f"({samtools_coverage}) >> {path_to_subtype_detection}/{coverage_summary}"
            os.system(creating_coverage_summary)
        for i in range(1, 11):
            consensus_seq = "cluster_B_" + str(i) + "_consensus"
            samtools_coverage = f"samtools coverage -r {consensus_seq}:3500-5500 {path_to_subtype_detection}/{alignment}.sorted.bam"
            creating_coverage_summary = f"({samtools_coverage}) >> {path_to_subtype_detection}/{coverage_summary}"
            os.system(creating_coverage_summary)

    # For each subtype the reference with the most reads aligned in that region is chosen.
    max_num_reads_A = 0
    max_num_reads_B = 0
    with open(f"{path_to_subtype_detection}/{coverage_summary}", "r") as fin:
        for line in fin:
            line = line.rstrip()
            if line[0] != "#":
                info = line.split("\t")
                subtype = info[0][8]
                num_reads = int(info[3])
                if subtype == "A":
                    if num_reads > max_num_reads_A:
                        max_num_reads_A = num_reads
                if subtype == "B":
                    if num_reads > max_num_reads_B:
                        max_num_reads_B = num_reads

    # Calculation of the ratios between the number of reads of the chosen subtype specific references.
    # If there is no reference with any reads mapped across all references of one subtype the value is set to 1 
    # to prevent an error that occurs when a number is divided by 0.
    if max_num_reads_A == 0:
        max_num_reads_A = 1
    if max_num_reads_B == 0:
        max_num_reads_B = 1
    A_B_ratio = max_num_reads_A / max_num_reads_B
    B_A_ratio = max_num_reads_B / max_num_reads_A

    # If one of the ratios between the maximum number of reads for each subtype is higher than 5 then the subtype detection is considered to be valid.
    if A_B_ratio > 5:
        subtype_A = True
    else:
        subtype_A = False
    if B_A_ratio > 5:
        subtype_B = True
    else:
        subtype_B = False

    if subtype_A == True and subtype_B == False:
        print("Subtype: A")
        final_subtype = "A"
    elif subtype_B == True and subtype_A == False:
        print("Subtype: B")
        final_subtype = "B"
    elif subtype_A == True and subtype_B == True:
        sys.exit("Error: Sample can be matched to subtype A and subtype B. This should not happen.")
    elif subtype_A == False and subtype_B == False:
        sys.exit("Error: Sample cannot be matched to any subtype. Not enough reads?")
    
    with open(output_dir + "/" + final_summary, "w") as fout:
        fout.write("Subtype: " + final_subtype + "\n")
    
    return(final_subtype)

def detect_duplication(subtype):

    # Read filtering
    os.system(f"artic guppyplex --skip-quality-check --min-length {amplicon_length[0]} --max-length {amplicon_length[1]} --directory {path_to_reads} --prefix demultiplexed")
    
    dup_alignment = f"{sample}_vs_subtype_{subtype}_with_duplications"
    minimap_command = f"minimap2 -ax map-ont -c {path_to_reference}/cluster_{subtype}_seqs_with_duplication.fasta {artic_dir}/demultiplexed_{barcode}.fastq > {duplication_dir}/{dup_alignment}.sam"
    convert_to_bam = f"samtools view -1 -S -b -F 256 {duplication_dir}/{dup_alignment}.sam > {duplication_dir}/{dup_alignment}.bam"
    samtools_sort = f"samtools sort {duplication_dir}/{dup_alignment}.bam -o {duplication_dir}/{dup_alignment}.sorted.bam"
    samtools_index = f"samtools index {duplication_dir}/{dup_alignment}.sorted.bam"
    os.system(f"{minimap_command} && {convert_to_bam} && {samtools_sort} && {samtools_index}")

    clusters_duplication_file = f"cluster_{subtype}_with_duplication.txt"
    dup_dict = {}
    with open(f"{path_to_reference}/{clusters_duplication_file}", "r") as fin3:
        for line3 in fin3:
            line_list = line3.rstrip().split()
            cluster_seq = line_list[0]
            dup_begin = int(line_list[1])
            dup_end = int(line_list[2])
            dup_dict[cluster_seq] = [dup_begin, dup_end]
            sam_duplication_region = f"{sample}_vs_duplication_region.sam"
            samtools_view_command = f"samtools view {duplication_dir}/{dup_alignment}.sorted.bam {cluster_seq}:{str(dup_begin - 50)}-{str(dup_end + 50)} >> {duplication_dir}/{sam_duplication_region}"
            os.system(samtools_view_command)

    duplication_summary = f"{sample}_deletion_ratios.txt"
    long_deletions = 0
    no_reads = 0
    with open(f"{duplication_dir}/{duplication_summary}", "w") as fout:
        fout.write(f"{sample} \n")
        with open(f"{duplication_dir}/{sam_duplication_region}", "r") as fin4:
                for line4 in fin4:
                    line4 = line4.rstrip().split()
                    alignment_type = line4[15]
                    if "P" in alignment_type:
                        aln_cluster_seq = line4[2]
                        cigar = line4[5]
                        end_of_cigar = len(cigar)
                        begin_of_cigar_part = 0
                        cigar_len = 0
                        aln_begin = int(line4[3])
                        new_long_deletions = 0
                        for i in range(0, end_of_cigar):
                            current_character = cigar[i]
                            if current_character.isdigit() == False:
                                cigar_part = cigar[begin_of_cigar_part:i + 1]
                                begin_of_cigar_part = i + 1
                                cigar_len += int(cigar_part[:-1])
                                if "D" in cigar_part:
                                    deletion_length = int(cigar_part[:-1])
                                    if subtype == "A":
                                        if deletion_length >= 60:
                                            new_long_deletions += 1
                                    elif subtype == "B":
                                        if deletion_length >= 50:
                                            new_long_deletions += 1
                        aln_end = aln_begin + cigar_len
                        if aln_begin <= dup_dict[aln_cluster_seq][0] and aln_end >= dup_dict[aln_cluster_seq][1]:
                            no_reads += 1
                            long_deletions += new_long_deletions

        duplication_ratio = long_deletions / no_reads
        fout.write(f"Number of reads: {str(no_reads)} \n")
        fout.write(f"Number of long deletions: {str(long_deletions)} \n")
        fout.write(f"Ratio: {str(duplication_ratio)}")
    if duplication_ratio > 0.7:
        reference_file = f"cluster_{subtype}_seqs_without_duplication.fasta"
    elif duplication_ratio < 0.7:
        reference_file = f"cluster_{subtype}_seqs_with_duplication.fasta"
    return reference_file

def reference_selection(ref_file):

    paf_output = f"{sample}_vs_{ref_file}.paf"
    minimap_command_ref = f"minimap2 -c -x map-ont {path_to_reference}/{ref_file} {path_to_reads}/*fastq* > {reference_selection_directory}/{paf_output}"
    os.system(minimap_command_ref)

    primary_alignment_dict = {}
    with open(f"{reference_selection_directory}/{paf_output}", "r") as fin:
        for line in fin:
            line = line.rstrip()
            line_list = line.split()
            alignment_tag = line_list[16][-1]
            cluster_seq = line_list[5]
            if alignment_tag == "P":
                if cluster_seq not in primary_alignment_dict:
                    primary_alignment_dict[cluster_seq] = 1
                else:
                    primary_alignment_dict[cluster_seq] += 1

    output_ref_select = f"{sample}_primary_alignments.txt"
    sorted_primary_alignment_dict = sorted(primary_alignment_dict.items(), key=lambda x:x[1])
    final_reference = sorted_primary_alignment_dict[-1][0]

    with open(f"{reference_selection_directory}/{output_ref_select}", "w") as fout:
        fout.write(sample + "\n")
        for item in sorted_primary_alignment_dict:
            fout.write(item[0] + "\t" + str(item[1]) + "\n")

    print("Reference: " + final_reference)
    with open(output_dir + "/" + final_summary, "a") as fout:
        fout.write("Reference: " + final_reference + "\n")

    # Convert reference to directory in which the correct primal scheme for that reference lies
    if "A" in final_reference:
        cluster_no = final_reference.split("_")[-2]
        version = "VA" + cluster_no
    if "B" in final_reference:
        cluster_no = final_reference.split("_")[-2]
        version = "VB" + cluster_no
    return version

def artic_minion(version):
    os.system(f"artic minion --no-longshot --medaka --medaka-model {medaka_model} --normalise 200 --threads 4 --scheme-directory {path_to_primer_scheme} --read-file demultiplexed_{barcode}.fastq {scheme_version}/{version} {sample}")

def nextclade(subtype, nextclade_output):
    nextclade_subtype = subtype.lower()
    consensus_seq = f"{sample}.consensus.fasta"
    os.system(f"nextclade run -d rsv_{nextclade_subtype} -O {nextclade_dir} -s={nextclade_output} {artic_dir}/{consensus_seq} --output-basename '{sample}_nextclade'")

    if nextclade_output == "all":
        nextclade_output = "tsv"
    with open(f"{nextclade_dir}/{sample}_nextclade.{nextclade_output}", "r") as fin:
        for line in fin:
            line = line.rstrip()
            if nextclade_output == "tsv" or nextclade_output == "csv":
                if nextclade_output == "tsv":
                    line_list = line.split("\t")
                elif nextclade_output == "csv":
                    line_list = line.split(";")
                if line_list[0] != "index":
                    clade = line_list[2]
                    g_clade = line_list[3]
            if nextclade_output == "json":
                if "\"clade\":" in line:
                    line_list = line.split()
                    clade = line_list[1][1:-2]
                if "\"G_clade\":" in line:
                    line_list = line.split()
                    g_clade = line_list[1][1:-2]

        if nextclade_output == "ndjson":
            clade = ""
            g_clade = ""
            with open(f"{nextclade_dir}/{sample}_nextclade.{nextclade_output}", "r") as fin:
                line_list = re.findall("\{(.*?)\}", line)
            for element in line_list:
                if "\"G_clade\"" in element:
                        new_list = element.split(":")
                        g_clade = new_list[1][1:-2]
                elif "clade" in element:
                        new_list = element.split(":")
                        clade = new_list[1][1:-2]

    with open(output_dir + "/final_summary.txt", "a") as fout:
        fout.write("Clade: " + clade + "\n")
        fout.write("G clade: " + g_clade)

    print("Clade: " + clade)
    print("G clade: " + g_clade)
    print("Pipeline finished. View results in the output file \"final_summary.txt\".")

def main():

    os.system(f"mkdir {path_to_subtype_detection}")
    determined_subtype = detect_subtype()

    os.system("mkdir " + artic_dir)
    os.system("mkdir " + duplication_dir)
    os.chdir(artic_dir)
    ref_file = detect_duplication(determined_subtype)

    os.system("mkdir " + reference_selection_directory)
    final_version = reference_selection(ref_file)

    os.chdir(artic_dir)
    artic_minion(final_version)

    os.system("mkdir " + nextclade_dir)
    nextclade(determined_subtype, nextclade_output)


parser = argparse.ArgumentParser()

path_to_python_file = os.path.abspath(os.path.dirname(__file__))
path_to_reference = path_to_python_file + "/references/"
path_to_primer_scheme =  path_to_python_file + "/primer_scheme/"
amplicon_length = "350, 900"
scheme_version = "RSV-2023"

parser.add_argument("-i", "--input", help = "Path to basecalled, demultiplexed fastq-files. It should end with the barcode directory (e.g. barcode15/)", required = True)
parser.add_argument("-s", "--sample", help = "Name of the sample", required = True)
parser.add_argument("-o", "--outputDir", help = "Output directory", required = True)
parser.add_argument("-m", "--medakaModel", help = "Medaka model that should be used for the artic pipeline (depends on basecaller used)", required = True)
parser.add_argument("-a", "--ampliconLength", help = "Minimum and maximum length of your amplicons comma separated (e.g. 350,900). Add 200 nt to your maximum length.", default = amplicon_length)
parser.add_argument("-p", "--schemeDir", help = "Path to primal scheme if the location of it was changed", default = path_to_primer_scheme)
parser.add_argument("-V", "--schemeVersion", help = "Name of your primer scheme version", default = scheme_version)
parser.add_argument("-r", "--refDir", help = "Path to directory containing the references if the location was changed", default = path_to_reference)
parser.add_argument("-n", "--nextcladeOutput", help = "Output file format of the nextclade results (tsv, csv, json, ndjson, all). Default: tsv", default = "tsv")


args = parser.parse_args()

path_to_reads = args.input
if path_to_reads[-1] == "/":
    path_to_reads = path_to_reads[:-1]
barcode = path_to_reads.split("/")[-1]
sample = args.sample
output_dir = args.outputDir
medaka_model = args.medakaModel
amplicon_length = args.ampliconLength.split(",")
path_to_primer_scheme = args.schemeDir
scheme_version = args.schemeVersion
if path_to_primer_scheme[-1] == "/":
    path_to_primer_scheme = path_to_primer_scheme[:-1]
path_to_reference = args.refDir
nextclade_output = args.nextcladeOutput
final_summary = "final_summary.txt"

path_to_subtype_detection = f"{output_dir}/{sample}_subtype_detection"
artic_dir = f"{output_dir}/{sample}_artic/"
duplication_dir = f"{output_dir}/{sample}_duplication/"
reference_selection_directory = f"{output_dir}/{sample}_reference_selection"
nextclade_dir = f"{output_dir}/{sample}_nextclade"


if __name__ == "__main__":
    main()
