#!/usr/bin/env python
""" This script runs the artic pipeline with Nanopore reads that have already been basecalled and demultiplexed. Variant calling is done by the experimental medaka pipeline.
It is specificially made for RSV samples as it determines the subtype and chooses the reference that matches the given reads the most."""
import os
import sys
import argparse
import re

def main():
    parser = argparse.ArgumentParser()

    path_to_python_file = os.path.abspath(os.path.dirname(__file__))
    # Path that contains all references. This includes a combined reference file containing all references (20), a reference file containing only references with the subtype A (10), a reference file
    # containing only references with the subtype B (10) and all reference files separately.
    path_to_reference = path_to_python_file + "/references/"
    # Path that contains the primer scheme "RSV-2023"
    path_to_primer_scheme =  path_to_python_file + "/primer_scheme/"

    parser.add_argument("-i", "--input", help = "Path to basecalled, demultiplexed fastq-files. It should end with the barcode directory (e.g. barcode15/)", required = True)
    parser.add_argument("-s", "--sample", help = "Name of the sample", required = True)
    parser.add_argument("-o", "--outputDir", help = "Output directory", required = True)
    parser.add_argument("-m", "--medakaModel", help = "Medaka model that should be used for the artic pipeline (depends on basecaller used)", required = True)
    parser.add_argument("-p", "--schemeDir", help = "Path to primal scheme if the location of it was changed", default = path_to_primer_scheme)
    parser.add_argument("-r", "--refDir", help = "Path to directory containing the references if the location was changed", default = path_to_reference)
    parser.add_argument("-n", "--nextcladeOutput", help = "Output file format of the nextclade results (tsv, csv, json, ndjson, all). Default: tsv", default = "tsv")
    args = parser.parse_args()

    path_to_reads = args.input
    sample = args.sample
    output_dir = args.outputDir
    medaka_model = args.medakaModel
    path_to_primer_scheme = args.schemeDir
    path_to_reference = args.refDir
    nextclade_output = args.nextcladeOutput


    # Subtype detection (subtype A or subtype B)
    # Directory for the output of the subtype detection
    path_to_subtype_detection = f"{output_dir}/{sample}_subtype_detection"
    os.system(f"mkdir {path_to_subtype_detection}")
    # The sequences in this merged reference file are consensus sequences of cluster sequences. The clusters were made using gisaid sequences without Ns that were sampled after 2014.
    reference_subtype = "all_consensus_references.fasta"
    coverage_summary = f"{sample}_coverage_summary.txt"
    alignment = f"{sample}_vs_{reference_subtype[:-6]}"

    # The reads are mapped to a combined fasta file of all reference sequences containing both subtypes.
    minimap_command = f"minimap2 -ax map-ont {path_to_reference}{reference_subtype} {path_to_reads}/*fastq* > {path_to_subtype_detection}/{alignment}.sam"
    os.system(minimap_command)

    # Only primary alignments are used for the subtype detection.
    primary_alignment = f"{alignment}_primary_only"
    with open(f"{path_to_subtype_detection}/{primary_alignment}.sam", "w") as fout:
        with open(f"{path_to_subtype_detection}/{alignment}.sam", "r") as fin:
            for line in fin:
                # Any secondary or inversed alignment is not written into the new file. Only the header and primary alignments will remain.
                  if "tp:A:S" not in line and "tp:A:I" not in line:
                        fout.write(line)

    # Conversion to bam file and creation of an index file which is needed for samtools coverage
    convert_to_bam = f"samtools view -1 -S -b {path_to_subtype_detection}/{primary_alignment}.sam > {path_to_subtype_detection}/{primary_alignment}.bam"
    sort = f"samtools sort {path_to_subtype_detection}/{primary_alignment}.bam -o {path_to_subtype_detection}/{primary_alignment}.sorted.bam"
    index = f"samtools index {path_to_subtype_detection}/{primary_alignment}.sorted.bam"
    command_chain = f"{convert_to_bam} && {sort} && {index}"
    os.system(command_chain)

    # Executing samtools coverage for every reference in the region 3400-5400. This region has been selected as it seems to have diverged between the subtypes.
    # The output for every reference is written into a combined file.
    with open(f"{path_to_subtype_detection}/{coverage_summary}", "w") as fout:
        for i in range(1, 11):
            consensus_seq = "cluster_A_" + str(i) + "_consensus"
            samtools_coverage = f"samtools coverage -r {consensus_seq}:3500-5500 {path_to_subtype_detection}/{primary_alignment}.sorted.bam"
            creating_coverage_summary = f"({samtools_coverage}) >> {path_to_subtype_detection}/{coverage_summary}"
            os.system(creating_coverage_summary)
        for i in range(1, 11):
            consensus_seq = "cluster_B_" + str(i) + "_consensus"
            samtools_coverage = f"samtools coverage -r {consensus_seq}:3500-5500 {path_to_subtype_detection}/{primary_alignment}.sorted.bam"
            creating_coverage_summary = f"({samtools_coverage}) >> {path_to_subtype_detection}/{coverage_summary}"
            os.system(creating_coverage_summary)

    # This ensures that the program continues only after the output of all references has been written into the coverage summary file.
    lines = 0
    while lines < 40:
        with open(f"{path_to_subtype_detection}/{coverage_summary}", "r") as fin3:
                lines = sum(1 for line3 in fin3)

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

    # If there is no reference with any reads mapped across all references of one subtype the value is set to 1 to prevent an error that occurs when a number is divided by 0.
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

    # The subtype and the final reference chosen for the artic pipeline are written into a file.
    final_summary = "final_summary.txt"
    with open(output_dir + "/" + final_summary, "w") as fout:
        fout.write("Subtype: " + final_subtype + "\n")



    # Selection of the reference that matches your reads best. The reference with the most primary alignments on the whole genome will be selected.
    # The output of the reference selection is written into a new directory.
    reference_selection_directory = f"{output_dir}/{sample}_reference_selection"
    os.system("mkdir " + reference_selection_directory)
    # Mapping the reads to a combined fasta file containing all references of the detected subtype.
    paf_output = f"{sample}_vs_cluster_{final_subtype}.paf"
    minimap_command_ref = f"minimap2 -c -x map-ont {path_to_reference}cluster_{final_subtype}_consensus_all.fasta {path_to_reads}/*fastq* > {reference_selection_directory}/{paf_output}"
    os.system(minimap_command_ref)

    # This counts the number of primary alignments for each reference.
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

    # Sorts by descending number of primary alignments. Therefore, the last element in the sorted directory is the reference with the most primary alignments.

    output_ref_select = f"{sample}_primary_alignments.txt"
    sorted_primary_alignment_dict = sorted(primary_alignment_dict.items(), key=lambda x:x[1])
    final_reference = sorted_primary_alignment_dict[-1][0]

    # Creates a file with the number of primary alignments for each reference
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

    # Running the artic pipeline in the conda environment
    artic_dir = f"{output_dir}/{sample}_artic"
    os.system("mkdir " + artic_dir)
    os.chdir(artic_dir)
    os.system(f"python3 {path_to_python_file}/artic_pipeline.py -i {path_to_reads} -m {medaka_model} -p {path_to_primer_scheme} -v {version} -s {sample} -o {output_dir}")
    os.chdir("../../")
    # Running nextclade
    nextclade_dir = f"{output_dir}/{sample}_nextclade"
    os.system("mkdir " + nextclade_dir)
    subtype = final_subtype.lower()
    consensus_seq = f"{sample}.consensus.fasta"

    os.system(f"nextclade run -d rsv_{subtype} -O {nextclade_dir} -s={nextclade_output} {artic_dir}/{consensus_seq} --output-basename '{sample}_nextclade'")


    # Writing out the clades into final_summary.txt
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

if __name__ == "__main__":
    main()
