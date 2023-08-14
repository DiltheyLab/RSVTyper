# RSVTyper
Genotyping RSV samples from nanopore sequencing data

# Installation

## Dependencies
- minimap2
- samtools >= 1.10
- nextclade 
- artic = 1.2.1

## Bioconda
This package can be installed via bioconda. All dependencies will be installed as well. 
`conda install -c bioconda rsv-typer`

# Usage 
To use the pipeline run the following command:
`rsv-typer -i /path/to/reads -s sample_name -o /path/to/output/directory -m medaka_model`
The input (-i) should consist of basecalled, demultiplexed nanopore sequencing reads that have passed the quality check. 
The sample (-s) determines the prefix of the output directories within the given output directory.
The medaka model (-m) should be chosen according to the version of the basecaller used. For further information on medaka models see:
https://github.com/nanoporetech/medaka
The output will consist of four different directories containing the output of the different steps in the pipeline and a summary file in which the subtype, the reference file used for the artic pipeline and the genotype can be found. 
