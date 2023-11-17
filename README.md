# Shotgun-metagenomic-pipeline

This repository hosts a pipeline for analysing shotgun metagenomic data. The pipeline is developed with Nextflow and is executed using Singularity containers.

# Running the pipeline

The pipeline is intended to be executed on a cluster or server. The pipeline provides configuration for a PBSPro cluster. A job submission script for running
the pipeline on a PBSPro cluster is available. The command for running the pipeline is:

`nextflow run metagen.nf`

Options:

`-c` [REQUIRED] path to the configuration file

`-profile` [REQUIRED] server or cluster. If running on a cluster enter cluster. If running on a server enter server

`--reads` [REQUIRED] Path to the input raw reads. Ensure that the path is encapsulated in double quotes ("")

`--blast_db` [REQUIRED] Full path to the BLAST database

`--deepARG` [REQUIRED] Full path to the deepARG database

`--quastref` [REQUIRED] Full path to the reference file

`--outdir` [OPTIONAL] directory where pipeline will write all output file, default name is "results"

# Requirements for running pipeline

Nextflow version 21+

Singularity version 3.5+

# Input

The pipeline takes raw paired-end reads in fastq or fasta format.

# Output

A directory named <results>, unless otherwise specified, will be generated in the current working directory. Notable output are:

BLAST tables

deepARG tables

Annotation files

Binning statistics reports

# Software algorithms

Software algorithms are available through Singularity containers. Definition files specify how the containers were developed and dependencies for each algorithm.
