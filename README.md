# Snakemake workflow: RNASeq
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.15.0-brightgreen.svg)](https://snakemake.bitbucket.io)

A Snakemake-based pipeline for **RNASeq data analysis**.

Starting from `fastq` files, the pipeline merges files from different `units` and perform reads `quality trimming`.

The pseudoaligner `kallisto` is used to estimate transcripts abundance, with resulting `.h5` files that can be imported into `DESeq2`for DE Analysis.

`STAR` 2-pass mapping is used for read alignment.

## Authors

* [Matteo Massidda](https://github.com/massiddamt), University of Sassari
* [Vincenzo Rallo](https://github.com/VincenzoRallo), Institute for Genetic and Biomedical Research (IRGB) - National Research Council (CNR)

## INSTRUCTIONS
Create a virtual environment with the command:
```commandline
mamba create -c bioconda -c conda-forge --name snakemake snakemake=6.15 snakedeploy
```
and activate it:
```commandline
conda activate snakemake
```
You can perform the pipeline deploy defining a directory `my_dest_dir` for analysis output and a pipeline tag for a specific version:
```bash
snakedeploy deploy-workflow https://github.com/GeneBANGS/RNASeq.git 
                    my_desd_dir 
                    --tag v1.0
```
To run the pipeline, go inside the deployed pipeline folder and use the command:
```bash
snakemake --use-conda -p --cores all
```