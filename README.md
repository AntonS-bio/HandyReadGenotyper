# HandyReadGenotyper
Tool for genotyping Oxford Nanopore amplicon sequencing data

## Purpose
This tool is intended to be used with a partner HandyAmpliconTool which helps to design primers for environmental surveilance. The partner tools designs primers for specific genotypes (based on SNPs unique to those genotypes) while minimising risk of amplyfing other DNA that might be present in the environment.

This tool, HandyReadGenotyper, has two modes: 
    1. Train a classification model using existing Nanopore data - this can be your own data, data from genomic repositories (ENA, NCBI, DDBJ) or a mix of these
    2. Classify ONT sequencing data using the model trained in (1)

## Setup
The easiest way to setup the tool is to use conda or mamba. If you are using macOS or Linux, you can install via command line. If you are using Windows 10+, the best option is to setup Windows Subsystem Linux (WSL) first which will give you access to Linux functionality from you Windows system. After that, you can use conda or mamba.

The best practice is to install packages in dedicated environment to avoid software conflicts. To create new environment and install HandyReadGenotyper into it use:
```
conda create --name  hrgENV -c bioconda -c conda-forge handyreadgenotyper
```
Once installed, use
```
conda activate hrgENV
```
to activate the environment. Now you are ready to use the HandyReadGenotyper. You would need to run activation command (but not create command) every time you start a new terminal. In both commands above "hrgENV" can be replaced with whatever you want to call the environment. 

You can also install the tool without creating new environment
```
conda install -c bioconda -c conda-forge handyreadgenotyper
```
but this is not recommended.

### Model training
```
usage: train.py -t  -r  -p  -v  -o  [-m] [-h] [-n] 

Classify reads in BAM file using existing model or train a model from bam files

options:
  -h, --help            show this help message and exit
  -t , --target_regions 
                        Bed file that specifies reference intervals to which reads where mapped
  -r , --reference      FASTA file with the reference to which reads where mapped
  -p , --positive_bams 
                        Directory with or list of NON-target BAM files and corresponding BAM index files (.bai)
  -n , --negative_bams 
                        Directory with or list of TARGET BAM files and corresponding BAM index files (.bai)
  -v , --vcf            VCF file containing positions that will be excluded from training as well as genotype defining SNPs (also excluded)
  -o , --output_dir     Directory for output files
  -m , --bams_matrix    Matrix with precalculated BAM matrices
```
### Reads classification (must have a trained model)
```
usage: classify.py -t  -r  -b  -m  -o  [-d] [-h]

Classify reads in BAM file using existing model or train a model from bam files

options:
  -h, --help            show this help message and exit
  -t , --target_regions 
                        Bed file that specifies reference intervals to which reads where mapped
  -r , --reference      FASTA file with the reference to which reads where mapped
  -b , --bam_files      Directory with, list of, or individual BAM and corresponding BAM index files (.bai)
  -d , --sample_descriptions 
                        File with sample descritions (tab delimited)
  -m , --model          Pickle (.pkl) file containing pretrained model. Model must be trained on same reference
  -o , --output_file    File to store classification results

```
