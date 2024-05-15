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



## Model training
If you have been given a model pickle file - you don't have to do this, so go to Read Classification section. If you haven't been given one, ask if anyone on your project already has it - chances are bioinformatician does.

The idea of training is that to distinguish target and non-target organism, HandyReadGenotyper needs to see a lot of examples of amplicon sequencing data from both. The examples of target organisms should be fairly easy to find - your own sequencing data should work. For non-target organism, you probably won't be able to generate sufficient data yourself, so you have to rely on public data. SRA database at NCBI is a good source. I recommend finding taxa for your target organism/serovar/strain and taking all available ONT sequencing data within two-three taxonomic layers above. For example, if you are looking at Salmonella, you may decide to take all ONT data from all taxa within Enterobacterales. Once raw ONT data is mapped to FASTQs these will be your negative BAMs. If in your test sequencing runs you often see an off-target amplification of same organism, add it to a negative set.

Once you have mapped the the sequencing data you've amassed to reference amplicons, you can train your model. The training process is quick, but generating BAM files for it can take several days. This is not related to HandyReadGenotyper, but to volume of data you need to process. However, training only needs to be done once and the "train" function itself below should take fewer than 10 minutes to run. It will generate a model file that you can share with colleagues to avoid them having to train the model.

```
usage: train -t  -r  -p  -v  -o  [-m] [-h] [-n] 

Classify reads in BAM file using existing model or train a model from bam files. You will get an error if you use "train.py" instead of "train"

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

## Reads classification (must have a trained model, see Model training above)
The purpose of classifier is to take all reads provided to it and identify which came from your target and which from non-target. The classifier will use extra information you provide to generate a nice HTML report file with results - you can see sample here in test_data/report.html. If you download this file it can be opened in any browser. 

**IMPORTANT all BAMs must be indexed using samtools index command.**

To run classifier you have to main options. First, if you have already mapped FASTQs to amplicons reference sequence, you can provide the report with either directory with BAM files, a single BAM file or a plain text file with list of BAM files (full address including directory) to classify. Examples below use supplied test data from test_data.


### I have BAMS files, I already mapped FASTQs to references

```
classify -t ./amplicons.bed -r ./amplicons.fna -b ./bams/ -m model.pkl -d metadata.tsv -c Description -g hierarchy.json -o report.html
```
Maps all BAMs in directory ./bams/ while this will only classify sample_1.bam
```
classify -t ./amplicons.bed -r ./amplicons.fna -b ./bams/sample_1.bam -m model.pkl -d metadata.tsv -c Description -g hierarchy.json -o report.html
```




### I only have FASTQ files, but no BAMs

If you don't have bams, you can use "-f" option to map FASTQs using minimap2 to the amplicon reference. The command only differs from above by this extra parameter. Not that "-b" is still present - it will be a directory to which the mapped BAMs will be saved. 
```
classify -t ./amplicons.bed -r ./amplicons.fna -b ./bams/ -f ./fastqs/ -m model.pkl -d metadata.tsv -c Description -g hierarchy.json -o report.html
```
If you look at structure of test_data/fastqs you will see it's a bit odd. It's very difficult to capture all possible ways the FASTQs will be provided, but the most likely one if a result of single run from ONT devices. When these were multiplexed (most likely use case) ONT device will create a directory for each barcode and will place multiple FASTQ files for this barcode into that directory. It will likely look like this:
```
Run_20
      ↳ barcode1
               ↳ FAX83461_pass_barcode1_503f1897_ffa628d1_1.fastq.gz
               ↳ FAX83461_pass_barcode1_503f1897_ffa628d1_2.fastq.gz
               ↳ FAX83461_pass_barcode1_503f1897_ffa628d1_3.fastq.gz
    
      ↳ barcode2
               ↳ FAX83461_pass_barcode2_503f1897_ffa628d1_1.fastq.gz
               ↳ FAX83461_pass_barcode2_503f1897_ffa628d1_2.fastq.gz
               ↳ FAX83461_pass_barcode2_503f1897_ffa628d1_3.fastq.gz
    
      ↳ barcode3
               ↳ FAX83461_pass_barcode3_503f1897_ffa628d1_1.fastq.gz
               ↳ FAX83461_pass_barcode3_503f1897_ffa628d1_2.fastq.gz
               ↳ FAX83461_pass_barcode3_503f1897_ffa628d1_3.fastq.gz
```
Normally, you'd need to merge the the files from each barcode before mapping, but HandyReadGenotype will do it for you.

**IMPORTANT if HandyReadGenotyper is doing the mapping, it will call each output BAM by the name with corresponding barcode (ex. barcode2.bam), this may overwrite some old BAMs**

Neither the metadata file (-d), nor genotypes hierarchy (-g) are required, but you probably already have them and they substantially enrich the classification report, so it's worth using them. The metadata file can be where you keep track of sequencing results, the genotypes hierarchy file was likely used when your PCR primers were generated. 

```
usage: classify -t  -r  -b  -m  -o  [-d] [-c] [--column_separator] [-g] [--cpu] [-h]

Classify reads in BAM file using existing model or train a model from bam files. You will get an error if you use "classify.py" instead of "classify"

options:
  -h, --help            show this help message and exit
  -t , --target_regions 
                        Bed file that specifies reference intervals to which reads where mapped
  -r , --reference      FASTA file with the reference to which reads where mapped
  -b , --bam_files      Directory with, list of, or individual BAM and corresponding BAM index files (.bai)
  -m , --model          Pickle (.pkl) file containing pretrained model. Model must be trained on same reference
  -o , --output_file    Name of HTML file to store classification results
  -d , --sample_descriptions 
                        File with sample descritions (tab delimited), first column must be the BAM file name without .bam
  -c , --description_column 
                        Column in sample descriptions file to use to augment samples descriptions
  --column_separator    The column separator for the sample descriptions file, default=tab
  -g , --genotypes_hierarchy
                        JSON file with genotypes hierarchy
  --cpus                Number of CPUs to use, default=1

```
