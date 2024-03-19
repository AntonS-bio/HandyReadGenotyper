from os.path import expanduser
from sys import path
import argparse
from typing import List
import pandas as pd
path.append( expanduser("~/HandyAmpliconTool/scripts/") )
from inputs_validation import ValidateFiles
from data_classes import Amplicon
from model_manager import ModelManager
from input_processing import InputProcessing

# bed_file=expanduser("~/HandyReadGenotyper/InputData/chr_primers_v4_merged.bed")
# fasta_file=expanduser("~/HandyReadGenotyper/InputData/amr_chr_primer_v4_merged.fna")
# read_meta_data=pd.read_csv(expanduser("~/HandyReadGenotyper/Old_InputData/fastqs/metadata.tsv"), sep="\t", index_col=4)
# valid_files=read_meta_data.index[read_meta_data["PrimerSet"]=="Genotype"]
print("\n")

parser = argparse.ArgumentParser(description='Classify reads in BAM file using existing model or train a model from bam files')
parser.add_argument('-t','--target_regions', metavar='', type=str,
                    help='Bed file that specifies reference intervals to which reads where mapped', required=True)
parser.add_argument('-r','--reference', metavar='', type=str,
                    help='FASTA file with the reference to which reads where mapped', required=True)
parser.add_argument('-p','--positive_bams', metavar='', type=str,
                    help='Directory with or list of NON-target BAM files and corresponding BAM index files (.bai)', required=True)
parser.add_argument('-n','--negative_bams', metavar='', type=str,
                    help='Directory with or list of TARGET BAM files and corresponding BAM index files (.bai)', required=False)
parser.add_argument('-v','--vcf', metavar='', type=str,
                    help='VCF file containing positions that will be excluded from training as well as genotype defining SNPs (also excluded)', required=True)
parser.add_argument('-o','--output_dir', metavar='', type=str,
                    help='Directory for output files', required=True)
parser.add_argument('-m', '--bams_matrix', metavar='', type=str, default="",
                    help='Matrix with precalculated BAM matrices ()', required=False)

try:
    args = parser.parse_args()
except:
    parser.print_help()
    exit(0)

input_processing=InputProcessing()
positive_bams=input_processing.get_bam_files( args.positive_bams )
negative_bams=input_processing.get_bam_files( args.negative_bams )
if len(positive_bams)==0 or len(negative_bams)==0:
    exit()

bed_file = expanduser(args.target_regions)
if not input_processing.file_exists(bed_file):
    exit(0)

fasta_file = expanduser(args.reference)
if not input_processing.file_exists(fasta_file):
    exit(0)

vcf_file=expanduser(args.vcf)
if not input_processing.file_exists(vcf_file):
    exit(0)

output_dir=args.output_dir
input_processing.make_dir(output_dir)
model_file=f'{output_dir}/model.pkl'
model_quality_file=f'{output_dir}/quality.tsv'

use_existing_bams_matrix=False
if args.bams_matrix!="":
    use_existing_bams_matrix=True
    matrix_file=expanduser(args.bams_matrix)
    if not input_processing.file_exists(matrix_file):
        exit(0)
else:
    matrix_file=f'{output_dir}/bams_matrix.pkl'


file_validator=ValidateFiles()
file_validator.validate_bed(bed_file)
file_validator.validate_fasta(fasta_file)
file_validator.contigs_in_fasta(bed_file, fasta_file )

target_regions: List[Amplicon] = []
with open(bed_file) as input_file:
    for line in input_file:
        target_regions.append(Amplicon.from_bed_line(line, fasta_file))


if True: #Model training
    model_manager=ModelManager(model_file)
    model_manager.model_evaluation_file=model_quality_file
    model_manager.load_genotype_snps(vcf_file)
    if use_existing_bams_matrix:
        model_manager.train_from_existing_matrices(target_regions=target_regions, matrices_file=matrix_file)
    else:
        model_manager.train_from_bams(target_regions, positive_bams, negative_bams, matrix_file)
    

