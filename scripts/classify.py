from os.path import expanduser
from sys import path
import argparse
from typing import List
from inputs_validation import ValidateFiles
from data_classes import Amplicon
from input_processing import InputProcessing
from model_manager import ModelManager


def main():

    parser = argparse.ArgumentParser(description='Classify reads in BAM file using existing model or train a model from bam files')
    parser.add_argument('-t','--target_regions', metavar='', type=str,
                        help='Bed file that specifies reference intervals to which reads where mapped', required=True)
    parser.add_argument('-r','--reference', metavar='', type=str,
                        help='FASTA file with the reference to which reads where mapped', required=True)
    parser.add_argument('-b','--bam_files', metavar='', type=str,
                        help='Directory with, list of, or individual BAM and corresponding BAM index files (.bai)', required=True)
    parser.add_argument('-d','--sample_descriptions', metavar='', type=str,
                        help='File with sample descritions (tab delimited)', required=False)
    parser.add_argument('-m','--model', metavar='', type=str,
                        help='Pickle (.pkl) file containing pretrained model. Model must be trained on same reference', required=True)
    parser.add_argument('-o','--output_file', metavar='', type=str,
                        help='File to store classification results', required=True)


    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        exit(0)

    input_processing=InputProcessing()
    file_to_classify=input_processing.get_bam_files( args.bam_files )
    if len(file_to_classify)==0:
        exit(0)

    bed_file, fasta_file, model_file = input_processing.get_ref_bed_model(args)
    if bed_file=="":
        exit(0)

    if not input_processing.check_address(args.output_file):
        exit(0)
    output_file=args.output_file

    file_validator=ValidateFiles()
    file_validator.validate_bed(bed_file)
    file_validator.validate_fasta(fasta_file)
    file_validator.contigs_in_fasta(bed_file, fasta_file )

    target_regions: List[Amplicon] = []
    with open(bed_file) as input_file:
        for line in input_file:
            target_regions.append(Amplicon.from_bed_line(line, fasta_file))

    model_manager=ModelManager(model_file)
    model_manager.classify_new_data(target_regions, file_to_classify, output_file)


if __name__=="__main__":
    main()