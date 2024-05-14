import argparse
from typing import List, Dict
from os.path import  expanduser
from inputs_validation import ValidateFiles
from data_classes import Amplicon
from input_processing import InputProcessing
from model_manager import ModelManager
from classifier_report import ClasssifierReport
from map import ReadMapper


def main():

    parser = argparse.ArgumentParser(description='Classify reads in BAM file using existing model or train a model from bam files')
    parser.add_argument('-t','--target_regions', metavar='', type=str,
                        help='Bed file that specifies reference intervals to which reads where mapped', required=True)
    parser.add_argument('-r','--reference', metavar='', type=str,
                        help='FASTA file with the reference to which reads where mapped', required=True)
    parser.add_argument('-b','--bams', metavar='', type=str,
                        help='Directory with, list of, or individual BAM and corresponding BAM index files .bai', required=True)
    parser.add_argument('-m','--model', metavar='', type=str,
                        help='Pickle .pkl file containing pretrained model. Model must be trained on same reference', required=True)
    parser.add_argument('-f','--fastqs', metavar='', type=str,
                        help='Directory with ONT run results or individual FASTQ file', required=False)
    parser.add_argument('-d','--sample_descriptions', metavar='', type=str,
                        help='File with sample descritions, tab delimited, first column must be the BAM file name without .bam', required=False)
    parser.add_argument('-c','--description_column', metavar='', type=str,
                        help='Column in sample description file to use to augmnet samples descriptions', required=False)
    parser.add_argument('-g','--genotypes_hierarchy', type=str,
                        help='JSON file with genotypes hierarchy', required=False)
    parser.add_argument('--column_separator', metavar='', type=str,
                        help='The column separator for the sample descriptions file', required=False, default="\t")

    parser.add_argument('--cpus', type=int,
                        help='Directory for output files', required=False, default=1)
    parser.add_argument('-o','--output_file', metavar='', type=str,
                        help='File to store classification results', required=True)


    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        exit(0)

    cpu_to_use=int(args.cpus)

    input_processing=InputProcessing()

    bed_file, fasta_file, model_file = input_processing.get_ref_bed_model(args)
    if bed_file=="":
        exit(0)

    file_validator=ValidateFiles()
    if not file_validator.validate_bed(bed_file):
        exit(0)
    if not file_validator.validate_fasta(fasta_file):
        exit(0)
    if not file_validator.contigs_in_fasta(bed_file, fasta_file ):
        exit(0)

    if not args.sample_descriptions is None: #this has to be checked early in case file does not exist
        metadata_file=args.sample_descriptions
        if not input_processing.file_exists(metadata_file):
            exit(0)

    #Map the raw reads or collect bams to classify
    if not args.fastqs is None:
        if not input_processing.file_exists(args.fastqs):
            exit(0)
        bams_dir=expanduser(args.bams)
        mapper=ReadMapper("~/HandyReadGenotyper/temp_data",
                        args.fastqs,
                        fasta_file, cpu_to_use)
        if not mapper.map(bams_dir):
            exit()
        file_to_classify=[f.bam_file for f in mapper.results]
    else: #use the bams provided
        file_to_classify=input_processing.get_bam_files( args.bams )
    if len(file_to_classify)==0:
        exit(0)

    #exit(0)

    if not input_processing.check_address(args.output_file):
        exit(0)
    output_file=args.output_file

    if not args.genotypes_hierarchy is None:
        if not input_processing.file_exists(args.genotypes_hierarchy):
            exit(1)

    sample_labels: Dict[str, str] = {}
    if args.description_column is None and not args.sample_descriptions is None:
        print("Description file was provided (-d), but the column to use was not (-c)")        
        exit(0)
    if not args.description_column is None and args.sample_descriptions is None:
        print("Description column was provided (-c), but the file was not (-d)")
        exit(0)
    elif not args.description_column is None and not args.sample_descriptions is None:
        metadata_column=args.description_column
        metadata_sep=args.column_separator
        if not input_processing.check_column_in_descriptions(metadata_file,metadata_column, metadata_sep):
            exit(1)
        if not input_processing.get_sample_labels(metadata_file,metadata_sep, metadata_column, file_to_classify, sample_labels):
            exit(1)


    target_regions: List[Amplicon] = []
    with open(bed_file) as input_file:
        for line in input_file:
            target_regions.append(Amplicon.from_bed_line(line, fasta_file))

    model_manager=ModelManager(model_file, cpu_to_use)
    results=model_manager.classify_new_data(target_regions, file_to_classify)
    if not args.fastqs is None:
        report=ClasssifierReport(output_file, model_file, args.genotypes_hierarchy, sample_labels, mapping_results=mapper.results)
    else:
        report=ClasssifierReport(output_file, model_file, args.genotypes_hierarchy, sample_labels)
    report.create_report(results)

if __name__=="__main__":
    main()