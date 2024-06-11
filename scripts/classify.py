import argparse
from typing import List, Dict
import pickle
from os.path import  expanduser, join, exists
from os import mkdir
from inputs_validation import ValidateFiles
from data_classes import Amplicon
from input_processing import InputProcessing
from model_manager import ModelManager
from classifier_report import ClasssifierReport
from read_classifier import Classifier
from map import ReadMapper
from shutil import which, rmtree
import uuid 

def generate_amplicons(model_file: str, fasta_file:str) -> List[Amplicon]:
    target_regions: List[Amplicon] = []
    with open(model_file, "rb") as input_model:
        model_manager: Dict[str, Classifier] =pickle.load(input_model)
    for name, model in model_manager.items():
        bed_line=f'{model.name}\t0\t{len(model.nucletoide_seq)}\t{model.name}'
        target_regions.append(Amplicon.from_bed_line(bed_line, fasta_file))
    return target_regions

def generate_ref_fasta(model_file: str, ref_file: str) -> None:
    with open(model_file, "rb") as input_model:
        model_manager: Dict[str, Classifier] =pickle.load(input_model)
    with open(ref_file, "w" ) as ref_output:
        for name, model in model_manager.items():
            ref_output.write(">"+str(model.name)+"\n")
            ref_output.write(str(model.nucletoide_seq)+"\n")

def main(temp_dir):

    parser = argparse.ArgumentParser(description='Classify reads in BAM file using existing model or train a model from bam files')
    parser.add_argument('-f','--fastqs', metavar='', type=str,
                        help='Directory with ONT run results or individual FASTQ file', required=False)
    parser.add_argument('-b','--bams', metavar='', type=str,
                        help='Directory with, list of, or individual BAM and corresponding BAM index files (.bai)', required=True)
    parser.add_argument('-m','--model', metavar='', type=str,
                        help='Pickle (.pkl) file containing pretrained model. Model must be trained on same reference', required=True)
    parser.add_argument('-d','--sample_descriptions', metavar='', type=str,
                        help='File with sample descriptions (tab delimited), first column must be the BAM file name without .bam', required=False)
    parser.add_argument('-c','--description_column', metavar='', type=str,
                        help='Column in sample description file to use to augmnet samples descriptions', required=False)
    parser.add_argument('--column_separator', metavar='', type=str,
                        help='The column separator for the sample descriptions file', required=False, default="\t")
    parser.add_argument('-g','--genotypes_hierarchy', type=str,
                        help='JSON file with genotypes hierarchy', required=False)
    parser.add_argument('--cpus', type=int,
                        help='Directory for output files', required=False, default=1)
    parser.add_argument('-o','--output_file', metavar='', type=str,
                        help='File to store classification results', required=True)
    parser.add_argument('--version', action='version', version='%(prog)s 0.1.15')

    mkdir(temp_dir)
    try:
        args = parser.parse_args()
    except:
        exit(0)

    #check minimap2 is installed
    if which("minimap2") is None:
        print(f'Missing minimap2 program. It is available via Bioconda.')
        exit(0)

    cpu_to_use=int(args.cpus)

    input_processing=InputProcessing()

    model_file = args.model


    fasta_file = join(temp_dir,"reference.fasta")  
    generate_ref_fasta(model_file, fasta_file)  
    generate_amplicons(model_file, fasta_file)

    if not args.sample_descriptions is None: #this has to be checked early in case file does not exist
        metadata_file=expanduser(args.sample_descriptions)
        if not input_processing.file_exists(metadata_file):
            exit(0)
        metadata_column, metadata_sep = input_processing.check_metadata(args, metadata_file)


    #Map the raw reads or collect bams to classify
    if not args.fastqs is None:
        if not input_processing.file_exists(args.fastqs):
            exit(0)
        bams_dir=expanduser(args.bams)
        mapper=ReadMapper(temp_dir,
                        args.fastqs,
                        model_file, fasta_file, cpu_to_use)
        if not mapper.map(bams_dir):
            exit()
        file_to_classify=[f.bam_file for f in mapper.results]
    else: #use the bams provided
        file_to_classify=input_processing.get_bam_files( args.bams )
    if len(file_to_classify)==0:
        print("No files to classify. Please check that you have specified correct FASTQs directory if using -f or existing BAMs if not using -f")
        exit(0)
    #exit(0)

    output_file=input_processing.check_output_dir(args.output_file)

    if not args.genotypes_hierarchy is None:
        if not input_processing.file_exists(args.genotypes_hierarchy):
            exit(1)

    sample_labels: Dict[str, str] = {}
    if not  args.sample_descriptions is None:
        if not input_processing.get_sample_labels(metadata_file,metadata_sep, metadata_column, file_to_classify, sample_labels):
            exit(1)

    target_regions: List[Amplicon] = generate_amplicons(model_file, fasta_file)

    model_manager=ModelManager(model_file, cpu_to_use)
    results=model_manager.classify_new_data(target_regions, file_to_classify)
    if not args.fastqs is None:
        report=ClasssifierReport(output_file, model_file, args.genotypes_hierarchy, sample_labels, mapping_results=mapper.results)
    else:
        report=ClasssifierReport(output_file, model_file, args.genotypes_hierarchy, sample_labels)
    report.create_report(results)
    rmtree(temp_dir)
    

if __name__=="__main__":
    temp_dir=expanduser( join("./",str(uuid.uuid4())) )
    try:
        main(temp_dir)
    except:
        if exists(temp_dir):
            rmtree(temp_dir)