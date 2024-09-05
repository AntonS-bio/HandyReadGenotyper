import argparse
from typing import List, Dict
import pickle
from os.path import  expanduser, join, exists
from os import mkdir
from shutil import which, rmtree
import uuid 
from data_classes import Amplicon
from input_processing import InputProcessing
from model_manager import ModelManager
from classifier_report import ClasssifierReport
from read_classifier import Classifier, ReadsMatrix
from map import ReadMapper

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

def classify(temp_dir):

    parser = argparse.ArgumentParser(description='Classify reads in BAM file using existing model or train a model from bam files')
    parser.add_argument('-f','--fastqs', type=str,
                        help='Directory with ONT run results or individual FASTQ file', required=False)
    parser.add_argument('-b','--bams',  type=str,
                        help='Directory with, list of, or individual BAM and corresponding BAM index files (.bai)', required=True)
    parser.add_argument('-m','--model', type=str,
                        help='Pickle (.pkl) file containing pretrained model. Model must be trained on same reference', required=True)
    parser.add_argument('-d','--sample_descriptions',  type=str,
                        help='File with sample descriptions (tab delimited), first column must be the BAM file name without .bam', required=False)
    parser.add_argument('-c','--description_column',  type=str,
                        help='Column in sample description file to use to augmnet samples descriptions', required=False)
    parser.add_argument('--column_separator',  type=str,
                        help='The column separator for the sample descriptions file', required=False, default="\t")
    parser.add_argument('-g','--genotypes_hierarchy', type=str, help='JSON file with genotypes hierarchy', required=False)
    parser.add_argument('--cpus', type=int, help='Number of CPU cores to use', required=False, default=1)
    parser.add_argument('-o','--output_file', type=str,
                        help='File to store classification results', required=True)
    parser.add_argument('-s', '--max_soft_clip', type=int,
                        help='Specifies maximum permitted soft-clip (length of read before first and last mapped bases) of mapped read', required=False, default=25)
    parser.add_argument('-l', '--max_read_len_diff', type=int,
                        help='Specifies maximum nucleotide difference between mapped portion of read and target amplicon', required=False, default=10)
    parser.add_argument('--organism_presence_cutoff', type=float,
                        help='Sample is reported as having target organism if at least this percentage  of non-transient amplicons have >10 reads. Values 0-100.', required=False, default=20.0)
   
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1.24')


    try:
        args = parser.parse_args()
    except Exception as e:
        parser.print_help()
        return

    mkdir(temp_dir)

    if args.organism_presence_cutoff<1 and args.organism_presence_cutoff!=0:
        print(f'Option "organism_presence_cutoff" has value {args.organism_presence_cutoff} which is below 1, but not 0. Valid values are 0 to 100.')
        return

    ReadsMatrix.permitted_read_soft_clip=int(args.max_soft_clip)
    ReadsMatrix.permitted_mapped_sequence_len_mismatch=int(args.max_read_len_diff)

    #check minimap2 is installed
    if which("minimap2") is None:
        print(f'Missing minimap2 program. It is available via Bioconda.')
        return

    cpu_to_use=int(args.cpus)

    input_processing=InputProcessing()

    model_file = args.model


    fasta_file = join(temp_dir,"reference.fasta")
    generate_ref_fasta(model_file, fasta_file)
    generate_amplicons(model_file, fasta_file)

    if not args.sample_descriptions is None: #this has to be checked early in case file does not exist
        metadata_file=expanduser(args.sample_descriptions)
        if not input_processing.file_exists(metadata_file):
            return
        metadata_column, metadata_sep = input_processing.check_metadata(args, metadata_file)


    #Map the raw reads or collect bams to classify
    if not args.fastqs is None:
        if not input_processing.file_exists(args.fastqs):
            return
        bams_dir=expanduser(args.bams)
        mapper=ReadMapper(temp_dir,
                        args.fastqs,
                        model_file, fasta_file, cpu_to_use)
        
        #check that sample labels are present in metadatafile
        sample_labels: Dict[str, str] = {}
        if not  args.sample_descriptions is None:
            files_to_check=mapper.output_names(bams_dir)
            if not input_processing.get_sample_labels(metadata_file,metadata_sep, metadata_column, files_to_check, sample_labels):
                return #terminate if some sample is not present in metadata.

        if not mapper.map(bams_dir):
            return
        file_to_classify=[f.bam_file for f in mapper.results]
    else: #use the bams provided
        file_to_classify=input_processing.get_bam_files( args.bams )
        sample_labels: Dict[str, str] = {}
        if not  args.sample_descriptions is None:
            files_to_check=mapper.output_names(bams_dir)
            if not input_processing.get_sample_labels(metadata_file,metadata_sep, metadata_column, file_to_classify, sample_labels):
                return #terminate if some sample is not present in metadata.

    if len(file_to_classify)==0:
        print("No files to classify. Please check that you have specified correct FASTQs directory if using -f or existing BAMs if not using -f")
        return
    

    output_file=input_processing.check_output_dir(args.output_file)

    if not args.genotypes_hierarchy is None:
        if not input_processing.file_exists(args.genotypes_hierarchy):
            return



    target_regions: List[Amplicon] = generate_amplicons(model_file, fasta_file)

    model_manager=ModelManager(model_file, cpu_to_use)
    results=model_manager.classify_new_data(target_regions, file_to_classify)

    # with open(expanduser("~/HandyReadGenotyper/temp/read_ids.tsv"),"w") as output:
    #     output.write("File\tAmplicon\tID\n")
    #     for result in results:
    #         for read_id in result.positive_read_ids:
    #             output.write(result.sample+"\t"+result.amplicon.name+"\t"+read_id+"\n")

    if not args.fastqs is None:
        report=ClasssifierReport(output_file, model_file, args.organism_presence_cutoff, args.genotypes_hierarchy, sample_labels, mapping_results=mapper.results)
    else:
        report=ClasssifierReport(output_file, model_file, args.organism_presence_cutoff, args.genotypes_hierarchy, sample_labels)
    report.create_report(results)
    rmtree(temp_dir)
    
def main():
    temp_dir=expanduser( join("./",str(uuid.uuid4())) )
    #try:
    classify(temp_dir)
    #except Exception as e:
    #    print(e)
    if exists(temp_dir): #clean up
        rmtree(temp_dir)

if __name__=="__main__":
    main()
