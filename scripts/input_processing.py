from typing import List, Tuple, Dict
from os.path import expanduser, isdir, exists, dirname, join
from os import listdir, mkdir
from pathlib import Path
import warnings

class InputProcessing:

    def __init__(self) -> None:
        pass

    def get_bam_files(self, bam_source: str) -> List[str]:
        bam_source=expanduser(bam_source)
        file_to_classify=[]
        if isdir(bam_source):  #Source is directory with BAM files
            file_to_classify=[join(bam_source,f) for f in listdir(bam_source)  if
                            f.split(".")[-1]=="bam"]
            if len(file_to_classify)==0:
                print(f'No files ending with .bam found in {bam_source}')
        elif bam_source.split(".")[-1]=="bam": #Source is single BAM file
            if exists(bam_source):
                file_to_classify.append(bam_source)
            else:
                print(f'BAM file {bam_source} does not exist. Maybe the location is incorrect?')
        else:
            #check that list of bams contains valid bams
            with open(bam_source) as bam_list:
                for line in bam_list:
                    tentative_bam=expanduser(line.strip())
                    if exists(tentative_bam):
                        file_to_classify.append(tentative_bam)
                    else:
                        print(f'BAM file {tentative_bam} does not exist. Maybe the location is incorrect?')
        return file_to_classify
    
    def check_column_in_descriptions(self, description_file: str, column: str, separator: str) -> bool:
        with open(description_file) as input_file:
            first_line=input_file.readline().strip().split(separator)
            if column in first_line:
                return True
            else:
                print(f'Column {column} not found in file {description_file} using column separator "{separator}"')
                return False

    def get_sample_labels(self, description_file: str, separator: str, label_column: str, bam_files: List[str], sample_labels: Dict[str, str]) -> bool:
        first_column_values:Dict[str, str]={}
        with open(description_file) as input_file:
            first_line = input_file.readline().strip().split(separator)
            column_index = first_line.index(label_column)
            for line in input_file:
                first_column_values[line.strip().split(separator)[0]]=line.strip().split(separator)[column_index]

        outcome=True
        for file in bam_files:
            if Path(file).stem not in first_column_values:
                warnings.warn(f'BAM derived sample name {Path(file).stem} not found in first column of metadata file.')
                outcome=False
            else:
                sample_labels[Path(file).stem]=first_column_values[Path(file).stem]
        return outcome

    def get_ref_bed_model(self, args) -> Tuple[str, str, str]:
        bed_file=expanduser(args.target_regions)
        fasta_file=expanduser(args.reference)
        model_file=expanduser(args.model)
        for file, name in zip( [bed_file, fasta_file, model_file], ["BED", "FASTA", "Model"]):
            if not exists(file):
                print(f'{name} file {file} does not exist.')
                return ("", "", "")
        return (bed_file, fasta_file, model_file)

    def file_exists(self, file_name: str) -> bool:
        if not exists(file_name):
            if isdir(file_name):
                print(f'Directory {file_name} does not exist. Check the address spelling.')
            else:
                print(f'File {file_name} does not exist. Check the address spelling.')
            return False
        else:
            return True

    def make_dir(self, dir_name: str) -> bool:
        if not exists(dir_name):
            mkdir(dir_name)
        else:
            print(f'Directory {dir_name} already exists.')
            return False
        return True
    
    def check_address(self, address:str) -> bool:
        if isdir(address):
            self.file_exists(dirname(address))
        else:
            return self.file_exists(dirname(address))
