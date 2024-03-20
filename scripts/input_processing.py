from typing import List, Tuple
from os.path import expanduser, isdir, exists, dirname
from os import listdir, mkdir

class InputProcessing:

    def __init__(self) -> None:
        pass

    def get_bam_files(self, bam_source: str) -> List[str]:
        bam_source=expanduser(bam_source)
        file_to_classify=[]
        if isdir(bam_source):  #Source is directory with BAM files
            file_to_classify=[bam_source+f for f in listdir(bam_source)  if
                            f.split(".")[-1]=="bam"]
            if len(file_to_classify)==0:
                print(f'No files ending with .bam found in {bam_source}')
        elif bam_source.split(".")[-1]=="bam": #Source is single BAM file
            if exists(bam_source):
                file_to_classify.append(bam_source)
            else:
                print(f'File {bam_source} does not exist. Maybe the location is incorrect?')
        else:
            #check that list of bams contains valid bams
            with open(bam_source) as bam_list:
                for line in bam_source:
                    tentative_bam=expanduser(line.strip())
                    if exists(tentative_bam):
                        file_to_classify.append(tentative_bam)
                    else:
                        print(f'File {tentative_bam} does not exist. Maybe the location is incorrect?')
        return file_to_classify
    
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
        directory=dirname(address)
        return self.file_exists(directory)
