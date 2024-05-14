### Checks all the inputs are valid before executing the program
### this avoids program crashing mid-way due to bad inputs.

from os.path import exists
from os import listdir
import warnings
from tqdm import tqdm
from typing import Set
from Bio import SeqIO

class ValidateFiles:

    def __init__(self) -> None:
        self._validated_files: Set[str]=set()
    
    @property
    def validated_files(self) -> Set[str]:
        return self._validated_files


    def validate_many(self, file_list, files_type) -> bool:
        if len(file_list)==0:
            raise ValueError(f'Cannot validate an empty list')
        with tqdm(total=len(file_list)) as progress_meter:
            for file in file_list:
                if files_type=="fasta":
                    self.validate_fasta(file)
                elif files_type=="vcf":
                    self.validate_vcf(file)
                elif files_type=="bed":
                    self.validate_bed(file)
                else:
                    raise ValueError(f'Unknown files type: {files_type}')
                progress_meter.update(1)
        return True

    def validate_bed(self, bed_file_name: str, **kwargs) -> bool:
        """Class for validate input bedfile
        :key min_col_number: minimum required number of column in bed file, default: 3, int
        """        
        if not exists(bed_file_name):
            print(f'File or directory {bed_file_name} does not exist')
            return False
        min_column_count=kwargs.get("min_col_number",3)
        with open(bed_file_name) as bed_file:
            for line_counter, line in enumerate(bed_file):
                if line.find("\t")==-1:
                    print(f'No tab-delimited found in file {bed_file_name} on line {line_counter}:\n {line}')
                    return False
                line_values=line.strip().split("\t")
                if len(line_values)<min_column_count:
                    print(f'Min number of tab-delimited columns is {min_column_count}, but only {len(line_values)} were found on line {line_counter}:\n {line} ')
                    return False
                if not line_values[1].isdigit() or not line_values[2].isdigit():
                    print(f'Expecting numeric values in column 1 and 2, but none numeric values (posibly decimal) values were found in {bed_file_name} on line {line_counter}:\n {line}')
                    return False
                if int(line_values[2])<=int(line_values[1]):
                    print(f'Value in column 2 must be greater than value in column 1. Check {bed_file_name} on line {line_counter}:\n {line}')
                    return False
            if 'line_counter' not in vars():
                print(f'{bed_file_name} is empty')
                return False
        self._validated_files.add(bed_file_name)
        return True

    def validate_fasta(self, fasta_file_name: str) -> bool:
        if not exists(fasta_file_name):
            print(f'File or directory {fasta_file_name} does not exist')
            return False
        with open(fasta_file_name) as fasta_file:
            first_fifty_char=fasta_file.readline()[0:50]
            if len(first_fifty_char)==0:
                print(f'Fasta file {fasta_file_name} is empty')
                return False
            if first_fifty_char[0]!=">":
                print(f'Fasta file must have ">" on first line in {fasta_file_name}\n {first_fifty_char}')
                return False
        self._validated_files.add(fasta_file_name)
        return True
    
    def fasta_has_dashes(self, fasta_file_name: str) -> bool:
        with open(fasta_file_name) as fasta_file:
            for line in fasta_file:
                if line[0]!=">":
                    if line.find("-")>-1:
                        warnings.warn(f'Fasta file {fasta_file_name} has "-". This is likely to cause problems')
                        return True
        return False
    
    def contigs_in_fasta(self, bed_file_name: str, fasta_file_name: str) -> bool:
        ### Assume that bed and fasta files have already been validated
        bed_contigs=set()
        with open(bed_file_name) as bed_file:
            for line in (bed_file):
                bed_contigs.add(line.split("\t")[0])

        for record in SeqIO.parse(fasta_file_name,"fasta"):
            if record.id in bed_contigs:
                bed_contigs.remove(record.id)
            if len(bed_contigs)==0:
                break
        if len(bed_contigs)!=0:
            missing_contigs="\n".join( list(bed_contigs)[0:min(len(bed_contigs),10)] )
            print(f'Some bed file {bed_file_name} contig IDs not found in fasta file \n Ex. {fasta_file_name} (first 10):\n {missing_contigs} ')
            return False
        return True
            
    def contigs_in_vcf(self, bed_file_name: str, vcf_file_name: str) -> bool:
        '''Checks that at least one bedfile contig is present in the VCF file
        Does NOT check if ALL bed file contigs are in the VCF file because some contigs
        may not have SNPs'''
        bed_contigs=set()
        with open(bed_file_name) as bed_file:
            for line in (bed_file):
                bed_contigs.add(line.split("\t")[0])
        if len(bed_contigs)==0:
            warnings.warn(f'Bed file {bed_file_name} is empty')
            return True
        
        vcf_contigs=set()
        with open(vcf_file_name) as vcf_file:
            for line in vcf_file:
                if line[0]!="#":
                    contig_id=line.split("\t")[0]
                    vcf_contigs.add(contig_id)
                    if contig_id in bed_contigs:
                        return True
        if len(vcf_contigs)==0:
            warnings.warn(f'VCF file {vcf_file_name} had no variants (i.e. no lines that do not start with #)')
            return True
        else:
            warnings.warn("\n"+f'None of the contigs in VCF file {vcf_file_name} are present in bedfile {bed_file_name}')
            return False

    def validate_vcf(self, vcfs_dir: str) -> None:
        vcf_files=[file for file in listdir(vcfs_dir) if file.split(".")[-1]=="vcf"]
        if len(vcf_files)==0:
            raise IOError(f'Directory {vcfs_dir} have no VCF files')
        vcf_file_name=vcfs_dir+vcf_files[0]
        if not exists(vcf_file_name):
            raise FileExistsError(f'File or directory {vcf_file_name} does not exist')
        vcf_file= open(vcf_file_name)
        first_two_char=vcf_file.readline()[0:2]
        if first_two_char!="##":
            vcf_file.close()
            raise ValueError(f'First line of VCF file {vcf_file_name} does not start with ##')
        for line in vcf_file:
            if len(line)>=6:
                if line[0:2]=="##":
                    pass #waiting to find line with #CHROM in it
                if line[0:6]=="#CHROM" and len(line.split("\t"))<=9:
                    vcf_file.close()
                    raise ValueError(f'The header VCF line (#CHROM...) should have 9 or more lines, but has fewer in file {vcf_file_name}')
        self._validated_files.add(vcf_file_name)                    
        return None

