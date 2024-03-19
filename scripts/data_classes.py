from typing import List, Tuple, Dict
import pandas as pd
import uuid
from io import TextIOWrapper
import copy
from json import load
from os.path import exists, expanduser
from multiprocessing import cpu_count

class SNP:
    def __init__(self, **kwargs) -> None:
        """Constructor

        :param ref_contig_id: The id of the contig which is the reference for this SNP, defaults to empty string
        :type ref_contig_id: str, optional

        :param ref_base: The nucleotide at this position on the reference contig, defaults to empty string
        :type ref_base: str, optional

        :param alt_base: The alternative nucleotide at this position on the reference contig, defaults to empty string
        :type alt_base: str, optional

        :param position: The nucleotide at this position on the reference contig, defaults to empty string
        :type position: int, optional
        
        :param passes_filters: Indicates if SNP passes some filters, defaults to False
        :type passes_filters: bool, optional
        """
        if kwargs.get("ref_contig_id","")!="":
            self._ref_contig_id=kwargs.get("ref_contig_id","")
        if kwargs.get("ref_base","")!="":
            self._ref_base=kwargs.get("ref_base","")
        if kwargs.get("alt_base","")!="":
            self._alt_base=kwargs.get("alt_base","")
        if kwargs.get("position","")!="":
            self._position=kwargs.get("position","")
        self._passes_filters=kwargs.get("passes_filters",False)

    @property
    def coordinate(self) -> Tuple[str, int]:
        return (self._ref_contig_id, self._position)

    @property
    def ref_contig_id(self) -> str:
        return self._ref_contig_id

    @ref_contig_id.setter
    def ref_contig_id(self, value: str):
        self._ref_contig_id = value

    @property
    def position(self) -> int:
        return self._position

    @position.setter
    def position(self, value: int):
        self._position = int(value)

    @property
    def ref_base(self) -> str:
        return self._ref_base

    @ref_base.setter
    def ref_base(self, value: str):
        self._ref_base = value

    @property
    def alt_base(self) -> str:
        return self._alt_base

    @alt_base.setter
    def alt_base(self, value: str):
        self._alt_base = value

    @property
    def sensitivity(self) -> float:
        return self._sensitivity

    @sensitivity.setter
    def sensitivity(self, value: float):
        self._sensitivity = float(value)

    @property
    def specificity(self) -> float:
        return self._specificity

    @specificity.setter
    def specificity(self, value: float):
        self._specificity = float(value)

    @property
    def passes_filters(self) -> bool:
        return self._passes_filters

    @passes_filters.setter
    def passes_filters(self, value: bool):
        self._passes_filters = value

    @property
    def is_genotype_snp(self) -> bool:
        return self._is_genotype_snp

    @is_genotype_snp.setter
    def is_genotype_snp(self, value: bool):
        self._is_genotype_snp = value
        self._is_species_snp=False
    
    @property
    def is_species_snp(self) -> bool:
        return self._is_species_snp

    @is_species_snp.setter
    def is_species_snp(self, value: bool):
        self._is_species_snp = value
        self._is_genotype_snp=False

    def to_file(self, file_handle: TextIOWrapper, **kwargs):
        """Constructor
        :param file_handle: File handle to which the SNP will be written
        :type file_handle: TextIOWrapper

        :param sep: Column separator to use, defaults to tab
        :type sep: str, optional
        """        
        sep=kwargs.get("sep","\t")
        name=f'{self._ref_base}/{self.alt_base}'
        if self._is_species_snp:
            name=name+"/species"
        elif self._is_genotype_snp:
            name=name+"/genotype"

        file_handle.write(sep.join( [str(f) for f in [self._ref_contig_id, self.position, self.position+1, name] ] )+"\n")

    def __eq__(self, other) -> bool:
        return self.coordinate==other.coordinate and self.alt_base==other.alt_base
    
    def __lt__(self, other):
        if self.ref_contig_id!=other.ref_contig_id:
            return sorted([self.ref_contig_id, other.ref_contig_id])[0]==self.ref_contig_id #sort by IDs alphabetically
        else:
            return self.position<other.position #sort by position if IDs are the same


    def __hash__(self):
        return hash( (self.ref_contig_id, self.position, self.alt_base) )

    def copy(self):
        """Creates deep copy of the SNP instance
        """
        return copy.deepcopy(self)

class Sample:
    """Represents a VCF sample and contains SNPs associated with it

    SNPs have to be loaded separately to avoid creating multitude of objects
    """
    def __init__(self, name: str, vcf_file: str) -> None:
        self._name=name
        self._snps=[]
        self._vcf_file=vcf_file
        self._uuid=str(uuid.uuid4())
        self._genotype=""

    
    @property
    def genotype(self) -> str:
        return self._genotype

    @genotype.setter
    def genotype(self, value: str):
        self._genotype = value

    @property
    def vcf_file(self) -> str:
        return self._vcf_file

    @property
    def snps(self) -> List[SNP]:
        return self._snps

    @property
    def name(self) -> str:
        return self._name
    
    @property
    def id(self) -> str:
        return self._uuid


class ReferenceSequence:
    def __init__(self, contig_id, seq_start, seq_end, sequence) -> None:
        self._refseq_id=contig_id
        self._ref_start=seq_start
        self._ref_end=seq_end
        self._sequence=sequence

    @classmethod
    def from_bed_line(cls, bed_line:str, ref_fasta_file: str):
        """Constructor using bedfile lines. 
        :param bed_line: String from bedfile, if the line has fourth column, this will be included in amplicon name
        :type bed_file: str

        :param ref_fasta_file: Path to fasta file on which the amplicon is based
        :type ref_fasta_file: str
        """
        bed_line_values = bed_line.strip().split("\t")
        contig_id, seq_start, seq_end=bed_line_values[0:3]
        seq_start=int(seq_start)
        seq_end=int(seq_end)
        with open(ref_fasta_file) as fasta_input:
            for record in SeqIO.parse(fasta_input,"fasta"):
                if record.id==contig_id:
                    if seq_end>len(str(record.seq)):
                        raise ValueError(f'Contig {contig_id} is shorter, {len(str(record.seq))}nt, than end position in the bedfile: {seq_end}')
                    new_refseq=cls(contig_id, seq_start, seq_end, str(record.seq[seq_start:seq_end]))
                    return new_refseq
        raise ValueError(f'Contig {contig_id} is not found in fasta file: {ref_fasta_file}')
    
    @property
    def refseq_id(self) -> str:
        return self._refseq_id

    @property
    def ref_start(self) -> int:
        return self._ref_start

    @property
    def ref_end(self) -> int:
        return self._ref_end
    
    @property
    def sequence(self) -> str:
        return self._sequence


class Amplicon:
    def __init__(self, name: str, seq: str) -> None:
        self._name: str=name
        self._seq: str=seq
        self._snps:List[SNP]=[]
        self._left_flanking_id=""
        self._right_flanking_id=""
        self._has_homologues=False
        self._uuid=str(uuid.uuid4())
        self._has_reference=False #indicates if the amplicon has associated reference sequence

    @classmethod
    def from_bed_line(cls, bed_line:str, ref_fasta_file: str):
        """Constructor using bedfile lines. 
        :param bed_line: String from bedfile, if the line has fourth column, this will be included in amplicon name
        :type bed_file: str

        :param ref_fasta_file: Path to fasta file on which the amplicon is based
        :type ref_fasta_file: str

        """
        bed_line_values = bed_line.strip().split("\t")
        ampl_chr, ampl_start, ampl_end=bed_line_values[0:3]
        if len(bed_line_values)>=4:
            name='_'.join( [str(f) for f in bed_line_values[0:4] ] )
        else:
            name='_'.join( [str(f) for f in bed_line_values[0:3] ] )

        new_refseq=ReferenceSequence.from_bed_line(bed_line,ref_fasta_file)
        new_amplicon=cls(name, new_refseq.sequence)
        new_amplicon.ref_seq=new_refseq
        return new_amplicon


    @property
    def ref_seq(self) -> ReferenceSequence:
        if not self._has_reference:
            raise ValueError(f'Amplicon {self._name} has no reference')
        return self._ref_seq

    @ref_seq.setter
    def ref_seq(self, value: ReferenceSequence):
        self._has_reference=True
        self._ref_seq=value

    @property
    def has_reference(self) -> bool:
        return self._has_reference

    @property
    def ref_contig(self) -> str:
        if not self._has_reference:
            raise ValueError(f'Amplicon {self._name} has no reference')
        return self._ref_seq.refseq_id

    @property
    def seq(self) -> str:
        if self._has_reference:
            return self.ref_seq.sequence
        else:
            return self._seq

    # @ref_contig.setter
    # def ref_contig(self, value: str):
    #     self._ref_contig = value

    # @property
    # def ref_start(self) -> int:
    #     if not self._has_reference:
    #         raise ValueError(f'Amplicon {self._name} has no reference')
    #     return self._ref_start

    # @ref_start.setter
    # def ref_start(self, value: int):
    #     self._ref_start = int(value)

    # @property
    # def ref_end(self) -> int:
    #     if not self._has_reference:
    #         raise ValueError(f'Amplicon {self._name} has no reference')
    #     return self._ref_end

    # @ref_end.setter
    # def ref_end(self, value: int):
    #     self._ref_end = int(value)

    @property
    def left_flanking_id(self) -> str:
        return self._left_flanking_id

    @left_flanking_id.setter
    def left_flanking_id(self, value: str):
        self._left_flanking_id = value

    @property
    def right_flanking_id(self) -> str:
        return self._right_flanking_id

    @right_flanking_id.setter
    def right_flanking_id(self, value: str):
        self._right_flanking_id = value

    @property
    def has_flanking(self) -> bool:
        return self._left_flanking_id!="" and self._right_flanking_id!=""

    @property
    def name(self) -> str:
        return self._name

    @property
    def id(self) -> str:
        return self._uuid

    @property
    def len(self) -> int:
        return len(self.seq)

    @property
    def has_homologues(self) -> bool:
        return self._has_homologues

    @has_homologues.setter
    def has_homologues(self, value: bool):
        self._has_homologues = value

    def snp_in_amplicon(self, snp:SNP) -> bool:
        if snp.ref_contig_id==self.ref_contig and \
        snp.position>=self.ref_seq.ref_start and snp.position<=self.ref_seq.ref_end:
            return True
        else:
            return False
        
    def coord_in_amplicon(self, coordinates:Tuple[str, int]) -> bool:
        if coordinates[0]==self.ref_contig and \
        coordinates[1]>=self.ref_seq.ref_start and coordinates[1]<=self.ref_seq.ref_end:
            return True
        else:
            return False

    @property
    def snps(self) -> List[SNP]:
        return self._snps

    @snps.setter
    def snps(self, value: List[SNP]):
        self._snps = value

    def __hash__(self):
        return hash(self.id)


class Genotype:

    def __init__(self, name: str) -> None:
        self._name=name
        self._subgenotypes:List[str]=[name] #every genotype has itself as subgenotypes
        self._alleles: Dict[SNP, str]={}
        self._allele_depths: Dict[SNP, int]={}
        self._amplicons: List[Amplicon]=[]


    @property
    def name(self) -> str:
        return self._name
    
    @property
    def subgenotypes(self) -> List[str]:
        return self._subgenotypes

    @subgenotypes.setter
    def subgenotypes(self, value: List[str]):
        self._subgenotypes = value

    @property
    def amplicons(self) -> List[Amplicon]:
        return self._amplicons

    @amplicons.setter
    def amplicons(self, value: List[Amplicon]):
        self._amplicons = value

    @property
    def defining_snps(self) -> List[SNP]:
        return list(self._alleles.keys())

    def get_genotype_allele(self, snp: SNP) -> str:
        if snp not in self._alleles:
            raise ValueError(f'SNP with coordinates {snp.coordinate} is not present among snps of genotype {self._name}')
        return self._alleles[snp]

    def get_genotype_allele_depth(self, snp: SNP) -> str:
        if snp not in self._allele_depths:
            raise ValueError(f'SNP with coordinates {snp.coordinate} is not present among snps of genotype {self._name}')
        return self._allele_depths[snp]

    def add_genotype_allele(self, snp: SNP, allele: str, depth: int):
        self._alleles[snp]=allele
        self._allele_depths[snp]=depth

    @property
    def defining_snp_coordinates(self) -> List[Tuple[str,int]]:
        return [f.coordinate for f in self.defining_snps if f.passes_filters]
    
class BlastResult:
    def __init__(self) -> None:
        pass

    @classmethod
    def from_blast_line(cls, blast_line:str):
        """Constructor using blast result output. 
        :param blast_line: String from blast output, blast must  use "-outfmt "6 delim=  qseqid qstart qend sseqid sstart send pident evalue qseq"
        :type blast_line: str

        """
        values=blast_line.split("\t")
        new_result=cls()
        new_result.qseqid=values[0]
        new_result.qstart=int(values[1])
        new_result.qend=int(values[2])
        new_result.sseqid=values[3]
        new_result.sstart=int(values[4])
        new_result.send=int(values[5])
        new_result.pident=float(values[6])
        new_result.evalue=float(values[7])
        new_result.qseq=values[8]
        return new_result 

    @property
    def q_hit_len(self) -> int:
        return len(self._qseq.replace("-",""))

    @property
    def qseq(self) -> str:
        return self._qseq

    @qseq.setter
    def qseq(self, value: str):
        self._qseq = value

    @property
    def evalue(self) -> float:
        return self._evalue

    @evalue.setter
    def evalue(self, value: float):
        self._evalue = float(value)

    @property
    def pident(self) -> float:
        return self._pident

    @pident.setter
    def pident(self, value: float):
        self._pident = float(value)

    @property
    def send(self) -> int:
        return self._send

    @send.setter
    def send(self, value: int):
        self._send = int(value)

    @property
    def sstart(self) -> int:
        return self._sstart

    @sstart.setter
    def sstart(self, value: int):
        self._sstart = int(value)

    @property
    def sseqid(self) -> str:
        return self._sseqid 

    @sseqid.setter
    def sseqid(self, value: str):
        self._sseqid = value

    @property
    def qseqid(self) -> str:
        return self._qseqid

    @qseqid.setter
    def qseqid(self, value: str):
        self._qseqid = value

    @property
    def qstart(self) -> int:
        return self._qstart

    @qstart.setter
    def qstart(self, value: int):
        self._qstart = int(value)

    @property
    def qend(self) -> int:
        return self._qend
    
    @qend.setter
    def qend(self, value: int):
        self._qend = int(value)

    @property
    def query_file_name(self) -> str:
        return self._query_file_name

    @query_file_name.setter
    def query_file_name(self, value: str):
        self._query_file_name = value

    @property
    def value(self) -> str:
        return f'{self._qseqid} {str(self._qstart)} {str(self._qend)} {self._sseqid} {str(self._sstart)} {str(self._send)}'
    
    @property
    def is_flipped(self) -> bool:
        """Indicates if query and subject sequences are same strand i.e. align --> -->/ <-- <-- (True) or --> <-- / <-- / --> (False)
        """
        return self.send<self.sstart
    
    def coordinates_match(self, blast_result) -> bool:
        for value in ["qstart", "qend", "sstart","send","pident","evalue","qseq"]:
            if getattr(self,value)!=getattr(blast_result,value):
                return False
        return True
        
class Primer:
    def __init__(self, seq: str, g_c: float, t_m: float, is_reverse) -> None:
        self._t_m=t_m
        self._seq=seq
        self._g_c=g_c
        self._ref_start=-1
        self._is_reverse=is_reverse

    @property
    def t_m(self) -> float:
        return self._t_m    

    @t_m.setter
    def t_m(self, value: float):
        self._t_m = float(value)

    @property
    def seq(self) -> str:
        return self._seq

    @seq.setter
    def seq(self, value: str):
        self._seq = value

    @property
    def g_c(self) -> float:
        return self._g_c

    @g_c.setter
    def g_c(self, value: float):
        self._g_c = float(value)

    @property
    def is_reverse(self) -> bool:
        return self._is_reverse

    @property
    def ref_start(self) -> int:
        return self._ref_start

    @ref_start.setter
    def ref_start(self, value: int):
        self._ref_start = int(value)

    @property
    def ref_end(self) -> int:
        return self._ref_start+len(self.seq)
    
    @property
    def length(self) -> int:
        return len(self.seq)

    def __hash__(self):
        return hash( (self.seq, self.ref_start, self.ref_end) )

class PrimerPair:
    def __init__(self, name: str, forward: Primer, reverse: Primer) -> None:
        self._name=name
        self._uuid=str(uuid.uuid4())
        self._forward=forward
        self._reverse=reverse
        self.ref_contig="Unknown"
        self.targets: List[str]=[]

    @property
    def forward(self) -> Primer:
        return self._forward

    @property
    def reverse(self) -> Primer:
        return self._reverse

    @property
    def ref_contig(self) -> str:
        return self._ref_contig

    @ref_contig.setter
    def ref_contig(self, value: str):
        self._ref_contig = value


    @property #Lower penalty is better
    def penalty(self) -> int:
        return self._penalty
    
    @penalty.setter
    def penalty(self, value: float):
        self._penalty = float(value)

    @property
    def targets(self) -> List[SNP]:
        return self._targets
    
    @targets.setter
    def targets(self, value: List[SNP]):
        self._targets = [f for f in value]

    @property
    def primers(self) -> List[Primer]:
        return [self._forward, self._reverse]

    @property
    def name(self) -> str:
        return self._name

    @property
    def uuid(self) -> str:
        return self._uuid

    @property
    def value(self) -> str:
        return '\t'.join([str(f) for f in [self.name,self.forward.ref_start, self.reverse.ref_end,
                          self.length, "NA", "NA", self.forward.seq, self.forward.t_m,
                          self.reverse.seq, self.reverse.t_m]] )

    @property
    def length(self) -> int:
        return self._reverse.ref_end - self._forward.ref_start

    def seq_in_pair(self, seq: str) -> bool:
        return self._forward.seq==seq or self._reverse.seq==seq
    
    def __eq__(self, __value: object) -> bool:
        return self.uuid==__value.uuid
    
    def to_string(self) -> str:
        return '\t'.join([str(f) for f in [self.name,
                                    '{0:.2f}'.format(self.penalty),
                                    self.ref_contig,
                                    self.forward.ref_start,
                                    self.reverse.ref_end,
                                    self.reverse.ref_end-self.forward.ref_start,
                                    self.forward.seq,
                                    '{0:.2f}'.format(self.forward.t_m),
                                    '{0:.2f}'.format(self.forward.g_c),
                                    self.reverse.seq,
                                    '{0:.2f}'.format(self.reverse.t_m),
                                    '{0:.2f}'.format(self.reverse.g_c)]])

class Genotypes:
    def __init__(self, **kwargs) -> None:
        """Constructor

        :param genotypes: List of Genotype objects, defaults to empty list
        :type genotypes: List[Genotype], optional
        """
        self._genotypes=kwargs.get("genotypes",[])

    @property
    def genotypes(self) -> List[Genotype]:
        return self._genotypes

    @genotypes.setter
    def genotypes(self, value: List[Genotype]):
        self._genotypes = list(value)

    def get_genotype(self, genotype_name: str) -> Genotype:
        for genotype in self._genotypes:
            if genotype.name==genotype_name:
                return genotype


    def all_snps_coord_sorted(self) -> List[Tuple[str, int]]:
        """Returns a list of unique contig + position pairs 
        sorted by contig and position
        """
        unique_contig_pos: List[Tuple[str, int]]=list(set([snp.coordinate for snps in self.genotypes for snp in snps.defining_snps]))
        unique_contig_pos=sorted(unique_contig_pos, key=lambda x: x)
        return unique_contig_pos

    def genotypes_to_snp_matrix(self) -> pd.DataFrame:
        """Converts list of genotypes into a matrix with each row a position on reference
        and each column a genotype. Two extra columns are Contig and Position
        The values indicate if the SNP passes some filter for that column (usually genotype)
        """
        if len(self._genotypes)==0:
            raise ValueError("The object has no genotypes in it.")
        unique_contig_pos: List[Tuple[str, int]]=self.all_snps_coord_sorted()
        gts: List[str]=[f.name for f in self.genotypes]
        output_df: pd.DataFrame=pd.DataFrame(columns=["Contig","Position"]+gts, index=range(0,len(unique_contig_pos))).fillna(False)
        output_df["Contig"]=[f[0] for f in unique_contig_pos]
        output_df["Position"]=[f[1] for f in unique_contig_pos]
        for gt in self.genotypes:
            for snp in gt.defining_snps:
                index=unique_contig_pos.index((snp.ref_contig_id, snp.position))
                output_df.loc[index, gt.name]=snp.passes_filters
        return output_df



    cpu_threads=1
    flank_len_to_check=-1 #this is intentional to avoid hiding this parameters,
    max_blast_length_diff=-1 #they must be defined in config file
    min_blast_identity=-1
    use_negative_genomes_subdir=False
    output_dir=""
    sensitivity_limit: float=-1.0
    specificity_limit: float=-1.0
    min_amplicon_length=200
    def __init__(self, file_name: str):
        file_name=expanduser(file_name)
        try:
            with open(file_name) as config_file:
                self._config_data = load(config_file)
                InputConfiguration.cpu_threads=min(  max(cpu_count()-1,1) , self._config_data["max_cpus"] )
                InputConfiguration.flank_len_to_check=self._config_data["analysis_parameters"]["flank_len_to_check"]
                InputConfiguration.max_blast_length_diff=self._config_data["analysis_parameters"]["max_blast_length_diff"]
                InputConfiguration.min_blast_identity=self._config_data["analysis_parameters"]["min_blast_identity"]
                InputConfiguration.use_negative_genomes_subdir=str.lower(self._config_data["input_directories"]["use_negative_genomes_subdir"])=="true"
                InputConfiguration.output_dir=expanduser(self._config_data["output_files"]["output_dir"])
                InputConfiguration.specificity_limit=self._config_data["analysis_parameters"]["snp_specificity"]/100
                InputConfiguration.sensitivity_limit=self._config_data["analysis_parameters"]["snp_sensitivity"]/100
                InputConfiguration.min_amplicon_length=self._config_data["analysis_parameters"]["min_amplicon_length"]
        except IOError as error:
            if not exists(config_file):
                raise IOError(f'Config file {config_file} does not exist') from error
            else:
                raise IOError(f'Error loading file {config_file}. It exits, but cannot be processed.') from error
        
    @property
    def config_data(self) -> Dict:
        return self._config_data

    @property
    def root_dir(self) -> str:
        return expanduser(self._config_data["root_dir"])

    @property
    def name_stubs(self) -> str:
        return self._config_data["name_stubs"]

    @property
    def reference_fasta(self) -> str:
        return expanduser(self._config_data["input_files"]["reference_fasta"])
    
    @property
    def repeats_bed_file(self) -> str:
        return expanduser(self._config_data["input_files"]["repeats_bed_file"])
    
    @property
    def hierarchy_file(self) -> str:
        return expanduser(self._config_data["input_files"]["hierarchy_file"])

    @property
    def meta_data_file(self) -> str:
        return expanduser(self._config_data["input_files"]["meta_data_file"])
    
    @property
    def vcf_dir(self) -> str:
        return expanduser(self._config_data["input_directories"]["vcf_dir"])

    @property
    def negative_genomes(self) -> str:
        return expanduser(self._config_data["input_directories"]["negative_genomes"])
    
    @property
    def temp_blast_db(self) -> str:
        return expanduser(self._config_data["input_directories"]["temp_blast_db"])
    
    @property
    def metadata_delim(self) -> str:
        return self._config_data["metadata_parameters"]["delimiter"]

    @property
    def genotype_column(self) -> str:
        return self._config_data["metadata_parameters"]["genotype_column"]
    

    @property
    def primer_opt_size(self) -> str:
        return self._config_data["primers_parameters"]["PRIMER_OPT_SIZE"]

    @property
    def primer_opt_tm(self) -> str:
        return self._config_data["primers_parameters"]["PRIMER_OPT_TM"]

    @property
    def primer_min_tm(self) -> str:
        return self._config_data["primers_parameters"]["PRIMER_MIN_TM"]

    @property
    def primer_max_tm(self) -> str:
        return self._config_data["primers_parameters"]["PRIMER_MAX_TM"]

       
    @property
    def gts_with_few_snps(self) -> str:
        return self._config_data["analysis_parameters"]["gts_with_few_snps"]

    @property
    def genotype_snps(self) -> str:
        return self.output_dir+self._config_data["output_files"]["genotype_snps"]
    
    @property
    def snps_bed(self) -> str:
        return self.output_dir+self._config_data["output_files"]["snps_bed"]
    
    @property
    def genotypes_data(self) -> str:
        return self.output_dir+self._config_data["output_files"]["genotypes_data"]
    
    @property
    def species_data(self) -> str:
        return self.output_dir+self._config_data["output_files"]["species_data"]
    
    @property
    def multi_gt_intervals(self) -> str:
        return self.output_dir+self._config_data["output_files"]["multi_gt_intervals"]
    
    @property
    def msa_dir(self) -> str:
        return self.output_dir+self._config_data["output_files"]["msa_dir"]
    
    @property
    def snps_vcf(self) -> str:
        return self.output_dir+self._config_data["output_files"]["snps_vcf"]