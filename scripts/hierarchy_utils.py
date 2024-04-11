from typing import Dict, List
import pandas as pd
from data_classes import Genotype

class HierarchyUtilities:
    """Class representing a hierarchy structure of the target organims.
    The data is supplied as tab-delimited file
    each row must contain (in first column) the target genotype
    and MAY also contain in subsequent column subgenotypes of target genotype
    There is no need to create rows for non-target genotypes
    If the genotype has not subgenotypes, the row would only have target genotype
    ex: 
    4.3.1 4.3.1.1.P1 4.3.1.2 4.3.1.2
    3.1.1
    """
    def __init__(self) -> None:
        self.genotype_hierarchy: Dict[str,Genotype]={}


    def load_hierarchy(self, filename: str) -> Dict[str,Genotype]:
        """Loads hierarchy data
        :param filename: Path to tab delimited file containing genotype hierarchy
        :type filename: str
        """
        self.genotype_hierarchy: Dict[str,Genotype]={}
        with open(filename) as input_file:
            for i, line in enumerate(input_file):
                values=line.strip().split("\t")
                values=[f.strip() for f in values] #remove any accidental spaces
                if values[0]=="":
                    raise ValueError(f'Line {i} in hierarchy file has empty first column: {line}')
                if values[0] in self.genotype_hierarchy:
                    raise ValueError(f'Genotype {values[0]} appears more than once in first column')
                genotype=Genotype(values[0])
                for value in values[1:]:
                    if value in genotype.subgenotypes:
                        raise ValueError(f'Genotype {value} appears more than once on line {i}: {line}')
                    genotype.subgenotypes.append(value)
                self.genotype_hierarchy[ values[0] ]=genotype
        return self.genotype_hierarchy

    _snp_data: pd.DataFrame=pd.DataFrame()
    _column_to_gt: List[str]
    _genotype_snps: pd.DataFrame

