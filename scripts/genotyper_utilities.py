from read_classifier import Classifier, ModelsData
from pickle import load, dump
from typing import Dict
from os.path import expanduser, exists
import argparse
import json

def get_reference_fasta(model_file: str) -> None:
    with open(model_file, "rb") as input_model:
        model_manager: Dict[str, Classifier] =load(input_model)
        for name, model in model_manager.items():
            print(">"+str(model.name))
            print(str(model.nucletoide_seq))

def get_snps(model_file: str) -> None:
    with open(model_file, "rb") as input_model:
        model_manager: Dict[str, Classifier] =load(input_model)
        for name, model in model_manager.items():
            for genotype_snp in model.genotype_snps:
                print(genotype_snp.to_string)

def load_hierarchy(model_file: str, heirarchy_file:str) -> None:
    with open(model_file, "rb") as input_model:
        model_manager: Dict[str, Classifier] =load(input_model)
        # for name, model in model_manager.items():
        #     for genotype_snp in model.genotype_snps:
        #         print(genotype_snp.to_string)                

def test(model_file: str) -> None:
    pass
    # with open(model_file, "rb") as input_model:
    #     model_manager: Dict[str, Classifier] =load(input_model)
    # genotypes={}
    # genotypes_file="/home/lshas17/HandyReadGenotyper/typhi_data_v3/hierarchy.json"
    # with open( expanduser(genotypes_file) ) as gt_data:
    #     genotypes = json.load(gt_data)
    # models_data: ModelsData = ModelsData()
    # #models_data.metadata["hierarchy"] = genotypes
    # for name, model in model_manager.items():
    #     models_data.classifiers[name]=model
    # with open("/home/lshas17/HandyReadGenotyper/models/paratyphi_A_v3.pkl", "wb") as output:
    #     dump(models_data,output)

    # transient_amplicons = [line.strip() for line in open(expanduser("~/HandyReadGenotyper/extra_citrobacter_DO_NOT_DELETE/transient_amplicon.txt"))]
    # print(transient_amplicons)
    # with open(model_file, "rb") as input_model:
    #     model_manager: Dict[str, Classifier] =load(input_model)
    #     for name, model in model_manager.items():
    #         if name in transient_amplicons:
    #             print(name)
    #             model._amplicon_transient=True
    #         else:
    #             model._amplicon_transient=False
    # with open(expanduser("~/HandyReadGenotyper/models/paratyphi_A_v1.pkl"), "wb") as output:
    #     dump(model_manager,output)

def main():

    parser = argparse.ArgumentParser(description='Various utilities for manipulating classification model')
    parser.add_argument('--reference', action="store_true",
                        help='Use to get reference sequences used in classification model. Must specify model (-m)', required=False)
    parser.add_argument('--snps', action="store_true",
                        help='Use to get SNPs in reference sequences used in genotyping. Must specify model (-m)', required=False)
    parser.add_argument('--test', action="store_true",
                        help='Used for testing only', required=False)
    parser.add_argument('-m','--model', metavar='', type=str,
                        help='Pickle (.pkl) file containing pretrained model.', required=False)

    try:
        args = parser.parse_args()
    except:
        parser.print_help()
        exit(0)
    if args.test:
        test(expanduser(args.model))
    if args.reference or args.snps:
        if args.model is None:
            print("To get reference sequence a model file must be specified!")
            exit(0)
        model_file=expanduser(args.model)
        if not exists(model_file):
            print(f'Model file {model_file} does not exist!')
            exit(0)
        if args.reference:
            get_reference_fasta(model_file)
        elif args.snps:
            get_snps(model_file)

if __name__=="__main__":
    main()