from Bio import Entrez
import time

# Set your email
Entrez.email = "your_email@example.com"

CLASSIFICATION=["Kingdom", "Phylum","Class","Order","Family","Genus","Species"]

def lower_classification(rank):
    for i, value in enumerate(CLASSIFICATION):
        if rank==value.lower():
            if value==CLASSIFICATION[-1]:
                return "Empty"
            else:
                return CLASSIFICATION[i+1]
        

def get_lower_taxa(tax_id, tax_level):
    search = Entrez.esearch(db="taxonomy", term=f"txid{tax_id}[Subtree] AND {tax_level}[Rank] ", retmode="xml", retmax=100)
    record = Entrez.read(search)
    id_list = record["IdList"]

    new_taxids=[]
    if len(id_list)>0:
        # Fetch details for each taxon
        print("id list len: " + str(len(id_list)))
        fetch = Entrez.efetch(db="taxonomy", id=",".join(id_list), retmode="xml")
        data = Entrez.read(fetch)
        for value in data:
            new_taxids.append({value['TaxId']})
    return new_taxids

def get_taxid_rank(tax_id):

    # get the rank of taxID 
    search = Entrez.efetch(db="taxonomy", id=f"{tax_id}", retmode="xml")
    record = Entrez.read(search)
    if len(record)>0:
        return record[0]["Rank"]

starting_taxa=91347
taxa_to_check=[starting_taxa]

while len(taxa_to_check)>0:
    taxa=taxa_to_check[0]
    taxa_to_check.remove(taxa)
    current_level=get_taxid_rank(taxa)
    if lower_classification(current_level)!="Empty":
        lower_taxa=get_lower_taxa(taxa, lower_classification(current_level))
    taxa_to_check=taxa_to_check+lower_taxa
    print(len(taxa_to_check))
    time.sleep(5)
