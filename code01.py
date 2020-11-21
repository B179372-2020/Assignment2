#!/usr/bin/python3
#coding utf-8

from Bio import Entrez
import subprocess
import os


def obtain_search_term():
 
    # user specify the protein family and the taxonomic group
    # store the name entered by user into "final_search_term"
        
    protein_family_name = input("Please enter the name of the protein family  eg. glucose-6-phosphatase:\n")
    taxonomic_group_name = input("Please enter the name of taxonomic group  eg. birds(Aves):\n")
    final_search_term = protein_family_name + " in " + taxonomic_group_name
    print("\nThe term you enter is:",final_search_term, "\n")

    tf = input("Do you want to change protein family or taxonomic group?(please enter 'yes' or 'no'):\n")
    if tf == 'yes':
            protein_family_name = input("Please enter the name of the protein family  eg. glucose-6-phosphatase:\n")
            taxonomic_group_name = input("Please enter the name of taxonomic group  eg. birds(Aves):\n")
            final_search_term = protein_family_name + " in " + taxonomic_group_name
            print("\nThe term you enter is:", final_search_term, "\n")
    if tf == 'no':
            print("\nWe are going to obtain the relevant protein sequence data..."+"\n")

    return protein_family_name,taxonomic_group_name

### protein_family_name,taxonomic_group_name =  obtain_search_term()
#print(protein_family_name)
#print(taxonomic_group_name)

### obtain the relevant protein sequence data, save to a file
# es_com = "esearch -db protein -query \""+taxonomic_group_name+" [organism] AND "+protein_family_name+" [protein]\"|efetch -format fasta > protein_seq.fa"

es_com = "esearch -db protein -query 'aves[organism] AND glucose-6-phosphatase[protein]' |efetch -format fasta > protein_seq.fa"
subprocess.call(es_com,shell=True)

#def _esearch(search_term,email,database):
    # Retrieve the data set specified by the user in the database
    # Return information that meets the criteria
#    Entrez.email = email      ### tell NCBI who you are
#    handle = Entrez.esearch(db = database, term = search_term)
#    search_results = Entrez.read(handle)
#    handle.close()
#    return search_results


#search_term = obtain_search_term()   ### glucose-6-phosphatase in birds!!!!!!!!!!!!!!!!!!!!!!
#print(search_term)
#search_term = "glucose-6-phosphatase in birds"
#email = "964145391@qq.com"     ### ask user for email!!!!!!!!!!!!!!!!!!!
#database = "protein"
#search_results = _esearch(search_term, email, database)

def _efetch(search_results,database):
    acc_list = search_results["IdList"]
    file_name = "protein_seq.fasta"
    with open (file_name,'w') as seq:
        for _id_ in acc_list:
            #print(_id_)
            handle = Entrez.efetch(db=database,id=_id_,rettype="fasta",retmode="text")
            #print(handle.read())
            seq.write(handle.read())  
#_efetch(search_results,database)


### alignment
os.system("clustalo -i protein_seq.fa --maxnumseq 250 -o ali.fa")
#subprocess.call()


