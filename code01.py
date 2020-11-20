#!/usr/bin/python3
#coding utf-8

from Bio import Entrez

def obtain_search_term():
        '''
        user specify the protein family and the taxonomic group
        store the name entered by user into "final_search_term"

        '''
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

        return final_search_term


def _esearch(search_term,email,database):
        '''
        Retrieve the data set specified by the user in the database
        Return information that meets the criteria

        '''
        Entrez.email = email
        net = Entrez.esearch(db = database, term = search_term, rettype = "fasta")
        search_results = Entrez.read(net)
        net.close()
        return search_results


search_term = obtain_search_term()
email = "964145391@qq.com"
database = "protein"
search_results = _esearch(search_term, email, database)
print(search_results)
