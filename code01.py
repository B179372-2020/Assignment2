#!/usr/bin/python3
#coding utf-8


import subprocess
import re
import sys



def obtain_search_term():
    '''
    User specify the protein family and the taxonomic group
    Confirm after input

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
            print("\nThe term you enter is:", final_search_term,
                  "\nWe are searching the relevant protein sequence data...\n")
    if tf == 'no':
            print("\nWe are searching the relevant protein sequence data...\n")

    return protein_family_name,taxonomic_group_name





def total_seq_number():
    '''
    Tell the user the total number of sequences

    '''
    with open("es_result.txt","r") as es:
        for line in es.readlines():
            line = line.strip()     ### Removing leading and trailing whitespace and tabs from a string
            line = line.strip("<")    ### Removing leading and trailing '<'
            if re.match(r'C',line):   ### Find the line of <Count> (now beginning with C)
                total_number = re.findall(r'\d+',line)    ### Grab the number of sequences in <Count>line 
                for count in total_number: 
                    ### Tell total sequence number to the user
                    print("There are "+count+" sequences in the dataset.") 
                    count = int(count)    ### Convert to int, compare sizes
    return count


 


def get_appropriate_total_num(count):
    '''
    Determine if the number of sequences is appropriate
    Give the user the option to continue or not continue with the current dataset

    Give advice to users: Greater than 1000 is not recommended
    If number > 1000: let the user choose whether or not to continue
    If number <= 0: number is wrong. Process ended.
    If 0<number<=1000: process continue.
    '''
    if count > 1000: 
        print("We recommend that the total number of sequences is less than 1000.\nThe current number is more than 1000.\n")
        yn=input("Do you want to continue? (please enter 'yes' or 'no'):\n")
        if yn == 'yes':
            print("Continuing...Please wait...Don't touch any buttons...\nWe are downloading protein sequences...\n")
        if yn == 'no':
            print("Process ended. Please restart the program.")
            sys.exit()

    elif count <= 0:
        ### If the number is less than or equal to 0, which is not logical, we will terminate the process.
        print("The current number is wrong.\nProcess ended.\nPlease restart the program and enter other terms to search.")
        sys.exit()
    else:
        ### Tell the user that the programme is continuing
        print("Continuing...Please wait...Don't touch any buttons...\nWe are downloading protein sequences...\n")





def similar_250_seq_and_plot():
    '''
    for PLOT:
    When total_seq_num > 250: we need to pick the 250 most similar protein sequences.
    Plot and show the graph.

    '''
    c = 0  ### Set a counter, count 250 sequences
    with open("blastoutput.out","r") as bout:
        for line in bout.readlines():   ### Read each line of the file
            line = line.strip()    ### Remove spaces and tabs at the beginning and end of each line
            if re.match(r'^\#(.*)',line):  ### Find heading lines
                pass
            else:
                if c <= 249: 
                    c = c+1
                    spline = line.split()    ### Split the line into a list
                    with open ("homo250.txt","a") as homo:
                        homo.write(str(spline[1])+"\n")  ### Extract the sequence name  
                else:
                    break
    
    ### get 250 protein sequences, save them to "pullseq_250.fa" for PLOT
    subprocess.call("/localdisk/data/BPSM/Assignment2/pullseq -i ali.fa -n homo250.txt > pullseq_250.fa", shell=True)            
    
    ### PLOT & show
    print("Plotting...")
    subprocess.call("plotcon -sequence pullseq_250.fa -winsize 5 -graph svg", shell=True)
    print("Show the graph...")
    subprocess.call("eog plotcon.svg", shell=True)





def motifs(count):
    '''
    Scan a protein sequence with motifs from the PROSITE database
    If total_seq_num <= 250: scan each protein sequence with motifs from the PROSITE database
    If total_seq_num > 250: pick 250 sequences that have the highest similarity to scan

    '''

    nameline = []
    name_list = []
    if count <= 250:
        with open("ali.fa","r") as alignment:
            for line in alignment.readlines():    ### Read each line of the file
                line = line.strip()    ### Remove spaces and tab at the beginning and end of each line
                if re.match(r'^\>(.*)',line):
                    line = line.strip(">")     ### Remove '>' at the begining
                    nameline = line.split()
                    filename = nameline[0]+".fasta"    ### Extract the sequence name                  
                    filename = filename.lower()
                    subprocess.call("patmatmotifs "+ filename,shell=True)        
                else:
                    pass
    else:
        ### count > 250
        with open ("homo250.txt","r") as pick:
            for line in pick.readlines():
                line = line.strip()
                filename = line.lower()
                subprocess.call("patmatmotifs "+ filename,shell=True)






def __main__():

    ### Call obtain_search_term(). Get the user's input.
    protein_family_name, taxonomic_group_name =  obtain_search_term()


    ### Display the result of esearch on the screen and save it to es_result.txt
    es = "esearch -db protein -query \""+taxonomic_group_name+" [organism] AND "+protein_family_name+" [protein]\""
    subprocess.call(es, shell=True)
    print("\n")
    es_com="esearch -db protein -query \""+taxonomic_group_name+" [organism] AND "+protein_family_name+" [protein]\">es_result.txt"
    subprocess.call(es_com, shell=True)


    ### Call total_seq_number(), the number shows on the screen
    ### Get the total number of sequences
    total = total_seq_number()


    ### Call get_appropriate_total_num()
    ### Determine if the number of sequences is appropriate
    ### Give advice to users: Greater than 1000 is not recommended
    ### If number > 1000, let the user choose whether or not to continue 
    get_appropriate_total_num(total)

    
    ### User chooses to continue
    ### Obtain the relevant protein sequence data, save to "protein_seq.fa"
    es_ef_com = "esearch -db protein -query \""+taxonomic_group_name+" [organism] AND "+protein_family_name+" [protein]\"|efetch -format fasta > protein_seq.fa"
    subprocess.call(es_ef_com, shell=True)


    ### ALIGN
    ### Align the original sequence. Get the aligned file "ali.fa".
    subprocess.call("clustalo -i protein_seq.fa -o ali.fa --force", shell=True)

    
    ### BLAST
    ### Make a database
    subprocess.call("makeblastdb -in protein_seq.fa -dbtype prot -out selfdb", shell=True)
    print("\n")
    ### Get the sequence that best fits the alignment
    ### Create a consensus sequence from a multiple alignment
    subprocess.call("cons -sequence ali.fa -outseq one_seq.fa", shell=True)
    ### Run blastp against selfdb database
    ### Output in the form of tables, save to the file blastoutput.out
    subprocess.call("blastp -db selfdb -query one_seq.fa -outfmt 7 > blastoutput.out", shell=True)
    print("BLAST finished. Output: blastoutput.out\n")


    ### PLOT
    ### PLOT conservation of a sequence alignment & show the graph 
    ### No more than 250 sequences can be used for PLOT
    ### If the total number of sequences is greater than 250, we have to pick the 250 that have the highest similarity
    if total <= 250: 
        ### When total_seq_num <= 250:
        ### Get a similarity plot of aligned sequences (plotcon.svg)
        print("Plotting...")
        subprocess.call("plotcon -sequence ali.fa -winsize 5 -graph svg", shell=True)
        ### show the graph
        print("Show the graph...")
        subprocess.call("eog plotcon.svg", shell=True)
        
    else:
        ### When total_seq_num > 250:
        ### According to the result of BLAST, pick 250 highest similarity sequences, save them to "pullseq_250.fa" for PLOT
        similar_250_seq_and_plot()
 
        
    ### MOTIFS 
    ### Read sequences(protein_seq.fa) and write them to individual files
    subprocess.call("seqretsplit -sequence protein_seq.fa -sformat fasta -osformat fasta",shell=True)
    ### Call motifs()
    ### Obtain name of each protein sequence
    ### Then run patmatmotifs: Scan a protein sequence with motifs from the PROSITE database
    motifs(total)


__main__()
