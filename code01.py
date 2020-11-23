#!/usr/bin/python3
#coding utf-8

from Bio import Entrez
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
    Give the user the option to continue or not continue with the current dataset

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
    Give advice to users: Greater than 1000 is not recommended
    If number > 1000, let the user choose whether or not to continue

    '''
    if count > 100:
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





def similar_250_seq_for_plot():
    '''
    for PLOT:
    When total_seq_num > 250
    We need to fine the 250 most similar protein sequences
    Plot and show the graph

    '''
    count_ = 0  ### Record the number of minus signs
    count_list = []  ### Store the number of minus signs per sequence
    nameline = []    ### Store the splited first line of FASTA format 
    name_list = []   ### Store the name of esch sequence
    
    ### Get the sequence name and number of minus signs, stored in two lists, respectively
    with open("ali.fa","r") as alignment:
        for line in alignment.readlines():   ### Read each line of the file
            line = line.strip()    ### Remove spaces and tabs at the beginning and end of each line
            if re.match(r'^\>(.*)',line):  ### Find the heading line
                line = line.strip(">")     ### Remove '>' at the begining 
                nameline = line.split()    ### Split the line into a list
                name_list.append(nameline[0])    ### Extract the sequence name
                count_list.append(count_)        ### Obtain the number of '-'
                count_ = 0     ### The counter returns to zero and is ready to count the next '-' number
            else:
                count_ = count_ + line.count("-")
    count_list.append(count_)   ### The '-' number of the last sequence was not added to the list in the loop, so it is added here
    del count_list[0]     ### The first element in count_list[] is invalid, so delete it
    #print(count_list,len(count_list),name_list,len(name_list))
    

    ### Merge the data from two lists into a dictionary   
    ### dictionary count_dict = {seq_name : number of minus signs}
    count_dict = {}
    for i in range(len(name_list)):
        key = name_list[i]
        value = count_list[i]
        count_dict[key] = value

    ### Dictionary sorted by value from small to large
    count_dict_sorted = sorted(count_dict.items(), key=lambda  kv:(kv[1],kv[0]))


    ### Find the most similar one's name for BLAST
    with open ("homo.txt","w") as f:   ### will get "pullseq_1.fa"
        f.write(count_dict_sorted[0][0])
    ### get the most similar protein sequence, save it to "pullseq_1.fa" for BLAST
    subprocess.call("/localdisk/data/BPSM/Assignment2/pullseq -i protein_seq.fa -n homo.txt > pullseq_1.fa", shell=True)


    ### Find the 250 most similar protein sequences' name
    for x in range(250):### 这里可以改数字
        with open ("homo250.txt","a") as fn:   ### will get "pullseq_250.fa"
            fn.write(count_dict_sorted[x][0]+"\n")
    #print(count_dict_sorted[0][0])  得到名字
    # 可以提取对齐文件嘛？？？
    ### get 250 protein sequences, save them to "pullseq_250.fa" for PLOT
    subprocess.call("/localdisk/data/BPSM/Assignment2/pullseq -i ali.fa -n homo250.txt > pullseq_250.fa", shell=True)            
    #subprocess.call("/localdisk/data/BPSM/Assignment2/pullseq -i protein_seq.fa -n homo.txt > pullseq_1.fa", shell=True)

    ### PLOT & show
    print("Plotting...")
    subprocess.call("plotcon -sequence pullseq_250.fa -winsize 5 -graph svg", shell=True)
    print("Show the graph...")
    subprocess.call("eog plotcon.svg", shell=True)




def one_seq_for_BLAST():
    '''
    When total_seq_num <= 250:
    Save the name of the most conservative seq to a txt file
    Get the most similar protein sequence, save it to "pullseq_1.fa" for BLAST
    '''
    count_ = 0  ### Record the number of minus signs
    count_list = []  ### Store the number of minus signs per sequence
    nameline = []    ### Store the splited first line of FASTA format
    name_list = []   ### Store the name of esch sequence

    ### Get the sequence name and number of minus signs, stored in two lists, respectively
    with open("ali.fa","r") as alignment:
        for line in alignment.readlines():   ### Read each line of the file
            line = line.strip()    ### Remove spaces and tabs at the beginning and end of each line
            if re.match(r'^\>(.*)',line):  ### Find the heading line
                line = line.strip(">")     ### Remove '>' at the begining
                nameline = line.split()    ### Split the line into a list
                name_list.append(nameline[0])    ### Extract the sequence name
                count_list.append(count_)        ### Obtain the number of '-'
                count_ = 0     ### The counter returns to zero and is ready to count the next '-' number
            else:
                count_ = count_ + line.count("-")
    count_list.append(count_)   ### The '-' number of the last sequence was not added to the list in the loop, so it is added here
    del count_list[0]     ### The first element in count_list[] is invalid, so delete it

    ### Get index numbers of most similar seq
    index_list = [i for i,x in enumerate(count_list) if x==min(count_list)]
    ### Save the name of the most conservative seq to "homo.txt"
    ### Just take the first sequence
    n = index_list[0]
    with open("homo.txt","w") as f:
        f.write(name_list[n])
    subprocess.call("/localdisk/data/BPSM/Assignment2/pullseq -i protein_seq.fa -n homo.txt > pullseq_1.fa", shell=True)
    

### Extrace similar sequences
#subprocess.call("/localdisk/data/BPSM/Assignment2/pullseq -i protein_seq.fa -n homo.txt > pullseq_1.fa", shell=True)

### Plot conservation of a sequence alignment
### Alignment again
#subprocess.call("clustalo -i pullseq_250.fa -o ali.fa --force", shell=True)



def motifs():
    '''
    Scan a protein sequence with motifs from the PROSITE database
    1. obtain name of each protein sequence
    2. run patmatmotifs : Scan a protein sequence with motifs from the PROSITE database
    '''

    nameline = []
    name_list = []

    with open("ali.fa","r") as alignment:
        for line in alignment.readlines():    ### Read each line of the file
            line = line.strip()    ### Remove spaces and tab at the beginning and end of each line
            if re.match(r'^\>(.*)',line):
                line = line.strip(">")     ### Remove '>' at the begining
                nameline = line.split()
                filename = nameline[0]+".fasta"    ### Extract the sequence name                  
                filename = filename.lower()
                #print(filename)
		
                subprocess.call("patmatmotifs "+ filename,shell=True)        
            else:
                pass



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

    ### PLOT
    ### PLOT conservation of a sequence alignment & show the graph AND find the most conservative sequence for BLAST
    ### No more than 250 sequences can be used for PLOT
    ### If the total number of sequences is greater than 250, we have to pick the 250 that have the highest similarity
    if total <= 250:
        ### When total_seq_num <= 250, call one_seq_for_BLAST()
        ### Save the name of the most conservative seq to a txt file
        ### Get the most similar protein sequence, save it to "pullseq_1.fa" for BLAST
        one_seq_for_BLAST() 

        ### get a similarity plot of aligned sequences (plotcon.svg)
        print("Plotting...")
        subprocess.call("plotcon -sequence ali.fa -winsize 5 -graph svg", shell=True)
        ### show the graph
        print("Show the graph...")
        subprocess.call("eog plotcon.svg", shell=True)
        
    else:
        ### When total_seq_num > 250, call find_similar_250_seq()
        ### Pick the 250 that have the highest similarity, save them to "pullseq_250.fa" for PLOT
        ### Get the most similar protein sequence, save it to "pullseq_1.fa" for BLAST
        find_similar_250_seq()


    ### BLAST
    ### Make a database
    subprocess.call("makeblastdb -in protein_seq.fa -dbtype prot -out selfdb", shell=True)
    ### Run blastp of our pullseq_1.fa against selfdb database
    ### Print out in the form of tables, save to the file blastoutput.out
    subprocess.call("blastp -db selfdb -query pullseq_1.fa -outfmt 7 > blastoutput.out", shell=True)
   

    ### MOTIFS 
    ### Read sequences(protein_seq.fa) and write them to individual files
    subprocess.call("seqretsplit -sequence protein_seq.fa -sformat fasta -osformat fasta",shell=True)
    ### Call motifs()
    ### Obtain name of each protein sequence
    ### Then run patmatmotifs: Scan a protein sequence with motifs from the PROSITE database
    motifs()


__main__()
