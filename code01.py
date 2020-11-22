#!/usr/bin/python3
#coding utf-8

from Bio import Entrez
import subprocess
import re
import sys

def obtain_search_term():
    '''
    User specify the protein family and the taxonomic group

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


### Assign user's input to variable
protein_family_name, taxonomic_group_name =  obtain_search_term()

### Display the result of esearch on the screen
es_com="esearch -db protein -query \""+taxonomic_group_name+" [organism] AND "+protein_family_name+" [protein]\" > es_result.txt"
subprocess.call(es_com, shell=True)
#es_com = "esearch -db protein -query birds [organism] AND glucose-6-phosphatase [protein]\" > es_result.txt"
#subprocess.call(es_com, shell=True)


#subprocess.call("blastp -db selfdb -query pullseq_pro_seq.fa -outfmt 7 > blastoutput.out",shell=True)


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
                total_seq_number = re.findall(r'\d+',line)    ### Grab the number of sequences in <Count>line 
                for num in total_seq_number: 
                    ### Tell total sequence number to the user
                    print("There are "+num+" sequences in the dataset.") 
                    num = int(num)    ### Convert to int, compare sizes
                    ### 这里可以返回count的数值
    return num

total_seq_number = total_seq_number()
#应该在主程序里写，决定要不要继续，继续就下载，不继续就回到input
if total_seq_number > 100:  ### 这里可以改
    print("We recommend that the total number of sequences is less than 1000.\nThe current number is more than 1000.\n")
    yn=input("Do you want to continue? (please enter 'yes' or 'no'):\n")
    if yn == 'yes':
        print("Continuing...Please wait...Don't touch any buttons...")
        print("We are downloading protein sequences...\n")
    if yn == 'no':
        print("Process ended. Please restart the program")
        sys.exit() 


### back to user input                  
#protein_family_name, taxonomic_group_name =  obtain_search_term()                



### Obtain the relevant protein sequence data, save to a file
es_ef_com = "esearch -db protein -query \""+taxonomic_group_name+" [organism] AND "+protein_family_name+" [protein]\"|efetch -format fasta > protein_seq.fa"
subprocess.call(es_ef_com, shell=True)

### Alignment
subprocess.call("clustalo -i protein_seq.fa -o ali.fa --force", shell=True)


def find_similar_250_seq():
    '''
    When total_seq_num > 250.
    We need to fine the 250 most similar protein sequences.

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
    ### Find the most similar one's name
    with open ("homo.txt","w") as f:   ### will get "pullseq_1.fa"
        f.write(count_dict_sorted[0][0])
    
    ### Find the 250 most similar protein sequences' name
    for x in range(100):### 这里可以改数字
        with open ("homo250.txt","a") as fn:   ### will get "pullseq_250.fa"
            fn.write(count_dict_sorted[x][0]+"\n")
    #print(count_dict_sorted[0][0])  得到名字大写的。。。.fasta
        
    ### Get index numbers of most similar seq
    #index_list = [i for i,x in enumerate(count_list) if x==min(count_list)]
    ### Save the name of the most conservative seq to a txt file 
    #for n in index_list:
    #    with open("homo.txt","a") as f:
    #        f.write(name_list[n]+"\n")
                
    #with open("homo.txt","w+") as f:
    #    f.write(name_list[0])

 
    ### 获得250个序列的名称，从原始文件提取出来250个序列，然后进行对齐，然后画图
    subprocess.call("/localdisk/data/BPSM/Assignment2/pullseq -i ali.fa -n homo250.txt > pullseq_250.fa", shell=True)
    subprocess.call("/localdisk/data/BPSM/Assignment2/pullseq -i protein_seq.fa -n homo.txt > pullseq_1.fa", shell=True)

find_similar_250_seq()
### Extrace similar sequences
#subprocess.call("/localdisk/data/BPSM/Assignment2/pullseq -i protein_seq.fa -n homo.txt > pullseq_1.fa", shell=True)


### make a database
subprocess.call("makeblastdb -in protein_seq.fa -dbtype prot -out selfdb", shell=True)


### Run blastp of our pullseq_pro_seq.fa against the selfdb database, using parameters, saved to the file blastoutput.out
subprocess.call("blastp -db selfdb -query pullseq_1.fa -outfmt 7 > blastoutput.out", shell=True)


### Plot conservation of a sequence alignment
### Alignment again
#subprocess.call("clustalo -i pullseq_250.fa -o ali.fa --force", shell=True)
### get a similarity plot of aligned sequences(plotcon.svg)
subprocess.call("plotcon -sequence pullseq_250.fa -winsize 5 -graph svg", shell=True)

### show the graph
subprocess.call("eog plotcon.svg", shell=True)


### read sequences(protein_seq.fa) and write them to individual files
subprocess.call("seqretsplit -sequence protein_seq.fa -sformat fasta -osformat fasta",shell=True)



### Scan a protein sequence with motifs from the PROSITE database
def motifs():
    '''
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

motifs()
