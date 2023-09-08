#Importamos los paquetes necesarios

from Bio import SeqIO
import os
import subprocess

#Creacción de directorios necesarios
if "all_16SrRNA_sequences" not in os.listdir():
    subprocess.run("mkdir all_16SrRNA_sequences", shell=True)
else:
    subprocess.run("rm -r all_16SrRNA_sequences", shell=True)
    subprocess.run("mkdir all_16SrRNA_sequences", shell=True)


#Creación de variables
multiple_rRNA_list = []
no_rRNA_list = []
rRNA_16S_multifasta = ""
rRNA_16S_dict = {}

#Listamos las secuencias contenidas en  "barrnap_results"
for barrnap_G in os.listdir("barrnap_results"):
    rRNA_16S = False
    tmp_RNA_desc = ""   
    tmp_RNA_seq = ""    
    rRNA_16S_num = 0    
    length_list = []    
    with open("barrnap_results/"+barrnap_G, "r") as f:  
        for record in SeqIO.parse(f, "fasta"):  
            if "16S" in record.description:  #si la str '16S' se encuentra en record.description...
                rRNA_16S_indiv = ">"+str(record.description)+"\n"+str(record.seq)+"\n"  

                #Escritura de las secuencias individuales
                with open("all_16SrRNA_sequences/"+barrnap_G.replace(".fa", "")+"_16SrRNA_"+str(rRNA_16S_num)+".fa", "w") as fo:
                    fo.write(rRNA_16S_indiv) 

                #Comprobación longitud de secuencias que previamente se ha confirmado la presencia de 16S
                if len(record.seq)>700: 
                    rRNA_16S = True 
                    length_list.append(len(record.seq))
                    rRNA_16S_num += 1
                    if tmp_RNA_seq == "":
                        tmp_RNA_seq = str(record.seq)
                        tmp_RNA_desc = str(record.description)
                    else:
                        if len(tmp_RNA_seq) < len(record.seq):
                            tmp_RNA_seq = str(record.seq)
                            tmp_RNA_desc = str(record.description)
#si hay mas de 1 rRNA_16S se añade a la correspondiente lista                            
    if rRNA_16S_num >= 1:
        multiple_rRNA_list.append(barrnap_G)

#Si la booleana es True (es decir si len(record.seq)>700...
                
    if rRNA_16S: 
        rRNA_16S_multifasta += ">"+barrnap_G.replace(".fa", "_")+tmp_RNA_desc+"\n"+tmp_RNA_seq+"\n"  
        rRNA_16S_dict[barrnap_G.replace(".fa", "")] = ">"+barrnap_G.replace(".fa", "_")+tmp_RNA_desc+"\n"+tmp_RNA_seq+"\n"
 #En caso de no haber rRNA 16S se añada a una lista diferente
    else:
        print("No 16S rRNA in:", barrnap_G)
        no_rRNA_list.append(barrnap_G.replace(".fa", ""))
#Agregacion de las secuencias al multifasta    
with open("deinococcaceae_longest_16S_rRNA.fna", "w") as f:
    f.write(rRNA_16S_multifasta)

if "Gmulti_16SrRNA_sequences" not in os.listdir():
    subprocess.run("mkdir Gmulti_16SrRNA_sequences", shell=True)
else:
    subprocess.run("rm -r Gmulti_16SrRNA_sequences", shell=True)
    subprocess.run("mkdir Gmulti_16SrRNA_sequences", shell=True)

#Renombramiento de las secuencias para mejor manejo
for rRNA in os.listdir("all_16SrRNA_sequences"): 
    G = "G_"+rRNA.split("_")[1]+".fa"
    if G in multiple_rRNA_list:
        subprocess.run("cp all_16SrRNA_sequences/"+rRNA+" Gmulti_16SrRNA_sequences/"+rRNA, shell=True)




if "Gmulti_16SrRNA_blast_results" not in os.listdir():
    subprocess.run("mkdir Gmulti_16SrRNA_blast_results", shell=True)
else:
    subprocess.run("rm -r Gmulti_16SrRNA_blast_results", shell=True)
    subprocess.run("mkdir Gmulti_16SrRNA_blast_results", shel#Si la booleana es True (es decir si len(record.seq)>700...
l=True)

#Se listas las secuencias obtenidas para compararlas a través de blastn con la base de datos concreta en busca de posibles alineamientos
command_part1 = 'ls Gmulti_16SrRNA_sequences/*.fa | parallel --verbose "blastn -query {} -db 16SrRNA_db/SILVA_138.1_SSURef_NR99_tax_silva -out Gmulti_16SrRNA_blast_results/{/.}.blast.out -max_target_seqs 5 -outfmt'
command_part2 = " '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle'"+str('"')
print(command_part1+command_part2)
subprocess.run(command_part1+command_part2, shell=True) 


genome_16Sgenera = {}
for blast_result in os.listdir("Gmulti_16SrRNA_blast_results"):
    with open("Gmulti_16SrRNA_blast_results/"+blast_result, "r") as f:
        content = f.readlines()    
        genus = content[0].split("/t")[-1].split(";")[-2]
        genome = blast_result.split("_")[1]  
        if genome not in list(genome_16Sgenera.keys()):
            genome_16Sgenera[genome] = [genus]
        else:
            if genus not in genome_16Sgenera[genome]:
                genome_16Sgenera[genome].append(genus)
contamination_list = []
for key in genome_16Sgenera:
    if len(genome_16Sgenera[key]) > 1:
        contamination_list.append("G_"+key)
        print(key, genome_16Sgenera[key])

clean_16rRNA_multifasta = ""

#Añadimos unicamente muestras con rRNA 16S y NO contaminadas
for G in rRNA_16S_dict:
    if G not in contamination_list and G not in no_rRNA_list:
        clean_16rRNA_multifasta += rRNA_16S_dict[G]
with open("Deinococcus_longest_16S_rRNA_clean.fna", "w") as f:
    f.write(clean_16rRNA_multifasta)

#Final stats
print("No 16S rRNA >700 bp list:", len(no_rRNA_list))   
print("Contaminated genomes list:", len(contamination_list), contamination_list)  #Muestras contaminadas

with open("Contaminated_genomes_list.txt", "w") as fo:
    fo.write("\n".join(contamination_list)+"\n")
print("no RNAs en:","\n", no_rRNA_list,"\n")
print(no_rRNA_list)

