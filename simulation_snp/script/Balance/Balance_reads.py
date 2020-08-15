#!/usr/bin/env python
# coding: utf-8

# In[2]:


import gzip
import argparse
import re
import random 

class SNP:
    def __init__(self,chr,position):
        self.chr = chr
        self.pos = int(position)
        self.vcf = [] 
        self.ref = []
    def add_read_id(self, str,read_id):
        if str == "vcf":
            self.vcf.append(read_id)
        elif str == "ref":
            self.ref.append(read_id) 
        
def File2SNPclass(file):
    all_snp = []
    vcf_handler = gzip.open(file,"rb+")
    for line in vcf_handler:
        line = line.decode("utf-8")
        line = line.strip()
        if not line.startswith("#"):
            info = line.split('\t')
            snp = SNP(info[0],int(info[1]))
            all_snp.append(snp)
            
    return all_snp    
    vcf_handler.close()

def Sam2region(line):
    cigar_hash = {'M':1,'I':0,'D':1,'N':1,'=':1,'X':1}
    read_id,flag,chromosome,start_position,tmp,cigar = line.strip().split('\t')[:6]
    start_position = int(start_position)
    factor = 1
    numbers = re.findall('[0-9]+',cigar)
    types = re.findall('[M|I|D|N|S|H|P|=|X]+',cigar)
    types2hash = [cigar_hash[i] for i in types]
    numbers = [int(i) for i in numbers]
    types2hash = [int(i) for i in types2hash]
    
    product = [ a*b for a,b in zip(numbers,types2hash)]
    addition = sum(product) * factor
    end_position = start_position + addition 
    if start_position > end_position:
        start_position, end_position = end_position, start_position 
    
    return read_id,chromosome,start_position,end_position

def addReadID_class(SNP_list,file_handler,read_length,samtype):
    file_handler = open(file_handler,"r+")
    
    #queue to store read id, start postion and chr
    #in order to append it to snp.ref or snp.vcf once it covers two snp. 
    read_id_queue  = []
    start_queue = []
    chr_queue = []
    
    line = file_handler.readline()
    
    for snp in SNP_list:
        for line in file_handler:
            if line.startswith("@"):
                continue 
            else: 
                read_id,chromosome,start,end = Sam2region(line)
                if chromosome == snp.chr and snp.pos - start >= 0: 
                    read_id_queue.append(read_id)
                    chr_queue.append(chromosome)
                    start_queue.append(start)
                    continue 
                elif chromosome == snp.chr  and snp.pos - start < 0:
                    read_id_queue.append(read_id)
                    chr_queue.append(chromosome)
                    start_queue.append(start)
                    
                    remove_index = []
                    for i in range(len(read_id_queue)):
                        if snp.pos - start_queue[i] > read_length:
                            remove_index.append(i)
                        elif snp.pos - start_queue[i]  <= read_length and snp.pos - start_queue[i] > 0 :
                            if samtype == 'ref':
                                if not read_id_queue[i] in snp.ref:
                                    snp.ref.append(read_id_queue[i])
                            elif samtype == 'vcf':
                                if not read_id_queue[i] in snp.vcf:
                                    snp.vcf.append(read_id_queue[i])
                                
                    j = len(remove_index)
                    read_id_queue = read_id_queue[j:]
                    chr_queue = chr_queue[j:]
                    start_queue = start_queue[j:]

                    break
                else: 
                    chr_queue = []
                    start_queue = []
                    read_id_queue = []
                
    file_handler.close()
    return SNP_list



def remove_ele(lst1, lst2):
    number = len(lst1) - len(lst2) 
    eles_rm = random.sample(lst1,number) 
    lst1 = [ele for ele in lst1 if ele not in eles_rm ]
    
    return lst1,eles_rm

def Balance(SNPs):
    id_rm_ref = []
    id_rm_vcf = []
     
    for snp in SNPs:
        for rm_id in id_rm_ref:
            if rm_id in snp.ref:
                snp.ref.remove(rm_id) 
            if rm_id in snp.vcf:
                snp.vcf.remove(rm_id) 
        len_ref = len(snp.ref)
        len_vcf = len(snp.vcf) 
        if len_ref > len_vcf:
            snp.ref,ele_rm = remove_ele(snp.ref, snp.vcf) 
            for ele in ele_rm:
                id_rm_ref.append(ele)
        elif len_ref < len_vcf:
            snp.vcf,ele_rm = remove_ele(snp.vcf,snp.ref)
            for ele in ele_rm:
                id_rm_vcf.append(ele)
        else:
            continue
    return SNPs     


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-refSam',"--ref",help = "address of simulated reference Sam file",
                   default = "/rhome/jzhan413/bigdata/proj/simulation_snp/dataset/simuFq/75bp/ref_genome/test.sam")
    ap.add_argument('-vcfSam',"--vcf",help = "address of simulated vcf Sam file",
                   default = "/rhome/jzhan413/bigdata/proj/simulation_snp/dataset/simuFq/75bp/vcf_genome/test.sam")
    ap.add_argument('-SNPFile',"--snps",help = "address of SNPs file",
                   default = "/rhome/jzhan413/bigdata/proj/simulation_snp/dataset/vcfdataset/vcf5perc/test.vcf.gz")
    ap.add_argument('-readLength',"--length",help = "read length",
                   default = 75 )
    args = ap.parse_args([])
    
    SNP_list = File2SNPclass(args.snps)
    SNP_list  = addReadID_class(SNP_list,args.ref,args.length,'ref')
    SNP_list = addReadID_class(SNP_list,args.vcf,args.length,'vcf')
    SNP_list = Balance(SNP_list)
    
    for snp in SNP_list:
        print(snp.chr,snp.pos,len(snp.ref),len(snp.vcf),snp.ref,snp.vcf)

        
if __name__=="__main__":
    main()


# In[ ]:




