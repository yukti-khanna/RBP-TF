import os,sys
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random
import math
#def remove_repeats(indict):
   # for key, value in indict.items():
   #     indict[key]=list(set(value))
    #return

def map_uniprot_gene(mapfile):
       unid_gene={}
       gene_unid={}
       f1=open(mapfile,'r')
       for line in f1:
           stripped=line.strip()
           splitted=stripped.split()
           unid_gene.setdefault(splitted[0],splitted[1])
           gene_unid.setdefault(splitted[1],splitted[0])
        
       f1.close()
       return unid_gene, gene_unid

def get_random_unids(mapfile, length):
    unid_list=[]
    ran_list=[]
    f1=open(mapfile,'r')
    for line in f1:
        stripped=line.strip()
        splitted=stripped.split()
        unid_list.append(splitted[0])
    ran_list=random.sample(unid_list, length)
    ran_list2=random.sample(unid_list, length)

    #for item in unid_list[randrange(0,length):randrange(0,length)]:
       # ran_list.append(item)
        #set seed
    return ran_list, ran_list2

def get_uni_ids(input_file):
    uniprot_list=[]
    f1=open(input_file, "r")
    for line in f1:
        stripped=line.strip()
        uniprot_list.append(stripped)
    f1.close()
    length=len(uniprot_list)
    return uniprot_list, length
        
def gene_to_unid(input_gene_file, gene_unid_dict):
    uniprot_list=[]
    gene_list=[]
    f1=open(input_gene_file, "r")
    for line in f1:
        stripped=line.strip()
        gene_list.append(stripped)
       
    f1.close()
    for gene in gene_list:
        uni_id=gene_unid_dict[gene]
        uniprot_list.append(uni_id)
    length=len(uniprot_list)    

    return uniprot_list, gene_list, length

def get_proteins_connectivity(conn_dict, unid_gene_dict,query_unids):

    #def remove_repeats(indict):
     #   for key, value in indict.items():
    #        indict[key]=list(set(value))
     #   return

    query_dict_unid={}
    query_dict_gene={}
    for query in query_unids:
        query_dict_unid.setdefault(query,conn_dict[query])
        #print (query_dict_unid)
        gene_query=unid_gene_dict[query]
        gene_conns=[unid_gene_dict[conn] for conn in conn_dict[query]]
        query_dict_gene.setdefault(gene_query,gene_conns)

   # query_dict_unid=remove_repeats(query_dict_unid)
    #query_dict_gene=remove_repeats(query_dict_gene)

    return query_dict_unid,query_dict_gene

   

def get_total_connectivity(conn_file, unid_gene_dict):
    
    def remove_repeats(indict):
        for key, value in indict.items():
            indict[key]=list(set(value))
        return indict

    conn_dict= {}
    conn_gene_dict= {}
    f1=open(conn_file, "r")
    for line in f1:
        stripped= line.strip()
        splitted= stripped.split()
        conn_dict.setdefault(splitted[0],[]).append(splitted[1])
        conn_dict.setdefault(splitted[1],[]).append(splitted[0])

        gene1=unid_gene_dict[splitted[0]]
        gene2=unid_gene_dict[splitted[1]]
        conn_gene_dict.setdefault(gene1,[]).append(gene2)
        conn_gene_dict.setdefault(gene2,[]).append(gene1)

    f1.close()

    conn_dict=remove_repeats(conn_dict)
    conn_gene_dict=remove_repeats(conn_gene_dict)

    return conn_dict,conn_gene_dict

def protein1_protein2_connectivity(conn_dict,unid_gene_dict,query1_unids, query2_unids):
    
    pp_dict_unid={}
    pp_dict_gene={}
    for query in query1_unids:
        p2=conn_dict[query]
        common=set(p2).intersection(query2_unids)
        #print (common)
        #print (p2)
        #if len(common)>1:
        pp_dict_unid.setdefault(query,common)
        gene_query=unid_gene_dict[query]
        gene_conns=[unid_gene_dict[conn] for conn in common]
        pp_dict_gene.setdefault(gene_query,gene_conns)
       
    return pp_dict_unid, pp_dict_gene

def get_len_array(unid_gene_dict, conn_dict):
    len_dict={}
    length=[]
    uniprot_list=[]
    for key in conn_dict.keys():
        uniprot_list.append(key)       
    for query in uniprot_list:
        #print (query)
        l1=len(conn_dict[query])
        length.append(l1)
        l=len(length)
        gene_query=unid_gene_dict[query]
        len_dict.setdefault(gene_query,l1)
        data=list(len_dict.items())
        len_array=np.array(data)
        len_list=np.array(length)
    #print (len(length))
    return len_array, len_list, l


if __name__=="__main__":
    list1='RBP_list.txt'
    list2='TF_list.txt'
    list3='RGG_list.txt'
    interactome_file=r"human_annotated_PPIs_unids_w_methods_evidence_nopred.txt"
    uniprot_gene_file=r'UNIPROTIDS_GENENAME.txt'
    outfile="out.txt"
    unid_gene_map,gene_unid_map=map_uniprot_gene(uniprot_gene_file)
    #print (unid_gene_map)
    #print (gene_unid_map)
    totcdict_un, totcdict_gene=get_total_connectivity(interactome_file,unid_gene_map)
    #totcdict_gene, totcdict_un=get_total_connectivity(interactome_file,gene_unid_map)
    #print (totcdict_un)
    query_list1, stupid_list1, len_list1=gene_to_unid(list1,gene_unid_map)
    query_list2, stupid_list2, len_list2=gene_to_unid(list2,gene_unid_map)
    query_list3, len_list3=get_uni_ids(list3)
    #print (query_list2)
    #print (query_list1)
    cdict_un,cdict=get_proteins_connectivity(totcdict_un,unid_gene_map,query_list1)
    #print(cdict)
    #print(cdict['SMN2'])
    pp_dict_un,pp_dict=protein1_protein2_connectivity(totcdict_un, unid_gene_map,query_list1,query_list2)
    
    #print (pp_dict)
    conn_array, conn_list, conn_len=get_len_array(unid_gene_map, cdict_un)
    pp_array, pp_list, len1=get_len_array(unid_gene_map, pp_dict_un)
    ratio1=np.divide(pp_list,conn_list)

    np.histogram(ratio1, bins=1)
    #plt.hist(ratio1, bins=200)
   # plt.savefig("pp_his_tf.pdf", dpi=300)

    random, random2=get_random_unids(uniprot_gene_file, len_list1)
    #print (random)
    #print (len_list1)
    #print (len(random))
    ran_dict_un,ran_dict=protein1_protein2_connectivity(totcdict_un, unid_gene_map,query_list1,random)
    ran_array, ran_list, len2=get_len_array(unid_gene_map, ran_dict_un)
    ratio2=np.divide(ran_list,conn_list)

    #plt.hist(ratio2, bins=200)
    #plt.savefig("random1_his_tf.pdf", dpi=300)
   
    ran2_dict_un,ran2_dict=protein1_protein2_connectivity(totcdict_un, unid_gene_map,query_list1,random2)
    ran2_array, ran2_list, len3=get_len_array(unid_gene_map, ran2_dict_un)
    ratio3=np.divide(ran2_list,conn_list)
    #plt.hist(ratio3, bins=200)
    #plt.legend()
    #plt.savefig("random2_his_tf.pdf", dpi=300)
    bins = np.linspace(0.01, 1, 200)

    plt.hist(ratio1, bins, alpha=0.5, label='TF')
    plt.hist(ratio2, bins, alpha=0.5, label='random set 1')
    plt.hist(ratio3, bins, alpha=0.5, label='random set 2')
    plt.legend(loc='upper right')
    plt.savefig("his_rbp.svg")
  

    compare=[ratio1,ratio2, ratio3]
    
    
    fig1, ax1 = plt.subplots()
    labels = ('TF', 'random set 1', 'random set 2') 
    fig1, ax1 = plt.subplots()
    ax1.boxplot(compare, showfliers=False)
    plt.xticks(np.arange(len(labels))+1,labels)
    plt.savefig("boxplot_rbp.svg", dpi=300)
    
    
 ##z test
    ## alpha=0.01
    ## 2 sided test
    ## z value for alpha/2 = -2.576,2.576
    a2=2.576
    
    mean1=np.mean(ratio1)
    mean2=np.mean(ratio2)
    mean3=np.mean(ratio3)
    std1=np.std(ratio1)
    std2=np.std(ratio2)
    std3=np.std(ratio3)
    print(len1, len2, len3)
    print(mean1, mean2, mean3)
    print(std1, std2, std3)
    
    
    #z_rbp_ran1=(mean2-mean1)/math.sqrt((std2^2/len2)+(std1^2/len1))
    m1=mean2-mean1
    stdmean2=std2**2/len2
    stdmean1=std1**2/len1
    sqrt1=math.sqrt(stdmean1+stdmean2)
    z_rbp_ran1=m1/sqrt1
    print (z_rbp_ran1)
    
    
    m2=mean3-mean1
    stdmean3=std3**2/len3
    sqrt2=math.sqrt(stdmean1+stdmean3)
    z_rbp_ran2=m2/sqrt1
    print (z_rbp_ran2)
                    
    m3=mean3-mean2
    sqrt3=math.sqrt(stdmean2+stdmean3)
    z_rbp_ran3=m3/sqrt3
    print (z_rbp_ran3)
    
    ## 99% confidence interval
    CI1_pos=m1+(a2*sqrt1)
    CI1_neg=m1-(a2*sqrt1)
    CI2_pos=m2+(a2*sqrt2)
    CI2_neg=m2-(a2*sqrt2)
    CI3_pos=m3+(a2*sqrt3)
    CI3_neg=m3-(a2*sqrt3)
    
    print ("[", CI1_neg, ", " , CI1_pos, "]")
    print ("[", CI2_neg, ", " , CI2_pos, "]")
    print ("[", CI3_neg, ", " , CI3_pos, "]")

