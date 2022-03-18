import os,sys
import matplotlib.pyplot as plt
import numpy as np


def remove_repeats(indict):
    for key, value in indict.items():
        indict[key]=list(set(value))
    return

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

    return uniprot_list, gene_list

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
    
    #def remove_repeats(indict):
     #   for key, value in indict.items():
      #      indict[key]=list(set(value))
       # return indict

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
        pp_dict_unid.setdefault(query,common)
        gene_query=unid_gene_dict[query]
        gene_conns=[unid_gene_dict[conn] for conn in common]
        pp_dict_gene.setdefault(gene_query,gene_conns)
       
    return pp_dict_unid, pp_dict_gene

def get_length_array(gene_list, conn_dict_gene):
    for query in gene_list:
        l1=len(conn_dict_gene[query])
        len_dict.setdefault(query,l1)
        data = list(len_tot_inter_dict.items())
        len_array = np.array(data)
    return len_array

if __name__=="__main__":
    list1='RBP_list.txt'
    list2='TF_list.txt'
    interactome_file=r"human_annotated_PPIs_unids_w_methods_evidence_nopred.txt"
    uniprot_gene_file=r'UNIPROTIDS_GENENAME.txt'
    outfile="out.txt"
    unid_gene_map,gene_unid_map=map_uniprot_gene(uniprot_gene_file)
    #print (unid_gene_map)
    #print (gene_unid_map)
    totcdict_un, totcdict_gene=get_total_connectivity(interactome_file,unid_gene_map)
    #totcdict_gene, totcdict_un=get_total_connectivity(interactome_file,gene_unid_map)
    #print (totcdict_un)
    query_list1, gene_list1=gene_to_unid(list1,gene_unid_map)
    query_list2, gene_list1=gene_to_unid(list2,gene_unid_map)
    #print (query_list2)
    #print (query_list1)
    cdict_un,cdict=get_proteins_connectivity(totcdict_un,unid_gene_map,query_list1)
    #print(cdict)
    pp_dict_un,pp_dict=protein1_protein2_connectivity(totcdict_un, unid_gene_map,query_list1,query_list2)
    #print (pp_dict)
    totc_len_array=get_length_array(cdict_gene, gene_list1)
    pp_len_array=get_length_array(pp_dict, gene_list1)
    print (totc_len_array)

