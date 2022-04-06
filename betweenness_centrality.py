import os,sys
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import random
import math
from collections import OrderedDict as od
from collections import Counter as ct
import operator
from random import sample


def two_sided_z_test(r1, r2):
    
    ## alpha=0.01
    ## 2 sided test
    ## z value for alpha/2 = -2.576,2.576
    a2=2.576 
    mean1=np.mean(r1)
    mean2=np.mean(r2)
    std1=np.std(r1)
    std2=np.std(r2)
    #print(len1, len2, len3)
    #print(mean1, mean2, mean3)
    #print(std1, std2, std3)
    m1=mean2-mean1
    stdmean2=std2**2/len2
    stdmean1=std1**2/len1
    sqrt1=math.sqrt(stdmean1+stdmean2)
    z_1_2=m1/sqrt1
    print (z_1_2)
    
    ## 99% confidence interval
    CI_pos=m1+(a2*sqrt1)
    CI_neg=m1-(a2*sqrt1)
      
    print ("Confidence interval=", "[", CI_neg, ", " , CI_pos, "]")
    if CI_pos>0>CI_neg or CI_pos<0<CI_neg:
        print ("two ratios are from different populations with significant difference")
    else:
        print ("two ratios are from the same population, no significant difference")
    
    return()

#def remove_repeats(indict):
   # for key, value in indict.items():
   #     indict[key]=list(set(value))
    #return   
    
def file_to_dict(any_file): #for a 2 column file
    any_dict={}
    f1=open(any_file,'r')
    for line in f1:
        stripped=line.strip()
        splitted=stripped.split()
        #print (splitted[1])
        any_dict.setdefault(splitted[0],splitted[1])
    f1.close()
    return any_dict

    
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

def get_interactor_distribution(conn_dict):
    uniprot_list=[]
    how_many_are={}
    which_are={}
    for key in conn_dict.keys():
        uniprot_list.append(key)       
    for query in uniprot_list:
        l1=len(conn_dict[query])
        keys= ["zero", "less_than_5", "between_5_10", "between_10_20", "more_than_20"]
        if l1==0:
            how_many_are.setdefault(keys[0],[]).append(l1)
            which_are.setdefault(keys[0],[]).append(query)
        if l1<=5:
            how_many_are.setdefault(keys[1],[]).append(l1)
            which_are.setdefault(keys[1],[]).append(query)
        elif (l1>5) and (l1<=10):
            which_are.setdefault(keys[2],[]).append(query)
            how_many_are.setdefault(keys[2],[]).append(l1)
        elif (l1>10) and (l1<=20):
            how_many_are.setdefault(keys[3],[]).append(l1)
            which_are.setdefault(keys[3],[]).append(query)
        elif (l1>20):
            how_many_are.setdefault(keys[4],[]).append(l1)
            which_are.setdefault(keys[4],[]).append(query)
            
    return how_many_are
                
        

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
    return len_array, len_list, l


def get_disorder_map(disorder_map, conn_dict, tot_conn_dict, unid_gene_dict):
    dis_dict={}
    disorder=[]
    #for query in conn_dict
    len_dict={}
    length=[]
    length2=[]
    uniprot_list=[]
    for key in conn_dict.keys():
        uniprot_list.append(key)       
    for query in uniprot_list:
        l1=len(conn_dict[query])
        l2=len(tot_conn_dict[query])
        length.append(l1)   
        length2.append(l2)
        prot_conn_list=np.array(length)
        tot_conn_list=np.array(length2)
        d1=float(disorder_map[query])
        #if d1==0.0:
            #print(query)
        disorder.append(d1)
        gene_query=unid_gene_dict[query]
        dis_dict.setdefault(gene_query, d1)
        dis_list=np.array(disorder)
    return dis_list, prot_conn_list, tot_conn_list

def common_counter(conn_dict):
    list_all=[]
    uniprot_list=[]
    number_list=[]
    for key in conn_dict.keys():
        uniprot_list.append(key)       
    for query in uniprot_list:
        arb=list(conn_dict[query])
        list_all=list_all + arb
    counter = ct(list_all)
    counter.most_common(1)
    print ("most common= ",counter.most_common(1))
    sort_c = dict( sorted(counter.items(), key=operator.itemgetter(1),reverse=True)) 
    c = od(sort_c)
    val_list=[]
    key_list=[]
    for query in c.keys():
        val_list.append(c[query])
        key_list.append(query)
    np.array(val_list)
    l=len(val_list)
    top_10={}
    key_at_index = list(c.keys())[int(l/10)]
    for element in key_list:
        val=c[element]
        top_10.setdefault(element,[]).append(val)
        if element==key_at_index:
            break
    top_10_list=list(top_10.keys())
    
    return c, top_10_list

def common_analysis(common_count):
    common_5=[]
    common_10=[]
    common_50=[]
    common_100=[]
    common_more_than_100=[]
    for query in common_count.keys():
        if 5>common_count[query]>=2:
            common_5.append(common_count[query])
        if 10>common_count[query]>=5:
            common_10.append(common_count[query])
        if 50>common_count[query]>=10:
            common_50.append(common_count[query])
        if 100>common_count[query]>=50:
            common_100.append(common_count[query])
        if common_count[query]>100:
            common_more_than_100.append(common_count[query])
     
    val_list= [[],[],[],[],[]]
    # val_list = np.empty((0, 4), int)        
    val_list[0]=common_5
    val_list[1]=common_10
    val_list[2]=common_50
    val_list[3]=common_100
    val_list[4]=common_more_than_100
    frequency=[['<5', '<10', '<50', '<100', '>100'],[]]
    
    for i in range (0,5):
            frequency[1].append(len(val_list[i]))
    
    return val_list, frequency

def get_network(conn_dict):
    edges_list=[]
    nodes_list=[]
    for key,value in conn_dict.items():
        nodes_list.append(key)
        for val in value:
            edge=(key,val)
            edges_list.append(edge)
            nodes_list.append(val)
            
    nodes_list=list(set(nodes_list))
    edges_list=list(set(edges_list))
    

    G = nx.Graph()
    G.add_nodes_from(nodes_list)
    G.add_edges_from(edges_list)
    

    #print(G.number_of_edges())
    #print(G.number_of_nodes())
    
    nx.is_connected(G)
    nx.number_connected_components(G)
    sg = (G.subgraph(c) for c in nx.connected_components(G))
    
   

    # remove low-degree nodes
    low_degree = [n for n, d in G.degree() if d < 10]
    G.remove_nodes_from(low_degree)

    # largest connected component
    components = nx.connected_components(G)
    largest_component = max(components, key=len)
    H = G.subgraph(largest_component)


    # compute centrality
    centrality = nx.betweenness_centrality(G, k=10, endpoints=True)
    #nx.degree_centrality

    # compute community structure
    lpc = nx.community.label_propagation_communities(G)
    community_index = {n: i for i, com in enumerate(lpc) for n in com}

    #### draw graph ####
    fig, ax = plt.subplots(figsize=(20, 15))
    pos = nx.spring_layout(G, k=0.15, seed=4572321)
    node_color = [community_index[n] for n in G]
    node_size = [v * 20000 for v in centrality.values()]
    nx.draw_networkx(
        G,
        pos=pos,
        with_labels=False,
        node_color=node_color,
        node_size=node_size,
        edge_color="gainsboro",
        alpha=0.4,
    )

    # Title/legend
    font = {"color": "k", "fontweight": "bold", "fontsize": 20}
    ax.set_title("TF interactome network", font)
    # Change font color for legend
    font["color"] = "r"

    ax.text(
        0.80,
        0.10,
        "node color = community structure",
        horizontalalignment="center",
        transform=ax.transAxes,
        fontdict=font,
    )
    ax.text(
        0.80,
        0.06,
        "node size = betweeness centrality",
        horizontalalignment="center",
        transform=ax.transAxes,
        fontdict=font,
    )

    # Resize figure for label readibility
    ax.margins(0.1, 0.05)
    fig.tight_layout()
    plt.axis("off")
    plt.savefig("centrality_plot_tf.pdf", dpi=300)
    
    # compute centrality
    #centrality = nx.betweenness_centrality(G, k=10, endpoints=True)
    sort_c = dict( sorted(centrality.items(), key=operator.itemgetter(1),reverse=True))
    c = od(sort_c)
    val_list=[]
    key_list=[]
    for query in c.keys():
        val_list.append(c[query])
        key_list.append(query)
    np.array(val_list)
    l=len(val_list)
    top_10={}
    key_at_index = list(c.keys())[20]
    for element in key_list:
        val=c[element]
        top_10.setdefault(element,[]).append(val)
        if element==key_at_index:
            break
    top_10_list=list(top_10.keys())
    print (top_10_list)

    #print (centrality)
    
    return (0)


if __name__=="__main__":
    #list1='RBP_list.txt'
    list1='TF_list.txt'
    list3='PTM_list.txt'
    disorder_file='proteome_disorder_coverage_30.txt'
    interactome_file=r"human_annotated_PPIs_unids_w_methods_evidence_nopred.txt"
    uniprot_gene_file=r'UNIPROTIDS_GENENAME.txt'
    outfile="out.txt"
    
    unid_gene_map,gene_unid_map=map_uniprot_gene(uniprot_gene_file)
    totcdict_un, totcdict_gene=get_total_connectivity(interactome_file,unid_gene_map)
    query_list1, len_list1=get_uni_ids(list1)
    #query_list1, stupid_list1, len_list2=gene_to_unid(list1,gene_unid_map)
    query_list3, len_list3=get_uni_ids(list3)

    cdict_un,cdict=get_proteins_connectivity(totcdict_un,unid_gene_map,query_list1)
    pp_dict_un,pp_dict=protein1_protein2_connectivity(totcdict_un, unid_gene_map,query_list1,query_list2)

    conn_array, conn_list, conn_len=get_len_array(unid_gene_map, cdict_un)
    
    
    get_network(cdict_un)
