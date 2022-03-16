import os,sys
import networkx as nx
import matplotlib.pyplot as plt


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

    return uniprot_list

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
        #print (p2)
        pp_dict_unid.setdefault(query,p2)
        gene_query=unid_gene_dict[query]
        gene_conns=[unid_gene_dict[conn] for conn in conn_dict[query]]
        pp_dict_gene.setdefault(gene_query,gene_conns)
        n=query2_unids.count(p2)
        #print (n)
        if n>0:
            print ("yes")
            pp_dict_unid.setdefault(query,p2)
            gene_query=unid_gene_dict[query]
            gene_conns=[unid_gene_dict[conn] for conn in conn_dict[query]]
            pp_dict_gene.setdefault(gene_query, gene_conns)
    return pp_dict_unid, pp_dict_gene
def get_network(conn_dict,fname,tag):
    fout1=open(fname+"_"+tag+"_vdegrees.out.txt","w")
    fout2=open(fname+"_"+tag+"_cliques_size.out.txt","w")
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
    
    #for entry in nodes_list:
     #   print(entry)
    
    G = nx.Graph()
    G.add_nodes_from(nodes_list)
    G.add_edges_from(edges_list)
    
    # export vertex degrees
    for node in G.nodes():
        fout1.write("%s\t%s\n"%(node,G.degree[node]))
    
    print(G.number_of_edges())
    print(G.number_of_nodes())
    
    # get connected subgraphs    
    sub_graphs = (G.subgraph(c) for c in nx.connected_components(G))
    print (sub_graphs)
    
    for i, sg in enumerate(sub_graphs):
        #print("subgraph {} has {} nodes".format(i, sg.number_of_nodes()))
        #
        #sg,pos=nx.spring_layout(sg)
        print (nx.spring_layout(sg))
        
        pos=nx.spring_layout(sg) # positions for all nodes
        print (pos)
        node_list=[node for node in sg.nodes()]
        node_labels={node:node for node in sg.nodes()}
        

        # nodes
        nx.draw_networkx_nodes(sg,pos,
                               nodelist=node_list,
                               node_color='r',
                               node_size=1000,
                               alpha=0.8)
        
        # edges
        nx.draw_networkx_edges(sg,pos,width=1.0,alpha=0.5)

                
        #labels=nx.draw_networkx_labels(sg,pos, font_size=14, font_color='k', font_family='sans-serif', font_weight='normal', alpha=1.0)
        
        nx.draw_networkx_labels(sg,pos,labels=node_labels,font_size=8,font_color='k',font_weight='bold')
        plt.show()
        #print("\tNodes:", sg.nodes(data=True))
        
        nx.draw(G)
        nx.draw(G,pos=nx.spring_layout(G)) # use spring layout
        
        for node in sg.nodes:
            print(node)
        print("\n\n")
        #print("\tEdges:", sg.edges())
    #print(i)
    
    #nx.draw(G, with_labels=True, font_weight='bold')
    #plt.show()
    fout1.close()
    fout2.close()

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
    query_list1=gene_to_unid(list1,gene_unid_map)
    query_list2=gene_to_unid(list2,gene_unid_map)
    #print (query_list2)
    #print (query_list1)
    cdict_un,cdict=get_proteins_connectivity(totcdict_un,unid_gene_map,query_list1)
    #print(cdict)
    pp_dict_un,pp_dict=protein1_protein2_connectivity(totcdict_un, unid_gene_map,query_list1,query_list2)
    #print (pp_dict)
    get_network(pp_dict,outfile,"set_w")


