import os,sys
import networkx as nx
import matplotlib.pyplot as plt

# read in uniprot IDs for query genes
def get_uniprot(input_file):
    uniprot_list=[]
    fin=open(input_file,"r")
    for line in fin:
        stripped=line.strip()
        splitted=stripped.split()
        uniprot_list.append(splitted[1])
    fin.close()
    
    return uniprot_list

# get total connectivity documented in the data file (conn_file)
def get_total_connectivity(conn_file,unid_gene_dict):
    
    def remove_repeats(indict):
        for key,value in indict.items():
            indict[key]=list(set(value))
            
        return indict
        
    conn_dict={} #initiate a connectivity dictionary
    conn_gene_dict={}
    fin=open(conn_file,"r")
    for line in fin:
        stripped=line.strip()
        splitted=stripped.split()
        conn_dict.setdefault(splitted[0],[]).append(splitted[1])
        conn_dict.setdefault(splitted[1],[]).append(splitted[0])
        
        gene1=unid_gene_dict[splitted[0]]
        gene2=unid_gene_dict[splitted[1]]
        conn_gene_dict.setdefault(gene1,[]).append(gene2)
        conn_gene_dict.setdefault(gene2,[]).append(gene1)
                        
    fin.close()
    
    conn_dict=remove_repeats(conn_dict)
    conn_gene_dict=remove_repeats(conn_gene_dict)
    
    return conn_dict,conn_gene_dict


def get_khops_conn(vertex_list,conn_dict,k=2):
    
    khops_conn={} # stores all connections from k-hops
    it=1
    while it<k:
        it+=1
        new_vertex_list=[] # empty for every iteration
        for vertex in vertex_list:
            vconns=conn_dict[vertex]
            khops_conn.setdefault(vertex,vconns)
            for vx in vconns:
                vxconns=conn_dict[vx]
                khops_conn.setdefault(vx,vxconns)
                new_vertex_list.extend(vxconns)
                
        vertex_list=list(set(new_vertex_list)) # update the vertex list
        
    return khops_conn          


def get_proteins_connectivity(conn_dict,unid_gene_dict,query_unids):

    def remove_repeats(indict):
        for key,value in indict.items():
            indict[key]=list(set(value))
        return indict
    
    query_dict_unid={} # subset conn_dict for query and its connections
    query_dict_gene={} 
    for query in query_unids:
        query_dict_unid.setdefault(query,conn_dict[query])
        
        gene_query=unid_gene_dict[query]
        gene_conns=[unid_gene_dict[conn] for conn in conn_dict[query]]
        query_dict_gene.setdefault(gene_query,gene_conns)
    
    query_dict_unid=remove_repeats(query_dict_unid)
    query_dict_gene=remove_repeats(query_dict_gene)
    
    return query_dict_unid,query_dict_gene


#get total connectivity of a set of proteins - including any other proteins not in the set
def get_proteins_connectivity_old(conn_file,asd_proteins):
    conn_dict={} #initiate a connectivity dictionary
    fin=open(conn_file,"r")
    for line in fin:
        stripped=line.strip()
        splitted=stripped.split()
        if (splitted[0] in asd_proteins) or (splitted[1] in asd_proteins): # if one of the interacting partners is ASD-risk gene
            conn_dict.setdefault(splitted[0],[]).append(splitted[1])
            conn_dict.setdefault(splitted[1],[]).append(splitted[0])
                        
    fin.close()
    
    for key,value in conn_dict.items():
        conn_dict[key]=list(set(value))
    
    return conn_dict

# get connectivity of the set of proteins to other proteins in the same set
def get_set_connectivity(conn_file,query_proteins):
    conn_dict={} #initiate a connectivity dictionary
    fin=open(conn_file,"r")
    for line in fin:
        stripped=line.strip()
        splitted=stripped.split()
        #print(splitted)
        if (splitted[0] in query_proteins) and (splitted[1] in query_proteins): # only if both proteins in query
            conn_dict.setdefault(splitted[0],[]).append(splitted[1])
            conn_dict.setdefault(splitted[1],[]).append(splitted[0])
                        
    fin.close()
    
    for key,value in conn_dict.items():
        conn_dict[key]=list(set(value))
    
    return conn_dict

#
def get_tot_degrees(conn_dict,fname,tag):
    fout1=open(fname+"_"+tag+"_vdegrees.out.txt","w")
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
    # get total vertex degrees
    for node in G.nodes():
        fout1.write("%s\t%s\n"%(node,G.degree[node]))
        
            
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
    
    for entry in nodes_list:
        print(entry)
    
    G = nx.Graph()
    G.add_nodes_from(nodes_list)
    G.add_edges_from(edges_list)
    
    # export vertex degrees
    for node in G.nodes():
        fout1.write("%s\t%s\n"%(node,G.degree[node]))
    
    print(G.number_of_edges())
    print(G.number_of_nodes())
    
    # get connected subgraphs    
    sub_graphs = nx.connected_component_subgraphs(G)
    
    for i, sg in enumerate(sub_graphs):
        print("subgraph {} has {} nodes".format(i, sg.number_of_nodes()))
        #
        #sg,pos=nx.spring_layout(sg)
        
        pos=nx.spring_layout(sg) # positions for all nodes
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
    
def map_uniprot_gene(mapfile):
    unid_gene={}
    fin=open(mapfile,'r')
    for line in fin:
        stripped=line.strip()
        splitted=stripped.split()
        unid_gene.setdefault(splitted[0],splitted[1])
                
    fin.close()
    
    return unid_gene
        
if __name__=="__main__":
    genelist='RBP_list.txt'#r'C:\Users\iva.pritisanac\Documents\toronto\projects\collaborations\jditlev\Protein_Sequences_for_PSC_scaffolds\SHAN2_SHAN3.txt'
    #r"C:\Users\iva.pritisanac\Documents\toronto\projects\idrs_disease\ASD_genelist_uniprotids.txt"
    interactome_list=r"human_annotated_PPIs_unids_w_methods_evidence_nopred.txt"
    
    uniprot_gene_file=r'UNIPROTIDS_GENENAME.txt'
    outfile='Hermann_out.txt'#r"C:\Users\iva.pritisanac\Documents\toronto\projects\collaborations\jditlev\Protein_Sequences_for_PSC_scaffolds\SHANK2_CAPRIN2"

    unid_gene_map=map_uniprot_gene(uniprot_gene_file)

    totcdict_un,totcdict_gene=get_total_connectivity(interactome_list,unid_gene_map)
    
    #get_tot_degrees(totcdict_un,outfile,"total")
    
    asdunids=get_uniprot(genelist)
    
    cdict_un,cdict_gene=get_proteins_connectivity(totcdict_un,unid_gene_map,asdunids) # ONLY BETWEEN PROTEINS IN THE LIST
    
    vals=[]
    for key,value in cdict_gene.items():
        vals.append(value)
        print('\n')
        print("QUERY:",key,len(value))
        print('\n')
        
        #for val in value:
            #print(val)     
    
    for i in range(len(vals)):
        
        for j in range(len(vals)):
            if i==j:
                continue
            set_overlap=list(set(vals[i]) & set(vals[j]))
            print("OVERLAP",len(set_overlap))
            #print(set_overlap)
            for entry in set_overlap:
                print(entry)
    #sys.exit(0)
    #get_network(cdict_gene,outfile,"set_w")
    acdict=get_set_connectivity(interactome_list,asdunids)
    #print(acdict)
    get_network(acdict,outfile,"inter_only")
 
    ## Hop across the graph from nodes in asdunids
    #cdict=get_khops_conn(asdunids,totcdict_un,2)
    #get_network(cdict,outfile,"set_w")
    
