#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 18:52:56 2023

@author: brooksantangelo
"""

import os
from urllib.request import urlretrieve
import pandas as pd
import re
from igraph import * 




def read_ocr_input(user_input_file):
	df = pd.read_csv(user_input_file, sep = "\t")
	return(df)



def interactive_search_wrapper(filename,g,known_mappings):
    df = read_ocr_input(filename)
    n = []
    if "genes" in filename:
        nodes = unique_nodes(df["ncbigene_id"].to_frame())
        nodes = ["http://www.ncbi.nlm.nih.gov/gene/" + str(n) for n in nodes]
        nodes_udpated = gene_to_protein(nodes, g)
        n.append(nodes_udpated)
        n = [item for items in n for item in items]
    else: 
        nodes = unique_nodes(df["word"].to_frame())
        for i in nodes:
            i_node,known_mappings = search_nodes(i,g,known_mappings)
            i_node_short = i_node.split('http://purl.obolibrary.org/obo/')[1].replace('_',':')
            
            n.append(i_node_short)
    
    return n,known_mappings

#Takes in list and converts gene uris to protein uris
def gene_to_protein(nodes,g):
    
    pr = []
    predicate = 'http://purl.obolibrary.org/obo/pr#has_gene_template'
    
    #looking for http://purl.obolibrary.org/obo/pr#has_gene_template rels, 
    #19521 has_gene_product and 19834 has_gene_template 
    for i in nodes:
        try:
            p = g.edgelist.loc[(g.edgelist['object'] == i) & (g.edgelist['predicate'] == predicate),'subject'].values[0]
            #p_lab = g.labels_all.loc[g.labels_all['entity_uri'] == p,'label'].values[0]
            p_short = p.split('http://purl.obolibrary.org/obo/')[1].replace('_',':')
            pr.append(p_short)
        except IndexError: continue
        
    return pr
    

def unique_nodes(examples):
	# get unique node values
	nodes = list(set(pd.melt(examples)["value"].values))
	return(nodes)


def map_input_to_nodes(node,kg):

	search_loop = True

	while(search_loop):
		print("User Search Node: ", node)
		found_nodes = find_node(node,kg)
		nrow = found_nodes.shape[0]
		if nrow == 0:
			print("No search terms returned")
			node = input("Please try another input term: ")
		else:
			search_loop = False	
	print("Found", nrow, "features in KG")

	return found_nodes,nrow

def find_node(node, kg, ontology = ""):
    nodes = kg.labels_all
    ### Check for exact matches first
    exact_matches = nodes[(nodes["label"].str.lower() == node.lower())|(nodes["entity_uri"].str.lower() == node.lower())][["label", "entity_uri"]]

    results = nodes[nodes["label"].str.contains(node,flags=re.IGNORECASE, na = False)|nodes["synonym"].str.contains(node,flags=re.IGNORECASE, na = False)|nodes["description/definition"].str.contains(node,flags=re.IGNORECASE, na = False)][["label", "entity_uri"]]
	#Remove exact matches from this df
    results = results[(~results.label.isin(exact_matches.label))]

    # sort results by ontology
    results = results.sort_values(['entity_uri'])

    #Concat both dfs so that exact matches are presented first
    all_results = pd.concat([exact_matches, results], axis=0)
    return(all_results)



# Check if search input is in the list of integer_ids
def node_in_search(found_nodes, user_input):
	if user_input in found_nodes[["label"]].values:
		return(True)
	else:
		return(False)

def manage_user_input(found_nodes,user_input,kg):
    
    if node_in_search(found_nodes,user_input):
        #Manage if there are 2 duplicate label names
        if len(found_nodes[found_nodes['label'] == user_input][['label','entity_uri']]) > 1:
            dup_node = True
            while(dup_node):
                user_id_input = input("Input node 'id': ")
                print(found_nodes[found_nodes['label'] == user_input]['entity_uri'].values.tolist())
                if user_id_input in found_nodes[found_nodes['label'] == user_input]['labels'].values.tolist():
					
                    node_id = kg.labels_all.loc[kg.labels_all['label'] == user_id_input,'entity_uri'].values[0]
                    bad_input = False
                    dup_node = False
                    
                else: 
                    print("Input id does not correspond with selected label.... try again")
        else:
            node_label = user_input
            node_id = kg.labels_all.loc[kg.labels_all['label'] == node_label,'entity_uri'].values[0]

            bad_input = False
            dup_node = False

    elif node_in_labels(kg,user_input):
        node_label= user_input
        node_id = kg.labels_all.loc[kg.labels_all['label'] == node_label,'entity_uri'].values[0]
        bad_input = False
    else:
        print("Input not in search results.... try again")
        node_label = ""
        bad_input = True


    return node_id,bad_input


def search_nodes(node,kg,known_mappings):
    print('Searching for node: ',node)
    
    vals_per_page = 20
    
    bad_input = True
    
    if node in list(known_mappings.keys()):
        node_label = known_mappings[node]
        return node_label
    
    else:
        found_nodes,nrow = map_input_to_nodes(node,kg)
        i = 1
        while(bad_input):	
            high = min(nrow,(i)*vals_per_page)	
            print(found_nodes.iloc[(i-1)*vals_per_page:high,].to_string())
            
            user_input = input("Input node 'label' or 'f' for the next " + str(vals_per_page) + " features, 'b' for the previous " + str(vals_per_page) + ", or 'u' to update the node search term: ")	
            if user_input == 'f':		
                if (nrow / i ) > vals_per_page:		
                    i+=1
            elif user_input == 'b':	
                if i > 1:
                    i-=1
            elif user_input == 'u':
                node = input("Input new node search term: ")
                found_nodes,nrow = map_input_to_nodes(node,kg)
                i = 1
            else:
                node_label,bad_input = manage_user_input(found_nodes,user_input,kg)
            
        known_mappings[node] = node_label
        return node_label,known_mappings

def process_pkl_files(triples_file,labels_file):
    #Edgelist uses IDs
    triples_df = pd.read_csv(triples_file,sep = '\t')
    triples_df.columns = triples_df.columns.str.strip().str.lower()
    triples_df = triples_df[['subject', 'predicate', 'object']]
    triples_df.replace({'<': ''}, regex=True, inplace=True)
    triples_df.replace({'>': ''}, regex=True, inplace=True)
    #triples_df.columns.str.lower()

    #labels = pd.read_csv(labels_file, sep = '\t', usecols = ['id','category','name','description'])  #,'iri','synonym'])
    labels = pd.read_csv(labels_file, sep = '\t', low_memory=False)
    labels.columns = labels.columns.str.strip().str.lower()
    labels = labels.rename(columns={'identifier':'entity_uri'})
    labels.replace({'<': ''}, regex=True, inplace=True)
    labels.replace({'>': ''}, regex=True, inplace=True)

    #Handle nodes with no label
    labels.loc[pd.isna(labels["label"]),'label'] = labels.loc[pd.isna(labels["label"]),'entity_uri']

    return triples_df,labels

class KnowledgeGraph:

    def __init__(self, edgelist, labels_all, igraph, igraph_nodes):
        self.edgelist = edgelist
        self.labels_all = labels_all
        self.igraph = igraph
        self.igraph_nodes = igraph_nodes


def create_igraph_graph(edgelist_df,labels):

    edgelist_df = edgelist_df[['subject', 'object', 'predicate']]

    g = Graph.DataFrame(edgelist_df,directed=True,use_vids=False)

    g_nodes = g.vs()['name']

    return g,g_nodes

def create_graph(triples_file,labels_file, kg_type):
    triples_df,labels = process_pkl_files(triples_file,labels_file)
        
    g_igraph,g_nodes_igraph = create_igraph_graph(triples_df,labels)
    # Created a PKL class instance
    pkl_graph = KnowledgeGraph(triples_df,labels,g_igraph,g_nodes_igraph)

    return pkl_graph


def get_graph():
    triples_list_file = '/Users/brooksantangelo/Documents/HunterLab/Cartoomics/PostRevisionUpdates/Inputs/pkl/PheKnowLator_v3.0.2_full_instance_relationsOnly_OWLNETS_Triples_Identifiers.txt'
    labels_file = '/Users/brooksantangelo/Documents/HunterLab/Cartoomics/PostRevisionUpdates/Inputs/pkl/PheKnowLator_v3.0.2_full_instance_relationsOnly_OWLNETS_NodeLabels.txt'
    
    kg_type = 'pkl'
    
    g = create_graph(triples_list_file,labels_file, kg_type)
    
    return g




def find_concept_annotations(pw,known_mappings):
    
    g = get_graph()
    
    folder = os.getcwd()+'/results/'

    concepts_found = []

    pfocr_id = pw.split("/")[-1].split(".")[0]
    ### Could add the Figure ID here. LG
    if not os.path.isdir(folder):
        # raise Exception('Missing folder input directory: ' + folder)
        # logging.error('Missing folder input directory: ' + folder)
        os.makedirs(folder)
    
    mentions = ["genes", "chemicals", "diseases"]
    for mention in mentions:
        print('known_mappings: ',known_mappings)
        url = "https://raw.githubusercontent.com/wikipathways/pfocr-database/main/_data/" + pfocr_id+ "-"+mention+".tsv"
        filename = folder + mention+".tsv"
        try:
            urlretrieve(url, filename)
            print("downloaded " + mention)
            c,known_mappings = interactive_search_wrapper(filename,g,known_mappings)
            concepts_found = concepts_found + c
            print(concepts_found)
        except:
            print("no content for " + mention)
            
    return concepts_found,known_mappings
    
    
    
def get_uri(value,labels_file):
    
    labels = pd.read_csv(labels_file, sep = '|', low_memory=False)
    labels.columns = labels.columns.str.strip().str.lower()

    
    id = labels.loc[labels['name'] == value,'id'].values[0]
        
    return id

def get_nodes_from_mech(filename,labels_file):
    
    df = pd.read_csv(filename,sep='|')
    
    l = df[['S','O']].values.tolist()
    l = list(set(sum(l, [])))
    l_ids = []
    for i in l:
        l_ids.append(get_uri(i,labels_file))
    
    return l_ids


pw_list = ['https://pfocr.wikipathways.org/figures/PMC6943888__40035_2019_179_Fig1_HTML.html',
           'https://pfocr.wikipathways.org/figures/PMC6365454__fnins-13-00025-g0002.html',
           'https://pfocr.wikipathways.org/figures/PMC6218675__fphys-09-01533-g0001.html',
           'https://pfocr.wikipathways.org/figures/PMC6143874__etm-16-04-3275-g02.html',
           'https://pfocr.wikipathways.org/figures/PMC6143874__etm-16-04-3275-g01.html',
           'https://pfocr.wikipathways.org/figures/PMC6121065__fncel-12-00273-g0001.html',
           'https://pfocr.wikipathways.org/figures/PMC6115613__fnmol-11-00244-g004.html',
           'https://pfocr.wikipathways.org/figures/PMC6115613__fnmol-11-00244-g003.html',
           'https://pfocr.wikipathways.org/figures/PMC6103457__BST-46-967-g0001.html',
           'https://pfocr.wikipathways.org/figures/PMC6087706__40169_2018_202_Fig2_HTML.html',
           'https://pfocr.wikipathways.org/figures/PMC6085908__BRB3-8-e00978-g001.html',
           'https://pfocr.wikipathways.org/figures/PMC6027339__biomedicines-06-00038-g005.html',
           'https://pfocr.wikipathways.org/figures/PMC5992185__41392_2018_15_Fig4_HTML.html',
           'https://pfocr.wikipathways.org/figures/PMC5923349__cancers-10-00094-g004.html',
           'https://pfocr.wikipathways.org/figures/PMC5923345__cancers-10-00090-g002.html',
           'https://pfocr.wikipathways.org/figures/PMC5923345__cancers-10-00090-g001.html',
           'https://pfocr.wikipathways.org/figures/PMC5883377__CN-16-151_F1B.html',
           'https://pfocr.wikipathways.org/figures/PMC5867295__fbioe-06-00018-g005.html',
           'https://pfocr.wikipathways.org/figures/PMC5854428__ppat.1006839.g001.html',
           'https://pfocr.wikipathways.org/figures/PMC5833345__41419_2017_66_Fig7_HTML.html',
           'https://pfocr.wikipathways.org/figures/PMC5822808__JIR2018-3696914.001.html',
           'https://pfocr.wikipathways.org/figures/PMC5819904__IJMM-41-03-1201-g03.html',
           'https://pfocr.wikipathways.org/figures/PMC5452883__ol-13-06-4047-g00.html',
           'https://pfocr.wikipathways.org/figures/PMC3514635__gr4.html']

#pw = 'https://pfocr.wikipathways.org/figures/PMC5095497__IMM-149-423-g007.html'
#pw_concepts = find_concept_annotations(pw)

    
 
   
mechs_labels_file = '/Users/brooksantangelo/Documents/Repositories/Microbiome_Disease_Models/nodes.csv'
mech_list = get_nodes_from_mech('/Users/brooksantangelo/Documents/Repositories/Microbiome_Disease_Models/ParkinsonsDisease/Salmonella_typhi_1.csv',mechs_labels_file)
print('mech_list')
print(mech_list)

known_mappings = {}
 
for pw in pw_list:
    
    print(pw)
    
    pw_concepts,known_mappings = find_concept_annotations(pw,known_mappings)

    print('pw_concepts')
    print(pw_concepts)

    overlapping_concepts = list(set(pw_concepts) & set(mech_list))
    print('overlapping_concepts')
    print(overlapping_concepts)


#Search through annotated mechanisms to see if any node pairs exist that also exist in one wikipathways diagram




