import pandas as pd
import csv
import py4cytoscape as p4c
from py4cytoscape import gen_node_color_map
from py4cytoscape import palette_color_brewer_d_RdBu
import os
from inputs import *


def process_labels_file(labels_file):

    labels = pd.read_csv("../"+labels_file, sep = '|', usecols = ['id','category','name','description','xref'])
    labels.columns = ['id','category', 'label','description/definition','xref']

    return labels


#subgraph_df is a dataframe with S, P, O headers and | delimited
def create_node_attributes(input_nodes_df,input_mech):
        
    full_list = []

    for i in range(len(input_mech)):
        #Only add subject and object columns, not the predicate
        for col in [0,2]:
            l = []
            node = input_mech.iloc[i,col]
            for j in range(len(input_nodes_df)):
                try:
                    cat = input_nodes_df.loc[input_nodes_df['label'].str.lower() == node.lower(),"category"].values[0]
                    att = cat.split('biolink:')[1]
                except IndexError: 
                    print('Node not included in labels file: ',node,'. Update labels file.')
                    exit()
            l.append(node)
            l.append(att)
            full_list.append(l)
    
    subgraph_attribute_df = pd.DataFrame(full_list,columns = ['Node','Attribute'])
    
    subgraph_attribute_df = subgraph_attribute_df.drop_duplicates(subset=['Node'])
    subgraph_attribute_df = subgraph_attribute_df.reset_index(drop=True)
    
    return subgraph_attribute_df

#subgraph_df is a dataframe with S, P, O headers and | delimited
def create_noa_file(input_mech_file,subgraph_attribute_df):

    orig_filename = input_mech_file.split('.csv')[0].split('/')[1]

    noa_file = "../"+input_mech_file.split('/')[0]+"/"+orig_filename+"_Subgraph_attributes.noa"

    subgraph_attribute_df.to_csv(noa_file,index=False,sep='|')

#subgraph_df is a dataframe with S, P, O headers and | delimited
def create_cytoscape_png(input_mech_file,input_mech,subgraph_attributes_df):

    orig_filename = input_mech_file.split('.csv')[0].split('/')[1]

    png_file = "../"+input_mech_file.split('/')[0]+"/"+orig_filename+"_Subgraph_Visualization.png"

    #Update column names for cytoscape
    input_mech.columns = ['source','interaction','target']
    subgraph_attributes_df.columns = ['id','group']

    p4c.create_network_from_data_frames(subgraph_attributes_df,input_mech,title='subgraph')

    #Ensure no network exists named subgraph in Cytoscape or you will have to manually override before it can be output
    p4c.set_visual_style('BioPAX_SIF',network='subgraph')

    p4c.set_node_color_mapping(**gen_node_color_map('group', mapping_type='d',style_name='BioPAX_SIF'))

    p4c.set_edge_label_mapping('interaction',style_name='BioPAX_SIF')

    p4c.set_edge_target_arrow_shape_default('ARROW',style_name='BioPAX_SIF')
    
    p4c.export_image(png_file,network='subgraph')

# Wrapper Function
def output_visualization(input_mech_file,input_nodes_df): #input_nodes_df,

    input_mech = pd.read_csv("../"+input_mech_file,sep='|')

    subgraph_attributes_df = create_node_attributes(input_nodes_df,input_mech)

    create_noa_file(input_mech_file,subgraph_attributes_df)

    create_cytoscape_png(input_mech_file,input_mech,subgraph_attributes_df)

def main():

    input_mech_file,labels_file = generate_arguments()

    labels_df = process_labels_file(labels_file)

    output_visualization(input_mech_file,labels_df)

if __name__ == '__main__':
    main()
