'''
Created on 29 nov 2019 
extacts the contactinfo from the a mpacts pics embryo.  The file needs to be enriched with 'contactid' (cfr mpacts_pics_enrich.py)
        
@author: wimth
'''

from VTK_utils import read_vtp_file, get_aggregate_data, get_data_array,get_mapper
from vtk.numpy_interface import dataset_adapter as dsa
import numpy as np
import pandas as pd
import sys,os,re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import networkx as nx

from param_XML import Param_xml
verbose=True


def read_parms():
    param_xml = Param_xml.get_param_xml(sys.argv, l_main_keys=['body'],verbose=True) 
    param_xml.add_prms_to_dict(l_keys=['extract_contactinfo'])
    param_xml.repr_prm()
    return param_xml.prm
    
def write_contact_csv(dd_cell1_cell2_contactinfo):
    l_cell1 = []
    l_cell2 = []
    l_contactarea = []
    for cell1,d_cell2_contactinfo in dd_cell1_cell2_contactinfo.items():
        for cell2, contactinfo in d_cell2_contactinfo.items():
            l_cell1.append(cell1)
            l_cell2.append(cell2)
            l_contactarea.append(contactinfo.area)

    if d_parentid_cellname:
        l_cellname1 = [d_parentid_cellname.get(i,'cellName not found') for i in l_cell1]
        l_cellname2 = [d_parentid_cellname.get(i,'cellName not found') for i in l_cell2]
        df = pd.DataFrame({'cell1':l_cell1, 'cell2':l_cell2,'name1':l_cellname1,'name2':l_cellname2, 'contactarea':l_contactarea})
    else:
        df = pd.DataFrame({'cell1':l_cell1, 'cell2':l_cell2, 'contactarea':l_contactarea})

    df.to_csv(str(prm['output_folder'] / "contactinfo.csv"),header=True,index=False)  

    return
    
def write_network_graph():
    G = nx.Graph()
    max_contactarea = 0
    l_contact_area = []
    for cell1,d_cell2_contactinfo in dd_cell1_cell2_contactinfo.items():
        for cell2, contactinfo in d_cell2_contactinfo.items():
                area = contactinfo.area
                G.add_edge(cell1, cell2, weight=area)
                l_contact_area.append(area)
                if area > max_contactarea:max_contactarea=area
                
    weights = [20*(G[u][v]['weight'] / max_contactarea) for u,v in G.edges()]
    
    pos = nx.spring_layout(G)  # positions for all nodes
    nx.draw_networkx_nodes(G, pos, node_size=700)
    # nx.draw_networkx_edges(G, pos,  width=(5* (l_contact_area / max_contactarea)))
    nx.draw_networkx_edges(G, pos,  width=weights)
    nx.draw_networkx_labels(G, pos, font_size=20, font_family="sans-serif")
    
    plt.axis("off")
    plt.savefig(prm['output_folder'] / "contactarea.png" )
    return

prm = read_parms()
prm['output_folder'].mkdir(parents=True,exist_ok=True)

input_file_vtp = prm['input_file'] if prm['input_file'] else prm['input_folder'] / 'selected_embryo.vtp'
poly_in = read_vtp_file(input_file_vtp)

#part2 : make contactarea network graph
#-step1) build datastructure containing all contact pair data dd_cell1_cell2_contactinfo
class DictOfDicts(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value
            
class Contactinfo:
    def __init__(self, l_idx_triangles,area):
        self.l_idx_triangles = l_idx_triangles
        self.area = area

#a_contact_id = dsa.WrapDataObject(poly_in).CellData['contactId'].__array__()
a_contact_id = get_data_array(poly_in, field_type="CELL", attribute='contactid',verbose=True)
#a_parentIndex = dsa.WrapDataObject(poly_in).CellData['parentIndex'].__array__()
a_parentIndex = get_data_array(poly_in, field_type="CELL", attribute='parentIndex',verbose=True)

d_parentid_cellname = get_mapper(poly_in,name_array_key = "parentIndex", name_array_value="cellName",field_type="CELL")

dd_cell1_cell2_contactinfo = DictOfDicts()  #dict of dicts : key 1= cell1, key2 = cell2 , value = contactinfo object
s_pairs = set([(i,j) for i,j in zip(a_parentIndex,a_contact_id)])

for t_pair in s_pairs:
    cell1,cell2 = t_pair
    l_idx_triangles = np.where( (a_parentIndex==cell1) & (a_contact_id==cell2) )
    contactarea = get_aggregate_data(poly_in, np.array(l_idx_triangles), l_aggregate_columns=['area'], func="sum")
    dd_cell1_cell2_contactinfo[cell1][cell2] = Contactinfo(l_idx_triangles,contactarea[0])

#-step2) write output
write_contact_csv(dd_cell1_cell2_contactinfo)
write_network_graph()





  