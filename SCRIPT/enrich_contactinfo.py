'''
Created on 29 nov 2019 
given a set of meshes determine which cells are making contact (not using trimesh, strictly vtp from mpacts)
<deprecated> : now contained in extract_contact_info.py
        
@author: wimth
'''


import vtk
from VTK_utils import extract_selection_vtp, read_vtp_file, write_vtp_file, write_stl_file, get_aggregate_data, add_array,enrich_embryo_with_contactID
from vtk.util.numpy_support import vtk_to_numpy
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs
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
    param_xml = Param_xml.get_param_xml(sys.argv,l_main_keys = ['body','enrich_contactinfo'],verbose=True) #param file must be passed as first argument
    
    input_folder  = param_xml.get_value('input_folder',['paths'])
    input_file_nb = param_xml.get_value('input_file_nb',['paths'])
    output_folder  = param_xml.get_value('output_folder',['paths'])
    
    return input_folder, input_file_nb ,output_folder
    
def write_contact_csv(dd_cell1_cell2_contactinfo):
    l_cell1 = []
    l_cell2 = []
    l_contactarea = []
    for cell1,d_cell2_contactinfo in dd_cell1_cell2_contactinfo.items():
        for cell2, contactinfo in d_cell2_contactinfo.items():
            l_cell1.append(cell1)
            l_cell2.append(cell2)
            l_contactarea.append(contactinfo.area)
            
    df = pd.DataFrame({'cell1':l_cell1, 'cell2':l_cell2, 'contactarea':l_contactarea})
    df.to_csv(str(output_folder / "contactinfo.csv"),header=True,index=False)  

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
    plt.savefig(output_folder / "contactarea.png" )
    return




input_folder, input_file_nb ,output_folder = read_parms()
output_folder.mkdir(parents=True,exist_ok=True)

# input_file_vtp = ""
# l_files = [(vtp_i,int(re.findall(r'\d+',vtp_i.stem)[0])) for vtp_i in input_folder.glob('*.vtp')]
# for t_vtp in sorted(l_files, key=lambda t: t[1]):  #sorted on file number ascending
# #for vtp_i in sorted((input_folder).glob('*.vtp')):
#     if input_file_nb:
#         if input_file_nb in re.findall(r'\d+',vtp_i.stem): 
#             input_file_vtp = vtp_i
#             break
#     input_file_vtp = vtp_i

#pick a VTP from folder to process, (or pick file directly via nb)
input_file_vtp = ""
l_files = [(vtp_i,int(re.findall(r'\d+',vtp_i.stem)[0])) for vtp_i in input_folder.glob('*.vtp')]
for t_vtp in sorted(l_files, key=lambda t: t[1]):  #sorted on file number ascending
    vtp_i, nb_vtp = t_vtp
    if input_file_nb:
        if int(input_file_nb) == int(nb_vtp):
            input_file_vtp = vtp_i
            break
input_file_vtp = vtp_i

#part1 : enrich embryo vtp with contactID
# - contact index is the index of the contacting triangle (from the other cell)
# - contactID is the parent index of the contacting cell
print('vtp file {0} will be used for enriching of contactinformation'.format(input_file_vtp)) if input_file_vtp else print('no vtp file found ! at {0}'.format(input_file_vtp))

poly_in = read_vtp_file(input_file_vtp)
poly_out = enrich_embryo_with_contactID(poly_in)


path_file = str(output_folder / (input_file_vtp.stem + "_enriched.vtp")) 
write_vtp_file(poly_out, path_file)
if verbose:print("output written : {0}".format(path_file))  

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

a_contact_id = dsa.WrapDataObject(poly_out).CellData['contactId'].__array__()
a_parentIndex = dsa.WrapDataObject(poly_out).CellData['parentIndex'].__array__()
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





  