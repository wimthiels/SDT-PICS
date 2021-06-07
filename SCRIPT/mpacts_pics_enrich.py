'''
Created on 29 dec 2020 
enrich mpacts-pics segmentation with contact info and cell id's
this is basically enrich_contactinfo.py with added cell id enrichment
@author: wimth
'''


import vtk
from VTK_utils import extract_selection_vtp, read_vtp_file, write_vtp_file, write_stl_file,  \
                     get_aggregate_data, add_array, get_data_array,enrich_embryo_with_contactID,get_parentIDs,extract_parentID_selection_from_VTP,poly_to_trimesh, \
                     add_array_with_mapper
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
from id_cells import identify

def read_parms():
	param_xml = Param_xml.get_param_xml(sys.argv,l_main_keys = ['body'],verbose=True) 
	param_xml.add_prms_to_dict(l_keys=['mpacts_pics_enrich'])
	param_xml.add_prms_to_dict(l_keys=['RAM'])
	param_xml.repr_prm()

	return param_xml.prm

def load_lineage_data(verbose=False):
	#read cell names
	df_cell_names = pd.read_csv(prm['input_file_lineage'],delimiter="\t")
	df_cell_names['x'] = prm['scaling_ZYX'][2]  *1.e-6 * df_cell_names['x']
	df_cell_names['y'] = prm['scaling_ZYX'][1] * 1.e-6 * df_cell_names['y']
	df_cell_names['z'] = prm['scaling_ZYX'][0] * 1.e-6 * df_cell_names['z']
	if verbose:print(f'df_cell_names = {df_cell_names}')

	return df_cell_names

def get_canonical_parentid(cell_name_search, df_canonical_parentid):
	cond = df_canonical_parentid['cell']==cell_name_search
	return 998 if df_canonical_parentid[cond].index.size==0 else df_canonical_parentid[cond].index[0]

def name_cells_using_lineage(poly,df_cell_names):
	""" provide a mapping between parentid and a cell name based on a named lineage tsv file"""
	d_parentid_VTP = extract_parentID_selection_from_VTP(poly_embryo,verbose=verbose)
	default_cellname = 'parentID'
	d_parentid_cellname = {k:"{}{}".format(default_cellname,str(k).zfill(3)) for k in d_parentid_VTP.keys()}
	for parentid_i,vtp_i in d_parentid_VTP.items():
		tmh_i = poly_to_trimesh(vtp_i)
		for row_i in df_cell_names.iterrows():
			point_xyz = row_i[1]['x'],row_i[1]['y'],row_i[1]['z']
			if tmh_i.contains(np.array([list(point_xyz)])):  #must have shape (1,3)
				if verbose:print('mesh {} contains point {}'.format(parentid_i, point_xyz))
				d_parentid_cellname[parentid_i]  = row_i[1]['cell']
			# else:
				# if verbose:print('mesh {} does not contain point {}'.format(parentid_i, point_xyz))
	print(d_parentid_cellname)

	return d_parentid_cellname

def write_mapper_cellname_csv():
	""" write a csv to log the mapping from cell name to parentid and canonical parentid"""
	df_parentid_cellname = pd.DataFrame.from_dict(d_parentid_cellname,orient='index')
	df_parentid_cellname.reset_index(level=0, inplace=True)
	df_parentid_cellname.columns = ['parentid','cell_name']

	df_cellname_canonical_parentid = pd.DataFrame.from_dict(d_cellname_canonical_parentid,orient='index')
	df_cellname_canonical_parentid.reset_index(level=0, inplace=True)
	df_cellname_canonical_parentid.columns = ['cell_name','canonical_parentid']

	df_merge = pd.merge(df_parentid_cellname,df_cellname_canonical_parentid,how='outer',left_on='cell_name',right_on='cell_name',indicator=True)
	df_merge.to_csv(str(prm['output_folder'] / "d_parentid_cellname.csv"),header=True,index=False)

	return

prm = read_parms()
prm['output_folder'].mkdir(parents=True,exist_ok=True)


# Part0: Select the embryo to be enriched
if prm['input_file']:
	input_file_vtp =  prm['input_folder'] / prm['input_file']
else:
	input_file_vtp = ""
	l_files = [(vtp_i,int(re.findall(r'\d+',vtp_i.stem)[0])) for vtp_i in prm['input_folder'].glob('Seeding_cells*[0-9].vtp')]
	for t_vtp in sorted(l_files, key=lambda t: t[1]):  #sorted on file number ascending
		vtp_i, nb_vtp = t_vtp
		if prm['input_file_nb']:
			if prm['input_file_nb'] == int(nb_vtp):
				input_file_vtp = vtp_i
				break
	input_file_vtp = vtp_i

print('vtp file {0} will be used for enriching'.format(input_file_vtp)) if input_file_vtp else print('no vtp file found ! at {0}'.format(input_file_vtp))
poly_embryo = read_vtp_file(input_file_vtp)



#part 2 enrich embryoVTP with cellname -> canonical parentid -> contactidCanonical
#naming_strategy='lineager'  #'micsla':strategy of michiel(needs debugging);  'fallback'=just use naming 'parentXX'; Ultimately this should be replaced with automatic naming
naming_strategy = prm.get('naming_strategy','lineager')
df_canonical_parentid =pd.DataFrame({})
d_parentid_canonical = {}
print(f'Automated naming via {naming_strategy} method:')
if naming_strategy=='micsla':
	a_centroids = np.genfromtxt(prm['input_file_centroids'], delimiter=',')
	err, l_names = identify(a_centroids)
	d_parentid_cellname={ix:x for ix,x in enumerate(l_names)}
	print("\nFound lineage match for {0} cells with score {1}".format(len(l_names), err))
	# a_parentIndex = get_data_array(poly, field_type="CELL", attribute='parentIndex')
	# a_parentIndex = np.where(a_parentIndex>len(l_names),len(l_names)+1,a_parentIndex)
	# l_names.append('Unidentified')
	# a_cellName = np.array([l_names[i] for i in a_parentIndex]).astype(str)
	# poly = add_array(poly,a_cellName, "cellName", field_type="CELL",dtype='str')
	
elif naming_strategy=='fallback':
	a_parentIndex = get_data_array(poly_embryo, field_type="CELL", attribute='parentIndex')
	d_parentid_cellname={i: ("parent" + str(i).zfill(3))for i in np.unique(a_parentIndex)}
	#a_cellName = np.array([("parent" + str(i).zfill(3)) for i in a_parentIndex]).astype(str)
	# poly = add_array(poly,a_cellName, "cellName", field_type="CELL",dtype='str')
elif naming_strategy=='lineager':
	df_cell_names = load_lineage_data(verbose=verbose)
	df_canonical_parentid = pd.read_csv(prm['input_file_canonical_tree'],delimiter="\t")
	d_parentid_cellname = name_cells_using_lineage(poly_embryo,df_cell_names)
add_array_with_mapper(poly_embryo,'parentIndex','cellName',d_parentid_cellname,default='Unidentified')

if not df_canonical_parentid.empty:
	d_cellname_canonical_parentid = {k:get_canonical_parentid(k,df_canonical_parentid) for k in d_parentid_cellname.values()}
	add_array_with_mapper(poly_embryo,'cellName','parentIndexCanonical',d_cellname_canonical_parentid,default=999)
	d_parentid_canonical={i:d_cellname_canonical_parentid[j] for i,j in d_parentid_cellname.items()}
	write_mapper_cellname_csv()
 

# enrich embryo VTP with 2 arrays : contactid and contactid_canonical

a_parentIndex  = get_data_array(poly_embryo, field_type="CELL", attribute='parentIndex')
a_contactIndex  = get_data_array(poly_embryo, field_type="CELL", attribute='contact_index')
a_contact_id = np.array([a_parentIndex[i] for i in a_contactIndex])
a_contact_id = np.where(a_contactIndex==0,a_parentIndex,a_contact_id)  #contactindex=0 means no contact => contact_id = parent_id in that case
poly_embryo = add_array(poly_embryo,a_contact_id, "contactid", field_type="CELL")
if d_parentid_canonical:
	add_array_with_mapper(poly_embryo,'contactid','contactidCanonical',d_parentid_canonical)

#part1 : enrich embryo vtp with contactID and canonical data
# - contact index is the index of the contacting triangle (from the other cell)
# - contactID is the parent index of the contacting cell
#d_parentID_map = {k:k for k in get_parentIDs(poly)} #fallback.  should eventually be replaced with a unique ID for every cell (name)
#poly = enrich_embryo_with_contactID(poly,d_parentID_map)

#part3 write enriched vtp
path_file = str(prm['output_folder'] / (input_file_vtp.stem + "_enriched.vtp")) 
write_vtp_file(poly_embryo, path_file)
if verbose:print("output written : {0}".format(path_file))  

