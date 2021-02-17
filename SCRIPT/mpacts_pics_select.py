'''
Created on 17jan 2021
selects the embryo from mpacts PiCS output and extracts the cells (in vtp and stl format).  stl's are scaled
This selection will be used further on (for SEG validation, geometry extraction etc)
		
@author: wimth
'''
from VTK_utils import write_vtp_file, write_stl_file, extract_parentID_selection_from_VTP,read_vtp_file,get_parentIDs,list_unique_values
import numpy as np
import re,os
from param_XML import Param_xml
import os,sys

def read_parms():
	param_xml = Param_xml.get_param_xml(sys.argv, l_main_keys=['body'],verbose=True) 
	param_xml.add_prms_to_dict(l_keys=['mpacts_pics_select'])
	param_xml.repr_prm()
	return param_xml.prm

prm = read_parms()
prm['output_folder'].mkdir(parents=True,exist_ok=True)

#pick a VTP to extract cells
input_file_vtp = ""
if prm['input_file']:
	input_file_vtp = input_file
else:
	input_file_vtp =  [vtp for vtp in prm['input_folder'].glob('*enriched.vtp')]
	if input_file_vtp:
		input_file_vtp = input_file_vtp[0]

print('vtp file {0} will be used for extraction cells'.format(input_file_vtp),flush=True) if input_file_vtp else print('no vtp file found !',flush=True)
poly_in = read_vtp_file(input_file_vtp)

# write VTP file per cell
d_parentID_VTP = extract_parentID_selection_from_VTP(poly_in,verbose=True)
for parent_ID_i,vtp_i in d_parentID_VTP.items():

	l_cell_name = list_unique_values(vtp_i,attribute='cellName', field_type="CELL",repr=True)
	if l_cell_name:
		cell_name = l_cell_name[0]
	else:
		cell_name = f"cell_parent{str(parent_ID_i).zfill(3)}" #fallback naming via parentid

	write_vtp_file(vtp_i, str(prm['output_folder'] / f"{cell_name}.vtp"))
	if prm['write_stl_files']:
		write_stl_file(vtp_i,str(prm['output_folder'] / f"{cell_name}.stl"),scale=1.0e6)
#write vtp embryo
write_vtp_file(poly_in, str(prm['output_folder'] / "selected_embryo.vtp") )