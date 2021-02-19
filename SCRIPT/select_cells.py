'''
Created on 29 may 2020 
selects the embryo from mpacts PiCS output and extracts the cells (in vtp and stl format).
This selection will be used further on (for SEG validation, geometry extraction etc)
It will also enrich the VTP if a mapping is provided (dds). It will assign a parentId_canonical and contactId_canonical
		
@author: wimth
'''
from VTK_utils import write_vtp_file, write_stl_file, extract_parentID_selection_from_VTP,read_vtp_file,enrich_embryo_with_contactID,get_parentIDs
from helper_functions import dds_get_value
from helper_pandas import excel_to_dds
import numpy as np
import re,os
from param_XML import Param_xml

import os,sys

def read_parms():
	param_xml = Param_xml.get_param_xml(sys.argv,l_main_keys = ['body','select_cells'],verbose=True) #param file must be passed as first argument
	
	input_folder  = param_xml.get_value('input_folder',['paths'])
	input_file_nb = param_xml.get_value('input_file_nb',['paths'])
	input_file = param_xml.get_value('input_file',['paths'])
	output_folder  = param_xml.get_value('output_folder',['paths'])
	write_stl_files  = param_xml.get_value('write_stl_files',['process_flow'])

	dds_file  = param_xml.get_value('dds_file',['body','RAM'],use_main_keys=False)
	repID  = param_xml.get_value('repID',['body','RAM'],use_main_keys=False)
	timestep  = param_xml.get_value('timestep',['body','RAM'],use_main_keys=False)
	if isinstance(timestep,list):
		timestep=timestep[0]

	
	return input_folder, input_file_nb , input_file, output_folder,write_stl_files,dds_file,repID,timestep

def update_mapping_for_cell_division(d_parentID_map,timestep,d_parentID_cellID):

	embID = dds_get_value(dds,'0',repID,"embID")
	d_cellID_parentID = {v:k for k,v in d_parentID_cellID.items()}

	for parentid_i, cellID_i in d_parentID_cellID.items():
		d_3= dds['3'].get((embID,cellID_i))
		if not d_3:
			continue
		if d_3['time_divided'] <= timestep:  #update the map with the daughter cell
			d_parentID_map[d_cellID_parentID[cellID_i]] = dds_get_value(dds,'4',d_3['daughter_cell'],'parentID')

	return d_parentID_map


#parms
input_folder, input_file_nb ,input_file, output_folder,write_stl_files,dds_file,repID,timestep = read_parms()
output_folder.mkdir(parents=True,exist_ok=True)

#pick a VTP from folder to process, (or pick file directly via nm)
input_file_vtp = ""
if input_file:
	input_file_vtp = input_file
else:
	l_files = [(vtp_i,int(re.findall(r'\d+',vtp_i.stem)[0])) for vtp_i in input_folder.glob('Seeding_cells*.vtp')]
	for t_vtp in sorted(l_files, key=lambda t: t[1]):  #sorted on file number ascending
		vtp_i, nb_vtp = t_vtp
		if input_file_nb:
			if int(input_file_nb) == int(nb_vtp):
				input_file_vtp = vtp_i
				break
	input_file_vtp = vtp_i
	
print('vtp file {0} will be used for extraction cells'.format(input_file_vtp),flush=True) if input_file_vtp else print('no vtp file found !',flush=True)
poly_in = read_vtp_file(input_file_vtp)

#enrich embryo-vtp with contactID, and, if mapping is provided, the canonical contactID
if dds_file and repID:
	dds = excel_to_dds(str(dds_file),sheet_name=None,d_index_cols={'1':[0],'3':[0]})
	a_parentIDs = get_parentIDs(poly_in)
	d_parentID_cellID = {k:dds_get_value(dds,'1',(repID,"cell_parent{0}".format(k)),'cellID') for k in a_parentIDs}
	d_parentID_map= {k:int(dds_get_value(dds,'4',d_parentID_cellID[k],'parentID',fallback=99)) for k in a_parentIDs}
	if timestep:  # update mapping to include cell division via tab 3 => embID	cellID	time_divided	daughter_cell
		d_parentID_map = update_mapping_for_cell_division(d_parentID_map, timestep,d_parentID_cellID)

	poly_enriched = enrich_embryo_with_contactID(poly_in,d_parentID_map)
else:
	poly_enriched = enrich_embryo_with_contactID(poly_in)

# write VTP file per cell
d_parentID_VTP = extract_parentID_selection_from_VTP(poly_enriched,verbose=True)

for parent_ID_i,vtp_i in d_parentID_VTP.items():
	write_vtp_file(vtp_i, str(output_folder / "cell_parent{0}.vtp".format(parent_ID_i)) )
	if write_stl_files:
		write_stl_file(vtp_i,str(output_folder / "cell_parent{0}.stl".format(parent_ID_i)),scale=1.0e6)

#write vtp embryo
write_vtp_file(poly_enriched, str(output_folder / "selected_embryo.vtp") )
write_vtp_file(poly_in, str(output_folder / "selected_embryo(not enriched).vtp") )

