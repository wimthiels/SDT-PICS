"""
Created on 29 may 2020
aggregate geometry data from all replicates
@author: wimth
"""
import numpy as np
import pandas as pd
from helper_pandas import excel_to_dds
import sys,re
from param_XML import Param_xml
from shutil import copyfile
from pathlib import Path


verbose=True

d_replicate_embid = {}
dd_replicate_label_cellID = {}
d_tcpart_time = {}

def read_parms():
	param_xml = Param_xml.get_param_xml(sys.argv, l_main_keys=['body', 'agg_geom_data'],verbose=True)  # param file must be passed as first argument

	input_dds = param_xml.get_value('input_dds', ['paths'])
	root_folder = param_xml.get_value('root_folder', ['paths'])
	regex_file_select = param_xml.get_value('regex_file_select', ['paths'])
	regex_file_select2 = param_xml.get_value('regex_file_select2', ['paths'])
	regex_file_select3 = param_xml.get_value('regex_file_select3', ['paths'])
	regex_file_select4 = param_xml.get_value('regex_file_select4', ['paths'])
	regex_file_select5 = param_xml.get_value('regex_file_select5', ['paths'])

	output_folder = param_xml.get_value('output_folder', ['paths'])

	l_select_cellID = param_xml.get_value('l_select_cellID')
	copy_and_rename_only = param_xml.get_value('copy_and_rename_only')
	prefix_folder = param_xml.get_value('prefix_folder')

	return input_dds, root_folder, regex_file_select, regex_file_select2, regex_file_select3, l_select_cellID,regex_file_select4, regex_file_select5,output_folder,copy_and_rename_only,prefix_folder

def get_replicate_and_tcpart(p_file):
	pattern = re.compile(r'(?:{})([0-9][0-9])([0-9][0-9])(?:/)'.format(prefix_folder))
	match=pattern.search(p_file)
	replicate = match.group(1)
	tcpart = match.group(2)
	#if verbose:print("{0} with replicate {1} and tcpart {2}".format(p_file,replicate,tcpart))

	return int(replicate), int(tcpart) # dds keeps reformatting to int...

def dds_get_value(label_ds,k1,k2):
	try:
		return dds[label_ds][k1][k2]
	except KeyError as ex:
		return 'not in dds : {0}'.format(k1)



def construct_vol_df(file_i):
	df_file_in = pd.read_csv(file_i,header=0,index_col=False,keep_default_na=True) 
	nb_rows, nb_columns = df_file_in.shape

	replicate,tcpart = get_replicate_and_tcpart(str(file_i))

	d_df = {}
	d_df['embID']  = [dds_get_value('0',replicate,'embID')] * nb_rows
	d_df['repID']  = [dds_get_value('0',replicate,'repID')] * nb_rows
	d_df['cellID'] = [dds_get_value('1',(replicate,label_i),'cellID') for label_i in df_file_in['name']]
	d_df['t'] = [dds_get_value('2',tcpart,'time')] * nb_rows
	d_df['vol'] = df_file_in['volume']
	d_df['surface_area'] = df_file_in['surface_area']
	d_df['sphericity'] = df_file_in['sphericity']
	d_df['extents_length'] = df_file_in['extents_length']
	d_df['extents_width'] = df_file_in['extents_width']
	d_df['extents_height'] = df_file_in['extents_height']
	d_df['nb_pix_last_z']  = [0] * nb_rows  #will be updated later
	d_df['orientation']  = [dds_get_value('5',dds_get_value('0',replicate,'repID'),'orientation')] * nb_rows
	d_df['delta_t_division'] = [int(dds_get_value('1',(replicate,label_i),'delta_t_division')) + int(dds_get_value('2',tcpart,'time')) - 1 for label_i in df_file_in['name']]
	[dds_get_value('1',(replicate,label_i),'delta_t_division')for label_i in df_file_in['name']]
	return pd.DataFrame.from_dict(d_df)

def update_nb_pix_last_z(file_i, df_agg_data_vol):
	df_file_in = pd.read_excel(file_i,sheet_name='nb_pix',header=0,index_col=False,keep_default_na=True) 
	replicate,tcpart = get_replicate_and_tcpart(str(file_i))
	cond1 = (df_agg_data_vol['repID'] == dds_get_value('0',replicate,'repID'))
	cond2 = (df_agg_data_vol['t'] == dds_get_value('2',tcpart,'time'))
	for _,row_i in df_file_in.iterrows():
		label_i = "cell_parent{0}".format(row_i['parentID'])
		cond3 = (df_agg_data_vol['cellID'] == dds_get_value('1',(replicate,label_i),'cellID'))
		df_agg_data_vol.loc[cond1 & cond2 & cond3,'nb_pix_last_z'] = row_i['nb_pixels']

	return

def get_ix_z_last_slice(file_i):
	df_file_in = pd.read_excel(file_i,sheet_name='ix_z',header=0,index_col=False,keep_default_na=True) 
	replicate,tcpart = get_replicate_and_tcpart(str(file_i))

	d_df = {}
	d_df['repID']  = [dds_get_value('0',replicate,'repID')] 
	d_df['t'] = [dds_get_value('2',tcpart,'time')] 
	d_df['ix_z_last_slice'] = df_file_in['ix_z']

	return pd.DataFrame.from_dict(d_df)

def construct_contact_df(file_i):
	df_file_in = pd.read_csv(file_i,header=0,index_col=False,keep_default_na=True) #cell1 (parentindex) , cell2, contactarea
	nb_rows, nb_columns = df_file_in.shape

	replicate,tcpart = get_replicate_and_tcpart(str(file_i))

	d_df = {}
	d_df['embID']  = [dds_get_value('0',replicate,'embID')] * nb_rows
	d_df['repID']  = [dds_get_value('0',replicate,'repID')] * nb_rows
	d_df['t'] = [dds_get_value('2',tcpart,'time')] * nb_rows
	d_df['cellID_1'] = [dds_get_value('1',(replicate,"cell_parent{0}".format(label_i)),'cellID') for label_i in df_file_in['cell1']]
	d_df['cellID_2'] = [dds_get_value('1',(replicate,"cell_parent{0}".format(label_i)),'cellID')  for label_i in df_file_in['cell2']]
	
	d_df['contactarea'] = df_file_in['contactarea']

	return pd.DataFrame.from_dict(d_df)

def construct_indent_volume_df(file_i):
	df_file_in = pd.read_excel(file_i,sheet_name='volume_cut',header=0,index_col=False,keep_default_na=True) #cell1 (parentindex) , cell2, contactarea
	nb_rows, nb_columns = df_file_in.shape

	replicate,tcpart = get_replicate_and_tcpart(str(file_i))

	d_df = {}
	d_df['embID']  = [dds_get_value('0',replicate,'embID')] * nb_rows
	d_df['repID']  = [dds_get_value('0',replicate,'repID')] * nb_rows
	d_df['t'] = [dds_get_value('2',tcpart,'time')] * nb_rows

	return pd.concat([pd.DataFrame(d_df),df_file_in],axis=1)

def blanc_cell_volumes_after_division(df_agg_data):

	for t_embID_cellID, d_value in dds['3'].items():
		cond1 = (df_agg_data.embID==t_embID_cellID[0])
		cond2 = (df_agg_data.cellID==t_embID_cellID[1])
		cond3 = (df_agg_data.t>=d_value['time_divided'])

		if not d_value['daughter_cell']:  #if no daughter cell is given, all data after cell division is blanked out
			df_agg_data.loc[cond1 & cond2 & cond3, 'vol'] = 'NaN'
		else:                             # replace parent cell with daughter label
			df_agg_data.loc[cond1 & cond2 & cond3, 'cellID'] = d_value['daughter_cell']


	return


def blanc_contactarea_after_division(df_agg_data):

	for t_embID_cellID, d_value in dds['3'].items():
		cond1 = (df_agg_data.embID==t_embID_cellID[0])
		cond2 = (df_agg_data.cellID_1==t_embID_cellID[1])
		cond3 = (df_agg_data.cellID_2==t_embID_cellID[1])
		cond4 = (df_agg_data.t>=d_value['time_divided'])

		if not d_value['daughter_cell']:  #if no daughter cell is given, all data after cell division is blanked out
			df_agg_data.loc[cond1 & (cond2 | cond3) & cond4, 'contactarea'] = 'NaN'
		else:                             # replace parent cell with daughter label
			df_agg_data.loc[cond1 & cond2 & cond4, 'cellID_1'] = d_value['daughter_cell']
			df_agg_data.loc[cond1 & cond3 & cond4, 'cellID_2'] = d_value['daughter_cell']

	return


def copy_selected_cells(file_i, replicate,tcpart):
	p_folder= file_i.parent
	for vtp_i in p_folder.glob('cell_parent*.vtp'):
		cellID = dds_get_value('1',(replicate,vtp_i.stem),'cellID') 
		if cellID in l_select_cellID:
			vtp_o = output_folder / "selected_embryos"/ "{0}_{1}_t{2}.vtp".format(cellID,dds_get_value('0',replicate,'repID'),dds_get_value('2',tcpart,'time'))
			copyfile(vtp_i, vtp_o)
			if verbose:print('{0} is copied'.format(vtp_i))

	return

#parms
input_dds, root_folder, regex_file_select ,regex_file_select2, regex_file_select3, l_select_cellID,regex_file_select4, regex_file_select5,output_folder,copy_and_rename_only, prefix_folder= read_parms()
output_folder.mkdir(parents=True,exist_ok=True)

#STEP 1 : read in datastructures for mapping
dds = excel_to_dds(str(input_dds),sheet_name=None,d_index_cols={'1':[0],'3':[0]})
#dds[0] : d_replicate_embID
#dds[1] : dd_replicate_cell-label_cellID  #multi-index because of d_index_cols={'1':[0]}
#dds[2] : d_tcpart_time
#dds[3] : d_embID_cellID_time-division : eg TL1  P2 3 = P2 is fully divided at t=3 #multi_index




#STEP0 : copy VTP files for generating a 4D movie in paraview
(output_folder / "selected_embryos").mkdir(parents=True,exist_ok=True)
for file_i in sorted(root_folder.glob(regex_file_select3)):
	try:
		replicate,tcpart = get_replicate_and_tcpart(str(file_i))
		file_o = output_folder / "selected_embryos"/ "selected_embryo_{0}_t{1}.vtp".format(dds_get_value('0',replicate,'repID'),dds_get_value('2',tcpart,'time'))
		copyfile(file_i, file_o)
		copy_selected_cells(file_i,replicate,tcpart)
		if verbose:print('{0} is copied'.format(file_i))
	except:
		print (file_i," is skipped ! check the file...")

if not copy_and_rename_only:

	writer = pd.ExcelWriter(output_folder/"agg_geom_data.xlsx")


	#STEP 2 : aggregate volume data
	l_df = []
	for file_i in sorted(root_folder.glob(regex_file_select)):
		try:
			l_df.append(construct_vol_df(file_i))
			if verbose:print('{0} is appended'.format(file_i))
		except:
			print (file_i," is skipped ! check the file...")

	df_agg_data_vol = pd.concat(l_df,ignore_index=True)
	blanc_cell_volumes_after_division(df_agg_data_vol)


	#STEP 3 : aggregate contact area data
	l_df = []
	for file_i in sorted(root_folder.glob(regex_file_select2)):
		try:
			l_df.append(construct_contact_df(file_i))
			if verbose:print('{0} is appended'.format(file_i))
		except:
			print (file_i," is skipped ! check the file...")

	df_agg_data_contactarea = pd.concat(l_df,ignore_index=True)
	blanc_contactarea_after_division(df_agg_data_contactarea)


	#STEP 4 : update last slice info on volume tab
	l_df = []
	for file_i in sorted(root_folder.glob(regex_file_select4)):
		try:
			update_nb_pix_last_z(file_i, df_agg_data_vol)
			l_df.append(get_ix_z_last_slice(file_i))
			if verbose:print('{0} :  number of pixels last Z is updated in the volume tab'.format(file_i))
		except:
			print (file_i," is skipped ! check the file...")
	df_ix_z_last_slice = pd.concat(l_df,ignore_index=True)

	#STEP 5 : Aggregate indentation volume data
	l_df = []
	for file_i in sorted(root_folder.glob(regex_file_select5)):
		try:
			l_df.append(construct_indent_volume_df(file_i))
			if verbose:print('{0} is appended'.format(file_i))
		except:
			print (file_i," is skipped ! check the file...")
	df_indent_volume = pd.concat(l_df,ignore_index=True)

	#write output excel
	with pd.ExcelWriter(output_folder/"agg_geom_data.xlsx") as writer:
	    df_agg_data_vol.to_excel(writer, sheet_name='cell_geometry',index=False)
	    df_agg_data_contactarea.to_excel(writer, sheet_name='contactarea',index=False)
	    df_ix_z_last_slice.to_excel(writer, sheet_name='ix_z_last_slice',index=False)
	    df_indent_volume.to_excel(writer, sheet_name='indent_volume',index=False)
	    if verbose:print("excel output : {0}".format(str(output_folder/"agg_geom_data.xlsx")))


