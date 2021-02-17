"""
Created on 29 may 2020
extract geometry data from an embryo
@author: wimth
"""
import numpy as np
import pandas as pd
import sys, os, re
from param_XML import Param_xml
from geom_classes import Cell

verbose=True
def read_parms():
	param_xml = Param_xml.get_param_xml(sys.argv, l_main_keys=['body', 'extract_geom_data'],verbose=True)  # param file must be passed as first argument

	input_folder = param_xml.get_value('input_folder', ['paths'])
	output_folder = param_xml.get_value('output_folder', ['paths'])

	return input_folder, output_folder


def init_ds_cells():
	for ix, stl_i in enumerate(sorted(input_folder.glob('*.stl'))):
		if verbose: print('Extracting geom data from stl ->  {0}'.format(stl_i.name))
		Cell(stl_i, ic_GT=False , cm=None)
	return 

def write_excel_vol_data():
	d_df = {}
	for cell_i in Cell.l_RES:
		d_df.setdefault('cell_name', []).append(cell_i.name)
		d_df.setdefault('volume', []).append(cell_i.volume)
		d_df.setdefault('surface_area', []).append(cell_i.surface_area)
		d_df.setdefault('sphericity', []).append(cell_i.sphericity)
		l,w,h = cell_i.extents
		d_df.setdefault('extents_length', []).append(l)
		d_df.setdefault('extents_width', []).append(w)
		d_df.setdefault('extents_height', []).append(h)

	df = pd.DataFrame.from_dict(d_df)
	df.to_csv(str(output_folder / "geom_info.csv"), index=False,mode='w')

	if verbose:print("excel with cell volumes written (geom_info.csv)")

	return


# Read parms
input_folder, output_folder = read_parms()
output_folder.mkdir(parents=True, exist_ok=True)
		
init_ds_cells()

write_excel_vol_data()