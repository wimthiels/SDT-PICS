"""
Created on 16 September 2020
extract indent volume
@author: wimth

datastructures:
d_cells:dict
--->cell:dict;(var = d_cell)
	--->'poly':poly,(poly of the input file)
		'parent_id':int
		'trimesh:Trimesh
		'contacts':dict;(var = d_contacts)
		--->parentid:dict; (var = d_contact ; parentid of contacting cell)
			--->'contact_area':poly; vtp of contact area of input file with contacting cell
				'border_points':poly; vtp of border points of contact area
				'cutting_plane':poly, vtp of the plane that does the cutting (normal plane, not inverted)
				'cut_plane': Trimesh; (trimesh of the cut with the normal plane)
				'cut_plane_inverted_normal': Trimesh; (trimesh of the cut with the inverted plane)
				'type_cut_plane': str; ('body' or 'indent' (added later, needed to give correct names to outputted polys))
				'type_cut_plane_inverted_normal':str; ('body' or 'indent' (added later, needed to give correct names to outputted polys))
				'pair_ctr':int;(a simple counter per cell-pair, for easy post processing and sorting)

"""
import numpy as np
import pandas as pd
import sys, os, re
import pprint;pp = pprint.PrettyPrinter(indent=4);ppr = pp.pprint
import copy

import trimesh
from pymeshfix import _meshfix
from param_XML import Param_xml


import vtk
from VTK_utils import extract_selection_vtp, read_vtp_file, write_vtp_file, \
				get_outside_points, get_unique_values,construct_plane, \
				write_stl_file,scale_poly,poly_to_trimesh
from vtk.util.numpy_support import vtk_to_numpy
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs



verbose=True
def read_parms():
	param_xml = Param_xml.get_param_xml(sys.argv, l_main_keys=['body', 'extract_indent_volume'],verbose=True)  # param file must be passed as first argument

	input_folder = param_xml.get_value('input_folder', ['paths'])
	l_input_vtp = param_xml.get_value('l_input_vtp', ['paths']) # this list must be thought of as a list of pairs
	output_folder = param_xml.get_value('output_folder', ['paths'])

	return l_input_vtp, input_folder,output_folder

def load_file_and_build_dict(p_cell,d_cells):

	file_error = False
	try:
		poly_in = read_vtp_file(p_cell)
		parentid = get_unique_values(poly_in,'parentIndexCanonical')[0] 
		poly_in = scale_poly(poly_in,1.e6) #trimesh fails in micron
		f_stl = str(output_folder / (p_cell.stem + ".stl"))
		write_stl_file(poly_in, f_stl)

		trimesh_in  = trimesh.load(f_stl, file_type="stl")
		if trimesh_in.is_volume:
			print("volume of input mesh {1} = {0}".format(trimesh_in.volume,p_cell.stem))

		d_cells[p_cell.stem] = {'poly':poly_in,'parentid':parentid,'trimesh':trimesh_in,'contacts':{}}  #d_contacts
	except:
		print('Error reading input file ({}), continuing with next cell-pair'.format(p_cell))
		file_error = True

	return file_error

def mesh_cut_with_plane(cell, parentid_contact, normal_plane, point_plane,invert_normal=False):

	if invert_normal:
		normal_plane *= -1

	trimesh_slice = trimesh.intersections.slice_mesh_plane(d_cells[cell]['trimesh'],normal_plane,point_plane,cap=True)

	#write and repair
	outputname = str(output_folder / '<temp>_{0}_cut({1}->{2})(plane_inverted={3}).stl'.format(cell,d_cells[cell]['parentid'],parentid_contact,invert_normal))
	_ = trimesh_slice.export(outputname)
	try:
		_meshfix.clean_from_file(outputname, outputname) #trimesh output is not watertight, so repair and reload from file
		trimesh_slice = trimesh.load(outputname, file_type="stl")
	except:
		print('_meshfix.clean_from_file crashed on {0}'.format(outputname))

	if not trimesh_slice.is_volume:
		print('warning: volume={0} is not watertight. volume measurement is unreliable'.format(outputname))
	else:
		if invert_normal:
			d_cells[cell]['contacts'][parentid_contact]['cut_plane'] = trimesh_slice
		else:
			d_cells[cell]['contacts'][parentid_contact]['cut_plane_inverted_normal'] = trimesh_slice

	return

def compose_df():
	""" Compose the dataframe based on d_cells
	per cell fill go over its contacts  
	each cell can have multiple contacts, and for each contact, 2 cuts are made, so 2 rows are added per contact"""

	l_d = []
	for cell, d_cell in d_cells.items():

		d_row_cell ={'input_mesh':cell,'parentid':d_cell['parentid']}

		for parent_id,d_contact in d_cell['contacts'].items():
			d_row = copy.deepcopy(d_row_cell)
			d_row['contactidCanonical']= parent_id
			d_row['pair_ctr']=d_contact['pair_ctr']
			if 'contact_area' in d_contact:
				d_row['contact_area']=poly_to_trimesh(d_contact['contact_area']).area
				d_row['invert_normal']='False'
				d_row['volume_cut(um³)']=d_contact['cut_plane'].volume
				d_row['type_cut']=d_contact['type_cut_plane']
				l_d.append(d_row)
				d_row_invert = copy.deepcopy(d_row)
				d_row_invert['invert_normal']='True'
				d_row_invert['volume_cut(um³)']=d_contact['cut_plane_inverted_normal'].volume
				d_row_invert['type_cut']=d_contact['type_cut_plane_inverted_normal']
				l_d.append(d_row_invert)
			else:  # no contact
				d_row['contact_area']=0
				d_row['invert_normal']='NA'
				d_row['volume_cut(um³)']=0
				d_row['type_cut']='indent'
				l_d.append(d_row)
				d_row_body = copy.deepcopy(d_row)
				d_row_body['volume_cut(um³)']=d_cell['trimesh'].volume
				d_row_body['type_cut']='body'
				l_d.append(d_row_body)

	return pd.DataFrame(l_d).sort_values(by=['pair_ctr','input_mesh'])


def add_type_cut(d_cells):
	""" go over d_cells and determine which of the two cuts is the body and which one is the indent volume.
	The smallest piece is always considered the indent"""

	for input_file_vtp,d_cell in d_cells.items():
		for _, d_contact in d_cell['contacts'].items():
			if 'cut_plane' in d_contact:
				if d_contact['cut_plane'] and d_contact['cut_plane_inverted_normal']:
					if d_contact['cut_plane'].volume > d_contact['cut_plane_inverted_normal'].volume:
						d_contact['type_cut_plane'] = 'body'
						d_contact['type_cut_plane_inverted_normal'] = 'indent'
					else:
						d_contact['type_cut_plane'] = 'indent'
						d_contact['type_cut_plane_inverted_normal'] = 'body'
			else:
				d_contact['type_cut_plane'] = 'NA'
				d_contact['type_cut_plane_inverted_normal'] = 'NA'


	return


def write_polys():
	print('... Writing polys')
	for cell, d_cell in d_cells.items():
		write_vtp_file(d_cell['poly'], str(output_folder / cell))
		for parentid_contact,d_contact in d_cell['contacts'].items():
			if d_contact.get('contact_area') is None:
				continue
			write_vtp_file(d_contact['contact_area'] , str(output_folder / '{0}_contact_area({1}->{2}).vtp'.format(cell,d_cell['parentid'],parentid_contact)))
			write_vtp_file(d_contact['cutting_plane'], str(output_folder / '{0}_cutting_plane({1}->{2}).vtp'.format(cell,d_cell['parentid'],parentid_contact)))
			write_vtp_file(d_contact['border_points'], str(output_folder / '{0}_border_points({1}->{2}).vtp'.format(cell,d_cell['parentid'],parentid_contact)))

			d_contact['cut_plane'].export(file_obj=str(output_folder / '{0}_cut_{3} ({1}->{2}).stl'.format(cell,d_cell['parentid'],parentid_contact,d_contact['type_cut_plane'])),
				file_type='stl')
			d_contact['cut_plane_inverted_normal'].export(file_obj=str(output_folder / '{0}_cut_{3} ({1}->{2}).stl'.format(cell,d_cell['parentid'],parentid_contact,d_contact['type_cut_plane_inverted_normal'])),
				file_type='stl')

	return


# Read parms
l_input_vtp, input_folder,output_folder = read_parms()
output_folder.mkdir(parents=True, exist_ok=True)

#Part 1 : build datastructure d_cells
d_cells = {}
for ix_cell in range (0, len(l_input_vtp),2):  #process in groups of two:
	input_vtp_1 = l_input_vtp[ix_cell]
	input_vtp_2 = l_input_vtp[ix_cell + 1]
	cell1 = (input_folder / input_vtp_1).stem
	cell2 = (input_folder / input_vtp_2).stem

	if cell1 not in d_cells:
		file_error = load_file_and_build_dict(input_folder / input_vtp_1,d_cells)
		if file_error:
			continue

	if cell2 not in d_cells:
		file_error = load_file_and_build_dict(input_folder / input_vtp_2,d_cells)
		if file_error:
			continue

	# get borderpoints (from both contactareas)
	a_border_points_total = np. empty((0,3))
	flag_no_overlap = False
	l_iter = [(cell1, d_cells[cell2]['parentid']),(cell2, d_cells[cell1]['parentid'])]
	for cell, parentid_contact in l_iter:
		d_cells[cell]['contacts'][parentid_contact] = {}  #d_contact
		d_cells[cell]['contacts'][parentid_contact]['pair_ctr'] = int(ix_cell/2)
		poly_contact_area = extract_selection_vtp(d_cells[cell]['poly'], query="contactidCanonical == {0}".format(parentid_contact) ,field_type="CELL")
		if not poly_contact_area:
			print('no overlap between cells detected.  Indent volume = 0')
			flag_no_overlap = True
			continue

		a_border_points,poly_border = get_outside_points(poly_contact_area)

		d_cells[cell]['contacts'][parentid_contact]['border_points'] = poly_border
		d_cells[cell]['contacts'][parentid_contact]['contact_area'] = poly_contact_area
		
		a_border_points_total = np.append(a_border_points_total,a_border_points,axis=0)

	if flag_no_overlap:
		continue
	else:
		# fit plane through points (=border points from both contactareas)
		point_plane,normal_plane = trimesh.points.plane_fit(a_border_points_total)
		poly_cutting_plane = construct_plane(point_plane,normal_plane)
		# cut meshes with cutting plane (2 cuts per cell)
		for cell, parentid_contact in l_iter:
			d_cells[cell]['contacts'][parentid_contact]['cutting_plane'] = poly_cutting_plane
			for invert_normal in [False,True]:
				mesh_cut_with_plane(cell,parentid_contact, normal_plane, point_plane,invert_normal=invert_normal)

	print('end cell1:{0}, cell 2:{1}'.format(cell1,cell2))
	ppr(d_cells)
add_type_cut(d_cells)

#Part 2 : write output
df = compose_df()
with pd.ExcelWriter(output_folder/"volume_mesh_plane_intersections.xlsx") as writer:
	df.to_excel(writer, sheet_name='volume_cut',index=False)
write_polys()

