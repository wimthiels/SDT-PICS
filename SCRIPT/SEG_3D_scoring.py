"""
Created on 11dec 2019
given 2 sets of stl's, 1 for the GT the other for the RES, calculate the 3D SEG score for the segmentation

CAVEAT : trimesh boolean operations do NOT work in micron !! so scale appropriate
@author: wimth
"""

import trimesh
from trimesh.collision import CollisionManager
# from trimesh.parent.Geometry import apply_scale
import numpy as np
import pandas as pd
import sys, os, re
from collections import Counter
from VTK_utils import write_stl_file, read_vtp_file, read_stl_file
from param_XML import Param_xml
from geom_classes import Cell, Contactinfo

verbose = True


def read_parms():
	param_xml = Param_xml.get_param_xml(sys.argv, l_main_keys=['body', 'SEG_3D_scoring'],verbose=True)  # param file must be passed as first argument

	input_folder_GT = param_xml.get_value('input_folder_GT', ['paths'])
	input_folder_RES = param_xml.get_value('input_folder_RES', ['paths'])
	output_folder = param_xml.get_value('output_folder', ['paths'])
	input_file_type = param_xml.get_value('input_file_type',['paths'])

	data = param_xml.get_value('l_error_intersection_trimesh')
	l_error_intersection_trimesh = []
	for i,k in zip(data[0::2], data[1::2]):
		l_error_intersection_trimesh.append(tuple((i,k)))

	scaling =  param_xml.get_value('scaling')

	return input_folder_GT, input_folder_RES, output_folder, l_error_intersection_trimesh, input_file_type, scaling

def init_output_folder():
	file_SEG3D.unlink() if file_SEG3D.exists() else print('') # not delete the whole dir (prevent delete disasters) 
	folder_intersection = output_folder / "intersection"
	folder_intersection.mkdir(parents=True, exist_ok=True)
	folder_RES_STL = output_folder / "RES_STL"
	folder_RES_STL.mkdir(parents=True, exist_ok=True)
	return folder_intersection,folder_RES_STL


def init_ds_cells(folder_cells, ic_GT):

	if input_file_type=="vtp":
		pattern  = "*.stl" if ic_GT else "cell*.stl"
	else:
		pattern  = "*.stl" 

	for ix, file_i in enumerate(sorted(folder_cells.glob(pattern))):
		if verbose: print('Adding {} ->  {}'.format("stl",file_i.name))
		cm = cm_GT if ic_GT else cm_RES
		Cell(file_i, ic_GT, cm,"stl")
	return

def write_excel_RES_data():
	d_df = {}
	for cell_i in Cell.l_RES:
		d_df.setdefault('RES_cell', []).append(cell_i.name)
		d_df.setdefault('volume_RES', []).append(cell_i.volume)

	df = pd.DataFrame.from_dict(d_df)
	df.to_csv(str(output_folder / "RES_geometry_details.csv"), index=False)

	print("RES_geometry_details.csv written")

	return

def write_excel_overlap_data():
	d_df = {}
	SEG3D_counter = 0
	volume_union_counter = 0
	volume_sum_counter = 0
	volume_intersection_counter = 0
	for cell_i in Cell.l_GT:
		d_df.setdefault('GT_cell', []).append(cell_i.name)
		cell_max_overlap, volume_max_overlap = cell_i.get_max_overlapping_cell_with_volume(based_on_SEG_score=True)
		if cell_max_overlap:
			d_df.setdefault('cell_max_overlap', []).append(cell_max_overlap.name)
			d_df.setdefault('volume_RES(mu_m3)', []).append(cell_max_overlap.volume)
			union_volume = cell_i.volume + cell_max_overlap.volume - volume_max_overlap
			sum_volume = cell_i.volume + cell_max_overlap.volume
		else:
			d_df.setdefault('cell_max_overlap', []).append("None")
			d_df.setdefault('volume_RES(mu_m3)', []).append(0)
			union_volume = cell_i.volume
			sum_volume = cell_i.volume

		d_df.setdefault('volume_intersect_vol', []).append(volume_max_overlap)
		volume_intersection_counter += volume_max_overlap
		d_df.setdefault('volume_GT(mu_m3)', []).append(cell_i.volume)
		# d_df.setdefault('volume_GT_pyvista', []).append(cell_i.volume)
		d_df.setdefault('volume_union', []).append(union_volume)
		volume_union_counter += union_volume
		volume_sum_counter += sum_volume
		if union_volume != 0:
			SEG3D = volume_max_overlap / union_volume
			dice = (2 * volume_max_overlap)/sum_volume
			d_df.setdefault('SEG3D', []).append(SEG3D)
			d_df.setdefault('dice', []).append(dice)
			SEG3D_counter += SEG3D
		else:
			print("warning : union volume supposedly zero. should not happen...{0}".format(cell_i.name))
			d_df.setdefault('SEG3D', []).append(0)
			d_df.setdefault('dice', []).append(0)

	# add aggregated extra row
	d_df.setdefault('GT_cell', []).append("SUM")
	d_df.setdefault('cell_max_overlap', []).append("N/A")
	d_df.setdefault('volume_intersect_vol', []).append(volume_intersection_counter)
	d_df.setdefault('volume_GT(mu_m3)', []).append("")
	# d_df.setdefault('volume_GT_pyvista', []).append("")
	d_df.setdefault('volume_RES(mu_m3)', []).append("")
	d_df.setdefault('volume_union', []).append(volume_union_counter)
	SEG3D = volume_intersection_counter / volume_union_counter
	d_df.setdefault('SEG3D', []).append(SEG3D)
	dice = (2* volume_intersection_counter)/ volume_sum_counter
	d_df.setdefault('dice', []).append(dice)

	df = pd.DataFrame.from_dict(d_df)
	df.to_csv(file_SEG3D, index=False)

	print("FINAL SEG3D score (averaged) of {0}".format(SEG3D_counter / len(Cell.l_GT)))

	return



# MAIN------------------------------------------------------------------------------------------------------------------------

# Read parms
print('reading parms')
input_folder_GT, input_folder_RES, output_folder,l_error_intersection_trimesh, input_file_type, scaling = read_parms()
file_SEG3D = output_folder / "SEG3D_score_details.csv"
folder_intersection,folder_RES_STL = init_output_folder()

if input_file_type=='vtp':

	# Part1 : create STL from RES
	#	trimesh cannot work with vtp directly, so convert to stl first  # print(trimesh.exchange.load.available_formats())
	for ix, vtp_path_i in enumerate(sorted(input_folder_RES.glob('cell*.vtp'))):
		poly = read_vtp_file(vtp_path_i)
		write_stl_file(poly, (folder_RES_STL / (vtp_path_i.stem + ".stl")))

	#	the STL's from ground truth are not scaled to micron (which is good, because trimesh fails in micron),  but now scale RES appropriately
	
	for ix, p_stl_i in enumerate(sorted(folder_RES_STL.glob('cell*.stl'))):
		tmh = trimesh.load(p_stl_i, file_type="stl")
		tmh.apply_scale(scaling)
		stl_ascii = trimesh.exchange.stl.export_stl_ascii(tmh)  # tmh.export writes a binary stl by default , problems with reloading this (so pass via file object (str))
		f = open(str(p_stl_i), 'w')
		f.write(stl_ascii)

elif input_file_type=='ply':  #convert to stl and scale
	for ix, ply_path_i in enumerate(sorted(input_folder_RES.glob('*.ply'))):
		tmh = trimesh.load(ply_path_i, file_type="ply")
		tmh.apply_scale(scaling)
		stl_ascii = trimesh.exchange.stl.export_stl_ascii(tmh)  # tmh.export writes a binary stl by default , problems with reloading this (so pass via file object (str))
		p_stl = folder_RES_STL / (ply_path_i.stem + ".stl")
		f = open(str(p_stl), 'w')
		f.write(stl_ascii)

print("Part2a: make 2 collision managers, one for GT , one for RES.  and build up datastructures of cells")
cm_GT = CollisionManager()
cm_RES = CollisionManager()


init_ds_cells(folder_RES_STL, ic_GT=False)
init_ds_cells(input_folder_GT, ic_GT=True)

print("# Part2b: already compose excel with the RES geometric data")
write_excel_RES_data()


# Part3 : find max contacting volume for each GT cell (always start from GT cells = similar to SEG scoring ).
#    1 : Use the collision detection to detect which cells overlap. 
#    2 : check the boolean difference volume.  
#    3 : pick the maximum for scoring

print('# Part3 : Checking collision pairs via trimesh...')
is_collision, s_t_colliding_pairs = cm_GT.in_collision_other(cm_RES, return_names=True, return_data=False)
if not is_collision:
	print("No collision between the cells, exitting")
	exit()

	
if verbose:print(s_t_colliding_pairs)
for t_colliding_pair in s_t_colliding_pairs:
	# print('--Checking collision pair ->  {0}'.format(t_colliding_pair))
	name_cell_1, name_cell_2 = t_colliding_pair
	cell_1 = Cell.d_name_cell[name_cell_1]
	cell_2 = Cell.d_name_cell[name_cell_2]
	print("start calculating intersection : {0}".format(t_colliding_pair),flush=True)
	trimesh_intersection = cell_1.trimesh.intersection(cell_2.trimesh)
	
	# diff.export(str(mesh_path.with_name(mesh_path.stem + '_no_overlap.stl')))  #export to same folder
	f_name = '{0}_{1}_intersection.stl'.format(cell_1.name, cell_2.name)
	f_path = str(folder_intersection / f_name)
	if t_colliding_pair in l_error_intersection_trimesh:
		print('{0} causes trimesh to crash, so is skipped (l_error_intersection_trimesh,)'.format(t_colliding_pair),flush=True)
		continue
	if trimesh_intersection.is_volume:
		cell_1.d_cell_contactinfo[cell_2] = Contactinfo(trimesh_intersection.volume)
	else:
		#Fallback strategy : random volume sampling
		print("intersection was not considered a volume. Fallback applied through volume sampling : ",flush=True)
		nb_sample_points=500  ##1000 is standard
		a_points_cell1 = trimesh.sample.volume_mesh(cell_1.trimesh, nb_sample_points)
		print("-->points sampled. calculating ratio inside<->outside..",flush=True)
		a_sign_dist = trimesh.proximity.signed_distance(cell_2.trimesh, a_points_cell1) #Points OUTSIDE the mesh will have NEGATIVE distance
		nb_inside = np.sum(a_sign_dist>=0,axis=0)
		score = nb_inside/a_sign_dist.shape[0]

		cell_1.d_cell_contactinfo[cell_2] = Contactinfo(score * cell_1.trimesh.volume)
		print("-->{0} points were sampled from cell_1, and {1} of those points fall inside cell_2 (= {3}%). So the intersection volume is estimated at {2}".
			format(a_sign_dist.shape[0],nb_inside,score * cell_1.trimesh.volume,score),flush=True)
		# if trimesh_union.is_volume:
		# 	print ("intersection did not work, but the union could be calculated and equals {0}".format(trimesh_union.volume))
		# else:
		# 	print("warning : boolean intersection between {0} and {1} is not considered as a volume by trimesh, so no overlap volume can be determined. Union volume also failed.".format(cell_1.name, cell_2.name))

#Part 3 extension : trimesh collision does NOT detect overlap if 1 mesh is fully contained in the other
# therefore we check if the centroid of a cell is contained in the other.  
# if this cell was not yet detected above, this mesh is considered to be fully contained in the other
l_centroid_RES = []
for cell_i in Cell.l_RES:
	l_centroid_RES.append(cell_i.centroid)

if any(x is None for x in l_centroid_RES):
	print(f"l_centroid_RES={l_centroid_RES}")
	print(f"The centroids of some cells were defined as None, this usually means that a mesh is rubbish (due to segmentation instabilities often).  Visually check the meshes.  Aborting the SEG scoring !")
	sys.exit()

print("l_centroid_RES=",l_centroid_RES)
for cell_i in Cell.l_GT:
	l_bool_contains = cell_i.trimesh.contains(np.array(l_centroid_RES))
	a_ix_contains = np.where(l_bool_contains)
	for ix_res in a_ix_contains[0].tolist():
		cell_RES = Cell.l_RES[ix_res]
		# print ("the centroid of RES-cell {0} is contained in the GT-cell {1} ".format(cell_RES.name, cell_i.name))
		if not cell_i.d_cell_contactinfo.get(cell_RES):
			print ("the centroid of RES-cell {0} is contained in the GT-cell {1} AND it was not yet detected. so fully contained presumably".format(cell_RES.name, cell_i.name))
			cell_i.d_cell_contactinfo[cell_RES] = Contactinfo(cell_RES.volume)
			# print ("the centroid of RES-cell {0} is contained in the GT-cell {1} AND it was not yet detected. so fully contained presumably".format(cell_RES.name, cell_i.name))
			# cell_i.d_cell_contactinfo[cell_RES] = Contactinfo(0)


# Part4 : write output SEG3D + details to excel
write_excel_overlap_data()

