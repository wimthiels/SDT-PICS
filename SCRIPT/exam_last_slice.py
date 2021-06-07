"""
Created on 29 may 2020
gather the pixels from the bottom slice in order to use this information for volume correction (more pixels on the bottom slice = more volume missed)
The z-index of the last slice will be the z-ix that was marked as still having signal (max signal > no-signal-threshold). However if a manually annotated exterior mask is created, this
can overrule this
@author: wimth
"""
import numpy as np
import pandas as pd
import sys,re
from param_XML import Param_xml
from helper_functions import get_3D_sphere_coordinates,slice_tif
from tifffile import TiffFile,imsave

verbose=True

def read_parms():
	param_xml = Param_xml.get_param_xml(sys.argv, l_main_keys=['body', 'exam_last_slice'],verbose=True)  # param file must be passed as first argument

	input_file_preprocessing = param_xml.get_value('input_file_preprocessing', ['paths'])
	input_file_sphere_meshing = param_xml.get_value('input_file_sphere_meshing', ['paths'])
	input_file_exterior_mask = param_xml.get_value('input_file_exterior_mask', ['paths'])
	output_folder = param_xml.get_value('output_folder', ['paths'])
	nb_stack = param_xml.get_value('nb_stack',['process_flow'])

	l_dim = param_xml.get_value('RAW_IMG_DIM_TZYX',['body','RAM'],use_main_keys=False)
	VOXEL_DEPTH,_,PIXEL_WIDTH = param_xml.get_value('scaling_ZYX',l_path_keys=['body','RAM'],use_main_keys=False)
	XY_VS_Z = VOXEL_DEPTH/PIXEL_WIDTH 


	return input_file_preprocessing,input_file_sphere_meshing, input_file_exterior_mask, output_folder,nb_stack,l_dim,XY_VS_Z

#parms
input_file_preprocessing,input_file_sphere_meshing, input_file_exterior_mask, output_folder,nb_stack,l_dim,XY_VS_Z = read_parms()
output_folder.mkdir(parents=True,exist_ok=True)


#STEP 1a : get ix last slice (from no signal filter)
d_nbstack_ixzlast={} 
d_nbstack_ixzlast['nb_stack'] = nb_stack
# d_nbstack_ixzlast['ix_z_no_signal'] = [l_dim[1]]
# d_nbstack_ixzlast['ix_z_exterior_mask'] = [l_dim[1]]
if input_file_preprocessing: 
	if input_file_preprocessing.suffix=='.xlsx':
		df_file_prepro = pd.read_excel(input_file_preprocessing,sheet_name='frangi_vesselness_scores')
	else:
		df_file_prepro = pd.read_csv(input_file_preprocessing,header=0,index_col=False,keep_default_na=True) 
	if not nb_stack:
		nb_stack = df_file_prepro['nb_stack'].unique()
	l_ix_z = [0] * len(nb_stack)
	for ix,t_i in enumerate(nb_stack):
		print (t_i)
		flag_signal_detected = False
		for index, row in df_file_prepro[df_file_prepro['nb_stack']==t_i].iterrows():
			if row['no_signal_filter']==0:
				flag_signal_detected=True
				l_ix_z[ix]=int(row['ix_z'])
			else:
				if flag_signal_detected:
					break
				else:
					continue
	d_nbstack_ixzlast['ix_z_no_signal'] = l_ix_z


#STEP 1b : get ix last slice (from exterior mask)
if input_file_exterior_mask:
	with TiffFile(input_file_exterior_mask) as tif:
		a_exterior_mask = tif.asarray()
	if len(a_exterior_mask.shape)>3:
		a_exterior_mask = a_exterior_mask[nb_stack[0]-1,:]
	a_pixel_per_slice = np.sum(a_exterior_mask,axis=(1,2))
	for ix in range(a_exterior_mask.shape[0]-1,0,-1):
		if a_pixel_per_slice[ix] > 0:
			d_nbstack_ixzlast['ix_z_exterior_mask']= [ix]
			break

if input_file_preprocessing and input_file_exterior_mask:
	d_nbstack_ixzlast['ix_z'] = min(d_nbstack_ixzlast['ix_z_no_signal'],d_nbstack_ixzlast['ix_z_exterior_mask'])
else:
	if input_file_preprocessing:
		d_nbstack_ixzlast['ix_z'] = d_nbstack_ixzlast['ix_z_no_signal']
	else:
		d_nbstack_ixzlast['ix_z'] = d_nbstack_ixzlast['ix_z_exterior_mask']
print(d_nbstack_ixzlast)

#STEP 2 : get nb pixels in last slice
df_file_spheres = pd.read_csv(input_file_sphere_meshing,header=0,index_col=False,keep_default_na=True) 
	#init
d_cellname_parentID = {}  #this mimics the way parentID is assigned, just sequentially, 0-indexed
d_cellname_cellID = {} 
for ix, cell_label_i in enumerate(df_file_spheres['cell_label'].unique()):
	d_cellname_parentID[cell_label_i] = ix
	d_cellname_cellID[cell_label_i] = int(re.findall(r'\d+',cell_label_i)[0])
a_stack = np.zeros(tuple(l_dim[1:]),'int16')

	#place spheres in stack
for index, row in df_file_spheres.iterrows():
	t_pixels = get_3D_sphere_coordinates(a_coordZYX=[row['z'],row['y'],row['x']],
		radius=row['radius'],limit_shape=(l_dim[1],l_dim[2],l_dim[3]),xy_vs_z=XY_VS_Z)  #[Z,Y,X] 
	a_stack[t_pixels] = d_cellname_cellID[row['cell_label']] #parentID has 0, so equal to background 

imsave(output_folder/"stack_cellID.tif", a_stack.astype('uint16'))

	#stats on last stack

a_bincount = np.bincount(a_stack[d_nbstack_ixzlast['ix_z'][0]].flatten())
a_bincount[0] = 0
d_last_stack = {}
d_last_stack['parentid'] = np.subtract(np.nonzero(a_bincount)[0],2)
d_last_stack['nb_pix_last_z'] = a_bincount[np.nonzero(a_bincount)[0]]


#write output
with pd.ExcelWriter(output_folder/"last_slice_info.xlsx") as writer:
	df_last_stack = pd.DataFrame.from_dict(d_last_stack)
	df_last_stack.to_excel(writer, sheet_name='nb_pix',index=False)
	df_nbstack_ixzlast = pd.DataFrame.from_dict(d_nbstack_ixzlast)
	df_nbstack_ixzlast.to_excel(writer, sheet_name='ix_z',index=False)
	if verbose:print("excel output : {0}".format(str(output_folder/"last_slice_info.xlsx")))
