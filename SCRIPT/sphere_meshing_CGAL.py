"""
Created on 08 January 2021
post processing of CGAL skin meshing : convert off-files into stl, remesh and repair
@author: wimth
"""

from pathlib import Path
import os,glob,sys
import trimesh
from pymeshfix import _meshfix
from pyacvd import clustering
import pyvista
from VTK_utils import write_stl_file
import pyvista
from param_XML import Param_xml
verbose=True

def read_parms():
	param_xml = Param_xml.get_param_xml(sys.argv, l_main_keys=['body'],verbose=True) 
	param_xml.add_prms_to_dict(l_keys=['sphere_meshing'])
	param_xml.repr_prm()
	return param_xml.prm


prm = read_parms()

d_path_mesh = {}
for path in prm['output_folder_meshes'].glob('*.off'):  #this script works on the output folder of the CGAL meshing

	if verbose:print(' post processing mesh : {}'.format(path.name),flush=True,end="")
	f_out = prm['output_folder_meshes'] / (path.stem  + ".stl")
	
	#convert off to stl
	d_path_mesh[path] = trimesh.load_mesh(str(path),process=False)
	if verbose:print(' -> converted to stl'.format(path.name),flush=True,end="")

	## remesh    
	mesh = pyvista.wrap(d_path_mesh[path])  #needs pyvista 0.27.x
	clus = clustering.Clustering(mesh)  #wrap gives VTK polydata
	clus.subdivide(3)
	clus.cluster(1200)
	remesh = clus.create_mesh()
	write_stl_file(remesh,f_out)
	if verbose:print(' -> remeshed'.format(path.name),flush=True,end="")
	
	## repair
	try:
		_meshfix.clean_from_file(str(f_out), str(f_out))   #repairing, just in case
		if verbose:print(' -> repaired'.format(path.name),flush=True,end="")
	except:
		print('something went wrong with _meshfix on {0}. continue'.format(f_name))
	if verbose:print("")
