#!/usr/bin/env python-mpacts
# param file must be passed as first argument

import mpacts.commands.misc.skinsurfacemesh as ssm
import mpacts.core.simulation as sim
import mpacts.core.arrays as ar
import mpacts.particles as prt
# from PyACVD import Clustering
from pyacvd import clustering
import mpacts.io.savefiles as write
import mpacts.io.filereaders as read
import mpacts.io.datasave as save
from mpacts.io.datasave import DataSaveCommand
from mpacts.geometrygenerators.meshobjects import TriangleMeshObject
from pymeshfix import _meshfix
import vtk, os, sys, glob
import numpy as np
from param_XML import Param_xml

import pyvista

os.system('hostname')




# files&folders
def read_param():
	param_xml = Param_xml.get_param_xml(sys.argv, l_main_keys=['body', 'sphere_meshing'],
										verbose=True)  # param file must be passed as first argument
	docker_run_ic = param_xml.get_value('docker_run_ic', ['paths'])
	output_folder = param_xml.get_value('output_folder', ['paths'], docker_run=docker_run_ic)
	if not output-folder:
		output_folder = param_xml.get_value('output_folder_meshes', ['paths'], docker_run=docker_run_ic) #backward compatibility 
	input_file_spheres = param_xml.get_value('input_file_spheres', ['paths'], docker_run=docker_run_ic)
	shrink_factor = param_xml.get_value('shrink_factor', ['parms'])
	threshold_radius_micron = param_xml.get_value('threshold_radius_micron', ['parms'])

	return output_folder, input_file_spheres, shrink_factor, threshold_radius_micron


# def init_output_folder():
#     try:
#         output_folder.mkdir(parents=True)  #exist_OK not available in pyhon2
#     except:
#         pass
#         # for f in list(output_folder.glob('*')):
#             # os.remove(str(f))
#     return

def remesh_file(inputname, outputname, N):
	
	_meshfix.clean_from_file(inputname, outputname)

	stlReader = vtk.vtkSTLReader()
	stlReader.SetFileName(outputname)
	stlReader.Update()
	mesh = stlReader.GetOutput()
	mesh = pyvista.PolyData(mesh)

	cobj = clustering.Clustering(mesh)
	cobj.subdivide(3)  # lijn toegevoegd tov python3 versie bart, want remeshing lukte niet cfr pyacvd documentatie
	cobj.cluster(N)
	remesh = cobj.create_mesh()

	w = vtk.vtkSTLWriter()
	w.SetFileName(outputname)
	w.SetInputData(remesh)
	w.Update()

	_meshfix.clean_from_file(outputname, outputname)

	ctrl, vil = read.read_stl(outputname)
	mesh = TriangleMeshObject()
	for xi in ctrl:
		mesh.add_vertex(tuple(xi))
	for ti in vil:
		mesh.add_triangle(tuple([int(el) for el in ti]))
	if mesh.is_cavity():
		vil = [el[::-1] for el in vil]
	write.save_stl(ctrl, vil, outputname)

	return


# MAIN######################################################################################################
output_folder, input_file_spheres, shrink_factor, threshold_radius_micron = read_param()
output_folder.mkdir(parents=True, exist_ok=True)
print("meshing started from input file: {0}".format(str(input_file_spheres)),flush=True)

# read inputexcel with cell / sphere data
data = np.recfromcsv(input_file_spheres)
data.sort(order='radius_micron')# ascending : smallest spheres first (skinsurface mesher seems to be more stable this way)
#print(data)
names = list(set(data['cell_label']))
# names = [i.astype('str') for i in set(data['cell_label'])] # datatype cell label = numpy.byte
print("found names {}".format(sorted(names)))

globsim = sim.Simulation("collection", timestep=1)
result = prt.ParticleContainer("result", prt.DeformableBody.compose((prt.Node, "nodes")
																	, (prt.DeformableRoundedTriangle, "triangles",
																	   "nodes")), globsim)

for label in sorted(names):
	print(" - start creating mesh of cell: '{}'".format(label.astype('str')), end="",flush=True)
	# pdb.set_trace()
	# sys.stdout = open(os.devnull, 'w')

	# PART1: skin surface meshing------------------------------
	mysim = sim.Simulation(label, timestep=1)

	spheres = prt.ParticleContainer("spheres", prt.Sphere0, mysim)

	cells = prt.ParticleContainer("cells", prt.DeformableBody.compose((prt.Node, "nodes")
																	  , (prt.DeformableRoundedTriangle, "triangles",
																		 "nodes")), mysim)

	glue_cmd = ssm.ComputeSphereSkinSurfaceMeshCommand("GlueSpheresCommand", mysim, pc=spheres
													   , nodes=cells('nodes')
													   , triangles=cells("triangles")
													   , shrink_factor=shrink_factor
													   , n_subd=0)

	mask = (data['cell_label'] == label) * (data['radius_micron'] > threshold_radius_micron)
	x = np.c_[data['x_micron'][mask], data['y_micron'][mask], data['z_micron'][mask]]
	r = data['radius_micron'][mask]

	# if r[0] < 1.e-4:
	# r *=1.e6;x *=1.e6  #if data is inputted in micron, pyACVD does not work properly (very coarse grained meshes).  so scaling required

	x = list(map(tuple, x))
	r = r.tolist()

	spheres.add_and_set_particles(x=x, r=r)
	try:
		print('starting skin mesh', end="",flush=True)
		mysim.run_n(1)  #the skin mesh can sometimes fail (Fatal Python error), without errorexception, script closes, no error warning (reason ?)
	except:
		print(' skin surface meshing of {0} failed, moving on to next cell'.format(label),flush=True)
		continue
	

	ctrl = cells('nodes')['x'].tolist()
	vil = cells('triangles')['vertexIndices'].tolist()

	# sys.stdout = sys.__stdout__

	filename = str(output_folder / 'skin_{}.stl'.format(label.astype('str')))
	write.save_stl(ctrl, vil, filename)
	print("...skin_surface meshing done..",end="")

	# PART2 : remeshing-----------------------------------
	print("...remeshing")
	stl_name = str(output_folder / 'cell_{}.stl'.format(label.astype('str')))
	try:
		remesh_file(filename, stl_name, 1200)
	except:
		print(' remeshing of {0} failed, moving on to next cell'.format(label))
		continue

	# sys.stdout = sys.__stdout__

	ctrl, vil = read.read_stl(stl_name)
	acell = result.add_particle()
	acell.nodes.add_and_set_particles(x=list(map(tuple, ctrl)))
	acell.triangles.add_and_set_particles(vertexIndices=list(map(tuple, vil)))

# sys.stdout = sys.__stdout__
globsim.run_n(1)
