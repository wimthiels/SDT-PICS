'''
Created on 15 May 2019
- adaption mpacts_Celegans.py = original script wim
- python 2.7 +  unix (binaries don't work under python 3, no binaries for windows)

@author: wimth,barts
'''

#IMPORT 
from PIL import Image
import numpy as np
import vtk
import VTK_utils
import matplotlib.pyplot as plt
from scipy import ndimage,signal
import os,shutil,glob, math,sys
import xmltodict
from tifffile import TiffFile,imshow,imsave
from scipy.ndimage import zoom
from param_XML import Param_xml
print("os.popen('hostname').read()= ", os.popen('hostname').read())

#conf.enable_fp_exceptions()

#IMPORT MPACTS
import mpacts.core.simulation as sim
from mpacts.core.valueproperties import Variable, VariableFunction
from mpacts.core.units import unit_registry as u
import mpacts.core.command as cmd
import mpacts.geometrygenerators.trianglegeometries as trigeo
from mpacts.particles.particlecontainer import ParticleContainer
import mpacts.particles as prt
from mpacts.contact.models.collision.linearforce.linearforce_fl_matrix import FlatLayerAdhesiveTensionFprimIntMatrix
import mpacts.commands.geometry.orientation as orient
import mpacts.geometrygenerators.trianglegeometries as trigeo
import mpacts.geometrygenerators.pointgenerators as poi
from scipy.ndimage import gaussian_filter
import mpacts.geometrygenerators.transformations as transfo
from mpacts.contact.models.continuum.immersedboundary import GridToCell
from mpacts.commands.onarrays.setvalue import SetValueCommand
import DEMutilities.enumerate_properties as ep
import DEMutilities.relaxcontacts as relax
from mpacts.core.arrays import create_array
# from mpacts.core.arrays.ArrayManager import get_array
import mpacts.core.arrays as arr
from mpacts.particles.specialcases import DeformableCell
import mpacts.contact.detectors.multigrid as cdm
from random import random
import mpacts.contact.matrix.conjugategradient as cg
import mpacts.commands.monitors.progress as pp
import mpacts.tools.random_seed as RS
import mpacts.geometrygenerators.polyhedron as polyhedron
from mpacts.io.datasave import DataSaveCommand
import mpacts.io.datareader as read
from mpacts.io.savefiles import save_stl
from mpacts.io.filereaders import read_stl
import mpacts.core.configuration as conf
from mpacts.tools.setcontactdatastorage import SetContactDataStorage

import mpacts.contact.models.misc.contactcounter as counter
from mpacts.contact.detectors.reuse_existing import ReuseExistingContactDetector
from mpacts.commands.geometry.nodenormals import ComputeNodeNormalsCommand
import convexhull
#ReuseExistingContactDetector


class VolumeLossPrinter( pp.PrinterBase ):
	def __init__( self, pc):
		self.pc = pc

	def PrintProgress( self ):
		retstr = 'dvol: '
		retstr += '%3.2f'%((np.sum(self.pc['volume'][:])-np.sum(self.pc['volume_eq'][:]))/np.sum(self.pc['volume_eq'][:])*100)
		retstr += r'%'


		return retstr

def pre_filter( image, sigma,dz_dx_ratio ):
	'''
	Gaussian smooth filter with standard deviation sigma. Implementation is numpy so very fast
	'''

	if sigma==0:
		print('No smoothing applied')
		return image
	image = np.array(image, dtype=float)
	image /= np.max(image)
	image *= 256
	
	l_sigma = np.array([sigma/dz_dx_ratio, sigma, sigma])  #even smoothing in all directions, so correct for z

	return np.array(gaussian_filter(image, l_sigma), dtype=float) / 256. * np.sqrt(2*np.pi * sigma**2)**3

def resample( data, to_zoom, dx, dy, dz ):
	'''
	Resampling data with given 'zoom' ratios. Used to make it more evenly distributed across dimensions
	'''
	dx /= to_zoom[2]
	dy /= to_zoom[1]
	dz /= to_zoom[0]
	print(to_zoom)
	if to_zoom == (1.0,1.0,1.0):
		new_data = data
		print(" - No resampling needed")
	else:
		new_data = zoom( data, to_zoom, order=1, prefilter=False )  #prefilter off. the spline filter messes up the image!
	print  (" - resampling from ", np.shape(data), "to", np.shape(new_data), "for data ratio", new_data.size/float(data.size)*100, "%")
	print  ("    - after resampling dimensions are scaled as follows dx=",dx,"dy=",dy,"dz=",dz)

	return new_data, dx,dy,dz


def filter_Signal( data, dx,dy,dz ,param_xml,calibration_run_ic=False):
	'''
	Ibase : pixel_strength (float)  (array of float)
	xb : list of tuples with coordinates of selected pixels
	'''
	print (" - Filtering and compressing the signal...")

	data = np.array(data, dtype=float)

	#Select only a percentage of strongest pixels
	if pixel_selector_type=='preprocessed':
		perc_pix_sel = param_xml.get_value('perc_pix_sel',['parms_pix']) 
	else:
		perc_pix_sel = 100
	nb_pix_sel = param_xml.get_value('nb_pix_sel',['parms_pix']) 
	n_keep = nb_pix_sel if nb_pix_sel>0 else int(perc_pix_sel / 100. * data.size)
	n_keep = min(n_keep,np.count_nonzero(data))  # don't pick up zero pixels
	xb = np.unravel_index(data.flatten().argsort()[ -n_keep :], data.shape) #xb=tuple of 3 arrays of pixel coordinates

	if param_xml.get_value('balance_z_signal',['parms_pix']):
		mid_z = int(np.round(data.shape[0]/2))
		print(' - rebalancing z-signal with z reference {0}'.format(mid_z))
		a_max = np.max(data,axis=(1,2))
		a_max = a_max / a_max[mid_z]
		data = data / a_max[:,np.newaxis,np.newaxis]  #rebalancing signal per slice

	Ibase = data[xb]

	xb=tuple(l*r for l,r in zip((xb[2],xb[1],xb[0]),(dx,dy,dz))) # tuple of 3 arrays in micron (in x,y,z order)
	xb=np.transpose(xb)  #array (count,3)
	xb=list(map(tuple,xb)) #list of 3D-tuples in micron (x,y,z)
	
	#Compress
	if pixel_selector_type=='preprocessed':
		def compressor_signal(x): # compresses data between interval [1,2]
			return np.exp(-((x_max - x) / ((x + 1.e-6) - x_min))) + 1.0 
		Ibase = np.log10(Ibase) #for a more even spread between min and max
		x_max = np.max(Ibase)  
		x_min = np.min(Ibase)
		compressor_signal_vec = np.vectorize(compressor_signal)
		Ibase = compressor_signal_vec(Ibase)
	else:
		if pixel_selector_type=='calibration':
			Ibase[Ibase > 0] = 1
			print("pixel_selector_type={}.  All pixels are set to 1".format(pixel_selector_type))
		else:
			print('pixel_selector_type={}.  values from pixel selector are kept as is'.format(pixel_selector_type))

	print (" -  pixel data is now between: ", np.min(Ibase), "and", np.max(Ibase), " with sum =",np.sum(Ibase), " and average = ", np.average(Ibase))
	print (" -  retained ", Ibase.size/float( data.size )*100,  "% pixels for a Basegrid size of ", Ibase.size)
	
	return Ibase,xb

def extract_parentID_selection_from_VTP(input_file_vtp,verbose=False):
	'''
	temporary => copy paste from VTK_Utils
	'''
	d_parentID_VTP={}
	poly_in = read_vtp_file(input_file_vtp)
	a_parent_id = np.unique(dsa.WrapDataObject(poly_in).CellData['parentIndex'].__array__())
	for parent_id in a_parent_id:
		if verbose:print('--processing cell with parentID {0}'.format(parent_id))
		poly_out = extract_selection_vtp(poly_in, query="parentIndex == {0}".format(parent_id) ,field_type="CELL")
		d_parentID_VTP[parent_id] = poly_out

	return d_parentID_VTP
	
def add_cell_to_particle_container(a_vertices,a_triangles,pc_cell,cd_kd,mysim):
	p = polyhedron.Polyhedron()
	p.vertices = list(map(tuple, a_vertices))
	p.triangles = list(map(tuple, a_triangles)) #vils = connectivity , shape = (nb_triangles,3), contains point indices
	con = polyhedron.Connectivity(p)
	Cell = pc_cell.add_particle()

	#verts = transfo.translate(p.vertices, tuple(xcell[i]))
	Cell.nodes.add_and_set_particles( x     = p.vertices
									, v     = (0., 0., 0.)
									)

	Cell.triangles.add_and_set_particles(vertexIndices   = p.triangles
										, F              = (0,0,0)
										, contact_area   = 0.
										, layer_width    = cd_kd
										)
	Cell.add_connectivity( mysim("loop_cmds/contact_cmds/CD_Springs_cells"), con.edgeCorners )

	return Cell

def print_simulation_parameters():


	print('Parameters used : baseVisco={0},pt_cortex={1},CM_cell_cell.w1={2}'.format(baseVisco.get_value(),pt_cortex.get_value(),CM_cell_cell.get('w1')))
	print("Internal_pressure  = {}".format(p_int.get_value()))
	print("kv_cell = {}".format(kv_cell.get_value()))
	print("cortical tension = {}".format(pt_cortex.get_value()))
	print("cell adhesion  = {}".format(adh_cell.get_value()))
	print("image_force_ratio  = {}".format(image_force_ratio.get_value()))
	print("k_layer  = {}".format(k_layer.get_value()))

	return

def add_pixels_to_simulation():
	#cells('triangles').VTKWriterCmd(executeInterval=1, select_all=True, directory=str(output_folder))  #execute interval large (only 1 or 2 vtp snapshots of initial relaxation)
	#Adding the pixels : (This is much much faster as it does no dependency checking...)
	Grid.resize(len( xb ))
	Grid['x'].set_array(xb)
	Grid['r'].set_array([0 for _ in xb])
	Grid['Signal'].set_array(list(Ibase))
	Grid['m'].set_array( [0. for _ in xb] )
	Grid['h'].set_array( [rval for _ in xb])
	print( "Introduced signal grid. Total of {} pixels".format( Grid.size() ))
	Grid.VTKWriterCmd(executeOnce=True, select_all=True, directory=str(output_folder))  #basegrid output

	return

def scale_hull(scale_factor,xyz_centroid_hull=[0,0,0]):  
	""" scale using center of hull as the origin"""
	xbf = np.array(hull('controlPoints/xBF'))
	for xi in range(0,3):
		xbf[:,xi] -= xyz_centroid_hull[xi]
	xbf *=  scale_factor  
	for xi in range(0,3):
		xbf[:,xi] += xyz_centroid_hull[xi]
	hull("controlPoints")['xBF'].set_array(list(map(tuple, xbf)))
	return

def add_hull_to_simulation():
	ctrl, vil, xyz_centroid = convexhull.giftwrap(cloud=np.array(cells('nodes/x')) * 1e6,
							ignore_half=True,
							quantile_gut=0.4, #original 0.4
							divide_remove=5,  #5 originally
							seedlen=75,
							debug=False)
	ctrl *= 1e-6
	xyz_centroid *= 1e-6
	vil = [el[::-1] for el in vil]  # Normals should point inwards
	thehull.controlPoints.add_and_set_particles(xBF=list(map(tuple, ctrl)))
	thehull.hull.add_and_set_particles(vertexIndices=list(map(tuple, vil)), layer_width=10 * erange.get_value())  # orig 5

	if param_xml.get_value('scale_hull_factor',['process_flow']) != 1:
		scale_hull(param_xml.get_value('scale_hull_factor',['process_flow']),xyz_centroid)  #wth: loose fit for smoothing (in conjunction with larger repulsion distance)
		
	print("Introduced and positioned hull")
	hull.VTKWriterCmd(executeInterval=.05, select_all=True, directory=str(output_folder))

	return

def update_parms_simulation(phase=2):

	p_int_mult.set_value(1.e3*param_xml.get_value('p_int_{}'.format(phase),['parms_cell_model'])/param_xml.get_value('p_int',['parms_cell_model']))
	kv_cell_mult.set_value(1.e3*param_xml.get_value('kv_cell_{}'.format(phase),['parms_cell_model'])/param_xml.get_value('kv_cell',['parms_cell_model']))
	pt_cortex_mult.set_value(1.e-3*param_xml.get_value('pt_cortex_{}'.format(phase),['parms_cell_model'])/param_xml.get_value('pt_cortex',['parms_cell_model']))
	adh_cell_mult.set_value(1e-1*param_xml.get_value('adh_cell_{}'.format(phase),['parms_cell_model'])/param_xml.get_value('adh_cell',['parms_cell_model']))

	if phase in [1]:
		k_layer_mult.set_value(1.e12)
	else:
		k_layer_mult.set_value(1.e7)


	return
   
if __name__ == "__main__":
	#files&folders
	import time;tic = time.time() 

	param_xml = Param_xml.get_param_xml(sys.argv,l_main_keys = ['body','mpacts_PiCS'],verbose=True)
	docker_run_ic =  param_xml.get_value('docker_run_ic',['paths'])
	input_file_pixel_selector= param_xml.get_value('input_file_pixel_selector',['paths'],docker_run=docker_run_ic)
	pixel_selector_type = param_xml.get_value('pixel_selector_type',['process_flow'],docker_run=docker_run_ic)
	output_folder = param_xml.get_value('output_folder',['paths'],docker_run=docker_run_ic)
	print (str(input_file_pixel_selector), " is used to extract pixels")

	with TiffFile(input_file_pixel_selector) as tif:
		a_pixel_selector=tif.asarray()
	if len(a_pixel_selector.shape)==4:
		a_pixel_selector=a_pixel_selector.reshape(a_pixel_selector.shape[1:])

	if output_folder.exists() and '012_mpacts_PiCS' in output_folder.parts: #safety check to prevent disastrous delete
		shutil.rmtree(output_folder, ignore_errors=False, onerror=None)
	output_folder.mkdir(parents=True,exist_ok=True)  #
	calibration_run_ic = param_xml.get_value('calibration_run_ic',['process_flow'])
	# for f in list(output_folder.glob('*')): # clear outputfolder
		# os.remove(str(f))
		
		
	SetContactDataStorage("Vector")
	#FILTER SIGNAL
	# dx,dy,dz = param_xml.get_value('dxyz',lookup=input_file_pixel_selector.name)
	dz,dy,dx = param_xml.get_value('scaling_ZYX',l_path_keys=['body','RAM'],use_main_keys=False)
	dx*=1e-6;dy*=1e-6;dz*=1e-6;
	print (" - pixels converted to micron using dx=",dx,"dy=",dy,"dz=",dz)

	a_pixel_selector = pre_filter(a_pixel_selector, sigma=param_xml.get_value('sigma_smoothing',['parms_pix']), dz_dx_ratio=dz/dx)                  #smoothing
	a_pixel_selector,dx,dy,dz = resample(a_pixel_selector, tuple(param_xml.get_value('zoom_factors',['parms_pix'])), dx, dy, dz) #resampling 
	Ibase,xb = filter_Signal(a_pixel_selector,dx,dy,dz,param_xml,calibration_run_ic)                                                             #filtering + scaling

	(nz,ny,nx)=a_pixel_selector.shape
	print (" - before correction=%d, after=%d" %((nz*(nx)*(ny)),Ibase.size))

	##SETTING UP SIMULATION-----------------------------------------------------------------------------------------

	dtg = 5e-3 #timestep
	mysim = sim.Simulation("Seeding", timestep=dtg)
	params = mysim('params') #Creating lib with settings

	image_force_ratio = param_xml.get_value('image_force_ratio',['parms_pix'])

	#calibrating the image force to the number of pixels
	image_force_ratio_mult =  Variable("image_force_ratio_mult", params, value=1) 
	if pixel_selector_type=='raw':
		image_force_ratio = Variable("image_force_ratio", params, value = image_force_ratio)#How much does image contribute compared to mechanics
	else:
		image_force_ratio_scaled = image_force_ratio/np.average(Ibase)
		print (" - The image_force_ratio from parameterfile = {0} and downscaled by the average pixel intensity to {1}".format(image_force_ratio,image_force_ratio_scaled ))
		print (" - total force of the pixels now equals ", np.sum(Ibase* image_force_ratio_scaled))
		image_force_ratio = VariableFunction("image_force_ratio", params, function="{}*$image_force_ratio_mult$".format(image_force_ratio_scaled))
		#image_force_ratio = Variable("image_force_ratio", params, value = image_force_ratio_scaled )#How much does image contribute compared to mechanics

	param_xml.l_main_keys.append('parms_cell_model')
	#Reference value for mechanical tension. This determines base mechanical state. Best to keep unchanged:
	
	pt_cortex_mult = Variable("pt_cortex_mult", params, value=1*u("nN/um"))  
	pt_cortex  = VariableFunction("pt_cortexn", params, function="{}*$pt_cortex_mult$".format(param_xml.get_value('pt_cortex'))) 
	#pt_cortex   = Variable( "pt_cortexn"      , params, value=  param_xml.get_value('pt_cortex')*u('nN/um'))
	function    = '$pt_cortexn$'
	pt_cortex_f   = VariableFunction("pt_cortex_f"  , params, function=function )

	kv_cell_mult = Variable("kv_cell_mult", params, value=1*u("kPa"))  
	kv_cell  = VariableFunction("kv_celln", params, function="{}*$kv_cell_mult$".format(param_xml.get_value('kv_cell'))) 
	#kv_cell     = Variable( "kv_celln"        , params, value=  param_xml.get_value('kv_cell')*u("kPa"))       #Volume conservation.

	p_int_mult = Variable("p_int_mult", params, value=1*u("kPa"))  
	p_int = VariableFunction("internal_pressure", params, function="{}*$p_int_mult$".format(param_xml.get_value('p_int'))) 
	# p_int     = Variable( "internal_pressure"        , params, value=  param_xml.get_value('internal_pressure')*u("kPa"))       #Balloon force

	h_eff       = Variable( "h_eff"           , params, value= param_xml.get_value('h_eff')*u('nm'))           #Adhesive range (smaller = stiffer, larger = soft and smooth)
	h_eff_pair  = Variable( "h_eff_pair"           , params, value= (param_xml.get_value('h_eff') / 1.3 )*u('nm'))  #Adhesive range (smaller = stiffer, larger = soft and smooth) previously 10
	#In the quadratic form, the maximal force will be at 1/sqrt(2)*h, so roughly 0.7*h
	h_smooth    = Variable( "smoothing_length", params, value=param_xml.get_value('h_smooth')*u('um'))         #Smoothing length of image pixels (big -> slow calculation)
	radius      = Variable( "radius"          , params, value=param_xml.get_value('radius')*u('um'))           #Average cell radius
	max_r       = Variable( "max_inv_curv"    , params, value=param_xml.get_value('max_r')*u('um'))
	opt_limit   = Variable( "cg_opt_limit"    , params, value=param_xml.get_value('opt_limit') )               #Tolerance for conjugate gradient solver (large=fast and imprecise)

	FrictionN   = Variable( "FrictionN"       , params, value= param_xml.get_value('FrictionN')*u('kPa*s/um')) #Normal contact friction
	FrictionT   = Variable( "FrictionT"       , params, value= param_xml.get_value('FrictionT')*u('kPa*s/um')) #Tangential contact friction
	viscosity   = Variable( "viscosity"       , params, value= param_xml.get_value('viscosity')*u('kPa*s'))    #Diagonal (liquid) viscosity: high = stable
	baseVisco   = Variable( "baseVisco"       , params, value= param_xml.get_value('baseVisco')*u('kPa*s'))    #In-plane cortex viscosity
	tc          = Variable( "tc"              , params, value= param_xml.get_value('tc')*u('um'))              #Cortex thickness (typically 250 nm)
	Vpoiss      = Variable( "Vpoiss"          , params, value= param_xml.get_value('Vpoiss'))                  #As inconsequential as parameters can be...

	function    = str(param_xml.get_value('gamma_liq')) + '*$viscosity$/$radius$'
	gamma_liq   = VariableFunction("gamma_liquid"  , params, function=function )

	function    = '$image_force_ratio$*$pt_cortexn$'
	k_signal    = VariableFunction('k_signal'      , params, function=function )            #important NOTE - k_signal should in principle be scaled based on the number of grid pixels.

	function    = str(param_xml.get_value('cd_kd'))     + '*$h_eff$'   
	cd_kd       = VariableFunction('cd_kd'         , params, function=function )

	adh_cell_mult = Variable("adh_cell_mult", params, value=1000*u("nN/um"))  #0.001
	function    = "{}*$adh_cell_mult$".format(param_xml.get_value('adh_cell'))  + '*$pt_cortexn$'  #multiply so adh_cell expresses a ratio of adhesion and surface tension. 
	adh_cell    = VariableFunction("adh_celln"     , params, function=function )
	
	adh_hull    = Variable("adh_hull", params, value=0.00002)  # Adhesion energy cell hull

	function    = str(param_xml.get_value('ka'))        + '*$pt_cortexn$'
	ka          = VariableFunction("ka"            , params, function=function )           #prevents triangles from collapsing or blowing up too much (no effect for small deformations, kicks in when triangles become too extreme)

	function    = '$baseVisco$'
	baseVisco_f   = VariableFunction("baseVisco_f"  , params, function=function )


	k_layer_mult = Variable("k_layer_mult", params, value=1e12)  
	k_layer = VariableFunction("k_layer", params, function="{}*$k_layer_mult$".format(1)) 


	del param_xml.l_main_keys[-1]

	#------------------------------------------------------#
	#Creating the base for the Signal

	Grid = ParticleContainer("BaseGrid", prt.Sphere, mysim)
	Grid.create_array('Scalar', 'Signal')
	Grid.create_array("Scalar", "h")

	#------------------------------------------------------#
	#Creating the base for the cells

	cells   = prt.ParticleContainer("cells", DeformableCell , mysim)
	create_array("Scalar",'edge_length',cells("triangles"))
	create_array("Scalar",'layer_width',cells("triangles"))
	create_array("Index",'reject',cells("triangles"))
	create_array("unsigned",'number_of_effective_contacts',cells("triangles"))
	create_array("int",'contact_index',cells("triangles"))

	cells.DeformableCellGeometryCmds(maximal_radius = max_r, enforce_positive = False)


	cells.DeformableCellInternalMechanicsCmds(kv                = kv_cell
											 , internal_pressure = p_int
											 , nu                = Vpoiss
											 , viscosity         = baseVisco_f
											 , surface_tension   = pt_cortex_f
											 , thickness         = tc
											 , indirect_surface_tension_model = False
											 , ka=0.
											 , ka_compress = ka
											 )
	#------------------------------------------------------#
	#Adding contact mechanics and interpolation stuff
	CM_cell_cell    = FlatLayerAdhesiveTensionFprimIntMatrix( pc1               = cells('triangles')
															, pc2               = cells('triangles')
															, w1                = adh_cell
															, w2                = adh_cell
															, gamma_normal      = FrictionN
															, gamma_tangential  = FrictionT
															, Fprim1            = cells('triangles')['F']
															, Fprim2            = cells('triangles')['F']
															, contact_range     = h_eff
															, implicitness      = 1.0
															, reject_angle      = np.pi/2.
															, k_layer=k_layer #1e3   #1e12 ##speedup repulsion  Layer stiffness. 
															)

	h0 = h_smooth.get_value()
	h = (h0,h0,h0)
	rval = min( h )
	CM_grid_cell= GridToCell(pc1=Grid, pc2=cells('nodes'), h=h,scale_factor=k_signal)

	#------------------------------------------------------#
	#Adding contact detectors

	CD_cell_cell = cdm.MultiGridContactDetector("Cell-Cell"
												, mysim
												, cmodel        = CM_cell_cell
												, update_every  = 5
												, keep_distance = cd_kd
												)

	CD_cell_grid = cdm.MultiGridContactDetector("Cell-grid"
											   , mysim
											   , cmodel        = CM_grid_cell
											   , update_every  = 25
											   , keep_distance = 1.799*h0 #Beyond this value, the function is less than 1/10th of maximal value
											   )

	#extra contactmodel for cell pair contact
	#-----------------------------------------

	#werkt niet met binaries : *** RuntimeError: Unknown type id: 'Deformable_RoundedTriangle Deformable_RoundedTriangle'
	CM_counter=counter.ContactIndexKeeper(pc1=cells('triangles'),pc2=cells('triangles'),contact_range=h_eff_pair,reject_angle = np.pi/2)

	CD_counter = ReuseExistingContactDetector("cell-Counter"
												, mysim
												, cmodel        = CM_counter
												, cd            = CD_cell_cell
												)


	#------------------------------------------------------#
	#Adding viscous friction
	cells("nodes").FrictionOnAreaCmd( area = cells('nodes')['area']
									, gamma_normal = gamma_liq
									, gamma_tangential = gamma_liq
									)

	#------------------------------------------------------#
	#Adding a GC solver for the system
	CG=cg.DefaultConjugateGradientSolver(sim=mysim, tolerance=opt_limit, reset_x=False)

	#------------------------------------------------------#
	#Stability command for a viscous model
	cells.MeshUniformizationCmds( relaxation_time = dtg*20, weight_with_area = False)

	#------------------------------------------------------#
	#Actually adding cells
	input_file_vtp = param_xml.get_value('input_file_vtp',['paths'],docker_run=docker_run_ic)
	if input_file_vtp:
		print("using vtp file as input for cells : {0}".format(input_file_vtp))
		from VTK_utils import get_point_array,read_vtp_file,extract_selection_vtp,get_point_array
		import vtk.numpy_interface.dataset_adapter as dsa
		d_parentID_VTP = extract_parentID_selection_from_VTP(str(input_file_vtp),verbose=False)
		for parentID,poly_cell in d_parentID_VTP.items(): 
			print(" - loading cell '{}' from vtp".format(parentID))
			wdo_cell  = dsa.WrapDataObject(poly_cell)
			add_cell_to_particle_container(wdo_cell.Points,wdo_cell.GetAttributes(vtk.vtkDataObject.CELL)['vertexIndices'],cells,cd_kd.get_value(),mysim)
	else:
		stlfiles = sorted(param_xml.get_value('input_folder_meshed_cells',['paths'],docker_run=docker_run_ic).glob('cell*.stl')) #Cells are saved as a collection of names stl files, in a separate folder.  
		SELECT_IDX_CELL = param_xml.get_value('SELECT_IDX_CELL',['process_flow'],docker_run=docker_run_ic)
		for i in range(len(stlfiles)):
			if (len(SELECT_IDX_CELL)==0) or (i in SELECT_IDX_CELL):
				print(" - loading cell '{}'".format(stlfiles[i]))
				ctrl,vils = read_stl( str(stlfiles[i]) )
				xs = np.array(ctrl)*1e-6 #BART: Convert from micron to SI units.
				add_cell_to_particle_container(xs,vils,cells,cd_kd.get_value(),mysim)

	print ("Done adding cells. Total of {} nodes and {} triangles".format( cells('nodes').size(),cells('triangles').size() ))

	if param_xml.get_value('fit_hull',['process_flow']):
		# Hull ----------------------------------------------------------------------------
		erange = Variable("contact_range", params, value=25 * u('nm'))  #wth previously 250, but expanded in conjuction with scaling of hull 
		hull = prt.ParticleContainer("hull", prt.RigidBody.compose((prt.RigidTriangle, "hull", "controlPoints")), mysim)
		tension = create_array("Scalar", "tension", hull("hull"))
		Fprim = create_array("Vector", "F", hull("hull"))
		create_array("Scalar", "layer_width", hull("hull"))
		ca = create_array("Scalar", "contact_area", hull("hull"))
		create_array("Scalar", "r", hull)

		ComputeNodeNormalsCommand("ComputeHullNodeNormals", mysim, pc=hull("hull"), nodes=hull("controlPoints"))
		hull("hull").TriangleNormalsAndAreaCmd(x=hull("controlPoints")['x'])
		SetValueCommand("ResetHullTension", mysim, value=0., array=tension)
		SetValueCommand("ResetHullContactArea", mysim, value=0., array=ca)
		SetValueCommand("ResetHullFPrim", mysim, value=(0., 0., 0.), array=Fprim)

		hull.SetContactMatrixDiagonalCmd(visc=1e8)  # hull may not migrate ;-)

		CM_cell_hull = FlatLayerAdhesiveTensionFprimIntMatrix(pc1=cells('triangles')
															  , pc2=hull('hull')
															  , w1=adh_hull
															  , w2=0
															  , gamma_normal=FrictionN.get_value() / 20
															  , gamma_tangential=FrictionT.get_value() / 20
															  , Fprim1=cells('triangles')['F']
															  , Fprim2=hull('hull')['F']
															  , contact_range=erange
															  , implicitness=1.0
															  , k_layer=1e12  #with : used to be 100e9  # flat layer adhesion (repulsion)
															  )  

		CD_cell_hull = cdm.MultiGridContactDetector("Cell-hull"
													, mysim
													, cmodel=CM_cell_hull
													, update_every=25  # orig 5
													, keep_distance=cd_kd.get_value())

		thehull = hull.add_particle()
		thehull.q = (1.0, 0, 0, 0)  # quaternions
		thehull.v = (0, 0, 0)  # did nothing
		thehull.r = radius.get_value() - 0.25 * erange.get_value()  # I have no idea why this is done
	#------------------------------------------------------#
	#some printers
	printerl = [pp.DefaultProgressPrinter(mysim) #simulation time, timestep, wall-time, and simulation-/wall-time fraction.
				, VolumeLossPrinter( cells)
				, pp.PropertyPrinter(CG('steps'), "CG ")
				]

	printer = pp.PrinterChain(printerl)
	pp.ProgressIndicator("PrintProgress"
						, mysim
						, printer         = printer
						, print_interval  = 3
						)

	#------------------------------------------------------#
	#RUNNING SIMULATION
	#------------------------------------------------------
	#set output frequency
	
	DataSaveCommand( "SaveDataCommand", mysim, executeInterval = 1.0, directory=str(output_folder) )  # mpacts data files...  execute interval smaller : more snapshots)
	cells('triangles').VTKWriterCmd(executeInterval=.05, select_all=True, directory=str(output_folder))  #vtp files


	#PHASEI------------
	if param_xml.get_value('fit_hull',['process_flow']):
		add_hull_to_simulation()
	print (">>>PhaseI:cell relaxation (without pixels); simulation time : {0}s".format(param_xml.get_value('simulation_time_1',['parms_sim'])))
	print_simulation_parameters()
	mysim.run_until(mysim.get_time() + param_xml.get_value('simulation_time_1',['parms_sim']))

	cells('triangles').VTKWriterCmd(executeInterval=.01, select_all=True, directory=str(output_folder))  #vtp files (default every 0.05s sim time)
	#PHASEII--------------
	add_pixels_to_simulation() 
	update_parms_simulation(phase=2)

	print (">>>PhaseII:cell expansion (with pixels); simulation time : {0}s".format(param_xml.get_value('simulation_time_2',['parms_sim'])))
	print_simulation_parameters()
	mysim.run_until(mysim.get_time() + param_xml.get_value('simulation_time_2',['parms_sim']))   


	#PHASEIII-------------------
	#add_pixels_to_simulation()
	update_parms_simulation(phase=3)
	#image_force_ratio_mult.set_value(0)
	if param_xml.get_value('fit_hull',['process_flow']):
		#mysim.remove_child(hull)  #problem later on : Multigrid still active
		scale_hull(0) # effectively shortcutting the hull


	print (">>>PhaseIII:Relaxation with pixels; simulation time : {0}s".format(param_xml.get_value('simulation_time_3',['parms_sim'])))
	print_simulation_parameters()
	mysim.run_until(mysim.get_time() + param_xml.get_value('simulation_time_3',['parms_sim']))

	print('>>>end of simulation')

	#save centroids
	np.savetxt(str(output_folder / "cell_centroids.csv"), np.array(cells('x')), delimiter=",")
	toc = time.time()
	print('runtime_mpacts_pics = ', toc-tic)