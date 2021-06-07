'''
Created on 4 Apr 2019

@author: wimth
'''
import numpy as np
from shutil import copyfile
from tifffile import TiffFile,imsave
import sys ,re,os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from skimage import img_as_ubyte
from skimage.color import grey2rgb
from skimage.draw.draw3d import ellipsoid
from skimage.draw import line
from math import ceil,sqrt
from param_XML import Param_xml
from skimage.morphology import remove_small_objects
from fileHandler import FileHandler
import tifffile
# import spam.helpers
#from pyevtk.hl import gridToVTK
import pyvista as pv
from VTK_utils import read_vtp_file, get_point_array, get_data_array
import vtk
from scipy.ndimage import zoom
from scipy import ndimage


# from pv.DataSetFilters import slice_along_axis

def display_top(snapshot, key_type='lineno', limit=3):
	import tracemalloc
	import linecache
	
	snapshot = snapshot.filter_traces((
		tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
		tracemalloc.Filter(False, "<unknown>"),
	))
	top_stats = snapshot.statistics(key_type)

	print("Top %s lines" % limit)
	for index, stat in enumerate(top_stats[:limit], 1):
		frame = stat.traceback[0]
		# replace "/path/to/module/file.py" with "module/file.py"
		filename = os.sep.join(frame.filename.split(os.sep)[-2:])
		print("#%s: %s:%s: %.1f KiB"
			  % (index, filename, frame.lineno, stat.size / 1024))
		line = linecache.getline(frame.filename, frame.lineno).strip()
		if line:
			print('    %s' % line)

	other = top_stats[limit:]
	if other:
		size = sum(stat.size for stat in other)
		print("%s other: %.1f KiB" % (len(other), size / 1024))
	total = sum(stat.size for stat in top_stats)
	print("Total allocated size: %.1f KiB" % (total / 1024))

def vtp_to_tif(vtp_file, t_tif_shape, path_output, l_zyx = [1e-6,0.129e-6,0.129e-6]):
	"""
	very inefficient : requires for-looping, drawing a line for every polyline
	"""
	l_xyz = l_zyx.copy()
	l_xyz.reverse()
	a_tif = np.zeros(t_tif_shape)
	pv_poly = pv.read(vtp_file)
	for ix_z in range(t_tif_shape[0]):
		pv_slice = pv_poly.slice(normal='z', origin=[ 0, 0, l_zyx[0] * (ix_z + 1)], generate_triangles=False, contour=False)
		for ix_cell in range(pv_slice.n_cells):
			a_points_cell= pv_slice.extract_cells(ix_cell).points
			a_points_cell_xyz = np.floor(a_points_cell / l_xyz).astype(int)
			l_points = line(a_points_cell_xyz[0][0], a_points_cell_xyz[0][1], a_points_cell_xyz[1][0], a_points_cell_xyz[1][1])  #line(r0, c0, r1, c1) 
			l_points_z_added = ( [ix_z] * len(l_points[0] )  , l_points[1], l_points[0] )  #back to zyx format
			a_tif[l_points_z_added] = 1
			
	imsave(str(path_output / "slices.tif" ),a_tif.astype('uint16'))

	return


def tif_to_paraview_tif(p_input_tif,output_folder,nb_stack=None):
	"""
	convert a tif into a tif that can be correctly vizualised in paraview (y-axis needs to be flipped = rotation around center line)
	in paraview one still needs to scale relatively (eg 0.129,0.129,1.1) to vizualize it correctly (per convention : do not use micron scaling)
	
	"""
	with TiffFile(str(p_input_tif)) as tif:a_raw = tif.asarray()

	if nb_stack:
		a_raw = a_raw[nb_stack[0]-1,:]
	a_paraview = np.zeros_like(a_raw)
	_,y_max,_ = a_paraview.shape 
	for i_y in range(0,y_max):
		a_paraview[:,i_y,:] = a_raw[:,-1*(i_y- y_max+1),:]
		  
	imsave(str(output_folder / (p_input_tif.stem + "_paraview.tif")) , a_paraview)

	return 

def vtp_to_vtm(vtp_file,t_tif_shape,path_output):
	"""
	vtp to vtm = vtk multiblock (will pick even slices across z, the exact z-locations are not under your control)
	"""
	pv_poly = pv.read(vtp_file)
	#ex
	slices = pv_poly.slice_along_axis(n=t_tif_shape[0], axis='z', tolerance=None, generate_triangles=False, contour=False, bounds=None, center=None)
	slices.save(str(path_output / "slices.vtm"))
	
	return
	
def numpy3D_to_vtk(data, path_file):
	"""
	https://www.kaggle.com/ironbar/saving-3d-files-to-vtr-for-later-visualization/
	
	save the 3d data to a .vtk file. 
	
	ea : werkt niet voor mij ? saved wel een vtr, maar paraview zegt dat hij deze niet kan lezen
	
	Parameters
	------------
	data : 3d np.array
		3d matrix that we want to visualize
	filepath : str
		where to save the vtk model, do not include vtk extension, it does automatically
	"""
	x = np.arange(data.shape[0]+1)
	y = np.arange(data.shape[1]+1)
	z = np.arange(data.shape[2]+1)
	
	gridToVTK(str(path_file.parent), x, y, z, cellData={path_file.stem:data.copy()})

# def tif_to_vtk(tif_path, vtk_path):
	# '''
	# cfr : https://ttk.gricad-pages.univ-grenoble-alpes.fr/spam/vtk.html
	# instead of importing a tiff into paraview (which flips the y-axis), turn it iot a vtk
	# => deprecate : does not work !  vtk file is empty
	# '''
	
	# greys = tifffile.imread(str(tif_path))
	# spam.helpers.writeStructuredVTK(cellData={'grey':greys}, fileName=str(vtk_path))
	
	# return 


def get_3D_sphere_coordinates(a_coordZYX=[0,0,0],radius = 10, limit_shape=(9999,9999,9999),xy_vs_z=1,verbose=False):
	'''
	get coordinates of a sphere on a certain locations, and within the perimeter set out by limits
	outputs an 3D array with the coordinates [Z,Y,X]
	:param a_coordZYX: where the center of the sphere must be placed
	:param radius: the radius in l_pixels (in the xy-plane, the z radius will be scale appropriately)
	:param limit_shape: the dimensions of the frame (Z,Y,X)
	:param verbose:
	
	:output [Z,Y,X]  representing the coordinates of the filled sphere 
	'''
	
	#get a filled sphere (in real size) as a boolean matrix (True=part of sphere)
	na_sphere = ellipsoid(max(1,round(radius/xy_vs_z)),radius,radius,levelset=False)
	if verbose:print('na_sphere', na_sphere,'shape', na_sphere.shape)
	
	# place this sphere in a X Y Z space defined by limits
	a_center_of_sphere = [int((i-1)/2) for i in na_sphere.shape]
	if verbose:print('a_center_of_sphere',a_center_of_sphere)
	a_delta_move = np.subtract(a_coordZYX,a_center_of_sphere)
	if verbose:print('a_delta_move',a_delta_move)
	Z,Y,X, = np.where(na_sphere)
	Z=np.add(Z,a_delta_move[0]);Y=np.add(Y,a_delta_move[1]);X=np.add(X,a_delta_move[2])
	
	z_lim, y_lim, x_lim = limit_shape
	sphere_coord = [coord for coord in zip(Z,Y,X) if (coord[0]>=0 and coord[1]>=0 and coord[2]>=0 and coord[0]<z_lim and coord[1]<y_lim and coord[2]<x_lim)]

	return tuple(tuple(map(int,i)) for i in zip(*sphere_coord))  #turn list of 3-tuples into a 3-tuple of tuples (Z,Y,X) and convert to int (used for indexing) (tuple->tuple->int)


def get_3D_sphere_coordinates_v2(a_coordZYX=[0,0,0],radius = 10, limit_shape=(9999,9999,9999),xy_vs_z=1,verbose=False):
	'''
	get coordinates of a sphere on a certain locations, and within the perimeter set out by limits
	outputs an 3D array with the coordinates [Z,Y,X]
	:adapted from ellipsoid method from scikit 
	=> improve performance
	'''
	a = radius
	b = radius
	c = radius
	spacing  = (xy_vs_z, 1., 1.)  #Spacing in (x, y, z) spatial dimensions.

	offset = np.r_[1, 1, 1] * np.r_[spacing]

	# Calculate limits, and ensure output volume is odd & symmetric
	low = np.ceil((- np.r_[a, b, c] - offset))
	high = np.floor((np.r_[a, b, c] + offset + 1))

	for dim in range(3):
		if (high[dim] - low[dim]) % 2 == 0:
			low[dim] -= 1
		num = np.arange(low[dim], high[dim], spacing[dim])
		if 0 not in num:
			low[dim] -= np.max(num[num < 0])

	# Generate (anisotropic) spatial grid
	x, y, z = np.mgrid[low[0]:high[0]:spacing[0],
					   low[1]:high[1]:spacing[1],
					   low[2]:high[2]:spacing[2]]
	
	arr = (((x / float(a)) ** 2 + (y / float(b)) ** 2 + (z / float(c)) ** 2) <= 1 ) #to be added = get coordinates (all int, around centerpoint and in the correct format (tuple (Z,Y,X)) +  apply filter for image dimensions
	
	return arr




def examine(na,na_name='no name np-array',output=sys.stdout,ix_print=0):  
	print('-----------EXAMINE ',na_name,'-------------------------------------------------->',file=output) 
	print(na_name,' has shape' ,na.shape,na.dtype,file=output)
	print('the max is ', np.max(na),', the min is ', np.min(na), 'the median is ', np.median(na), 'the mean is ', np.mean(na),file=output)
	print('histogram =>',file=output) 
	[print ('count:{0} for bin:{1}'.format(*i),file=output) for i in zip(*np.histogram(na))]
	print(str(ix_print),'th value is ',file=output)
	print(na[ix_print,...],file=output)
	print('<-----------------------------------------------------------------------------------',file=output)
	return

def examine_img(img,name='no name given'):
	print('examine img name={0} :'.format(name))
	print('shape:' ,img.shape)
	print('dtype:' ,img.dtype)
	print('bincount membrane: ',np.bincount(img.flatten()))
	print("Histogram of DT :{2} values{2} {0} {2} with bin_edges {2}{1}".format(np.histogram(img)[0],np.histogram(img)[1],"\n"))
	return

def slice_tif(na_tif=np.zeros((0)),ix_slice_z=None,ix_slice_time=None,
			  RGB=False,
			  ix_keep_channel=None,
			  filehandler=None,
			  save_selection=False,
			  write_excel_slice=False,
			  trim_1D=True,
			  verbose=True):
	'''
	slice a np-array
	outer dimensions=1 will be erased 
	
	:param ix_slice_z: eg [6,7] ; set to None if you don't want to select (27 zpos)
	:param ix_slice_time:eg [6,7] ;set to None if you don't want to select (18 timesteps)
	:param na_tif: input tif, if not given, then you have to give load info
	:param RGB: 
	:param ix_keep_channel : only use for TZCYX format , eg (21, 25, 2, 768, 768).  will remove the channel dimension
	:param d_load_info: if given the tif will be loaded
	:param d_save_info: if given the output tif will be saved. leave filename '' for automatic filename creation
	:param write_excel_slice: if true, then the slice will be written to excel
	
	:teturn the na_slice
	'''
	
	if not na_tif.any(): 
		if verbose:print('No input data given (or all zeros), so pass on load info via a filehandler object')
		na_tif = filehandler.load_tif(verbose=verbose)
	
	if ix_keep_channel:
		na_tif=na_tif[:,:,ix_keep_channel,:,:]  #channel dimension will automatically be removed
		
	#Select part of tiff file
	if ix_slice_z is not None:
		if len(ix_slice_z)==2:
			ix_slice_z=list(range(ix_slice_z[0],ix_slice_z[1]))

	if ix_slice_time is not None:
		if len(ix_slice_time)==2:
			ix_slice_time=list(range(ix_slice_time[0],ix_slice_time[1]))

	if len(na_tif.shape)==4 or (len(na_tif.shape)==5 and RGB):
		if ix_slice_z  and ix_slice_time:
			na_tif_slice = na_tif[ix_slice_time,ix_slice_z,:,:]
		elif ix_slice_time:
			na_tif_slice = na_tif[ix_slice_time,...]
		elif ix_slice_z:
			na_tif_slice = na_tif[:,ix_slice_z,:,:]
		else:
			if verbose:print('no select')
			na_tif_slice=na_tif
	elif len(na_tif.shape)==3 or (len(na_tif.shape)==4 and RGB):
		if ix_slice_z:
			na_tif_slice = na_tif[ix_slice_z,:,:]
			if na_tif_slice.shape[0]==1:
				na_tif_slice = na_tif_slice.reshape(na_tif_slice.shape[1],na_tif_slice.shape[2]) 
		else:
			if verbose:print('no select')
			na_tif_slice=na_tif
		
	else:
		print('no selection : dimension is not 3 or 4')
	
	#remove top dimensions of length 1
	if trim_1D:
		for i in range(5,1,-1):
			if len(na_tif_slice.shape)==i:
				if na_tif_slice.shape[0]==1:
					na_tif_slice = na_tif_slice.reshape(*na_tif_slice.shape[1:i])
				elif na_tif_slice.shape[0]==0:
					print('Warning : first dimension = 0, check your slice range !')
	
	if save_selection:
		#determine suffixes for output filenames
		if not ix_slice_z:
			z_suffix = 'z=All'
		else:
			z_suffix = 'z=' + str(ix_slice_z[0]) + 'to' + str(ix_slice_z[-1])
			
		if not ix_slice_time:
			time_suffix = 't=All'
		else:
			time_suffix = 't=' + str(ix_slice_time[0]) + 'to' + str(ix_slice_time[-1])   
	
		#write selection to output tiff
		if filehandler.d_save_info:
			if not filehandler.d_save_info['f_name']:
				filehandler.d_save_info['f_name'] = '{0}_{1}#{2}'.format(z_suffix,time_suffix,filehandler.d_load_info['f_name'])
			filehandler.save_data(na_tif_slice,verbose=verbose)
			
			
	if write_excel_slice:  #write slice to excel (for now, only to root directory)
		filehandler.save_data(na_tif_slice,file_ext='xlsx',verbose=verbose)

	if verbose:print('file sliced to shape={0}'.format(na_tif_slice.shape))   
	
	return na_tif_slice


def convert_16bit_to_RGB(a_img_16bit,clip_not_scale=True,d_save_info=None):

	if a_img_16bit.dtype!='int16':
		a_img_16bit = a_img_16bit.astype('int16')

	a_img_16bit = a_img_16bit - np.min(a_img_16bit)  #subtract mnimimum
	
	if clip_not_scale:
		a_img_8bit=img_as_ubyte(np.where(a_img_16bit>255,255,a_img_16bit))
	else:
		a_img_8bit =img_as_ubyte(a_img_16bit)

	#a_4D_RGB = np.asarray([grey2rgb(i,alpha=None) for i in a_img_8bit]) #scikit1702
	a_4D_RGB = np.asarray([grey2rgb(i) for i in a_img_8bit])
	
	if d_save_info:
		save_data(a_4D_RGB,**d_save_info,file_ext='tif',RGB=True,resolution='uint8',channel_reshape=False)
	
	
	return a_4D_RGB 

def plot3D_stack_with_overlay_1axes(ax,na_3D, overlay_image,
									intensity_threshold,trim_overlay,
									overlay_type,
									a_shapes,
									view_init,
									shape_by_shape,
									number=1,
									label_exterior_cell=1,
									verbose=True):
	z,y,x = na_3D.shape
	
	if shape_by_shape:
		ax.get_xaxis().set_visible(False)
		ax.get_yaxis().set_visible(False)
		ax.get_xaxis().set_ticks([])
		ax.get_yaxis().set_ticks([])
		ax.set_zlim(z, 0) # slice 27 on the bottom
		#ax.get_zaxis().set_visible(False)
	else:
		ax.set_xlabel("x-axis")
		ax.set_ylabel("y-axis")
		ax.set_zlabel("z-axis")
		
	ax.set_xlim(0, x)  
	ax.set_ylim(0, y)  
	ax.set_zlim(z, 0) # slice 27 on the bottom
	   

	#scaling the axis not supported in 3D for the moment.  this oneliner solves this 
	#https://stackoverflow.com/questions/30223161/matplotlib-mplot3d-how-to-increase-the-size-of-an-axis-stretch-in-a-3d-plo/30315313
	Z_DIM,Y_DIM,X_DIM = na_3D.shape
	pixel_width = 0.1294751   #check fiji properties
	voxel_depth = 1.1000000
	xy_vs_z = voxel_depth/pixel_width #1 increase in z corresponds to 5 increases in x or y in absolute distance (= back of the envelope calculation)
	Z_STRETCH = 3 # for plotting purposes it can be handy to stretch the z-dimension so you are better able to see the slices. 1=true size ; 5 = max otherwise it will not fit the page anymore
	ax.get_proj = lambda: np.dot(Axes3D.get_proj(ax), np.diag([X_DIM/Y_DIM,1, Z_DIM/Y_DIM * xy_vs_z * Z_STRETCH,1]))
	
	
	if view_init:ax.view_init(elev=view_init[0],azim=view_init[1])
	
	Z1,Y1,X1 = np.where(na_3D>intensity_threshold)
	if verbose:
		if not len(Z1):print('warning : no points selected with intensity threshold', intensity_threshold, 'the maximum of the stack is ',np.max(na_3D))
	ax.plot3D(X1,Y1,Z1,**d_plot3D_red)

	#show a levelset / overlag image over the original
	if len(overlay_image):
		if overlay_type=='level_set':
			print(overlay_image.dtype)
			if trim_overlay: 
				trim_overlay(na_3D, overlay_image)
			Z2,Y2,X2 = np.where(overlay_image<0) #interior of a level set is below zero
		else:
			Z2,Y2,X2 = np.where(overlay_image!=0)      

			
		ax.plot3D(X2,Y2,Z2,**d_plot3D_green) 

	#draw spheres inside the cell
	if len(a_shapes):
		for ix_stack, shape in enumerate(a_shapes):
			if shape.__class__.__name__=='Sphere':
				ax.plot3D(*shape.l_pixels[::-1],**d_plot3D_sphere,color=a_colors[ix_stack%len(a_colors)])
			if shape.__class__.__name__=='Cell':
				if shape.l_pixels:
					ax.plot3D(*shape.l_pixels[::-1],**d_plot3D_sphere,color=a_colors[ix_stack%len(a_colors)])
					if shape.label==label_exterior_cell:
						s_title = str(shape.label) + '(exterior)'
						ax.set_title(s_title)
					else:
						ax.set_title(shape.label)
			
#     if len(a_spheres):
#         for ix_stack, d_sphere in enumerate(a_spheres):
#             ax.plot3D(*get_3D_sphere_coordinates(d_sphere['centre'], d_sphere['radius'], na_3D.shape, verbose=False)[::-1],
#                       **d_plot3D_sphere,color=a_colors[ix_stack%len(a_colors)])
	
	
	if shape_by_shape:
		pass
	else:
		title =  'intensity_threshold={0};view_init={1}'.format(intensity_threshold,view_init)
		ax.set_title(title, fontsize=12)
	
	if view_init: d_save_info['f_name']= d_save_info['f_name'] + '_elev=' + str(view_init[0]) +'_azim=' + str(view_init[1])
	
	return
	
def plot3D_stack_with_overlay(na_3D,overlay_image=[],
							  intensity_threshold=0,
							  trim_overlay=True,
							  overlay_type='level_set',
							  a_shapes=[],
							  d_save_info=None,
							  view_init=None,
							  shape_by_shape=False,
							  label_exterior_cell=1,
							  verbose=True):
	'''
	
	:param na_3D: cell l_pixels, mandatory. set intensity treshold to '999999' to only show overlay
	:param overlay_image: an overlay as a level set, so floating values, with interior > 0
	:param intensity_threshold: the cut off pixelintensity for the cell l_pixels
	:param trim_overlay: cut the overlay to the same size as na_3D.  leave on, unless you know the sizes are the same
	:param d_save_info: if filled in , image will be saved
	'''
	
	
	dim_axes_grid = ceil(sqrt(len(a_shapes)))
	if (not shape_by_shape) or (dim_axes_grid==1):
		fig = plt.figure(figsize=(20, 20))
		ax = fig.add_subplot(111, projection='3d') # by doing this you get access to the 3D toolkits functions , which contains plot3D for example
		plot3D_stack_with_overlay_1axes(ax,number=1,na_3D=na_3D,overlay_image=overlay_image,intensity_threshold=intensity_threshold,
										trim_overlay=trim_overlay,overlay_type=overlay_type,shape_by_shape=shape_by_shape,
										a_shapes=a_shapes,view_init=view_init,label_exterior_cell=label_exterior_cell)
	else:

		fig, axes = plt.subplots(dim_axes_grid, dim_axes_grid, 
								 figsize=(30, 30),
								 subplot_kw=dict(projection='3d'))
		a_axes = axes.flatten()
		if verbose:print('Starting to plot cell by cell.  this might take a while....')
		for ix_stack in range(0,len(a_shapes)):
			plot3D_stack_with_overlay_1axes(ax=a_axes[ix_stack],number=ix_stack,na_3D=na_3D,overlay_image=[],intensity_threshold=intensity_threshold,
										trim_overlay=trim_overlay,overlay_type=overlay_type,shape_by_shape=shape_by_shape,
										a_shapes=[a_shapes[ix_stack]],view_init=view_init,label_exterior_cell=label_exterior_cell)

			if verbose:print('cell{0} is plotted'.format(ix_stack),end="")
	if d_save_info:save_data(fig,**d_save_info,file_ext='png',verbose=verbose)
		
	print("")
	return fig

def get_color_values(fun_cd='1'):
	'''
	fun_cd = '1' : get a list of reserved 16 bit fiji colors, used for manual annotation
	
	d_reserved_colors : are the reserved colors that can be used for manual annotation
	->  the first value is the Fiji integer that is used for this color, this will be filtered
	->  the second value (tuple) is the RGBA value, this is just for reference here (eg. d_reserved_colors.get('orange')[1][0:3] to get RGB)
	@@deprecated : turns out that Fiji changes the number it assigns to a specific RGB :(
	'''
	d_reserved_colors = {
				'orange':[326,(255,128,0,255)],
				'orange+':[282,(204,102,0,255)],
				'orange++':[260,(153,76,0,255)],
				'orange+++':[247,(130,40,0,255)],
				'orange-':[356,(255,153,51,255)],
				'orange--':[487,(255,178,102,255)],
				'orange---':[505,(255,204,153,255)],
				
				'yellow':[459,(255,255,0,255)],
				'yellow+':[417,(204,204,0,255)],
				'yellow++':[356,(153,153,0,255)],
				'yellow+++':[260,(102,102,0,255)],
				'yellow-':[477,(255,255,51,255)],
				'yellow--':[492,(255,255,102,255)],
				'yellow---':[505,(255,255,153,255)],
				
				'white':[542,(255,255,255,255)],
				# 'green':(0,255,0,255),  # is mapped to zero in Fiji, can not be used
				# 'green+':(0,204,0,255),
				# 'green++':(0,153,0,255),
				# 'green+++':(0,102,0,255),
				# 'green-':[477,(51,255,51,255)], # is mapped to yellow values in Fiji, can not be used
				# 'green--':(102,255,102,255),
				# 'green---':(153,255,153,255),
				
				# 'skyblue': (0,255,255,255),
				# 'skyblue+': [70,(0,204,204,255)],
				'skyblue++': [54,(0,153,153,255)],
				'skyblue+++': [37,(0,102,102,255)],
				# 'skyblue-': (51,255,255,255),
				# 'skyblue--': (102,255,255,255),
				# 'skyblue---': (153,255,255,255),
				
				'blue':[87,(0,0,255,255)],     
				# 'blue+':(0,0,204,255),
				# 'blue++':(0,0,153,255),
				# 'blue+++':(0,0,102,255),
				'blue-':[100,(51,51,255,255)],
				'blue--':[115,(102,102,255,255)],
				# 'blue---':(153,153,255,255),       
			   
			   'purple':[122,(127,0,255,255)],
			   # 'purple+':(102,0,204,255),
			   'purple++':[74,(76,0,153,255)],
			   'purple+++':[50,(51,0,102,255)],
			   'purple-':[129,(153,51,255,255)],
			   'purple--':[140,(178,102,255,255)],
			   # 'purple---':(204,153,255,255),
			   
			   'fuchia':[146,(255,0,255,255)],
			   'fuchia+':[157,(204,0,204,255)],
			   'fuchia++':[186,(153,0,153,255)],
			   'fuchia+++':[63,(102,0,102,255)],
			   'fuchia-':[146,(255,51,255,255)],
			   # 'fuchia--':(255,102,255,255),
			   # 'fuchia---':(255,153,255,255),
			   
			   'magenta':[205,(255,0,127,255)],
			   # 'magenta+':(204,0,102,255),
			   'magenta++':[207,(153,0,76,255)],
			   'magenta+++':[212,(102,0,51,255)],
			   'magenta-':[192,(255,51,153,255)],
			   'magenta--':[175,(255,102,178,255)],
			   'magenta---':[520,(255,153,204,255)],
			   
			   'gray':[179,(128,128,128,255)],
			   'gray+':[56,(96,96,96,255)],
			   'gray++':[19,(64,64,64,255)],
			   'gray+++':[8,(32,32,32,255)],
			   'gray-':[507,(160,160,160,255)],
			   'gray--':[516,(192,192,192,255)],
			   'gray---':[525,(224,224,224,255)],
			   
			   # 'lime':(0,255,0,255),
			   'cyan':[70,(0,204,204,255)],
			   # 'black':(0,0,0,255),
			   
			   'off-white':[529,(230,230,230,255)],
			   'red':[256,(255,0,0,255)],
			   # 'red vague':(255,0,0,15),
			   'silver':[516,(192,192,192,255)],
			   # 'gray':(128,128,128,255),
			   'maroon':[234,(128,0,0,255)],
			   'olive':[277,(128,128,0,255)],
			   'teal':[45,(0,128,128,255)],
			   # 'navy':(0,0,128,255),
			   'beige':[525,(245,245,220,255)],
			   # 'light pink':(255,182,193,255),
			   'deep pink':[192,(255,20,147,255)],
			   'blue violet':[131,(138,43,226,255)],
			   'chocolate':[286,(210,105,30,255)],
			   'salmon':[326,(250,128,114,255)],
			   'gold':[433,(255,215,0,255)],
			   'dark sea green':[503,(143,188,143,255)],
			   'plum':[525,(221,160,221,255)],
			   'lawn green':[459,(124,252,0,255)]
			   }
	
	l_colors = []
	if fun_cd == 1:
		l_colors = [i[0] for i in d_reserved_colors.values()]
		
	return l_colors
	
def filter_tif_values(a_img,l_filter_values=[], filter_reserved_fiji=False,invers=False,remove_small_objects=True):
	'''
	:keep only the values in l_filter_values, all other values will be set to zero

	'''
	
	# filter the values
	if filter_reserved_fiji: #do not use input list, filter all the reserved colors
		l_filter_values = get_color_values(fun_cd=1)
	
	if invers:
		a_img_filtered = np.where(np.isin(a_img,l_filter_values), 0, a_img)
	else:
		a_img_filtered = np.where(np.isin(a_img,l_filter_values), a_img, 0)
	
	
	#try to remove small isolated pixels, check per slice
	examine(a_img_filtered,'a_img_filtered')
	if remove_small_objects:
		a_img_filter2 = np.zeros_like(a_img_filtered)
		for ix_z,a_z_i in enumerate(a_img_filtered):
			a_bool = a_z_i > 0
			a_img_filter2 [ix_z,...]= remove_small_objects(a_bool, min_size=10, connectivity=1, in_place=False)
		
		examine(a_img_filter2,'a_img_filter2')
	else:
		a_img_filter2 = a_img_filtered
	
	return np.where(a_img_filter2,a_img_filtered,0)

def replace_tif_values(a_img,ix_z,old_value,new_value,constrain=None):
	'''
	replaces values in a tif image
	z_ix = index of z_slice.  if = 'all' then apply to whole stack
	constrain = to indicate a special action instead of just a simple replace all operation
		'unconnected' = will only replace values that are not connected = loose pixels
	'''
	
	if constrain=='unconnected':
		a_proc = a_img if ix_z == 'all' else a_img[int(ix_z)]
		a_bool = a_proc==old_value
		a_bool_proc= remove_small_objects(a_bool, min_size=2, connectivity=1, in_place=False)
		a_proc[a_bool != a_bool_proc] = new_value #want to know the pixels that are removed, so the inverse
			
		return a_img

	if ix_z == 'all':
		a_img[a_img==old_value] = new_value
	else:
		a_z = a_img[int(ix_z)]
		a_z[a_z==old_value]=new_value
		
	return a_img
	
def make_image_sequence(input_file,output_folder,prefix='man_seg_000_'):
	'''
	the fiji 'image sequence' can lead to the error :'Incorrect count for "ColorMap" when using Tiffreader, therefore made my own
   
	'''
	filehandler = FileHandler()
	filehandler.d_load_info['f_name'] = str(input_file)
	filehandler.d_save_info['save_dir'] = str(output_folder)
	a_img = filehandler.load_tif()
	for ix_z,a_slice in enumerate(a_img):
		suffix = "{:03d}".format(ix_z)
		filehandler.d_save_info['f_name'] = "{0}{1}".format(prefix,suffix)
		filehandler.save_data(a_slice)
		
	return

def calculate_CTC_SEG_slice_detail(a_gt,a_res,output_folder,ix_z="",verbose=False,output=sys.stdout,d_d_gt_res_overlap={}):
	'''
	calculate SEG score for 2 slices, with subcalculations
	'''
	if output_folder:
		# f_name =  output_folder / "CTC_SEG_detail_calculations.txt"
		# f_output = open( f_name  ,"w" )
		fh = FileHandler()
		fh.d_save_info['save_dir'] = str(output_folder)

	
	a_gt_labels = np.unique(a_gt)
	a_gt_labels =np.delete(a_gt_labels, [0])
	if verbose: print('ground truth labels are ', a_gt_labels,file=output )

	
	SEG_counter = 0
	SEG_sum = 0
	jaccard_sum = 0

	for gt_label in a_gt_labels:
		a_bool_gt = a_gt==gt_label
		nb_gt_label_pix = np.sum(a_bool_gt)
		a_res_overlap = np.where(a_bool_gt,a_res,0)
		a_bincount_overlap = np.bincount(a_res_overlap.flatten())
		a_bincount_overlap[0] = 0
		a_res_labels_overlap = np.unique(a_res_overlap)
		a_res_labels_overlap =np.delete(a_res_labels_overlap, [0])
		if verbose: print('  ground truth label {0} overlaps with res labels {1} '.format(gt_label, a_res_labels_overlap),file=output )
		
		max_overlap_label = 0
		max_overlap_ratio = 0
		max_jaccard       = 0

		d_GT_volumes = d_d_gt_res_overlap.setdefault("GT_volumes",{})  #d_d_gt_res_overlap :  key = gt-label, value = dict (key = RES label, value = nb of overlapping pixels with gt-label)
		nb_pixel_self = d_GT_volumes.setdefault(gt_label, 0)
		d_GT_volumes[gt_label] = nb_pixel_self + np.sum(a_bool_gt)

		d_res_overlap = d_d_gt_res_overlap.setdefault(gt_label,{}) 
		for res_label in a_res_labels_overlap:
			nb_res_label_pix = np.sum(a_res==res_label)
			nb_overlap_pix = a_bincount_overlap[res_label]
			jaccard = (nb_overlap_pix ) / (nb_gt_label_pix + nb_res_label_pix - nb_overlap_pix)
			overlap_ratio = nb_overlap_pix / nb_gt_label_pix
			if overlap_ratio > max_overlap_ratio:
				max_overlap_label = res_label
				max_overlap_ratio = overlap_ratio
				max_jaccard = jaccard
				
			rec = "    GT_label {0} ({4} pixels)  overlaps with RES_label {1}({5} pixels) with {6} pixels overlap.  \n    Jaccard = {3} with ratio (nb pix overlap/nb pix GT) = {2}\n".format(
				gt_label,
				res_label,
				(nb_overlap_pix/nb_gt_label_pix),
				jaccard,
				nb_gt_label_pix,
				nb_res_label_pix,
				nb_overlap_pix
				)

			nb_pixel_overlap = d_res_overlap.setdefault(res_label, 0) 
			d_res_overlap[res_label] += nb_overlap_pix

			if verbose:print(rec,file=output)
			
		SEG_score = 0
		 
		if max_overlap_ratio < 0.5:
			SEG_score = 0
			rec = "  SEG_score = 0 because max_overlap_ratio < 0.5 !"
		else:
			SEG_score = max_jaccard
			rec = "  SEG_score = {0} for gt_label {1} overlapping res_label {2} ".format(SEG_score,gt_label,max_overlap_label)
				
		SEG_counter += 1
		SEG_sum += SEG_score
		jaccard_sum +=max_jaccard
			
		if verbose:print(rec,file=output)

		#write visual    
		if output_folder:
			# f_output.write(rec)
			fh.d_save_info['f_name'] = "GT_{0}(z{5})_vs_RES_{1}_SEG={2}|jac={3}|overlap={4}".format(
				gt_label, 
				max_overlap_label,
				round(SEG_score,2),
				round(max_jaccard,2),
				round(max_overlap_ratio,2),
				str(ix_z).zfill(3)
				)
				
			a_output = np.zeros_like(a_gt)
			# a_output [a_gt == gt_label] = gt_label
			# a_output [a_res == max_overlap_label] = max_overlap_label
			# a_output [(a_gt == gt_label) & (a_res == max_overlap_label)] = np.round(SEG_score*100)
			a_output [a_gt == gt_label] = 1
			a_output [a_res == max_overlap_label] = 2
			a_output [(a_gt == gt_label) & (a_res == max_overlap_label)] = 3
			fh.save_data(a_output,verbose=verbose)
		
	
	# if output_folder:f_output.close()
	a_res_labels = np.unique(a_res)
	a_res_labels =np.delete(a_res_labels, [0])

	if list(a_gt_labels): # a slice with no GT is skipped (=SEG score rules)
		for res_label in a_res_labels:
			d_RES_volumes = d_d_gt_res_overlap.setdefault("RES_volumes",{})  #d_d_gt_res_overlap :  key = gt-label, value = dict (key = RES label, value = nb of overlapping pixels with gt-label)
			nb_pixel_self = d_RES_volumes.setdefault(res_label, 0)
			d_RES_volumes[res_label] = nb_pixel_self + np.sum(a_res==res_label)
			if res_label==2:
				print(np.sum(a_res==res_label),'for a total of ', d_RES_volumes[res_label] )
	
	return SEG_sum,jaccard_sum,SEG_counter,d_d_gt_res_overlap
	
	
def hessian(x,spacing=[1,1,1],verbose=False):
	"""
	Calculate the hessian matrix with finite differences
	Parameters:
	   - x : ndarray
	Returns:
	   an array of shape (x.dim, x.ndim) + x.shape
	   where the array[i, j, ...] corresponds to the second derivative x_ij
	"""
	x = x.astype('float64')   #important : calculations are completely wrong using uint16 !!  (gradients of 4000 eg)
	x_grad = np.gradient(x,*spacing)
	
	if verbose:
		a_grad = np.array(x_grad)
		#examine(a_grad,"a_grad is de eerste gradient")
		print('spacing = ',spacing)
		print('the 3 gradients (first gradient) for point 10,200,100 (',x[10,200,100],') are a_grad[:,10,200,100] \n',a_grad[:,10,200,100])
		print('the neighbouring points (XY plane) are x[10,199:202,99:102] =\n',x[10,199:202,99:102])
		print('the neighbouring points (ZY plane) are x[9:12,199:202,100] =\n',x[9:12,199:202,100])
		x_grad_test = np.gradient(x[10,199:202,99:102])
		a_grad = np.array(x_grad_test)
		print('gradient van XYplane x[10,199:202,99:102] = \n=', a_grad)

		# print('the neighbouring points in z are 10,200,100 are x[9,200,100] =',x[9,200,100])
		# print('the neighbouring points in z are 10,200,100 are x[10,200,100] =',x[10,200,100])
		# print('the neighbouring points in z are 10,200,100 are x[11,200,100] =',x[11,200,100])
		
	hessian = np.empty((x.ndim, x.ndim) + x.shape, dtype=x.dtype) 
	for k, grad_k in enumerate(x_grad):
		# iterate over dimensions
		# apply gradient again to every component of the first derivative.
		tmp_grad = np.gradient(grad_k,*spacing) 
		for l, grad_kl in enumerate(tmp_grad):
			hessian[k, l, :, :] = grad_kl
	return hessian
	
def extract_hessian_info(a_hessian,verbose=False):
	"""
	given a hessian for a point (z,y,x), extract 4 values of interest
	1,2,3 : the 3 eigenvalues from large to small (in absolute value)
	4 : the cos(z-angle) of the eigenvector associated with the largest eigenvalue
	"""
	a_hessian_info = np.zeros((4,))
	lambdas, v = np.linalg.eig(a_hessian)
	lambdas = np.absolute(lambdas)
	v_max = v[:,np.argmax(lambdas)] 
	a_hessian_info[0:3] = sorted(lambdas, reverse=True)
	a_hessian_info[3] = np.absolute(np.dot(v_max,np.array([1,0,0])))  #cosine of the z-angle (without sign)
	
	return a_hessian_info
	
def extract_hessian_info_stack(a_stack,spacing=[1,1,1],verbose=False):
	"""
	for all points in stack calculate 4 values of interest
	1,2,3 : the 3 eigenvalues from large to small (in absolute value)
	4 : the cos(z-angle) of the eigenvector associated with the largest eigenvalue
	"""
	z,y,x = a_stack.shape
	a_hessian = hessian(a_stack,spacing=spacing,verbose=verbose)
	a_hessian_info = np.zeros((4,z,y,x))
	for z_i in range(0,z):
		for y_i in range(0,y):
			for x_i in range(0,x):
				a_hessian_info[:,z_i,y_i,x_i] = extract_hessian_info(a_hessian[:,:, z_i, y_i, x_i],verbose=False)
  
	return a_hessian_info
	
	
def score_hessian_info(a_hessian_info,func="simple",verbose=False):
	"""
	for all points in stack calculate 4 values of interest
	1,2,3 : the 3 eigenvalues from large to small (in absolute value)
	4 : the cos(z-angle) of the eigenvector associated with the largest eigenvalue
	"""
	a_lambda1,a_lambda2,a_lambda3,a_z_angle = a_hessian_info
	np.clip(a_lambda2,0.01,65000)  #remove zeros
	np.clip(a_lambda3,0.01,65000)  #remove zeros
	
	if func=="simple":
		a_score = a_z_angle * a_lambda1/(a_lambda2*a_lambda3)
	if func=="frangi":
		a_Rb = a_lambda1/(a_lambda2*a_lambda3)
		a_Ra = a_lambda2 / a_lambda3
		a_S =  np.sqrt(np.square(a_lambda1) + np.square(a_lambda2) + np.square(a_lambda3) )
		a_score = (1- np.exp(np.square(a_Ra))) * np.exp(np.square(a_Rb)) * (1 - np.exp(-np.square(a_S)))
		
	return a_score

def interpollate_image(a_img,spacing_ZYX = [1.1,0.1294,0.1294],XY_downsampling_ratio = 1,verbose=False):
	"""
	interpollates an image so that the resulting relative spacing  = (1,1,1)
	default will cause upsampling (dimension with highest resolution remains unchanged).  XY_downsampling_ratio can be used to downscale XY however
	"""
	spacing_ZYX = [i / min(spacing_ZYX) for i in spacing_ZYX]  #normalize e.g. (10,1,1) will upsample z by factor 10
	l_zoom = [spacing_ZYX[0] , spacing_ZYX[1] * XY_downsampling_ratio , spacing_ZYX[2] * XY_downsampling_ratio]
	a_img_zoom = zoom( a_img, l_zoom , order=1, prefilter=False ) 
	
	if verbose: print (" - resampling from ", a_img.shape, "to", a_img_zoom.shape, "using spacing_ZYX=",spacing_ZYX)


	return a_img_zoom

def zoom_data(nd_array, spacing=None,target_shape = None, keep_XY_ratio=True, target_res = None, curr_res = None):
	"""
	specify a target shape OR a target resolution in micron
	use 'spacing' to target spacing of the z-axis specifically
	"""

	if len(nd_array.shape) ==4: #handle a timelapse recursively
		l_array_zoom= []
		for ix,stack in enumerate(nd_array):
			array_zoom, [z_zoom,y_zoom,x_zoom], spacing = zoom_data(stack,spacing=spacing,target_shape = target_shape, keep_XY_ratio=keep_XY_ratio, target_res = target_res, curr_res = curr_res)
			l_array_zoom.append(array_zoom)
		return np.array(l_array_zoom).astype('int16'), [z_zoom,y_zoom,x_zoom], spacing
		
	elif spacing is not None:
		print ("start interpollating the data to expand z-range: ")
		# xy_down = 1/4
		xy_down = 1
		print  (" - to compress size x and y are downscaled with a factor of {0}".format(xy_down))
		z_zoom = spacing[0]
		y_zoom = spacing[1] * xy_down
		x_zoom = spacing[2] * xy_down
		spacing=[1,1*xy_down,1*xy_down]
		print  (" - zoom used {0} with spacing {1} ".format([z_zoom,y_zoom,x_zoom], spacing))
	elif target_shape is not None:
		z,y,x = nd_array.shape
		z_zoom = target_shape[0] / z
		y_zoom = target_shape[1] / y
		if keep_XY_ratio:
			x_zoom = y_zoom
		else:
			x_zoom = target_res[2] / x
		print  (" - zoom used {0} with target shape {1} ".format([z_zoom,y_zoom,x_zoom], target_shape))
	elif target_res is not None:
		if curr_res is None:
			print('current resolution ZYX in microns is needed as input')
			return None
		z_zoom,y_zoom,x_zoom = [i/j for i,j in zip(curr_res,target_res)]
		print  (" - zoom used {0} with target resolution {1} ".format([z_zoom,y_zoom,x_zoom], target_res))
	else:
		print("Specify target resolution or shape !")
		return None
	
	nd_array_zoom = ndimage.zoom( nd_array, [z_zoom,y_zoom,x_zoom], order=1, prefilter=False ) 
	
	print  (" - resampling from ", nd_array.shape, "to", nd_array_zoom.shape)
	
	return nd_array_zoom, [z_zoom,y_zoom,x_zoom], spacing

def dds_get_value(dds,label_ds,k1,k2,fallback=None):
	#labels_ds = name of dict = name of excel tab
	# k1 = the key of the dict. in case of multiindex this is a tuple ('key1,'key2')
	# k2 = the name of the value (=the name of the column header)
	try:
		return dds[label_ds][k1][k2]
	except KeyError as ex:
		if fallback:
			return fallback
		else:
			return 'not in dds : {0}'.format(k1)

def agg_files(input_folder,output_folder,file_regex,tc_regex=None,l_suffix=[],verbose=False):
	""" 
	the input_folder will be searched for files complying to the file_regex
	these files will be copied to the output, but the name of these files will get a suffix
	this suffix will be automatically composed by searching the path name with the tc_regex
	Each capturing group will be used in the output file name preceded with the l_suffix label
	e.g. 
	if tc_regex = 'fmn(\d\d)(\d\d)' and l_suffix = ['R','t'] and fmn0109 is found in the path, 
	then the output file name will get _R01_t09 as suffix
	if no tc_regex is given a simple counter will be used (0-indexed)
	"""
	for ix_path,path_i in enumerate(input_folder.rglob(file_regex)):
		if verbose:print('found file : ',path)
		file_name  = path_i.stem
		if tc_regex:
			pattern = re.compile(tc_regex)
			match=pattern.search(path_i)
			if match: # match will be None if nothing is found
				for ix in range(pattern.groups):
					file_name += "_{}{}".format(l_suffix[ix],match.group(ix+1))
			else:
				file_name += "_{}".format(ix_path)
		else:
			file_name += "_{}".format(ix_path)

		copyfile(path_i, output_folder / file_name)

	return










