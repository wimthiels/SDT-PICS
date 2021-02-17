'''
Created on 30 nov 2019 
extracts a selection of a vtp file based on a query
		
@author: wimth
'''
import vtk
import vtk.numpy_interface.dataset_adapter as dsa
import vtk.numpy_interface.algorithms as algos
import numpy as np
import os
import re
import pyvista as pv
import trimesh
import pprint
pp = pprint.PrettyPrinter(indent=4)



def read_vtp_file(vtp_file):

	vtp_reader = vtk.vtkXMLPolyDataReader()
	vtp_reader.SetFileName(str(vtp_file))
	vtp_reader.Update()

	return vtp_reader.GetOutput()   #vtkPolyData

def copy_poly_minimal(poly_in,scale=1):
	""" a copy of poly_in only keeping the minimal data (needed for mesh construction, STL generation etc), With optional scaling"""
	wdo_in = dsa.WrapDataObject(poly_in)
	wdo_out = dsa.WrapDataObject(vtk.vtkPolyData()) # we need a new poly, because updating values through dsa API not possible

	wdo_out.Points = wdo_in.Points * scale

	wdo_out.GetAttributes(vtk.vtkDataObject.POINT).append(wdo_in.PointData['x'] * scale ,'x')
	wdo_out.GetAttributes(vtk.vtkDataObject.POINT).append(wdo_in.PointData['localIndex'] ,'localIndex')
	wdo_out.GetAttributes(vtk.vtkDataObject.POINT).append(wdo_in.PointData['normal'] ,'normal')

	wdo_out.GetAttributes(vtk.vtkDataObject.CELL).append(wdo_in.CellData['x'] * scale ,'x')
	wdo_out.GetAttributes(vtk.vtkDataObject.CELL).append(wdo_in.CellData['localIndex'] ,'localIndex')
	wdo_out.GetAttributes(vtk.vtkDataObject.CELL).append(wdo_in.CellData['normal'] ,'normal')
	wdo_out.GetAttributes(vtk.vtkDataObject.CELL).append(wdo_in.CellData['vertexIndices'] ,'vertexIndices')
	#all other attributes are left unscaled, not important for stl conversion

	#Construct Polygons (Cells) 
	#wdo_out.Polygons    = wdo_in.Polygons[a_ix_point] => not supported by wrapper, so use native VTK
	vtkCellArray = vtk.vtkCellArray()
	for triangle_i in wdo_in.CellData['vertexIndices'].__array__():
		vtkCellArray.InsertNextCell (3, triangle_i)
	wdo_out.VTKObject.SetPolys(vtkCellArray)    

	
	return wdo_out.VTKObject

	
def write_stl_file(poly_in, path_file,scale=1):

	if scale!=1:
		poly_in = copy_poly_minimal(poly_in,scale=scale)

	vtkSTLwriter = vtk.vtkSTLWriter()
	vtkSTLwriter.SetFileTypeToASCII ()
	vtkSTLwriter.SetInputData(poly_in)
	vtkSTLwriter.SetFileName(str(path_file))
	
	return vtkSTLwriter.Write()
	
def write_vtp_file(polyout, path_file):

	vtk_writer = vtk.vtkXMLPolyDataWriter()
	vtk_writer.SetDataMode(vtk.vtkXMLWriter.Ascii)
	vtk_writer.SetByteOrderToBigEndian()
	vtk_writer.SetFileName(path_file)
	vtk_writer.SetInputData(polyout)  
		
	return vtk_writer.Write()


def read_stl_file(p_stl):

	reader = vtk.vtkSTLReader()
	reader.SetFileName(str(p_stl))
	reader.Update()
	poly = reader.GetOutput()
	
	return poly

	
def set_attributetype(field_type):

	if field_type == "CELL":
		attributeType = vtk.vtkDataObject.CELL
	elif field_type == "POINT":
		attributeType = vtk.vtkDataObject.POINT
	elif field_type == "ROW":
		attributeType = vtk.vtkDataObject.ROW
	else:
		raise RuntimeError ("Unsupported field attributeType".format(field_type))
		
	return attributeType
	
def _create_id_array(dataobject, attributeType):
	"""Returns a VTKArray or VTKCompositeDataArray for the ids"""
	if not dataobject:
		raise RuntimeError ("dataobject cannot be None")
	if dataobject.IsA("vtkCompositeDataSet"):
		ids = []
		for ds in dataobject:
			ids.append(_create_id_array(ds, attributeType))
		return dsa.VTKCompositeDataArray(ids)
	else:
		return dsa.VTKArray(\
				np.arange(dataobject.GetNumberOfElements(attributeType)))
				
def maskarray_is_valid(maskArray):
	"""Validates that the maskArray is either a VTKArray or a
	VTKCompositeDataArrays or a NoneArray other returns false."""
	return maskArray is dsa.NoneArray or \
		isinstance(maskArray, dsa.VTKArray) or \
		isinstance(maskArray, dsa.VTKCompositeDataArray)
		
def compute(wdo_in, expression, ns=None):
	#  build the locals environment used to eval the expression.
	mylocals = {}
	if ns:
		mylocals.update(ns)
		
	try:
		mylocals["points"] = wdo_in.Points  # the points are added in case selections is based on the points
	except AttributeError: pass
	
	return eval(expression, globals(), mylocals)  #boolean mask based on query
	
def get_arrays(attribs):
	"""Returns a 'dict' referring to arrays in dsa.DataSetAttributes
	"""
	arrays = {}
	for key in attribs.keys():
		# varname = paraview.make_name_valid(key)
		arrays[key] = attribs[key]


def get_point_array(poly,l_dim=None):
	""" get back point array = shape(nb_cells, connnectivity, 3) #output will be transformed to zyx format
	if l_dim is filled in, the coordinates get transformed into pixel coordinates (eg l_dim = [1e-6,0.129e-6,0.129e-6])
	these coordinates are indexes, because flooring is applied (e.g. all pixels with z-location between 0 and 1e-6 will get z=0)
	"""
	try:
		pvw = pv.wrap(poly)
	except:
		pvw = poly
	
	conn = pvw.extract_cells(0).n_points
	a_points = np.zeros((pvw.n_cells,conn,3))
	for i in range(pvw.n_cells):
		points_cell= pvw.extract_cells(i).points
		for j in range(conn):
			a_points[i,j,:] = points_cell[j]
	
	a_points = a_points[:,:,[2,1,0]] # from xyz to zyx order
	if l_dim:
		a_points= np.floor(a_points / np.array(l_dim)).astype(int)  #turns metric coordinates into pixel coordinates

	return a_points
	
def get_data_array(poly, field_type="CELL", attribute='parentIndex',verbose=True):
	""" get back a point or cell data array """
	if field_type == "CELL":
		data_array = dsa.WrapDataObject(poly).CellData[attribute]
	elif field_type == "POINT":
		data_array = dsa.WrapDataObject(poly).PointData[attribute]

	if isinstance(data_array,vtk.numpy_interface.dataset_adapter.VTKNoneArray):
		if verbose:print("{} is no attribute in this vtk object (field_type={})".format(attribute,field_type))
		return None

	if data_array.__str__().startswith('vtkString'): #
		data_array = pv.convert_string_array(data_array)

	return data_array
	
	
def add_array(poly, a_added, name_array, field_type="CELL",dtype='int'):

	if a_added.dtype.type is np.str_:
	#if dtype =='str':
		a_added = pv.convert_string_array(a_added)
		a_added.SetName(name_array)
		if field_type == "CELL":
			poly.GetCellData().AddArray(a_added);
		elif field_type == "POINT":
			poly.GetPointData().AddArray(a_added);

	else:
		wdo  = dsa.WrapDataObject(poly)
		if field_type == "CELL":
			wdo.GetAttributes(vtk.vtkDataObject.CELL).append(a_added,name_array)
		elif field_type == "POINT":
			wdo.GetAttributes(vtk.vtkDataObject.POINT).append(a_added,name_array)
			
	return poly

def add_array_with_mapper(poly,name_array_source, name_array_dest,mapper={},field_type="CELL",default=999):
	""" takes an array (name_array_source) and maps  these values to a new array (name_array_dest) using a dict as mapper"""
	a_source = get_data_array(poly, field_type=field_type, attribute=name_array_source)
	a_dest = np.vectorize(mapper.get)(a_source,default)
	poly = add_array(poly,a_dest, name_array_dest, field_type=field_type)

	return a_dest


def extract_selection_vtp(poly_in, query=None, l_cell_idx = None, field_type="CELL", inverse_selection=False, verbose=False):
	"""Returns a vtk polydata object as a subselection of an input polydata object
	Some words about indexes : 
	Cells
		- CellId              = index (of CellData and Polys (=connectivity and offset arrays)) = PK
		- vtkOriginalCellIds  = index of original vtp (never changed, even after consecutive selections)
		- SelectionCellIds    = index of previous selection(changes every selection step, allows you to go 1 level up)
		- vertexIndices       = related point indices (=redundant = same as connectivity) = FK
	Points  
		- PointId             = index (of PointData and Points) = PK
		- vtkOriginalPointIds = index of original vtp (never changed, even after consecutive selections) 
		- SelectionPointIds   = index of previous selection(changes every selection step, allows you to go 1 level up)
		
	naming chosen as to comply with the paraview naming
		
	"""
	a_ix_cell = np.array([])
	a_ix_point = np.array([])
	wdo_in  = dsa.WrapDataObject(poly_in)
	wdo_out = dsa.WrapDataObject(vtk.vtkPolyData()) 
	   
	if query:
		attributeType =  set_attributetype(field_type)
		if verbose: print(wdo_in.GetAttributes(attributeType))
		# d_attr_data = get_arrays(inputs[0].GetAttributes(attributeType))
		attribs = wdo_in.GetAttributes(attributeType)
		
		d_attr_data = {key:attribs[key] for key in attribs.keys()} #dict : key = 'attributename' , value = array of data 
		if not "id" in d_attr_data.keys() and re.search(r'\bid\b', query): # add "id" array if the query string refers to id.
			d_attr_data["id"] = _create_id_array(wdo_in, attributeType)
		try:
			maskArray = compute(wdo_in, query, ns=d_attr_data)   # SELECTION :returns boolean mask 
		except:
			print ("Error: Failed to evaluate Expression '%s'. "\
				"The following exception stack should provide additional developer "\
				"specific information. This typically implies a malformed "\
				"expression. Verify that the expression is valid.\n", query)
			raise
			
		if not maskarray_is_valid(maskArray):
			raise RuntimeError("Expression '%s' did not produce a valid mask array. The value "\
				"produced is of the type '{0}'. This typically implies a malformed "\
				"expression. Verify that the expression is valid. {1}".format(query, type(maskArray)))
		
		if inverse_selection:
			maskArray = algos.logical_not(maskArray)

		
		nonzero_indices =  algos.flatnonzero(maskArray)  # returns indices of boolean mask 

		if field_type=="CELL":
			a_ix_cell = np.array(nonzero_indices)
		else:
			print("only cell based selections are supported right now.")
			return
	else:
		a_ix_cell = np.array(l_cell_idx)
	
	if not isinstance(a_ix_cell,np.ndarray) or a_ix_cell.size == 0:
		print('warning : nothing selected.  None value will be returned')
		return None
	
	if field_type=="CELL":
	
		#STEP1 : Replace CellData
		nb_arrays = wdo_in.GetCellData().GetNumberOfArrays()
		if verbose:print("{0} arrays are present in {1}Data".format(nb_arrays,field_type))
		
		
		for i in range(nb_arrays):
			
			vtk_array        = wdo_in.GetAttributes(vtk.vtkDataObject.CELL).GetArray(i)
			attr_name        = wdo_in.GetAttributes(vtk.vtkDataObject.CELL).GetArrayName(i)
			if attr_name=="SelectionCellIds":
				continue
			
			if a_ix_cell.size==0:
				vtk_array_select = np.array([])
			else:
				if isinstance(vtk_array.GetValue(0),str):
					vtk_array_select = np.array([vtk_array.GetValue(i) for i in a_ix_cell])  #indexing not possible on vtkStringArray
				else:
					vtk_array_select = vtk_array[a_ix_cell]
			
			if attr_name=="vertexIndices":
				a_vtkOriginalPointIds = vtk_array_select.__array__()  #not stored in output_poly, only used further down
				continue
		   
			if verbose:print("{0}),{1},{2} ==> {3}".format(i, attr_name, vtk_array.size, vtk_array_select.size))
			
			# wdo_out.GetAttributes(vtk.vtkDataObject.CELL).append(vtk_array_select,attr_name)
			if isinstance(vtk_array_select[0],str):
				add_array(wdo_out.VTKObject, vtk_array_select, attr_name, field_type="CELL",dtype='str')
			else:
				wdo_out.GetAttributes(vtk.vtkDataObject.CELL).append(vtk_array_select,attr_name)
				if verbose:print(wdo_out.GetAttributes(vtk.vtkDataObject.CELL)[attr_name],"compared to input : \n", wdo_in.GetAttributes(vtk.vtkDataObject.CELL)[attr_name])
					
		wdo_out.GetAttributes(vtk.vtkDataObject.CELL).append(a_ix_cell,"SelectionCellIds")  #backup selectionIds to easily refer to the selection 1 level up
		if isinstance(wdo_out.CellData['vtkOriginalCellIds'], dsa.VTKNoneArray):  # at first selection, this column is newly added
			wdo_out.GetAttributes(vtk.vtkDataObject.CELL).append(a_ix_cell,"vtkOriginalCellIds") 
		
		#STEP2 : Get points to be selected based on the cell selection
		a_ix_point = np.unique(a_vtkOriginalPointIds) #unique gives 1D SORTED ascending array
		d_oldPointId_newPointID = { old:new for new,old in enumerate(a_ix_point) }  
		
		#STEP3 : Copy PointData
		nb_arrays = wdo_in.GetPointData().GetNumberOfArrays()
		if verbose:print("{0} arrays are present in CellData".format(nb_arrays))
		for i in range(nb_arrays):
			vtk_array        = wdo_in.GetAttributes(vtk.vtkDataObject.POINT).GetArray(i)
			attr_name        = wdo_in.GetAttributes(vtk.vtkDataObject.POINT).GetArrayName(i)
			if attr_name == "SelectionPointIds":
				continue
				
			if a_ix_point.size==0:
				vtk_array_select = np.array([])
			else:   
				vtk_array_select = vtk_array[a_ix_point]
			if verbose:print("{0}),{1},{2} ==> {3}".format(i, attr_name, vtk_array.size, vtk_array_select.size))

			wdo_out.GetAttributes(vtk.vtkDataObject.POINT).append(vtk_array_select,attr_name)   
			
		wdo_out.GetAttributes(vtk.vtkDataObject.POINT).append(a_ix_point,"SelectionPointIds")  #backup original ids as extra column
		if isinstance(wdo_out.PointData['vtkOriginalPointIds'], dsa.VTKNoneArray): # at first selection, this column is newly added
			wdo_out.GetAttributes(vtk.vtkDataObject.POINT).append(a_ix_point,"vtkOriginalPointIds")
		
		#STEP4: Copy Points
		if a_ix_point.size!=0: wdo_out.Points = wdo_in.Points[a_ix_point]

		#STEP5: Construct Polygons (Cells) 
		#wdo_out.Polygons    = wdo_in.Polygons[a_ix_point] => not supported by wrapper, so use native VTK
		vtkCellArray = vtk.vtkCellArray()
		vertexIndices_new = []
		for p1,p2,p3 in a_vtkOriginalPointIds:
			l_new_triangle = [d_oldPointId_newPointID[p1],d_oldPointId_newPointID[p2],d_oldPointId_newPointID[p3]]
			vtkCellArray.InsertNextCell (3, l_new_triangle) 
			vertexIndices_new.append(l_new_triangle)
			
		wdo_out.VTKObject.SetPolys(vtkCellArray)
		
		#STEP6: update CellData vertexIndices
		wdo_out.GetAttributes(vtk.vtkDataObject.CELL).append(np.array(vertexIndices_new),"vertexIndices")
		
	if a_ix_point.size!=0:
		return wdo_out.VTKObject 
	else:
		return None
	
	
def get_aggregate_data(poly,a_cell_id,l_aggregate_columns,func="sum"):
	wdo = dsa.WrapDataObject(poly)
	l_output = []
	
	if func=="sum":
		for ix, column in enumerate(l_aggregate_columns):
			l_output.append( np.sum(wdo.CellData[l_aggregate_columns[ix]][a_cell_id].__array__()) )
	else:
		print('func not supported')

	return l_output


def extract_parentID_selection_from_VTP(poly_in,verbose=False):
	d_parentID_VTP={}
	a_parent_id = np.unique(dsa.WrapDataObject(poly_in).CellData['parentIndex'].__array__())
	for parent_id in a_parent_id:
		if verbose:print('--processing cell with parentID {0}'.format(parent_id),flush=True)
		poly_out = extract_selection_vtp(poly_in, query="parentIndex == {0}".format(parent_id) ,field_type="CELL")
		d_parentID_VTP[parent_id] = poly_out

	return d_parentID_VTP

def enrich_embryo_with_contactID(poly_in,d_parentID_map=None):
	"""
	d_parentID maps the standard parentID into a canonical parent ID

	"""
	a_parentIndex = dsa.WrapDataObject(poly_in).CellData['parentIndex'].__array__()
	a_contactIndex = dsa.WrapDataObject(poly_in).CellData['contact_index'].__array__()
	a_contact_id = np.array([a_parentIndex[i] for i in a_contactIndex])
	a_contact_id = np.where(a_contactIndex==0,a_parentIndex,a_contact_id)  #contactindex=0 means no contact => contact_id = parent_id in that case
	poly_out = add_array(poly_in,a_contact_id, "contactId", field_type="CELL")

	if d_parentID_map:
		a_parentID_canonical = np.array([d_parentID_map[x] for x in a_parentIndex])
		a_contactID_canonical = np.array([d_parentID_map[x] for x in a_contact_id])
		poly_out = add_array(poly_out,a_parentID_canonical, "parentId_canonical", field_type="CELL")
		poly_out = add_array(poly_out,a_contactID_canonical, "contactId_canonical", field_type="CELL")


	return poly_out


def get_parentIDs(poly_in):
	a_parentIndex = dsa.WrapDataObject(poly_in).CellData['parentIndex'].__array__()
	return np.sort(np.unique(a_parentIndex))

def get_unique_values(poly_in,attrib='parentIndex'):
	a_attrib= dsa.WrapDataObject(poly_in).CellData[attrib].__array__()
	return np.sort(np.unique(a_attrib))



def convert_to_line_mesh(poly_in):
	"""returns same poly mesh, but with cell_type = lines instead of triangles"""
	extractEdges = vtk.vtkExtractEdges()
	extractEdges.SetInputData(poly_in)
	extractEdges.Update()
	poly_line = extractEdges.GetOutput()  
	wdo_line = dsa.WrapDataObject(poly_line)
	#a_edges =  poly_line.GetLines()  #vtkCellArray of lines
	return poly_line


def get_outside_points(poly_in):
	"""returns outside points of a triangulated mesh. 
	These are points associated with lines that are only part of one triangle"""

	#step1 : build ds that counts occurence of lines (=2-sets) in triangles
	from collections import defaultdict
	d_line_count=defaultdict(int)
	wdo_in = dsa.WrapDataObject(poly_in)
	a_connectivity = wdo_in.CellData['vertexIndices'].__array__() #=the connectivity array
	for triangle_i in a_connectivity:
		triangle_i.sort()
		line1 = [triangle_i[0],triangle_i[1]]
		line2 = [triangle_i[0],triangle_i[2]]
		line3 = [triangle_i[1],triangle_i[2]]
		for line in [line1,line2,line3]:
			d_line_count[tuple(line)] += 1

	#step2 : get indices of outside points
	s_outside_points = set()
	for line_i, count in d_line_count.items():
		if count==1:
			for point_i in line_i:
				s_outside_points.add(point_i)

	#step3 : get xyz of outside points + poly
	a_points_xyz = wdo_in.Points.__array__()[np.array(list(s_outside_points))]
	wdo_border_points = dsa.WrapDataObject(vtk.vtkPolyData())
	wdo_border_points.SetPoints(a_points_xyz)


	return a_points_xyz,wdo_border_points.VTKObject


def construct_plane(point_plane,normal_plane):
	vtk_plane = vtk.vtkPlaneSource()
	vtk_plane.SetCenter(point_plane[0],point_plane[1],point_plane[2])
	vtk_plane.SetNormal(normal_plane[0],normal_plane[1],normal_plane[2])
	vtk_plane.Update()
	return vtk_plane.GetOutput()


def scale_poly(poly_in,scale=1.e6):
	'''the shape of the poly is scaled, all other data remains unaffected'''
	wdo_in = dsa.WrapDataObject(poly_in)
	a_points_xyz = wdo_in.Points.__array__()
	wdo_in.SetPoints(a_points_xyz * scale )

	return wdo_in.VTKObject

def poly_to_trimesh(poly_in):
    """convert a vtk poly to a trimesh object"""
    
    wdo_in = dsa.WrapDataObject(poly_in)

    return trimesh.Trimesh(vertices=wdo_in.Points,faces=wdo_in.CellData['vertexIndices'])

def list_attributes(poly,repr=True):
	""" returns the attribute names in a list (point and cell attributes)"""
	wdo = dsa.WrapDataObject(poly)
	l_point_attr = [wdo.GetAttributes(vtk.vtkDataObject.POINT).GetArrayName(i) for i in range(wdo.GetPointData().GetNumberOfArrays())]
	l_cell_attr = [wdo.GetAttributes(vtk.vtkDataObject.CELL).GetArrayName(i) for i in range(wdo.GetCellData().GetNumberOfArrays())]

	if repr:
		print('POINT attributes -> \n',l_point_attr)
		print('CELL attributes -> \n' ,l_cell_attr)

	return l_point_attr,l_cell_attr

def list_unique_values(poly,attribute='cellName', field_type="CELL",repr=True):
	""" returns the attribute names in a list (point and cell attributes)"""
	wdo = dsa.WrapDataObject(poly)
	data_array = get_data_array(poly, field_type=field_type, attribute=attribute)
	if data_array is None:
		return []
	else:
		a_unique = np.unique(data_array)
		print(f"Unique values in {attribute} are {a_unique}")
		return a_unique


def get_mapper(poly,name_array_key, name_array_value,field_type="CELL"):
	"""finds the all combinations between two arrays and provides a mapper between them"""
	a_key = get_data_array(poly, field_type=field_type, attribute=name_array_key)
	if a_key is None:
		return None

	a_value = get_data_array(poly, field_type=field_type, attribute=name_array_value)
	if a_value is None:
		return None

	s_pairs = set([(i,j) for i,j in zip(a_key,a_value)])
	mapper = {i:j for i,j in s_pairs}

	return mapper