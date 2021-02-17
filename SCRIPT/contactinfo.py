'''
Created on 29 nov 2019 
given a set of meshes determine which cells are making contact and determine the contact_area
        
@author: wimth
'''


import vtk
from VTK_utils import extract_selection_vtp, read_vtp_file, write_vtp_file, write_stl_file, get_aggregate_data, add_array
from vtk.util.numpy_support import vtk_to_numpy
from vtk.numpy_interface import dataset_adapter as dsa
from vtk.numpy_interface import algorithms as algs

import trimesh
from trimesh.collision import CollisionManager
import numpy as np
import pandas as pd
import sys,os,re
from collections import Counter

from param_XML import Param_xml

verbose = True

class Contactinfo:
    def __init__(self, cell, area, area_eq, area_flat, contact_area, d_idx_l_idx_triangles, l_intersection_point):
        self.cell = cell
        self.area = area
        self.area_eq = area_eq
        self.area_flat = area_flat
        self.contact_area = contact_area
        self.d_idx_l_idx_triangles = d_idx_l_idx_triangles  #dict of lists : key = idx of triangle of cell 1 ; value = list of indices of contacting triangles of cell 2
        self.l_intersection_point = l_intersection_point
        self.poly_collision = extract_selection_vtp(poly_in = d_cellname_poly[cell], l_cell_idx = [*d_idx_l_idx_triangles])  #the poly that comes out trimesh collision manager (still contains gaps)
        self.s_points = None
        self.l_triangles_idx_hole = []
        self.l_triangles_popped = []


def read_parms():
    param_xml = Param_xml.get_param_xml(sys.argv,l_main_keys = ['body','contactinfo'],verbose=True) #param file must be passed as first argument
    
    input_folder= param_xml.get_value('input_folder',['paths'])
    output_folder = param_xml.get_value('output_folder',['paths'])

    return input_folder,output_folder
    
class DictOfDicts(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def write_summary_csv(dd_cell1_cell2_contactinfo, output_path):        
    l_columns = ["cell_1",
               "cell_2",
               "area",
               "area_eq",
               "area_flat",
               "contact_area",
               "nb_triangles",
               "ix_triangles"
            ]

    cell_1=[]
    cell_2=[]
    area=[]
    area_eq=[]
    area_flat=[]
    contact_area=[]
    nb_triangles=[]
    ix_triangles=[]
    
    for cell_1_i, d_cell2_contactinfo in dd_cell1_cell2_contactinfo.items():
        for cell_2_i, contactinfo in d_cell2_contactinfo.items():
            cell_1.append(cell_1_i)
            cell_2.append(cell_2_i)
            area.append(contactinfo.area)
            area_eq.append(contactinfo.area_eq)
            area_flat.append(contactinfo.area_flat)
            contact_area.append(contactinfo.contact_area)
            nb_triangles.append(len(contactinfo.d_idx_l_idx_triangles.keys()))
            ix_triangles.append([contactinfo.d_idx_l_idx_triangles.keys()])
            
    df = pd.DataFrame({l_columns[0]:cell_1,
                       l_columns[1]:cell_2,
                       l_columns[2]:area,
                       l_columns[3]:area_eq,
                       l_columns[4]:area_flat,
                       l_columns[5]:contact_area,
                       l_columns[6]:nb_triangles,
                       l_columns[7]:ix_triangles
                        })
                        
    df.to_csv(output_path,columns=l_columns,header=True,index=False)    

    return
    
    
def enrich_poly_cell_with_contacting_angles(cell):
    if verbose: print('enrich_poly_cell_with_contacting_angles ({0})'.format(cell))
    
    wdo_1 = dsa.WrapDataObject(d_cellname_poly[cell])
    a_max_angle = np.zeros_like(wdo_1.CellData['parentIndex'].__array__())

    d_idx_l_angle = {} #dict of list of angles; key=index of triangle of cell 1; value : list of angles with the triangles of cell 2
    
    for cell2_i, contactinfo_i in dd_cell1_cell2_contactinfo[cell].items():
        wdo_2 = dsa.WrapDataObject(d_cellname_poly[cell2_i])
        for idx_triangle_1_i, l_idx_triangles_2_i in contactinfo_i.d_idx_l_idx_triangles.items():
            a_normal_1 = wdo_1.CellData['normal'][idx_triangle_1_i].__array__()
            for idx_triangle_2_i in l_idx_triangles_2_i:
                a_normal_2 = wdo_2.CellData['normal'][idx_triangle_2_i].__array__()
                angle = np.degrees(np.arccos(np.clip(np.dot(a_normal_1,((a_normal_2 * -1))),-1.0,1.0)))  #1 normal is flipped to make it inward facing 
                # if angle > 90:angle = angle - (2 * (angle-90)) # we are only interested in angles between 0 and 90.
                d_idx_l_angle.setdefault(idx_triangle_1_i,[]).append(angle)  #normals are normalized, clipping as safety measure against rounding errors
                   
            a_max_angle[idx_triangle_1_i] = max(d_idx_l_angle[idx_triangle_1_i])

    
    d_cellname_poly[cell] = add_array(wdo_1.VTKObject,a_max_angle, "maxAngle", field_type="CELL")

    return 
    
def remove_shared_triangles_between_contactareas(cell):
    if verbose: print('remove_shared_triangles_between_contactareas ({0})'.format(cell))
    
    #build up datastructure d_idx_triangle_l_contactarea
    wdo_1 = dsa.WrapDataObject(d_cellname_poly[cell])
    d_idx_triangle_l_contactarea = {}
    for cell2_i, contactinfo_i in dd_cell1_cell2_contactinfo[cell].items():
        for idx_triangle_1_i in contactinfo_i.d_idx_l_idx_triangles.keys():
            d_idx_triangle_l_contactarea.setdefault(idx_triangle_1_i,[]).append(cell2_i) 
    
    #remove shared triangles from contactareas
    l_triangles_popped = []
    for idx_triangle_i, l_contactarea_i in d_idx_triangle_l_contactarea.items():
        if len(l_contactarea_i)>1:
            l_triangles_popped.append(idx_triangle_i)
            for contactarea_i in l_contactarea_i:
                dd_cell1_cell2_contactinfo[cell][contactarea_i].d_idx_l_idx_triangles.pop(idx_triangle_i)
                dd_cell1_cell2_contactinfo[cell][contactarea_i].l_triangles_popped.append(idx_triangle_i)
                
                if verbose:print('triangle {0} has been removed from contactarea {2} of cell {1} because it was shared with another contactarea from the same cell'.format(idx_triangle_i,cell,contactarea_i))
                
    if verbose:
        poly_popped = extract_selection_vtp(poly_in = d_cellname_poly[cell], l_cell_idx = l_triangles_popped)
        write_vtp_file(poly_popped, str(output_folder / "popped_{0}.vtp".format(cell) ) )
        
    return
    
def write_vtp_files():
    for cell_1_i, d_cell2_contactinfo in dd_cell1_cell2_contactinfo.items():
        for cell_2_i, contactinfo in d_cell2_contactinfo.items():
            poly_collision = extract_selection_vtp(poly_in = d_cellname_poly[cell_1_i], l_cell_idx = [*contactinfo.d_idx_l_idx_triangles] + contactinfo.l_triangles_popped) 
            write_vtp_file(poly_collision, str(output_folder / "{0}_contact_{1}.vtp".format(cell_1_i, cell_2_i) ) )  #with holes
            
            poly_holes_filled = extract_selection_vtp(poly_in = d_cellname_poly[cell_1_i], l_cell_idx = [*contactinfo.d_idx_l_idx_triangles] + contactinfo.l_triangles_idx_hole + contactinfo.l_triangles_popped)
            write_vtp_file(poly_holes_filled, str(output_folder / "{0}_contact_filled{1}.vtp".format(cell_1_i, cell_2_i) ) )   #without holes
            
        write_vtp_file(d_cellname_poly[cell_1_i], str(output_folder / "{0}_enriched.vtp".format(cell_1_i) ) )
    
    return


def get_assigned_area(vi1,vi2,vi3):
    """assign a triangle to the area for which it shares the most points (with a 2 point minimum) """
    l_contactarea = []
    for vi in [vi1,vi2,vi3]:l_contactarea += d_vi_arealabel.get(vi,[]) 
    counter_most_common = Counter(l_contactarea).most_common(1)    # {'most_common', counter)}
    if not counter_most_common:return []
    return counter_most_common [0][0] if counter_most_common[0][1]> 1 else []

    
#MAIN
l_selection_cells = ['cell_parent2']  #'all' will process all cells
l_selection_cells = ['all']  #'all' will process all cells

input_folder,output_folder = read_parms()
(output_folder / "STL_converted").mkdir(parents=True,exist_ok=True)


#Part1: find out contact pairs via Trimesh collision manager/
#trimesh cannot work with vtp directly, so convert to stl first  # print(trimesh.exchange.load.available_formats())
 
d_cellname_poly = {}
for ix,vtp_path_i in enumerate(sorted(input_folder.glob('*.vtp'))): 
    poly = read_vtp_file(vtp_path_i)
    d_cellname_poly[vtp_path_i.stem] = poly
    write_stl_file(poly, (output_folder / "STL_converted" / (vtp_path_i.stem + ".stl") ))
    
    
cm = CollisionManager()
d_cellname_trimeshStl = {}
for ix,stl_i in enumerate(sorted((output_folder / "STL_converted").glob('cell*.stl'))): 
    print('Adding stl ->  {0}'.format(stl_i.name))
    f = open(str(stl_i),'r')
    trimesh_stl = trimesh.load(f, file_type="stl")
    d_cellname_trimeshStl[stl_i.stem] = trimesh_stl
    cm.add_object(stl_i.stem, trimesh_stl)

is_collision, s_t_colliding_pairs = cm.in_collision_internal(return_names=True, return_data=False)

print("cm.in_collision_internal = {0} \nwith colliding pairs {1} ".format(is_collision, s_t_colliding_pairs))

if not is_collision:
    print("No collision between the cells, exitting")
    exit()

dd_cell1_cell2_contactinfo = DictOfDicts() #dict of dicts : key 1= cell1, key2 = cell2 , value = contactinfo object
l_vtp_columns = ['area','area_flat','area_eq','contact_area' ]

for t_colliding_pair in s_t_colliding_pairs:
    print('--Checking collision pair ->  {0}'.format(t_colliding_pair))
    cell_1,cell_2= t_colliding_pair
    
    if l_selection_cells[0] != 'all' and cell_1 not in l_selection_cells:continue
    
    cm_pair = CollisionManager()
    cm_pair.add_object( cell_1, d_cellname_trimeshStl[cell_1] )
    cm_pair.add_object( cell_2, d_cellname_trimeshStl[cell_2] )
    
    _, l_contactdata = cm_pair.in_collision_internal(return_names=False, return_data=True)
    
    l_intersection_point =  []
    d_idx_l_idx_triangles_1 = {}  #dict of lists : key = idx of triangle of cell 1 ; value = list of indices of contacting triangles of cell 2
    d_idx_l_idx_triangles_2 = {}
    for contactdata in l_contactdata:
        idx_triangle_1 = contactdata.index(cell_1) #stl face number matches the cellId in the vtp-file = triangle index
        idx_triangle_2 = contactdata.index(cell_2)
        l_intersection_point.append(contactdata.point)
        d_idx_l_idx_triangles_1.setdefault(idx_triangle_1,[]).append(idx_triangle_2)
        d_idx_l_idx_triangles_2.setdefault(idx_triangle_2,[]).append(idx_triangle_1)
        # print("{3}->{4} : face0 {0} intersects with face1 {1} at point {2}".format(contactdata.index(cell_1), contactdata.index(cell_2), contactdata.point,cell_1,cell_2))
    
    area_1, area_eq_1, area_flat_1, contact_area_1 = get_aggregate_data(d_cellname_poly[cell_1], np.array([*d_idx_l_idx_triangles_1]), l_aggregate_columns=l_vtp_columns, func="sum")
    dd_cell1_cell2_contactinfo[cell_1][cell_2] =  Contactinfo(cell_1, area_1, area_eq_1, area_flat_1, contact_area_1, d_idx_l_idx_triangles_1, l_intersection_point) 

    area_2, area_eq_2, area_flat_2, contact_area_2 = get_aggregate_data(d_cellname_poly[cell_2], np.array([*d_idx_l_idx_triangles_2]), l_aggregate_columns=l_vtp_columns, func="sum")
    dd_cell1_cell2_contactinfo[cell_2][cell_1] =  Contactinfo(cell_2, area_2, area_eq_2, area_flat_2, contact_area_2, d_idx_l_idx_triangles_2, l_intersection_point) 
    
#part 2 enrich cell VTP's with contacting angles
for cell_i in d_cellname_poly.keys(): enrich_poly_cell_with_contacting_angles(cell_i)

#part 3 remove shared triangles between contactareas
for cell_i in d_cellname_poly.keys(): remove_shared_triangles_between_contactareas(cell_i)

#part 4 : fill in holes of contactareas (all indices used in this step, refer to indices at cell-level)
for cell_i in dd_cell1_cell2_contactinfo.keys():
    if l_selection_cells[0] != 'all' and cell_i not in l_selection_cells:
        print('skipping cell {0}'.format(cell_i))
        continue
    
    d_vi_arealabel = {} #key  vertex index (=pointID);  value = the label of the contactarea (=name of contacting cell)
    #step 1 : get all contactpoints based on contactarea vtp 
    poly_totalcontactarea = extract_selection_vtp(poly_in = d_cellname_poly[cell_i], 
                                                query = "contact_area > 0" , 
                                                field_type="CELL", inverse_selection=False, verbose=False)
    wdo_total_contactarea = dsa.WrapDataObject(poly_totalcontactarea)
   
    a_contactIndex = np.copy(dsa.WrapDataObject(d_cellname_poly[cell_i]).CellData['parentIndex'].__array__()) #a_contactIndex = the array that will be used to enrich the cellVTP with contacting cellids, initialize with parentID
 
    s_idx_triangles_unassigned = set(wdo_total_contactarea.CellData['SelectionCellIds'].__array__())  #SelectionCellIds are the indices of 1 level up (cell level)
    
    write_vtp_file(poly_totalcontactarea, str(output_folder / "poly_totalcontactarea{0}.vtp".format(cell_i) ) )   #TEMP

    #step2 : gather pointsets of the different contactareas of the cell, and label them
    for cell_2_i, contactinfo in dd_cell1_cell2_contactinfo[cell_i].items():
        poly_collision = extract_selection_vtp(poly_in = d_cellname_poly[cell_i], l_cell_idx = [*contactinfo.d_idx_l_idx_triangles]) 
        wdo_contact_area =  dsa.WrapDataObject(contactinfo.poly_collision)
        contactinfo.s_points = set(list(wdo_contact_area.PointData['SelectionPointIds'].__array__()))
        
        # d_vi_arealabel.update({k:cell_2_i for k in contactinfo.s_points}) #  construct index to link a point index with contactarea (contactarea has cell2 as a label)
        for k in contactinfo.s_points: d_vi_arealabel.setdefault(k,[]).append(cell_2_i) #a point can be associated with multiple contactareas (in general)
        
        len_before = len(s_idx_triangles_unassigned)
        a_ix_triangles_contactarea = wdo_contact_area.CellData['SelectionCellIds'].__array__()
        s_idx_triangles_unassigned -= set(a_ix_triangles_contactarea) #these triangles are already assigned, take them out of the set
        a_contactIndex[a_ix_triangles_contactarea] = re.search(r'\d+', cell_2_i).group()  #extract number from cell
        if verbose:print('{0} triangles were already assigned and removed from the set'.format(len_before - len(s_idx_triangles_unassigned)))
        write_vtp_file(extract_selection_vtp(poly_in = d_cellname_poly[cell_i], l_cell_idx = list(s_idx_triangles_unassigned)), str(output_folder / "Missing triangles cell{0}.vtp".format(cell_i) ) )   #TEMP
        
    #step3 : assign the unassigned triangles (the ones in the holes, not part of the trimesh collision set).  Use the index we build up, adjacent triangles will be clustered to a contactarea at every iteration
    wdo_cell_i =  dsa.WrapDataObject(d_cellname_poly[cell_i])
    ic_update_happened = True
    while s_idx_triangles_unassigned and ic_update_happened:  #keep iterating until every triangle is assigned
        ic_update_happened = False
        s_iter = s_idx_triangles_unassigned.copy()
        for idx_triangle_i in s_iter:
            vi1,vi2,vi3 = wdo_cell_i.CellData['vertexIndices'][idx_triangle_i]
            assigned_area = get_assigned_area(vi1,vi2,vi3)
            if not assigned_area: 
                continue
            else:  #triangle gets assigned
                contactinfo = dd_cell1_cell2_contactinfo[cell_i][assigned_area]
                contactinfo.l_triangles_idx_hole.append(idx_triangle_i)
                contactinfo.s_points.update(set([vi1,vi2,vi3]))
                for k in [vi1,vi2,vi3]: d_vi_arealabel.setdefault(k,[]).append(assigned_area)
                s_idx_triangles_unassigned.discard(idx_triangle_i)
                a_contactIndex[idx_triangle_i] = re.search(r'\d+', assigned_area).group() 
              
                ic_update_happened = True
                if verbose:print ('Triangle {0} of cell {1} is clustered with the contactarea of cell {2}'.format( idx_triangle_i, cell_i, assigned_area))
    
    if s_idx_triangles_unassigned:
        print ('Warning : These triangles for cell {0} could not be assigned to an area : {1}'.format(cell_i, [i for i in s_idx_triangles_unassigned]))
    else:
        print ('All triangles for cell {0} could  be assigned to an area'.format(cell_i))
    
    d_cellname_poly[cell_i] = add_array(d_cellname_poly[cell_i],a_contactIndex, "contactIndex", field_type="CELL")  #enrich cell VTP with contactIndex
    
    
    
write_vtp_files()
write_summary_csv(dd_cell1_cell2_contactinfo, (output_folder / "contactinfo.csv"))

  