import trimesh
from trimesh.collision import CollisionManager
import math
# from trimesh.parent.Geometry import apply_scale
import sys, os, re
from collections import Counter
from VTK_utils import write_stl_file, read_vtp_file, read_stl_file
# import pyvista as pv
# import tetgen
from pymeshfix import _meshfix


class Cell:
	l_GT = []
	l_RES = []
	d_name_cell = {}

	def __init__(self, p_cell, ic_GT, cm, f_type='stl'):
		print('creating Cell {0}'.format(p_cell.stem))
		self.name = p_cell.stem
		self.ic_GT = ic_GT

		_meshfix.clean_from_file(str(p_cell), str(p_cell))  # CRUCIAL ! some stl are not closed, needed for volume and
		# correct boolean operation
		if f_type=="ply":
			self.ply = open(str(p_cell), 'r')
		else:
			self.stl = open(str(p_cell), 'r')
		# self.poly = read_stl_file(p_stl)
		# self.pyvista = pv.wrap(self.poly)
		# self.tet = tetgen.TetGen(self.pyvista)
		
		if f_type=="ply":
			self.trimesh = trimesh.load(self.ply, file_type="ply")
		else:
			self.trimesh = trimesh.load(self.stl, file_type="stl")

		self.d_cell_contactinfo = {}
		self.surface_area = self.trimesh.area
		self.extents= self.trimesh.extents #The length, width, and height of the axis aligned bounding box of the mesh.
		self.centroid = self.trimesh.centroid
		self.principal_inertia_components = self.trimesh.principal_inertia_components
		self.scale = self.trimesh.scale

		if self.trimesh.is_volume:
			self.volume = self.trimesh.volume
			self.centroid = self.trimesh.centroid
			self.sphericity = ((math.pi**(1/3))*((6.0e-6*self.trimesh.volume)**(2/3)))/(self.surface_area*1.e-6) # sphericity = ratio of surface area of a sphere that has the same volume as the object to the surface area of the particle itself.
			self.center_mass = self.trimesh.center_mass
		else:
			self.volume = None
			self.centroid = None
			self.sphericity = None
			self.center_mass = None
			print('warning : the cell {0} is not a proper volume according to trimesh'.format(self.name))

		self.cm = cm
		if cm:
			self.cm.add_object(self.name, self.trimesh)

		Cell.l_GT.append(self) if self.ic_GT else Cell.l_RES.append(self)

		Cell.d_name_cell[self.name] = self

	def get_max_overlapping_cell_with_volume(self,based_on_SEG_score=False):

		def get_SEG_score (cell1,cell2,contactinfo):

			if not (cell1.volume and cell2.volume):
				return 0
				
			union_volume = cell1.volume + cell2.volume - contactinfo.overlapvolume
			return contactinfo.overlapvolume / union_volume if union_volume > 0 else 0

		cell_max_overlap = None
		max_SEG_score = 0
		volume_max_overlap = 0
		for cell_i, contactinfo_i in self.d_cell_contactinfo.items():
			if based_on_SEG_score:
				SEG_i = get_SEG_score (self,cell_i,contactinfo_i)
				if SEG_i > max_SEG_score:
					volume_max_overlap = contactinfo_i.overlapvolume
					cell_max_overlap = cell_i
					max_SEG_score = SEG_i
			else:
				if contactinfo_i.overlapvolume > volume_max_overlap:
					volume_max_overlap = contactinfo_i.overlapvolume
					cell_max_overlap = cell_i
		return cell_max_overlap, volume_max_overlap


class Contactinfo:
	def __init__(self, overlapvolume):
		self.overlapvolume = overlapvolume