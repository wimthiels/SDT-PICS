'''
Created on 29 Jul 2018 
sync testupdata

DATAMODEL : 
-----------
there is a hierarchy of 4 classes that compose an actual cluster :
1) sphere : the basic building block.
2) cluster : a grouping of overlapping spheres.  has an ID (a counter), not a label
3) Cell : 1 or more clusters within the same stack that are assigned the same label. seed_with_cells_from_prev_stack
4) cell4D  : a cell that is tracked throughout time.  This is the actual cell in the biological sense

other classes
* Label : conceptually, the label is linked to a cell4D. But because this Label object is used throughout the build up of a cell, it is present on every level (except sphere)
* Featuretracker : this class is for reporting, it keeps track of data over the whole time dimension

Remark about RI: 
- Because of the bottom up construction of a cell, the data model does not have referential integrity DURING the construction, but it should be consistent AFTER processing. 
eg 
- you can have a Cluster that is not yet part of a Cell.
- you can have a Cell that is not yet part of a Cell4D

the chronology of 1 stack processing looks like this 
-> step 1 : create Spheres that group into Clusters
-> step 2 : label the Clusters
-> step 3 : use the labels to aggregate the Clusters into Cells
-> step 4 : apply cell division on the Cells
-> step 5 : update Cell4D objects

INDEXES
------- 
besides the obvious links via the object attributes some extra class indices are available for performant access : 

Cluster.d_ctr_cluster = {} #reverse index
Cluster.l_curr_stack = []  #only used in the labelling process of the clusters.  this will be ordered in decreasing size and can be processed sequentially

Cell4D.l_cell4D =[] #not deleted over time (only wise if cell4D is not connected to other objects)
Stack.curr_stack=[]  #only 2 stacks are available (to avoid memory problems)
Stack.prev_stack=[]

Label.d_ctr_label   #reverse index

within each stack you have access to cell-indices
self.cells = []
self.d_label_cell = {}  # queries with keys = stack and label)

    
        
@author: wimth
'''
from fileHandler import FileHandler
from helper_functions import get_3D_sphere_coordinates,plot3D_stack_with_overlay,examine

import numpy as np
import sys,glob
import pandas as pd
import os

DBG = False

                    
from scipy.ndimage.morphology import distance_transform_edt
#from scipy.spatial import distance
from skimage import img_as_ubyte
from math import floor,pi



class Stack:
    ctr = 0
    timestep_last_stack = 0
    
    curr_stack=[]
    prev_stack=[]
    l_RGBA_clusterview=[]  #only used when not saving the clusterview (not used when working with a temp save folder)
    
    tot_TP = 0
    tot_FN = 0 
    tot_FP = 0
    set_validation_names = set()
    d_name_lineage_ctr= {}

    def __init__(self, img,timestep,feat_tracker,name='no name stack'):
        Stack.ctr += 1
        self.ctr = Stack.ctr
        self.timestep = timestep
        self.name = name
        
        self.img = img   #format [Z,Y,X], masked array of the original image (=must be boolean image : True=membrane)
        self.exterior_mask = []     #a boolean np.array
        
        self.mask = np.zeros(self.img.shape,'bool')  #used as a working canvas for calculating overlap
        
        self.dist_transf = []           #the original DT
        self.max_DT = 0                #radius of the biggest sphere that can be placed
        self.dist_transf_carved = []   #the DT that will be carved out by spheres
        
        self.clusterview = np.zeros(img.shape,'int16') #will be filled with the counters of the clusters (NOT the cluster-id's)
        self.seedview = np.array([]) #this is a snapshot of the clusterview after the seeding (before extra spheres get added)
        # self.minseedview = np.zeros(img.shape,'int16')  #this is a seed clusterview where every seed has a minimum allowed size(min-seed-radius)
        self.labelview = np.zeros(img.shape,'int16') 
    
        self.feat_tracker = feat_tracker
        
        self.cells = []
        self.d_label_cell = {}  #Index (queries with keys = stack and label)
        
        self.l_spheres_filtered = [] #these spheres are logged to show on the clusterview (will not be linked to a cluster)
        
        self.nb_FN = 0
        self.nb_TP = 0
        self.nb_FP = 0
        
        self.z_min_ix_signal = 0
        self.z_max_ix_signal = Param.Z_DIM - 1
        
        
        self.move_to_next_stack()
        
    @staticmethod
    def init_class_variables():

        Stack.ctr = 0
        Stack.timestep_last_stack = 0
    
        Stack.curr_stack=[]
        Stack.prev_stack=[]
        Stack.l_RGBA_clusterview=[]  #only used when not saving the clusterview (not used when working with a temp save folder)
        
        Stack.tot_TP = 0
        Stack.tot_FN = 0 
        Stack.tot_FP = 0
        Stack.set_validation_names = set()
        
        Cell4D.init_class_variables()
        Cluster.init_class_variables()
        Sphere.init_class_variables()
        Label.init_class_variables()
        
        return
    
        
        
    def get_clusters_for_label(self,label):
        if label.is_real():
            cell = self.d_label_cell.get(label)
            if cell: 
                return cell.get_clusters()
        else:
            l_clusters = []
            for cluster_i in Cluster.d_ctr_cluster.values():
                if not cluster_i.stack == self:
                    continue
                if cluster_i.label == label:
                    l_clusters.append(cluster_i)
            return l_clusters
                    
        return None

    
    def get_clusters(self):
        l_clusters = []
        for cell_i in self.cells:
            l_clusters += cell_i.get_clusters()
        
        return l_clusters
 

    
    def move_to_next_stack(self):
        '''
        do the memory housekeeping that results from shifting stacks (only 2 stacks are kept in memory)
        '''
        #clean up class memory zones (not cell4D)
        if self.prev_stack:
            l_cluster_2Bdeleted = []
            for cluster_i in Cluster.d_ctr_cluster.values():
                if cluster_i.stack == Stack.prev_stack:
                    l_cluster_2Bdeleted.append(cluster_i)           
            for cluster_i in l_cluster_2Bdeleted:       
                    Cluster.d_ctr_cluster.pop(cluster_i.ctr,None)
            
        Cluster.l_curr_stack = []
        
        
        
        #shift timestep
        Stack.prev_stack = Stack.curr_stack
        Stack.curr_stack = self
        
        
        return
    
    
    def add_exterior_mask(self, img_exterior,mask=True):
        '''
        makes a boolean stack where the interior is set to True, and the exterior to False
        
        :param img_exterior: input image that marks the exterior
        :param mask: if the input is a 'mask', it means it marks the thing to hide (=exterior in this case) as False
        '''
        if img_exterior.dtype != 'bool':
            exterior_mask_bool = np.zeros(img_exterior.shape,'bool')
            if mask:
                exterior_mask_bool[img_exterior>0] = True
            else:
                exterior_mask_bool[img_exterior==0] = True
            
        self.exterior_mask = img_exterior
        
        self.clusterview[self.exterior_mask == 0]  = Cluster.CTR_MASK  
        
        return
        
    
            
    def add_distance_transform(self,verbose=True):
        '''
        adds a DT for a stack. but now apply the mask afterwards and also do not use slabs
        remark : the DT measures the euclidean distance to a zero/false value.   So we need to invert the membrane image.  We take into account the actual size of the voxels
        '''
        if verbose:print("###make distance transform for stack{0}".format(self.timestep))
        
        def put_y_membranes_bars_on_blank_slices():
            
            def draw_y_membranes_on_slice(z_i):
                for y_i in range(0,Param.Y_DIM,Param.MIN_SPHERE_RADIUS):
                        self.img[z_i,y_i,...]=1
                
                return 1
         
         
            z_dim = self.img.shape[0]
                
            for z_i in range(z_dim):   #top to bottom
                if np.sum(self.img[z_i,...])==0:
                    draw_y_membranes_on_slice(z_i)
                else:
                    self.z_min_ix_signal = z_i -1
                    break        
                
            for z_i in range(z_dim-1,-1,-1):   #bottom to top
                if np.sum(self.img[z_i,...])==0:
                    draw_y_membranes_on_slice(z_i)
                else:
                    self.z_max_ix_signal = z_i - 1 
                    break 
            
            return
        
        
        def use_blank_img_planes_as_exterior_mask():
            z_dim = self.img.shape[0]
            for z_i in range(z_dim):   #top to bottom
                if np.sum(self.img[z_i,...])==0:
                    self.dist_transf[z_i,...] = 0
                else:
                    break
            for z_i in range(z_dim-1,-1,-1): #bottom to top
                if np.sum(self.img[z_i,...])==0:
                    self.dist_transf[z_i,...] = 0
                else:
                    break
                
    #             for y_i in range(Param.Y_DIM):   #front to back
    #                 if np.sum(self.img[:,y_i,...])==0:
    #                     self.dist_transf[:,y_i,...] = 0
    #                 else:
    #                     break
    #             for y_i in range(Param.Y_DIM-1,-1,-1): #back to front
    #                 if np.sum(self.img[:,y_i,...])==0:
    #                     self.dist_transf[:,y_i,...] = 0
    #                 else:
    #                     break
    #                 
    #                 
    #             for x_i in range(Param.X_DIM):   #left to right
    #                 if np.sum(self.img[:,:,x_i])==0:
    #                     self.dist_transf[:,:,x_i] = 0
    #                 else:
    #                     break
    #             for x_i in range(Param.X_DIM-1,-1,-1): #right to left
    #                 if np.sum(self.img[:,:,x_i])==0:
    #                     self.dist_transf[:,:,x_i] = 0
    #                 else:
    #                     break
        
            return


        if Param.MASK_AS_MEMBRANE:
            if len(self.exterior_mask):
                self.img[self.exterior_mask==0] = True
        
        #sandwich the image with 2 empty slices (these will be filled with y-membranes)
        self.img=np.insert(self.img ,0,False,axis=0)   #ztop
        self.img=np.insert(self.img ,self.img.shape[0],False,axis=0) #zbottom


        #         if not len(self.exterior_mask):
        put_y_membranes_bars_on_blank_slices()
        
        self.dist_transf = distance_transform_edt(np.invert(self.img), sampling = [Param.XY_VS_Z,1,1], 
                                       return_distances=True, return_indices=False) [1:-1,:,:]
                                       
        self.img = self.img[1:-1,:,:]
                
        
        #apply the exterior mask
        if len(self.exterior_mask):
            self.dist_transf[self.exterior_mask==0] = 0
            self.img[self.exterior_mask==0] = False
        else:
            use_blank_img_planes_as_exterior_mask()
        
        self.dist_transf_carved = self.dist_transf.copy()
        self.max_DT = round(np.max(self.dist_transf),0)  
        return   

        
    
    def get_overlap_maxDT_threshold(self,sphere,overlapping_cluster=None):  
        '''
        An H-max threshold will be calculated.  This H-max threshold indicates a radius length (expressed in nb of pixels).  Think of shrinking the radius of the sphere in order for it to be able
        to pass to the overlapping cluster.  If you need to shrink the radius more than the H-max threshold, the sphere will not be joined with the overlapping cluster
        The standard H-max (=h_max_relative_max) is set at 1/3 of the size of the sphere. (if Param_XML.FRAGMENTATION_LEVEL=2/3)

        extension : if overlapping cluster is given this can affect the threshold by incorporating zspan
           
        :param min_radius: the radius of the current sphere

        '''

        def get_zspan(sphere1,sphere2):
            '''
            z-span is the z-distance between top and bottom z of 2 spheres
            '''
            if sphere1.centre[0] > sphere2.centre[0]:
                sphere_min = sphere2
                sphere_max = sphere1
            else:
                sphere_min = sphere1
                sphere_max = sphere2

            return (sphere_max.centre[0] + (sphere_max.radius/Param.XY_VS_Z)) - (sphere_min.centre[0] - (sphere_min.radius/Param.XY_VS_Z)) 
 

        if overlapping_cluster:
            z_span = get_zspan(sphere,overlapping_cluster.spheres[0])
            if z_span > Param.z_span_max:
                #  print("zspan criterium satisfied zspan = {2} with zmin = {0}, zmax = {1}".format(sphere.centre[0] ,overlapping_cluster.spheres[0].centre[0],z_span))
                if sphere.radius / self.max_DT > Param.radius_zspan: 
                    print("FRAGMENTATION_LEVEL raised to 1 : radius sphere = {0}, radius cluster = {1} and z_span = {2}".format(sphere.radius,overlapping_cluster.spheres[0].radius,z_span))
                    return sphere.radius * 1 #FRAGMENTATION_LEVEL = 1

        return  sphere.radius * Param.FRAGMENTATION_LEVEL
       
    
    def seed_with_cells_from_prev_stack(self,verbose=True):
        if verbose:print("###seeding with cells from the previous stack, stack=stack{0}".format(self.timestep)) 
        
        for cell_i in Stack.prev_stack.cells:
            if cell_i.label in Label.L_INDICATOR:
                continue
            self.add_sphere_on_max_dt(seeding_cell=cell_i, verbose=verbose)
        
        self.seedview = self.clusterview.copy()
        
        return

    def seed_with_validation_file(self, df_nuclei, validation_start_time=1, verbose=True):
        if verbose: print("###seeding with cells from validation file, stack=stack{0}".format(self.timestep))

        for _, nucleus in df_nuclei[df_nuclei.time == self.timestep + validation_start_time-1].iterrows():
            # Beware : Lineager is 1-indexed ! so x,y and z should always be decremented for indexing (cfr seeding and validating)
            self.add_sphere_on_max_dt(
                GT_seedZYX=[int(floor(nucleus.z) - 1), int(floor(nucleus.y) - 1), int(floor(nucleus.x) - 1)],
                GT_radius=nucleus.radius,
                verbose=verbose)

        self.seedview = self.clusterview.copy()
    
    def reset_seeds(self,func="all",verbose=True):
        if func=="no_reset":return  
        if not self.seedview.any(): return
        if verbose:print("###resetting seeds, stack=stack{0}".format(self.timestep))
        
        
        if func=='all':
            self.clusterview = np.where(self.seedview != 0,self.seedview,self.clusterview)
        if func=='wipeout':  #only restore the wiped out seeds
            a_seed_labels    = np.unique(self.seedview)
            a_seed_wipedout  = a_seed_labels[np.in1d(a_seed_labels, self.clusterview, invert=True)] #np.in1d returns boolean array with missing values = True
            for cluster_ctr_i in a_seed_wipedout:
                print('restoring seed cluster {0} from wipeout'.format(cluster_ctr_i))
                self.clusterview = np.where(self.seedview==cluster_ctr_i,cluster_ctr_i,self.clusterview)
        
        return
    
    
    def put_spheres_in_distance_transform(self,verbose=True):
        if verbose:print("");print("###put spheres on maximum of the DT of stack{0}".format(self.timestep))
        
        ctr= 0  
        while True:
            if ctr>=Param.MAX_NB_SPHERES :break
            progress_bar=True if ctr%5==0 else False
            sphere_added = self.add_sphere_on_max_dt(progress_bar=progress_bar)
            if not sphere_added: break
            current_radius = sphere_added.radius
            if current_radius <Param.MIN_SPHERE_RADIUS:break
            ctr += 1
            
        
    #         if len(self.exterior_mask):
    #             self.clusterview[self.exterior_mask==0] = Cluster.CTR_MASK  #clusterview gets trimmed off, the spheres remain intact
        
        print("")  
        return    
    
    def get_nb_overlap_pix_with_mask(self,sphere):
        
        if len(self.exterior_mask)==0:
            return 0
        
        bincount = np.bincount(self.exterior_mask[tuple(sphere.l_pixels)].flatten())
        
        return bincount[0]
    
    def get_cluster_ctr_of_max_overlap(self,sphere):
        '''
        the cluster counter of the sphere that maximally overlaps with the current sphere is returned
        '''
        bincount = np.bincount(self.clusterview[tuple(sphere.l_pixels)].flatten())  
        bincount[0] = 0 #not interested in overlap with zero region, so set zero count equal to zero
        
        nb_overlap_pix_mask=0
        
        if len(bincount)>Cluster.CTR_MASK:
            nb_overlap_pix_mask = bincount[Cluster.CTR_MASK]
            bincount[Cluster.CTR_MASK] = 0 #not interested in overlap with mask, so set zero count equal to zero
        
        return np.argmax(bincount),nb_overlap_pix_mask,bincount

    def add_sphere_on_max_dt(self, seeding_cell=None,
                             GT_seedZYX=None,
                             GT_radius=None,
                             verbose=True,
                             progress_bar=True):
        '''
        -a sphere is placed inside the DT at the current max position 
        -for every sphere it is determined if it is the seed of a new cluster, or that is should be added to an existing cluster (based on an overlap threshold)
        -the DT will be carved out (sphere set to zero), and the cluster view will be updated with the cluster counter of the cluster to which the sphere is assigned
        -At this point, cluster are NOT yet labelled (=get Label.UNASSIGNED_CTR). Except in 2 cases  
            -  Label.EXTERIOR_CTR if the sphere touches the exterior frame
            -  during the seeding process : the seeding cluster comes from the previous stack.   the max finding will be restricted to pixels inside this cluster and the label will already be inherited

        Args:
            GT_radius (float): radius of GT seed sphere in micron.  seed sphere will not be smaller than this value
        '''
        if progress_bar:print('|',end="",flush=True)
        
        
        def create_sphere_on_max_DT():
            '''
            The biggest possible sphere possible is created and placed in the (carved out) DT.  
            -When seeding, only the seeding volume is considered
            -if the maximum sphere is smaller then the minimum cluster radius, no sphere is created
            '''
            if seeding_cell: #only consider the seeding volume 
                self.mask.fill(False);
                self.mask[tuple(seeding_cell.get_seeding_volume())] = True;
                a_dist_transf= np.where(self.mask, self.dist_transf_carved,0)
            elif GT_seedZYX: #for a GT seed we want a sphere on that exact location, only if DT is zero we settle for a tiny sphere as seeding volume
                self.mask.fill(False)
                if self.dist_transf_carved[tuple(GT_seedZYX)] == 0:
                    print('GT logic : GT seeds lands on DT value = 0 , so a tiny seeding volume will be taken (small sphere of min-seed-radius',GT_seedZYX)
                    t_tiny_seeding_volume = get_3D_sphere_coordinates(a_coordZYX=GT_seedZYX,
                                                  radius=Param.MIN_SEED_RADIUS,
                                                  limit_shape=(Param.Z_DIM,Param.Y_DIM,Param.X_DIM),
                                                  xy_vs_z=Param.XY_VS_Z);
                                                  
                    self.mask[t_tiny_seeding_volume] = True;
                else:
                    self.mask[tuple(GT_seedZYX)] = True; #seeding volume of 1 pixel
                a_dist_transf= np.where(self.mask, self.dist_transf_carved,0)
            else:
                a_dist_transf = self.dist_transf_carved
                
            #get coordinates of the maximum value of the carved DT
            #max_dt = round(np.max(a_dist_transf),0)   
            max_dt = round(np.max(a_dist_transf),2)
            if max_dt < Param.MIN_SPHERE_RADIUS and not GT_seedZYX: 
                print('-',end="",flush=True)
                if seeding_cell:
                    self.feat_tracker.add_event(19,[seeding_cell.label,max_dt])
                return None

            if GT_seedZYX:
                max_dt_zyx = GT_seedZYX
                # max_dt = min(max_dt,Param.MIN_SEED_RADIUS)  ##GT seeding, limit movement because GT is correct position!
                max_dt = max(GT_radius,
                             max_dt,
                             Param.MIN_SEED_RADIUS)  # GT seeding, too small spheres will get engulfed, so seed larger
                max_dt = 25 # unequal seeds can cause problems
            elif seeding_cell:
                ##max_dt = min(max_dt,Param.MIN_SEED_RADIUS)  ##normal seeding, to be assessed  limit movement good idea ?? : No worsen results...
                max_dt_zyx = list(np.unravel_index(np.argmax(a_dist_transf), a_dist_transf.shape))
            else:
                max_dt_zyx = list(np.unravel_index(np.argmax(a_dist_transf), a_dist_transf.shape))
                
            
        
            #create sphere on that max value
            sphere = Sphere(centre_ZYX=max_dt_zyx,
                         radius=max_dt,
                         stack=Stack.curr_stack)
             
            #shrink radius if overlapping with min seed radius of previous cells, prevent wipe outs
            # if GT_seedZYX or seeding_cell:
            #     if np.sum(self.minseedview[tuple(sphere.l_pixels)]) > 0:
            #         print('information : a seed still overlaps with another sphere (minseedview) . radius seed last added =',max_dt)

                    # if GT_seedZYX:
                    #     print('a seeding sphere is shrunken because it comes to close to a previous seed (minseedview check)',GT_seedZYX ,seeding_cell,"old radius=",max_dt," new_radius=",Param.MIN_SEED_RADIUS)
                    #     sphere = Sphere(centre_ZYX=max_dt_zyx, radius=Param.MIN_SEED_RADIUS, stack=Stack.curr_stack)
             
                         
            
            return sphere
        
        
        
        
        def create_new_cluster_with_sphere(sphere,max_cluster_dcopy_overlap=None):
            
            cluster = Cluster(sphere,self,seeding_cell=seeding_cell,GT_seedZYX=GT_seedZYX)  #use sphere to create a new cluster  

            if max_cluster_dcopy_overlap:
                sphere.d_clusterctr_overlap[max_cluster_ctr] = max_cluster_dcopy_overlap  
                cluster.update_overlap_pixels(max_cluster_ctr,max_cluster_dcopy_overlap)   
            return cluster
            
        def add_sphere_to_existing_cluster():
            cluster = Cluster.d_ctr_cluster[max_cluster_ctr]
            cluster.add_sphere(sphere)
            sphere.overlap_max  = bincount[max_cluster_ctr]
            return cluster
        
        def get_overlap_info_of_current_sphere_with_cluster(cluster_ctr):
            '''
            the masked_clusterview only looks at a certain region of the clusterview (so a_masked_clusterview must be initialized properly eg. by delineating the current sphere.)
            From this view the overlap info is gathered wrt the other clusters
            The result is given back as a filled in copybook
            :param cluster_ctr:
            '''
            cluster = Cluster.d_ctr_cluster[cluster_ctr]
            Cluster.dcopy_overlap['overlap_pix'] = np.where(a_masked_clusterview==cluster_ctr)  #the coordinates of the overlapping pixels of a cluster with the current sphere
            Cluster.dcopy_overlap['overlap_max_DT'] = np.max(self.dist_transf[cluster.dcopy_overlap['overlap_pix']])  #the maximum DT value in the overlap zone
          
            return Cluster.dcopy_overlap.copy()
        
        def get_overlap_with_other_clusters():
            nonlocal max_cluster_ctr
            bincount[max_cluster_ctr] = 0 #the highest overlap is with the own cluster, only interested in the other clusters
            max_cluster_ctr = np.argmax(bincount)
      
            while max_cluster_ctr != 0:
                cluster.update_overlap_pixels(max_cluster_ctr,get_overlap_info_of_current_sphere_with_cluster(max_cluster_ctr))
                bincount[max_cluster_ctr] = 0
                max_cluster_ctr = np.argmax(bincount)
            
            return 
        
        def filter_sphere():
                                
            if GT_seedZYX:
                return False #a sphere put there by GT seeding should never be filtered
                
            if Param.FILTER_SPHERE_SEED_TOP_OR_BOTTOM:
                if sphere.centre[0] in [self.z_min_ix_signal,self.z_max_ix_signal]:
                    sphere.reason_cd_filtering = '4'
                    return True # a sphere with its centre on top or bottom is highly unlikely
                
                
            if max_cluster_ctr!=0: #if there is overlap, there will never be filtering
                return False
            
            if Param.FILTER_SPHERE_TOUCH_EXTERIOR_FRAME:
                if sphere.flag_exterior_frame:
                    sphere.reason_cd_filtering = '1'
                    if DBG:print('reason_cd_filtering =1', sphere.ctr, sphere.centre[0])
                    return True  #a sphere that overlaps with the frame and does not overlap with any other cluster is filtered
            
            
            if Param.FILTER_SPHERE_TOO_MUCH_OVERLAP_WITH_MASK:
                if nb_overlap_pix_mask/len(sphere.l_pixels[0]) > Param.FILTER_SPHERE_OVERLAP_MASK_PC:
                    sphere.reason_cd_filtering = '2'
                    if DBG:print('reason_cd_filtering =2', sphere.ctr)
                    return True  #a sphere that overlaps heavily with the mask and does not overlap with any other cluster is filtered
            
            if Param.FILTER_SPHERE_SMALL_RADIUS:
                if Stack.prev_stack:
                    if sphere.radius <= Param.MIN_SEED_RADIUS:
                        sphere.reason_cd_filtering = '3'
                        if DBG:print('reason_cd_filtering =3', sphere.ctr)
                        return True  # small cells with no overlap should be filtered
                else:
                    if sphere.radius <= Param.MIN_CELL_RADIUS:
                        sphere.reason_cd_filtering = '3'
                        return True  # small cells with no overlap should be filtered (more stringent first stack)

            
            return False
        
        def threshold_cluster_creation():

            #test the conditions
            if GT_seedZYX:#if the ground truth is provided, then SpheresDT must obey
                return True
                
            if max_cluster_ctr==0:
                return True #a sphere with no overlap to an existing cluster should give rise to a new cluster, no other option
            
            if Stack.prev_stack:
                if sphere.radius <= Param.MIN_SEED_RADIUS: 
                    return False # the seeding sphere should be big enough to start a new cluster
            else:
                if sphere.radius <= Param.MIN_CELL_RADIUS: 
                    return False # the seeding sphere should be big enough to start a new cluster.  in the first stack this will lead to a cell creation
                
            if seeding_cell:
                return True #the conditions below do not apply anymore for a seeding cluster
       
            if (nb_overlap_pix_mask/len(sphere.l_pixels[0]) > Param.CLUSTER_THRESHOLD_OVERLAP_MASK_PC or sphere.flag_exterior_frame) and max_cluster_ctr!=0:
                return False #a sphere that overlaps heavily with the exterior but does have some overlap with a cluster should never start a new cluster
            

            if  max_cluster_dcopy_overlap['overlap_max_DT'] < self.get_overlap_maxDT_threshold(sphere,Cluster.d_ctr_cluster.get(max_cluster_ctr)):
                return True #a sphere that does not overlap well enough with another cluster (based on DTmax of overlapzone), should start a cluster of its own

            if seeding_cell:
                return True
            
            
            return False
        
        #create sphere
        sphere = create_sphere_on_max_DT()
        if not sphere: return
        
        #gather sphere info (overlap info with cluster of highest overlap)                                                   
        self.mask.fill(False);self.mask[tuple(sphere.l_pixels)] = True;a_masked_clusterview = np.where(self.mask, self.clusterview,0)
        max_cluster_ctr,nb_overlap_pix_mask, bincount = self.get_cluster_ctr_of_max_overlap(sphere)
        max_cluster_dcopy_overlap=None
        if max_cluster_ctr!=0:
            max_cluster_dcopy_overlap = get_overlap_info_of_current_sphere_with_cluster(max_cluster_ctr)
        
        #filter sphere
        if filter_sphere():
            self.dist_transf_carved[tuple(sphere.l_pixels)] = Cluster.CTR_BACKGROUND
            self.l_spheres_filtered.append(sphere)
            return sphere
        
        #create cluster (threshold)
        if threshold_cluster_creation():
            cluster = create_new_cluster_with_sphere(sphere,max_cluster_dcopy_overlap)
        else:
            cluster = add_sphere_to_existing_cluster()

        if cluster:    
            get_overlap_with_other_clusters() 

        #update stack views   
        if cluster: 
            self.clusterview[tuple(sphere.l_pixels)] = cluster.ctr
            # self.minseedview[get_3D_sphere_coordinates(a_coordZYX=sphere.centre,
            #                                            radius=min(sphere.radius,Param.MIN_SEED_RADIUS),
            #                                            limit_shape=(Param.Z_DIM,Param.Y_DIM,Param.X_DIM),
            #                                            xy_vs_z=Param.XY_VS_Z)] = cluster.ctr
        else:
            self.clusterview[tuple(sphere.l_pixels)] = Cluster.CTR_BACKGROUND
        self.dist_transf_carved[tuple(sphere.l_pixels)] = 0
            
        return sphere
    
    
    def compose_label_view(self):
        '''
        labelview is created (based on clusterview)
        '''

        for label_i in self.d_label_cell.keys(): 
            l_cluster_ctr = [cluster_i.ctr for cluster_i in self.get_clusters_for_label(label_i)]
            self.labelview[np.isin(self.clusterview,l_cluster_ctr)] = label_i.ctr  
        return
    
    
        
    def label_clusters_of_first_stack(self,ic_GT_corr,verbose=True):
        if verbose:print("###label clusters for stack {0} based on current stack information".format(self.timestep))
        
        def threshold_cluster_fusion_first_stack():
            
            if ic_GT_corr:  #in case GT is delivered via a validation file
                if overlap_cluster.GT_seeded:
                    return False  #a GT cluster should never be merged with another cluster (removes FN)
                else:
                    return True   #a non_GT cluster should always be merged (removes FP)
            
            if overlap_cluster.get_radius() < Param.MIN_CELL_RADIUS:
                return True # a cluster to small to be a seed should always clump to a bigger cluster

            if d_overlap['overlap_max_DT'] > self.get_overlap_maxDT_threshold(overlap_cluster.spheres[0],cluster_i):

                for overlap_ctr_reverse, d_overlap_reverse in overlap_cluster.d_clusterctr_overlap.items():
                    if d_overlap_reverse['overlap_max_DT'] > d_overlap['overlap_max_DT']:
                        if Cluster.d_ctr_cluster[overlap_ctr_reverse].get_seed_radius() > overlap_cluster.get_seed_radius():
                            print ('the big cluster {0} overlaps sufficiently with the smaller {1}, BUT the smaller has an even bigger overlap with another cluster {2} which is bigger than itself'.format(cluster_i.ctr,overlap_cluster_ctr,overlap_ctr_reverse))
                            return False # the big cluster overlaps sufficiently with the smaller, BUT the smaller has an even bigger overlap with another cluster which is bigger than itself

                return True # if the overlap between the two clusters is big enough, they should clump together
            
            
            
            
            return False
            
        
        l_processed_clusters = []

        #will be ordered from big to small (cluster.get_seed_radius()). 
        #A big cluster will always search for smaller cluster to fuse with, never the other way around (by design). 
        #once a cluster is fused it will stay fused (no stealing)
        for cluster_i in Cluster.l_curr_stack:  
            if verbose:print('first stack handling cluster_i {0} (radius={1})'.format(cluster_i.ctr,cluster_i.get_seed_radius()))
            
            
            if ic_GT_corr and not cluster_i.GT_seeded: 
                continue #if GT is delivered than a cluster that is not GT-seeded should get a new label
            
            if Param.THRESHOLD_CELL_OVERLAP_TOO_MUCH_WITH_MASK and not ic_GT_corr:
                if cluster_i.overlaps_too_much_with_mask(self): 
                    cluster_i.label = Label.UNASSIGNED
                    l_processed_clusters.append(cluster_i)
                    continue
            
            if cluster_i in l_processed_clusters: 
                print('---cluster_i already processed, no new label will be drawn. but still checking overlapping cluster for label passing')
            else:
                cluster_i.label = Label(timestep=self.timestep,label_parent_ctr=Label.PARENT_ZEUS_CTR)  #draw new label with no parent
        #             if DBG:print('new label is created {0}'.format(cluster_i.get_color()))
            
            for overlap_cluster_ctr, d_overlap in cluster_i.d_clusterctr_overlap.items():  #check all the overlapping clusters
                overlap_cluster = Cluster.d_ctr_cluster[overlap_cluster_ctr]
                if verbose:print('---first stack handling overlap_cluster_ctr {0} (radius={1})'.format(overlap_cluster_ctr,overlap_cluster.get_seed_radius()))
                if overlap_cluster in l_processed_clusters: 
                    if verbose:print('--- overlap_cluster already processed, move on the next overlap cluster')
                    continue
                if threshold_cluster_fusion_first_stack():
                    if verbose:print('---fusing, move on to next cluster ...')
                    overlap_cluster.label = cluster_i.label
                    l_processed_clusters.append(overlap_cluster)
            #<-- 1 overlapping cluster for 1 cluster handled 
            
            l_processed_clusters.append(cluster_i)
            self.feat_tracker.add_event(1,[cluster_i])
        #<--1 cluster handled  
        #if verbose:print(['{0}_{1}'.format(i.ctr,i.label.color) for i in Cluster.l_curr_stack])
        return
    
    def label_clusters_via_prev_stack_information(self,verbose=True):
        if verbose:print("###label clusters for stack {0} based on previous stack information".format(self.timestep),end='')
        
        def assign_max_overlapping_label(cluster):
            bincount = np.bincount(Stack.prev_stack.labelview[cluster.get_l_pixels()].flatten())
            bincount[0] = 0  #get rid of overlap with background
            label_ctr = np.argmax(bincount)
            if label_ctr ==0:
                cluster.set_label(label=Label.UNASSIGNED,trace='no_match') #no overlap is not normal, we can filter these out later if we want
            else:
                label=Label.d_ctr_label[label_ctr]
                cluster.set_label(label=label,trace='labelview')
                
            return
        
        Stack.prev_stack.compose_label_view()  #compose labelview of the previous stack
        
        for cluster in Cluster.l_curr_stack:   
            if cluster.seeding_cell_label and not cluster.GT_seeded:
                cluster.set_label(cluster.seeding_cell_label, trace='seeded') #seeding label always takes precedence
            else:
                assign_max_overlapping_label(cluster) #go over every cluster of current stack and inherit label based on maximum overlap (of cluster)
        
            
        
        #if verbose:print(['{0}_{1}'.format(i.ctr,i.label.color) for i in Cluster.l_curr_stack])
        return
    

    def filter_clusters(self,ic_GT_corr=False,verbose=True):
        if verbose:print(">>>Filter clusters<<<",end="")
        

        
        ctr_filter = 0
        for cluster_i in Cluster.l_curr_stack:
            # cluster that failed to get a label from the previous have moved too much and should be removed
            if cluster_i.label == Label.UNASSIGNED:
                ctr_filter +=1
                cluster_i.set_label(Label.FILTERED,trace='unassigned')
                
            #never filter a cluster seeded by validation
            if cluster_i.GT_seeded:
                continue  
                
            # a cluster with its seed on top or bottom is highly unlikely (given the shell for one)
            if Param.FILTER_CLUSTER_SEED_TOP_OR_BOTTOM:
                if cluster_i.centre[0] in [self.z_min_ix_signal,self.z_max_ix_signal]:
                    ctr_filter +=1
                    cluster_i.set_label(Label.FILTERED,trace='z top or bottom')
                

            # a cluster that touches the exterior frame should be filtered 
            if Param.FILTER_CLUSTER_TOUCH_EXTERIOR_FRAME:
                if cluster_i.flag_exterior_frame:
                    ctr_filter +=1
                    cluster_i.set_label(Label.FILTERED,trace='touches exterior frame')
                
        if verbose:print("->",ctr_filter,"filtered")
        return
        
    def create_cells_based_on_cluster_labels(self,verbose=True): 
        if verbose:print("###create cells for stack {0} ".format(self.timestep), end='')
        
        def sort_clusters_on_label_seeded_and_seed_radius():
            l_t = []
            for cluster_i in Cluster.l_curr_stack:
                if cluster_i.seeding_cell_label or cluster_i.GT_seeded:
                    ic_seeded = 1
                else:
                    ic_seeded = 0
                l_t.append((cluster_i,cluster_i.label.ctr,-ic_seeded,-cluster_i.get_seed_radius())) #(cluster(key),ctr(ASC),ic_seeded(DESC),volume(DESC))
            l_t = sorted(l_t, key=lambda l_t: (l_t[1],l_t[2],l_t[3]))
            return [t[0] for t in l_t]
        
        nb_cells_created = 0
        l_clusters_sorted_on_label_seeded_and_volume = sort_clusters_on_label_seeded_and_seed_radius()
        label_old = ''
        cell_current = ''
        
        for cluster_i in l_clusters_sorted_on_label_seeded_and_volume:
            label_new = cluster_i.label.ctr
            if label_old==label_new:
                cell_current.add_cluster(cluster_i)
            else:
                cell_current = Cell(label=cluster_i.label,
                                    cluster_seed=cluster_i,
                                    stack=self)
                nb_cells_created +=1

                    
            label_old=label_new
            
        if verbose:print("->",nb_cells_created,"created",end='')
        #if verbose:print(['{0} from {1}'.format(i.label.color,i.cluster_seed.ctr) for i in self.cells])
        return
        
    #     def create_cells(self, verbose=True):
    #         if verbose:print("###create cells for stack {0} ".format(self.timestep), end='')
    #         
    #         #mark seeding clusters
    #         
    #         
    #         #create a cell for each seeding cluster 
    #         nb_cells_created = 0
    #         for cluster_i in Cluster.l_curr_stack: #create cells using the seed cluster
    #             if cluster_i.label in Label.L_INDICATOR:continue #only create cells out of 'real' clusters (that passed the filter)
    #             if cluster_i.ic_use_for_seeding: 
    #                 Cell(label=cluster_i.label,
    #                      cluster_seed=cluster_i,
    #                      stack=self)
    #                 #if DBG:print('{0} cell created from seed cluster_{1}'.format(cluster_i.get_color(),cluster_i.ctr))
    #                 if cluster_i.ctr in [203,223]:
    #                     print('break')
    #                 nb_cells_created +=1
    #         
    #         #add the other clusters to the newly created cells        
    #         for cluster_i in Cluster.l_curr_stack: 
    #             if cluster_i.label in Label.L_INDICATOR:continue 
    #             if not cluster_i.ic_use_for_seeding: 
    #                 if DBG:print('cluster(non seed) {0} is going to be added to label {1}'.format(cluster_i.ctr,cluster_i.get_color()))
    #                 label_cell = self.d_label_cell.get(cluster_i.label)
    #                 if label_cell:
    #                     cell = self.d_label_cell[cluster_i.label]
    #                     cell.add_cluster(cluster_i)
    #                 else:
    #                     Cell(label=cluster_i.label,
    #                      cluster_seed=cluster_i,
    #                      stack=self)
    #                 #if DBG:print('{0} cell created from seed cluster_{1}'.format(cluster_i.get_color(),cluster_i.ctr))
    #                     nb_cells_created +=1 
    #         
    #         if verbose:print("->",nb_cells_created,"created",end='')
    #         if verbose:print(['{0} from {1}'.format(i.label.color,i.cluster_seed.ctr) for i in self.cells])
    #         return
    #     
    def get_matching_cluster_from_stack_based_on_max_overlap(self,cluster):
            '''
            get the matching cluster based on maximum overlap
            if no overlap is found, None is returned
            '''
            bincount = np.bincount(self.clusterview[cluster.l_pixels].flatten())
            bincount[0] = 0
    
            return Cluster.d_ctr_cluster.get(np.argmax(bincount),None)
    
    def create_cells_through_cell_divisions(self,verbose=True):
        '''
        go over all cells and check if the cell should undergo cell division.  
        at first, the daughter cluster will get the 'candidate_new_cell' indication.  But if the same pattern occurs in the next timestep, a true cell division is executed (=confirmation process)
        the seed for the new cell will always be the biggest non seed cluster of the current cell.   Likewise, only the biggest non seed cluster can have the candidate new cell indication
        '''
        if verbose:print("###handle cell divisions for stack {0} ".format(self.timestep), end='')
        
        
        
        l_new_labels = []
        for cell_i in self.cells:
            if cell_i.label in Label.L_INDICATOR:
                continue
            if cell_i.check_cell_division():
                cluster_match_prev_stack = Stack.prev_stack.get_matching_cluster_from_stack_based_on_max_overlap(cell_i.cluster_largest_unseeded)
                cell_match_prev_stack  = cluster_match_prev_stack.get_cell()
                if cluster_match_prev_stack.candidate_new_cell:  
                    new_label=Label(timestep=Stack.prev_stack.timestep,label_parent_ctr=cell_i.label.ctr,cell=None,cell4D=None,ctr='')
                    new_cell_prev_stack = cell_match_prev_stack.execute_cell_division(label_forced=new_label,true_division_time=True)  #division in the previous stack
                    new_cell_curr_stack = cell_i.execute_cell_division(label_forced=new_label,true_division_time=False)                 #division in the current stack
                    new_cell_curr_stack.cluster_seed.seed_new_cell=False #was true in the previous stack, not anymore in the current stack
                    l_new_labels.append(new_label)
                    
                    # cell_match_prev_stack.label.update_lineage_name_cell_division() #extend lineage_name with d1
                    # cell_match_prev_stack.name_lineage = cell_match_prev_stack.label.name_lineage #take snapshot 
                    # cell_i.name_lineage = cell_match_prev_stack.label.name_lineage #take snapshot
                else:
                    cell_i.cluster_largest_unseeded.candidate_new_cell = True   
                    
        if verbose:print("({0} cell divisions {1})".format(len(l_new_labels),[i.color for i in l_new_labels])) 
        return
        
    def create_cells_through_cell_divisions_GT_corr(self,verbose=True):
        '''
        only divide cells that contain more than 1 GT_seeded cluster.  can even be more than 2 per cell
        we need to do some ad-hoc manipulations to fit the GT_cells nicely into the lineage tree (so lineager_feeder will correctly interpret it)
        # GT seeding does not trigger cell divisions in the previous stack.  it just does an cell division in the current stack , without the confirmation process
        # example : a normal cell division in the cluster summary looks likes this
        # t   color cluster_label         parent_label 
        # 1	blue	    4	    seed	    99	
        # 2	blue	    4	    seed	    99	        DIVIDED
        # 2	yellow+++	46	    seed	    4	        NEW CELL
        # so we need to first inherit cluster_label and parent_label from previous stack (as the default)
        # then go over the newly created cells and correct them : assign prev cluster label to the parent_label, and mark them as a 'new cell'
        '''
        if verbose:print("###handle cell divisions GT_seeded for stack {0} ".format(self.timestep), end='')
        
        l_new_labels = []
        l_cells_to_process = self.cells
        l_cells_set_aside = []
        while l_cells_to_process:
            for cell_i in l_cells_to_process:
                if not cell_i.cluster_largest_unseeded:continue
                if cell_i.cluster_seed.GT_seeded:
                    if cell_i.cluster_largest_unseeded.GT_seeded:
                        cluster_match_prev_stack = Stack.prev_stack.get_matching_cluster_from_stack_based_on_max_overlap(cell_i.cluster_largest_unseeded)
                        if not cluster_match_prev_stack:print('debug : cluster_match_prev_stack.get_cell() returned null ! not ok, cell division is cancelled. cell repr= ',cell_i.repr());continue
                        cell_match_prev_stack  = cluster_match_prev_stack.get_cell()
                        
                        
                        new_label=Label(timestep=Stack.prev_stack.timestep, label_parent_ctr=cell_match_prev_stack.label.ctr, cell=None, cell4D=None, ctr='')
                        if verbose:print('TEMP GT logic : forced cell division (cluster_largest_unseeded.ctr {2}) because of cluster_seed.ctr {1}->new label {0}, old label {3}'.format(new_label.ctr,cell_i.cluster_seed.ctr,cell_i.cluster_largest_unseeded.ctr,cell_i.label.ctr))
                        new_cell_curr_stack = cell_i.execute_cell_division(label_forced=new_label,true_division_time=True)                 #division in the current stack //GT does not trigger the prev stack
                        l_new_labels.append(new_label)

                        
                    #check for GT-seeds in the other clusters.  if found, promote to largest_unseeded and put on list to process
                    for ix,other_cluster_i in enumerate(cell_i.l_clusters_other):
                        print('in loop 1 cell_i.label={0}   other_cluster_i.GT_seeded={1}   other_cluster_i.ctr={2} {3}'.format(cell_i.label.color,other_cluster_i.GT_seeded, other_cluster_i.ctr,ix))
                        if other_cluster_i.GT_seeded:
                            
                            if other_cluster_i != cell_i.cluster_largest_unseeded:
                                print('swithching other to largest unseeded')
                                t = cell_i.cluster_largest_unseeded
                                cell_i.cluster_largest_unseeded = other_cluster_i
                                cell_i.l_clusters_other[ix] = t
       
                            l_cells_set_aside.append(cell_i)

                            continue
            l_cells_to_process = l_cells_set_aside
            l_cells_set_aside = []

        
        if verbose:
            gt_count = 0
            cell_count = 0
            for cell_i in self.cells:
                cell_count += 1
                if cell_i.cluster_seed.GT_seeded:
                    gt_count +=1
                    print("{0} has a gt seeded cluster_seed {1}=ok".format(cell_i.label.color,cell_i.cluster_seed.ctr))
                    for cluster_i in cell_i.get_clusters('no_cluster_seed'):
                        if cluster_i==cell_i.cluster_seed:continue
                        if cluster_i.GT_seeded:
                            print("warning , still a cell present with a GT cluster non seed",cluster_i.ctr)
                else:
                    print("warning: this cell {0} has no GT seed".format(cell_i.label.color))
            print("there a {0} cells of which {1} containing a GT seed".format(cell_count,gt_count))
                
            if self.timestep==3:cell_i.repr()
        if verbose:print("({0} cell divisions {1})".format(len(l_new_labels),[i.color for i in l_new_labels])) 
        return
        
    
    def filter_cells(self,ic_GT_corr=False,verbose=True):
        if verbose:print(">>>Filter cells<<<", end='')
        
        nb_filtered_cells = 0
        if not self.prev_stack:  #only first stack . filter on radius of pixel volume
            for cell_i in self.cells:
                if cell_i.cluster_seed.GT_seeded:continue
                
                if cell_i.get_radius() < Param.MIN_CELL_RADIUS:  #this should be higher because these cells will already be filtered out at cluster level, but for now
                    cell_i.replace_label(new_label=Label.FILTERED,trace='cell-filter')
                    nb_filtered_cells +=1
        
        if ic_GT_corr:
            for cell_i in self.cells:
                if not cell_i.GT_seeded:
                    cell_i.replace_label(new_label=Label.FILTERED,trace='cell-filter')
                    nb_filtered_cells +=1
                    if verbose:print("(GT-logic :",cell_i.label.color, "cell was filtered because the cluster_seed was not gt seeded)")
                    
        if verbose:print("(",nb_filtered_cells, "cells filtered)")
        return


    def create_cell4D(self,verbose=True):
        if verbose:print("###create cell4D for stack {0}".format(self.timestep),end='')
        
        l_cell4D_created = []
        #update the cluster information
        l_existing_labels = [i.label for i in Cell4D.l_cell4D]
        for cell_i in self.cells:
    #             if cell_i.label in Label.L_INDICATOR:
    #                 continue
            if cell_i.label not in l_existing_labels:
                cell4d = Cell4D(cell_i,self.timestep)  #create cell4D
                l_cell4D_created.append(cell4d)

        if verbose:print("({0} cell4D created {1})".format(len(l_cell4D_created), [i.get_color() for i in l_cell4D_created] ))
        return
    
    def update_cell4D(self,verbose=True):
        if verbose:print("###update cell4D info for stack {0}".format(self.timestep))
    
        for cell_i in self.cells:
            if cell_i.cell4D:
                cell_i.cell4D.incorporate_cell_data(cell_i) #update some info based on current stack
                cell_i.cell4D.lifetime +=1
       
        return
    

    
    def filter_cell4D(self,verbose=True):
        if verbose:print(">>>Filter cell4D<<<", end='')
        
        nb_filtered_cell4d=0
        if self.prev_stack: #not applicable for the first stack
            for cell4D_i in Cell4D.l_cell4D:
               ## if cell4D_i.initiating_cell.cluster_seed.GT_seeded:continue
                if (cell4D_i.lifetime < 1) and cell4D_i.has_died_in_this_timestep():  
                    cell4D_i.replace_label(new_label=Label.FILTERED,trace='cell4D-filter')
        
        if verbose:print("(",nb_filtered_cell4d, "cell4D filtered)")
        return

 
    def validate_clusters(self,df_nuclei,validation_start_time=1,verbose=True):
        '''
        this will validate clusters (and by cascade, also cell and cell4) based on a ground truth excel.
        beware : x, y and z in validation file are 1-indexed.  when manually added they are integers, otherwise they can be float (and will be floored)
        4 possibilities : 
        1) VAlIDATED  (TP): the ground truth nuclei falls inside an existing, unvalidated cluster
        2) OVERSEGMENTATION (FP): the detected cluster is not present in the ground truth 
        3) UNDERSEGMENTATION (FN): a detected cluster covers more than 1 ground truth cluster.  an extra fictive sphere will be created (small red ball)
        4) cluster NOT DETECTED (FN): a ground truth nuclei falls in a region where no cluster is detected.  an extra fictive sphere will be created (small red ball)
        :param df_nuclei: the dataframe with the ground truth nuclei
        '''
        if verbose:print("###validating the clusters based on the nuclei input (stack{0})".format(self.timestep),end='')
        
        def create_validation_error_cluster(s_desc):
            #validation clusters are not put in the clusterview so they will not interfere with the cluster detection algorithm;  they will show up in the tiff and excel
            validation_error_sphere  = Sphere(centre_ZYX=[floor(centre.z)-1,centre.y-1,centre.x-1],
                                           radius=4,
                                           stack=self)
            validation_error_cluster = Cluster(validation_error_sphere,stack=self)
            validation_error_cluster.set_label(Label.VALIDATION, trace='creation')
            validation_error_cluster.centre= [floor(centre.z)-1,centre.y-1,centre.x-1] 
            
            if s_desc=='UNDERSEGMENTATION':
                cluster_landed.validation_error = s_desc # the landed cluster gets marked with undersegmentation, not the validation sphere
                validation_error_cluster.set_val_info(centre.cell, [floor(centre.z)-1,centre.y-1,centre.x-1] , '')
            else:
                validation_error_cluster.set_val_info(centre.cell, [floor(centre.z)-1,centre.y-1,centre.x-1] , s_desc)
            
            self.nb_FN +=1
            
            return
    
        if not len(df_nuclei): return      
                                                                           
        for _,centre in df_nuclei[df_nuclei.time==self.timestep + validation_start_time - 1].iterrows():
            Stack.set_validation_names.add(centre.cell)
            try:
                ctr_cluster_landed = self.clusterview[int(floor(centre.z)-1),int(centre.y-1),int(centre.x-1)]
            except:
                print("warning : validation centre {0} lands outside the image. skipping...".format(centre))
                continue
            cluster_landed = Cluster.d_ctr_cluster.get(ctr_cluster_landed,None)
            
            if cluster_landed and cluster_landed.label not in Label.L_INDICATOR:
                if cluster_landed.cell.validation_name: #already validated
                    create_validation_error_cluster('UNDERSEGMENTATION')
                else:
                    cluster_landed.set_val_info(centre.cell, [floor(centre.z)-1,centre.y-1,centre.x-1] , "")
                    self.nb_TP +=1
            else:
                create_validation_error_cluster('CELL NOT DETECTED') 
        #<<<end validation centres loop
                
        for cell_i in self.cells:
            if not cell_i.validation_name and cell_i.label not in Label.L_INDICATOR:
                cell_i.cluster_seed.set_val_info("","",'OVERSEGMENTATION') #the cluster seed will get marked with oversegmentation
                self.nb_FP +=1
                
        #update the total stack stats
        Stack.tot_TP += self.nb_TP
        Stack.tot_FN += self.nb_FN
        Stack.tot_FP += self.nb_FP

        if verbose:print('...(TP={0} FN={1} FP={2})'.format(self.nb_TP,self.nb_FN,self.nb_FP))
        return
     
    
    def write_excel_summary_3D(self,filehandler,initialize=False, verbose=True):
        
        if verbose:print('###cluster_summary_over_all_stacks (append)')
        filehandler.d_save_info['f_name'] = 'cluster_summary'
        filehandler.store_f_name('cluster_summary',os.path.join(filehandler.get_save_location(),'cluster_summary.csv'))
        if initialize:
            filehandler.init_f_name('cluster_summary') #prevent appending to old file
            return
        
        ls_columns = [
            "stack_nb",
            "color",
            "cluster_label",
            "cluster_counter",
            "cell_rank",
            "label_parent",
            "division",
            "no division reason",
           
            "validation_error",
            "validation_name",
            "validated_center",
            "z",
            "y",
            "x",
            "label_hist",
            "radius_seed",
            "radius_cluster",
            "pixels_count",
            
            "name_lineage",
            "ctr_lineage",
            
            "overlapping_clusters_ctr",
           
            "overlap_max_DT",  #the maximum DT value in the overlap zone (=the overlap you observe)
            "DT_max_threshold", #the value that is compared to overlap_max_DT (=the threshold you calculate)
            "ic_divide_cell",        #will be = True is overlap_max_DT < DT_max_threshold (=the final decision)
            
            "seeded_by_GT",
                             
            "spheres(ctr,centre,radius)"
               ]
        
        #put filler in empty ls_columns
        
        #intialize ls_columns
        color=[]
        cluster_labels = []
        cluster_counter = []
        
        cell_rank=[]
        parent_label = []
        stack_nb = []
        division = []
        no_division_reason=[]
        
        validation_error = []
        validated_name= []
        validated_centre = []
        
        
        
        z = []
        y = []
        x = []
        
        label_history = []
       
        radius_seed =[]
        radius_cluster=[]
        pixels_count = []
        
        name_lineage = []
        ctr_lineage =[]
        
        overlapping_clusters_ctr= []
        overlap_max_DT=[]
        
        DT_max_threshold=[]
        ic_divide_cell=[]
        
        seeded_by_GT = []
    
        spheres = []
        
        #compose ls_columns
        def compose_sorted_cluster_list():
            
            l_cluster_sort = []
            for cluster_i in Cluster.d_ctr_cluster.values():
                if cluster_i.stack !=self:
                    continue #filter : only clusters of this stack
                l_cluster_sort.append((cluster_i,cluster_i.label.ctr,-cluster_i.get_volume())) #sorting on [label (asc),volume(desc)]
            
            return [i for i,_,_ in sorted(l_cluster_sort, key=lambda l_cluster_sort: (l_cluster_sort[1],l_cluster_sort[2]))]  
        
        
        for cluster in compose_sorted_cluster_list():
            
            # cluster_labels.append(cluster.label.ctr)
            cluster_labels.append("{:02d}".format(cluster.label.ctr))
            
            if cluster.label==Label.EXTERIOR_CTR:
                color.append(Label.COLOR_EXTERIOR)
            else:
                color.append(cluster.label.color)
            
            if cluster.cell:
                if cluster.cell.cluster_seed==cluster: 
                    cell_rank.append('seed')
                elif cluster.cell.cluster_largest_unseeded==cluster:
                    cell_rank.append('largest unseeded')
                else:
                    cell_rank.append('')
            else:
                cell_rank.append('')
            
            parent_label.append(cluster.label.label_parent_ctr)
            stack_nb.append(cluster.stack.timestep)
            
            if cluster.seed_new_cell:division.append('NEW CELL')
            elif cluster.candidate_new_cell:division.append('CANDIDATE')
            elif cluster.seed_divided:division.append('DIVIDED')
            else: division.append('')
            
            no_division_reason.append(Cluster.d_reason_code_no_division_desc.get(cluster.reason_cd_no_division,'No description'))
            
            validation_error.append(cluster.validation_error)
            validated_name.append(cluster.validation_name)
            validated_centre.append(cluster.validation_centre)
            
            cluster_counter.append(cluster.ctr)
            
            z.append(cluster.centre[0] + 1)# from index to number
            y.append(cluster.centre[1] + 1)
            x.append(cluster.centre[2] + 1)
            
            
            
            label_history.append(cluster.label_history)
            
            radius_seed.append(cluster.spheres[0].radius)
            radius_cluster.append(cluster.get_radius())
            
            pixels_count.append(len(cluster.l_pixels[0]))
            
            if not cluster.cell:
                name_lineage.append("")
                ctr_lineage.append("")
            else:
                name_lineage.append(cluster.label.name_lineage)
            
                if cluster.label.name_lineage not in Stack.d_name_lineage_ctr:   #just give a unique number to every name
                    Stack.d_name_lineage_ctr[cluster.label.name_lineage]=len(Stack.d_name_lineage_ctr) + 1
                ctr_lineage.append(Stack.d_name_lineage_ctr[cluster.label.name_lineage])
            
            overlapping_clusters_ctr.append([i for i in cluster.d_clusterctr_overlap.keys()])
            
            #report the overlap info for clusters involved in cell division
            if cluster.is_largest_non_seed_cluster() or cluster.seed_new_cell:
                cluster_seed = cluster.cell.cluster_seed
                overlap_max_DT.append(cluster.d_clusterctr_overlap.get(cluster_seed.ctr)['overlap_max_DT'] if cluster.d_clusterctr_overlap.get(cluster_seed.ctr) else 0)
                min_radius = min(cluster.spheres[0].radius, cluster_seed.spheres[0].radius)
                DT_max_threshold.append(self.get_overlap_maxDT_threshold(cluster.spheres[0],cluster_seed))
                if  overlap_max_DT[-1] < DT_max_threshold[-1]:       
                    ic_divide_cell.append('True')
                else:
                    ic_divide_cell.append('False')
            else:
                overlap_max_DT.append('')
                DT_max_threshold.append('')
                ic_divide_cell.append('')
               
            seeded_by_GT.append(cluster.GT_seeded)
            
            spheres.append([(i.ctr,i.centre,i.radius) for i in cluster.spheres])
            
        #<--end cluster loop (all clusters of this stack)
      
        df = pd.DataFrame({ls_columns[0]:stack_nb,
                           
                           ls_columns[1]:color,
                           ls_columns[2]:cluster_labels,
                           ls_columns[3]:cluster_counter,
                           
                           ls_columns[4]:cell_rank,
                           
                           ls_columns[5]:parent_label,
                           
                           ls_columns[6]:division,    
                           ls_columns[7]:no_division_reason,
                           
                           ls_columns[8]:validation_error,
                           ls_columns[9]:validated_name,
                           ls_columns[10]:validated_centre,
                           
                           ls_columns[11]:z,
                           ls_columns[12]:y,
                           ls_columns[13]:x,
                           
                           ls_columns[14]:label_history,                         
                           ls_columns[15]:radius_seed,
                           ls_columns[16]:radius_cluster,
                           ls_columns[17]:pixels_count,  
                           ls_columns[18]:name_lineage,
                           ls_columns[19]:ctr_lineage,
                                                    
                           ls_columns[20]:overlapping_clusters_ctr,
                           ls_columns[21]:overlap_max_DT,
                           ls_columns[22]:DT_max_threshold,
                           ls_columns[23]:ic_divide_cell,
                           ls_columns[24]:seeded_by_GT,
                           
                           ls_columns[25]:spheres
                       })
        
        
        #save the data
        filehandler.set_save_info('nb_stack','')
        filehandler.save_data(data=df,csv_columns=ls_columns,file_ext='csv',verbose=False)
        filehandler.reset_save_info()
        
        return
    
    def write_excel_validation(self,filehandler,print_total=False,verbose=True):
        
        if verbose:print('###write excel validation (append stack)')
        filehandler.d_save_info['f_name'] = 'validation_cells'
        
        ls_columns = ["timestep",
                   "TP",
                   "FN",
                   "FP",
               ]
        
        #put filler in empty ls_columns
        
        #intialize ls_columns
        timestep = []
        TP=[]
        FN = []
        FP = []
        
        
        #append the data
        timestep.append(self.timestep);TP.append(self.nb_TP);FN.append(self.nb_FN);FP.append(self.nb_FP)
        
        if print_total:
            timestep.append('TOTAL');TP.append(Stack.tot_TP);FN.append(Stack.tot_FN);FP.append(Stack.tot_FP)
            
            TP_cell4D= 0
            FN_cell4d= 0
            FP_cell4d = 0
            set_match_w_cell4D = set()
            for cell4D_i in Cell4D.l_cell4D:
                if cell4D_i.label in Label.L_INDICATOR:
                    continue
                if len(cell4D_i.l_validation_names)>0:
                    for validation_name_i in cell4D_i.l_validation_names:
                        set_match_w_cell4D.add(validation_name_i)
                    TP_cell4D +=1 
                else:
                    FP_cell4d +=1
                
            set_diff= Stack.set_validation_names.difference(set_match_w_cell4D)
            FN_cell4d = len(set_diff)
                
            timestep.append('');TP.append('');FN.append('');FP.append('')
            timestep.append('cell4D');TP.append(TP_cell4D);FN.append(FN_cell4d);FP.append(FP_cell4d)
            timestep.append('detail->');TP.append('');FN.append('');FP.append('')
            
            for cell4D_i in Cell4D.l_cell4D:
                if cell4D_i.label in Label.L_INDICATOR:
                    continue
                timestep.append(cell4D_i.label.ctr)
                TP.append(cell4D_i.label.color)
                if len(cell4D_i.l_validation_names)>0:
                    FN.append(cell4D_i.l_validation_names)
                else:
                    FN.append('NO MATCH')
                FP.append('')

            for validation_name_i in set_diff:
                    timestep.append('NO MATCH')
                    TP.append('')
                    FN.append(validation_name_i)
                    FP.append('')                

        #compose dataframe
        df = pd.DataFrame({ls_columns[0]:timestep,
                           ls_columns[1]:TP,
                           ls_columns[2]:FN,
                           ls_columns[3]:FP
                        })
        
        
        #save the data

        filehandler.set_save_info('nb_stack','')
        filehandler.save_data(data=df,csv_columns=ls_columns,file_ext='csv',verbose=False)
        filehandler.reset_save_info()
  
        
        return
              
    def construct_RGBA_clusterview(self,filehandler=None,temp_folder=False,
                                show_cells_only=False,df_validation=None,validation_start_time=1,verbose=True):
        
        if verbose:print("###constructing the RGBAview for stack {0}".format(self.timestep))
        
        def print_cluster_id_RGB(print_string,cluster):
            #getting starting position upper left corner
            box_size,_= Param.a_empty.shape
            z_cluster,y_cluster,x_cluster = cluster.centre
            x_left_upper = max(0,x_cluster - round(box_size/2))
            y_left_upper = max(0,y_cluster - round(box_size/2))
            length_total = len(print_string) * box_size
            if x_left_upper + length_total > Param.X_DIM:
                x_left_upper -= x_left_upper + length_total - Param.X_DIM +1 
            if y_left_upper + box_size > Param.Y_DIM:
                y_left_upper -= y_left_upper + box_size - Param.Y_DIM +1
                
            #printing the letters of the string onto the image
            x,y,z = x_left_upper,y_left_upper,z_cluster
            for letter_i in print_string:
                a_print_area = RGBA_clusterview[z,y:y+box_size,x:x+ box_size,:]
                a_mask = Param.d_char_typeletter.get(letter_i,Param.a__)
                if a_print_area.shape[0:2]==a_mask.shape:
                    a_print_area[np.where(a_mask==1)] = Label.D_COLOR_RGBA['black'] 
                x = x + box_size
            
            return 
            
           
        if temp_folder:
            filehandler.extend_save_info(extra_dir_1 = Param.NAME_TEMP_FOLDER)
        filehandler.d_save_info['f_name']='RGBA_clusterview'
        filehandler.d_save_info['nb_stack']=self.timestep
        
        z,y,x = self.clusterview.shape
        RGBA_clusterview = np.zeros((z,y,x,4),'uint8')                                      
        RGBA_clusterview.fill(255)                                                         #white background
        
        Z,Y,X = np.where(self.img==True)                                                    #superimposing the black regions (=  membrane)
        RGBA_clusterview[Z,Y,X,:] = Label.D_COLOR_RGBA['black'] 
        
        
        for sphere in self.l_spheres_filtered:                                                      #superimposing the filtered spheres
            if show_cells_only:
                continue
            Z,Y,X = sphere.l_pixels
            RGBA_clusterview[Z,Y,X,:] = Label.D_COLOR_RGBA[Label.COLOR_FILTERED_SPHERE] 
            z,y,x = sphere.centre
            if sphere.reason_cd_filtering == '1':
                RGBA_clusterview[z,y,x,:] = Label.D_COLOR_RGBA[Label.COLOR_FILTERED_SPHERE_CENTRE_EXTERIOR] 
            elif sphere.reason_cd_filtering == '2':
                RGBA_clusterview[z,y,x,:] = Label.D_COLOR_RGBA[Label.COLOR_FILTERED_SPHERE_CENTRE_OVERLAP_MASK] 
            elif sphere.reason_cd_filtering == '3':
                RGBA_clusterview[z,y,x,:] = Label.D_COLOR_RGBA[Label.COLOR_FILTERED_SPHERE_CENTRE_MIN_SEED_RADIUS] 
            elif sphere.reason_cd_filtering == '4':
                RGBA_clusterview[z,y,x,:] = Label.D_COLOR_RGBA[Label.COLOR_FILTERED_SPHERE_CENTRE_Z_TOP_BOTTOM] 

                
              
        for cluster in self.get_clusters():    
                                                                                        #superimposing the clusters
            Z,Y,X = cluster.l_pixels
            if cluster.label==Label.EXTERIOR:
                RGBA_clusterview[Z,Y,X,:]=np.asarray(Label.D_COLOR_RGBA.get(Label.COLOR_EXTERIOR))
            else: 
                RGBA_clusterview[Z,Y,X,:]=np.asarray(Label.D_COLOR_RGBA.get(cluster.label.color,((255,255,255,255))))
            
            if show_cells_only:
                continue
            #if cluster.label == Label.VALIDATION:continue  
                                                   
            
            for ix,sphere in enumerate(cluster.spheres):          #superimposing the spheres (dots)
                z,y,x = sphere.centre
                RGBA_clusterview[z,y,x,:] = Label.D_COLOR_RGBA['black']  #a spheres gets a dot by default
                            
            z,y,x = cluster.centre
            color_centre_cross='black' # superimposing the clusters (crosses)
            if cluster.is_seed():#default black cross for clusters, #size shrinks with cluster importance (seed-> largest non seeds -> others)
                size_cross = 5  
            elif cluster.is_largest_non_seed_cluster():
                size_cross = 4  
            else:
                size_cross = 3
                
            if cluster.seed_new_cell or cluster.candidate_new_cell:
                print_string = 'O' #big square
                print_cluster_id_RGB(print_string,cluster)
                color_centre_cross='green'  #green cross for clusters involved in cell division
                size_cross = 6  # a new cell gets an extra big cross
            if cluster.candidate_new_cell:
                print_string = 'o' #small circle
                print_cluster_id_RGB(print_string,cluster)
                color_centre_cross='green'  #green cross for clusters involved in cell division

                    

            x_min_clip = max(x-(size_cross-1),0)
            x_max_clip = min(x+size_cross,Param.X_DIM-1)
            y_min_clip = max(y-(size_cross-1),0)
            y_max_clip = min(y+size_cross,Param.Y_DIM-1)
            RGBA_clusterview[z,y,x_min_clip:x_max_clip,:] = Label.D_COLOR_RGBA[color_centre_cross]  
            RGBA_clusterview[z,y_min_clip:y_max_clip,x,:] = Label.D_COLOR_RGBA[color_centre_cross]  
            if cluster.seeding_cell_label:
                RGBA_clusterview[z,y,x,:] = Label.D_COLOR_RGBA['white']  #cluster that was seeded gets a white dot
    #             if cluster.validation_error == 'OVERSEGMENTATION':          #underline with red line in case of oversegmentation
    #                 RGBA_clusterview[z,y_max_clip,x_min_clip:x_max_clip,:] = Label.D_COLOR_RGBA['red']
                            
        #<---cluster loop 1
        
        
        
            
        if not show_cells_only:
            for cluster in self.get_clusters():  
                for d_overlap_cluster in cluster.d_clusterctr_overlap.values():    #superimposing the overlap l_pixels
                    if len(d_overlap_cluster['overlap_pix']):
                        Z,Y,X = d_overlap_cluster['overlap_pix']
                        RGBA_clusterview[Z,Y,X,:] = Label.D_COLOR_RGBA['black_light']     
                #<--end of overlapping clusters
            #<--cluster loop 2  
                    
                    
            if df_validation is not None:                                                                         #superimposing the cluster nuclei for validation
                size_sword = 3
                for _, centre in df_validation[(df_validation.time + validation_start_time - 1)==self.timestep].iterrows():
                    RGBA_clusterview[floor(centre.z)-1,floor(centre.y),floor(centre.x)-size_sword:floor(centre.x)+round(centre.radius),:] = Label.D_COLOR_RGBA['white']  
                    RGBA_clusterview[floor(centre.z)-1,floor(centre.y)-size_sword:floor(centre.y)+size_sword+1,floor(centre.x),:] = Label.D_COLOR_RGBA['white'] 
         
          
                                                                                               
            for cluster_i in self.get_clusters():                                                         #superimposing label and cluster ctr numbers
                print_string = '  ' + str(cluster_i.label.ctr) + '_' + str(cluster_i.ctr)
                print_cluster_id_RGB(print_string,cluster_i)
            
            if cluster_i.validation_error == 'OVERSEGMENTATION':  
                z,y,x = cluster_i.centre
                size_cross = 6 
                x_min_clip = max(x-(size_cross-1),0)
                x_max_clip = min(x+size_cross,Param.X_DIM-1)
                y_max_clip = min(y+size_cross,Param.Y_DIM-1)        #underline with red line in case of oversegmentation
                RGBA_clusterview[z,y_max_clip,x_min_clip:x_max_clip,:] = Label.D_COLOR_RGBA['red']
            #<--cluster loop 4   
        
        #validation errors and cell/label numbers
            for cluster in self.get_clusters_for_label(Label.VALIDATION):  
                Z,Y,X = cluster.get_l_pixels()
                RGBA_clusterview[Z,Y,X,:]=np.asarray(Label.D_COLOR_RGBA.get(cluster.label.color,((255,255,255,255))))  #superimposing the validation clusters
            #<--cluster loop 3 
        
        if filehandler:
            filehandler.save_data(RGBA_clusterview,file_ext='tif',resolution='uint8',use_nb_stack=True,verbose=False)  
        else:
            Stack.l_RGBA_clusterview.append(RGBA_clusterview)  
        
        #Also store clusterview (used for tracking CTC)
        filehandler.d_save_info['f_name']='raw_clusterview'
        filehandler.save_data(self.clusterview.astype('uint16'),file_ext='tif',resolution='uint16',use_nb_stack=True,verbose=False)  
        
        save_cell_view=True
        if save_cell_view:
            self.compose_label_view()
            filehandler.d_save_info['f_name']='labelview'
            filehandler.save_data(self.labelview.astype('uint16'),file_ext='tif',resolution='uint16',use_nb_stack=True,verbose=False)  
        
        if temp_folder:
            filehandler.pop_save_info()
            
        return
    
    
    
    @staticmethod
    def save_4D_RGBA_clusterview(filehandler,nb_channels=1,extra_channel=None,verbose=False):
        
        if verbose:print("###saving the RGB cluster view (fiji)")
        
        filehandler.d_save_info['f_name']= 'RGB_cellview'
        
        #concatenate tif files from temp folder
        filehandler.extend_save_info(extra_dir_1 = Param.NAME_TEMP_FOLDER)
        temp_folder_select = filehandler.get_save_location() +  '*RGBA*.tif'
        l_files = sorted(glob.glob(temp_folder_select))
        l_4D_RBGA_clusterview = []     
        
        filehandler_empty = FileHandler()
        for path in l_files: 
            filehandler_empty.d_load_info['f_name']=path
            l_4D_RBGA_clusterview.append(filehandler_empty.load_tif(add_extension=False,verbose=False))
        
        a_4D_RGBA_clusterview = np.asarray(l_4D_RBGA_clusterview)
        
        #create channeldimension,  CTZYXS format
        t,z,y,x,rgba = a_4D_RGBA_clusterview.shape
        a_4D_RGBA_clusterview = a_4D_RGBA_clusterview.reshape((nb_channels,int(t/nb_channels),z,y,x,rgba)) 
        
        #add extra channel (use for raw data for example)
        if len(extra_channel):
            if extra_channel.shape[-1] ==3:
                extra_channel = np.insert(extra_channel,3,255,axis=4)
            t,z,y,x,rgba = extra_channel.shape
            if str(extra_channel.dtype) not in ('uint8'):extra_channel=img_as_ubyte(extra_channel) 
            extra_channel = extra_channel.reshape((1,t,z,y,x,rgba)) 
            a_4D_RGBA_clusterview = np.append(a_4D_RGBA_clusterview ,extra_channel,axis=0)
        
        #save data
        filehandler.pop_save_info()
        if len(a_4D_RGBA_clusterview):
            filehandler.save_data(a_4D_RGBA_clusterview,file_ext='tif',RGB=True,resolution='uint8',verbose=False)
        
        return
    

    @staticmethod
    def agg_temp_viz_files(filehandler,search_string='clusterview',blank_filtered=True,verbose=False):
        
        if verbose:print("###saving the {0} (fiji)".format(search_string))
        
        filehandler.d_save_info['f_name']= search_string 
        filehandler.store_f_name(search_string,os.path.join(filehandler.get_save_location(),'{0}.tif'.format(search_string)))
        
        #concatenate tif files from temp folder
        filehandler.extend_save_info(extra_dir_1 = Param.NAME_TEMP_FOLDER)
        temp_folder_select = filehandler.get_save_location() +  '*{0}*.tif'.format(search_string)
        l_4D_view = []     
        filehandler_empty = FileHandler()
        for path in sorted(glob.glob(temp_folder_select)): 
            filehandler_empty.d_load_info['f_name']=path
            l_4D_view.append(filehandler_empty.load_tif(add_extension=False,verbose=False))

        z,y,x = l_4D_view[0].shape
        t = len(l_4D_view)
        a_4D_view = np.asarray(l_4D_view).reshape((t,z,y,x))

        if blank_filtered:
            a_4D_view = np.where((a_4D_view==Label.FILTERED_CTR),0,a_4D_view)

        #save data
        filehandler.pop_save_info()
        if len(a_4D_view):
            filehandler.save_data(a_4D_view.astype('uint16'),file_ext='tif',RGB=False,verbose=verbose)

        return    
    
       
    def save_3Dplot(self,d_save_info):
        
        plot3D_stack_with_overlay(na_3D=self.img,
                                  intensity_threshold=None,
                                  overlay_image=self.exterior_mask,
                                  overlay_type='binary',
                                  a_shapes=self.get_clusters(),
                                  d_save_info=d_save_info,
                                  view_init = None,
                                  label_exterior_cell=Label.EXTERIOR_CTR,
                                  shape_by_shape=True)
        
  
    def check_consistency_datamodel(self,verbose=True):
        if verbose:print("...checking consistency of the datamodel for stack{0}...".format(self.timestep))
        
        nb_inconsistencies = 0
        #self.cells contains all cells from this stack.  => no check, must be true
        
        #Cell checks
        l_label = []
        for cell_i in self.cells:
            if cell_i.label:
                if cell_i.label in [Label.EXTERIOR,Label.UNASSIGNED]: #at this point, these labels should not be present anymore
                    print('Inconsistency : cell {0} (with seed cluster {1}) has an invalid label [Label.EXTERIOR,Label.UNASSIGNED]'.format(cell_i.get_color(),cell_i.cluster_seed.ctr))
                    nb_inconsistencies +=1;continue
                if cell_i.label in Label.L_INDICATOR:
                    continue #only reals label need to be RI consistent
                if cell_i.label in l_label:
                    print('Inconsistency : cell {0} (with seed cluster {1}) has the same label as another cell in the stack'.format(cell_i.get_color(),cell_i.cluster_seed.ctr))
                    nb_inconsistencies +=1;continue
                l_label.append(cell_i.label)
            else:
                print('Inconsistency : cell(with seed cluster {0} has no label'.format(cell_i.cluster_seed.ctr))
                nb_inconsistencies +=1;continue
            
            if not cell_i.stack:
                print('Inconsistency : cell(with seed cluster {0} it not attached to a stack')
                nb_inconsistencies +=1;continue
                
            if cell_i.cell4D:
                if cell_i.label != cell_i.cell4D.label:
                    print('Inconsistency : cell {0} (with seed cluster {1}) is connected to a cell4D with label {2} '.
                          format(cell_i.get_color(),cell_i.cluster_seed.ctr,cell_i.cell4D.label.color))
                    nb_inconsistencies +=1;continue
            else:
                print('Inconsistency : cell {0} (with seed cluster {1}) is not connected to a cell4D'.format(cell_i.get_color(),cell_i.cluster_seed.ctr))
                nb_inconsistencies +=1;continue
            
            l_clusters = cell_i.get_clusters()
            if l_clusters:
                for cluster_i in l_clusters:
                    if cluster_i.label != cell_i.label:
                        print('Inconsistency : cell {0} (with seed cluster {1}) has a different label than cluster {2}_{3}'.
                              format(cell_i.get_color(),cell_i.cluster_seed.ctr,cluster_i.get_color(),cluster_i.ctr))
                        nb_inconsistencies +=1;continue
                    if cluster_i != cell_i.cluster_seed:
                        if cluster_i.seed_new_cell:
                            print('Inconsistency : cell {0} (with seed cluster {1}) has non-seed cluster {2}_{3} that is marked as a seed of new cell (cluster.seed_new_cell)'.
                              format(cell_i.get_color(),cell_i.cluster_seed.ctr,cluster_i.get_color(),cluster_i.ctr))
                        if cluster_i.seed_divided:
                            print('Inconsistency : cell {0} (with seed cluster {1}) has non-seed cluster {2}_{3} that is marked as seed of divided cell (cluster.seed_divided)'.
                              format(cell_i.get_color(),cell_i.cluster_seed.ctr,cluster_i.get_color(),cluster_i.ctr))
            else:
                print('Inconsistency : cell {0} (with seed cluster {1}) has no clusters'.format(cell_i.get_color(),cell_i.cluster_seed.ctr))
                nb_inconsistencies +=1;continue
            
            if cell_i.l_clusters_other:
                if not cell_i.cluster_largest_unseeded:
                    print('Inconsistency : cell {0} (with seed cluster {1}) has a list of other clusters {2} but no largest unseeded cluster'.
                          format(cell_i.get_color(),cell_i.cluster_seed.ctr,[i.label.ctr for i in cell_i.l_clusters_other]))
                    nb_inconsistencies +=1;continue
                else:
                    for cluster_other_i in cell_i.l_clusters_other:
                        if cluster_other_i.get_seed_radius() > cell_i.cluster_largest_unseeded.get_seed_radius():
                            print('Inconsistency : cell {0} (with seed cluster {1}) has a cluster {2}(seed radius {4})in its other list, which is bigger than the largest unseeded cluster {3}(seed radius {5})'.
                                  format(cell_i.get_color(),
                                         cell_i.cluster_seed.ctr,
                                         cluster_other_i.ctr,
                                         cell_i.cluster_largest_unseeded.ctr,
                                         cluster_other_i.get_seed_radius(),
                                         cell_i.cluster_largest_unseeded.get_seed_radius()))
                            nb_inconsistencies +=1;continue
               
            
        #cluster check
        for cluster_i in Cluster.d_ctr_cluster.values():
            if not cluster_i.stack:   
                print('Inconsistency : cluster {0}  is not attached to a stack'.format(cluster_i.ctr)) 
                nb_inconsistencies +=1;continue
            if cluster_i.stack != self:
                continue
            
            if not cluster_i.cell:
                label = cluster_i.label
                if label:
                    if label not in Label.L_INDICATOR:
                        print('Inconsistency : cluster {0} (with label {1}) is not attached to a cell'.format(cluster_i.ctr,cluster_i.get_color())) 
                        nb_inconsistencies +=1;continue
                else:
                    print('Inconsistency : cluster {0} is not attached to a cell, and it also has no label'.format(cluster_i.ctr)) 
                    nb_inconsistencies +=1;continue
        
            for sphere_i in cluster_i.spheres:
                if not sphere_i.cluster:
                    print('Inconsistency : sphere{0} has no reciprocal link to its cluster {1}_{2}'.format(sphere_i.ctr,cluster_i.ctr,cluster_i.get_color()))
                    nb_inconsistencies +=1;continue
            
            if Cluster.d_ctr_cluster.get(cluster_i.ctr) != cluster_i:
                print('Inconsistency : cluster {0}_{1} => Cluster.d_ctr_cluster.get(cluster_i.ctr) != cluster_i'.format(cluster_i.ctr,cluster_i.get_color())) 
                nb_inconsistencies +=1;continue
            
            if Stack.curr_stack == self:
                if cluster_i not in Cluster.l_curr_stack:
                    print('Inconsistency : cluster {0}_{1} not in Cluster.l_curr_stack'.format(cluster_i.ctr,cluster_i.get_color())) 
                    nb_inconsistencies +=1;continue

            
            if cluster_i.candidate_new_cell:
                if cluster_i.cell.cluster_largest_unseeded!=cluster_i:
                    print('Inconsistency : cluster {0}_{1} is a candidate new cell, but not a cluster_largest_unseeded'.format(cluster_i.ctr,cluster_i.get_color())) 
                    nb_inconsistencies +=1;continue
                    
                    
        # #lineage tree check
        d_new_cell ={}      #key = parent_label_ctr ; value = new cell
        d_divided_cell = {} #key = label_ctr        ; value = divided cell
        for cell_i in self.cells:
            if cell_i.cluster_seed.seed_new_cell:
                parent_label_ctr = cell_i.cluster_seed.label.label_parent_ctr
                other_new_cell = d_new_cell.get(parent_label_ctr)
                if other_new_cell:
                    print('Inconsistency : cell {0} (with seed cluster {1}) is a new cell, but has the same parent label {4} as another new cell {2} (with seed cluster {3}).  double cell division ??'.
                    format(cell_i.get_color(),cell_i.cluster_seed.ctr,other_new_cell.get_color(),other_new_cell.cluster_seed.ctr,parent_label_ctr))
                    nb_inconsistencies +=1;continue
                else:
                    d_new_cell[parent_label_ctr] = cell_i
            elif cell_i.cluster_seed.seed_divided:
                d_divided_cell[cell_i.cluster_seed.label.ctr] = cell_i
        
        for parent_label_ctr_i, cell_i in d_new_cell.items():
            if not d_divided_cell.get(parent_label_ctr_i):
                print('Inconsistency : cell {0} (with seed cluster {1}) is a new cell, but the parent_label of this cell {2} is not found as a divided cell in this stack.  parentless division ??'.
                format(cell_i.get_color(),cell_i.cluster_seed.ctr,parent_label_ctr_i))
                nb_inconsistencies +=1;continue
                
        for label_ctr_i, cell_i in d_divided_cell.items():
            if not d_new_cell.get(parent_label_ctr_i):
                print('Inconsistency : cell {0} (with seed cluster {1}) is a divided cell, but its label {2} is not found as a parent label of a divided cell. fake division ??'.
                format(cell_i.get_color(),cell_i.cluster_seed.ctr,label_ctr_i))
                nb_inconsistencies +=1;continue
                
                    
        return nb_inconsistencies


    def update_lineage_names_cell_division(self,verbose=True):
        #a full update (complete stack) after the consistency checks. lineage names are only used to output to cluster_summary so they do not need to be maintained up-to-date during processing all the time
        #the entry point is a new cell (a cell division gives rise to a 'divided' cell and a 'new cell', the new cell draws a new label, the divided cell keeps the old label, but the lineage name is extended with d1
        if verbose:print('...update_lineage_names...')
        for cell_i in self.cells:
            if cell_i.cluster_seed.seed_new_cell:
            
                label_new_cell     = cell_i.label
                label_divided_cell = Label.d_ctr_label[label_new_cell.label_parent_ctr]

                label_new_cell.name_lineage      = label_divided_cell.name_lineage + 'd2'
                label_divided_cell.name_lineage += 'd1'
    
        return
            
        
    @staticmethod
    def check_memory_allocation(verbose=True):
        if verbose:print('...checking for memory leaks...')
        for cluster_i in Cluster.d_ctr_cluster.values():
            if not cluster_i.stack in [Stack.curr_stack,Stack.prev_stack]:
                print('Memory leak: cluster {0}_{1} in Cluster.d_ctr_cluster is connected to a stack other than the current or previous. timestep stack={2}'.
                      format(cluster_i.ctr,cluster_i.get_color(),cluster_i.stack.timestep)) 
                
        return
        
    def repr(self,detail_cells=False,output=sys.stdout):
            
        print('{0}info on Stack {1}{2}'.format
              ('*',self.name,'-->'),file=output)
        print('{0}img.shape ={1}'.format('\t', self.img.shape),file=output)
        print('{0}stack ={1}'.format('\t', self.timestep),file=output)
    #         print('{0}{1} spheres attached to stack'.format('\t',self.ctr_spheres ),file=output)
    #         print('{0}spheres (added in this order) ='.format('\t'),['{0}{1}'.format(i.centre,'\n') for i in self.spheres],file=output)

        
        exterior_counter=0;interior_counter=0
        for cluster in self.get_clusters():
            if cluster.label==Label.EXTERIOR_CTR:
                exterior_counter +=1
            else:
                interior_counter +=1
        print('{0}-->{1} interior clusters, {2} exterior clusters'.format('\t',interior_counter,exterior_counter),file=output)
            
        if len(self.exterior_mask):
            print('{0}exterior_mask present {1}'.format('\t',self.exterior_mask.shape),file=output)
        else:
            print('{0}no exterior_mask attached'.format('\t'),file=output)
        if len(self.dist_transf_carved):
            print('{0}distance transform present {1}'.format('\t',self.dist_transf_carved.shape),file=output)
        else:
            print('{0}no distance transform attached'.format('\t'),file=output)

        if detail_cells:
            for cluster in self.get_clusters():
                print('detailed info about cluster ', cluster.ctr, ': ','\n ----------------------------',file=output)
                cluster.repr(output)
        
        return output

          
    @classmethod
    def how_many(cls):
        print("We have {:d} stack.".format(cls.ctr))   
        
        
class Sphere:
    ctr = 0


    def __init__(self, centre_ZYX,radius,stack):
        Sphere.ctr +=1
        self.ctr = Sphere.ctr
        self.centre = centre_ZYX 
        self.radius = radius
        self.cluster = None
        self.stack = stack
        
        self.reason_cd_filtering = ''            
        
        self.l_pixels = get_3D_sphere_coordinates(a_coordZYX=centre_ZYX,
                                                  radius=radius,
                                                  limit_shape=(Param.Z_DIM,Param.Y_DIM,Param.X_DIM),
                                                  xy_vs_z=Param.XY_VS_Z)  #[Z,Y,X] 
        #self.perim = trace_perimeter_coord(self.l_pixels)                               #[Z,Y,X] 

        self.flag_exterior_frame = self.check_exterior_frame_flag()
        self.overlap_max = 0
        
        self.d_clusterctr_overlap = {}  
        
    @staticmethod
    def init_class_variables():
        Sphere.ctr = 0
        return
        
    def check_exterior_frame_flag(self):
        '''
        checks if the sphere overlaps with the frame of the picture.  If so the sphere will be marked as a exterior sphere
        if you use an exterior_mask : this will always be false
        '''
        
        #check if it touches the exterior frame
        if ((self.centre[1] - self.radius < 0) or 
            (self.centre[2] - self.radius < 0) or 
            (self.centre[0] - (self.radius * Param.PIXEL_WIDTH/Param.VOXEL_DEPTH) < 0) or 
            (self.centre[1] + self.radius > Stack.curr_stack.img.shape[1]) or
            (self.centre[2] + self.radius > Stack.curr_stack.img.shape[2]) or
            (self.centre[0] + (self.radius * Param.PIXEL_WIDTH/Param.VOXEL_DEPTH) > Stack.curr_stack.img.shape[0])
            ):
            return True
     
        return False

    def repr(self):
        print('Info on Sphere (ctr={0})-->'.format(self.ctr))
        print('{0}centre={1}'.format('\t',self.centre))
        print('{0}radius={1}'.format('\t',self.radius))
        print('{0}cluster={1}'.format('\t',self.cluster))
        print('{0}stack ={1}'.format('\t',Stack.curr_stack))

        
        
    @classmethod
    def how_many(cls):
        print("We have {:d} spheres.".format(cls.ctr))
       
        
class Cluster:
    
    CTR_BACKGROUND=0
    CTR_MASK = 1 #do not use negative ctr's, will mess up bincount
    ctr = 1 #the mask already occupies 1  

    dcopy_overlap = {'overlap_pix':['Z','Y','X'],      #copybook that will keep track of overlap info with a certain cluster (cfr value in d_clusterctr_overlap)
                     'overlap_max_DT':0       #the maximum of the DT in the overlap zone 
                     }
    dcopy_overlap_unk = {'overlap_pix':[[],[],[]],
                         'overlap_max_DT':0
                         }
    
    dcopy_label_info = {'ctr':'',
                        'color':'',
                        'cell':None
                        }
        
    d_reason_code_no_division_desc  = {'':'not considered',
                                       '0':'considered',
                                       '1':'no match in prev stack',
                                       '2':'not large enough',
                                       '3':'too much overlap with seed',
                                       '4':'radii of dividing cells too different'
                                        }
    d_ctr_cluster = {}
    
    l_curr_stack = []  #just used for the labelling of the clusters.  this will be ordered in decreasing size and can be processed sequentially
    
    def __init__(self, sphere, stack, seeding_cell=None,GT_seedZYX=None):
        self.label_history = [] # list of tuples (label, history_trace)
        
        if seeding_cell:
            self.seeding_cell_label = seeding_cell.label #only keeping label will increase performance
        else:
            self.seeding_cell_label = None
        #self.ic_use_for_seeding=False #will be determined later on, after labelling clustering (when creating cells)
        
        Cluster.ctr += 1
        self.ctr = Cluster.ctr
        
        
        self.spheres = [sphere] 
        self.centre = sphere.centre
        self.l_pixels = []   #[Z,Y,X]
        self.radius = 0 #this will be the radius calculated from the pixel volume (so if the cluster would be a perfect sphere what would be it's radius)
        self.volume = 0 #expressed in nb pixels
        self.stack = stack #
        
        self.validation_name = ''
        self.validation_centre = ''
        self.validation_error = ''

        
        self.d_clusterctr_overlap = {}   #key : a cluster counter   ;  value :dcopy_overlap (=the overlap info with that cluster )
        
        self.candidate_new_cell=False #marks a possible, unconfirmed, cluster division
        self.seed_new_cell = False #marks a confirmed cluster division
        self.seed_divided = False #marks a confirmed cluster division
        self.reason_cd_no_division=''  #used for tracking cell division, for information only (cfr d_reason_code_no_division_desc )
        
        self.cell = None
        self.cell_ctr_CTC = 0  #this is a cell counter that's
        self.label = None
        
        sphere.cluster=self
        Cluster.d_ctr_cluster[self.ctr] = self
        Cluster.l_curr_stack.append(self)
        
        self.flag_exterior_frame=sphere.flag_exterior_frame
        
        self.GT_seeded = True if GT_seedZYX else False
        
        
    @staticmethod
    def init_class_variables():
        
        Cluster.ctr = 1 #the mask already occupies 1  
        Cluster.dcopy_overlap = {'overlap_pix':['Z','Y','X'],      #copybook that will keep track of overlap info with a certain cluster (cfr value in d_clusterctr_overlap)
                     'overlap_max_DT':0       #the maximum of the DT in the overlap zone 
                     }
        Cluster.d_ctr_cluster = {}
        Cluster.l_curr_stack = [] 
        
        return
    
    @staticmethod
    def print_stats(info=""):
        nb_total = 0
        nb_GT_seeded = 0
        nb_seeded = 0
        nb_filtered = 0
        
        for cluster_i in Cluster.l_curr_stack:
            nb_total += 1
            if cluster_i.label == Label.FILTERED:nb_filtered += 1
            if cluster_i.GT_seeded:nb_GT_seeded += 1
            if cluster_i.seeding_cell_label: nb_seeded += 1
            
        print('___Cluster stats ({4}) : {0} TOTAL, of which {1} GT_seeded, {2} seeded, {3} filtered (current stack)'.format(nb_total,nb_GT_seeded,nb_seeded,nb_filtered,info))
        
        return 
        
    def get_radius(self):
        if self.radius==0:
            self.radius = (self.get_volume() * (3/4*pi))**(1/3)
        
        return self.radius
            
    def get_seed_radius(self):
        return self.spheres[0].radius


    def get_cell(self):
        if not self.cell:
            self.cell = self.stack.d_label_cell.get(self.label)
        
        return self.cell
           
    def get_volume(self):
        if self.volume==0:
            if not self.l_pixels:
                self.gather_pixels()
            self.volume=len(self.l_pixels[0])
        
        return self.volume
    
    def get_l_pixels(self):
        
        if not self.l_pixels:
            self.gather_pixels()
        
        return self.l_pixels
    
  
        
    def add_sphere(self,sphere):
        self.spheres.append(sphere)
        if sphere.flag_exterior_frame:
            self.flag_exterior_frame = True
        sphere.cluster=self

        
    #        
    #     def delete_cluster(self):
    #         del self.stack.d_number_cluster[self.ctr]
    #         Cluster.ctr -= 1
         
    def merge_with_cluster(self,merge_cluster):
        
        for sphere in self.spheres: 
            merge_cluster.spheres.append(sphere)
            sphere.cluster = merge_cluster
        
        if sphere.flag_exterior_frame: merge_cluster.flag_exterior_frame=True
        
        for stack in Stack.curr_stack.stack:
            stack.clusterview[stack.clusterview==self.ctr] = merge_cluster.ctr
            stacks:stack.remove_cluster(self)
        
        print('cluster', self.ctr, ' merged with cluster ', merge_cluster.ctr, 'now delete it') 
          
        self.delete_cluster()
        
        return
     

    def gather_pixels(self,hide_exterior=True): 
        temp_img = np.zeros(self.stack.img.shape,'bool')
        for sphere in self.spheres: temp_img[tuple(sphere.l_pixels)] = True
        
        if hide_exterior:
            temp_img[self.stack.exterior_mask==False]=False  #cut of exterior mask
            self.l_pixels = np.where(temp_img==True)
        return 
        
    def is_seed(self):
        if not self.cell:
            return False
        
        if self.cell.cluster_seed==self:
            return True
        
        return False
    
    def is_largest_non_seed_cluster(self):
        
        if not self.cell:
            return False
        
        if self.cell.cluster_largest_unseeded==self:
            return True
        
        return False
    
 
    
    def update_overlap_pixels(self,cluster_ctr_overlap,d_overlap_sphere):
        '''
        updates the overlap info for the cluster(self) with another cluster (cluster_ctr_overlap) given the overlap info of a new sphere
        :param cluster_ctr_overlap: the cluster counter for which the overlap info must be updated
        :param d_overlap_sphere: the overlap info of a sphere that will be integrated into the overlap info of the cluster (self)
        '''
        cluster_1 = self
        cluster_2 = Cluster.d_ctr_cluster[cluster_ctr_overlap]  #symmetric relationship
        
        
        d_overlap_cluster= cluster_1.d_clusterctr_overlap.get(cluster_2)  #retrieve the overlap info for cluster_ctr_overlap
        if d_overlap_cluster:  #appending the sphere l_pixels to the cluster
            Z,Y,X = d_overlap_cluster['overlap_pix'] 
            Z_append,Y_append,X_append = d_overlap_sphere['overlap_pix']
            Z=np.append(Z,Z_append);Y=np.append(Y,Y_append);X=np.append(X,X_append)
            d_overlap_cluster['overlap_pix'] = [Z,Y,X]
            
            if d_overlap_sphere['overlap_max_DT'] > d_overlap_cluster['overlap_max_DT']:
                d_overlap_cluster['overlap_max_DT'] = d_overlap_sphere['overlap_max_DT'] 
                
            cluster_1.d_clusterctr_overlap[cluster_2.ctr] = d_overlap_cluster
            cluster_2.d_clusterctr_overlap[cluster_1.ctr] = d_overlap_cluster
                
        else:
            cluster_1.d_clusterctr_overlap[cluster_2.ctr]= d_overlap_sphere
            cluster_2.d_clusterctr_overlap[cluster_1.ctr]= d_overlap_sphere
  
        return 
    
    def overlaps_too_much_with_mask(self,stack):
        
        seed_sphere  = self.spheres[0]
        if len(seed_sphere.l_pixels[0])==0:
            return False
        
        nb_overlap_pix_mask = stack.get_nb_overlap_pix_with_mask(seed_sphere)
        if nb_overlap_pix_mask/len(seed_sphere.l_pixels[0]) > Param.FILTER_SPHERE_OVERLAP_MASK_PC:
            return True

        return False
   
    def collides_exterior_frame(self):
        for sphere in self.spheres:
            if sphere.flag_exterior_frame:return True
        return False
    
    def set_new_label(self,label_parent):
        '''
        a cluster gets a new label = a number/counter
        the label_parent can be forced, or if none given, it will be set to Label.PARENT_UNKNOWN_CTR
        :param label_parent: the label_parent assigned
        '''

        label_new = Label(timestep=self.stack.timestep,label_parent=label_parent)
        
        self.set_label(label_new.ctr, trace = 'NEW',label_parent=label_parent)
        return  
    
    def set_val_info(self,val_name, val_centre, val_error):
        
        self.validation_name=val_name
        self.validation_centre=val_centre
        self.validation_error=val_error
        
        if self.cell:
            self.cell.set_val_info(val_name, val_centre, val_error)
        
        return
    
    def get_color(self):
        
        return self.label.color if self.label else 'NoColor'
    
    def get_label(self):
        return self.cell.label.ctr
        
    def set_label(self,label,trace=''):
        self.label = label
        self.label_history.append((label.ctr,trace))  
        
    #         if self.label in Label.L_INDICATOR:
    #             self.ic_use_for_seeding=False
            

    #     def set_candidate_new_cluster(self,TrueOrFalse):  
    #         
    #         if self.label in  [Label.EXTERIOR_CTR,Label.VALIDATION]:
    #             print('Warning : an exterior or validation cluster was ordered to become a candidate new cluster. request rejected')
    #             return
    #         
    #         self.candidate_new_cell=TrueOrFalse
    #         if TrueOrFalse:
    #             self.seeding_power=0
    #         else:
    #             self.seeding_power=1      
    #         
    #         return
            
    def repr(self,output=sys.stdout):
        
        print('Info on cluster (ctr={0},label={1})'.format(self.ctr,self.label),end='',file=output)
        print('cluster centre at ', self.centre,file=output)
        if self.label:
            print('the label is {0}'.format(self.label.color),file=output)
        print('it is validated by {0} with centre {1}'.format(self.validation_name,self.validation_centre),file=output)
        print('this cluster is built with {0} spheres'.format(len(self.spheres)),file=output)
        print('this cluster has radius {0}'.format(self.get_radius()),file=output)
        print('{0}spheres (ctr,centre,radius) ='.format('-->'),['{0}_{1}_{2}|'.format(i.ctr,i.centre,i.radius) for i in self.spheres],file=output)
        print('candidate_new_cell=',self.candidate_new_cell,file=output)
        print('seed_new_cell=',self.seed_new_cell,file=output)
        print('seed_divided=',self.seed_divided,file=output)
        print('reason_cd_no_division=',self.reason_cd_no_division,file=output)
        print('GT_seeded=',self.GT_seeded,file=output)
        if self.cell:
            if self.cell.label:
                print('cell=',self.cell," with label ",self.cell.label.ctr,file=output)
            else:
                print('cell=',self.cell, ' with no label attached to the cell',file=output)

        for cluster_overlap_i, d_overlap_max_dt_i in self.d_clusterctr_overlap.items():
            print('overlaps with cluster {0} with overlap max_dt of {1}'.format(cluster_overlap_i,d_overlap_max_dt_i['overlap_max_DT']),file=output)
        
        print('',file=output)
        
        return output
               
        

class Cell:

    def __init__(self,label,cluster_seed,stack):
    
        self.stack = stack
        
        self.cluster_seed = cluster_seed
        self.cluster_largest_unseeded = None
        self.l_clusters_other = []
        
        cluster_seed.cell = self
        
        self.set_label(label) #obligatory to first have a label
        
        self.cell4D = self.label.cell4D #for cell births this is only consistent after cell4D creation (= when cell4D links to a label)
        
        self.stack.cells.append(self)
        
    
        self.stack.d_label_cell[self.label] = self
        #self.stack.add_cell_to_label_index(self)
        
        self.radius = 0 #radius of pixel volume
        self.volume = 0 #pixel volume
        
        self.validation_name = ''
        self.validation_centre = ''
        self.validation_error = ''
       
        
    def set_label(self,label):
        for cluster_i in self.get_clusters():
            cluster_i.label = label
            
        self.label = label
        
        #         if label in Label.L_INDICATOR:
        #             self.cluster_seed.ic_use_for_seeding=False
            
        return
    
    def set_val_info(self,val_name,val_centre,val_error):
        
        self.validation_name=val_name
        self.validation_centre=val_centre
        self.validation_error=val_error
        
        if self.cell4D:
            self.cell4D.set_val_info(val_name,val_centre,val_error)  # will always overwrite with the validation of latest stack
        
        return
        
    
    def repr(self):
        print('Info on Cell(name_lineage={0})-->'.format(self.name_lineage))
        if self.label:
            print('{0}label={1}'.format('\t',self.label.color))
        else:
            print('{0}has no label'.format('\t'))
        print('{0}cluster_seed_ctr={1} (GT_seeded={2})'.format('\t',self.cluster_seed.ctr,self.cluster_seed.GT_seeded))
        if self.cluster_largest_unseeded:
            print('{0}cluster_largest_unseeded_ctr={1} (GT_seeded={2})'.format('\t',self.cluster_largest_unseeded.ctr,self.cluster_largest_unseeded.GT_seeded))
        else:
            print('{0}has no largest unseeded cluster'.format('\t'))
        if self.l_clusters_other:
            print('{0}other clusters={1}'.format('\t',[i.ctr for i in self.l_clusters_other]))
        else:
            print('{0}has no other clusters'.format('\t'))
        print('{0}radius={1}'.format('\t',self.get_radius()))
        print('{0}volume={1}'.format('\t',self.get_volume()))
        print('{0}seeding_volume={1}'.format('\t',self.get_seeding_volume()))
        print('{0}stack ={1}'.format('\t',self.stack))

        return
    
    def get_clusters(self,mode=''):
        l_clusters=[]
        
        if mode!='no_cluster_seed':
            if self.cluster_seed:
                l_clusters.append(self.cluster_seed)
            else:
                if DBG:print('inconsistency : no seed cluster for cell {0}!'.format(self.ctr))
                
        if self.cluster_largest_unseeded:
            l_clusters.append(self.cluster_largest_unseeded)
        else:
            if DBG and self.l_clusters_other:print('inconsistency :cluster_largest_unseeded is None, but other cells are filled in ! {0}!'.format(self.ctr))
        if self.l_clusters_other:
            for cluster_i in self.l_clusters_other:
                l_clusters.append(cluster_i)
            
        return l_clusters
    
    def add_cluster(self,cluster):
        
        
        if not self.cluster_largest_unseeded:
            self.cluster_largest_unseeded = cluster
            
        else:    
            #if cluster.get_volume() > self.cluster_largest_unseeded.get_volume():
            if cluster.get_seed_radius() > self.cluster_largest_unseeded.get_seed_radius():
                self.l_clusters_other.append(self.cluster_largest_unseeded)
                self.cluster_largest_unseeded = cluster
            else:
                self.l_clusters_other.append(cluster)
                
        self.stack.labelview[cluster.get_l_pixels()] = self.label.ctr
        cluster.label = self.label
        cluster.cell = self
        #cluster.ic_use_for_seeding=False
    
        return
    
    def gather_pixels(self):
            
        temp_img = np.zeros(self.stack.img.shape,'bool')
        for cluster_i in self.get_clusters(): 
            temp_img[cluster_i.get_l_pixels()] = True
        temp_img[self.stack.exterior_mask==False]=False  #cut of exterior mask
        l_pixels = np.where(temp_img==True)
        
        
        self.volume = len(l_pixels[0])
        self.radius = (self.volume * (3/4*pi))**(1/3)
        
        return 
        
    def replace_label(self,new_label,trace=''):
        old_label = self.label
        self.label = new_label
        for cluster_i in self.get_clusters(): 
            cluster_i.set_label(new_label,trace=trace)
            
        self.stack.d_label_cell[self.label] = self
        del self.stack.d_label_cell[old_label]

        return
        
       
    def get_seeding_volume(self,func="disk"):
        '''
        get the seeding volume for a cell that will be used for seeding
        The seeding sphere is used.  this will limit the range of the cell as a whole
        update : a seeding disk is used, constrain the z movement by seeding_z_range (up or down)
        '''
        #l_pixels_sphere = self.cluster_seed.spheres[0].l_pixels  Z,Y,X
     
        if func=="disk":
            a_pixels_sphere = np.array(self.cluster_seed.spheres[0].l_pixels)
            z_centre = self.cluster_seed.centre[0]
            l_pixels_disk = list(a_pixels_sphere[:,np.where(abs(a_pixels_sphere[0]-z_centre) <Param.SEEDING_Z_RANGE + 1)[0]])
        elif func=="min_seed_disk":  ##still a disk , but x,y range is limited to minimum seed sphere
            a_pixels_sphere = np.array(get_3D_sphere_coordinates(a_coordZYX=self.cluster_seed.centre,
                                                  radius=Param.MIN_SEED_RADIUS,
                                                  limit_shape=(Param.Z_DIM,Param.Y_DIM,Param.X_DIM),
                                                  xy_vs_z=Param.XY_VS_Z))  #get_3D_sphere_coordinates returns tuple (Z,Y,X)
            
            z_centre = self.cluster_seed.centre[0]
            l_pixels_disk = list(a_pixels_sphere[:,np.where(abs(a_pixels_sphere[0]-z_centre) <2)[0]])
            
        
        
        return l_pixels_disk
    
    def get_color(self):       
        return self.label.color if self.label else 'NoColor'
    
    def get_volume(self):
        """gives back volume in nb of pixels"""
        if self.volume == 0:
            self.gather_pixels()
        
        return self.volume 

    def get_volume_um(self):
        """gives back volume in um3"""
        if self.volume == 0:
            self.gather_pixels()
        
        return self.volume * Param.PIXEL_WIDTH * Param.VOXEL_DEPTH * Param.PIXEL_WIDTH 
    
    def get_radius(self):
        
        if self.radius == 0:
            self.gather_pixels()
        
        return self.radius
        
    
    def check_cell_division(self):
        '''
        check if a cell is eligible for cell division.  This is mostly based on information of the current stack (=looking at overlap),but also the previous stack
        (check if a matching cluster can be found)
        the reason is tracked for ease of debugging
        '''
        
        verbose=False
        if not self.cluster_largest_unseeded:
            return False #only 1 cluster, so no cell division possible
        
        self.cluster_largest_unseeded.reason_cd_no_division = '0' # consider this cluster for cell division
        

                   
        
        cluster_match_prev_stack = Stack.prev_stack.get_matching_cluster_from_stack_based_on_max_overlap(self.cluster_largest_unseeded)
        if not cluster_match_prev_stack:
            self.cluster_largest_unseeded.reason_cd_no_division = '1' 
            return False #if the possible daughter cell does not have a matching cluster in the previous stack (=unlikely), no cell division should take place
        
        if self.cluster_largest_unseeded.get_radius() < Param.MIN_CELL_RADIUS: 
            self.cluster_largest_unseeded.reason_cd_no_division = '2' 
            if verbose:print('temp : the largest non seed cluster {2} is not large enough to be considered a new seed cluster : {0} < {1}'.format(self.cluster_largest_unseeded.get_radius(),Param.MIN_CELL_RADIUS,self.cluster_largest_unseeded.ctr))
            return False #the largest non seed cluster is not large enough to be considered a new seed cluster## MIN_CELL_RADIUS = minium cell size of a new cell
        
        
        if Param.THRESHOLD_CELL_OVERLAP_TOO_MUCH_WITH_MASK:
            if self.cluster_largest_unseeded.overlaps_too_much_with_mask(self.stack):
                if verbose:print('temp : cluster_largest_unseeded.overlaps_too_much_with_mask=TRUE')
                return False #a cluster that overlaps too heavily with the exterior should never become a cell
        
        
        d_overlap= self.cluster_seed.d_clusterctr_overlap.get(self.cluster_largest_unseeded.ctr)
        if not d_overlap:
            return True  #if the 2 clusters have no overlap then off course they are candidate for cell division
        
        if self.cluster_largest_unseeded.spheres[0].radius < self.cluster_seed.spheres[0].radius:
            threshold_dt = self.stack.get_overlap_maxDT_threshold(self.cluster_largest_unseeded.spheres[0])
        else:
            threshold_dt = self.stack.get_overlap_maxDT_threshold(self.cluster_seed.spheres[0])
        # min_radius = min(self.cluster_largest_unseeded.spheres[0].radius,self.cluster_seed.spheres[0].radius)
        # if d_overlap['overlap_max_DT'] > self.stack.get_overlap_maxDT_threshold(min_radius): 
        if d_overlap['overlap_max_DT'] > threshold_dt:
            self.cluster_largest_unseeded.reason_cd_no_division = '3' 
            if verbose:print('the largest non seed cluster {0} overlaps too much with the seed : {1} > {2}'.format(self.cluster_largest_unseeded.ctr,d_overlap['overlap_max_DT'],threshold_dt))
            return False #the seed cluster and largest non seed cluster have too much overlap to be a candidate for cell division
        

        if self.cluster_seed.get_radius() > self.cluster_largest_unseeded.get_radius():
            radius_ratio = self.cluster_largest_unseeded.get_radius() /  self.cluster_seed.get_radius()
        else:
            radius_ratio = self.cluster_seed.get_radius() / self.cluster_largest_unseeded.get_radius() 
        if verbose:print('the largest non seed cluster {0} has ratio {1}<{2} ?radius seed={3} //radius daughter= {4} ('.format(self.cluster_largest_unseeded.ctr,radius_ratio,Param.THRESHOLD_CELL_DIVISION_MIN_RADIUS_RATIO,self.cluster_seed.get_radius(),self.cluster_largest_unseeded.get_radius()))
        if radius_ratio < Param.THRESHOLD_CELL_DIVISION_MIN_RADIUS_RATIO:
            if verbose:print('the largest non seed cluster {0} is not equal in size of seed :ratio {1}<{2}'.format(self.cluster_largest_unseeded.ctr,radius_ratio,Param.THRESHOLD_CELL_DIVISION_MIN_RADIUS_RATIO))
            self.cluster_largest_unseeded.reason_cd_no_division = '4' 
            return False #the 2 dividing clusters should have roughly the same volume


        
        if verbose:
            print('check_cell_division =TRUE:the largest non seed cluster {0} does not overlap too much with the seed : {1} < {2}.(radius = {3})'.
                format(self.cluster_largest_unseeded.ctr, d_overlap['overlap_max_DT'], threshold_dt,self.cluster_largest_unseeded.get_radius()))
        

        return True
        
    def execute_cell_division(self,label_forced=None,true_division_time=True):
        '''
        the largest non seed cluster will be marked as a condidate new Cell.  
        If this cluster already have that label in the previous stack, then cluster division will take place
        This involves creating a new Cell for the largest non seed cluster, and also dividing the other clusters of the Cell based on maximum overlap.  
        as a fallback they will be assigned to the original Cell
        
        for the current stack, we force the label to be the same as in the previous stack by using th label_forced
        
        the label view will be updated accordingly
        
        true_division_time indicates if this timestep is the actual time of cell division.  because of the confirmation process a cell is divided in both previous stack and current stack
        but only the previous stack is the true time of division
        '''
        def check_merge_with_new_cell():
            overlap_info_w_new_seed= cluster_i.d_clusterctr_overlap.get(cell_new.cluster_seed.ctr)
            if not overlap_info_w_new_seed:
                return False
            overlap_pix_w_new_seed=overlap_info_w_new_seed['overlap_pix']
            
            overlap_info_w_old_seed= cluster_i.d_clusterctr_overlap.get(self.cluster_seed.ctr)
            if not overlap_info_w_old_seed:
                return True
            overlap_pix_w_old_seed=overlap_info_w_old_seed['overlap_pix']
            
            if len(overlap_pix_w_new_seed) > len(overlap_pix_w_old_seed):
                return True
            
            return False
        
        #creating the new cell
        cell_new = Cell(label=label_forced,
                        cluster_seed=self.cluster_largest_unseeded,
                        stack=self.stack)
        
        #removing the largest unseeded cluster of the original cell, it has promoted to a new cell (=becomes a seed)
        if true_division_time:
            self.cluster_largest_unseeded.seed_new_cell=True
            self.cluster_seed.seed_divided=True
        self.cluster_largest_unseeded.candidate_new_cell=False
        self.cluster_largest_unseeded=None
        
        
        #redividing the other clusters 
        l_clusters_to_redivide = self.l_clusters_other
        if l_clusters_to_redivide:
            self.l_clusters_other=[]
            for cluster_i in l_clusters_to_redivide:
                if check_merge_with_new_cell():
                    cell_new.add_cluster(cluster_i)
                else:
                    self.add_cluster(cluster_i)

        self.gather_pixels() # recalculates new volume
                
        return cell_new

        
class Cell4D:
    ctr = 0
    l_cell4D = [] #not deleted over timestep (only wise if cell4D is not connected to other objects)
    set_validation_names_assigned = set() # set of validation names that are assigned to cell4D objects. If every cell division changes the names of the 2 daughtercells, this does not need initialization each stack, if 1 of the daughters inherits the mothers name, this should be cleaned up every timestep
    
    
    def __init__(self,initiating_cell,timestep):
        self.label = initiating_cell.label #the label plays the role of a cell4D counter
        
        self.create_timestep = timestep #timestep when the object is created
        if timestep==1:
            self.birth_timestep= timestep #
            self.lifetime = 0
        else:
            self.birth_timestep= timestep - 1 #cell confirmation lags 1 timestep
            self.lifetime = 1
            
        self.max_volume = initiating_cell.get_volume()
        self.max_radius = initiating_cell.get_radius()
        initiating_cell.cell4D = self
        #self.initiating_cell = initiating_cell  #bad idea ! keep Cell4D light, because preserved across the timelapse
        self.initiating_cell_ctr = initiating_cell.cluster_seed.ctr

        Cell4D.l_cell4D.append(self)

        
        self.label.set_cell4D(self)
        
        self.l_validation_names=[]
        self.validation_centre=''
        self.validation_error=''
        

    @staticmethod
    def init_class_variables():
        Cell4D.ctr = 0
        Cell4D.l_cell4D = []
        Cell4D.set_validation_names_assigned = set() 
        
        return
    
    def incorporate_cell_data(self,cell):
        cell4D = cell.cell4D
        if cell.cell4D.max_volume < cell.get_volume(): 
            cell4D.max_volume = cell.get_volume()
            cell4D.max_radius = cell.get_radius()

    
    #     @staticmethod
    #     def age():
    #         for cell4D_i in Cell4D.l_cell4D:
    #             cell4D_i.lifetime +=1
    #             
            
    def has_died_in_this_timestep(self):
        
        if not Stack.prev_stack:
            return False #a cell cannot die in the first stack
        
        if Stack.prev_stack.d_label_cell.get(self.label) and not Stack.curr_stack.d_label_cell.get(self.label):
            return True  #the label was present in the previous stack, but not present in the current stack anymore
        
        return False
    
    
    def replace_label(self,new_label,trace=''):
        
        cell_curr_stack = Stack.curr_stack.d_label_cell.get(self.label)
        if cell_curr_stack:
            cell_curr_stack.replace_label(new_label,trace=trace)
            
        cell_prev_stack = Stack.prev_stack.d_label_cell.get(self.label)
        if cell_prev_stack:
            cell_prev_stack.replace_label(new_label,trace=trace)
        
        self.label=new_label
        return
        
    def get_color(self):       
        return self.label.color if self.label else 'NoColor'
    
    def repr(self,output=sys.stdout):
        print('Info on cell4D (label={0})-->'.format(self.label),file=output)
        print('{0}label={1}'.format('\t',self.label.color),file=output)
        print('{0}birth_timestep={1}'.format('\t',self.birth_timestep),file=output)
        print('{0}lifetime={1}'.format('\t',self.lifetime),file=output)
        print('{0}max_volume={1}'.format('\t',self.max_volume),file=output)       

    
    def set_val_info(self,val_name,val_centre,val_error,overwrite_info=False):
        
    #         if not overwrite_info and self.validation_name:
    #             return  #initial assigned info will not be overwritten 
        
        def check_validation_name():
            
            def other_daughter_cell_is_already_linked_to_cell4D():
                
                if not(len(val_name) > 1 and val_name[-2]=='d'):
                    return False # the validation cell is not a daughter cell to begin with, so no problem
                    
                s_stem_name = val_name[:-1]
                for validation_name_i in self.l_validation_names:
                    if validation_name_i[:-1] == s_stem_name:
                        return True
                return False
            
            if val_name in Cell4D.set_validation_names_assigned:
                return False # once a validation cell is assigned to a Cell4D object, it sticks, cannot be changed
            
    #             if other_daughter_cell_is_already_linked_to_cell4D():
    #                 return False # a cell4D cannot be assigned to 2 daughter cells at the same time
                                  
            return True
        
        
        if not val_name:  #if val_name is blank, just update the error info
            self.validation_error=val_error
            return
            
        if check_validation_name():
            self.l_validation_names.append(val_name)
            Cell4D.set_validation_names_assigned.add(val_name)
            self.validation_centre=val_centre
            self.validation_error=val_error
            
        return        
        
        
class Label():
    UNASSIGNED_CTR = 0
    EXTERIOR_CTR=1
    FILTERED_CTR = 9998
    VALIDATION_CTR = 9999
    
    UNASSIGNED = None
    EXTERIOR=None
    FILTERED=None
    VALIDATION=None
    L_INDICATOR=[]

    
    PARENT_UNKNOWN_CTR=0
    PARENT_ZEUS_CTR = 99
    
    D_LABEL_COLOR = {1:'salmon',2:'lime',3:'light pink',4:'magenta',5:'yellow',6:'light sky blue',7:'cyan',
                 8:'blue',9:'olive',10:'purple',11:'green',12:'orange',13:'maroon',14:'blue violet',15:'chocolate',
                 16:'teal',17:'navy',18:'beige',19:'dark sea green',20:'gold'}
    
    
    L_MAIN_COLORS = ['orange','yellow','green','skyblue','blue','purple','fuchia','magenta','gray']
    COLOR_EXTERIOR = 'white'
    COLOR_VALIDATION= 'red'
    COLOR_UNASSIGNED = 'white'
    COLOR_OVERLAP = 'black'
    
    COLOR_FILTERED_SPHERE = COLOR_EXTERIOR
    COLOR_FILTERED_SPHERE_CENTRE_EXTERIOR = 'red'
    COLOR_FILTERED_SPHERE_CENTRE_MIN_SEED_RADIUS = 'blue' #caveat : in the first stack MIN_CELL_RADIUS instead of MIN_SEED_RADIUS 
    COLOR_FILTERED_SPHERE_CENTRE_OVERLAP_MASK = 'green'
    COLOR_FILTERED_SPHERE_CENTRE_Z_TOP_BOTTOM = 'yellow'
  
    
    #https://www.rapidtables.com/web/color/RGB_Color.html
    D_COLOR_RGBA = {'orange':(255,128,0,255),
                    'orange+':(204,102,0,255),
                    'orange++':(153,76,0,255),
                    'orange+++':(102,50,0,255),
                    'orange-':(255,153,51,255),
                    'orange--':(255,178,102,255),
                    'orange---':(255,204,153,255),
                    
                    'yellow':(255,255,0,255),
                    'yellow+':(204,204,0,255),
                    'yellow++':(153,153,0,255),
                    'yellow+++':(102,102,0,255),
                    'yellow-':(255,255,51,255),
                    'yellow--':(255,255,102,255),
                    'yellow---':(255,255,153,255),
                    
                    'green':(0,255,0,255),
                    'green+':(0,204,0,255),
                    'green++':(0,153,0,255),
                    'green+++':(0,102,0,255),
                    'green-':(51,255,51,255),
                    'green--':(102,255,102,255),
                    'green---':(153,255,153,255),
                    
                    'skyblue': (0,255,255,255),
                    'skyblue+': (0,204,204,255),
                    'skyblue++': (0,153,153,255),
                    'skyblue+++': (0,102,102,255),
                    'skyblue-': (51,255,255,255),
                    'skyblue--': (102,255,255,255),
                    'skyblue---': (153,255,255,255),
                    
                    'blue':(0,0,255,255),     
                    'blue+':(0,0,204,255),
                    'blue++':(0,0,153,255),
                    'blue+++':(0,0,102,255),
                    'blue-':(51,51,255,255),
                    'blue--':(102,102,255,255),
                    'blue---':(153,153,255,255),       
                   
                   'purple':(127,0,255,255),
                   'purple+':(102,0,204,255),
                   'purple++':(76,0,153,255),
                   'purple+++':(51,0,102,255),
                   'purple-':(153,51,255,255),
                   'purple--':(178,102,255,255),
                   'purple---':(204,153,255,255),
                   
                   'fuchia':(255,0,255,255),
                   'fuchia+':(204,0,204,255),
                   'fuchia++':(153,0,153,255),
                   'fuchia+++':(102,0,102,255),
                   'fuchia-':(255,51,255,255),
                   'fuchia--':(255,102,255,255),
                   'fuchia---':(255,153,255,255),
                   
                   'magenta':(255,0,127,255),
                   'magenta+':(204,0,102,255),
                   'magenta++':(153,0,76,255),
                   'magenta+++':(102,0,51,255),
                   'magenta-':(255,51,153,255),
                   'magenta--':(255,102,178,255),
                   'magenta---':(255,153,204,255),
                   
                   'gray':(128,128,128,255),
                   'gray+':(96,96,96,255),
                   'gray++':(64,64,64,255),
                   'gray+++':(32,32,32,255),
                   'gray-':(160,160,160,255),
                   'gray--':(192,192,192,255),
                   'gray---':(224,224,224,255),
                   
                   'lime':(0,255,0,255),
                   'cyan':(0,204,204,255),
                   'black':(0,0,0,255),
                   'black_light':(50,50,50,255),
                   'white':(255,255,255,255),
                   'off-white':(230,230,230,255),
                   'red':(255,0,0,255),
                   'red vague':(255,0,0,15),
                   'silver':(192,192,192,255),
                   'gray':(128,128,128,255),
                   'maroon':(128,0,0,255),
                   'olive':(128,128,0,255),
                   'teal':(0,128,128,255),
                   'navy':(0,0,128,255),
                   'beige':(245,245,220,255),
                   'light pink':(255,182,193,255),
                   'deep pink':(255,20,147,255),
                   'blue violet':(138,43,226,255),
                   'chocolate':(210,105,30,255),
                   'salmon':(250,128,114,255),
                   'gold':(255,215,0,255),
                   'dark sea green':(143,188,143,255),
                   'plum':(221,160,221,255),
                   'lawn green':(124,252,0,255)
                   }
    
    ctr=0
    d_ctr_label = {}
    
    l_labels_ctr = [] #not deleted over timestep (only wise if Label is not connected to other objects)
    
    d_name_lineage_ctr = {} #this will only be used to give a unique number to every name_lineage, used when outputtin cluster_summary
    
    def __init__(self,timestep,label_parent_ctr='',cell=None,cell4D=None,ctr=''):
        '''
        A label has a counter, which is just a unique ID
        A label is associated with a certain color, if colors run out, more than 1 label can have the same color (only impacts visualization)
        A label is linked to a parent_label through the counter of the parent_label
        A label is mostly linked 1-to-1 to a cell4D.  
               Still, a label is not exactly a cell4D, because of the technical aspect of building up cells. a Cell4D can only be created at the very end, while labels must already by 
               assigned to clusters in the beginning
               also the exterior has a label, and validation cells.. these are not Cell4D
        A label is linked to the name lineage, but this lineage is extended over time ! This is because the seed cell is retained after cell division (instead of having 2 new daughtercells) 
        to track the log you need to check the cluster_summary.
        '''
        if ctr:
            self.ctr=ctr 
        else:
            Label.ctr +=1
            self.ctr = Label.ctr
            
        self.color = self.__get_color_for_label()
        
        self.label_parent_ctr = label_parent_ctr if label_parent_ctr else Label.PARENT_UNKNOWN_CTR
        
        self.birth_timestep= timestep
        
        self.cell4D = cell4D  #is light, so ok

        Label.d_ctr_label[self.ctr] = self
        Label.l_labels_ctr.append(self.ctr)
        
        if label_parent_ctr in ['',Label.PARENT_ZEUS_CTR]:
            self.name_lineage = str(self.ctr)
        else:
            self.name_lineage = '' # the lineagenames after cell division will only be filled in before outputting to csv
        
    @staticmethod
    def init_class_variables():
        Label.ctr=0
        Label.d_ctr_label = {}
        Label.l_labels_ctr = [] 
        return   
     
    # def update_lineage_name_cell_division(self):
        # self.name_lineage += 'd1'
         
        
        # return
        

    
    @staticmethod
    def create_indicator_labels():
        Label.UNASSIGNED = Label(timestep=0,label_parent_ctr='',cell=None,cell4D=None,ctr=Label.UNASSIGNED_CTR)
        Label.EXTERIOR= Label(timestep=0,label_parent_ctr='',cell=None,cell4D=None,ctr=Label.EXTERIOR_CTR)
        Label.FILTERED = Label(timestep=0,label_parent_ctr='',cell=None,cell4D=None,ctr=Label.FILTERED_CTR)
        Label.VALIDATION = Label(timestep=0,label_parent_ctr='',cell=None,cell4D=None,ctr=Label.VALIDATION_CTR)
        Label.L_INDICATOR = [Label.VALIDATION,Label.FILTERED,Label.EXTERIOR,Label.UNASSIGNED]
        
        return
        
    def is_real(self):
        if self in Label.L_INDICATOR:
            return False

        return True
    
    def set_cell4D(self,cell4D):
        
        #label is the link between a cell4D and its composing cells (=different timepoints), but can only be updated when cell4D links to the label
        cell = Stack.curr_stack.d_label_cell.get(self)
        if cell:
            cell.cell4D=cell4D
        if Stack.prev_stack:
            cell = Stack.prev_stack.d_label_cell.get(self)
            if cell:
                cell.cell4D=cell4D
  
        self.cell4D=cell4D
        
        return
        
    
    def __get_color_for_label(self):
        if self.ctr==Label.EXTERIOR_CTR:
            return Label.COLOR_EXTERIOR
        if self.ctr==Label.VALIDATION_CTR:
            return Label.COLOR_VALIDATION
        if self.ctr==Label.UNASSIGNED_CTR:
            return Label.COLOR_UNASSIGNED

        
        nb_cell_colors = len(Label.L_MAIN_COLORS)
        color_pick = Label.L_MAIN_COLORS[int(self.ctr)%nb_cell_colors]
        color_dim_level = np.floor(int(self.ctr)/nb_cell_colors)
        while color_dim_level >=0:
            if color_dim_level==0:
                return color_pick
            elif color_dim_level==1:
                return color_pick + '++'
            elif color_dim_level==2:
                return color_pick + '--'
            elif color_dim_level==3:
                return color_pick + '+' 
            elif color_dim_level==4:
                return color_pick + '-'
            elif color_dim_level==5:
                return color_pick + '+++' 
            elif color_dim_level==6:
                return color_pick + '---' 
            else:
                color_dim_level -= 7
                



class Featuretracker:
    
    d_feature_ix = {'seed_centre_z':0,
                       'seed_centre_y':1,
                       'seed_centre_x':2,
                       'nb_pixels':3,
                       'name_lineage':4,
                       'volume':5
                       }   
    
    def __init__(self,nb_timesteps):
        self.nb_timesteps = nb_timesteps
        self.tracker = np.zeros((0,self.nb_timesteps,len(Featuretracker.d_feature_ix)),dtype='str') #first dimension = labels (labels function as indexes), timesteps are nb not indexes
        self.d_label_events = {}
        
        
    def add_event(self,nb_event,agents):
        
        def add_log_for_label(label,text):
            s_text = '{0} {1} {2}'.format(s_time,label.color,text)
            l_chron = self.d_label_events.get(label,[])
            l_chron.append(s_text)
            self.d_label_events[label] = l_chron
            return
        
        s_time = 'In timestep {0} '.format(Stack.curr_stack.timestep)
        
        if nb_event == 1:
            add_log_for_label(agents[0].label,
                              'was created (no checks were done)'
                              )
            
        elif nb_event == 17: #ft.add_event(17,[cluster,seeding_label])
            add_log_for_label(agents[0].label,
                              '_{0} inherits the label({1}) from the previous stack through seeding'.
                              format(agents[0].ctr,agents[1])
                              )

        elif nb_event == 19: # self.feat_tracker.add_event(19,[seeding_label,max_dt])
            add_log_for_label(agents[0],
                              '_{0} could not be seeded to the next timestep, because the max radius sphere that could be found( {1} )was smaller than the MIN_cluster_VOLUME.'.
                              format(agents[0].color,agents[1])
                              )
        else:
            print('event not handled')       
        return
        
    def expand_label_dimension(self,nb_label):
        if nb_label > self.tracker.shape[0] - 1:
            self.tracker = np.append(self.tracker,np.zeros((nb_label-(self.tracker.shape[0] - 1),self.tracker.shape[1],self.tracker.shape[2])),axis=0)
            
    def update_features_v2(self,stack,verbose=True):
        if verbose:print('updating feature tracker v2')
        
        labels_processed=[]
        for cell_i in stack.cells:
            if not cell_i.label.is_real():
                continue
            if cell_i.label.ctr > (self.tracker.shape[0] - 1):
                self.expand_label_dimension(cell_i.label.ctr)
            self.tracker[cell_i.label.ctr,stack.timestep - 1,:] = cell_i.cluster_seed.centre + [cell_i.get_volume()] + [cell_i.label.name_lineage] + [cell_i.get_volume_um()] #updates done
            labels_processed.append(cell_i.label)        
        
        return
            
    def update_features(self,stack,verbose=True):
        if verbose:print('updating feature tracker')
        
        labels_processed=[]
        for cluster in stack.get_clusters():
            if not cluster.cell.cluster_seed==cluster:
                continue # only the seed cluster will be represented in the featuretracker ! 
            if not cluster.label.is_real():
                continue
            if cluster.label.ctr > (self.tracker.shape[0] - 1):
                self.expand_label_dimension(cluster.label.ctr)
            self.tracker[cluster.label.ctr,stack.timestep - 1,:] = cluster.centre + [len(cluster.l_pixels[0])] + [len(cluster.spheres)]  #updates done
            labels_processed.append(cluster.label)
         
    def write_excel_summary_4D_v2(self,featuretracker,verbose=True):
        
        if verbose:print("###write_excel_summary_4D")
        featuretracker.d_save_info['f_name']='summary4D'
        
        ls_columns = ["cluster_label",   #this array determines the order of the columns (besides the names)
                   "color",
                   "label_parent",
                   "birth_stack_nb",
                   "lifetime",
                   "validation_chronology",
                   "validation_names",
                   "validation_error"]

        d_df = {}  # the dict that will be used to create the dataframe   
        
        
            
        #initialize the ls_columns
        l_label = []
        l_colors = []
        l_label_parent = []
        l_birth_stack_nb = []
        l_lifetime= []
        l_validation_chronology = []
        l_validation_names = []
        l_validation_error = []
        
        #compose the ls_columns
        for cell4d_i in Cell4D.l_cell4D:
            if not cell4d_i.label.is_real():
                continue
            l_label.append(cell4d_i.label.ctr)
            l_colors.append(cell4d_i.label.color)
            l_label_parent.append(cell4d_i.label.label_parent_ctr)
            l_birth_stack_nb.append(cell4d_i.birth_timestep)
            l_lifetime.append(cell4d_i.lifetime)
            l_validation_chronology.append(self.d_label_events.get(cell4d_i.label,'no chronology recorded'))
            l_validation_names.append(cell4d_i.l_validation_names)
            l_validation_error.append(cell4d_i.validation_error)

        
        d_df[ls_columns[0]] = l_label
        d_df[ls_columns[1]] = l_colors
        d_df[ls_columns[2]] = l_label_parent
        d_df[ls_columns[3]] = l_birth_stack_nb 
        d_df[ls_columns[4]] = l_lifetime  
    #         d_df[ls_columns[5]] = l_max_cluster_size 
        d_df[ls_columns[5]] = l_validation_chronology
        d_df[ls_columns[6]] = l_validation_names
        d_df[ls_columns[7]] = l_validation_error
            
        #parse the featuretracker data into the df columns
        ls_nb_timesteps = [str(i) for i in range(1, 1+self.nb_timesteps)]
        for s_feature,ix_feat in Featuretracker.d_feature_ix.items():
            ls_append_columns = ["{1}_{0}".format(s_feature,i) for i in ls_nb_timesteps]   
            ls_columns += ls_append_columns  #a list of column names for the df
            for ix_timestep, s_column in enumerate(ls_append_columns):
                d_df[s_column] = self.tracker[l_label,ix_timestep,ix_feat]  #extract the data from the featuretracker
                
  
        #save the data
        featuretracker.d_save_info['nb_stack'] = ''
        featuretracker.save_data(data=pd.DataFrame(d_df),csv_columns=ls_columns,file_ext='csv',verbose=False)
          
        return   
          
           
class Param():
    NAME_TEMP_FOLDER = 'temp'
    
    a_empty=  np.array([ [0,0,0,0,0] \
                        ,[0,0,0,0,0] \
                        ,[0,0,0,0,0] \
                        ,[0,0,0,0,0] \
                        ,[0,0,0,0,0] \
                        ])
    
    a_0 = np.array([ [0,1,1,1,0] \
                    ,[1,0,0,0,1] \
                    ,[1,0,0,0,1] \
                    ,[1,0,0,0,1] \
                    ,[0,1,1,1,0] \
                    ])
    
    a_1 = np.array([ [0,0,1,0,0] \
                    ,[0,1,1,0,0] \
                    ,[0,0,1,0,0] \
                    ,[0,0,1,0,0] \
                    ,[0,0,1,0,0] \
                    ])
    
    a_2 = np.array([ [0,1,1,1,1] \
                    ,[0,0,0,0,1] \
                    ,[1,1,1,1,1] \
                    ,[1,0,0,0,0] \
                    ,[1,1,1,1,1] \
                    ])
    
    a_3 = np.array([ [0,1,1,1,0] \
                    ,[0,0,0,0,1] \
                    ,[0,0,1,1,0] \
                    ,[0,0,0,0,1] \
                    ,[0,1,1,1,0] \
                    ])
    
    a_4 = np.array([ [0,1,0,0,1] \
                    ,[0,1,0,0,1] \
                    ,[0,1,1,1,1] \
                    ,[0,0,0,0,1] \
                    ,[0,0,0,0,1] \
                    ])
    
    a_5 = np.array([ [1,1,1,1,1] \
                    ,[1,0,0,0,0] \
                    ,[0,1,1,1,1] \
                    ,[0,0,0,1,1] \
                    ,[1,1,1,1,1] \
                    ])
    
    a_6 = np.array([ [0,0,1,0,0] \
                    ,[0,1,0,0,0] \
                    ,[0,1,1,1,0] \
                    ,[0,1,0,1,0] \
                    ,[0,1,1,1,0] \
                    ])
    
    a_7 = np.array([ [0,1,1,1,1] \
                    ,[0,0,0,1,0] \
                    ,[0,0,1,0,0] \
                    ,[0,1,0,0,0] \
                    ,[1,0,0,0,0] \
                    ])
    
    a_8 = np.array([ [0,1,1,1,0] \
                    ,[1,0,0,0,1] \
                    ,[0,1,1,1,0] \
                    ,[1,0,0,0,1] \
                    ,[0,1,1,1,0] \
                    ])
    
    
    a_9 = np.array([ [0,1,1,1,0] \
                    ,[0,1,0,1,0] \
                    ,[0,1,1,1,0] \
                    ,[0,0,0,1,0] \
                    ,[0,1,1,1,0] \
                    ])
    
    a__ = np.array([ [0,0,0,0,0] \
                    ,[0,0,0,0,0] \
                    ,[0,0,0,0,0] \
                    ,[0,0,0,0,0] \
                    ,[0,1,1,1,0] \
                    ])
    
    a_o = np.array([ [0,0,0,0,0] \
                    ,[0,1,1,1,0] \
                    ,[0,1,0,1,0] \
                    ,[0,1,1,1,0] \
                    ,[0,0,0,0,0] \
                    ])
    
    a_O = np.array([ [1,1,1,1,1] \
                    ,[1,0,0,0,1] \
                    ,[1,0,0,0,1] \
                    ,[1,0,0,0,1] \
                    ,[1,1,1,1,1] \
                    ])
    
    d_char_typeletter = {'0':a_0,
                         '1':a_1,
                         '2':a_2,
                         '3':a_3,
                         '4':a_4,
                         '5':a_5,
                         '6':a_6,
                         '7':a_7,
                         '8':a_8,
                         '9':a_9,
                         '_':a__,
                         ' ':a_empty,
                         'o':a_o,
                         'O':a_O
                         }
        
    
    def __init__(self,param_xml,filehandler):
         
        verbose=False
        param_xml.l_main_keys = ['body','spheresDT','parms']
        Param.FRAGMENTATION_LEVEL = param_xml.get_value('FRAGMENTATION_LEVEL')
        Param.z_span_max = param_xml.get_value('z_span_max')
        Param.radius_zspan = param_xml.get_value('radius_zspan')
        Param.FILTER_SPHERE_OVERLAP_MASK_PC = param_xml.get_value('FILTER_SPHERE_OVERLAP_MASK_PC')
        Param.CLUSTER_THRESHOLD_OVERLAP_MASK_PC = param_xml.get_value('CLUSTER_THRESHOLD_OVERLAP_MASK_PC')
        Param.THRESHOLD_CELL_OVERLAP_MASK_PC = param_xml.get_value('THRESHOLD_CELL_OVERLAP_MASK_PC')
        Param.MASK_AS_MEMBRANE = param_xml.get_value('MASK_AS_MEMBRANE')
        
        
        Param.MIN_SPHERE_RADIUS = param_xml.get_value('MIN_SPHERE_RADIUS')
        Param.MIN_SEED_RADIUS = param_xml.get_value('MIN_SEED_RADIUS')
        Param.MIN_CELL_RADIUS = param_xml.get_value('MIN_CELL_RADIUS')
        Param.SEEDING_Z_RANGE = param_xml.get_value('SEEDING_Z_RANGE')
        Param.THRESHOLD_CELL_DIVISION_MIN_RADIUS_RATIO = param_xml.get_value('THRESHOLD_CELL_DIVISION_MIN_RADIUS_RATIO')
        Param.MAX_NB_SPHERES = param_xml.get_value('MAX_NB_SPHERES')
        Param.D_ALGO_SWITCHES = param_xml.get_value("ALGO_SWITCHES",dtype='dict')
        

        param_xml.l_main_keys = ['body','spheresDT','parms','ALGO_SWITCHES']
        Param.FILTER_SPHERE_TOUCH_EXTERIOR_FRAME = param_xml.get_value('FILTER_SPHERE_TOUCH_EXTERIOR_FRAME')
        Param.FILTER_SPHERE_TOO_MUCH_OVERLAP_WITH_MASK = param_xml.get_value('FILTER_SPHERE_TOO_MUCH_OVERLAP_WITH_MASK')
        Param.FILTER_SPHERE_SMALL_RADIUS = param_xml.get_value('FILTER_SPHERE_SMALL_RADIUS')
        Param.FILTER_SPHERE_SEED_TOP_OR_BOTTOM = param_xml.get_value('FILTER_SPHERE_SEED_TOP_OR_BOTTOM')
        Param.FILTER_CLUSTER_SEED_TOP_OR_BOTTOM = param_xml.get_value('FILTER_CLUSTER_SEED_TOP_OR_BOTTOM')
        Param.FILTER_CLUSTER_TOUCH_EXTERIOR_FRAME = param_xml.get_value('FILTER_CLUSTER_TOUCH_EXTERIOR_FRAME')
        Param.THRESHOLD_CELL_OVERLAP_TOO_MUCH_WITH_MASK = param_xml.get_value('THRESHOLD_CELL_OVERLAP_TOO_MUCH_WITH_MASK')

        Param.VOXEL_DEPTH,_,Param.PIXEL_WIDTH = param_xml.get_value('scaling_ZYX',l_path_keys=['body','RAM'],use_main_keys=False)
        Param.XY_VS_Z = Param.VOXEL_DEPTH/Param.PIXEL_WIDTH 
        if verbose:print('temp sdt uses voxel depth of {0} and a pixel depth of {1}'.format(Param.PIXEL_WIDTH,Param.VOXEL_DEPTH))

        
          
        #to be filled dynamically
        Param.X_DIM = 0
        Param.Y_DIM = 0 
        Param.Z_DIM = 0 
        Param.T_DIM = 0 
        


    @staticmethod
    def set_img_dimensions(img_4D):
        if len(img_4D.shape)==3:
            Param.Z_DIM,Param.Y_DIM,Param.X_DIM=img_4D.shape
            Param.T_DIM = 1
        else:
            Param.T_DIM,Param.Z_DIM,Param.Y_DIM,Param.X_DIM=img_4D.shape
        return
    