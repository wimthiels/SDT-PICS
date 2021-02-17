#!/usr/bin/python3.7

import pandas as pd
from collections import namedtuple
import os
from helper_functions import get_3D_sphere_coordinates
import numpy as np


def generate_CTC_TRA_files(param_xml,filehandler):

    def read_parms(param_xml):
        param_xml.l_main_keys = ['body','CTC_TRA']
        input_file_csv = param_xml.get_value('input_file_csv',['paths'])
        if not input_file_csv:
            input_file_csv = filehandler.get_f_name('cluster_summary')
        input_file_clusterview = param_xml.get_value('input_file_clusterview',['paths'])
        if not input_file_clusterview:
            input_file_clusterview = filehandler.get_f_name('clusterview')
        INPUT_FILE_VALIDATION  = param_xml.get_value('INPUT_FILE_VALIDATION',['paths'])
        CTC_TRA_EXE  = param_xml.get_value('CTC_TRA_EXE',['paths'])
        
        GT_mode = param_xml.get_value('GT_mode',['flowcontrol'])
        
        return input_file_csv, input_file_clusterview, INPUT_FILE_VALIDATION,CTC_TRA_EXE , GT_mode
        
    def read_input_files():
        full_path = input_file_csv if input_file_csv else filehandler.get_f_name('cluster_summary')
        df_ipt = pd.read_csv(full_path,delimiter=',')
        
        filehandler.d_load_info['load_dir'] = ''
        filehandler.d_load_info['f_name'] = input_file_clusterview if input_file_clusterview else filehandler.get_f_name('clusterview')
        
        a_clusterview = filehandler.load_tif()
        
        return df_ipt,a_clusterview

    def read_validation_nuclei():
        print("#load the Cell nuclei for validation")
        df_val = pd.read_csv(INPUT_FILE_VALIDATION ,delimiter='\t')
        df_val["x"] = pd.to_numeric(df_val["x"],downcast='integer')
        df_val["y"] = pd.to_numeric(df_val["y"],downcast='integer')
        df_val["z"] = pd.to_numeric(df_val["z"],downcast='integer')
            
        return df_val


    #SDT_MAIN################################################################# 
    input_file_csv, input_file_clusterview, INPUT_FILE_VALIDATION ,CTC_TRA_EXE, GT_mode = read_parms(param_xml)
    filehandler.extend_save_info(extra_dir_1='006_CTC_TRA',from_root=True,take_snapshot_after=True)
    filehandler.extend_save_info(extra_dir_1='01_RES')
    df_ipt, a_clusterview = read_input_files()
    
    #Change the clusterview to a tracking_clusterview
    LABEL_BACKGROUND=0
    a_clusterview[a_clusterview==1]=LABEL_BACKGROUND #set background to zero
    a_ctc_mask = np.zeros_like(a_clusterview)
    LineageData=namedtuple('LineageData',['starttime','endtime','link_parent'])  # zero-based time indexing
    d_lineagectr_lineagedata = {} #every item here will be a record in the res.track.txt (lineagectr is the cell counter used for tracking=cell_ctr_CTC)

    if input_file_csv:
        #FIRST USE THE CLUSTERVIEW EXCEL TO GET THE CTC tracking files (RES)
        #-------------------------------------------------------------------
        d_name_lineage_ctr = {} #key=lineage_name, value = lineage_counter.  needed to retrieve link-parent
        
        time_range = range(df_ipt['stack_nb'].min()-1,df_ipt['stack_nb'].max())  #time for CTC must be zero-indexed
        for nb_time in time_range:
            print(nb_time)
            a_lineagectr = np.zeros_like(a_clusterview[nb_time],dtype='uint16')
            for _,row_cluster in df_ipt[df_ipt.stack_nb==nb_time+1].iterrows():    #row_cluster is row from the cluster_summary

                if row_cluster.cluster_label in [9998,9999]:
                    #a_clusterview[a_clusterview==int(row_cluster.cluster_counter)]= LABEL_BACKGROUND
                    continue
                if int(row_cluster.ctr_lineage) in d_lineagectr_lineagedata:
                    nt_lineage = d_lineagectr_lineagedata[int(row_cluster.ctr_lineage)]
                    nt_lineage= nt_lineage._replace(endtime=int(row_cluster.stack_nb) -1)
                    d_lineagectr_lineagedata[int(row_cluster.ctr_lineage)]= nt_lineage
                else:
                    if isinstance(row_cluster.name_lineage, (int,float)):
                        link_parent = 0
                    elif len(row_cluster.name_lineage) >2 and row_cluster.name_lineage[-2]=='d': # get parent = the cell without d1/d2 suffix
                        link_parent = d_name_lineage_ctr.get(row_cluster.name_lineage[:-2],0)
                    else:
                        link_parent = 0
                            
                    nt_lineage = LineageData(starttime=int(row_cluster.stack_nb) -1, 
                                             endtime=int(row_cluster.stack_nb) -1, 
                                             link_parent=int(link_parent))
                    
        
                    d_lineagectr_lineagedata[int(row_cluster.ctr_lineage)] = nt_lineage
                    d_name_lineage_ctr[row_cluster.name_lineage] = row_cluster.ctr_lineage
                
                a_lineagectr[np.where(a_clusterview[nb_time]==int(row_cluster.cluster_counter))]= int(row_cluster.ctr_lineage)  #make tiff
                #a_clusterview[a_clusterview==int(row_cluster.cluster_counter)]=int(row_cluster.ctr_lineage)  #make tiff
             
            filehandler.d_save_info['f_name'] = 'mask' + str(nb_time).zfill(3)
            filehandler.save_data(a_lineagectr,resolution='uint16',clip_not_scale=False,verbose=False)
        
        full_path = os.path.join(filehandler.get_save_location(),"res_track.txt")
        print(full_path, "items", len(d_lineagectr_lineagedata))
        with open(full_path, 'w',newline='') as f:
             
            for lineage_ctr_i, lineage_data_i in d_lineagectr_lineagedata.items():
                rec_track = str(lineage_ctr_i) + " " + "{0}\n".format(' '.join(map(str, lineage_data_i)))
                print(rec_track)
                f.write(rec_track)
             
        print('the end : excel')

    #SECOND,  USE THE LINEAGER FILES TO GET THE CTC tracking files (GT)
    #-------------------------------------------------------------------

    #the saved lineage file is made quite confusing. to build the lineage you must traverse the links, timestep per timestep, the self counter changes constantly
    # a cell division is detectable is child2 is filled in.  I will use the first self number as the cell counter.  

    def create_new_cell_track(LINEAGE_CTR):
       
        LINEAGE_CTR += 1    # draw a new number for the cell
        d_self_lineager_ctr_curr[row_i.self] = LINEAGE_CTR  #
        
        parent_cell = d_self_lineager_ctr_prev.get(row_i.parent,0)
        
        nt_lineage = LineageData(starttime=nb_time, 
                                             endtime=nb_time,
                                             link_parent=parent_cell)
        d_lineagectr_lineagedata[LINEAGE_CTR] = nt_lineage
        
        d_lineagectr_celldesc[LINEAGE_CTR] = row_i.name
        
        return LINEAGE_CTR

    def extend_end_time_cell_track():
        lineage_ctr_cell = d_self_lineager_ctr_prev[row_i.parent]
        nt_lineage = d_lineagectr_lineagedata[lineage_ctr_cell]
        nt_lineage= nt_lineage._replace(endtime=nb_time)
        d_lineagectr_lineagedata[lineage_ctr_cell]= nt_lineage

        return lineage_ctr_cell

    d_self_lineager_ctr_curr = {}  #this dict will link every self counter to the actual lineager counter. only needed for looking 1 timestep in the past
    d_self_lineager_ctr_prev = {} 
    s_cell_division_curr = set()
    s_cell_division_prev = set()
     
    Z_DIM,Y_DIM,X_DIM = a_clusterview.shape[1:]
    filehandler.extend_save_info(extra_dir_1='01_GT',extra_dir_2='TRA',reset_to_snapshot=True,take_snapshot_after=True)
    
    # f_name = filehandler.get_f_name('IMG_RAW_FILE')
    # pixel_width ,_,voxel_depth = param_xml.get_value('dxyz',l_path_keys=['body','spheresDT','parms'],use_main_keys=False,lookup=f_name)
    voxel_depth,_,pixel_width= param_xml.get_value('scaling_ZYX',l_path_keys=['body','RAM'],use_main_keys=False)
    XY_VS_Z = voxel_depth/pixel_width

    d_lineagectr_lineagedata = {} #this is the data for the CTC tracking text file
    LINEAGE_CTR = 0
    d_lineagectr_celldesc = {}

    if GT_mode:
        df_val = read_validation_nuclei()
        time_range = range(df_val['time'].min()-1,df_val['time'].max())
        for nb_time in time_range:
            a_GT_mask = np.zeros((Z_DIM,Y_DIM,X_DIM))
            print(nb_time)
            for _,row_i in df_val[df_val.time==nb_time+1].iterrows():
                
                #determine the cell-nb for this cell (either new number, of inherit from previous stack)
                if (row_i.parent==-1) or (row_i.parent in s_cell_division_prev) or (not d_self_lineager_ctr_prev.get(row_i.parent)):
                    LINEAGE_CTR = create_new_cell_track(LINEAGE_CTR)
                    lineage_ctr_cell=LINEAGE_CTR    
                else:
                    lineage_ctr_cell= extend_end_time_cell_track()
                
                d_self_lineager_ctr_curr[row_i.self] = lineage_ctr_cell
                
                if row_i.child2!=-1:  #remember cell division
                    s_cell_division_curr.add(row_i.self)
                    
                # update the image
                centre_ZYX = [int(round(row_i.z)),row_i.y,row_i.x]
                l_pixels = get_3D_sphere_coordinates(a_coordZYX=centre_ZYX,
                                                      radius=row_i.radius,
                                                      limit_shape=(Z_DIM,Y_DIM,X_DIM),
                                                      xy_vs_z=XY_VS_Z)  #[Z,Y,X] 
                
                a_GT_mask[l_pixels]=lineage_ctr_cell 
                
            # save image for this timestep
            filehandler.d_save_info['f_name'] = 'man_track' + str(nb_time).zfill(3)
            filehandler.save_data(a_GT_mask.astype('uint16'),resolution='uint16',clip_not_scale=False,verbose=False)
            
            # shift to next timestep
            d_self_lineager_ctr_prev = d_self_lineager_ctr_curr.copy()
            d_self_lineager_ctr_curr = {}
            s_cell_division_prev = s_cell_division_curr.copy()
            s_cell_division_curr = set()
            
        # write lineage text file in CTC format
        full_path = os.path.join(filehandler.get_save_location(),"man_track.txt")
        print(full_path, "items", len(d_lineagectr_lineagedata))
        with open(full_path, 'w',newline='') as f:
        
            for lineage_ctr_i, lineage_data_i in d_lineagectr_lineagedata.items():
                rec_track = str(lineage_ctr_i) + " " + "{0}\n".format(' '.join(map(str, lineage_data_i)))
                print(rec_track)
                f.write(rec_track)
        
        
        #write file as extra information (link cell number with actual name)
        full_path = os.path.join(filehandler.get_save_location(),"link_file.txt")
        with open(full_path, 'w') as f:
            for i,j in d_lineagectr_celldesc.items():
                s_rec = str(i) + " = " + str(j)
                f.write(s_rec)
               
    print('the end , GT')
    
    #THIRD : EXECUTE CTC TRA score
    filehandler.extend_save_info(extra_dir_1='006_CTC_TRA',from_root=True)
    
    cmd_TRA =  str(CTC_TRA_EXE) + " " + filehandler.get_save_location() + ' 01'
    print('DEBUG cmd_TRA =', cmd_TRA)
    os.system(cmd_TRA)
    print('DEBUG after executing TRA')
    

            