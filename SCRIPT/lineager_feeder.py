'''
Created on 4 Apr 2019
interface for the Lineager application (WormDeLux)

lineage tree is built up using the 
* cluster labels (own cluster_label and parent_label)
* cluster.division indication to mark cell division (marked 'NEW CELL')
* cell_rank : to mark the seed clusters (only seed clusters get an entry)

al this information comes from the cluster_summary report from spheresDT

example : a cell division in the cluster summary looks likes this
t   color cluster_label         parent_label 
1	blue	    4	    seed	    99	
2	blue	    4	    seed	    99	        DIVIDED
2	yellow+++	46	    seed	    4	        NEW CELL

@author: wimth
'''


import pandas as pd
from collections import namedtuple
import os
from datamodel import Label
from helper_functions import slice_tif


def generate_lineager_input_files(param_xml,filehandler):
        
    def use_cluster_as_representative_of_cell():  ##deprecated
        if row_cluster.cluster_label in [Label.EXTERIOR,Label.VALIDATION,Label.UNASSIGNED]:
            return False
        
        
        if row_cluster.cluster_label not in d_label_lineagedata: return True #first cluster will always get an entry, but maybe overwritten later on
        
        prev_ic_labelled_by_seeding = d_label_seeded.get(row_cluster.cluster_label)
        if prev_ic_labelled_by_seeding: return False  #the previous cluster came first and was seeded, this always wins
        
        if 'seeded' in row_cluster.label_hist:return True#the new cluster is seeded and the old one not, so the new cluster wins
     
            
        prev_z,z_competitor_cluster = d_label_deltaZ[row_cluster.cluster_label]  #both unseeded, let the delta z decide which is best
        if abs(prev_z - row_cluster.z) < abs(prev_z - z_competitor_cluster): return True
            
        return False
        
    def read_parms(param_xml):
        param_xml.l_main_keys = ['body','lineager_feeder']
        INPUT_FILE = param_xml.get_value('INPUT_FILE',['paths'])
        MAKE_IMAGE_SEQUENCE = param_xml.get_value('MAKE_IMAGE_SEQUENCE',['flowcontrol'])
        offset_lineage_files = param_xml.get_value('offset_lineage_files',['flowcontrol'])
        if offset_lineage_files==-1:
            offset_lineage_files = param_xml.get_value('l_stack_number',['body','preprocessing','flowcontrol'],use_main_keys=False)[0] - 1
        param_xml.l_main_keys = ['body','MAIN']
        img_raw_file = param_xml.get_value('img_raw_file',['paths'])

        return INPUT_FILE, MAKE_IMAGE_SEQUENCE,img_raw_file,offset_lineage_files
    
    def read_input_file():
        full_path = INPUT_FILE if INPUT_FILE else filehandler.get_f_name('cluster_summary')
        df_ipt = pd.read_csv(full_path,delimiter=',')
    
        return df_ipt
    
    def write_lineage_file_prev():
        filehandler.d_save_info['f_name'] = 't{0}-nuclei'.format(str(nb_time_prev+offset_lineage_files).rjust(3,'0'))
        filehandler.save_data(d_label_lineagedata_prev,file_ext='lineage',verbose=False)
        print('|',end='')
        return
    
    def make_image_sequence():
        print('  - making image sequence')
        
        # filehandler.d_load_info['load_dir'] = IMG_RAW_DIR
        filehandler.d_load_info['f_name'] = str(img_raw_file)
        a_raw = filehandler.load_tif()
        if len(a_raw.shape)==3:
            z,_,_ = a_raw.shape
            t=1
        else:
            t,z,_,_ = a_raw.shape
        
        filehandler.extend_save_info(extra_dir_1='Image_Sequence',reset_to_snapshot=True)
        
        for t_i in range(t):
            for z_i in range(z):
                a_slice = slice_tif(a_raw,ix_slice_time=[t_i,t_i+1],ix_slice_z=[z_i,z_i+1],verbose=False)
                filehandler.d_save_info['f_name'] = img_raw_file.stem + '_t' +  str(t_i+1).zfill(3) + '_z' + str(z_i+1).zfill(3)
                filehandler.save_data(a_slice,verbose=False)
                print('|',end='')

        filehandler.reset_save_info()
        
        return
        
    #MAIN#################################################################   
    
    INPUT_FILE, MAKE_IMAGE_SEQUENCE,img_raw_file,offset_lineage_files = read_parms(param_xml)
    filehandler.extend_save_info(extra_dir_1='005_Lineager_feeder',from_root=True,take_snapshot_after=True)
    filehandler.extend_save_info(extra_dir_1='lineage_files')
    
    df_ipt = read_input_file()
    
    LineageData=namedtuple('LineageData',['id','unk','link_past','link_future_1','link_future_2','x','y','z','diameter','name'])
    d_label_lineagedata = {};d_label_lineagedata_prev = {};d_label_deltaZ={};d_label_seeded={}
    nb_time_prev = 0

    
    time_range = range(df_ipt['stack_nb'].min(),df_ipt['stack_nb'].max()+1)
    for nb_time in time_range:
        for _,row_cluster in df_ipt[df_ipt.stack_nb==nb_time].iterrows():
    #         if use_cluster_as_representative_of_cell():  
            if row_cluster.cluster_label in [Label.EXTERIOR_CTR,Label.VALIDATION_CTR,Label.UNASSIGNED_CTR,Label.FILTERED_CTR]:
                continue
            if row_cluster.cell_rank=='seed':  
                nt_lineage = LineageData(id=row_cluster.cluster_label, #i use the same id's in every stack ! so not unique over stacks
                                         unk=1,
                                         link_past=row_cluster.cluster_label,
                                         link_future_1=row_cluster.cluster_label,
                                         link_future_2=-1,
                                         x=row_cluster.x,
                                         y=row_cluster.y,
                                         z=row_cluster.z,
                                         diameter=8.5,
                                         name='{0}({1}_{2})'.format(row_cluster.name_lineage,row_cluster.cluster_label,row_cluster.color)
                                         )
                
                if nb_time==1:
                    nt_lineage=nt_lineage._replace(link_past=-1)
                    
                            
    #             if row_cluster.division=='DIVIDED':  # a new cell must correct 2 records, current and previous
    #                 nt_lineage=nt_lineage._replace(link_past=row_cluster.name_lineage[0:-2]) #correct the past link of the current record (=remove 'd1')
    #                 if d_label_lineagedata_prev:
    #                     nt_lineage_prev = d_label_lineagedata_prev.get(row_cluster.label_parent)
    #                     if nt_lineage_prev:
    #                         nt_lineage_prev=nt_lineage_prev._replace(link_future_1=row_cluster.name_lineage) #add the second future link to the previous record
    #                         d_label_lineagedata_prev[row_cluster.label_parent]=nt_lineage_prev
    #                                                 
                                                                 
                if row_cluster.division=='NEW CELL':  # a new cell must correct 2 records, current and previous
                    nt_lineage=nt_lineage._replace(link_past=row_cluster.label_parent) #correct the past link of the current record (=remove 'd2')
                    if d_label_lineagedata_prev:
                        nt_lineage_prev = d_label_lineagedata_prev.get(row_cluster.label_parent)
                        if nt_lineage_prev:
    #                         nt_lineage_prev=nt_lineage_prev._replace(link_future_1=row_cluster.name_lineage[0:-2] + 'd1') #add the first future link to the previous record
                            nt_lineage_prev=nt_lineage_prev._replace(link_future_2=row_cluster.cluster_label) #add the second future link to the previous record
                            d_label_lineagedata_prev[row_cluster.label_parent]=nt_lineage_prev
    
    
                
                                       
                d_label_lineagedata[row_cluster.cluster_label] = nt_lineage # a new entry or overruling of the previous cluster entry
                
    #             if row_cluster.cluster_label in d_label_lineagedata_prev:   #updating the smallest delta z
    #                 nt_lineage_prev = d_label_lineagedata_prev.get(row_cluster.cluster_label)
    #                 d_label_deltaZ[row_cluster.cluster_label] = [nt_lineage_prev.z,nt_lineage.z]
    #             else:
    #                  d_label_deltaZ[row_cluster.cluster_label] = [nt_lineage.z,nt_lineage.z]
    #             
    #             
    #             d_label_seeded[row_cluster.cluster_label] = 'seeded' in row_cluster.label_hist
            else:
                continue
            
        if d_label_lineagedata_prev:
            write_lineage_file_prev()
        
        nb_time_prev = nb_time
        d_label_lineagedata_prev = d_label_lineagedata   
        d_label_lineagedata = {}
        d_label_deltaZ={}
    
    write_lineage_file_prev()

    if MAKE_IMAGE_SEQUENCE:make_image_sequence()
    
    print('END Lineager Feeder')



