'''
Created on 27 Mar 2019

@author: wimth
'''

import os
import numpy as np
from tifffile import TiffFile,imsave
from skimage import img_as_uint,img_as_ubyte
import pandas as pd
from matplotlib.pyplot import savefig
from pathlib import Path

class FileHandler():
    '''
    helper class : loading and saving data
    '''

    
    d_locations = {'1':r'E:\Downloads\@tiffhandling',
               'download':r'E:\Downloads\@tiffhandling',
               '2':r'C:\Users\wimth\OneDrive\Eclipse_MyFiles\Thesis',
               'onedrive':r'C:\Users\wimth\OneDrive\Eclipse_MyFiles\Thesis',
               '3':r'C:\Users\wimth\Documents\THESIS\baseimages',
               'base':r'C:\Users\wimth\Documents\THESIS\baseimages',
               '4':r'C:\Users\wimth\Documents\THESIS\worm-img-TL1',
               'rob':r'C:\Users\wimth\Documents\THESIS\worm-img-TL1',
               'mpacts_local':'/mnt/c/Users/wimth/Documents/THESIS/MPACTS/C_elegans_Membrane_signal_for_toy_model',
               'mpacts_vsc':'/data/leuven/317/vsc31741/mpacts_input'}

    d_f_names = {'1':'raw4D',
                 '2':'raw4D_cutcorner',  #pick this for the first timelapsee
                 '3':'raw3D_timesquash',
                 '4':'raw3D_timesquash_cutcorner',
                 '5':'4D_membrane_contour',  #francesca tif
                 '6':'3D_membrane_contour_1',
                 '7':'3Dmask_simple',   # a simple ellipsoid, not deformed
                 '8':'3Dmask_simple_perimeter',   
                 '9':'raw4D_stack1_BinaryThresholdImageFilter_Out1', #straight from race , the first stack
                 '10':'exterior_overlay_3DtimesquashChanVese',  #see morphacwe2D
                 '11':'raw4D_RGB',  #the raw4D in RGB form
                 '12':'TL2_crop_1ch_emb1of3' #at baseimages, 1 cell from the second timelapse 
                 }
    
    d_load_info_cleared = {'load_dir':'',
               'sub_dir_1':'',
               'sub_dir_2':'',
               'f_name':''}
    
    imagej_metadata = """ImageJ=1.52e
    images={nr_images}
    frames={nr_channels}
    slices={nr_slices}
    frames={nr_frames}
    hyperstack=true
    mode=gray
    loop=false"""

    imagej_metadata2 = """ImageJ=1.52e
    images={nr_images}
    slices={nr_slices}
    frames={nr_frames}
    hyperstack=true
    mode=gray
    loop=false"""




    def __init__(self,filehandler_template=None):
        '''
        Constructor
        '''
        
        if filehandler_template:
            self.d_load_info=filehandler_template.d_load_info
            self.d_save_info=filehandler_template.d_save_info
            self.d_save_info_root=filehandler_template.d_save_info_root
            self.d_save_info_snapshot = filehandler_template.d_save_info_snapshot
            self.d_load_info_snapshot = filehandler_template.d_load_info_snapshot
            
        else:
            self.d_load_info = {'load_dir':'',
                   'sub_dir_1':'',
                   'sub_dir_2':'',
                   'f_name':''}
    
    
            self.d_save_info = {'save_dir':'',
                           'sub_dir_1':'',
                           'sub_dir_2':'',
                           'sub_dir_3':'',
                           'sub_dir_4':'',
                           'sub_dir_5':'',
                           'sub_dir_6':'',
                           'pre_f_name':'',   #prefix of the file name
                           'nb_stack':'',
                           'f_name':'default'}
        
            self.d_save_info_root = self.d_save_info.copy()
            self.d_save_info_snapshot = self.d_save_info.copy()
            self.d_load_info_snapshot = self.d_load_info.copy()
        
        self.d_name_f_name = {} # a dict for storing variable file names
        self.d_loaded_files = {}  #a dict for files (to prevent reloading the same image)
    
    
    def set_save_info(self,key,value):
        self.d_save_info_snapshot = self.d_save_info.copy()
        self.d_save_info[key]=value
        
        return
    
    def take_snapshot(self):
        self.d_save_info_snapshot = self.d_save_info.copy()
        self.d_load_info_snapshot = self.d_load_info.copy()
        return
    
    def store_f_name(self,key,value):
        self.d_name_f_name[key] = value
        return
    
    def get_f_name(self,key):
        return self.d_name_f_name[key]

    def init_f_name(self,key):
        file = Path(self.get_f_name(key))
        file.unlink() if file.exists() else print('') 
        return 
        
    
    def set_save_info_root(self):
        self.d_save_info_root = self.d_save_info.copy()
        return
    
    def set_save_info_to_root(self):
        self.d_save_info = self.d_save_info_root.copy()
        return

    def pop_save_info(self,set_root=False):
        if self.d_save_info['sub_dir_6']:
            self.d_save_info['sub_dir_6']=''
        
        elif self.d_save_info['sub_dir_5']:
            self.d_save_info['sub_dir_5']=''
           
        
        elif self.d_save_info['sub_dir_4']:
            self.d_save_info['sub_dir_4']=''
            
        
        elif self.d_save_info['sub_dir_3']:
            self.d_save_info['sub_dir_3']=''
            
        
        elif self.d_save_info['sub_dir_2']:
            self.d_save_info['sub_dir_2']=''
            
        
        elif self.d_save_info['sub_dir_1']:
            self.d_save_info['sub_dir_1']=''
            
        if set_root:
            self.d_save_info_root = self.d_save_info.copy()
        
        return
        


    def extend_save_info(self,extra_dir_1='',extra_dir_2='',f_name='',from_root=False,take_snapshot_before=False,take_snapshot_after=False,set_root=False,reset_to_snapshot=False):
        
        def add_extra_dir(extra_dir_x):
            if extra_dir_x:
                if self.d_save_info['sub_dir_1']:
                    if self.d_save_info['sub_dir_2']:
                        if self.d_save_info['sub_dir_3']:
                            if self.d_save_info['sub_dir_4']:
                                if self.d_save_info['sub_dir_5']:
                                    self.d_save_info['sub_dir_6'] = extra_dir_x
                                else:
                                    self.d_save_info['sub_dir_5'] = extra_dir_x
                            else:
                                self.d_save_info['sub_dir_4'] = extra_dir_x
                        else:
                            self.d_save_info['sub_dir_3'] = extra_dir_x
                    else:
                        self.d_save_info['sub_dir_2'] = extra_dir_x
                else:
                    self.d_save_info['sub_dir_1'] = extra_dir_x
            else:
                return
        
            return
        
        if take_snapshot_before:self.d_save_info_snapshot=self.d_save_info.copy()
        
        if from_root:
            self.d_save_info = self.d_save_info_root.copy()
        if reset_to_snapshot:
            self.reset_save_info()
        
        add_extra_dir(extra_dir_1)
        add_extra_dir(extra_dir_2)
        
        if f_name:
            self.d_save_info['f_name']=f_name
        
        if set_root:
            self.set_save_info_root()
            
        if take_snapshot_after:self.d_save_info_snapshot=self.d_save_info.copy()
        
        return
        
    def clear_save_dir(self):
        save_dir = self.get_save_location()
        if save_dir:
            for file in Path(save_dir).iterdir(): 
                os.remove(file)

        return
        
    def reset_save_info(self):
        self.d_save_info = self.d_save_info_snapshot.copy()
        return
    
    def reset_to_snapshot(self):
        self.d_save_info = self.d_save_info_snapshot.copy()
        self.d_load_info = self.d_load_info_snapshot.copy()
        return
        
        
    def load_tif(self,add_extension=False,storage_name = '',use_mem_map=False,input_shape=None,verbose=True):
        '''
        load a tif into an array
        dimensions of length 1 will be removed
        :param verbose:
        '''
        
        if storage_name and storage_name in self.d_loaded_files.keys(): return self.d_loaded_files[storage_name]  #no reloading of files
        
        load_dir = self.d_load_info['load_dir']
        sub_dir_1 = self.d_load_info['sub_dir_1']
        sub_dir_2 = self.d_load_info['sub_dir_2']
        f_name = self.d_load_info['f_name']
        
        load_dir=FileHandler.d_locations.get(load_dir,load_dir)
        f_name=FileHandler.d_locations.get(f_name,f_name)
        
        f_name =  '{0}.tif'.format(f_name) if add_extension else '{0}'.format(f_name)
            
        full_path=os.path.join(load_dir,sub_dir_1,sub_dir_2, f_name)
        
        if use_mem_map:
            a_tif = np.memmap(full_path, dtype='uint16', mode='r', input_shape=input_shape)
        else:  
            with TiffFile(full_path) as tif:
                a_tif = tif.asarray()
        
        if storage_name:
            self.d_loaded_files[storage_name] = a_tif
        
        if verbose:print("file loaded->{0}, with input_shape={1}{4} from location{2}{3}{5}".format(f_name,a_tif.shape,'\n',full_path,a_tif.dtype,'\n')) 
        return a_tif


    def save_data(self,data,pre_f_name='', f_name='default',
                  verbose=True,file_ext='tif', RGB=False, 
                  resolution='uint16',clip_not_scale=False,
                  use_nb_stack=False,add_ext=True,
                  csv_columns=None, fiji=False):
        '''
        save data as a tif file, excel file or image in a certain location 
        the data can be a numpy array (tif or excel) or a figure (png)
        (usage : save_data(data,**d_save_info))
        (extra:fill in a file_ext not in the list to just get back the save directory)
        
        :param data: tiff to be saved
        :param f_extension: tif, (implicit assumptions -> dimensions are ordered according to CTZYXR (C=channel , R = RGB value)
                            xlsx, png, txt, csv, lineage,stl
        :param RGB: a tif will be written with RGB.  because these are always tifs for visualization, these will be automatically converted to 
                    fiji format (=tif but with dimensions flipped (TZCYXR)
        resolution : default the data will be converted to a uint16
        '''
        save_dir = self.d_save_info['save_dir']
        sub_dir_1 = self.d_save_info['sub_dir_1']
        sub_dir_2 = self.d_save_info['sub_dir_2']
        sub_dir_3 = self.d_save_info['sub_dir_3']
        sub_dir_4 = self.d_save_info['sub_dir_4']
        sub_dir_5 = self.d_save_info['sub_dir_5']
        sub_dir_6 = self.d_save_info['sub_dir_6']
        pre_f_name=self.d_save_info['pre_f_name']
        if use_nb_stack:
            nb_stack=self.d_save_info['nb_stack']
        else:
            nb_stack=''
        f_name = self.d_save_info['f_name']
       
        
        save_dir=FileHandler.d_locations.get(save_dir,save_dir)
        f_name=FileHandler.d_locations.get(f_name,f_name)                      
    
        if sub_dir_1:
            directory = os.path.join(save_dir,sub_dir_1)
            if not os.path.exists(directory):
                if verbose:print('creating a subdir1', directory)
                os.makedirs(directory)
                
        if sub_dir_2:
            directory = os.path.join(save_dir,sub_dir_1,sub_dir_2)
            if not os.path.exists(directory):
                if verbose:print('creating a subdir2', directory)
                os.makedirs(directory)
                
        if sub_dir_3:
            directory = os.path.join(save_dir,sub_dir_1,sub_dir_2, sub_dir_3)
            if not os.path.exists(directory):
                if verbose:print('creating a subdir3', directory)
                os.makedirs(directory)   
                
        if sub_dir_4:
            directory = os.path.join(save_dir,sub_dir_1,sub_dir_2, sub_dir_3,sub_dir_4)
            if not os.path.exists(directory):
                if verbose:print('creating a subdir4', directory)
                os.makedirs(directory)      
                
        if sub_dir_5:
            directory = os.path.join(save_dir,sub_dir_1,sub_dir_2, sub_dir_3,sub_dir_4,sub_dir_5)
            if not os.path.exists(directory):
                if verbose:print('creating a subdir5', directory)
                os.makedirs(directory) 
                
        if sub_dir_6:
            directory = os.path.join(save_dir,sub_dir_1,sub_dir_2, sub_dir_3,sub_dir_4,sub_dir_5,sub_dir_6)
            if not os.path.exists(directory):
                if verbose:print('creating a subdir6', directory)
                os.makedirs(directory)    
    
        if nb_stack:     
            full_f_name = pre_f_name + 'ST' + str(nb_stack).zfill(2) +'_' + f_name
        else:
            full_f_name = pre_f_name + f_name
        
        if file_ext=='tif':
            #comment : use np.abs(np_array).astype('uint16') : to get no scaling and clipping 
            #if data.dtype not in ('uint16','uint8','int'):data=img_as_uint(data)  #WARNING if values exist above max, img_as_uint returns all zeros !!
            if resolution in ['float64','float32']:
                data=data.astype('float32') #for now change both to float32 which is compatible with fiji (float64 can be saved as bigtiff however)
            if resolution=='uint16':
                if data.dtype in ['float64','float32']:
                    data=data.astype(int,casting='unsafe')
                if str(data.dtype) not in ('uint16'):
                    data=img_as_uint(data) #Negative input values will be clipped. Positive values are scaled between 0 and 255.
            elif resolution=='uint8':
                if data.dtype  in ['float64','float32']:
                    data=data.astype(int,casting='unsafe') #Images of type float must be between -1 and 1.for img_as_ubyte
                if str(data.dtype) not in ('uint8'):
                    if clip_not_scale:
                        data=img_as_ubyte(np.where(data>255,255,data))
                    else:
                        data=img_as_ubyte(data) 
                            
            if add_ext:
                full_path=os.path.join(save_dir,sub_dir_1,sub_dir_2, sub_dir_3,sub_dir_4,sub_dir_5,sub_dir_6,full_f_name + '.tif')
            else:
                full_path=os.path.join(save_dir,sub_dir_1,sub_dir_2, sub_dir_3,sub_dir_4,sub_dir_5,sub_dir_6,full_f_name)
        
            
            if RGB:
                if len(data.shape)==6:
                    c,t,z,y,x,rgb = data.shape
                    data= np.moveaxis(data, [0, 1, 2], [2, 0, 1])  #fiji always interprets TZCYXR (while input has shape CTZYXR)
                    nr_channels = c 
                    nr_slices = z
                    nr_frames = t
                    hyperstack = 'true' 
                    s_axes = 'TZCYXS'
                if len(data.shape)==5:   #4D RGB file
                    t,z,y,x,rgb = data.shape 
                    nr_channels = 1
                    data = data.reshape(t,z,nr_channels,y,x,rgb)  # put a 1 dimensional channel in there to comply to fiji format
                    nr_slices = z
                    nr_frames = t
                    hyperstack = 'true'
                    s_axes = 'TZCYXS'
                if len(data.shape)==4:   #3D RGB file
                    z,y,x,rgb = data.shape 
                    nr_channels = 1
                    nr_frames = 1
                    data = data.reshape(z,nr_channels,y,x,rgb)  # put a 1 dimensional channel in there (
                    nr_slices = z
                    hyperstack = 'false'  
                    s_axes = 'TZCYXS'
                              
                nr_images = nr_channels *nr_slices * nr_frames
                metadata = FileHandler.imagej_metadata.format(nr_images=nr_images,nr_slices=nr_slices,nr_channels=nr_channels,nr_frames=nr_frames,hyperstack=hyperstack)
                imsave(full_path, data,imagej=metadata,metadata={'axes': s_axes})  
            else:
                if len(data.shape)==4:   #4D non-RGB file
                    t,z,y,x = data.shape 
                    if fiji:
                        data = data.reshape(t,z,nr_channels,y,x)  # put a 1 dimensional channel in there to comply to fiji format
                    metadata = FileHandler.imagej_metadata.format(nr_images=z * t,nr_slices=z,nr_frames=t,nr_channels=1,hyperstack=True)
                    imsave(full_path, data,imagej=metadata,metadata={'axes': 'TZYX'}) #
                else:
                    imsave(full_path, data)
            
            if verbose:print("*TIF file written->{0}, with shape={1}{4} to location{2}{3}{5}".format(f_name,data.shape,'\n',full_path,data.dtype,'\n')) 
            
        elif file_ext == 'xlsx':
            a_2D = data.reshape(data.shape[-2],data.shape[-1])
            df = pd.DataFrame (a_2D)
            full_path=os.path.join(save_dir,sub_dir_1,sub_dir_2, sub_dir_3,sub_dir_4,sub_dir_5,sub_dir_6,full_f_name + '.xlsx')
            df.to_excel(full_path, index=False)
            if verbose:print("excel file written->{0}.xlsx, with shape{1} to location{2}{3}{4}".format(f_name,a_2D.shape,'\n',full_path,'\n')) 
        
        elif file_ext == 'png':
            full_path=os.path.join(save_dir,sub_dir_1,sub_dir_2, sub_dir_3,sub_dir_4,sub_dir_5,sub_dir_6,full_f_name + '.png')
            savefig(full_path)
            if verbose:print("*PNG file written->{0}.png, to location{1}{2}{3}".format(f_name,'\n',full_path,'\n')) 
        
        elif file_ext == 'txt':
            full_path=os.path.join(save_dir,sub_dir_1,sub_dir_2, sub_dir_3,sub_dir_4,sub_dir_5,sub_dir_6,full_f_name + '.txt')
            with open(full_path, "a") as text_file:
                print(data, file=text_file)
            if verbose:print("*TXT file written->{0}.txt, to location{1}{2}{3}".format(f_name,'\n',full_path,'\n'))    
        
        elif file_ext == 'csv':
            full_path=os.path.join(save_dir,sub_dir_1,sub_dir_2, sub_dir_3,sub_dir_4,sub_dir_5,sub_dir_6,full_f_name + '.csv')
            if os.path.isfile(full_path) :
                with open(full_path, 'a') as f:
                    data.to_csv(f,columns=csv_columns,header=False,index=False)
                    if verbose:print("*CSV file appended>{0}.png, to location{1}{2}{3}".format(f_name,'\n',full_path,'\n')) 
            else:
                data.to_csv(full_path,columns=csv_columns,index=False)
                if verbose:print("*CSV file written->{0}.png, to location{1}{2}{3}".format(f_name,'\n',full_path,'\n')) 
        
        elif file_ext == 'lineage':
            full_path=os.path.join(save_dir,sub_dir_1,sub_dir_2, sub_dir_3,sub_dir_4,sub_dir_5,sub_dir_6,full_f_name)
            with open(full_path, 'w',newline='') as f:
                for nt_lineage in data.values():
                    f.write("{0}\n".format(','.join(map(str, nt_lineage))))
            if verbose:print("*lineage file written->{0}, to location{1}{2}{3}".format(f_name,'\n',full_path,'\n')) 
            
        elif file_ext == 'stl':  #from stl import mesh
            full_path=os.path.join(save_dir,sub_dir_1,sub_dir_2, sub_dir_3,sub_dir_4,sub_dir_5,sub_dir_6,full_f_name + '.stl')
            data.save(full_path)
            if verbose:print("*stl-file written->{0}, to location{1}{2}{3}".format(f_name,'\n',full_path,'\n')) 
            
        return os.path.join(save_dir,sub_dir_1,sub_dir_2,sub_dir_3,sub_dir_4,sub_dir_5,sub_dir_6)
    
    
    def get_save_location(self):
        return self.save_data(data=None,file_ext='xxx')
        
    def get_load_location(self):
        return os.path.join(self.d_load_info['load_dir'],self.d_load_info['sub_dir_1'],self.d_load_info['sub_dir_2'])
    
    def get_root_save_location(self):
        temp = self.d_save_info
        self.set_save_info_to_root()
        path= self.save_data(data=None,file_ext='xxx')
        self.d_save_info=temp
        return path
        
    def __repr__(self):
        for k in self.__dict__.keys(): print(k,"=",self.__dict__[k])
        return