'''
Created on 27 Mar 2019

@author: wimth  
'''
from fileHandler import FileHandler
from param_XML import Param_xml
from preProcess import preprocess_image
from spheres_DT import spheres_DT
from lineager_feeder import generate_lineager_input_files
from CTC_Tracking import generate_CTC_TRA_files
from mpacts_input_generator import mpacts_input_generator
from fig_z_trajectory import fig_z_trajectory
import os,datetime,sys
import shutil
import pprint
from pathlib import Path
pp = pprint.PrettyPrinter(indent=4)

os.system('hostname')




d_execution_blocks = {'1':'preprocessing',
                      '2':'spheresDT',
                      '3':'lineager_feeder',
                      '4':'CTC_TRA',
                      '5':'mpacts_input_generator',
                      '6':'fig_z_trajectory'
                      }

def check_params_and_data():
    
    return True




if __name__ == '__main__':
    
    
    def read_parms():
           
        l_execution_blocks = param_xml.get_value('l_execution_blocks', ['flowcontrol']) 
        ix_restart = param_xml.get_value('ix_restart', ['flowcontrol']) 
        l_execution_blocks = l_execution_blocks if (not ix_restart or ix_restart==0) else l_execution_blocks[ix_restart:]
        
        img_raw_file = param_xml.get_value('img_raw_file', ['paths'])
        img_exterior_outline = param_xml.get_value('img_exterior_outline', ['paths'])

        IMG_RAW_DIR = str(img_raw_file.parent)
        IMG_RAW_FILE= img_raw_file.name
        OUTPUT_FOLDER_ROOT = param_xml.get_value('OUTPUT_FOLDER_ROOT', ['paths'])
        ic_timestamp_subfolder = param_xml.get_value('ic_timestamp_subfolder', ['paths'])
        
        return l_execution_blocks, IMG_RAW_DIR,IMG_RAW_FILE,OUTPUT_FOLDER_ROOT,ic_timestamp_subfolder,img_exterior_outline
    
    def initialize_filehandler():
        
        filehandler = FileHandler()
        filehandler.d_save_info['save_dir']= OUTPUT_FOLDER_ROOT
        #print("DEBUG1: filehandler {0},  OUTPUT_FOLDER_ROOT {1}, filehandler.get_save_location() {2}".format(filehandler.__repr__(), OUTPUT_FOLDER_ROOT,filehandler.get_save_location()))
        
        if ic_timestamp_subfolder:
            if IMG_RAW_FILE:
                filehandler.d_save_info['sub_dir_1']=str(datetime.datetime.now()).replace(":","_").replace(" ","_") + "_" + IMG_RAW_FILE[0:-4]
            else:
                filehandler.d_save_info['sub_dir_1']=str(datetime.datetime.now()).replace(":","_").replace(" ","_") 
                
        filehandler.set_save_info_root()
        #print("DEBUG2: filehandler {0},  OUTPUT_FOLDER_ROOT {1}, filehandler.get_save_location() {2}".format(filehandler.__repr__(), OUTPUT_FOLDER_ROOT,filehandler.get_save_location()))
        filehandler.d_load_info['load_dir']=IMG_RAW_DIR
        filehandler.take_snapshot()
        #print("DEBUG3: filehandler {0},  OUTPUT_FOLDER_ROOT {1}, filehandler.get_save_location() {2}".format(filehandler.__repr__(), OUTPUT_FOLDER_ROOT,filehandler.get_save_location()))
        print('SDT_MAIN.py:all info will be written to : ',filehandler.get_save_location())
        
        return filehandler
    
    def update_filehandler_with_f_name(filename_raw):
        if not IMG_RAW_FILE or ic_timestamp_subfolder: #if input file name is not filled in , create a subfolder for each image processed
            filehandler.extend_save_info(extra_dir_1=filename_raw[0:-4],from_root=True,set_root=True)
        filehandler.store_f_name('IMG_RAW_FILE', filename_raw)
        filehandler.d_load_info['f_name']=filename_raw

        filehandler.store_f_name('img_exterior_outline', img_exterior_outline)


        filehandler.take_snapshot()
        
        return filehandler
    
    def load_raw_image():
        filehandler.d_load_info['load_dir'] = param_xml.get_value('IMG_RAW_DIR',['paths'],  dtype='string')
        filehandler.d_load_info['f_name'] = param_xml.get_value('IMG_RAW_FILE',['paths'],  dtype='string')
        return filehandler.load_tif()
    

    print('STEP0-SETUP------------------------------------------------------------------------------------------------------------')
                
    param_xml = Param_xml.get_param_xml(sys.argv,l_main_keys = ['body','MAIN'],verbose=True)
    # file_param, param_xml,param_xml.l_main_keys = read_param_xml_file()
    l_execution_blocks,IMG_RAW_DIR,IMG_RAW_FILE,OUTPUT_FOLDER_ROOT,ic_timestamp_subfolder,img_exterior_outline = read_parms()
    
    filehandler = initialize_filehandler()
    os.makedirs(filehandler.get_save_location(),exist_ok=True)
    #shutil.copyfile(str(param_xml.file), os.path.join(filehandler.get_save_location(),param_xml.file.name)) #backup param file
    
    print('checks done, proceeding...') if check_params_and_data() else print('checks failed. stopping executing')
    

    for filename_raw in Path(IMG_RAW_DIR).glob('**/*.tif'):
    # for filename_raw in os.listdir(IMG_RAW_DIR):
    #     if filename_raw.endswith(".tif") or filename_raw==IMG_RAW_FILE:
            
        if IMG_RAW_FILE and filename_raw.name != IMG_RAW_FILE:continue  #if input file is given , only process this one file
        
        update_filehandler_with_f_name(filename_raw.name)
        
        for execution_block in l_execution_blocks:
            if d_execution_blocks.get(execution_block)=='preprocessing':
                print('STEP1-PREPROCESSING---------------------------------------------------------------------------------------------------')
                preprocess_image(param_xml,filehandler)
                
            elif d_execution_blocks.get(execution_block)=='spheresDT':
                print('STEP2-SPHERES_DT------------------------------------------------------------------------------------------------------')
                spheres_DT(param_xml,filehandler)
                
            elif d_execution_blocks.get(execution_block)=='lineager_feeder':
                print('STEP3-LINEAGER_FEEDER------------------------------------------------------------------------------------------------------')
                generate_lineager_input_files(param_xml, filehandler)
                
            elif d_execution_blocks.get(execution_block)=='CTC_TRA':
                print('STEP4-CTC_TRA------------------------------------------------------------------------------------------------------')
                generate_CTC_TRA_files(param_xml, filehandler)
                
            elif d_execution_blocks.get(execution_block)=='mpacts_input_generator':
                print('STEP5-mpacts_input_generator-----------------------------------------------------------------------------------------------------')
                mpacts_input_generator(param_xml, filehandler)
                
            elif d_execution_blocks.get(execution_block)=='fig_z_trajectory':
                print('STEP6-fig_z_trajectory-----------------------------------------------------------------------------------------------------')
                fig_z_trajectory(param_xml)
                
            else:
                print('execution block not recognized')
            filehandler.set_save_info_to_root()
        #<--next executionstep
        filehandler.pop_save_info(set_root=True) #remove filename from root
    #<--next tif
    
    # if os.popen('hostname').read().startswith('DESKTOP'):ipdb.pm()
    
