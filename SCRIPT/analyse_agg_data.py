"""
Created on 29 January 2021
This is a python wrapper around several R modules that will analyse the aggregated data obtained from different testcases  
this module will enable parameter resolving and passing it on the R modules 
It will also create the directory structure for all the files (each R module will get its own submap)  
For now the same parameters will be passed on to each module to keep it generic  (input file and outputfolder
@author: wimth
"""

import pandas as pd
import sys, os, re
from param_XML import Param_xml
from pathlib import Path
from shutil import copy
import subprocess
import logging
verbose=True

def read_parms():
	param_xml = Param_xml.get_param_xml(sys.argv, l_main_keys=['body'],verbose=True) 
	param_xml.add_prms_to_dict(l_keys=['analyse_agg_data'])
	param_xml.repr_prm()
	return param_xml.prm

# MAIN-----------------------------------------------
prm = read_parms()
prm['output_folder'].mkdir(parents=True,exist_ok=True)
logging.basicConfig(filename=(prm['output_folder'] / "log.txt"), filemode='w', format='%(name)s - %(levelname)s - %(message)s',level=logging.INFO)


for R_module in prm['l_R_modules']:
	# Define command and arguments
	command ='Rscript'
	path2script = R_module #'path/to your script/max.R'

	
	output_folder_script = prm['output_folder']/Path(R_module).stem
	output_folder_script.mkdir(parents=True,exist_ok=True)  #sub folder for each R_script

	# Variable number of args in a list
	args = [prm['input_file'],str(output_folder_script)+os.sep]

	# Build subprocess command
	cmd = [command, path2script] + args

	# check_output will run the command and store to result
	print(f" Executing cmd = {cmd}",flush=True)
	R_output = subprocess.check_output(cmd, universal_newlines=True) 
	logging.info(R_output)
	print("Flushing collected R_output==>",R_output,flush=True)
	print(f" End of Execution of cmd = {cmd}",flush=True)