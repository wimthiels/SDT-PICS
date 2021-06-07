'''
Created on 28 Mar 2019

@author: wimth
'''

import numpy as np        
from helper_functions import slice_tif,plot3D_stack_with_overlay,convert_16bit_to_RGB
from datamodel import Stack,Featuretracker,Cell4D,Label,Param,Cluster
from io import StringIO
import time
import pandas as pd
import os
import shutil
#import linecache
import tracemalloc
from skimage.exposure import adjust_gamma
#tracking memory allocation  : https://stackoverflow.com/questions/552744/how-do-i-profile-memory-usage-in-python



def spheres_DT(param_xml,filehandler):
	#>>>>>FUNCTIONS>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	
	# def display_top(snapshot, key_type='lineno', limit=3):
	# 	if prm['verbose']:print("###tracemalloc report for SpheresDT.py")
	# 	snapshot = snapshot.filter_traces((
	# 		tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
	# 		tracemalloc.Filter(False, "<unknown>"),
	# 	))
	# 	top_stats = snapshot.statistics(key_type)

	# 	print("Top %s lines" % limit)
	# 	for index, stat in enumerate(top_stats[:limit], 1):
	# 		frame = stat.traceback[0]
	# 		# replace "/path/to/module/file.py" with "module/file.py"
	# 		filename = os.sep.join(frame.filename.split(os.sep)[-2:])
	# 		print("#%s: %s:%s: %.1f KiB"
	# 			  % (index, filename, frame.lineno, stat.size / 1024))
	# 		line = linecache.getline(frame.filename, frame.lineno).strip()
	# 		if line:
	# 			print('    %s' % line)

	# 	other = top_stats[limit:]
	# 	if other:
	# 		size = sum(stat.size for stat in other)
	# 		print("%s other: %.1f KiB" % (len(other), size / 1024))
	# 	total = sum(stat.size for stat in top_stats)
	# 	print("Total allocated size: %.1f KiB" % (total / 1024))
	

	def initialize_class_variables():
		
		if prm['verbose']:print("#initialize_class_variables>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")

		Param(param_xml,filehandler=filehandler)
		Stack.l_RGBA_cellview=[]
		Label.create_indicator_labels()
		
		return
	
	 
	def append_summary_of_stack_txt_repr(stack):
		if not stack:return
		if prm['verbose']:print("###saving a txt-summary of the stack repr (append)")
		
		output = StringIO()
		stack.repr(detail_cells=True,output=output)
	
		filehandler.d_save_info['f_name']='debug_logging_repr'
		filehandler.save_data(output.getvalue(),file_ext='txt',verbose=False)
		
		return
	
	def append_summary_of_4Dcell_txt_repr():
		if prm['verbose']:print("###saving a txt-summary of the 4D repr (append)")
		
		output = StringIO()
		for cell4D in Cell4D.l_cell4D:
			cell4D.repr(output=output)
	
		filehandler.d_save_info['f_name']='summary'
		filehandler.save_data(output.getvalue(),file_ext='txt',verbose=False)
		
		return
	
	def append_labelling_chronology_to_txt_log():
		if prm['verbose']:print("###saving the labelling chronology to the txt summary (append)")
		
		output = StringIO()
		for label,l_labelling_log in feat_tracker.d_label_events.items():
			print('',file=output)
			print('Start labelling_log for ',label,' **********',file=output)
			for ix,event in enumerate(l_labelling_log):
				print(str(ix+1),'->',event,file=output)
	
		filehandler.d_save_info['f_name']='summary'
		filehandler.save_data(output.getvalue(),file_ext='txt',verbose=False)
		
		return
	
	def add_exterior_mask():
		if not prm['USE_EXTERIOR_OVERLAY']:return
		if prm['verbose']:print("###add a mask to shield off the exterior")
		
		filehandler.d_load_info['load_dir'] = prm['DIR_EXTERIOR_MASK'] if prm['DIR_EXTERIOR_MASK'] else os.path.join(filehandler.get_root_save_location(),'002_preprocessing')
		filehandler.d_load_info['f_name'] = prm['FILE_EXTERIOR_MASK'] if prm['FILE_EXTERIOR_MASK'] else 'exterior_mask.tif'

		stack.add_exterior_mask(slice_tif(filehandler=filehandler,ix_slice_time=[nb_stack-1,nb_stack],verbose=False),mask=True)
		
		return
	
	
	def save_dist_transf_carved():
		if not prm['SAVE_EXTRA_OUTPUT']: return
			
		filehandler.d_save_info['f_name'] = 'distanceTransform_withmask';
		filehandler.save_data(stack.dist_transf_carved,verbose=False)
		
		return
	
	def create_stack():
		if prm['verbose']:print("##create a stack ----------------------------------------------------------------------------timestep{0}".format(nb_stack))
		img = slice_tif(na_tif=a_4D_membrane,ix_slice_time=[nb_stack-1,nb_stack],verbose=False)
		img_bin = np.zeros(img.shape,'bool')
		img_bin[img > 0] = True 
		stack = Stack(img = img_bin,timestep=nb_stack,feat_tracker=feat_tracker,name=filehandler.d_load_info['f_name'])
		
		return stack
	
	def save_extra_output():
		filehandler.d_save_info['sub_dir_3']= ''
		if not prm['SAVE_EXTRA_OUTPUT']:return
		
		if prm['verbose']:print("###saving OUTPUT for stack {0}".format(nb_stack))
		if prm['verbose']:print("####plot the img data of stack{0}".format(nb_stack))
		filehandler.d_save_info['f_name'] = 'fransceca'
		#plot3D_stack_with_overlay(stack.img, d_save_info=filehandler.d_save_info,intensity_threshold=INTENSITY_THRESHOLD)
		
		if prm['verbose']:print("####save the distance transform after carving out of spheres")
		np.where(stack.dist_transf_carved<0,stack.dist_transf_carved*(-1),stack.dist_transf_carved)
		filehandler.d_save_info['f_name'] = 'stack.dist_transf_carved';
		filehandler.save_data(stack.dist_transf_carved,verbose=False)
		
	#     if verbose:print('####make a 3Dplot of the cells detected in stack{0}'.format(nb_stack))
	#     filehandler.d_save_info['f_name'] = 'stack.save_3Dplot'
	#     stack.save_3Dplot(filehandler.d_save_info)
			
		return
	
	def save_compare_tifs():
		if prm['verbose']:print('####save some stacks for comparison (img and membrane contour')
		
		filehandler.d_load_info['load_dir']= 'base'
		filehandler.d_load_info['f_name'] = '4D_membrane_contour'
		filehandler.d_save_info['f_name'] = 'Compare_3D_membrane_contour' + '_ST' + str(nb_stack)
		slice_tif(ix_slice_time=[nb_stack-1,nb_stack],filehandler=filehandler,verbose=False)
		
		filehandler.d_load_info['f_name'] = 'raw4D'
		filehandler.d_save_info['f_name'] = 'Compare_raw3D' + '_ST' + str(nb_stack)
		slice_tif(ix_slice_time=[nb_stack-1,nb_stack],filehandler=filehandler,verbose=False)
		
		return
	
	def read_validation_nuclei():
		if not prm['VALIDATE_CELLS']:return None
		if prm['verbose']:print("#load the Cell nuclei for validation")
	#     if FILE_PRESET in ['TL1_with_preprocessing','TL1_no_preprocessing']:
	#     df_val = pd.read_csv(Param_XML.FILE_VALIDATION,delimiter=',')
	# #     else:
		df_val = pd.read_csv(prm['FILE_VALIDATION'],delimiter='\t')
		df_val["x"] = pd.to_numeric(df_val["x"],downcast='integer')
		df_val["y"] = pd.to_numeric(df_val["y"],downcast='integer')
		df_val["z"] = pd.to_numeric(df_val["z"],downcast='integer')
		
		return df_val
	
	def finish_last_stack(stack):
		if prm['verbose']:print('##finishing last stack')
		
		if prm['VALIDATE_CELLS']:
			stack.validate_clusters(df_nuclei,verbose=prm['verbose'])
		feat_tracker.update_features_v2(stack,verbose=prm['verbose'])
		
		stack.construct_RGBA_clusterview(filehandler,temp_folder=True,
										show_cells_only=prm['SHOW_CELLS_ONLY'],
										df_validation = df_nuclei,
										validation_start_time = prm['validation_start_time'],
										verbose=prm['verbose'])
		
		stack.update_lineage_names_cell_division()
		stack.write_excel_summary_3D(filehandler,verbose=prm['verbose'])
		stack.write_excel_validation(filehandler,print_total=True,verbose=prm['verbose'])
		
		if prm['DBG']:append_summary_of_stack_txt_repr(stack)
		
		return 
	
	
	def get_extra_channel():
		
		if prm['DIR_EXTRA_CHANNEL']:
			filehandler.d_load_info['load_dir'] = prm['DIR_EXTRA_CHANNEL'] 
			filehandler.d_load_info['f_name'] = prm['FILE_EXTRA_CHANNEL']
			filehandler.d_save_info['pre_f_name'] = ''
			a_extra_channel = filehandler.load_tif()
		else:
			filehandler.d_load_info['load_dir'] = prm['IMG_RAW_DIR']
			filehandler.d_load_info['f_name'] = prm['IMG_RAW_FILE']
			filehandler.d_save_info['pre_f_name'] = ''
			a_extra_channel = convert_16bit_to_RGB(filehandler.load_tif(storage_name='IMG_RAW_FILE'),clip_not_scale=True)

		if len(a_extra_channel.shape)==4:  #turn stack into timelapse
			z,y,x,rgba = a_extra_channel.shape
			a_extra_channel = a_extra_channel.reshape((1,z,y,x,rgba)) 
		
		
		
		if prm['l_stack_raw'] == 'all':
			a_slice_tif =  slice_tif(na_tif=a_extra_channel,ix_slice_time= [prm['lnb_stack'][0] - 1, prm['lnb_stack'][-1]],verbose=False,RGB=True,trim_1D=False)
		else:
			a_slice_tif =  slice_tif(na_tif=a_extra_channel,ix_slice_time= [i-1 for i in prm['l_stack_raw']],verbose=False,RGB=True,trim_1D=False)

		make_brighter=True #causes memory crash for big timelapses,so go stack by stack
		if make_brighter:
			for ix,i in enumerate(a_slice_tif):
				a_slice_tif[ix,:] = adjust_gamma(i, gamma=0.5)

		return a_slice_tif
	
	def load_extra_channel():
		
		if not prm['DIR_EXTRA_CHANNEL']:
			return None
		
		filehandler.d_load_info['load_dir'] = prm['DIR_EXTRA_CHANNEL'] 
		filehandler.d_load_info['f_name'] = prm['FILE_EXTRA_CHANNEL']
		filehandler.d_save_info['pre_f_name'] = ''
		return slice_tif(ix_slice_time= [prm['lnb_stack'][0] - 1, prm['lnb_stack'][-1]], filehandler=filehandler,verbose=False,RGB=True,trim_1D=False)
		
	
	def read_parms(param_xml):
		
		prm = {}
		param_xml.l_main_keys = ['body','spheresDT']
		
		prm['DIR_MEMBRANE']  = param_xml.get_value('DIR_MEMBRANE',['paths'])
		prm['FILE_MEMBRANE'] = param_xml.get_value('FILE_MEMBRANE',['paths'])
		prm['DIR_EXTERIOR_MASK'] = param_xml.get_value('DIR_EXTERIOR_MASK',['paths'])
		prm['FILE_EXTERIOR_MASK'] = param_xml.get_value('FILE_EXTERIOR_MASK',['paths'])
		prm['FILE_VALIDATION'] = param_xml.get_value('FILE_VALIDATION',['paths'])
		prm['DIR_EXTRA_CHANNEL'] = param_xml.get_value('DIR_EXTRA_CHANNEL',['paths'])
		prm['FILE_EXTRA_CHANNEL']= param_xml.get_value('FILE_EXTRA_CHANNEL',['paths'])
		prm['lnb_stack'] = param_xml.get_value('lnb_stack',['flowcontrol'])
		prm['lnb_GT_stack'] = param_xml.get_value('lnb_GT_stack',['flowcontrol'])
		prm['validation_start_time'] = param_xml.get_value('validation_start_time',['flowcontrol'])
		prm['reset_seeds'] = param_xml.get_value('reset_seeds',['flowcontrol'])
		
		param_xml.l_main_keys = ['body','MAIN']
		prm['img_raw_file'] = param_xml.get_value('img_raw_file',['paths'])
		prm['IMG_RAW_DIR']  = filehandler.d_load_info['load_dir'] if filehandler.d_load_info['load_dir'] else str(img_raw_file.parent)
		prm['IMG_RAW_FILE'] = filehandler.d_load_info['f_name'] if filehandler.d_load_info['f_name'] else img_raw_file.name
		
		param_xml.l_main_keys = ['body','preprocessing']
		prm['l_stack_raw'] = param_xml.get_value('l_stack_number',['flowcontrol'])

		if not prm['FILE_VALIDATION']:
			try:
				path_to_file= os.path.join(prm['IMG_RAW_DIR'],'validation_files',filehandler.d_name_f_name['IMG_RAW_FILE'][0:-4]+'.csv')
				fh = open(path_to_file, 'r')
				prm['FILE_VALIDATION'] = path_to_file
				# Store configuration file values
			except FileNotFoundError:
				pass # Keep preset values
		
		return prm
	
	def set_process_parameters_unparameterized(prm):

		prm['verbose']=False
		prm['SAVE_COMPARE_TIFS'] = False
		prm['SAVE_EXTRA_OUTPUT'] = False
		prm['INTENSITY_THRESHOLD'] = 100
		prm['USE_EXTERIOR_OVERLAY'] = True
		
		prm['VALIDATE_CELLS'] = True if prm['FILE_VALIDATION'] else False
		prm['SHOW_CELLS_ONLY'] = False if prm['VALIDATE_CELLS'] else True

		prm['DBG'] = False
		
		return 
	
	def load_membrane_data():
		
		filehandler.d_load_info['load_dir'] = prm['DIR_MEMBRANE'] if prm['DIR_MEMBRANE'] else os.path.join(filehandler.get_root_save_location(),'002_preprocessing')
		filehandler.d_load_info['f_name'] = prm['FILE_MEMBRANE'] if prm['FILE_MEMBRANE'] else 'membranes.tif'
		a_4D = filehandler.load_tif(verbose=prm['verbose'])
		if len(a_4D.shape) == 3:
			a_4D = a_4D.reshape(1,*a_4D.shape)

		return a_4D
	
	def remove_temp_files():
		filehandler.extend_save_info(extra_dir_1='003_spheresDT',extra_dir_2='temp',from_root=True)
		shutil.rmtree(filehandler.get_save_location())
		
		return
	
	############################################################################################################################
	#tracemalloc.start()
	tic = time.time() 
	#set parms
	prm = read_parms(param_xml)
	set_process_parameters_unparameterized(prm)
	
	#fill in general load and save info
	filehandler.extend_save_info(extra_dir_1='003_spheresDT',from_root=True)
	filehandler.clear_save_dir()

	initialize_class_variables()
	
	a_4D_membrane = load_membrane_data()
	Param.set_img_dimensions(a_4D_membrane)
	
	nb_inconsistencies = 0
	
	df_nuclei =read_validation_nuclei()
	
	feat_tracker = Featuretracker(nb_timesteps = a_4D_membrane.shape[0])
	
	if prm['lnb_stack']=='all':prm['lnb_stack'] = range(1,len(a_4D_membrane)+1) 
	if prm['lnb_GT_stack']=='all':prm['lnb_GT_stack'] = range(1,len(a_4D_membrane)+1)  

	
	for nb_stack in prm['lnb_stack']:
		print("Processing stack {}".format(nb_stack))
		#read/initialize--------------------- 
		start_time = time.time()
		filehandler.d_save_info['nb_stack'] = nb_stack
		stack = create_stack()
		prev_stack = Stack.prev_stack
		add_exterior_mask()
		stack.add_distance_transform(verbose=prm['verbose'])
		if prm['SAVE_EXTRA_OUTPUT']:save_dist_transf_carved()
		ic_GT_corr= True if (nb_stack + prm['validation_start_time']-1) in prm['lnb_GT_stack'] else False
		
		#SPHERES/CLUSTERS--------------------- 
		  #seeding
		if ic_GT_corr:
			stack.seed_with_validation_file(df_nuclei,prm['validation_start_time'],verbose=prm['verbose'])
		elif prev_stack:
			stack.seed_with_cells_from_prev_stack(verbose=prm['verbose'])

		  #clustering
		stack.put_spheres_in_distance_transform(verbose=prm['verbose'])
		
		stack.reset_seeds(func=prm['reset_seeds'])
		
		  #labeling
		if not prev_stack:
			stack.label_clusters_of_first_stack(ic_GT_corr, verbose=prm['verbose'])
		else:
			stack.label_clusters_via_prev_stack_information(verbose=prm['verbose'])
		
		  #filtering
		stack.filter_clusters(ic_GT_corr,verbose=prm['verbose'])
		Cluster.print_stats("after filtering")
		
		#CELLS---------------------------------    
		stack.create_cells_based_on_cluster_labels(verbose=prm['verbose'])
		
		if prev_stack:
			if ic_GT_corr:
				stack.create_cells_through_cell_divisions_GT_corr(verbose=prm['verbose'])
				Cluster.print_stats("after GT cell division creation")
			else:
				stack.create_cells_through_cell_divisions(verbose=prm['verbose'])
			
		stack.filter_cells(verbose=prm['verbose'])
		
		#CELL4D--------------------------------   
		stack.create_cell4D(verbose=prm['verbose'])
		stack.update_cell4D(verbose=prm['verbose'])
		stack.filter_cell4D(verbose=prm['verbose'])
		
		
		#Finalize previous stack-------------------------------   
		if prev_stack:
			#check datamodel before writing output
			nb_inconsistencies += prev_stack.check_consistency_datamodel(verbose=prm['verbose'])
			
			#validation
			if prm['VALIDATE_CELLS']:
				prev_stack.validate_clusters(df_nuclei,prm['validation_start_time'],verbose=prm['verbose'])
				prev_stack.write_excel_validation(filehandler)
				
			#determine the lineage_names for outputting
			prev_stack.update_lineage_names_cell_division(verbose=prm['verbose'])
			
			#output
			prev_stack.construct_RGBA_clusterview(filehandler,temp_folder=True,
				show_cells_only=prm['SHOW_CELLS_ONLY'],df_validation = df_nuclei,
				validation_start_time = prm['validation_start_time'],
				verbose=prm['verbose'])
			prev_stack.write_excel_summary_3D(filehandler,verbose=prm['verbose'])
			feat_tracker.update_features_v2(prev_stack,verbose=prm['verbose'])
			if prm['DBG']:append_summary_of_stack_txt_repr(prev_stack)
		else:
			stack.write_excel_summary_3D(filehandler,initialize=True,verbose=prm['verbose'])
			
		if prm['SAVE_EXTRA_OUTPUT']:save_extra_output()
		if prm['SAVE_COMPARE_TIFS']:save_compare_tifs()                                                                                                                                                                                                                                          
	
		if prm['verbose']:print('elapsed time for stack =', time.time() - start_time)
		
		#CHECKS
		Stack.check_memory_allocation(verbose=prm['verbose'])
		# snapshot = tracemalloc.take_snapshot()
		# display_top(snapshot,limit=20)
	#<<<end of stack loop
	
	finish_last_stack(Stack.curr_stack)
	
	if prm['DBG']:
		append_summary_of_4Dcell_txt_repr()
		append_labelling_chronology_to_txt_log()
	feat_tracker.write_excel_summary_4D_v2(filehandler)
	Stack.save_4D_RGBA_clusterview(filehandler,nb_channels=1,extra_channel=get_extra_channel())
	#Stack.save_4D_RGBA_clusterview(filehandler,nb_channels=1)
	Stack.agg_temp_viz_files(filehandler,search_string='raw_clusterview')
	Stack.agg_temp_viz_files(filehandler,search_string='labelview')
	#Stack.save_4D_clusterview(filehandler)
	
	Stack.init_class_variables()
	print('{0} datamodel inconsistencies found across all stacks'.format(nb_inconsistencies))
	
	remove_temp_files()
	
	# snapshot = tracemalloc.take_snapshot()
	# display_top(snapshot,limit=10)
	toc = time.time()
	print('runtime_spheresDT = ', toc-tic)

	return