"""
Created on 01 july 2020
aggregate geometry data from all replicates
@author: wimth
"""
import numpy as np
import pandas as pd
from helper_pandas import excel_to_dds
import sys,re
from param_XML import Param_xml
from shutil import copyfile
from pathlib import Path


verbose=True

def dds_get_value(label_ds,k1,k2):
	try:
		return dds[label_ds][k1][k2]
	except KeyError as ex:
		return 'not in dds : {0}'.format(k1)

def read_parms():
	param_xml = Param_xml.get_param_xml(sys.argv, l_main_keys=['body', 'agg_seg_data'],verbose=True)  # param file must be passed as first argument

	input_dds = param_xml.get_value('input_dds', ['paths'])
	root_folder = param_xml.get_value('root_folder', ['paths'])
	regex_file_select = param_xml.get_value('regex_file_select', ['paths'])
	regex_file_select2 = param_xml.get_value('regex_file_select2', ['paths'])
	regex_file_select3 = param_xml.get_value('regex_file_select3', ['paths'])
	regex_file_select4 = param_xml.get_value('regex_file_select4', ['paths'])
	regex_file_select5 = param_xml.get_value('regex_file_select5', ['paths'])
	output_folder = param_xml.get_value('output_folder', ['paths'])


	return input_dds, root_folder, regex_file_select,regex_file_select2, regex_file_select3, regex_file_select4, regex_file_select5,output_folder




def process_SEG_file(file_i,tcID,tc_info):
	df_cell = pd.read_csv(file_i,header=0,index_col=False,keep_default_na=True) 
	cond1 = df_cell['GT_cell'] != "SUM"
	df_cell = df_cell[cond1]  #all the details of all cells
	df_cell.insert(0,'tcID',tcID)

	d_seg = {}
	d_seg['tcID'] = [tcID]
	d_seg['type'] = tc_info['score']
	d_seg['score'] =  df_cell["SEG3D"].mean()

	d_dice = {}
	d_dice['tcID'] = [tcID]
	d_dice['type'] = 'dice'
	d_dice['score'] =  df_cell["dice"].mean()

	return [df_cell, pd.DataFrame(d_seg),pd.DataFrame(d_dice)]

def process_SEG2D_file(file_i,tcID,tc_info):
	try:
		df_cell = pd.read_csv(file_i,header=0,index_col=False,keep_default_na=True) 
	except Exception as e:
		print("No SEG2D validationfile...",end="")
		return [None,None,None]
	
	cond1 = df_cell['GT_label'] != "SUM"
	df_cell = df_cell[cond1]  #all the details of all cells
	df_cell.insert(0,'tcID',tcID)

	d_seg2D = {}
	d_seg2D['tcID'] = [tcID]
	d_seg2D['type'] = 'SEG2D'
	d_seg2D['score'] =  df_cell["SEG_score"].mean()

	d_SEG2D_corr = {}
	d_SEG2D_corr['tcID'] = [tcID]
	d_SEG2D_corr['type'] = 'SEG2D_bg_ignored'
	d_SEG2D_corr['score'] =  df_cell["SEG_score_background_ignored"].mean()

	return [df_cell, pd.DataFrame(d_seg2D),pd.DataFrame(d_SEG2D_corr)]

def process_tracking_file(file_i,tcID,tc_info):
	df_tra_detail = pd.read_csv(file_i,header=0,index_col=False,keep_default_na=True) 
	row_last = df_tra_detail[df_tra_detail['timestep'] == 'TOTAL'].index[0] 
	df_tra_detail = df_tra_detail[0:row_last]
	d_dtype =  {name : int for name in df_tra_detail.columns}
	df_tra_detail = df_tra_detail.astype(d_dtype)
	df_tra_detail.insert(0,'tcID',tcID)

	d_tra = {}
	d_tra['tcID'] = [tcID]  * 3
	
	tp_sum = df_tra_detail['TP'].sum()
	fn_sum = df_tra_detail['FN'].sum()
	fp_sum = df_tra_detail['FP'].sum()
	d_tra.setdefault('type', []).append('accuracy')
	d_tra.setdefault('score', []).append(tp_sum/(tp_sum + fn_sum + fp_sum))
	d_tra.setdefault('type', []).append('sensitivity')
	d_tra.setdefault('score', []).append(tp_sum/(tp_sum + fp_sum))
	d_tra.setdefault('type', []).append('precision')
	d_tra.setdefault('score', []).append(tp_sum/(tp_sum + fn_sum))

	return [df_tra_detail, pd.DataFrame(d_tra)]

def process_tracking_metrics(file_i,tcID,tc_info):
	df_cell4D = pd.read_csv(file_i,header=0,index_col=False,keep_default_na=True) 
	df_cell4D.insert(0,'tcID',tcID)
	df_cell4D.insert(1,'cell_stage',dds['tc'][tcID]['cell_stage'])

	return [pd.DataFrame(df_cell4D)]

def process_TRA_file(file_i,tcID,tc_info):
	s_file = open(file_i).read()
	d_tra = {}
	d_tra['tcID'] = [tcID]
	d_tra['type'] = tc_info['score']
	d_tra['score'] = float(re.search(r"(?:TRA measure: )(\d*.\d*)",s_file).group(1)) #TRA measure: 0.990486

	return pd.DataFrame(d_tra)

def enrich_SEG_REL_score(df_agg_scores):
	d_seg_rel = {}
	for index, row in df_agg_scores.iterrows():
		if row['type'] != 'SEG3D':
			continue
		tcID_rel = dds['tc'][row['tcID']]['tcID_relative']
		if tcID_rel:
			cond = df_agg_scores['tcID'] ==tcID_rel
			cond2 = df_agg_scores['type'] =='SEG3D'
			
			SEG_relative = df_agg_scores[cond & cond2 ].get('score')
			if not SEG_relative.empty:
				d_seg_rel.setdefault('tcID',[]).append(row['tcID'])
				d_seg_rel.setdefault('type',[]).append('SEG_REL')
				d_seg_rel.setdefault('score',[]).append((row['score'] / float(SEG_relative)))

	return pd.concat([df_agg_scores,pd.DataFrame(d_seg_rel)],ignore_index=True)

#parms
input_dds, root_folder,regex_file_select,regex_file_select2, regex_file_select3 , regex_file_select4,regex_file_select5, output_folder = read_parms()
output_folder.mkdir(parents=True,exist_ok=True)

#STEP 1 : read in datastructures for mapping
dds = excel_to_dds(str(input_dds),sheet_name=None) #dds['tc'] : d_tcID_info

#STEP 2 : aggregate segmentation score data
writer = pd.ExcelWriter(output_folder/"agg_SEG_TRA_data.xlsx")
l_df_agg_scores = []
l_df_seg_detail = []
l_df_seg2D_detail = []
l_df_tracking_detail = []
l_df_tracking_metrics = []
l_d_tcID_info = []
d_error = {}
for tcID_i, tc_info_i in dds['tc'].items():
	if tc_info_i['aggregate'] != 'yes':
		continue
	try:
		if tc_info_i['score'] == 'SEG3D':
			#SEG and SEG_REL score
			file_i = root_folder / regex_file_select.replace("<tcID>",tcID_i)
			l_df = process_SEG_file(file_i,tcID_i,tc_info_i)
			l_df_seg_detail.append(l_df[0])
			l_df_agg_scores.append(l_df[1])
			l_df_agg_scores.append(l_df[2])
			# #SEG2D and SEG2D_bg_corr
			# file_i_2 = root_folder / regex_file_select5.replace("<tcID>",tcID_i)
			# l_df = process_SEG2D_file(file_i_2,tcID_i,tc_info_i)
			# if l_df:
			# 	l_df_seg2D_detail.append(l_df[0])
			# 	l_df_agg_scores.append(l_df[1])
			# 	l_df_agg_scores.append(l_df[2])

		elif tc_info_i['score'] == 'TRA':
			#tracking detail
			file_i = root_folder / regex_file_select2.replace("<tcID>",tcID_i)
			l_df = process_tracking_file(file_i,tcID_i,tc_info_i)
			l_df_tracking_detail.append(l_df[0])
			l_df_agg_scores.append(l_df[1])
			#TRA score
			file_i = root_folder / regex_file_select3.replace("<tcID>",tcID_i)
			df_TRA = process_TRA_file(file_i,tcID_i,tc_info_i)
			l_df_agg_scores.append(df_TRA)
			#tracking extra metrics
			file_i = root_folder / regex_file_select4.replace("<tcID>",tcID_i)
			l_df = process_tracking_metrics(file_i,tcID_i,tc_info_i)
			l_df_tracking_metrics.append(l_df[0])
			

		if verbose:print('{0} is appended'.format(file_i))
	except:
		d_error.setdefault('tc_error', []).append(tcID_i)
		print (file_i," is skipped ! check the file...")

df_seg_detail = pd.concat(l_df_seg_detail,ignore_index=True)
df_seg2D_detail = pd.concat(l_df_seg_detail,ignore_index=True)
df_tracking_detail = pd.concat(l_df_tracking_detail,ignore_index=True)
df_tracking_metrics = pd.concat(l_df_tracking_metrics,ignore_index=True)
df_agg_scores = pd.concat(l_df_agg_scores,ignore_index=True)
df_error = pd.DataFrame(d_error)
for key in dds['tc'].keys():
	dds['tc'][key]['tcID'] = key 
df_tcinfo = pd.DataFrame.from_dict(dds['tc'],orient='index')

#STEP 3 : calculate SEG_REL
df_agg_scores = enrich_SEG_REL_score(df_agg_scores)


#STEP4 write output excel
with pd.ExcelWriter(output_folder/"agg_seg_data.xlsx") as writer:
	df_agg_scores.to_excel(writer, sheet_name='agg_scores',index=False)
	df_tracking_detail.to_excel(writer, sheet_name='tracking detail',index=False)
	df_tracking_metrics.to_excel(writer, sheet_name='tracking metrics',index=False)
	df_seg_detail.to_excel(writer, sheet_name='SEG detail',index=False)
	df_seg2D_detail.to_excel(writer, sheet_name='SEG2D detail',index=False)
	df_tcinfo.to_excel(writer, sheet_name='tc_info',index=False)
	df_error.to_excel(writer, sheet_name='errors',index=False)
	if verbose:print("excel output : {0}".format(str(output_folder/"agg_seg_data.xlsx")))

print('END')