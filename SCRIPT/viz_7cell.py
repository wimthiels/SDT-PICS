'''
Created on 22 July 2020 
Vizualize results from 7cell geometry data (besides the R-graphs)
		
@author: wimth
'''


import numpy as np
import pandas as pd
import sys,os,re
import seaborn as sns
import matplotlib.pyplot as plt
from param_XML import Param_xml
verbose=True

def read_parms():
	param_xml = Param_xml.get_param_xml(sys.argv,l_main_keys = ['body','viz_7cell'],verbose=True) #param file must be passed as first argument
	
	input_file_excel  = param_xml.get_value('input_file_excel',['paths'])
	output_folder  = param_xml.get_value('output_folder',['paths'])
	
	return input_file_excel,output_folder

def compose_ca_matrix(df_cell_ca,select_by = 'first',t=1,values='contactarea_mean'):
	""" 
	select_by
			'time' : pick a certain timestep t. cells not present at that time will get a dummy entry
			'first' : for every cell the first timepoint it appears is selected. 
			'last':  for every cell the last timepoint it apppears is selected
	"""
	#select the data 
	if select_by == 'time':
		rowcond1 = (df_cell_ca['t']==t)
		df_data = df_cell_ca[rowcond1]

		# missing cells are added for consistency between graphs
		for missing_cellid_i in np.setdiff1d(pd.unique(df_cell_ca['cellID_1']),pd.unique(df_data['cellID_1'])):
			extra_row = df_data.iloc[1].copy()
			extra_row['cellID_1']=missing_cellid_i
			extra_row['cellID_2']=missing_cellid_i
			extra_row[values]=np.nan
			df_data = df_data.append(extra_row)

	elif select_by in ['first','last']:
		if select_by == 'last':
			df_cell_ca.sort_values(by='t',axis=0,ascending=False,inplace=True,na_position='last')
		else:
			df_cell_ca.sort_values(by='t',axis=0,ascending=True,inplace=True,na_position='last')
		cond_exclude = (df_cell_ca['contactarea_mean'].isnull())
		ser_cellID_1_first = df_cell_ca[~cond_exclude].groupby(['cellID_1','cellID_2']).first()['t']  #Series, group by field becomes the index

		df_data = pd.DataFrame()
		for index_i, t_i in ser_cellID_1_first.iteritems():
			rowcond1= (df_cell_ca['t']==t_i)
			rowcond2= (df_cell_ca['cellID_1']==index_i[0])
			rowcond3= (df_cell_ca['cellID_2']==index_i[1])
			rowcond = rowcond1 & rowcond2 & rowcond3
			print(df_cell_ca[rowcond])
			df_data = df_data.append(df_cell_ca[rowcond])

	df_matrix = df_data.pivot(index='cellID_1', columns='cellID_2', values=values)
	return df_matrix

def strip_axis_matplotlib(ax):
		""" strips tick labels, axis labels (not ticks) from 1 AxesSubplot.  in matplotlib
		This cannot strip colorbars (you need to hold on to the colorbar when creating it and then do "cbar.set_ticks([])"")"""

		ax.xaxis.set_ticklabels([]) #  removes ticklabels
		ax.yaxis.set_ticklabels([])

		ax.xaxis.set_label_text("") #  removes axis label text
		ax.yaxis.set_label_text("")

		if hasattr(ax, 'zaxis'): 
			ax.zaxis.set_ticklabels([])
			ax.zaxis.set_label_text("")
			
		legend = ax.get_legend()
		if legend:
			legend.remove()

		#ax.set_visible(False)
		
		return ax


def create_heatmap(df_matrix, title='Contactarea between cells', filename = 'heatmap_ca.png',center=None,cmap='Reds'):


	mask = np.zeros_like(df_matrix, dtype=np.bool)
	mask[np.triu_indices_from(mask)]= True

	#get max
	a_helper = np.array(df_matrix)
	np.fill_diagonal(a_helper, 0)
	vmax = np.nanmax(a_helper)
	vmin = np.nanmin(a_helper)

	f, ax = plt.subplots(figsize=(11, 15)) 

	heatmap = sns.heatmap(df_matrix, 
						mask = mask,
						square = True,
						linewidths = .5,
						cmap = cmap,
						cbar_kws = {'shrink': 0.6}, 
						# cbar_kws = {'shrink': .4, 
						# 					"ticks" : [-1, -.5, 0, 0.5, 1]},
						vmin = vmin, 
						vmax = vmax,
						annot = True,
						center=center,
						fmt='.0f',
						annot_kws = {"size": "xx-large"})
	#add the column names as labels
	#ax.set_title(title)
	ax.set_title("")
	ax.set_yticklabels(df_matrix.columns, rotation = 0)
	ax.set_xticklabels(df_matrix.columns)
	ax.set_ylabel('') 
	ax.set_xlabel('')
	for tick in ax.xaxis.get_major_ticks():
		tick.label.set_fontsize('xx-large') 
		tick.label.set_weight('bold') 
		# specify integer or one of preset strings, e.g.
		#tick.label.set_fontsize('x-small') 
		tick.label.set_rotation('vertical')
	for tick in ax.yaxis.get_major_ticks():
		tick.label.set_fontsize('xx-large') 
		tick.label.set_weight('bold') 
	# ax.set_visible(False)
	sns.set_style({'xtick.bottom': True}, {'ytick.left': True})
	for tick in heatmap.figure.axes[-1].yaxis.get_major_ticks():
		tick.label.set_fontsize(20) #does not work (last axis = legend)

	p_filename = output_folder / filename
	heatmap.get_figure().savefig(str(p_filename), bbox_inches="tight")

	strip_axis_matplotlib(heatmap.figure.axes[0]) #strip main figure
	strip_axis_matplotlib(heatmap.figure.axes[-1]) #strip legend
	s_filename = str(p_filename.parent / p_filename.stem) + "(stripped).pdf"
	heatmap.get_figure().savefig(s_filename, bbox_inches="tight")

	return

#parms
input_file_excel,output_folder = read_parms()
output_folder.mkdir(parents=True,exist_ok=True)

#get in right format
df_cell_ca  = pd.read_excel(input_file_excel, sheet_name ='df_ca_agg')
df_cell_ca.rename(columns={'nb_stack_analysis': 't','name1':'cellID_1','name2':'cellID_2'}, inplace=True) # needed in case of newest R version input
# df_cell_ca['contactarea_mean'] = df_cell_ca['contactarea_mean'] * 1.e12 # convert m² to um²

l_timepoints = [1,5]
for select_by_time in [True,False]:
	l_df=[]
	if select_by_time:
		for t_i in l_timepoints:
			df_matrix = compose_ca_matrix(df_cell_ca,select_by='time',t=t_i)
			# viz
			create_heatmap(df_matrix,
										title='Contact area between cells (timepoint {0})'.format(t_i),
										filename="1-1_heatmap contactarea timepoint{0}.png".format(t_i))
			l_df.append(df_matrix)
	else:
			for select_i in ['first','last']:
				df_matrix = compose_ca_matrix(df_cell_ca,select_by=select_i)
				# viz
				create_heatmap(df_matrix,
											title='Contact area between cells ({0} appearance of cell-cell contact)'.format(select_i),
											filename="2-1_heatmap contactarea ({0} appearance cell-cell contact).png".format(select_i))
				l_df.append(df_matrix)

	#create delta viz
	df_delta = l_df[1].fillna(0) - l_df[0].fillna(0)
	if select_by_time:
		s_title = 'Contact area change between time {0} and {1}'.format(l_timepoints[0],l_timepoints[1])
		s_filename = "1-2_heatmap contactarea delta(between time {0} and {1}).png".format(l_timepoints[0],l_timepoints[1])
	else:
		s_title = 'Contact area change between first and last appearance of cell-cell contacts'
		s_filename = "2-2_heatmap contactarea delta (first to last appearance cell).png"

	create_heatmap(df_delta,
					title=s_title,
					filename=s_filename,
					cmap='seismic',
					center=0)


# INFO Parameters 
#Makes each cell square-shaped.
# square = True,
# #Set width of the lines that will divide each cell to .5
# linewidths = .5,
# #Map data values to the coolwarm color space
# cmap = 'coolwarm',
# #Shrink the legend size and label tick marks at [-1, -.5, 0, 0.5, 1]
# cbar_kws = {'shrink': .4, ‘ticks’ : [-1, -.5, 0, 0.5, 1]},
# #Set min value for color bar
# vmin = -1, 
# #Set max value for color bar
# vmax = 1,
# #Turn on annotations for the correlation values
# annot = True,
# #Set annotations to size 12
# annot_kws = {“size”: 12})
# #Add column names to the x labels 
# ax.set_xticklabels(corr_matrix.columns)
# #Add column names to the y labels and rotate text to 0 degrees
# ax.set_yticklabels(corr_matrix.columns, rotation = 0)
# #Show tickmarks on bottom and left of heatmap
# sns.set_style({'xtick.bottom': True}, {'ytick.left': True})




