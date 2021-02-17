'''
Created on  02 July 2020 
Vizualisation of Results : tracking and segmentation from agg_seg_data

@author: wimth
'''
import seaborn as sns
import pandas as pd 
from pathlib import Path
from param_XML import Param_xml
import matplotlib
import matplotlib.pyplot as plt
import os,sys,re



def read_parms():
	param_xml = Param_xml.get_param_xml(sys.argv,l_main_keys = ['body','viz_results2'],verbose=True) #param file must be passed as first argument

	d_prm = {}
	d_prm['f_agg_seg_data'] = param_xml.get_value('f_agg_seg_data',['paths'])
	d_prm['output_folder'] = param_xml.get_value('output_folder',['paths'])
	
	return d_prm

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

	# ax.xaxis.set_visible(False)  # removes axis labels, tick labels and ticks
	# ax.yaxis.set_visible(False)

	return ax

def strip_seaborn_plot(catplot):
	catplot._legend.remove()
	catplot.set_xticklabels([])
	catplot.set_yticklabels([])

	return catplot


def add_sortkey_and_tcinfo(df_agg_scores,df_tc_info):

	df_agg_scores = df_agg_scores.merge(d_df['tc_info'],how='left',on='tcID')

	l_sortkey = []
	for i in df_agg_scores['cell_stage']:
		sortkey = re.search("\d*",i).group(0)
		sortkey = int(sortkey) if sortkey else 999
		l_sortkey.append(sortkey)
	df_agg_scores['sortkey'] = l_sortkey
	df_agg_scores = df_agg_scores.sort_values(by='sortkey')

	df_agg_scores.to_csv(d_prm['output_folder'] / 'agg_scores_w_tc_ID.csv',index=False)

	return df_agg_scores

def plot_results_tracking(df_tra):

	cond = df_tra['type'].isin(['TRA','sensitivity','precision'])
	df_tra = df_tra[cond]

	# FIGURE1A : point plot for results tracking (all 3 metrics)
	#fig, ax = plt.subplots(figsize=(10,6))
	catplot = sns.catplot(x="cell_stage", y="score_x",hue='type',kind="point", data=df_tra,legend=False) #facetplot
	catplot.add_legend(title='', loc='center',fontsize='xx-small')
	catplot.set(ylim=(0, 1.0))
	catplot.set_xlabels('')
	catplot.set_ylabels('')
	l_cell_stage = [re.findall(r'\d+', s)[0] for s in df_tra['cell_stage'].unique()]
	catplot.set_xticklabels(l_cell_stage ,size = 'x-small',weight='normal')
	catplot.set_yticklabels([0.0, 0.2,0.4,0.6,0.8,1] ,size = 'x-small',weight='normal')
	#grid.set(xticks=np.arange(5), yticks=[-3, 3],xlim=(-.5, 4.5), ylim=(-3.5, 3.5)) example set
	catplot.savefig(str(d_prm['output_folder'] / 'results_tracking.png'))

	catplot = strip_seaborn_plot(catplot)
	catplot.savefig(str(d_prm['output_folder'] / 'results_tracking(stripped).pdf'))

	# FIGURE1B: same results in a barplot
	fig, ax_barplot = plt.subplots(figsize=(3.2,2.4),dpi=350) #1inch = 2.54 cm (3.2,2.4) = 1/2 pagewidth; dpi=350=publication ready
	ax_barplot = sns.barplot(x="cell_stage", y="score_x", hue='type',hue_order=['TRA','sensitivity','precision'],palette='summer',data=df_tra)
	ax_barplot.legend(title='', loc='lower left',fontsize='xx-small',framealpha=1)
	ax_barplot.set_ylabel('') 
	ax_barplot.set_xlabel('')
	ax_barplot.set(ylim=(0, 1.0))
	ax_barplot.set_xticklabels([2,7,26,48] ,size = 'x-small',weight = 'normal')
	ax_barplot.set_yticklabels([0.0, 0.2,0.4,0.6,0.8,1] ,size = 'x-small',weight = 'normal')
	
	fig.savefig(str(d_prm['output_folder'] / 'results_tracking(barplot).pdf'))
	ax_barplot = strip_axis_matplotlib(ax_barplot)
	fig.savefig(str(d_prm['output_folder'] / 'results_tracking(barplot)(stripped).pdf'))


	#FIGURE2: same data both shown as combo chart (bar + line)
	#problem : i want the same y-ticks on both sides
	fig, ax1 = plt.subplots(figsize=(3.2,2.4),dpi=350)
	color = 'tab:green'
	#bar plot creation
	ax1.set_title('results tracking', fontsize=16)
	ax1.set_xlabel('embryo', fontsize=16)
	ax1.set_ylabel('TRA', fontsize=16,color=color)
	ax1 = sns.barplot(x='cell_stage', y='score_x', data = df_tra[df_tra['type']=='TRA'], palette='summer')
	ax1.tick_params(axis='y')
	#specify we want to share the same x-axis
	ax2 = ax1.twinx()
	
	color = 'tab:red'
	ax2.set_ylabel('%', fontsize=16,color=color)
	# df_tidy2 = pd.melt(df_tra, id_vars=['ID','nb_cells_t0'], value_vars=['sensitivity','precision'],
	# 		 var_name='metric', value_name='score')
	# df_tidy2 = df_tidy2.sort_values(by='nb_cells_t0')
	# ax2 = sns.catplot(x="ID", y="score",hue="metric",kind="point", data=df_tra_tidy2); #comboplot does not work for catplot (overwrites other axis)
	# ax2 = sns.stripplot(x="ID", y="score",hue="metric",data=df_tidy2,jitter=False,size=20)
	ax2 = sns.pointplot(x="cell_stage", y="score_x",hue="type",data=df_tra[df_tra['type'].isin(['sensitivity','precision'])],color=color)
	ax2.tick_params(axis='y') 
	# ax2.set_yticklabels(color=color)
	#show plot
	# plt.show()
	plt.savefig(d_prm['output_folder'] / 'results_tracking_combo.png')


	return


def plot_results_tracking_detail(df_agg_scores,df_tracking_detail):
	#enrich with nb_cell_t0 needed as a sort key
	df_tra = df_tracking_detail.merge(df_agg_scores[df_agg_scores['type']=='TRA'],how='left',on='tcID')
	df_tra.sort_values(by=['sortkey','timestep'],inplace=True)
	df_tra_long = pd.melt(df_tra, id_vars=['tcID','timestep'], value_vars=['TP','FN','FP'], var_name='type', value_name='nb_cells')
	df_tra_long['time'] = df_tra_long['timestep'] * 3  #delta_t = 3 minutes

	#all together on a strip 
	colors = ["grass green", "red", "orange"]
	sns.set(style="ticks", color_codes=True)
	sns_grid = sns.FacetGrid(df_tra_long, col="tcID", col_wrap=2,hue="type",palette=sns.xkcd_palette(colors),height=2.4,aspect=1.3) #aspect * height gives the width of each facet in inches.
	# sns_grid.map(plt.scatter, "t", "nb_cells", edgecolor="w")
	#sns_grid.map(plt.plot, "timestep", "nb_cells",linewidth=3,alpha=0.7,xticks=[0,5,10,15,20], yticks=[0, 25, 50,75,100])
	sns_grid.map(plt.plot, "time", "nb_cells",linewidth=3,alpha=0.7, marker="o")
	sns_grid.set(xticks=[0,15,30,45,60],yticks=[0,25,50,75,100])
	axes = sns_grid.axes.flatten()
	l_cell_stage = [re.findall(r'\d+', s)[0] for s in df_tra['cell_stage'].unique()]
	sns_grid.add_legend(title="")
	for ix,cell_stage in enumerate(l_cell_stage):
		#axes[ix].set_title(cell_stage)
		axes[ix].set_title("")  # blank and use inkscape
		for tick in axes[ix].xaxis.get_major_ticks():
			tick.label.set_fontsize(14)
		for tick in axes[ix].yaxis.get_major_ticks():
			tick.label.set_fontsize(14)

	sns_grid.savefig(str(d_prm['output_folder'] / 'results_tracking_detail.pdf'))

	sns_grid.set_xlabels('',size=16)  # blank instead of name : use inkscape
	sns_grid.set_ylabels('',size=16)
	#sns_grid.set(xticks=[],yticks=[])
	#
	sns_grid.savefig(str(d_prm['output_folder'] / 'results_tracking_detail(stripped).pdf'))


	# sns_grid = sns_grid.map_dataframe(timeplot, "t", "nb_cells")

	#each tc separately
	for tcID_i in df_tra.tcID.unique():
		fig, ax = plt.subplots(figsize=(3.2,2.4),dpi=350)
		cond = df_tra['tcID'] == tcID_i
		df_data = df_tra[cond][['timestep','TP','FN','FP']].apply(pd.to_numeric)
		df_data.set_index(keys=['timestep'],inplace=True)
		ax = sns.lineplot(data=df_data, palette="tab10", linewidth=7)
		ax.set_xlim([0,25])
		ax.set_ylim([0,100])
		plt.ylabel('nb cells', fontsize=18)
		plt.xlabel('t', fontsize=18)
		ax.tick_params(axis='both', which='major', labelsize=18)
		ax.get_legend().remove()
		#ax.legend(loc='upper right', shadow=False, fontsize='x-large')
		plt.savefig(d_prm['output_folder'] / 'detail_TRA_{0}.png'.format(tcID_i))


		ax = strip_axis_matplotlib(ax)
		plt.savefig(d_prm['output_folder'] / 'detail_TRA_{0}(stripped).pdf'.format(tcID_i))

	return

def plot_results_tracking_metrics(df_tracking_metrics):

	def melt_column_range_and_add_t(s_column,df):
		l_col_names = [col for col in df.columns if s_column in col]
		df_long =  pd.melt(df, id_vars=['cell_stage','cluster_label'], value_vars=l_col_names)
		df_long['t'] = [int(re.findall(r'\d+', s)[0]) for s in df_long['variable']]

		df_long.rename(columns={"value": s_column},inplace=True)
		df_long.drop(["variable"],inplace=True,axis=1)
		df_long.set_index(['cell_stage', 'cluster_label','t'],inplace=True)

		return df_long


	#First get into long format-----------------------------------------
	# 	step 1 : melt 3 column ranges into 3 df
	df_tracking_metrics.reset_index(level=0, inplace=True)  #turn index cell_stage into column
	df_pixel_count   = melt_column_range_and_add_t('nb_pixels',df_tracking_metrics)
	df_seed_centre_z = melt_column_range_and_add_t('seed_centre_z',df_tracking_metrics)
	df_name_lineage = melt_column_range_and_add_t('name_lineage',df_tracking_metrics)

	# 	step 2 : join everything back together to get a long df with key columns = [cell_stage,name_lineage,z] and variables= nb_pixels and seed_centre_z
	df_join = df_name_lineage.join([df_pixel_count,df_seed_centre_z])  #joining on index = cell_stage, cluster_label
	df_join.reset_index(inplace=True)
	df_long =  pd.melt(df_join, id_vars=['cell_stage', 'name_lineage','t'], value_vars=['nb_pixels','seed_centre_z'])
	#	step 3 : clean up , filter
	cond =  (df_long['name_lineage'].isnull() | (df_long['value'].isnull()) )
	df_long.drop(df_long[cond].index,inplace=True)
	df_long = df_long.astype({"value":int,"name_lineage":str})
	cond = (df_long['name_lineage'] =='0') | (df_long['name_lineage'] =='0.0')
	df_long.drop(df_long[cond].index,inplace=True)

	df_long.sort_values(by=['cell_stage','name_lineage','t'],inplace=True)


	#	step 3 : get the reference values in order to calculate the deltas
	cond1 = (df_long['variable'] == 'seed_centre_z')
	cond2 = ~(df_long['value'].isnull() | (df_long['value'] == 0))
	df_long_sel = df_long[cond1 & cond2]
	df_long_sel = df_long_sel.sort_values(by=['cell_stage','name_lineage','t'])
	ser_seed_centre_z_first = df_long[cond1 & cond2].groupby(['cell_stage','name_lineage']).first()['value'] 
	cond3 = (df_long['variable'] == 'nb_pixels')
	ser_nb_pixels_mean = df_long[cond3 & cond2].groupby(['cell_stage','name_lineage']).mean()['value']
	ser_nb_pixels_std = df_long[cond3 & cond2].groupby(['cell_stage','name_lineage']).std()['value']

	def calc_delta(row):
		if row['variable']=='seed_centre_z':
			return row['value'] - ser_seed_centre_z_first[row['cell_stage'],row['name_lineage']]
			return 
		if row['variable']=='nb_pixels':
			return row['value'] - ser_nb_pixels_mean[row['cell_stage'],row['name_lineage']]
		return row['value']

	for index, row in df_long.iterrows():
		df_long.at[index,'value'] = calc_delta(row)

	#	step 4 : replace naming
	#df_seg["type"].replace({"SEG_REL": "JAC_REL", "SEG3D": "JAC", "dice": "Dice"}, inplace=True) # swithing to JAC naming

	with pd.ExcelWriter(d_prm['output_folder'] / 'tracking_metrics_long.xlsx') as writer:
		df_long.to_excel(writer, sheet_name='df_long',index=True)
		ser_nb_pixels_mean.to_excel(writer, sheet_name='vol_mean',index=True)
		ser_nb_pixels_std.to_excel(writer, sheet_name='vol_std',index=True)

	#Next,make the graphs----------------------------------------

	def save_lineplot(variable,y_label):
		fig, ax = plt.subplots(figsize=(3.2,2.4),dpi=350)
		sns.set(style="ticks")
		cond = (df_long['variable'] == variable)
		ax = sns.lineplot(x="t", y="value",
			units="name_lineage",estimator=None,
			hue="cell_stage",hue_order=['2_cell','7_cell','26_cell','48_cell'],
			data=df_long[cond])
		plt.setp([ax.get_children()],alpha=.05)
		ax.set_ylabel(y_label)
		plt.savefig(d_prm['output_folder'] / 'detail_TRA_metrics_lineplot_{0}.png'.format(variable))

		return

	def save_relplot():
		# 	Plot the lines on two facets
		# 		data should be in long format
		# 		col = the 2 palletes left and right (one for delta vol, one for delta z)
		# 		size = differentiates size of line (not needed in our case
		# 		dfhue = differentiate color based on cell_stage
		fig, ax = plt.subplots(figsize=(3.2,2.4),dpi=350)
		palette = dict(zip(df_long.cell_stage.unique(), sns.color_palette("rocket_r", 6)))
		ax =sns.relplot(x="t", y="value",
			hue="cell_stage", col="variable", palette=palette,
			height=2.4, aspect=1.3, facet_kws=dict(sharex=False,sharey=False),
			kind="line", legend="full", data=df_long)

		plt.savefig(d_prm['output_folder'] / 'detail_TRA_metrics_relplot.png')

		return

	save_lineplot('seed_centre_z',y_label='delta z')
	save_lineplot('nb_pixels',y_label='delta cell volume around mean')
	save_relplot()
	return

def plot_results_segmentation(df_seg):
	cond = df_seg['type'].isin(['SEG3D','dice','SEG_REL'])
	cond2 = ~df_seg['tcID_relative'].isnull()
	df_seg = df_seg[cond & cond2]

	fig, ax_barplot = plt.subplots(figsize=(3.2,2.4),dpi=350)
	df_seg["type"].replace({"SEG_REL": "JAC_REL", "SEG3D": "JAC", "dice": "Dice"}, inplace=True) # swithing to JAC naming
	#df_seg.rename(columns={'SEG3D':'JAC','SEG_REL':'JAC_REL'}, inplace=True)  
	ax_barplot = sns.barplot(x="sortkey", y="score_x", hue='type',hue_order=['JAC','Dice','JAC_REL'],palette='summer',data=df_seg,errwidth=1, capsize=0) # ,capsize=0.1
	ax_barplot.legend(title='', loc='lower left',fontsize='xx-small',framealpha=1)
	ax_barplot.set_ylabel('') 
	ax_barplot.set_xlabel('')
	ax_barplot.set(ylim=(0, 1.0))
	ax_barplot.set_xticklabels([2,7,26,48] ,size = 'x-small',weight = 'normal')
	ax_barplot.set_yticklabels([0.0, 0.2,0.4,0.6,0.8,1] ,size = 'x-small',weight = 'normal')

	
	fig.savefig(str(d_prm['output_folder'] / 'results_segmentation.pdf'))

	ax_barplot = strip_axis_matplotlib(ax_barplot)
	fig.savefig(str(d_prm['output_folder'] / 'results_segmentation(stripped).pdf'))


	return


d_prm = read_parms()

#read data
d_df = pd.read_excel(d_prm['f_agg_seg_data'],sheet_name=None,header=0,index_col=0,dtype=object)

#add sortkey
df_agg_scores = add_sortkey_and_tcinfo(d_df['agg_scores'],d_df['tc_info'])

#plotting
plot_results_tracking(df_agg_scores)
plot_results_tracking_detail(df_agg_scores,d_df['tracking detail'])
plot_results_tracking_metrics(d_df['tracking metrics'])
plot_results_segmentation(df_agg_scores)



print('Finished...viz_results2')


