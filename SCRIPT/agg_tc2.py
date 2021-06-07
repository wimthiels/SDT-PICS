"""
Created on 22 January 2021
aggregation of testcases
generic way to aggregate data over several testcases
requirement : the testcase label must be in the file name path !
the program works as follows : 
- via regex a the files are selected (1 regex per file starting from a root folder)
- the testcase label is extracted from the file path
- if the data is an excel file this data will result in a tab in an excel file preceded by the testcase label
- all the data from the testcase itself (metadata), will be extracted from the param xml file for that testcase. this data will be written to the tc_xml tab
@author: wimth
"""

import pandas as pd
import sys, os, re
from param_XML import Param_xml
from pathlib import Path
from shutil import copy
import logging
verbose=True


class Testcase:
	d_label_instance = {}
	df_error =  pd.DataFrame(columns=['file','tc_label','regex_file_ix','regex','message'])
	
	def __init__(self, label):
		#print(f'creating cell label {label}')
		self.label = label
		self.xml_data = {}
		self.files = {}  #ey is the index of the regex list, value is path file name
		Testcase.d_label_instance[label] = self
		
	def get_param_xml_file(self):
		"first xml found will be considered"
		for path_i in self.files.values():
			if path_i.suffix=='.xml':
				return path_i

		print(f'No param xml file found for tc {self.label}')
		return None
	
	@staticmethod
	def add_error(file='',tc_label="",regex_file_ix="",regex="",message="no message specified",verbose=False):
		new_row = pd.DataFrame( [[file,tc_label,regex_file_ix,regex,message]] ,columns=Testcase.df_error.columns) #you need the right column names !
		Testcase.df_error=Testcase.df_error.append(new_row,ignore_index=True)
		if verbose:print(message)
		return

		
		
	def repr(self,show_xml=True,logging=None):
		if show_xml:
			message = f" {self.label} with {len(self.files.items())} files ({self.files.keys()})  and xml_data={self.xml_data}"
		else:
			message = f" {self.label} with {len(self.files.items())} files ({self.files.keys()})"
		print(message)
		if logging:
			logging.info(message)
	
		
	@staticmethod
	def repr_all(logging=None):
		print(f"{len(Testcase.d_label_instance.keys())} instances in Testcase class")
		for label, instance_i in sorted(Testcase.d_label_instance.items()):
			instance_i.repr(show_xml=False,logging=logging)
			
# METHODS-----------------------------------------------

def read_parms():
	param_xml = Param_xml.get_param_xml(sys.argv, l_main_keys=['body'],verbose=True) 
	param_xml.add_prms_to_dict(l_keys=['agg_tc'])
	param_xml.repr_prm()
	return param_xml.prm

def create_testcase(path,tcpattern): 
	match=tcpattern.search(str(path))
	if match: # match will be None if nothing is found
		tc_label = match.group(0)
		tc = Testcase.d_label_instance.get(tc_label)
		if not tc:
			tc = Testcase(tc_label)
	else:
		Testcase.add_error(file=path,regex=tcpattern,message=f"testcase label not found for path  {path}.  This file is skipped therefore ",verbose=True)
		return None
			
	return tc

def fill_datastructure():
	tcpattern = re.compile(prm['regex_tc'])
	for ix_regex,regex_i in enumerate(prm['l_regex_files']):
		print(f"Building up ds -> Handling regex = {regex_i}")
		for ix_path,path in enumerate(prm['root_folder'].rglob(regex_i)):
			tc = create_testcase(path,tcpattern)
			if not tc:continue
			tc.files[ix_regex] = path

	return

def enrich_xml_testcase_metadata():
	for tc_i in Testcase.d_label_instance.values():
		#print(tc_i.get_param_xml_file())
		tc_param_xml =  Param_xml.get_param_xml(['_',tc_i.get_param_xml_file()],l_main_keys = ['body'],verbose=False)
		for xml_key in prm['l_xml_keys_tc_data']:
			l_keys = xml_key.split('/')
			l_keys = [i for i in l_keys if i!='.']  #remove .
			search_field = l_keys[-1]
			key_fields = [] if len(l_keys)<2 else l_keys[:-1]
			try:
				tc_i.xml_data[search_field] = tc_param_xml.get_value(search_field,l_path_keys=key_fields)
			except:
				tc_i.xml_data[search_field]  = 'SEARCH FIELD NOT FOUND !'
	return

def construct_tc_xml(file,tc):
	"""the param xml of the testcase will be used to extract metadata """
	#print(tc.get_param_xml_file())
	d_df = {"tc_label":tc.label}
	for xml_field,value in tc.xml_data.items():
		d_df[xml_field] =value
		
	return pd.DataFrame(d_df, index=[0]) #if using all scalar values you must pass an index
	# return pd.DataFrame(d_df)

def construct_tc_xml_deprecated(file,tc):
	"""the param xml of the testcase will be used to extract metadata """
	#print(tc.get_param_xml_file())
	d_df = {"tc_label":tc.label}
	tc_param_xml =  Param_xml.get_param_xml([_,file],l_main_keys = ['body'],verbose=False)

	for xml_key in prm['l_xml_keys_tc_data']:
		l_keys = xml_key.split('/')
		l_keys = [i for i in l_keys if i!='.']  #remove .
		search_field = l_keys[-1]
		key_fields = [] if len(l_keys)<2 else l_keys[:-1]
		
		d_df[search_field] =tc_param_xml.get_value(search_field,l_path_keys=key_fields)
		
	return pd.DataFrame(d_df)

def construct_csv_tab(file,tc):
	"""simple copy with the tclabel as an extra column at the start"""
	df_file_in = pd.read_csv(file,header=0,index_col=False,keep_default_na=True) 
	df_file_in.insert(0, "tc_label", tc.label)
	return df_file_in

def construct_xlsx_tab(file,tc):
	"""simple copy with the tclabel as an extra column at the start"""
	df_file_in = pd.read_excel(file,header=0,index_col=False,keep_default_na=True) 
	df_file_in.insert(0, "tc_label", tc.label)
	return df_file_in
	
def copy_vtp_file(file,tc):
	file_extension = ""
	for ix, parm_i in enumerate(prm['l_ext_parms']):
		extension = tc.xml_data[parm_i]
		extension = str(extension) if not isinstance(extension,list) else str(extension[0])
		file_extension += f"_{prm['l_ext_letters'][ix]}{extension.zfill(2)}"
	#new_file_name = f"{file.stem}_R{str(tc.xml_data['repID']).zfill(2)}_t{str(tc.xml_data['nb_stack_analysis'][0]).zfill(2)}"
	#copy(file, (prm['output_folder'] / file.stem / (new_file_name + file.suffix)))
	print('vtp written to ', (prm['output_folder'] / file.stem / (file.stem + file_extension + file.suffix)),flush=True)
	copy(file, (prm['output_folder'] / file.stem / (file.stem + file_extension + file.suffix)))
	
	return


# MAIN-----------------------------------------------
prm = read_parms()
prm['output_folder'].mkdir(parents=True,exist_ok=True)
logging.basicConfig(filename=(prm['output_folder'] / "log.txt"), filemode='w', format='%(name)s - %(levelname)s - %(message)s',level=logging.INFO)

fill_datastructure()
enrich_xml_testcase_metadata()
if verbose:
	print("Testinfra Class filled")
	Testcase.repr_all(logging)

file_suffix = ""
d_tab_dfagg = {}
for ix_regex,regex_i in enumerate(prm['l_regex_files']):
	l_df = []
	tab_name = ""
	for ix_file,(tclabel_i, tc_i) in enumerate(sorted(Testcase.d_label_instance.items())):
		#print(f"handling file-{ix_regex} for tc={tclabel_i}")
		file = tc_i.files.get(ix_regex)
		if not file:
			Testcase.add_error(tc_label=tclabel_i,regex_file_ix=ix_regex,regex=regex_i,message=f'tc={tclabel_i} does not have a file with index-{ix_regex}, moving on with next testcase..',verbose=True)
			continue
		
		if not tab_name:
			file_suffix = file.suffix
			tab_name = file.stem if file_suffix!='.xml' else "tc_xml"
			print(f"These {ix_regex}-files are treated as {file_suffix}-files and will get tabname= {tab_name}")
			
		try:
			if file_suffix=='.xml':
				l_df.append(construct_tc_xml(file,tc_i))
			elif file_suffix=='.csv':
				l_df.append(construct_csv_tab(file,tc_i))
			elif file_suffix=='.xlsx':  #only first tab is extracted
				l_df.append(construct_xlsx_tab(file,tc_i))
			elif file_suffix in ['.vtp','.tif']:
				if not (prm['output_folder'] / file.stem).exists():
					(prm['output_folder'] / file.stem).mkdir(parents=True)
				copy_vtp_file(file,tc_i)
			else:
				Testcase.add_error(tc_label=tclabel_i,regex_file_ix=ix_regex,regex=regex_i,message=f'appending file {file} for tc={tclabel_i} was not possible because the file extension is not supported (yet)',verbose=True)
		except Exception as e:
			Testcase.add_error(tc_label=tclabel_i,regex_file_ix=ix_regex,regex=regex_i,message=f'appending file {file} for tc={tclabel_i} resulted in an error={e} Check the file...',verbose=True)
	
	if l_df:
		d_tab_dfagg[tab_name] = pd.concat(l_df,ignore_index=True) #store away
	
	
with pd.ExcelWriter(prm['output_folder']/"agg_tc.xlsx",mode='w') as writer:
	
	for tab_name_i, df_agg in d_tab_dfagg.items():
		df_agg.sort_values(by='tc_label',ascending=True,inplace=True) 
		df_agg.to_excel(writer, sheet_name=tab_name_i,index=False)
		
	Testcase.df_error.sort_values(by='tc_label',ascending=True,inplace=True) 
	Testcase.df_error.to_excel(writer, sheet_name="errors",index=False)


print("Aggregated data can be found in \n {}".format(prm['output_folder']/"agg_tc.xlsx"))


#append extra excel
if prm['append_file1']:
	excel1 = prm['output_folder']/"agg_tc.xlsx"
	excel2 = prm['append_file1']

	d_sheet_df = {}
	for sheet_name_i in pd.ExcelFile(excel1).sheet_names: 
		print(sheet_name_i)
		try:
			d_sheet_df[sheet_name_i] = pd.concat([pd.read_excel(excel1,sheet_name=sheet_name_i),pd.read_excel(excel2,sheet_name=sheet_name_i)])
		except Exception as e:
			print(f'{sheet_name_i} caused error {e}. Copying data from excel1 only')
			d_sheet_df[sheet_name_i] = pd.read_excel(excel1,sheet_name=sheet_name_i)
			pass

	with pd.ExcelWriter(prm['output_folder']/"agg_tc_appended.xlsx",mode='w') as writer:
		for sheet_name_i, df_i in d_sheet_df.items():
			df_i.to_excel(writer, sheet_name=sheet_name_i,index=False)
	print("Appended data can be found in \n {}".format(prm['output_folder']/"agg_tc_appended.xlsx"))

