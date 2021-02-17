'''
Created on 27 Mar 2019

@author: wimth
'''
import xmltodict
import os,re
from pathlib import Path

class Param_xml():
	
	def __init__(self,file,l_main_keys):
		
		self.file = Path(file)
		self.l_main_keys = l_main_keys
		
		with open(file) as fd:
			self.d_prms= xmltodict.parse(fd.read()) #full dict unresolved parms (raw)
		self.prm = {} #dict for resolved parms 


			
	def get_value(self,s_parm,l_path_keys=[], use_main_keys=True, lookup=None,dtype=None,all=None, docker_run=False):
		
		def resolve_path_name(value):
			if "<ROOT_OUTPUT>" in value:
				OUTPUT_FOLDER_ROOT = self.get_value("OUTPUT_FOLDER_ROOT",l_path_keys=['body','MAIN','paths'], use_main_keys=False,dtype='string')
				if OUTPUT_FOLDER_ROOT != "None":
					value = value.replace('<ROOT_OUTPUT>',str(OUTPUT_FOLDER_ROOT))
			
			if "<SUB_FOLDER>" in value:
				OUTPUT_FOLDER_ROOT_SUB = self.get_value("OUTPUT_FOLDER_ROOT_SUB",l_path_keys=['body','MAIN','paths'], use_main_keys=False,dtype='string')
				if OUTPUT_FOLDER_ROOT_SUB!="None":
					value = value.replace('<SUB_FOLDER>',OUTPUT_FOLDER_ROOT_SUB)

			if "<runID>" in value:
				runID = self.get_value("runID",l_path_keys=['body','LOG'], use_main_keys=False,dtype='string')
				if runID!="None":
					value = value.replace('<runID>',runID)

			if value.startswith('..'):
				value = value.replace('..',str(Path(os.path.realpath(__file__) ).parents[1])) #location of this script

			if value.startswith('.'):
				value = value.replace('..',str(Path(os.path.realpath(__file__) ).parents[0]))

			return os.path.normpath(value)

		def fill_in_root_folder(value):
			OUTPUT_FOLDER_ROOT = self.get_value("OUTPUT_FOLDER_ROOT",l_path_keys=['body','MAIN','paths'], use_main_keys=False,dtype='string')
			value = value.replace('<ROOT_OUTPUT>/','')

			if "<SUB_FOLDER>" in value:
				OUTPUT_FOLDER_ROOT_SUB = self.get_value("OUTPUT_FOLDER_ROOT_SUB",l_path_keys=['body','MAIN','paths'], use_main_keys=False,dtype='string')
				value = value.replace('<SUB_FOLDER>/','')
				if OUTPUT_FOLDER_ROOT_SUB!="None":
					return os.path.join(OUTPUT_FOLDER_ROOT,OUTPUT_FOLDER_ROOT_SUB,value)
						
			return os.path.join(OUTPUT_FOLDER_ROOT,value)
		

		def compose_file_path():
			file_nm = d_param_sub[s_parm]['@value']
			s_folder = str(l_path_keys[-1])
			folder_path = self.get_value(s_folder,l_path_keys=l_path_keys[:-1], use_main_keys=False)

			if not isinstance(folder_path, Path):
				folder_path = Path(folder_path)
			  
			return folder_path / file_nm

			
		s_parm = s_parm.strip()
		if isinstance(l_path_keys, str):
			l_path_keys = [l_path_keys]
		
		if use_main_keys:
			l_path_keys = self.l_main_keys + l_path_keys
				
		
		d_param_sub = self.retrieve_element(l_path_keys)
		
		#@deprecated : is replaced by +
		# ic_link=True                                                                 # apply links in iteration
		# while ic_link:
			# d_param_sub = retrieve_element(l_path_keys)
			# ic_link=False
			# if not d_param_sub:return
			# value = d_param_sub[s_parm].get('@value')
			# if value=='link' or "+" in value:
				# l_path_keys = d_param_sub[s_parm].get('link').split('+')
				# if l_path_keys and len(l_path_keys)>1:   
					# s_parm = l_path_keys[-1]
					# l_path_keys = l_path_keys[0:-1]    
					# ic_link=True
		
		
		
		
		if dtype=='dict':
			return d_param_sub[s_parm]
				
		elif s_parm.startswith('file'): # 'file' triggers recursion to retrieve folder name ! dtype path is implicit
			file_path = str(compose_file_path())

			if docker_run:
				file_path = file_path.replace(self.d_prms["body"]["RAM"]["script_folder"]["@value"],self.d_prms["body"]["RAM"]["docker_script_folder"]["default"])
				file_path = file_path.replace(self.d_prms["body"]["RAM"]["data_folder"]["@value"],self.d_prms["body"]["RAM"]["docker_data_folder"]["default"])

			return Path(file_path)
			
		else:
			value = d_param_sub[s_parm]['@value']
			if value == 'link': print('temp : stop using link , convert to +')
			if value == "default":
				value = d_param_sub[s_parm]['default']
		

		if value and "+" in value:  # + triggers recursion !
			l_path_keys=value.split('+')
			return self.get_value(l_path_keys[-1],l_path_keys=l_path_keys[:-1], use_main_keys=False, lookup=lookup, docker_run=docker_run)
			
		if lookup:                                                                  # look up a value
			value = d_param_sub[s_parm].get(lookup)
			if not value:
				value = d_param_sub[s_parm].get('default')
			if not value:
				value = d_param_sub[s_parm].get('@value')  #just as backup         
		
		dtype = d_param_sub[s_parm].get('dtype',dtype)                                #apply dtype
		if isinstance(dtype,dict):
			dtype=dtype['@value'] #in the rare cases we have updated dtype, we retrieve the update here (the old entry can be retrieved with '#text')
		
		ic_list=True if (value and ';' in value) else False
		if not dtype or dtype=='string':
			value = value.split(';') if ic_list else str(value)  #if dtype is not mentioned in param file, string is assumed
		if dtype=='float':
			value = list(map(float,value.split(';'))) if ic_list else float(value) 
		if dtype=='int':
			value = list(map(int,value.split(';'))) if ic_list else int(value)
		if dtype=='list;int':
			if value=='all':
				pass
			else:
				try:
					value = list(map(int,value.split(';'))) if ic_list else [int(value)]
				except:
					value = []
		if dtype=='list;float':
			if value=='all':
				pass
			else:
				try:
					value = list(map(float,value.split(';'))) if ic_list else [float(value)]
				except:
					value = []
		if dtype=='list;string':
			value = list(map(str,value.split(';'))) if ic_list else [str(value)]
		if dtype=='bool':
			if "false" in value.lower() or value=='0':
				value= False
			else:
				value = True 
		if dtype=='tuple;int':
			value = tuple(map(int,value.split(';')))
		if dtype=='tuple':
			value = tuple(value.split(';')) 
				
		if dtype=='path':
			if value:
				value = resolve_path_name(value)
				# if value.startswith('<ROOT_OUTPUT>'):
				#     value = fill_in_root_folder(value)
				
				if docker_run:
					value = value.replace(self.d_prms["body"]["RAM"]["script_folder"]["@value"],self.d_prms["body"]["RAM"]["docker_script_folder"]["default"])
					value = value.replace(self.d_prms["body"]["RAM"]["data_folder"]["@value"],self.d_prms["body"]["RAM"]["docker_data_folder"]["default"])

				value = Path(value)

		if value == 'all':
			if not all:
				pass
			else:
				value = range(1,all)
			
		
		return value  #default= string
		
	def __repr__(self):
		print('file=',str(self.file))
		print('l_main_keys=',self.l_main_keys)
		print('self.d_prms=',self.d_prms)
		return
	
	
	def store_value(self,s_parm,value_to_store,l_path_keys=["body","RAM"],verbose=True):
		d_parm = self.d_prms
		for i in l_path_keys:
			d_parm = d_parm[i]
		if not d_parm.get(s_parm):print('error : the key {0} is not present in the PARAM.xml file'.format(l_path_keys.append(s_parm)));return
		d_parm[s_parm] = str(value_to_store)
		
		if verbose:print ('value {0} stored at {1}'.format(format(str(value_to_store),l_path_keys.append(s_parm))))
		return
		
	@staticmethod
	def get_param_xml(sys_argv, l_main_keys=[],add_to_dict=False,verbose=False):
		'''
		in order of priority
		1) if a Param_xml object is passed as the FIRST argument use this
		2) if a path is passed as the FIRST argument use this to create a new Param_xml object
		3) if no arguments are given, load the PARAM.xml file in the current working directory
		'''
		if len(sys_argv)>1:
			if isinstance(sys_argv[1],Param_xml): 
				param_xml = sys_argv[1]
			else:    
				param_xml = Param_xml(file=sys_argv[1],l_main_keys=l_main_keys)
		else:
			file_name='PARAMS.xml'
			param_xml = Param_xml(file=os.path.join(os.getcwd(),file_name),l_main_keys=l_main_keys)


		if add_to_dict:
			param_xml.add_prms_to_dict()
		
		if verbose:print('param_xml_object from location {0} is used'.format(param_xml.file))      
		return param_xml       

	def retrieve_element(self, l_path_keys):
		if len(l_path_keys)==0:
			d_param_sub = self.d_prms
		elif len(l_path_keys)==1:
			d_param_sub = self.d_prms[l_path_keys[0]]
		elif len(l_path_keys)==2:
			d_param_sub = self.d_prms[l_path_keys[0]][l_path_keys[1]]
		elif len(l_path_keys)==3:
			d_param_sub = self.d_prms[l_path_keys[0]][l_path_keys[1]][l_path_keys[2]]
		elif len(l_path_keys)==4:
			d_param_sub = self.d_prms[l_path_keys[0]][l_path_keys[1]][l_path_keys[2]][l_path_keys[3]]
		elif len(l_path_keys)==5:
			d_param_sub = self.d_prms[l_path_keys[0]][l_path_keys[1]][l_path_keys[2]][l_path_keys[3]][l_path_keys[4]]
		else:
			print('ERROR : exceeded max depth for a xml path',l_path_keys)
			return None
		return d_param_sub
		
	def add_prms_to_dict(self,l_keys=[]):
		""" resolve all the parameters in a param-object and add the parameters to a dict
		a parameter is identifed if it has a key= '@value'
		"""
		def add_prms_to_dict_recursion(d_iter,prm={},key='',keys=[]):
			if  '@value' in d_iter:
				#print('key={},keys={}'.format(key,keys))
				prm[key] = self.get_value(key, keys)
			else:
				for ix,key_i in enumerate(d_iter.keys()):
					if key and ix==0:
						keys.append(key)
					add_prms_to_dict_recursion(d_iter[key_i],prm,key_i, keys)
				if keys:
					keys.pop()
			return prm

		d_sel = self.retrieve_element(self.l_main_keys + l_keys)

		d_prm = add_prms_to_dict_recursion(d_sel,keys=l_keys)
		self.prm.update(d_prm)
	
		return self.prm

	def repr_prm(self):
		import pprint
		pp = pprint.PrettyPrinter(indent=4)
		pp.pprint(self.prm)
		return