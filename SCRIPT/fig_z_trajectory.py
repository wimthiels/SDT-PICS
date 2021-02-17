'''
Created on 22 Jul 2018

@author: wimth
'''

import pandas as pd
import sys
import pprint as pp
import copy
from param_XML import Param_xml


def fig_z_trajectory(param_xml):
    def read_parms():
        param_xml = Param_xml.get_param_xml(sys.argv,l_main_keys = ['body','fig_z_trajectory'],verbose=True) #param file must be passed as first argument
        input_file_csv= param_xml.get_value('input_file_csv',['paths'])
        output_file_csv= param_xml.get_value('output_file_csv',['paths'])
        
        return input_file_csv,output_file_csv 
        
    def add_cell_division_marker(ix, z):
        nonlocal nb_division
        nb_division += 1
        label = 'birth_{0}'.format(nb_division)
        l_z = [EMPTY] * (max_z_ix +1)
        l_z[ix] = z
        d_result[label] = l_z
        print('cell division was created (label lab_i, ix, z)', label , lab_i, ix, z)

        return 
        
    def add_cell_death_marker(ix, z):
        nonlocal nb_death
        nb_death += 1
        label = 'death_{0}'.format(nb_death)
        l_z = [EMPTY] * (max_z_ix +1)
        l_z[ix] = z
        d_result[label] = l_z
        print('cell death was created (label lab_i, ix, z)',label,  lab_i, ix, z)

        return
        
    input_file_csv,output_file_csv  = read_parms()

    df_ipt = pd.read_csv(input_file_csv)
    d_lab = {}  #  key = label , value = list with z coordinates
    d_lab_name = {}
    d_lab_parent = {}
    columns = list(df_ipt)
    for _,row_i in df_ipt.iterrows():
        l_z = []
        for column_i in columns:
            
            if column_i.endswith('_seed_centre_z'):
                l_z.append(int(row_i[column_i]))
            d_lab[str(row_i.cluster_label).zfill(3)]=l_z
            d_lab_name[str(row_i.cluster_label).zfill(3)] = row_i.color
            d_lab_parent[str(row_i.cluster_label).zfill(3)] = str(row_i.label_parent).zfill(3)

    # pp.pprint(d_lab)
    # pp.pprint(d_lab_name)
    # pp.pprint(d_lab_parent)

    d_result = copy.deepcopy(d_lab)
    EMPTY= '#N/A'
    nb_division = 0
    nb_death = 0

    #update z-parentlink
    for lab_i, l_z_i in d_lab.items():
        max_z_ix = len(l_z_i) - 1
        ic_birth = False
        #if d_lab_parent[lab_i]=='99':continue
        for ix,z_ii in enumerate(l_z_i):
            if z_ii != 0: 
                ic_birth= True
                continue
            if ix== max_z_ix:
                d_result[lab_i][ix] = EMPTY
                continue
            if l_z_i[ix+1]==0:
                d_result[lab_i][ix] = EMPTY
                if ic_birth:
                    add_cell_death_marker(ix-1,l_z_i[ix-1],nb_death)
                    ic_birth=False
                continue
            if d_lab_parent[lab_i]=='99': print('error, row has zeros but no parent', lab_i, l_z_i)
            d_result[lab_i][ix] = d_lab[d_lab_parent[lab_i]][ix]
            add_cell_division_marker(ix,d_lab[d_lab_parent[lab_i]][ix])
            

            
    #print('d_lab');pp.pprint(d_lab)      
    print('d_RESULT')
    for i,j in d_result.items():print(i,'=>',j)   
       
    df = pd.DataFrame.from_dict(d_result,orient='index').sort_index()  #index orientatie is meest logische voor een simpele dict
    df.to_csv(output_file_csv,index=True)
        
    return