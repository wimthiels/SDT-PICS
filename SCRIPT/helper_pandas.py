# -*- coding: utf-8 -*-
"""
Spyder Editor
******************************************************************************
creating panda dataframes : cdf chapter 8.2
******************************************************************************
DataFrame is a 2-dimensional labeled data structure with columns of potentially different types. You can think of it
like a spreadsheet or SQL table, or a dict of Series objects. It is generally the most commonly used pandas object.
Like Series, DataFrame accepts many different kinds of input:
• Dict of 1D ndarrays, lists, dicts, or Series
• 2-D numpy.ndarray
• Structured or record ndarray
• A Series
• Another DataFrame
Along with the data, you can optionally pass index (row labels) and columns (column labels) arguments. If you pass
an index and / or columns, you are guaranteeing the index and / or columns of the resulting DataFrame. Thus, a dict
of Series plus a specific index will discard all data not matching up to the passed index.
If axis labels are not passed, they will be constructed from the input data based on common sense rules
"""
#help("itertools")
import itertools
import pandas as pd
import numpy as np

def excel_to_dds(excel_name,sheet_name=None,d_index_cols={},verbose=True):
    """ 
    turn an excel into a dict with datastructures
    access as follows : 
    dds['tabname']['index']['columnname value']
    by default the index will be constructed from the first column of the excel. in case of a multiindex specify this via d_index_cols 
    eg d_index_cols{'tabname2':[0]} => this will construct a multiindex by appending column 1 (the first column will always be part of the index)  'index' must then be specified with a tuple ! eg (2,'cell_label_0')
    sheet_name= None will select all sheets, but you can also state sheets =  [0,1,2] to select only the first three tabs

    """

    d_df = pd.read_excel(excel_name,sheet_name=sheet_name,header=0,index_col=0,dtype=object)  

    if d_index_cols: #mutli indexing is needed for a certain tab
        for tabname_i, index_cols_i in d_index_cols.items():
            df = d_df[tabname_i]
            df.set_index(keys=[df.columns[i] for i in index_cols_i],append=True,inplace=True)  

    #  now turn the dict of df into a dict of dicts
    dds = {}  # dictionary of datastructures
    for tabname_i, df_i in d_df.items():
        dds[tabname_i] = df_i.to_dict(orient="index")

    return dds


    #sheet_name None will process all tabes, 2 refers to third sheet, string as name is also possible
    #header=0 takes first row as header
    #dtype=object : Use object to preserve data as stored in Excel, otherwise autmatic conversion (eg '01' will become number 1) ea : werkt niet...
    # remember pandas is 0-indexed
    # a list of index_col will generate a multi index (a dict of dicts)

    d_ds = df.to_dict(orient="index")




    return d_dds


if __name__ == "__main__":
    ## SIMPLE DICT
       # The “orientation” of the data. If the keys of the passed dict should be the columns of the resulting DataFrame, pass ‘columns’ (default). 
       #Otherwise if the keys should be rows, pass ‘index’.
    word2index={}
    df = pd.DataFrame.from_dict(word2index,orient='index')  #index orientatie is meest logische voor een simpele dict
       
    ##DICT OF DICTS
    # By column =>the dict is sorted on key and then each dict will become a column, 
        ##with the columnname = key of outerdict
        ##index = key of innerdict (columns without that index will get NaN)
        ##remark : get the same effect by using a series with index instead of innerdict eg pd.Series([1., 2., 3.], index=['a', 'b', 'c']
        ## control the order and selection via index and column parameter !!! cfr df2
    d_d_1 = {'a':{'x':1,'y':2}
            ,'c':{'x':3,'y':4,'z':5}
            ,'b':{'x':6,'y':7}}

    df = pd.DataFrame(d_d_1)
    df2 = pd.DataFrame(d_d_1,index=['x','z'],columns=['c','a']) ##selecting and sorting!!
    df3 = df = pd.DataFrame.from_dict(d_d_1,orient='columns') #completely the same result
    df3 = df = pd.DataFrame.from_dict(d_d_1) #completely the same result
    df4 = pd.DataFrame.from_dict(d_d_1,orient='index')  ## by row !!!!so transpose result #1 problem : you have no 'columns'-argument, so no control over this classmethod DataFrame.from_dict(data, orient='columns', dtype=None)

    ##DICT OF LISTS
    #By column => same as dict of dicts, but 
        ##lists have to be the same length (because lists don't have indexes...)
        ## the indexes will be range (so numbers)
    d_l_1 = {'a':[1,2]
            ,'c':[3,4]
            ,'b':[6,7]}

    df = pd.DataFrame(d_l_1)

    ##LIST OF DICTS
        ## by row => every dict is a row, 
        ##      index = a counter, 
        ##      columnnames = keys of dict (automatically sorted)
    l_d_1 = [{'x':1,'a':2}
            ,{'x':3,'a':4,'z':5}
            ,{'x':6,'a':7}]

    df = pd.DataFrame(l_d_1)


    ##LIST OF ITEMS
        ##by colummn: the tuple acts as the inner dict (cfr dict of dicts) (first element = key, second = the columnvalues),
        ##      but now this doesn't get sorted (because no outerdict with keys to sort on), and the index = a counter
        ##    remark : This can be useful for constructing a DataFrame with the columns in a particular order
        ##              without having to pass an explicit list of columns
        ##    you can transpose but then you need to give the columnnames of course (no counter for columnnames allowed)
    df= pd.DataFrame.from_dict(dict([('A', [1, 2, 3]), ('B', [4, 5, 6])]))

    df2= pd.DataFrame.from_dict(dict([('A', [1, 2, 3]), ('B', [4, 5, 6])]), 
                                  orient='index',columns=['one', 'two', 'three'])




    ###SELECTING VALUES
    ####################
    # df_genre_moods.Reverent                               #select a column by name
    # df_genre_moods['Reverent']                            #select a column by name
    # df_genre_moods[0:3]                                   #select rows by index
    # df_genre_moods[0:1]                                   #select 1 row by index
    # df_genre_moods['Blues':'Blues']                       #select 1 row by name
    # df_genre_moods.loc['Blues']                           #select 1 row by name
    # df_genre_moods.loc[:,'Reverent']                      #select 1 column by name
    # df_genre_moods.loc['Blues','Reverent']                #select 1 value by name
    # df_genre_moods.loc['Blues',['Reverent','Sad'] ]       #select 2 values by name
    # df_genre_moods.loc['Blues':'Rap',['Reverent','Sad'] ] #select a block by name
    # df_genre_moods.iloc[2]                                #select a row by index
    # df_genre_moods.iloc[2:5,4:8]                          #select a block by index
    # df_genre_moods.iloc[:,4:8]                            #select columns by index
    # df_genre_moods.iloc[[1,3,8],[4,7]]                    #select a selection by index
    # df_genre_moods.iloc[5,8]                              #select a value by index
    # df_genre_moods.iat[5,8]                               #select a value by index

    # df_genre_moods[df_genre_moods.Rustic>0]               #select rows by condition on a column
    # df_genre_moods[df_genre_moods>0]                      #select rows by condition on all values (NaN for not selected)
    # df_genre_moods[df_genre_moods.Rustic.isin([0.5,0.3])] #select rows by condition on a column




    #################

    # by row (from a dictionary)
    # rows = tuple (from an iterable)
    # columnnames = list
    # rowindex = counter (=automatic)
    def expand_grid(data_dict):

        rows = itertools.product(*data_dict.values()) 
        for row in rows:
            print(row)
        return pd.DataFrame.from_records(rows, columns=data_dict.keys())

    df = expand_grid(
            {'height': [60, 70],
             'weight': [100, 140, 180],
             'sex': ['Male', 'Female']})
    		 
    		

    ###input excel with datastructures as different tabs into pandas, and then into dicts
    #https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_excel.html?highlight=excel#pandas.read_excel
    f_ds = "/home/wth/Downloads/testinfra/INPUT/M030a/ds_agg_geom_data.xlsx"

    df = pd.read_excel(f_ds,sheet_name=1,header=0,index_col=[0,1],dtype=object)  
    #sheet_name None will process all tabes, 2 refers to third sheet, string as name is also possible
    #header=0 takes first row as header
    #dtype=object : Use object to preserve data as stored in Excel, otherwise autmatic conversion (eg '01' will become number 1)
    # remember pandas is 0-indexed
    # a list of index_col will generate a multi index (a dict of dicts)

    d_ds = df.to_dict(orient="index")
    # access the data as follows : d_ds[(1, 'cell_parent_0')]['cellID']
    #     so beware : the key is a tuple !! the value is a dict with the key = column name ! (so u can use multiple values as well = flexible
    #     so datastructure = d_t_d (dict of tuple to dict))

    #all tabs read simutaneously, leading to a dict of dataframes with keys = tabnames 
    dd_tab_df = pd.read_excel(f_ds,sheet_name=None,header=0,index_col=0,dtype=object) #dict of dataframes ! because of sheet-name=None
    df2_tab1 = df2['repl_embID']
    df2_tab2 = df2['repl_lab_cellID']  # this should be a multiindex, so we have to expand the index with the first column 
    df2_tab2_multi_index  = set_index(keys=df2_tab2.columns[0],append=True,inplace=False)  
    df2_tab3 = df2['repl_lab_cellID']
    d_repl_embID = df2_tab1.to_dict(orient="index")
    d_repl_lab_cellID = df2_tab2_multi_index.to_dict(orient="index")
    d_repl_lab_cellID = df2_tab1.to_dict(orient="index")
    #now c
