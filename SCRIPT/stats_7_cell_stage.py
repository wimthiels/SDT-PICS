"""
Created on 11 April 2020
a script that does some descriptive statistics on geometry data of the 7 cell stage (volume)
also some hypothesis testing to check for significant volume differences between cells (E vs MS, and the AB cells)

@author: wimth
"""

import os
from plotnine import ggplot, geom_point, geom_boxplot, stat_smooth, facet_wrap,save_as_pdf_pages,aes,scale_x_discrete#https://plotnine.readthedocs.io/en/stable/index.html
import pandas as pd  
from pathlib import Path
import scipy.stats as stats
import researchpy as rp
import statsmodels.api as sm
from statsmodels.formula.api import ols




if __name__ == "__main__":
    output_folder = Path('/home/wth/Downloads/testinfra/OUTPUT/M027')
    output_folder.mkdir(parents=True,exist_ok=True) 
    input_file = Path('/home/wth/Downloads/testinfra/INPUT/M027/volume_7_cell_stage.csv')
    f = open(str( output_folder / "stats_analysis_7_cell_stage.txt") ,'w')
    #read in volume measurments into a df (from /home/wth/Downloads/SYNC/E_MS_symmetry_data.xlsx)
    df = pd.read_csv(input_file,index_col=0)
    df = df.stack().to_frame()  #stack give a dataSeries(with multi index)
    df.reset_index(inplace=True) #turns multiindex into columns
    df.columns=['cellID','embryo','volume']
    print(df,file=f)

    #part 1 : descriptive statistics 
    gg_box_volume = ggplot(df) +   \
        geom_boxplot(mapping = aes(x='factor(cellID)', y='volume')) + \
        scale_x_discrete(name='month') 
    gg_box_volume.save(filename='boxplot_volume.pdf',path=str(output_folder),verbose=False)

    # Gettin summary statistics
    print(rp.summary_cont(df['volume'].groupby(df['cellID'])),file=f)
    # save_as_pdf_pages(gg_box_volume,file_name = 'boxplot_volume',path=output_folder,verbose=True)

    #PARTII HYPOTHESIS TESTING   (https://pythonfordatascience.org/anova-python/)
    # difference in volume AB cells ?
    print("HYPOTHESIS TEST : ALL CELLS-------------------------------------------------------",file=f)
    model_allcells_volume = ols('volume ~ C(cellID)', data=df).fit()
    print(model_allcells_volume.summary(),file=f)
    aov_table_all_cells = sm.stats.anova_lm(model_allcells_volume, typ=2)
    print(aov_table_all_cells,file=f)
    
    # df.loc[(df['cellID']=='ABal') | (df['cellID']=='ABar') |(df['cellID']=='ABpl') | (df['cellID']=='ABpr')  ]
    print("HYPOTHESIS TEST : AB CELLS--------------------------------------------------------",file=f)
    df_AB = df.loc[df['cellID'].str.startswith('AB')]
    model_AB_volume = ols('volume ~ C(cellID)', data=df_AB).fit()
    print(model_AB_volume.summary(),file=f)
    aov_table_AB = sm.stats.anova_lm(model_AB_volume, typ=2)
    print(aov_table_AB,file=f)

    # df_EMS = df.loc[~df['cellID'].str.startswith('AB')]
    print("HYPOTHESIS TEST : EMS ------------------------------------------------------------",file=f)
    df_EMS = df.loc[(df['cellID']=='E') | (df['cellID']=='MS') ]
    model_EMS_volume = ols('volume ~ C(cellID)', data=df_EMS).fit()
    print(model_EMS_volume.summary(),file=f)
    aov_table_EM = sm.stats.anova_lm(model_EMS_volume, typ=2)
    print(aov_table_EM,file=f)

#LOG---------------------------------------------------------------------------------------
#    cellID       embryo       volume
# 0      P2          TL1  3880.124887
# 1      P2  TL2_emb2_t2  3399.959019
# 2      P2       TL3_t2  3807.286568
# 3      P2  TL9_emb2_t8  3111.135547
# 4    ABar          TL1  3477.016211
# 5    ABar  TL2_emb2_t2  3265.999599
# 6    ABar       TL3_t2  3568.638810
# 7    ABar  TL9_emb2_t8  2513.198513
# 8    ABpr          TL1  3449.559332
# 9    ABpr  TL2_emb2_t2  2864.767532
# 10   ABpr       TL3_t2  2932.807394
# 11   ABpr  TL9_emb2_t8  2393.291320
# 12     MS          TL1  3140.363867
# 13     MS  TL2_emb2_t2  2330.607235
# 14     MS       TL3_t2  3242.188675
# 15     MS  TL9_emb2_t8  2544.086222
# 16      E          TL1  2848.317161
# 17      E  TL2_emb2_t2  2298.304885
# 18      E       TL3_t2  2641.500746
# 19      E  TL9_emb2_t8  1979.587506
# 20   ABal          TL1  3003.280165
# 21   ABal  TL2_emb2_t2  3349.843176
# 22   ABal       TL3_t2  3546.252932
# 23   ABal  TL9_emb2_t8  3323.746160
# 24   ABpl          TL1  3442.904375
# 25   ABpl  TL2_emb2_t2  2840.974834
# 26   ABpl       TL3_t2  2930.733579
# 27   ABpl  TL9_emb2_t8  3064.991466


#         N         Mean          SD          SE    95% Conf.     Interval
# cellID                                                                  
# ABal    4  3305.780608  224.794519  112.397260  3051.401663  3560.159554
# ABar    4  3206.213283  479.072146  239.536073  2664.092116  3748.334451
# ABpl    4  3069.901063  265.160569  132.580284  2769.843633  3369.958494
# ABpr    4  2910.106395  432.310938  216.155469  2420.900515  3399.312274
# E       4  2441.927575  382.695598  191.347799  2008.866738  2874.988411
# MS      4  2814.311500  445.862662  222.931331  2309.770401  3318.852598
# P2      4  3549.626505  360.690341  180.345170  3141.466961  3957.786049
# HYPOTHESIS TEST : ALL CELLS------------
#                             OLS Regression Results                            
# ==============================================================================
# Dep. Variable:                 volume   R-squared:                       0.508
# Model:                            OLS   Adj. R-squared:                  0.368
# Method:                 Least Squares   F-statistic:                     3.615
# Date:                Sat, 11 Apr 2020   Prob (F-statistic):             0.0128
# Time:                        17:49:07   Log-Likelihood:                -202.05
# No. Observations:                  28   AIC:                             418.1
# Df Residuals:                      21   BIC:                             427.4
# Df Model:                           6                                         
# Covariance Type:            nonrobust                                         
# =====================================================================================
#                         coef    std err          t      P>|t|      [0.025      0.975]
# -------------------------------------------------------------------------------------
# Intercept          3305.7806    190.168     17.383      0.000    2910.304    3701.257
# C(cellID)[T.ABar]   -99.5673    268.939     -0.370      0.715    -658.856     459.721
# C(cellID)[T.ABpl]  -235.8795    268.939     -0.877      0.390    -795.168     323.409
# C(cellID)[T.ABpr]  -395.6742    268.939     -1.471      0.156    -954.963     163.614
# C(cellID)[T.E]     -863.8530    268.939     -3.212      0.004   -1423.142    -304.564
# C(cellID)[T.MS]    -491.4691    268.939     -1.827      0.082   -1050.758      67.820
# C(cellID)[T.P2]     243.8459    268.939      0.907      0.375    -315.443     803.135
# ==============================================================================
# Omnibus:                        2.319   Durbin-Watson:                   3.140
# Prob(Omnibus):                  0.314   Jarque-Bera (JB):                1.473
# Skew:                          -0.304   Prob(JB):                        0.479
# Kurtosis:                       2.055   Cond. No.                         7.87
# ==============================================================================

# Warnings:
# [1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#                  sum_sq    df         F    PR(>F)
# C(cellID)  3.137328e+06   6.0  3.614698  0.012761
# Residual   3.037778e+06  21.0       NaN       NaN
# HYPOTHESIS TEST : AB CELLS------------
# /home/wth/anaconda3/lib/python3.7/site-packages/scipy/stats/stats.py:1450: UserWarning: kurtosistest only valid for n>=20 ... continuing anyway, n=16
#   "anyway, n=%i" % int(n))
#                             OLS Regression Results                            
# ==============================================================================
# Dep. Variable:                 volume   R-squared:                       0.180
# Model:                            OLS   Adj. R-squared:                 -0.025
# Method:                 Least Squares   F-statistic:                    0.8783
# Date:                Sat, 11 Apr 2020   Prob (F-statistic):              0.480
# Time:                        17:49:07   Log-Likelihood:                -114.86
# No. Observations:                  16   AIC:                             237.7
# Df Residuals:                      12   BIC:                             240.8
# Df Model:                           3                                         
# Covariance Type:            nonrobust                                         
# =====================================================================================
#                         coef    std err          t      P>|t|      [0.025      0.975]
# -------------------------------------------------------------------------------------
# Intercept          3305.7806    183.243     18.040      0.000    2906.529    3705.032
# C(cellID)[T.ABar]   -99.5673    259.144     -0.384      0.708    -664.194     465.059
# C(cellID)[T.ABpl]  -235.8795    259.144     -0.910      0.381    -800.506     328.747
# C(cellID)[T.ABpr]  -395.6742    259.144     -1.527      0.153    -960.301     168.952
# ==============================================================================
# Omnibus:                        0.912   Durbin-Watson:                   2.591
# Prob(Omnibus):                  0.634   Jarque-Bera (JB):                0.554
# Skew:                          -0.435   Prob(JB):                        0.758
# Kurtosis:                       2.732   Cond. No.                         4.79
# ==============================================================================

# Warnings:
# [1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#                  sum_sq    df         F    PR(>F)
# C(cellID)  3.539055e+05   3.0  0.878321  0.479573
# Residual   1.611737e+06  12.0       NaN       NaN
# HYPOTHESIS TEST : EMS ------------
# /home/wth/anaconda3/lib/python3.7/site-packages/scipy/stats/stats.py:1450: UserWarning: kurtosistest only valid for n>=20 ... continuing anyway, n=8
#   "anyway, n=%i" % int(n))
#                             OLS Regression Results                            
# ==============================================================================
# Dep. Variable:                 volume   R-squared:                       0.211
# Model:                            OLS   Adj. R-squared:                  0.080
# Method:                 Least Squares   F-statistic:                     1.607
# Date:                Sat, 11 Apr 2020   Prob (F-statistic):              0.252
# Time:                        17:49:07   Log-Likelihood:                -58.436
# No. Observations:                   8   AIC:                             120.9
# Df Residuals:                       6   BIC:                             121.0
# Df Model:                           1                                         
# Covariance Type:            nonrobust                                         
# ===================================================================================
#                       coef    std err          t      P>|t|      [0.025      0.975]
# -----------------------------------------------------------------------------------
# Intercept        2441.9276    207.741     11.755      0.000    1933.605    2950.251
# C(cellID)[T.MS]   372.3839    293.790      1.268      0.252    -346.493    1091.261
# ==============================================================================
# Omnibus:                        3.393   Durbin-Watson:                   3.177
# Prob(Omnibus):                  0.183   Jarque-Bera (JB):                0.938
# Skew:                          -0.125   Prob(JB):                        0.626
# Kurtosis:                       1.341   Cond. No.                         2.62
# ==============================================================================

# Warnings:
# [1] Standard Errors assume that the covariance matrix of the errors is correctly specified.
#                  sum_sq   df         F    PR(>F)
# C(cellID)  2.773396e+05  1.0  1.606604  0.251939
# Residual   1.035748e+06  6.0       NaN       NaN
