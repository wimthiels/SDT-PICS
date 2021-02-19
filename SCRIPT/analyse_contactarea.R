print("CONTACT AREA ANALYSIS#################################################################################################################")
#print('# 0-Fetch command line arguments------------') 
myArgs <- commandArgs(trailingOnly = TRUE)

#<<<MANUAL INPUT : REPLACE DATA_FOLDER TO LOCAL PATH !>>>>>
DATA_FOLDER <- "/home/wth/Downloads/testinfra/"
SCRIPT_FOLDER <- "/home/wth/Downloads/SYNC/"

if (length(myArgs)==0) {
  print("Input arguments filled with default values:")
  
  input_file <- paste0(DATA_FOLDER,"OUTPUT/7cellAgg/agg_tc.xlsx")
  output_folder <- paste0(DATA_FOLDER,"OUTPUT/7cellAgg/")

} else {
  input_file <- myArgs[1]
  output_folder <- myArgs[2]
}
print(paste0("Input file: ", input_file))
print(paste0("output_folder: ", output_folder))

# 0- load packages ---------
print(">>> 0 - Load packages <<< ")
packages <- c("ggplot2", "dplyr","nlme","lme4","lattice","writexl","stargazer","Rmisc","stringr")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages()))) 
}

library(ggplot2)
library(lme4)
library(lattice)
library(magrittr)
library(Matrix)
library(MASS)
#library(car)
library(pbkrtest)
library(stats)
library(readxl)
#library(nlme)
library(lattice)
#library(doBy)
library(nlme)
library(reshape2)
library(writexl)
library(Rmisc)
library(dplyr)
library(devtools)
library(lmerTest)
#library(rstudioapi) 
library(stringr) 

source(paste0(SCRIPT_FOLDER,"helper_functions.R") )
# source(paste0(dirname(getSourceEditorContext()$path) ,"/helper_functions.R"))


# Conventions : 
# - values expressed in um²
# - area = 0 means the cells exist together at that timestep, but do not make contact
# - area = NA means one of the cells does not yet exists at that timestep
# - mean : when averaging a NA-value will be ignored (not considered = 0)
print(">>> 0 - read data / enrich / select<<< ")
# 0 - read data / enrich / select -------------
input_data_ca_full <-read_excel(input_file,sheet = "contactinfo",na="NaN")
# join tc-data
tc_xml<-read_excel(input_file,sheet = "tc_xml",na="NaN")
input_data_ca_full <-merge(x = input_data_ca_full, y = tc_xml, by = c("tc_label"), all.x = TRUE)


## (filter, convert to um², pad with zeros, apply NA) ------
# filter
l_cell_name <- c('ABpl','ABal','ABpr','ABar','E','MS','P2','C','P3')
df_ca <- subset(input_data_ca_full,(input_data_ca_full$name1 %in% l_cell_name) & (input_data_ca_full$name2 %in% l_cell_name) )

# convert to um²
df_ca$contactarea <- df_ca$contactarea * 1.0e12 
max_ca = max(df_ca$contactarea,na.rm=T)
#max_ca = 4.e-10

# get mean contactarea (per cell pair) and drop the original columns
# checked : in extreme cases the delta is 30um², but usually quite close
get_reverse_ca = function(x, output) {
  cond1 <- df_ca$repID == x['repID']
  cond2 <- df_ca$nb_stack_analysis == x['nb_stack_analysis']
  cond3 <- df_ca$name1 == x['name2']
  cond4 <- df_ca$name2 == x['name1']
  df_ca[cond1&cond2&cond3&cond4,]$contactarea[1]
}

df_ca$contactarea_rev <- apply(df_ca, 1, get_reverse_ca)
df_ca$contactarea_mean <- rowMeans(cbind(df_ca$contactarea, df_ca$contactarea_rev), na.rm=FALSE)
df_ca <- df_ca[c('repID','nb_stack_analysis','name1','name2','contactarea_mean')]

# pad with zero values for all combinations
pad_zero_ca_values = function(df_ca) {
  #first create a df with all zero values over all cell1 / cell 2 / t combinations
  df_ca_pad <- as.data.frame(expand.grid(repID=sort(unique(df_ca$repID)),
                                         nb_stack_analysis= sort(unique(df_ca$nb_stack_analysis)), 
                                         name1=sort(unique(df_ca$name1)),
                                         name2=sort(unique(df_ca$name2))))
  df_ca_pad$contactarea_mean <- 0
  #second : add original data via left outer join
  df_ca_pad <-merge(x = df_ca_pad, y = df_ca, by = c("repID","name1","name2","nb_stack_analysis"), all.x = TRUE)
  df_ca_pad$contactarea_mean <- apply(df_ca_pad[c("contactarea_mean.x","contactarea_mean.y")], 1, max, na.rm=TRUE)
  df_ca_pad <- subset(df_ca_pad, select = -c(contactarea_mean.x,contactarea_mean.y) )
  df_ca_pad
}


# set NA values :when cell does not exist in a certain timestep : set to NA
#first find cell birth and death of every cell
get_lifespan<- function(df_ca){
  df_lifespan <- df_ca %>% 
    group_by(repID,name1) %>%
    summarize(time_birth = min(nb_stack_analysis),time_death= max(nb_stack_analysis))
  df_lifespan
}

overlay_NA_ca<- function(df_ca_pad,df_lifespan) {
  #if a cell does not exist in a certain timestep then set contactarea to NA (instead of zero)
  for( i in rownames(df_lifespan) ) {
    row <- df_lifespan[i,]
    cond1 = df_ca_pad$repID==as.character(row['repID'])
    cond2 = df_ca_pad$name1==as.character(row['name1'])
    cond3 = df_ca_pad$name2==as.character(row['name1'])
    cond4 = df_ca_pad$nb_stack_analysis<as.numeric(row['time_birth'])
    cond5 = df_ca_pad$nb_stack_analysis>as.numeric(row['time_death'])
    cond = cond1 & (cond2 | cond3) & (cond4 | cond5)
    # assign(df_ca_pad[cond,'contactarea_mean'],value=NA,inherits = TRUE)
    df_ca_pad[cond,'contactarea_mean'] <- NA
  }
  df_ca_pad
}
# by(df_lifespan, 1:nrow(df_lifespan), function(row) {
#   cond1 = df_ca_pad$repID==as.character(row['repID'])
#   cond2 = df_ca_pad$name1==as.character(row['name1'])
#   cond3 = df_ca_pad$name2==as.character(row['name1'])
#   cond4 = df_ca_pad$t<as.numeric(row['time_birth'])
#   cond5 = df_ca_pad$t>as.numeric(row['time_death'])
#   cond = cond1 & (cond2 | cond3) & (cond4 | cond5)
#   # assign(df_ca_pad[cond,'contactarea_mean'],value=NA,inherits = TRUE)
#   df_ca_pad[cond,'contactarea_mean'] <- NA  #this does not change df_ca_pad
#}

df_ca_pad<- pad_zero_ca_values(df_ca)
df_lifespan <- get_lifespan(df_ca)
df_ca_pad <- overlay_NA_ca(df_ca_pad,df_lifespan)
df_ca <- df_ca_pad


# add column ca_perc_emb = express a percentage of total contactarea per embryo
l_total_ca_emb <-by(df_ca$contactarea_mean,df_ca$repID,sum,na.rm=TRUE)
l_total_ca_emb <- l_total_ca_emb / (5*2) #5 timesteps, and symmetric 

df_ca$ca_total_emb<- unlist(lapply(df_ca$repID, function(x) l_total_ca_emb[x]))
df_ca$ca_perc_emb <- unlist(lapply(split(df_ca, seq(nrow(df_ca))), function(x) x$contactarea_mean * 100 / l_total_ca_emb[x$repID]))
df_ca <- subset(df_ca, select = -c(ca_total_emb) )
# check <-by(df_ca$ca_perc_emb, df_ca$repID, sum, na.rm=TRUE)
# check #must add up to 100 * 5 * 2 = 1000 per embryo 
# check <-sum(subset(df_ca,repID=='emb_1' & t==2,select=c(ca_perc_emb))) 
# check # must be around 200 (because of symmetric relation)
summary(df_ca)

#aggregate data over repID 
# remember: NA are skipped, not taken up in the mean. this means for example
# if C appears at t=2 for 1 replicate, the C contacts will be filled in for t=2,
# but these values are entirely determined by that 1 replicate
df_ca_agg<- df_ca %>%
  group_by(name1,name2,nb_stack_analysis) %>% 
  summarise(contactarea_mean=mean(contactarea_mean,na.rm = TRUE),
            ca_perc_emb=mean(ca_perc_emb,na.rm = TRUE))
summary(df_ca_agg)

#saving data (for viz in python)
write_xlsx(list(df_ca = df_ca,df_ca_agg = df_ca_agg,df_lifespan=df_lifespan), 
           paste0(output_folder,"R_7cell.xlsx",sep=''))



# 3.2 explore the data-----------------
print("# plotting all contactareas of all cells over time:")
for (cellname_i in l_cell_name) {
  s_title =  sprintf("%s contacts",cellname_i)
  s_f_name1 = sprintf("<contactarea %s>OLS_facetwrap(cellID).jpg",cellname_i)
  s_f_name2 = sprintf("<contactarea %s>OLS_facetgrid(cellID-repID).jpg",cellname_i)
  ca_celli_all <- subset(input_data_ca_full,(input_data_ca_full$name1 == cellname_i) & (input_data_ca_full$name2 %in% l_cell_name) )
  max_ca = max(ca_celli_all$contactarea,na.rm=T)
  ca_celli_plot<- ggplot(ca_celli_all, aes(x=nb_stack_analysis, y=contactarea)) + geom_point() + ylim(0, max_ca*1.1) + xlim(1,5)
  ca_celli_facet <-ca_celli_plot +  facet_wrap(~name2) + stat_smooth(method="lm") + ggtitle(s_title)
  ca_celli_facet <- ca_celli_facet + aes(color=name2)
  ca_celli_facet
  ggsave(paste0(output_folder,s_f_name1))
  ca_celli_facetgrid <-ca_celli_plot +  facet_grid(repID~name2) + stat_smooth(method="lm") + ggtitle(s_title)
  ca_celli_facetgrid
  ggsave(paste0(output_folder,s_f_name2))
}

# make 1 plot with all contactareas, aggregated over repID, normalized
ca_cells_plot<- ggplot(df_ca_agg, aes(x=nb_stack_analysis, y=ca_perc_emb)) + geom_point() + ylim(0, NA) + xlim(1,5)
ca_cells_facetgrid <-ca_cells_plot +  facet_grid(name1~name2) + stat_smooth(method="lm") + ggtitle("All cell-cell contacts over time (averaged over replicates)")
ca_cells_facetgrid




