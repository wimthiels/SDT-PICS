print("CONTACT AREA ANALYSIS#################################################################################################################")

#<<<MANUAL INPUT : REPLACE DATA_FOLDER TO LOCAL PATH !>>>>>
DATA_FOLDER <- "/home/wth/Downloads/testinfra/"
SCRIPT_FOLDER <- "/home/wth/Downloads/SYNC/"

source(paste0(SCRIPT_FOLDER,"fetch_command_line_arguments.R"))
source(paste0(SCRIPT_FOLDER,"helper_functions.R"))
source(paste0(SCRIPT_FOLDER,"set_input_parameters.R"))
output_folder <- paste0(output_folder,"/analyse_contactarea/")
source(paste0(SCRIPT_FOLDER,"load_R_packages_standard.R"))

# 0- load packages ---------
print(">>> 0 - Load extra packages <<< ")
packages <- c("nlme","lattice","MASS")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages()))) 
}
library(lattice)
library(Matrix)
library(MASS)
library(pbkrtest)
library(nlme)
library(lmerTest)

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

#fix types
input_data_ca_full$repID <- as.factor(input_data_ca_full$repID)
input_data_ca_full$nb_stack_EMS_division <- as.integer(input_data_ca_full$nb_stack_EMS_division)
input_data_ca_full$res_t <- as.integer(input_data_ca_full$res_t)
input_data_ca_full$name1<- as.factor(input_data_ca_full$name1)
input_data_ca_full$name2<- as.factor(input_data_ca_full$name2)

#calculate time (=minutes from EMS division)
input_data_ca_full['t'] <- (input_data_ca_full$nb_stack_raw  - input_data_ca_full$nb_stack_EMS_division)  *  input_data_ca_full$res_t/60   #expressed in minutes,t=0 should be a reference time (e.g. EMS division)

## (filter, convert to um², pad with zeros, apply NA) ------
# select
if (!exists("l_cell_name") | anyNA(l_cell_name)  | length(l_cell_name)==0 ){
  l_cell_name <- unique(input_data_ca_full['name1'])[[1]]
  #l_cell_name <- c('ABal','ABar','ABpl','ABpr','E','MS','P2','C','P3') #also sets sort order
}
if (!exists("l_emb_ID") | anyNA(l_emb_ID) | length(l_emb_ID)==0){
  l_emb_ID <- unique(input_data_ca_full['embID'])[[1]]
  #l_emb_ID <- c('emb_1','emb_2','emb_3','emb_4', 'emb_5','emb_6','emb_7','emb_8','emb_9')
}
if (!exists("l_rep_ID") | anyNA(l_rep_ID)  | length(l_rep_ID)==0){
  l_rep_ID <- unique(input_data_ca_full['repID'])[[1]]
}
if (!exists("l_time_interval") | anyNA(l_time_interval)  | length(l_time_interval)==0){
  l_time_interval <- list(min(input_data_ca_full['t'])[[1]],max(input_data_ca_full['t'])[[1]])
}

print('Data will be filtered by  : ')
print(l_cell_name)
print(l_emb_ID)
print(l_rep_ID)
print(l_time_interval)

df_ca <- subset(input_data_ca_full,(input_data_ca_full$name1 %in% l_cell_name) &
                  (input_data_ca_full$name2 %in% l_cell_name) &
                  (input_data_ca_full$embID %in% l_emb_ID) &
                  (input_data_ca_full$repID %in% l_rep_ID) &
                  (input_data_ca_full$t >= as.integer(l_time_interval[[1]]) )&  
                  (input_data_ca_full$t <= as.integer(l_time_interval[[2]]))
                )

print(paste('The full dataframe contains ',nrow(input_data_ca_full),' rows. The filtered data contains ',nrow(df_ca),' rows.'))
print('summary of df_ca after selection, before padding : ')
print(summary(df_ca))
# convert to um²
df_ca$contactarea <- df_ca$contactarea * 1.0e12 
max_ca = max(df_ca$contactarea,na.rm=T)
#max_ca = 4.e-10

#set flags
if (length(unique(df_ca['category'])[[1]]) > 1){
  flag_single_category = FALSE
} else {
  flag_single_category = TRUE
}


# get mean contactarea (per cell pair) and drop the original columns
# checked : in extreme cases the delta is 30um², but usually quite close
get_reverse_ca = function(x,result) {
  cond1 <- df_ca$repID == x['repID']
  cond2 <- df_ca$tc_label == x['tc_label']  #cond2 <- df_ca$t == x['t'] does not work ?!?
  cond3 <- df_ca$name1 == x['name2']
  cond4 <- df_ca$name2 == x['name1']
  df_ca[cond1&cond2&cond3&cond4,]$contactarea[1]
}

df_ca$contactarea_rev <- apply(df_ca, 1, get_reverse_ca)
df_ca$contactarea_mean <- rowMeans(cbind(df_ca$contactarea, df_ca$contactarea_rev), na.rm=FALSE)


#sanity check : lets calculate the total surface area of each cell and plot it
# it should stay more or less constant, and should not differ widely between replicates
df_ca_cell<- df_ca %>%
  group_by(name1,repID,t) %>% 
  summarise(contactarea_sum=sum(contactarea_mean,na.rm = TRUE))
df_ca_cell <- df_ca_cell[df_ca_cell$contactarea_sum!=0,]
sanity_plot <- ggplot(df_ca_cell, aes(x=t, y=contactarea_sum, color=name1)) +geom_line() +
  labs(title="Sanity check : total surface area per cell over time")
sanity_plot
ggsave(paste0(output_folder,"_Sanity-check_total_surface area per cell over time.jpg"))

df_ca_embryo<- df_ca %>%
  group_by(repID,t) %>% 
  summarise(contactarea_sum=sum(contactarea_mean,na.rm = TRUE))
df_ca_embryo<- df_ca_embryo[df_ca_embryo$contactarea_sum!=0,]
sanity_plot2 <- ggplot(df_ca_embryo, aes(x=t, y=contactarea_sum, color=repID)) +geom_point() +
  labs(title="Sanity check : total surface area per embryo (only selected cells !) over time")
sanity_plot2
ggsave(paste0(output_folder,"_Sanity-check_total_surface embryo-only selected cells - per cell over time.jpg"))

df_ca <- df_ca[c('repID','t','name1','name2','contactarea_mean')]


# pad with zero values for all combinations
pad_zero_ca_values = function(df_ca) {
  #first create a df with all zero values over all cell1 / cell 2 / t combinations
  df_ca_pad <- as.data.frame(expand.grid(repID=sort(unique(df_ca$repID)),
                                         t= sort(unique(df_ca$t)), 
                                         name1=sort(unique(df_ca$name1)),
                                         name2=sort(unique(df_ca$name2))))
  df_ca_pad$contactarea_mean <- 0
  #second : add original data via left outer join
  df_ca_pad <-merge(x = df_ca_pad, y = df_ca, by = c("repID","name1","name2","t"), all.x = TRUE)
  df_ca_pad$contactarea_mean <- apply(df_ca_pad[c("contactarea_mean.x","contactarea_mean.y")], 1, max, na.rm=TRUE)
  df_ca_pad <- subset(df_ca_pad, select = -c(contactarea_mean.x,contactarea_mean.y) )
  df_ca_pad
}


# set NA values :when cell does not exist in a certain timestep : set to NA
#first find cell birth and death of every cell
get_lifespan<- function(df_ca){
  df_lifespan <- df_ca %>% 
    group_by(repID,name1) %>%
    summarize(time_birth = min(t),time_death= max(t))
  df_lifespan
}

overlay_NA_ca<- function(df_ca_pad,df_lifespan) {
  #if a cell does not exist in a certain timestep then set contactarea to NA (instead of zero)
  for( i in rownames(df_lifespan) ) {
    row <- df_lifespan[i,]
    cond1 = df_ca_pad$repID==as.character(row['repID'])
    cond2 = df_ca_pad$name1==as.character(row['name1'])
    cond3 = df_ca_pad$name2==as.character(row['name1'])
    cond4 = df_ca_pad$t<as.numeric(row['time_birth'])
    cond5 = df_ca_pad$t>as.numeric(row['time_death'])
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
print('summary of df_ca after selection, AFTER padding (check for NAs in contacarea !): ')
summary(df_ca)



#aggregate data over repID 
# remember: NA are skipped, not taken up in the mean. this means for example
# if C appears at t=2 for 1 replicate, the C contacts will be filled in for t=2,
# but these values are entirely determined by that 1 replicate
df_ca_agg<- df_ca %>%
  group_by(name1,name2,t) %>% 
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
  s_f_name1 = sprintf("_contactarea %s_OLS_facetwrap(cellname).jpg",cellname_i)
  s_f_name2 = sprintf("_contactarea %s_OLS_facetgrid(repid-cellname).jpg",cellname_i)
  ca_celli_all <- subset(input_data_ca_full,(input_data_ca_full$name1 == cellname_i) 
                         & (input_data_ca_full$name2 %in% l_cell_name) 
                         & (input_data_ca_full$name2 != cellname_i))
  max_ca = max(ca_celli_all$contactarea,na.rm=T)
  ca_celli_plot<- ggplot(ca_celli_all, aes(x=t, y=contactarea)) + geom_point(alpha=0.5) + ylim(0, max_ca*1.1) 
  ca_celli_facet <-ca_celli_plot +  facet_wrap(~name2) + ggtitle(s_title) #+ stat_smooth(method="lm") 
  ca_celli_facet <- ca_celli_facet + aes(color=repID)
  ca_celli_facet
  ggsave(paste0(output_folder,s_f_name1))
  ca_celli_facetgrid <-ca_celli_plot +  facet_grid(repID~name2) + ggtitle(s_title) # + stat_smooth(method="lm") 
  ca_celli_facetgrid
  ggsave(paste0(output_folder,s_f_name2))
  if (!flag_single_category){
    s_f_name3 = sprintf("_contactarea %s_OLS_facetgrid(category-cellname).jpg",cellname_i)
    ca_celli_facetgrid <-ca_celli_plot +  facet_grid(category~name2) + ggtitle(s_title)   + stat_smooth(method="lm")
    ca_celli_facetgrid 
    ggsave(paste0(output_folder,s_f_name3))
    
    s_f_name4 = sprintf("_contactarea %s_OLS_facetwrap(cellname-colour-category).jpg",cellname_i)
    ca_celli_facetgrid <-ca_celli_plot +  facet_wrap(~name2) + stat_smooth(method="lm") + ggtitle(s_title)
    ca_celli_facetgrid <- ca_celli_facetgrid + aes(color=category,alpha=0.5)
    ca_celli_facetgrid 
    ggsave(paste0(output_folder,s_f_name4))
    
  }
}


# make 1 plot with all contactareas, aggregated over repID, normalized
ca_cells_plot<- ggplot(df_ca_agg, aes(x=t, y=ca_perc_emb)) + geom_point() + ylim(0, NA) + xlim(1,5)
ca_cells_facetgrid <-ca_cells_plot +  facet_grid(name1~name2) + stat_smooth(method="lm") + ggtitle("All cell-cell contacts over time (averaged over replicates)")
ca_cells_facetgrid




