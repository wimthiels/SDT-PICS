print("SPHERICITY ANALYSIS#################################################################################################################")
#print('# 0-Fetch command line arguments------------') 
myArgs <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(myArgs)==0) {
  # input_file <- "/home/wth/Downloads/testinfra/OUTPUT/fmnAgg/agg_tc.xlsx"
  # output_folder <- "/home/wth/Downloads/testinfra/OUTPUT/fmnAna"
  print("Input arguments filled with default values:")
  input_file <- "/home/wth/Downloads/testinfra/OUTPUT/7cellAgg/agg_tc.xlsx"
  input_file_fmn <- "/home/wth/Downloads/testinfra/OUTPUT/fmnAgg/agg_tc.xlsx"
  output_folder <- "/home/wth/Downloads/testinfra/OUTPUT/7cellAna/analyse_sphericity/"
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
# library(rstudioapi) 
library(stringr) 

source("/home/wth/Downloads/SYNC/helper_functions.R") #todo: replace this via parm passing
# source(paste0(dirname(getSourceEditorContext()$path) ,"/helper_functions.R")) 

# 0 - read data / enrich / select -------------
print(">>> 0 - read data / enrich / select<<< ")
geom_info_full <-read_excel(input_file,sheet = "geom_info",na="NaN")
print('geom_info_full (raw) before enrichment -> ')
str(geom_info_full)
if (exists("input_file_fmn")){
  geom_info_fmn <-read_excel(input_file_fmn,sheet = "geom_info",na="NaN")
  geom_info_full2 <- rbind(geom_info_full,geom_info_fmn)
}
# join tc-data
tc_xml<-read_excel(input_file,sheet = "tc_xml",na="NaN")
geom_info_full <-merge(x = geom_info_full, y = tc_xml, by = c("tc_label"), all.x = TRUE)

#fix types
geom_info_full$res_t <- as.numeric(geom_info_full$res_t)
geom_info_full$repID <- as.character(geom_info_full$repID)
geom_info_full$nb_stack_EMS_division <- as.numeric(geom_info_full$nb_stack_EMS_division)


#enrich
#we define relative to EMS-division 
geom_info_full['t'] <- (geom_info_full$nb_stack - 1 - geom_info_full$nb_stack_EMS_division)  *  geom_info_full$res_t/60   #expressed in minutes,t=0 should be a reference time (e.g. EMS division)
if ("nb_pix_last_z" %in% colnames(geom_info_full)){} else {geom_info_full$nb_pix_last_z <- 0 }

# To get the time to the next division we need manual input 
      #future update : derive automatically instead of manual input
file_time_until_division <- '/home/wth/Downloads/SYNC/dds_7cell_geom.xlsx'
df_time_until_division <-read_excel(file_time_until_division ,sheet = "1",na="NaN")
if (exists("df_time_until_division")){
  df_time_until_division$repID<-str_remove(df_time_until_division$replicate, "^0+") # some ad-hoc parsing to match repID
  names(df_time_until_division)[names(df_time_until_division) == "cellID"] <- "cell_name"
  df_time_until_division = subset(df_time_until_division,select = c(cell_name,repID,delta_t_division))
  df_time_until_division$delta_t_division <- as.numeric(df_time_until_division$delta_t_division)
  geom_info_full <-merge(x = geom_info_full, y = df_time_until_division, by = c("repID","cell_name"), all.x = TRUE)
  geom_info_full$delta_t_division <- as.numeric(geom_info_full$delta_t_division)
  geom_info_full$phase_division <- 'after division'
  geom_info_full$t_division <- (geom_info_full$delta_t_division * 3 ) + geom_info_full$t #3minute time intervals
  geom_info_full$phase_division[geom_info_full$t_division < 0] <- 'before division'
  
}else{
  print('To get the time to the next division we need manual input. No file provided !')
  geom_info_full$phase_division <- 'No division info provided'
  geom_info_full$t_division <- 0
}




# select
l_cell_name_all <- unique(geom_info_full['cell_name'])[[1]]
l_cell_name <- c('ABal','ABar','ABpl','ABpr','E','MS','P2','C','P3') #also sets sort order
l_emb_ID <- c('emb_1','emb_2','emb_3','emb_4', 'emb_5','emb_6','emb_7')
l_rep_ID <- unique(geom_info_full['repID'])[[1]]

geom_info<- subset(geom_info_full,geom_info_full$cell_name %in% l_cell_name &
                     geom_info_full$embID %in% l_emb_ID)

attach(geom_info)
str(geom_info)

## 1 -  exploratory analysis (fig S3)--------
print('## 1 - exploratory analysis (fig S3)----')
summary(geom_info)

#   OLS for every repID/cell_name combination (sphericity)
sphericityplot<- ggplot(geom_info, aes(x=t, y=sphericity)) +geom_point() + 
  ylim(min(sphericity,na.rm = T), max(sphericity,na.rm = T))

sphericityplot<-sphericityplot +  facet_grid(repID~cell_name) + stat_smooth(method="lm")
sphericityplot <- sphericityplot  + theme(axis.title.x=element_text(size=16),
                                          axis.title.y=element_text(size=16),
                                          legend.title = element_blank())
#sphericityplot<- sphericityplot + scale_x_continuous(name="t", breaks=seq(-10,10,10),limits=c(-15, 15))
sphericityplot
ggsave(paste0(output_folder,"<sphericity>OLS_facet_grid(repID-cell_name).jpg"),dpi='retina')

#   OLS for every repID/cell_name combination (bounding box volume)
geom_info$bounding_box_vol <- geom_info$extents_length * geom_info$extents_width * geom_info$extents_height
bounding_box_plot<- ggplot(geom_info, aes(x=t, y=bounding_box_vol)) +geom_point() + 
  ylim(min(geom_info$bounding_box_vol,na.rm = T), max(geom_info$bounding_box_vol,na.rm = T))
bounding_box_plot<-bounding_box_plot +  facet_grid(repID~cell_name) + stat_smooth(method="lm")
bounding_box_plot
ggsave(paste0(output_folder,"<bounding box>OLS_facet_grid(repID-cell_name).jpg"))

#   fig S3 - sphericityplot OLS (fig5c) for every repID/cell_name combination (sphericity vs time to cell division)----
sphericityplot<- ggplot(geom_info, aes(x=t, y=sphericity,colour=phase_division)) +geom_point() + 
  ylim(min(sphericity,na.rm = T), max(sphericity,na.rm = T))
sphericityplot<- sphericityplot + scale_x_continuous(name="\u0394t cell division (min)", breaks=seq(-10,10,10),limits=c(-15, 15))
sphericityplot<-sphericityplot +  facet_grid(repID~cell_name) + stat_smooth(method="lm",se=FALSE)
sphericityplot <-sphericityplot + labs(y="sphericity",title="")
sphericityplot <- sphericityplot  + theme(axis.title.x=element_text(size=16),
                                          axis.title.y=element_text(size=16),
                                          legend.title = element_blank(),
                                          legend.position="bottom",legend.direction = "horizontal")
sphericityplot
ggsave(paste0(output_folder,"<sphericity>OLS_facet_grid[repID-cell_name]_vs_delta_t_division.jpg"),dpi='retina')


# Same graph, but now aggregated over replicates and sync on time of dvision
cell_spher_prediv <- geom_info %>% group_by(t_division,cell_name) %>% 
  filter(phase_division == 'before division') %>% 
  summarise(sphericity =mean(sphericity,na.rm=TRUE)) 
cell_spher_prediv$phase_division <- 'before division'
cell_spher_postdiv <- geom_info %>% group_by(t_division,cell_name) %>% 
  filter(phase_division == 'after division') %>% 
  summarise(sphericity =mean(sphericity,na.rm=TRUE)) 
cell_spher_postdiv$phase_division <- 'after division'
cell_spher <- rbind(cell_spher_prediv,cell_spher_postdiv)


sphericityplot<- ggplot(cell_spher, 
                        aes(x=t_division, y=sphericity, 
                            group=phase_division, color = cell_name,
                        )) 
sphericityplot <-sphericityplot + geom_point() + ylim(0.7, NA) +    ylim(min(sphericity,na.rm = T), max(sphericity,na.rm = T))
# sphericityplot<-sphericityplot + stat_smooth(method="lm", se=FALSE,size=3)

sphericityplot<-sphericityplot + 
  stat_smooth(data=cell_spher_prediv,method="lm", se=FALSE,
              inherit.aes=FALSE,
              aes(x=t_division,y=sphericity,
                  colour=cell_name),size=2) 
sphericityplot<-sphericityplot + 
  geom_smooth(data=cell_spher_postdiv,method="lm", se=FALSE,
              inherit.aes=FALSE,
              aes(x=t_division,y=sphericity,
                  colour=cell_name),size=2) 

sphericityplot<-  sphericityplot + geom_vline(xintercept=0)

sphericityplot <-sphericityplot  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                         panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                         axis.title.x=element_blank(),axis.title.y=element_blank(),
                                         axis.text.x= element_text(face="bold", size=16),
                                         axis.text.y= element_text(face="bold", size=16))
# sphericityplot <- sphericityplot  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))

sphericityplot <- sphericityplot + theme(legend.title = element_blank(),
                                         legend.position="bottom",legend.direction = "horizontal",
                                         legend.text = element_text(size = "large"))
sphericityplot
ggsave(paste0(output_folder,"<sphericity>delta_t_division_turning point.pdf"),device = "pdf")

sphericityplot_stripped <- strip_plot(sphericityplot)
sphericityplot_stripped
ggsave(paste0(output_folder,"<sphericity>delta_t_division_turning point(stripped).pdf"),device = "pdf")


#2.2 fitting ----
# model ABpl sphericity
geom_info_ABpl <- geom_info[geom_info$cell_name == 'ABpl',]
lmer_cell_spher_ABpl <- lmer(geom_info_ABpl$sphericity~ 1 + geom_info_ABpl$t  + (1 | geom_info_ABpl$repID))
lmer_cell_spher_ABpl.p <- get_p_values_lmer(lmer_cell_spher_ABpl)
lmer_cell_spher_ABpl.p
if (lmer_cell_spher_ABpl.p["geom_info_ABpl$t","p-value"]<0.05){
  print("<RESULT> Sphericity of ABpl significantly changes over time")
} else {
  print("<RESULT> Sphericity of ABpl does NOT significantly change over time")
}
summary(lmer_cell_spher_ABpl)
