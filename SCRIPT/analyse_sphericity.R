print("SPHERICITY ANALYSIS#################################################################################################################")
#<<<MANUAL INPUT : REPLACE DATA_FOLDER TO LOCAL PATH !>>>>>
DATA_FOLDER <- "/home/wth/Downloads/testinfra/"
SCRIPT_FOLDER <- "/home/wth/Downloads/SYNC/"

source(paste0(SCRIPT_FOLDER,"fetch_command_line_arguments.R"))
source(paste0(SCRIPT_FOLDER,"helper_functions.R"))
source(paste0(SCRIPT_FOLDER,"set_input_parameters.R"))
output_folder <- paste0(output_folder,"/analyse_sphericity/")
source(paste0(SCRIPT_FOLDER,"load_R_packages_standard.R"))

# 0- load packages ---------
print(">>> 0 - Load packages <<< ")
packages <- c("lattice","MASS","pbkrtest","writexl")
if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages()))) 
}
library(lattice)
library(MASS)
library(pbkrtest)
library(lattice)
library(lmerTest)
library(writexl)

# 0 - read data / enrich / select -------------
print(">>> 0 - read data / enrich / select<<< ")
geom_info_full <-read_excel(input_file,sheet = "geom_info",na="NaN")
print('geom_info_full (raw) before enrichment -> ')
str(geom_info_full)

# join tc-data
tc_xml<-read_excel(input_file,sheet = "tc_xml",na="NaN")
geom_info_full <-merge(x = geom_info_full, y = tc_xml, by = c("tc_label"), all.x = TRUE)

#fix types
geom_info_full$res_t <- as.numeric(geom_info_full$res_t)
geom_info_full$repID <- as.character(geom_info_full$repID)
geom_info_full$nb_stack_EMS_division <- as.numeric(geom_info_full$nb_stack_EMS_division)


#enrich
#we define relative to EMS-division 
geom_info_full['t'] <- (geom_info_full$nb_stack_raw  - geom_info_full$nb_stack_EMS_division)  *  geom_info_full$res_t/60   #expressed in minutes,t=0 should be a reference time (e.g. EMS division)
if ("nb_pix_last_z" %in% colnames(geom_info_full)){} else {geom_info_full$nb_pix_last_z <- 0 }

# To get the time to the next division we need manual input 
      #future update : derive automatically instead of manual input
result <- tryCatch({
  #file_time_until_division <- '/home/wth/Downloads/SYNC/dds_7cell_geom.xlsx'
  file_time_until_division <- '/home/wth/Downloads/SYNC/dummy.xlsx'
  df_time_until_division <-read_excel(file_time_until_division ,sheet = "1",na="NaN")
  df_time_until_division$repID<-str_remove(df_time_until_division$replicate, "^0+") # some ad-hoc parsing to match repID
  names(df_time_until_division)[names(df_time_until_division) == "cellID"] <- "cell_name"
  df_time_until_division = subset(df_time_until_division,select = c(cell_name,repID,delta_t_division))
  df_time_until_division$delta_t_division <- as.numeric(df_time_until_division$delta_t_division)
  geom_info_full <-merge(x = geom_info_full, y = df_time_until_division, by = c("repID","cell_name"), all.x = TRUE)
  geom_info_full$delta_t_division <- as.numeric(geom_info_full$delta_t_division)
  geom_info_full$phase_division <- 'after division'
  geom_info_full$t_division <- (geom_info_full$delta_t_division * 3 ) + geom_info_full$t #3minute time intervals
  geom_info_full$phase_division[geom_info_full$t_division < 0] <- 'before division'
  return(geom_info_full)
}, error = function(err) {
  print('To get the time to the next division we need manual input. No file provided !')
  print(err)
})
if(!("t_division" %in% colnames(geom_info_full))){
  geom_info_full$phase_division <- 'No division info provided'
  geom_info_full$t_division <- 0
  }


write.csv(geom_info_full,paste0(output_folder,"geom_info_full_enriched.csv"), row.names = FALSE)


source(paste0(SCRIPT_FOLDER,"set_selection_filters.R"))

geom_info<- subset(geom_info_full,
                   geom_info_full$cell_name %in% l_cell_name  &
                     geom_info_full$embID %in% l_emb_ID       & 
                     geom_info_full$repID %in% l_rep_ID       &
                     geom_info_full$t >= as.integer(l_time_interval[[1]]) &  
                     geom_info_full$t <= as.integer(l_time_interval[[2]])
                    )
# geom_info<- subset(geom_info_full,geom_info_full$cell_name %in% l_cell_name &
#                      geom_info_full$embID %in% l_emb_ID)

attach(geom_info)
str(geom_info)

# set process flags
if (length(unique(geom_info['t'])[[1]]) > 1){
  flag_single_timepoint = FALSE
} else {
  flag_single_timepoint = TRUE
}

if (length(unique(geom_info['category'])[[1]]) > 1){
  flag_single_category = FALSE
} else {
  flag_single_category = TRUE
}

## 1 -  exploratory analysis (fig S3)--------
print('## 1 - exploratory analysis (fig S3)----')
summary(geom_info)

#   OLS for every repID/cell_name combination (sphericity)
sphericityplot<- ggplot(geom_info, aes(x=t, y=sphericity)) +geom_point(alpha=0.5) + 
  ylim(min(sphericity,na.rm = T), max(sphericity,na.rm = T))

sphericityplot<-sphericityplot +  facet_grid(repID~cell_name) + stat_smooth(method="lm")
sphericityplot <- sphericityplot  + theme(axis.title.x=element_text(size=16),
                                          axis.title.y=element_text(size=16),
                                          legend.title = element_blank())
sphericityplot <- sphericityplot+ scale_x_continuous(labels = scales::number_format(accuracy = 1)) #or 0.01
#sphericityplot<- sphericityplot + scale_x_continuous(name="t", breaks=seq(-10,10,10),limits=c(-15, 15))
sphericityplot
ggsave(paste0(output_folder,"_sphericity_OLS_facet_grid(repID-cell_name).jpg"),dpi='retina')

#   OLS for every repID/cell_name combination (bounding box volume)
geom_info$bounding_box_vol <- geom_info$extents_length * geom_info$extents_width * geom_info$extents_height
bounding_box_plot<- ggplot(geom_info, aes(x=t, y=bounding_box_vol)) +geom_point() + 
  ylim(min(geom_info$bounding_box_vol,na.rm = T), max(geom_info$bounding_box_vol,na.rm = T))
bounding_box_plot<-bounding_box_plot +  facet_grid(repID~cell_name) #+ stat_smooth(method="lm")
bounding_box_plot <-bounding_box_plot + scale_x_continuous(labels = scales::number_format(accuracy = 1)) #or 0.01
bounding_box_plot
#ggsave(paste0(output_folder,"<bounding box>OLS_facet_grid(repID-cell_name).jpg"))

#   fig S3 - sphericityplot OLS (fig5c) for every repID/cell_name combination (sphericity vs time to cell division)----
sphericityplot<- ggplot(geom_info, aes(x=t, y=sphericity,colour=phase_division)) +geom_point(alpha=0.5) + 
  ylim(min(sphericity,na.rm = T), max(sphericity,na.rm = T))
sphericityplot<- sphericityplot + scale_x_continuous(name="\u0394t cell division (min)") #, breaks=seq(-10,10,10),limits=c(-15, 15))
sphericityplot<-sphericityplot +  facet_grid(repID~cell_name) #+ stat_smooth(method="lm",se=FALSE)
sphericityplot <-sphericityplot + labs(y="sphericity",title="")
sphericityplot <- sphericityplot  + theme(axis.title.x=element_text(size=16),
                                          axis.title.y=element_text(size=16),
                                          legend.title = element_blank(),
                                          legend.position="bottom",legend.direction = "horizontal")
sphericityplot
ggsave(paste0(output_folder,"_sphericity_OLS_facet_grid[repID-cell_name]_vs_delta_t_division.jpg"),dpi='retina')

#aggregated over replicates
sphericityplot<- ggplot(geom_info, aes(x=t, y=sphericity,colour=repID)) +geom_point(alpha=0.5) + 
  ylim(min(sphericity,na.rm = T), max(sphericity,na.rm = T))
sphericityplot<- sphericityplot + scale_x_continuous(name="\u0394t cell division EMS (min)")#, breaks=seq(-10,10,10),limits=c(-15, 15))
sphericityplot<-sphericityplot +  facet_wrap(~cell_name) #+ stat_smooth(method="lm",se=FALSE)
sphericityplot <- sphericityplot + geom_smooth(method = "loess",se = FALSE)
sphericityplot
ggsave(paste0(output_folder,"_sphericity_OLS_facet_wrap_cell_name.jpg"),dpi='retina')

#aggregated over category
if (!flag_single_category){
  sphericityplot<- ggplot(geom_info, aes(x=t, y=sphericity,colour=category)) +geom_point(alpha=0.5) + 
    ylim(min(sphericity,na.rm = T), max(sphericity,na.rm = T))
  sphericityplot<- sphericityplot + scale_x_continuous(name="\u0394t cell division EMS (min)")
  sphericityplot<-sphericityplot +  facet_wrap(~cell_name)# + stat_smooth(method="lm",se=FALSE)
  sphericityplot <- sphericityplot + geom_smooth(method = "loess")
  #sphericityplot <- sphericityplot + geom_jitter()
  sphericityplot
  ggsave(paste0(output_folder,"_sphericity_OLS_facet_wrap_cell_name-colour-category.jpg"),dpi='retina')
}

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
ggsave(paste0(output_folder,"_sphericity_delta_t_division_turning point.pdf"),device = "pdf")

sphericityplot_stripped <- strip_plot(sphericityplot)
sphericityplot_stripped
#ggsave(paste0(output_folder,"_sphericity_delta_t_division_turning point(stripped).pdf"),device = "pdf")


#2.2 fitting ----

# E  sphericity significantly different between categories ?


if (!flag_single_category){
  cell_name = 'E'
  cond1 = geom_info$cell_name == cell_name
  cond2 = geom_info$t < 10  #!!!!!  cut out time when E starts rounding again (building up to division)
  sphericity_test <- geom_info[cond1 & cond2,]
  
  source('sphericity_test_category.R')

}

# MS  sphericity significantly different between categories ?
if (!flag_single_category){
  cell_name = 'MS'
  cond1 = geom_info$cell_name == cell_name
  cond2 = geom_info$t < 10  #!!!!!  cut out time when E starts rounding again (building up to division)
  sphericity_test <- geom_info[cond1 & cond2,]
  
  source('sphericity_test_category.R')

}


# ABpl sphericity significantly different between categories ?
if (!flag_single_category){
  cell_name = 'ABpl'
  cond1 = geom_info$cell_name == cell_name
  cond2 = geom_info$t < 10  #!!!!!  cut out time when E starts rounding again (building up to division)
  sphericity_test <- geom_info[cond1 & cond2,]
  
  source('sphericity_test_category.R')
  
}


