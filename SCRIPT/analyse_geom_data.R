print("VOLUME ANALYSIS#################################################################################################################")


#<<<MANUAL INPUT : REPLACE SCRIPT AND DATA FOLDER TO LOCAL PATH !>>>>>
SCRIPT_FOLDER <- "/home/wth/Downloads/SYNC/"
DATA_FOLDER <- "/home/wth/Downloads/testinfra/"

source(paste0(SCRIPT_FOLDER,"fetch_command_line_arguments.R"))
source(paste0(SCRIPT_FOLDER,"helper_functions.R"))
source(paste0(SCRIPT_FOLDER,"set_input_parameters.R"))
output_folder <- paste0(output_folder,"/analyse_geom_data/")
source(paste0(SCRIPT_FOLDER,"load_R_packages_standard.R"))


# 0- load packages ---------
print(">>> 0 - Load extra packages <<< ")
require(remotes)
packages <- c("nlme","lme4","lattice","devtools")
if (length(setdiff(packages, rownames(installed.packages())))> 0) {
  install.packages(setdiff(packages, rownames(installed.packages())),repos = c("https://CRAN.R-project.org","https://cran.r-project.org") )
}
library(lme4)
library(lattice)
library(magrittr)
library(Matrix)
library(MASS)
library(lattice)
library(nlme)
library(lmerTest)
library(pbkrtest)


# 0 - read data / enrich / select -------------
print(">>> 0 - read data / enrich / select<<< ")
geom_info_full <-read_excel(input_file,sheet = "geom_info",na="NaN")
print('geom_info_full (raw) before enrichment -> ')
str(geom_info_full)

# join tc-data
tc_xml<-read_excel(input_file,sheet = "tc_xml",na="NaN")
geom_info_full <-merge(x = geom_info_full, y = tc_xml, by = c("tc_label"), all.x = TRUE)

#join parentids
d_parentid_cellname<-read_excel(input_file,sheet = "d_parentid_cellname",na="NaN")
if (sum(is.na(d_parentid_cellname$cell_name))>0) {
  print(paste('<DATA CORRUPT> : d_parentid_cellname has ',sum(!is.na(d_parentid_cellname$cell_name)), 'NA values. Merge will fail. Check the file !'))
  print(str(d_parentid_cellname))
} else{
  geom_info_full <-merge(x = geom_info_full, y = d_parentid_cellname, by = c("tc_label","cell_name"), all.x = TRUE)
  }


# join nb pixels last slice (if available)
result <- tryCatch({
  last_slice_info<-read_excel(input_file,sheet = "last_slice_info",na="NaN")
}, error = function(err) {
  print('No last info available. Default 0 is presumed')
})
if (exists("last_slice_info")){
  if (sum(is.na(last_slice_info$parentid))>0) {
    print(paste('<DATA CORRUPT> : last_slice_info has ',sum(!is.na(last_slice_info$parentid)), 'NA values. Merge will fail. Check the file !'))
    print(str(last_slice_info))
  } else{
    geom_info_full <-merge(x = geom_info_full, y = last_slice_info, by = c("tc_label","parentid"), all.x = TRUE)
    geom_info_full['nb_pix_last_z'][is.na(geom_info_full['nb_pix_last_z'])] <- 0
  }
}

#fix types
geom_info_full$res_t <- as.numeric(geom_info_full$res_t)
geom_info_full$repID <- as.character(geom_info_full$repID)
geom_info_full$nb_stack_EMS_division <- as.numeric(geom_info_full$nb_stack_EMS_division)

#enrich
geom_info_full['t'] <- (geom_info_full$nb_stack_raw - geom_info_full$nb_stack_EMS_division)  *  geom_info_full$res_t/60   #expressed in minutes,t=0 should be a reference time (e.g. EMS division)
geom_info_full$t <- as.numeric(geom_info_full$t)
if ("nb_pix_last_z" %in% colnames(geom_info_full)){} else {geom_info_full$nb_pix_last_z <- 0 }

write.csv(geom_info_full,paste0(output_folder,"geom_info_full_enriched.csv"), row.names = FALSE)

source(paste0(SCRIPT_FOLDER,"set_selection_filters.R"))

geom_info<- subset(geom_info_full,
                   geom_info_full$cell_name %in% l_cell_name &
                     geom_info_full$embID %in% l_emb_ID      & 
                    # geom_info_full$repID %in% l_rep_ID      &
                     geom_info_full$t >= as.integer(l_time_interval[[1]]) &  
                     geom_info_full$t <= as.integer(l_time_interval[[2]])
                   )


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

attach(geom_info)
str(geom_info)

print(">>> 1. exploratory analysis <<<")
# 1 exploratory analysis----------

## Spaghettiplot (mean volume per cell_name over time) ####
#print (geom_info_full)
if (!flag_single_timepoint){
par(mfrow=c(1,1))
interaction.plot(geom_info$t,geom_info$cell_name,geom_info$volume, xlab="time(min)", ylab="volume (mu_m^3)", legend=T)
mean_volumes<-aggregate(geom_info$volume, by=list(geom_info$t), FUN=mean)
mean_volumes$Group.1 <-mean_volumes$Group.1
lines(mean_volumes,col='red',lwd=2) 
#ggsave(paste0(output_folder,"_volume_spaghetti_plot.jpg"))
}

print('-> ordinary least squares')
## ordinary least squares -------
# OLS for every cell_name
#   ignoring batch effect (aggrating over cell_name)
volumeplot<- ggplot(geom_info, aes(x=t, y=volume, factor(cell_name), color = repID)) +geom_point(alpha=0.5) + ylim(0, NA)
volumeplot2<-volumeplot +  facet_wrap(~cell_name,ncol=7) #+ stat_smooth(method="lm", se=FALSE)
#volumeplot2 <- volumeplot2 + theme(axis.text.x = element_text(face="bold", color="#993333", size=14, angle=45))
volumeplot2 <- volumeplot2 + scale_x_continuous(name="\u0394t cell division EMS (min)",labels = scales::number_format(accuracy = 1)) #or 0.01
volumeplot2
ggsave(paste0(output_folder,"_volume_OLS_facetwrap(cell_name).jpg"))

#   ignoring cell_name effect (aggrating over cell_name)
volumeplot<- ggplot(geom_info, aes(x=t, y=volume, factor(repID), color = cell_name)) +geom_point(alpha=0.5) + ylim(0, NA)
volumeplot3<-volumeplot +  facet_wrap(~repID,ncol=5) + stat_smooth(method="lm",se=FALSE)
volumeplot3 <- volumeplot3 + scale_x_continuous(name="\u0394t cell division EMS (min)",labels = scales::number_format(accuracy = 1)) #or 0.01
volumeplot3
ggsave(paste0(output_folder,"_volume_OLS_facetwrap(repID).jpg"))

#   OLS embryo size over time per replicate
#   (rem)
df_emb_vol <- aggregate(geom_info_full$volume, by=list(geom_info_full$t,geom_info_full$repID), FUN=sum)
names(df_emb_vol)<-cbind('t','repID','emb_vol')
str(df_emb_vol)
volumeplot<- ggplot(df_emb_vol, aes(x=t, y=emb_vol, color=repID)) +geom_point() + ylim(min(df_emb_vol$emb_vol,na.rm = T), max(df_emb_vol$emb_vol,na.rm = T))
volumeplot4<-volumeplot  + stat_smooth(method="lm",se=FALSE)
volumeplot4 <- volumeplot4 +  ggtitle("total embryo volume over time")
volumeplot4
ggsave(paste0(output_folder,"_volume_embryosize total(repID).jpg"))
#l_emb_vol <-by(geom_info$vol,geom_info$t,sum,na.rm=TRUE)

#   OLS for every repID/cell_name combination (volume)
volumeplot<- ggplot(geom_info, aes(x=t, y=volume)) +geom_point() + 
  ylim(min(volume,na.rm = T), max(volume,na.rm = T))
volumeplot5<-volumeplot +  facet_grid(repID~cell_name) + stat_smooth(method="lm")
volumeplot5 <- volumeplot5 + scale_x_continuous(name="\u0394t cell division EMS (min)",labels = scales::number_format(accuracy = 1)) #or 0.01
volumeplot5
ggsave(paste0(output_folder,"_volume_OLS_facet_grid(repID-cellname).jpg"))

#   OLS for every category/cell_name combination (volume)
if (!flag_single_category){
  volumeplot<- ggplot(geom_info, aes(x=t, y=volume, color=repID)) +geom_point(alpha=0.5) + 
    ylim(min(volume,na.rm = T), max(volume,na.rm = T))
  volumeplot6<-volumeplot +  facet_grid(category~cell_name) + stat_smooth(method="lm",se=FALSE)
  volumeplot6 <- volumeplot6 + scale_x_continuous(name="\u0394t cell division EMS (min)",labels = scales::number_format(accuracy = 1)) #or 0.01
  volumeplot6
  ggsave(paste0(output_folder,"_volume_OLS_facet_grid(category-cellname).jpg"))
  
  volumeplot<- ggplot(geom_info, aes(x=t, y=volume,color=category)) +geom_point(alpha=0.5) + 
    ylim(min(volume,na.rm = T), max(volume,na.rm = T))
  volumeplot7<-volumeplot +  facet_wrap(~cell_name) + stat_smooth(method="lm",se=FALSE)
  volumeplot7 <- volumeplot7 + scale_x_continuous(name="\u0394t cell division EMS (min)",labels = scales::number_format(accuracy = 1)) #or 0.01
  volumeplot7
  ggsave(paste0(output_folder,"_volume_OLS_facet_wrap(-cellname-color=category).jpg"))
  
}



print('>>> 2 Fitting the fixed effects model (initial fit, before simplifications) <<<')
#2  Fitting the data  ####### ###################################################################
if (flag_single_category){
  lmer_cell_vol <- lmer(volume~ 1 + cell_name + t + nb_pix_last_z + (1 | repID),data=geom_info)
} else {
  lmer_cell_vol <- lmer(volume~ 1 + cell_name + category + t + nb_pix_last_z + (1 | repID),data=geom_info)
}
lmer_cell_vol.p.KR <- get_p_values_lmer(lmer_cell_vol)
lmer_cell_vol.p.KR

print('---> Checking assumptions of model')
#check assumptions ---------------

#CHECK is there correlation between residuals ? (ie is it a time-series)
# Rbook crawley p645 "Time Series Analysis in Mixed-Effects Models"
# different packages (not so obvious in lme4): "https://bbolker.github.io/mixedmodels-misc/notes/corr_braindump.html"
#lme produces exactly the same ouput as lmer, but offers time correlation modeling
result <- tryCatch({
  print(paste('packageVersion(nlme): ',packageVersion("nlme")))
  lme_cell_vol<-lme(volume~1 + cell_name +t    ,data=geom_info,random=~ 1  | repID)
  summary(lme_cell_vol)
  ACFplot <- plot(ACF(lme_cell_vol,maxLag = 4),alpha=0.05) # no time correlation ACFplot
  ACFplot
  #ggsave(paste0(output_folder,"<check>ACF plot time correlation of residuals).jpg"))
  plot(lme_cell_vol,resid(.,type="p")~fitted(.)|repID)
  #ggsave(paste0(output_folder,"<check>residual plot.jpg"))
  qqnorm(lme_cell_vol,~resid(.)|repID)  
  #ggsave(paste0(output_folder,"<check>qqplot.jpg"))
}, error = function(err) {
  print('Correlation between residuals could not be determined')
})


print('-----> Hypothesis testing with model :  Group E and MS ?')
#group E and MS ?--------
#   recommend method for LMMs = conditional F-tests with df correction cfr Bolker

geom_info$cell_name_group1 <- geom_info$cell_name
geom_info$cell_name_group1 [geom_info$cell_name_group1  %in% c('E','MS')] <- 'E/MS'


lmer_cell_vol <- lmer(volume~ 1 + cell_name + t + nb_pix_last_z + (1 | repID) ,
                      data = geom_info,
                      REML = FALSE)
lmer_cell_vol_EMS <- lmer(volume~ 1 + cell_name_group1 + t + nb_pix_last_z + (1 | repID) ,
                          data = geom_info,
                          REML = FALSE)
summary(lmer_cell_vol)
summary(lmer_cell_vol_EMS)

#likelihood ratio test: is general model significantly better than EMS model ?
# anova(null.model, model1, test = “Chisquare”) # null is simple, constraint model
anov_EMS_test = anova(lmer_cell_vol_EMS, lmer_cell_vol,test=TRUE)
anov_EMS_test
if (anov_EMS_test$`Pr(>Chisq)`[2]<0.05){
  print("<RESULT> E is significantly different in volume from MS (null model with E/MS combined is rejected)")
} else {
  print("<RESULT> E is NOT significantly different in volume from MS (null model with E/MS combined can not be rejected)")
}

print('----> Hypothesis testing with model :  Group all AB cells ?')
#groups all AB cells ? ----

geom_info$cell_name_group2  <- geom_info$cell_name
geom_info$cell_name_group2 [geom_info$cell_name_group2  %in% c('ABal','ABar','ABpl','ABpr')] <- 'AB'
lmer_cell_vol_AB <- lmer(volume~ 1 + cell_name_group2+ t + nb_pix_last_z  + (1 | repID),data=geom_info )
summary(lmer_cell_vol_AB)

#likelihood ratio test
anov_AB_test = anova(lmer_cell_vol_AB,lmer_cell_vol,test=TRUE)
anov_AB_test
if (anov_AB_test$`Pr(>Chisq)`[2]<0.05){
  print("<RESULT> AB cells are significantly different in volume (null model with all AB cells combined is rejected)")
} else {
  print("<RESULT> AB cells are NOT significantly different in volume (null model with all AB cells combined can not be rejected)")
}

print('----> #groups [ABar, ABal,ABpl] and [ABpr] cells ?')
#groups [ABar, ABal,ABpl] and [ABpr] cells -------
geom_info$cell_name_group3  <- geom_info$cell_name
geom_info$cell_name_group3 [geom_info$cell_name_group3  %in% c('ABal','ABar','ABpl')] <- 'AB3'

attach(geom_info)
lmer_cell_vol_AB3 <- lmer(volume~ 1 + cell_name_group3+  t + nb_pix_last_z + (1 | repID) ,data=geom_info)
summary(lmer_cell_vol_AB3)
lmer_cell_vol_AB3.p <- get_p_values_lmer(lmer_cell_vol_AB3)
lmer_cell_vol_AB3.p
#likelihood ratio test:
anov_AB3_test = anova(lmer_cell_vol_AB3,lmer_cell_vol,test=TRUE)
anov_AB3_test
if (anov_AB3_test$`Pr(>Chisq)`[2]<0.05){
  print("<RESULT> (ABal;ABar;ABpl) cells are significantly different in volume (null model with (ABal;ABar;ABpl) cells combined is rejected)")
} else {
  print("<RESULT> (ABal;ABar;ABpl)  cells are NOT significantly different in volume (null model with (ABal;ABar;ABpl) cells combined can not be rejected)")
}

print('---> Hypothesis testing with model :  is the random effect (=embryo size) significant ?')
#random effect (embryo size) ---------

lmer_cell_vol_norandeff <- glm(volume~ 1 + cell_name + t,data=geom_info)
summary(lmer_cell_vol_norandeff)

# variance can not be < 0 so we can not use anova
logratio_statistic <- -2 * (logLik(lmer_cell_vol_norandeff,REML=TRUE) - logLik(lmer_cell_vol,REML=TRUE))
P_chi_df1 <-as.numeric(pchisq(logratio_statistic, df=1,lower.tail = FALSE))
P_chi_df2 <-as.numeric(pchisq(logratio_statistic, df=2,lower.tail = FALSE))
P_null_model <- 1/2 * (P_chi_df1 + P_chi_df2)
P_null_model
if (P_null_model<0.05){
  print("<RESULT> The random effect (=embryo size) is sigificant")
} else {
  print("<RESULT> The random effect (=embryo size) is NOT sigificant")
}


if (!flag_single_category){
  print('---> Hypothesis testing with model :  is the category effect significant ?')
  lmer_cell_nocategory <- lmer(volume~ 1 + cell_name +  t + nb_pix_last_z + (1 | repID) ,data=geom_info)
  summary(lmer_cell_nocategory)

  # variance can not be < 0 so we can not use anova
  logratio_statistic <- -2 * (logLik(lmer_cell_nocategory,REML=TRUE) - logLik(lmer_cell_vol,REML=TRUE))
  P_chi_df1 <-as.numeric(pchisq(logratio_statistic, df=1,lower.tail = FALSE))
  P_chi_df2 <-as.numeric(pchisq(logratio_statistic, df=2,lower.tail = FALSE))
  P_null_model <- 1/2 * (P_chi_df1 + P_chi_df2)
  P_null_model
  if (P_null_model<0.05){
    print(paste("<RESULT> The category effect is sigificant (p-null-model=", P_null_model,")"))
  } else {
    print(paste("<RESULT> The category effect is NOT sigificant (p-null-model=", P_null_model,")"))
  }
}


print('<<< 3- Fitting the FINAL MODEL >>>')
#3 Final model (volume) ###############################

if (all(geom_info$nb_pix_last_z == 0)) {
  if (flag_single_timepoint) {
    cell_vol.lmer.final <-
      lmer(volume ~ 1 + cell_name  + (1 | repID), data = geom_info)
  } else {
    cell_vol.lmer.final <-
      lmer(volume ~ 1 + cell_name  + t + (1 | repID), data = geom_info)
  }
}  else {
  if (flag_single_timepoint) {
    cell_vol.lmer.final <-
      lmer(volume ~ 1 + cell_name + nb_pix_last_z + (1 |
                                                       repID), data = geom_info)
  } else {
    cell_vol.lmer.final <-
      lmer(volume ~ 1 + cell_name + nb_pix_last_z + t + (1 |
                                                           repID), data = geom_info)
  }
}
print("--> final model (volume) : ")
summary(cell_vol.lmer.final)

print('--> add corrected volume based on fit')
if (all(geom_info$nb_pix_last_z==0)){
geom_info$vol_corr<-geom_info$volume
}else{
  geom_info$vol_corr<-geom_info$volume - (coef(summary(cell_vol.lmer.final))['nb_pix_last_z','Estimate'] * geom_info$nb_pix_last_z)
}

print('--> visual check of residuals')
plot(cell_vol.lmer.final) 
qqnorm(resid(cell_vol.lmer.final))
qqline(resid(cell_vol.lmer.final)) 

print('-->  CI (fixed and random effects) viz (fig5a)')
newdat2 <- expand.grid(
  cell_name=l_cell_name,
  nb_pix_last_z=0,
  t=1,
  volume = 0
)
newdat2$volume<-predict(cell_vol.lmer.final,newdat2,re.form=NA)
mm2 <- model.matrix(object=terms(cell_vol.lmer.final),data=newdat2)
pvar2 <- diag(mm2 %*% tcrossprod(vcov(cell_vol.lmer.final),mm2)) # tcrossproduct = vcov(cell_vol.lmer.final) %*% t(mm2).
# identical to : pvar2 <- diag(mm2 %*% vcov(cell_vol.lmer.final) %*% t(mm2))
#pvar2 = variance of fixed effects 
#   the model gives the std. dev, BUT these estimates are defined versus a reference level (ABal),
#   so the variance of estimates off the individual cell must be calculated by incorporating the covariance
#   of the cell with the reference and also the covariance of all other cells on the ref. level
#   this is the variance of FIXED effects only
tvar2 <- pvar2+VarCorr(cell_vol.lmer.final)$repID[1]  # adds estimated variance of random intercept effect of repID (=std Dev ** 2)
cmult <- 1.96 
newdat2 <- data.frame(
  newdat2
  , plo = newdat2$vol-cmult*sqrt(pvar2)  #low limit 0.95 CI of fixed effect
  , phi = newdat2$vol+cmult*sqrt(pvar2)  #high limit 0.95 CI of fixed effect
  , tlo = newdat2$vol-cmult*sqrt(tvar2)  #low limit 0.95 CI of fixed effect + random effect
  , thi = newdat2$vol+cmult*sqrt(tvar2) #high limit 0.95 CI of fixed effect + random effect
)

g0 <- ggplot(newdat2, aes(x=cell_name, y=volume))+geom_point()
g0 <- g0 + geom_pointrange(aes(ymin = plo, ymax = phi))+
  labs(title="CI based on fixed-effects uncertainty ONLY")
g0 <- g0  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"))
g0
#ggsave(paste0(output_folder,'_volume_CI based on fixed-effects uncertainty(CI).pdf'),device = "pdf")
write.csv(newdat2,paste0(output_folder,"final_model_volume_estimations_newdat2.csv"), row.names = FALSE)

g0 <- ggplot(newdat2, aes(x=cell_name, y=volume))+geom_point()
g0 + geom_pointrange(aes(ymin = tlo, ymax =thi))+
  labs(title="CI based on fixed-effect + random effects uncertainty")
g0 <- g0 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
#ggsave(paste0(output_folder,'_volume_CI based on fixed-effect + random effects uncertainty.pdf'),device = "pdf")



print('---> add estimates per repID')
newdat3 <- expand.grid(
  repID = l_rep_ID,
  cell_name=l_cell_name,
  nb_pix_last_z=0,
  t=1,
  volume = 0
)
newdat3$repID <- as.character(newdat3$repID)

newdat3$volume<-predict(cell_vol.lmer.final,newdat3)
g0 <- ggplot(newdat2, aes(x=cell_name, y=volume))+geom_point()
g0 <- g0 + geom_pointrange(aes(ymin = plo, ymax = phi),size=1)
#g0 <- g0+ labs(title="CI based on fixed-effects uncertainty ONLY")
g0 <- g0  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                  axis.title.x=element_blank(),axis.title.y=element_blank(),
                  axis.text.x= element_text(face="bold", size=16),
                  axis.text.y= element_text(face="bold", size=16))
# g0 <-t heme(axis.title.x=element_blank(),

g0
g1 <- g0 + geom_line(data=newdat3, 
                     aes(x=cell_name, y=volume,group=repID,colour=repID))
g1  <- g1 + labs(title="model prediction cell volumes",hjust=0.5) +
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()
g1 
ggsave(paste0(output_folder,'_volume_model prediction cell volumes.pdf'),device = "pdf")


print('---> overlay with real data (averaged over time)')
geom_info_agg <- geom_info %>%
  group_by(category,repID,cell_name) %>%
  summarize(volume = mean(volume, na.rm=TRUE),
            vol_corr = mean(vol_corr, na.rm=TRUE))

geom_info_agg_complete <- geom_info %>%
  group_by(cell_name) %>%
  summarize(volume = mean(volume, na.rm=TRUE),
            vol_corr = mean(vol_corr, na.rm=TRUE))
print('Aggregate volumes over all replicates and all timepoints:->');print(geom_info_agg_complete)


g2 <- g0 + geom_line(data=geom_info_agg, 
                     aes(x=cell_name, y=volume,group=repID,colour=repID),
                     linetype = "dashed")
g2  <- g2 + labs(title="volume measurements (time averaged, uncorrected)",hjust=0.5) +
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()
g2 <- g2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
g2
#ggsave(paste0(output_folder,'_volume_measurements (time averaged, uncorrected).pdf'),device = "pdf")


g3 <- ggplot(newdat3, aes(x=cell_name, y=volume,group=repID,colour=repID))+geom_line()
g3 <- g3 + geom_line(data=geom_info_agg, 
                     aes(x=cell_name, y=volume,group=repID,colour=repID),
                     linetype = "dashed")
g3  <- g3 + labs(title="volume measurements and model predictions",hjust=0.5) +
  theme(plot.title = element_text(hjust = 0.5))
g3 <- g3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
g3
ggsave(paste0(output_folder,'_volume_measurements and model predictions.pdf'),device = "pdf")

#fig5a-------
print('---> fig5a- _volume_measurements (time averaged, corrected)-')
if (flag_single_category){
  g4 <- g0 + geom_line(data=geom_info_agg, size=1,
                         aes(x=cell_name, y=vol_corr,group=repID,colour=repID),
                               linetype = "dashed")
} else {
  g4 <- g0 + geom_line(data=geom_info_agg, size=1,
                       aes(x=cell_name, y=vol_corr,group=repID,colour=category),
                       linetype = "dashed")
  
  }



# g4  <- g4 + labs(title="volume measurements (time averaged, corrected)",hjust=0.5) +
#   theme(plot.title = element_text(hjust = 0.5))
g4 <- g4 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
g4 <- g4 + theme(legend.title = element_blank(),
                 legend.position="bottom",legend.direction = "horizontal",
                 legend.text = element_text(size = "large"),
                 axis.text.x = element_text(angle=90,size=12))
g4


write.csv(geom_info_agg,paste0(output_folder,"g4_geom_info_agg.csv"), row.names = FALSE)
ggsave(paste0(output_folder,'_volume_measurements (time averaged, corrected).pdf'),device = "pdf")

g4_stripped_full <- strip_plot(g4,strip='full')
g4_stripped_full
#ggsave(paste0(output_folder,'_volume_measurements (time averaged, corrected)(stripped_full).pdf'),device = "pdf")

g4_stripped <- strip_plot(g4)
g4_stripped
#ggsave(paste0(output_folder,'_volume_measurements (time averaged, corrected)(stripped).pdf'),device = "pdf")

#Final
#Corrected volume viz  (fig S1)------------
print('---> Corrected volume viz  (fig S1 : facet grids) (final results)------------')

#   OLS for every repID/cell_name combination (volume corr + vol)
df_melt_vol <- melt(geom_info,
                    id.vars=c("repID","cell_name","nb_stack_mpacts","t","category"),
                    measure.vars = c("volume","vol_corr"))
volumeplot<- ggplot(df_melt_vol, aes(x=t, y=value,colour=variable)) +geom_point() + 
  ylim(min(df_melt_vol$value,na.rm = T), max(df_melt_vol$value,na.rm = T))
volumeplot5<-volumeplot +  facet_grid(repID~cell_name) + stat_smooth(method="lm")
volumeplot5 <- volumeplot5 + scale_x_continuous(name="\u0394t cell division EMS (min)",labels = scales::number_format(accuracy = 1)) #or 0.01
volumeplot5
ggsave(paste0(output_folder,"_volumeCorr+volume_OLS_facet_grid(repID-cell_name).jpg"))

#   OLS for every repID/cell_name combination (only volume corr) = fig S1
cond <-df_melt_vol['variable'] == "vol_corr"
volumeplot<- ggplot(df_melt_vol[cond,], aes(x=t, y=value)) +
  geom_point() + 
  ylim(min(df_melt_vol$value,na.rm = T), max(df_melt_vol$value,na.rm = T))
volumeplot6<-volumeplot +  facet_grid(repID~cell_name) + stat_smooth(method="lm",se=FALSE)
volumeplot6<-volumeplot6 + labs(y="volume (\u03bcm³)", x = "\u0394t cell division EMS (min)",title="")
volumeplot6<-volumeplot6 + scale_x_continuous(breaks=c(0,5,10))
volumeplot6<-volumeplot6 + theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16))
volumeplot6
ggsave(paste0(output_folder,"_volumeCorr_OLS_facet_grid(repID-cell_name)_S1.jpg"),dpi="retina")
#ggsave(paste0(output_folder,"<volumeCorr>OLS_facet_grid(repID-cell_name)_S1.pdf"),dpi="retina")


## ordinary least squares -------
# OLS for every cell_name
volumeplot<- ggplot(df_melt_vol[cond,], aes(x=t, y=value,colour=repID)) +geom_point(alpha=0.5) + 
  ylim(min(df_melt_vol$value,na.rm = T), max(df_melt_vol$value,na.rm = T))
volumeplot5b <-volumeplot +  facet_wrap(~cell_name,ncol=7) 
volumeplot5b<-volumeplot5b + labs(y="volume_corr (\u03bcm³)", x = "\u0394t cell division EMS (min)",title="")
volumeplot5b<-volumeplot5b + scale_x_continuous(breaks=c(0,5,10))
volumeplot5b <- volumeplot5b + geom_smooth(method = "loess",se = FALSE)
volumeplot5b <-volumeplot5b + theme(axis.title.x=element_text(size=16),axis.title.y=element_text(size=16))

volumeplot5b
ggsave(paste0(output_folder,"_volumeCorr+volume_OLS_facet_wrap(-cell_name).jpg"))

volumeplot<- ggplot(df_melt_vol[cond,], aes(x=t, y=value,color=category)) +geom_point(alpha=0.5) + 
  ylim(min(volume,na.rm = T), max(volume,na.rm = T))
volumeplot7b<-volumeplot +  facet_wrap(~cell_name) + stat_smooth(method="lm",se=FALSE)
volumeplot7b<-volumeplot7b + labs(y="volume_corr (\u03bcm³)", x = "\u0394t cell division EMS (min)",title="")
#volumeplot7b <- volumeplot7b + scale_x_continuous(name="\u0394t cell division EMS (min)",labels = scales::number_format(accuracy = 1)) #or 0.01
volumeplot7b
ggsave(paste0(output_folder,"_volumeCorr_OLS_facet_wrap(-cellname-color=category).jpg"))




## Viz volumeratios (fig 5b)--------------
print('-->  Viz volumeratios (fig 5b) (ratio plots) --------------')
#  with ratios, the random effect of embryosize is neutralized, however
# correlation over time is not taken into account (not an hierarchical model)
#  vol_corr is picked so corrected for z-positioning
#  time aspect must be modelled

#functions
get_vol_ratios <- function(geom_info=df_ipt, cell1="E", cell2='MS') {
  n_ratio <- numeric(length = 0)
  l_rep <- c()
  l_t <- numeric(length = 0)
  l_z_orientation <- c()
  df_iter = expand.grid(unique(geom_info$repID),unique(geom_info$t))
  for (row in 1:nrow(df_iter)){
    rep_i = df_iter[row,1]
    t_i = df_iter[row,2]
    z_orientation_i = geom_info$z_orientation[geom_info$t==t_i & geom_info$repID==rep_i & geom_info$cell_name==cell1] 

    ratio <- geom_info$vol_corr[geom_info$t==t_i & geom_info$repID==rep_i & geom_info$cell_name==cell1] /
      geom_info$vol_corr[geom_info$t==t_i & geom_info$repID==rep_i & geom_info$cell_name==cell2] 
    
    if (length(ratio)>1){
      print(paste0('<DATA INCONSISTENCY>!!There is a duplicate cell pair in the data ! this is not possible.  Check this combination =>','cell1=',cell1,' cell2=',cell2,' time(minutes)=',t_i, ' replicate=',rep_i))
      print('    -->Probably caused by 2 parentids linked to the same cell name (cfr d_parentid_cellname).  You have to correct the lineage file in that case')
      ratio <- ratio[1] #only use the first
    }
    if (length(ratio)>0 & is.numeric(ratio)) {
      n_ratio <- append(n_ratio, ratio )
      l_rep <- append(l_rep,rep_i)
      l_t <- append(l_t,t_i)
      l_z_orientation <- append(l_z_orientation,z_orientation_i)
    }
    
  }
  
  df_ratio= data.frame(repID=l_rep,t=l_t,ratioID=paste0(cell1,"/",cell2),ratio=n_ratio,z_orientation=l_z_orientation)
  
}


l_cell_pairs = list(c('ABpr','ABpl'),c('E','MS'),c('P3','C')) # INPUT VALUES HERE
df_ratios = {}
df_ci = {}
l_category<- unique(geom_info_full['category'])[[1]]
for (category_i in l_category) {
  print(category_i)
  geom_info_category <- geom_info[geom_info['category']==category_i,]
  for (cell_pair_i in l_cell_pairs) {
    print(cell_pair_i)
    df_ratio = get_vol_ratios(geom_info_category,cell1 = cell_pair_i[1],cell2 = cell_pair_i[2])
    df_ratio$category <- category_i
    df_ratios = rbind(df_ratios,df_ratio)
    df_ratio <- df_ratio[!is.na(df_ratio$ratio),]
    n_ci = CI(df_ratio$ratio, ci=0.95)
    df_ci<-rbind(df_ci,
                 data.frame(category=category_i,ratioID=df_ratio$ratioID[1],
                            upper=n_ci['upper'],mean=n_ci['mean'],lower=n_ci['lower']))
    
  }
}

print(df_ratios)
write.csv(df_ratios,paste0(output_folder,"cell_volume_ratios.csv"), row.names = FALSE)
if (any(is.na(df_ratios$ratio))) {
  print ("<WARNING> the rows below contain NA for volume-ratio.  These NA values were removed for the ratio calculation, but check if volume was correctly determined for these cells")
  print(df_ratios[is.na(df_ratios$ratio),])
  df_ratios <- df_ratios[!is.na(df_ratios$ratio),]
}
print(df_ci)
# viz
if (flag_single_category){
  ratioplot <- ggplot(df_ci, aes(x=ratioID, y=mean, group=1) )+
    geom_errorbar(width=.1, size=2,aes(ymin=lower, ymax=upper)) 
} else {
  ratioplot <- ggplot(df_ci, aes(x=ratioID, y=mean, group=1,colour=category) )+
    geom_errorbar(width=.1, size=2,aes(ymin=lower, ymax=upper)) 
}

#+ geom_col(width = 0.1,position = position_dodge(width=0.7)) 
#ratioplot <- ratioplot  +   ylim(0.4,1.0)

#element_text(family = NULL, face = NULL, colour = NULL, size = NULL, hjust = NULL, vjust = NULL, angle = NULL, lineheight = NULL, color = NULL, margin = NULL, debug = NULL)
ratioplot <-  ratioplot +  geom_hline(yintercept=1,color = "red",size=2)
ratioplot <- ratioplot  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                axis.title.x=element_blank(),axis.title.y=element_blank(),
                                axis.text.x= element_text(face="plain", size=16),
                                axis.text.y= element_text(face="plain", size=16))

ratioplot
ggsave(paste0(output_folder ,'_volume_ volume ratio.pdf'), device = "pdf")

ratioplot_stripped <- ratioplot

ratioplot_stripped <- strip_plot(ratioplot)
ratioplot_stripped
#ggsave(paste0(output_folder ,'_volume_ volume ratio(stripped).pdf'), device = "pdf")

# Calculate ratios with hierarchical model  ----
print(" Calculate ratios with hierarchical model")
#The previous ratioplot did not take the time correlation into consideration
# so now we will create the same ratioplot but using an hierarchical model
# this will affect the CI which are now calculated for the fixed effects like we did for volume

#     first approach : separate out the data and fit model for each (deprecated !)
#     problem : the CI are too conservative (not correct)
df_ratios$category <- as.factor(df_ratios$category)
l_rows <- list()
for (category_i in l_category) {
  print(category_i)
  df_ratio_category <- df_ratios[df_ratios['category']==category_i,]
  l_ratioID<- unique(df_ratio_category['ratioID'])[[1]]
  for (ratioID_i in l_ratioID) {
    print(ratioID_i)
    df_ratio_category_i <- df_ratio_category[df_ratio_category['ratioID']==ratioID_i,]
    
    # fit model
    vol_ratio.lmer <-lmer(ratio ~ 1   + t + (1 | repID), data = df_ratio_category_i)

    #get CI
    newdat4 <- expand.grid(t=1,ratio = 0)
    newdat4$ratio <- predict(vol_ratio.lmer,newdat4 ,re.form=NA)
    mm4 <- model.matrix(object=terms(vol_ratio.lmer),data=newdat4)
    pvar4 <- diag(mm4 %*% tcrossprod(vcov(vol_ratio.lmer),mm4))
    tvar4 <- pvar4+VarCorr(vol_ratio.lmer)$repID[1]
    cmult <- 1.96 
    newdat4 <- data.frame(
      newdat4
      , plo = newdat4$ratio-cmult*sqrt(pvar4)  #low limit 0.95 CI of fixed effect
      , phi = newdat4$ratio+cmult*sqrt(pvar4)  #high limit 0.95 CI of fixed effect
      , tlo = newdat4$ratio-cmult*sqrt(tvar4)  #low limit 0.95 CI of fixed effect + random effect
      , thi = newdat4$ratio+cmult*sqrt(tvar4) #high limit 0.95 CI of fixed effect + random effect
    )
    newdat4$category <- category_i
    newdat4$ratioID <- ratioID_i
    l_rows <- append(l_rows,list(newdat4))
  }
}
df_ratio_lmer <- do.call(rbind, l_rows)
write.csv(df_ratio_lmer,paste0(output_folder,"cell_volume_ratios_lmer.csv"), row.names = FALSE)

# second approach : fit entire model and calculate CI via bootstrap (preferred)
cond_ABprl  <- df_ratios['ratioID']=='ABpr/ABpl'
df_ratios_ABprl <- df_ratios[cond_ABprl,]
vol_ratio_full.lmer <-lmer(ratio ~ 1  + category + (1 | repID), data = df_ratios[cond_ABprl,])

print('--> visual check of residuals')
plot(vol_ratio_full.lmer) 
qqnorm(resid(vol_ratio_full.lmer))
qqline(resid(vol_ratio_full.lmer)) 
vol_ratio_full.lme<-lme(ratio~1  + category ,data=df_ratios[cond_ABprl,],random=~ 1  | repID)
ACFplot <- plot(ACF(vol_ratio_full.lme,maxLag = 5),alpha=0.05) # no time correlation residuals ACFplot
ACFplot

library(boot)
set.seed(42)
CI_lmer_compressed.fun <- function(data, idx)
{
  df <- data[idx, ]
  
  vol_ratio_full.lmer <-lmer(ratio ~ 1  + category + (1 | repID), data = df)
  newdat4 <- expand.grid(category="compressed",ratio = 0)
  predict(vol_ratio_full.lmer,newdat4 ,re.form=NA)
}
CI_lmer_uncompressed.fun <- function(data, idx)
{
  df <- data[idx, ]
  
  vol_ratio_full.lmer <-lmer(ratio ~ 1  + category + (1 | repID), data = df)
  newdat4 <- expand.grid(category="uncompressed",ratio = 0)
  predict(vol_ratio_full.lmer,newdat4 ,re.form=NA)
}
df_ratio_lmer_bootstrap  <- data.frame (ratio  = numeric(0),
                                        plo=numeric(0),
                                        phi=numeric(0),
                                        category=c(),
                                        ratioID=c())
for (ratioID_i in l_ratioID) {
  cond  <- df_ratios['ratioID']==ratioID_i
  df_ratios_i <- df_ratios[cond,]
  print(paste(nrow(df_ratios_i), ' is the nb of rows selected for bootstrapping'))
  print(paste('...bootstrapping ratios',ratioID_i,' compressed...'))
  bootstrap <- boot(df_ratios_i, CI_lmer_compressed.fun , R = 1000)
  #bootstrap
  #plot(bootstrap)
  boot_ci <- boot.ci(boot.out = bootstrap, type = c("norm")) #type = c("norm", "basic","perc", "bca")
  df_ratio_lmer_bootstrap <- rbind(df_ratio_lmer_bootstrap,
                                   data.frame(ratio  = bootstrap$t0,
                                              plo=boot_ci$normal[2],
                                              phi=boot_ci$normal[3],
                                              category=c('compressed'),
                                              ratioID=ratioID_i))
  
  print(paste('...bootstrapping ratios',ratioID_i,' uncompressed...'))
  bootstrap <- boot(df_ratios_i, CI_lmer_uncompressed.fun , R = 1000)
  boot_ci <- boot.ci(boot.out = bootstrap, type = c("norm")) #type = c("norm", "basic","perc", "bca")
  df_ratio_lmer_bootstrap <- rbind(df_ratio_lmer_bootstrap,
                                   data.frame(ratio  = bootstrap$t0,
                                                      plo=boot_ci$normal[2],
                                                      phi=boot_ci$normal[3],
                                                      category=c('uncompressed'),
                                                      ratioID=ratioID_i))
}






if (flag_single_category){
  ratioplot <- ggplot(df_ratio_lmer_bootstrap, aes(x=ratioID, y=ratio, group=1) )+
    geom_errorbar(width=.1, size=2,aes(ymin=plo, ymax=phi)) 
} else {
  ratioplot <- ggplot(df_ratio_lmer_bootstrap, aes(x=ratioID, y=ratio, group=1,colour=category) )+
    geom_errorbar(width=.1, size=2,aes(ymin=plo, ymax=phi)) 
}
ratioplot <-  ratioplot +  geom_hline(yintercept=1,color = "red",size=2)
ratioplot <- ratioplot  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                axis.title.x=element_blank(),axis.title.y = element_text(size=22),
                                axis.text.x= element_text(face="plain", size=22),
                                axis.text.y= element_text(face="plain", size=22))
ratioplot <- ratioplot + theme(legend.title = element_blank(),
                 # legend.position="bottom",legend.direction = "horizontal",
                 legend.text = element_text(size = 22)) +  ylab("Volume ratio")
                 #axis.text.x = element_text(angle=90,size=12))



ratioplot 
write.csv(df_ratio_lmer_bootstrap,paste0(output_folder,"cell_volume_ratios_lmer_bootstrap.csv"), row.names = FALSE)
ggsave(paste0(output_folder ,'_volume_ volume ratio_lmer.pdf'))

# ratios hyphothesis testing : 
# Test if z_orientation is significant (with random effects)
cond_ABprl <- df_ratios['ratioID']=='ABpr/ABpl'
null_model <- lmer(ratio ~ 1   + category + (1 | repID), data = df_ratios[cond_ABprl,] , REML = FALSE)
alt_model <- lmer(ratio ~ 1    + category + z_orientation + (1 | repID), data = df_ratios[cond_ABprl,], REML = FALSE)

#likelihood ratio test:
# anova(null.model, model1, test = “Chisquare”) # null is simple, constraint model
anov_test = anova(null_model, alt_model,test=TRUE)
anov_test
if (anov_test$`Pr(>Chisq)`[2]<0.05){
  print("<RESULT> z orientation has a significant effect on ABpr/l ratio (null model without z orientation is rejected)")
} else {
  print("<RESULT> z orientation has NO significant effect on ABpr/l ratio (null model without z orientation is not rejected)")
}

# test if category is significant (with random effects)
null_model <- lmer(ratio ~ 1   + t  + z_orientation + (1 | repID), data = df_ratios[cond_ABprl,] , REML = FALSE)
alt_model <- lmer(ratio ~ 1   + t + category  +  z_orientation  + (1 | repID), data = df_ratios[cond_ABprl,], REML = FALSE)

#likelihood ratio test:
# anova(null.model, model1, test = “Chisquare”) # null is simple, constraint model
anov_test = anova(null_model, alt_model,test=TRUE)
anov_test
if (anov_test$`Pr(>Chisq)`[2]<0.05){
  print("<RESULT> category has a significant effect on ABpr/l ratio (null model without z orienation is rejected)")
} else {
  print("<RESULT> category has NO significant effect on ABpr/l ratio (null model without z orienation is not rejected)")
}

print('<<< 4  Observations  >>>>')
#4  Observations ------------

print(' -->  Observations : embryo volumes')
## 4.1 embryo volumes (random effects) -----------
#  (only the b's , so without fixed effect)
cell_vol.lmer.re=ranef(cell_vol.lmer.final)
cell_vol.lmer.fe=fixef(cell_vol.lmer.final)
df_rep_vol <- cell_vol.lmer.re$repID
rownames(df_rep_vol)
ix_col = 1
for(col_i in names(cell_vol.lmer.fe)){
  print(col_i)
  # rand_eff = cell_vol.lmer.re
  if (ix_col ==1) {
    df_rep_vol[ix_col + 1] = unlist(cell_vol.lmer.re) + (cell_vol.lmer.fe[1])
  }
  else {
    if (col_i=='t' | col_i=='nb_pix_last_z') {}
    else {
      df_rep_vol[ix_col + 1] = unlist(cell_vol.lmer.re) + (cell_vol.lmer.fe[1] + cell_vol.lmer.fe[ix_col])
    }}
  
  ix_col = ix_col + 1
}
names(df_rep_vol)[1] <- 'random_effect'
names(df_rep_vol)[2:ncol(df_rep_vol)] <- names(cell_vol.lmer.fe)
df_rep_vol
write.csv(df_rep_vol,paste0(output_folder,"4_df_rep_vol_observations.csv"), row.names = FALSE)

print('# Stacked barplot with multiple groups')
library(tidyr)
library(reshape2)
df_rep_vol$repID <- row.names(df_rep_vol)
subset(df_rep_vol, select = -c(random_effect) )
df_rep_vol_long <-melt(subset(df_rep_vol, select = -c(random_effect) ), id.vars=c("repID"),variable.name="cell_name",value.name="volume")
ggplot(data=df_rep_vol_long, aes(x=repID, y=volume,fill=cell_name)) +
  geom_bar(stat="identity")
ggsave(paste0(output_folder,"_volume_Stacked barplot(embryo size).jpg"))


print(' --- > Observations : ratios of daughtercells')
#4.2ratios of daughtercells -----

#in [12] Fickentsher : MS/E = 1.2, C / P3 = 1.8
df_rep_vol$MS_E_ratio <- df_rep_vol$cell_nameMS / df_rep_vol$cell_nameE
print(paste0("<RESULT> The mean MS / E volume ratio is ", mean(df_rep_vol$MS_E_ratio),"(Fickentsher claims around 1.2)"))
df_rep_vol$C_P3_ratio <- df_rep_vol$cell_nameC / df_rep_vol$cell_nameP3
print(paste0("<RESULT> The mean C / P3 volume ratio is ", mean(df_rep_vol$C_P3_ratio),"(Fickentsher claims around 1.8)"))
df_rep_vol$ABpl_ABpr_ratio <- df_rep_vol$cell_nameABpl / df_rep_vol$cell_nameABpr
print(paste0("<RESULT> The mean ABpl / ABpr volume ratio is ", mean(df_rep_vol$ABpl_ABpr_ratio)))

