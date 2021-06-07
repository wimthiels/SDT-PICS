#helper functions

library(stringr)

strip_plot <-function (plot,strip="axis") {
  library(ggplot2)
  if (strip=="full"){
    plot_stripped <- plot + theme_void()  #removes labels, but also axes
    return(plot_stripped)
  }
  plot_stripped <- plot + theme(axis.title.x=element_blank(),
                                axis.text.x=element_blank(),
                                axis.title.y=element_blank(),
                                axis.text.y=element_blank(),
                                legend.text =  element_blank(), #removes legend labels
                                legend.background = element_rect(fill = "white"),
                                legend.key = element_rect(fill = "white", color = NA)
                                #axis.ticks.x=element_blank(),
                                #axis.ticks.y=element_blank()
  )
  return(plot_stripped)
}


get_p_values_lmer <-function (lmer) {
  #kenward rogers is used to get degrees of freedom for the t-values
  print(paste('packageVersion(pbkrtest): ',packageVersion("pbkrtest")))
  #lmer.df.KR <- pbkrtest::get_ddf_Lb(lmer, fixef(lmer))
  lmer.df.KR <- pbkrtest::get_Lb_ddf(lmer, fixef(lmer))
  lmer.coef=coef(summary(lmer))
  lmer.p.KR <- cbind(lmer.coef,'p-value'=2 * (1 - pt(abs(lmer.coef[,3]), lmer.df.KR)))
  lmer.p.KR 
}

parse_as_Rlist <- function(l_python) {
  # input is a python list as a character string. some parsing to turn it into
  # an R list (vector of strings)
  remove_these <- c("\\'","\\[","\\]"," ")
  for (i in 1:length(remove_these)) {
    l_python <- str_remove_all(l_python, remove_these[i])
  }
  l_R <- strsplit(l_python, ",",fixed=TRUE)[[1]] 
  l_R
  
}


