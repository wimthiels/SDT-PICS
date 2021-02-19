#helper functions
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
  lmer.df.KR <- get_ddf_Lb(lmer, fixef(lmer))
  lmer.coef=coef(summary(lmer))
  lmer.p.KR <- cbind(lmer.coef,'p-value'=2 * (1 - pt(abs(lmer.coef[,3]), lmer.df.KR)))
  lmer.p.KR 
}