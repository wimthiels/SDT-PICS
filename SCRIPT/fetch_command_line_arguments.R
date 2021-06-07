myArgs <- commandArgs(trailingOnly = TRUE)
if (length(myArgs)==0){
  myArgs <- c('/home/wth/Downloads/testinfra/OUTPUT/7cellAna/param.csv')
  # myArgs <- c('/home/wth/Downloads/testinfra/ISOLATED/Cshaper/CShaper/Data/active_contour_effect_on_volume/original(manual_exterior)/param_original.csv')
  # myArgs <- c('/home/wth/Downloads/testinfra/INPUT/uncomp/param.csv')
}
print(paste('loading parameters from file ', myArgs))
