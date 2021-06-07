#set input parameters
if (length(myArgs)==0) {
  print("Input arguments filled with default values:")
  input_file <- paste0(DATA_FOLDER,"OUTPUT/7cellAgg/agg_tc.xlsx")
  output_folder <- paste0(DATA_FOLDER,"OUTPUT/7cellAna/")
} else {
  prm <- read.csv(myArgs[1],header=FALSE,row.names=1,stringsAsFactors = FALSE)
  print(paste("prm=",prm))
  input_file <- prm['input_file',]
  output_folder <- prm['output_folder',]
  l_cell_name <- parse_as_Rlist(prm['l_cell_name',])
  l_rep_ID <- parse_as_Rlist(prm['l_rep_ID',])
  l_emb_ID <- parse_as_Rlist(prm['l_emb_ID',])
  l_time_interval <- parse_as_Rlist(prm['l_time_interval',])
}
print(paste0("Input file: ", input_file))
print(paste0("output_folder: ", output_folder))
print(l_time_interval)
