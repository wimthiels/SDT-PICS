# select
if (!exists("l_cell_name") | anyNA(l_cell_name)  | length(l_cell_name)==0){
  l_cell_name <- unique(geom_info_full['cell_name'])[[1]]
  l_cell_name <- c('ABal','ABar','ABpl','ABpr','E','MS','P2','C','P3') #also sets sort order
}
if (!exists("l_emb_ID") | anyNA(l_emb_ID) | length(l_emb_ID)==0){
  l_emb_ID <- unique(geom_info_full['embID'])[[1]]
  #l_emb_ID <- c('emb_1','emb_2','emb_3','emb_4', 'emb_5','emb_6','emb_7','emb_8','emb_9')
}
if (!exists("l_rep_ID") | anyNA(l_rep_ID)  | length(l_rep_ID)==0){
  l_rep_ID <- unique(geom_info_full['repID'])[[1]]
}

if (!exists("l_time_interval") | anyNA(l_time_interval)  | length(l_time_interval)==0){
  l_time_interval <- list(min(geom_info_full['t'])[[1]],max(geom_info_full['t'])[[1]])
}


print('Data will be filtered by  : ')
print('cell_name')
print(l_cell_name)
print('emb_ID')
print(l_emb_ID)
print('rep_ID:')
print(l_rep_ID)
print('time_interval:')
print(l_time_interval)