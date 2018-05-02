dir="/home/labs/amit/diklag/qc_paper/"
setwd(dir)
source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")
library(GEOquery)
library(Seurat)
require(XML)

###########################################
##### getting mouce cell atlas meta data

gse = getGSEDataTables("GSE108097")
samples = gse@header$sample_id
l = list()
for(sample in samples){
  l[sample] = getGEO(sample,GSEMatrix = FALSE)
}

r = sapply(l, function(y) data.frame(lapply(y@header, function(x) t(x)),stringsAsFactors = FALSE))
for(i in 1:length(r)){
  fields = c("characteristics","relation")
  for(field in fields){
    cols_ch = grep(field,colnames(r[[i]]))
    if(length(cols_ch) !=0){
      col_repl = apply(r[[i]][cols_ch],1, function(x) sapply(strsplit(x, ": "), "[[", 1))
      names(r[[i]])[cols_ch] = as.character(trimws(col_repl))
      cell_repl = trimws(apply(r[[i]][cols_ch],1, function(x) sapply(strsplit(x, ": "), "[[", 2)))
      if(length(cols_ch) ==1){
        r[[i]][cols_ch[1]] = cell_repl
      }
      else{
        for(j in 1:length(cols_ch)){
          r[[i]][cols_ch[j]] = cell_repl[j,]
        } 
      }
    }
  }
}

cols = unique(unlist(lapply(r, function(x) names(x))))
df = df <- data.frame(matrix(ncol = length(cols), nrow = 0))
colnames(df) = cols
for(i in 1:length(r)){
  current_cols = intersect(names(r[[i]]),cols)
  print(r[[i]][current_cols] )
  df[names(r)[i],current_cols]=r[[i]][current_cols]
}
meta_data = df
write.csv(meta_data,file="GSE108097_RAW/meta_data.csv")

###########################################

### generate seurat
non_info_cols = grep("channel_count|data_row_count|supplementary|data_processing|contact|extract_protocol|taxid|library|molecule",colnames(meta_data))
meta = meta_data[,-non_info_cols]
gsm_files = list.files("GSE108097_RAW/")
setwd("GSE108097_RAW/")
s_list = list()
for(gsm in rownames(meta_data)){
  if(gsm == 'GSM2906480'){
    next
  }
  pos = grep(gsm,gsm_files)
  path = gsm_files[pos]
  print(gsm)
  print(path)
  md = meta[gsm,]
  umi_tab = read.table(path)
  md = merge(colnames(umi_tab),md,all.x=TRUE)
  rownames(md) = md$x
  md = md[,-grep("x",colnames(md))]
  if(is.na(s_list[])){
    s_list[[i]] = CreateSeuratObject(raw.data = umi_tab, project = amp_batch,meta.data = meta_data)
  }
  else{
    new.data = CreateSeuratObject(raw.data = umi_tab, project = amp_batch, meta.data = meta_data)
    s_list[[i]] =  MergeSeurat(object1 = s_list[[i]], object2 = new.data,project = as.character(organs[i]))
  }
}


######## mammary gland tables
non_info_cols = grep("channel_count|data_row_count|supplementary|data_processing|contact|extract_protocol|taxid|library|molecule",colnames(meta_data))
meta = meta_data[,-non_info_cols]
exps = rownames(meta)[grep("MammaryGland",meta$source_name_ch1)]
md = meta[exps,]
dir = "/home/labs/amit/diklag/qc_paper/GSE108097_RAW"
s_list = get_s_list(exps,md,dir)
dir.create("mouse_atlas")
dir.create("mouse_atlas/saved_work")
save(s_list,file="mouse_atlas/saved_work/seurat_list_mammary_gland.RData")
for(i in 1:length(s_list)){
  s_list[[i]] = FilterCells(s_list[[i]], subset.names = "nGene", low.thresholds = 100, high.thresholds = Inf)
  s_list[[i]] = NormalizeData(s_list[[i]])
  s_list[[i]] = ScaleData(s_list[[i]], display.progress = F)
  s_list[[i]] = FindVariableGenes(s_list[[i]], do.plot = F)
}

save(s_list,file="mouse_atlas/saved_work/seurat_list_mammary_gland_normalized.RData")

hvg_l = lapply(s_list,function(x) head(rownames(x@hvg.info),2000))
gene.use = unique(unlist(hvg_l))
for(i in 1:length(s_list)){
  gene.use = intersect(gene.use, rownames(s_list[[i]]@scale.data))
}



###### function for SEURAT list
get_s_list <- function(table_list,meta_table,dir){
  s_list = list()
  projects_list = list.files(dir)
  if(substr(dir,start=nchar(dir),stop=nchar(dir))!= "/"){
    dir = paste0(dir,"/")
  }
  for(project in table_list){
    pos = grep(project,projects_list)
    path = projects_list[pos]
    umi_tab = read.table(paste0(dir,path))
    m_data = merge(colnames(umi_tab),meta_table[project,],all.x=TRUE)
    rownames(m_data) = m_data$x
    m_data = m_data[,-grep("x",colnames(m_data))]
    s_list[project] = CreateSeuratObject(raw.data = umi_tab, project = project,meta.data = m_data)
  }
  s_list
}
