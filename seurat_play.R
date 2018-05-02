dir="/home/labs/amit/diklag/qc_paper/"
setwd(dir)
scdb="/home/labs/amit/eyald/sc_pipeline/scdb_hisat_mouse"
wells_cells=read.delim("/home/labs/amit/eyald/sc_pipeline/scdb_hisat_mouse/annotations/wells_cells.txt")
amp_batches_file = paste0(scdb,"/annotations/amp_batches.txt")

amp_file = read.delim(amp_batches_file)

scheme=read.csv("files/scheme.csv",row.names = 1)
colors = unique(scheme[,c("general.cell.type","general.color")])
cols = colors$general.color
names(cols)=colors$general.cell.type

data = read.csv("files/merged_immune_db.csv",row.names = 1)

data$general.cell.type = scheme[as.character(data$cell_type),"general.cell.type"]
abs = unique(data$Amp_batch_ID)
organs = unique(data$tissue)

library(Seurat)
library(dplyr)


s_list = list()
for(i in 1:length(organs)){
  s_list[[i]] = NA
  abs = unique(data[data$tissue == organs[i],]$Amp_batch_ID)
  for(amp_batch in abs){
    path = paste0(scdb,"/output/umi.tab/",amp_batch,".txt")
    if(!file.exists(path)){
      print(paste0(path," not exsits!"))
      next
    }
    umi_tab = read.delim(path)
    meta_data =  data.frame(well_id = colnames(umi_tab), Amp_batch_ID = amp_batch)
    meta_data$tissue = organs[i]
    meta_data = merge(meta_data,amp_file,all.x = TRUE,by="Amp_batch_ID")
    rownames(meta_data) = meta_data$well_id
    meta_data = merge(meta_data,data[,c("cell_type","general.cell.type","Replicate")],by=0,all.x=TRUE)
    rownames(meta_data) = meta_data$well_id
    meta_data = meta_data[-c(1,3)]
    if(is.na(s_list[[i]])){
      s_list[[i]] = CreateSeuratObject(raw.data = umi_tab, project = amp_batch,meta.data = meta_data)
    }
    else{
      new.data = CreateSeuratObject(raw.data = umi_tab, project = amp_batch, meta.data = meta_data)
      s_list[[i]] =  MergeSeurat(object1 = s_list[[i]], object2 = new.data,project = as.character(organs[i]))
    }
  }
}

old_list = s_list
save(s_list,file="saved_work/seurat_list.RData")

for(i in 1:length(organs)){
  s_list[[i]] = FilterCells(s_list[[i]], subset.names = "nGene", low.thresholds = 100, high.thresholds = Inf)
  s_list[[i]] = NormalizeData(s_list[[i]])
  s_list[[i]] = ScaleData(s_list[[i]], display.progress = F)
  s_list[[i]] = FindVariableGenes(s_list[[i]], do.plot = F)
}

save(s_list,file="saved_work/seurat_list_normalized.RData")

load("saved_work/seurat_list_normalized.RData")
organs = unlist(lapply(s_list, function(x) x@project.name))

g.1 <- head(rownames(s_list[[1]]@hvg.info), 2000)
g.2 <- head(rownames(s_list[[2]]@hvg.info), 2000)
g.3 <-head(rownames(s_list[[3]]@hvg.info), 2000)
genes.use <- unique(c(g.1, g.2, g.3))
for(i in 1:length(organs)){
  genes.use = intersect(genes.use, rownames(s_list[[i]]@scale.data))
}

data.combined <- RunMultiCCA(object.list = s_list, genes.use = genes.use, num.ccs = 30)

save(data.combined,file="saved_work/data.combined.RData")


p1 <- DimPlot(object = data.combined, reduction.use = "cca", group.by = "tissue", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = data.combined, features.plot = "CC1", group.by = "tissue", 
              do.return = TRUE)
plot_grid(p1, p2)

PrintDim(object = data.combined, reduction.type = "cca", dims.print = 1:3, 
         genes.print = 20)

p3 <- MetageneBicorPlot(data.combined, grouping.var = "tissue", dims.eval = 1:30, 
                        display.progress = FALSE)

load("saved_work/after_alignment.RData")
aln = immune.combined

aln <- RunTSNE(aln, reduction.use = "cca.aligned", dims.use = 1:20, 
                           do.fast = T)
aln <- FindClusters(aln, reduction.type = "cca.aligned", 
                                resolution = 0.6, dims.use = 1:20)

# Visualization
p1 <- TSNEPlot(aln, do.return = T, pt.size = 0.5, group.by = "tissue")
p2 <- TSNEPlot(aln, do.label = T, do.return = T, pt.size = 0.5)
plot_grid(p1, p2)

marker_list = list()
for(i in unique(aln@meta.data$res.0.6)){
  marker_list[i] = FindConservedMarkers(aln, ident.1 = as.integer(i), grouping.var = "tissue", 
                                  print.bar = FALSE)
}



g = ggplot(subset(combined.data@meta.data,!is.na(organ)), aes(x=nUMI, y=nGene)) + geom_point()
g + facet_grid(.~organ)

g = ggplot(subset(filtered.data@meta.data,!is.na(general.cell.type)), aes(x=nUMI, y=nGene)) + geom_point()
g + facet_grid(general.cell.type~.)

