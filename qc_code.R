dir="/home/labs/amit/diklag/qc_paper/"
setwd(dir)
scdb="/home/labs/amit/eyald/sc_pipeline/scdb_hisat_mouse"
wells_cells=read.delim("/home/labs/amit/eyald/sc_pipeline/scdb_hisat_mouse/annotations/wells_cells.txt")
df = read.csv("summary/cells_table.csv",row.names = 1)
ds = merge(df,wells_cells[,c("Well_ID","Amp_batch_ID")],by="Well_ID")
rownames(ds)=ds$Well_ID
ds=ds[-1]
ds$umis_count=0
abs = unique(ds$Amp_batch_ID)
umis = data.frame()
for(amp_batch in abs){
  path = paste0(scdb,"/output/umi.tab/",amp_batch,".txt")
  if(!file.exists(path)){
    print(paste0(path," not exsits!"))
    next
  }
  umi_tab = read.delim(path)
  #c = colSums(umi_tab)
  #ds[names(c),"umis_count"]=as.integer(as.character(c))
  if(length(umis) == 0){
    umis = umi_tab
  }
  else{
    umis = cbind(umis,umi_tab)
  }
}


db = read.csv("cells_table_counts.csv")
index_fn = read.delim("files/index_fn.txt")

merged_db = merge(db,index_fn, by.x = "Amp_batch_ID",by.y = "amp_batch")
rownames(merged_db) = merged_db$X
merged_db = merged_db[-c(2,6)]
colnames(merged_db)[5] = "umi_count"
write.csv(merged_db, file = "files/merged_immune_db.csv")


db = read.csv("files/merged_immune_db.csv")

ord = order(db$umi_count,decreasing = TRUE)
scheme=read.csv("files/scheme.csv",row.names = 1)
colors = unique(scheme[,c("general.cell.type","general.color")])
cols = colors$general.color
names(cols)=colors$general.cell.type

library(ggplot2)
library(scales)
good_cells=rownames(db[ord,])[-1]
data = db
data$general.cell.type = as.factor(as.character(scheme[as.character(data$cell_type),"general.cell.type"]))
types = levels(data$general.cell.type)
g = ggplot(subset(data,!is.na(general.cell.type)), aes(x=general.cell.type, y=umi_count, fill=general.cell.type)) + scale_x_discrete(limits=types) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                                                                                                                                                   labels = trans_format("log10", math_format(10^.x)))

p<- g + geom_boxplot()+scale_fill_manual(values=as.character(cols[types])) + facet_grid(tissue ~.) + 
  theme(text =element_text(size=18),legend.text = element_text(size=rel(3)),axis.text = element_text(size=rel(1))) + guides(fill = guide_legend(ncol = 2))
p 

png("saved work/umi_cont_vs_cell_type.png",width = 1440,height = 1051)
p
dev.off()

cell_ord = order(data$cell_type)
g = ggplot(subset(data,!is.na(cell_type)), aes(x=general.cell.type, y=umi_count, fill=general.cell.type)) + scale_x_discrete(limits=types) + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
p<- g + geom_boxplot()+scale_fill_manual(values=as.character(cols[types])) + facet_grid(tissue ~.) + 
  theme(text =element_text(size=18),legend.text = element_text(size=rel(3)),axis.text = element_text(size=rel(1))) + guides(fill = guide_legend(ncol = 2))
p 
