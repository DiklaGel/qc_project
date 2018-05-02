library(devtools)
library(tgstat)
library(tgconfig)
load_all("/home/labs/amit/diklag/sc_packages/scrdb_old/")

setwd("/home/labs/amit/diklag/qc_paper")
scdb="/home/labs/amit/eyald/sc_pipeline/scdb_hisat_mouse"
umi_dir = "/home/labs/amit/eyald/sc_pipeline/scdb_hisat_mouse/output/umi.tab"
amp_batches = "/home/labs/amit/eyald/sc_pipeline/scdb_hisat_mouse/annotations/amp_batches.txt"
index_fn = "/home/labs/amit/diklag/qc_paper/files/index_fn.txt"
batches_stats = read.delim(index_fn)
#batches = unique(batches_stats[batches_stats$tissue == "spleen" & batches_stats$area %in% c("B","T"),]$amp_batch)
gene_markers = as.vector(read.delim("files/genes.csv",header = FALSE)[[1]])

data = read.csv(file = "files/merged_immune_db.csv",row.names = 1)
batches = unique(data$Amp_batch_ID)
data[,gene_markers]=NA
for(amp_batch in batches){
  path = paste0(scdb,"/output/umi.tab/",amp_batch,".txt")
  if(!file.exists(path)){
    print(paste0(path," not exsits!"))
    next
  }
  umi_tab = read.delim(path)
  cells = intersect(colnames(umi_tab),rownames(data))
  data[cells,gene_markers] = t(umi_tab[gene_markers,cells])
}

genes_shades = colorRampPalette(c("white", "orange", "tomato","mediumorchid4", "midnightblue"))(1000)

data$cell_type = ordered(data$cell_type,levels=cell_types)
cell_ord = order(data$cell_type)
umis = t(data[,gene_markers])
umis_n = sweep(umis,2,colSums(umis),"/") * 500
foc = log(1 + 7 * umis_n)
cls = cumsum(table(data$cell_type)) / length(rownames(data))
cell_types = c("B","B_CD21_hi","B_gc","ILC","CD4" ,"T","iTreg","Treg","CD8","CD8_lo", "CD_hl","CD8_activated","NK","Baso","Mast","pDC",
               "DC","PC","Mon","Mono_Ly6c_lo","Mono_Ly6c_hi","eMac","Mac","AM","Neut" )

mm = melt(data[c("cell_type",gene_markers)], id=c("cell_type"))
ggplot(mm)+geom_boxplot(aes(x=paste(variable,factor.col,sep="_"), y=value))

ggplot(mm)+geom_boxplot(aes(x=variable, y=value))+facet_grid(.~cell_type)


png("saved work/heatmap.png", height = 1500, width = 2475)
par(fig = c(0,1,0,0.9), mar = c(4,10,0,0))
image(t(foc[gene_markers, cell_ord]), col = genes_shades, axes = F)
mtext(gene_markers, side = 2, las = 2, at = (1 - seq_along(gene_markers)) / (1 - length(gene_markers)))
abline(v = cls, lty = 2)
mtext(names(cls), side = 1, at = rowMeans(cbind(as.vector(c(0,cls[-length(cls)])),cls)))
dev.off()

modules = read.csv("/home/labs/amit/diklag/singleCell/lung_20170910/bitbucket/lung_development/files/modules_big.csv", stringsAsFactors = F)
cc_genes = modules[ modules$annotation == "CC", "gene"]
ribo_genes = modules[ modules$annotation == "Ribo", "gene"]
bad_genes = unique(c(cc_genes, ribo_genes, "Malat1", "7SK", "Xist", "mmu-mir-689-2"))

mars_batches = batches_stats[batches_stats$amp_batch %in% batches,]

index = batches_stats[batches_stats$amp_batch %in% batches,"amp_batch"]

lung_sc_mat = scm_read_scmat_mars(base_dir = umi_dir,
  mars_batches = index ,
  batch_meta = batches_stats,
  min_umis_n = 200)

sc_mat = sc_pipe_clean(index_fn = index_fn,
                       batch_meta_attr = "amp_batch",
                       base_dir= umi_dir,
                       mars_batches = batches,
                       outdir = "tg_comb_spleen",
                       filt_amb_on_clusts = T,
                       filt_outliers_on_clusts = T,
                       remove_outliers_before_clustering = F,
                       mark_blacklist_terms = bad_genes,
                       amb_epsilon = 0.03)

load(file = "saved work/sc_mat_spleen.RData")
