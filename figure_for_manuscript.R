# 2022-03-29
# figures for manuscript

# volcano plot
# DEGs GO terms
# Heatmap
# GSEA GO terms

# colors to use: upregulated; darkseagreen3, downregulated; deeppink2
if (Sys.info()['sysname'] == "Darwin"){
  setwd("Dropbox/karlssonlab/projects/DBP/")
} else {
  setwd("F:/Dropbox/karlssonlab/projects/DBP")
}

library(dplyr)
library(DESeq2)
library(tximport)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(clusterProfiler)
library(org.Mm.eg.db)

# functions
## generate a neat table for results
make_resultstable <- function(num, obj = dds) {
  name <- resultsNames(obj)[num]
  res <- results(obj, name = name, alpha = 0.05, lfcThreshold = 0) # alpha = FDR correction
  
  d <- data.frame(res)
  d$ensembl_gene_id <- rownames(d)
  d$gene_name <- ens[match(d$ensembl_gene_id, ens$ensembl_gene_id), "gene_name"] %>% toupper
  
  d$sign <- "ns"
  d$sign[d$padj <= 0.05 & d$log2FoldChange > 0] <- "up"
  d$sign[d$padj <= 0.05 & d$log2FoldChange < 0] <- "down"
  rownames(d) <- NULL
  
  d
}

# ================================= Code ===================================== #
# load metadata
metadata <- read.csv("data/liver/liver_dpb_mouse_metadata.csv")

metadata$exposure <- factor(metadata$exposure, # adjust factor levels order
                            levels = c("zero_dose", "low_dose", "high_dose"))
metadata$diet <- factor(metadata$diet, levels=c("normal_fat", "high_fat"))
metadata$group <- factor(paste0(metadata$exposure, ".", metadata$diet),
                         levels=c("zero_dose.normal_fat", "zero_dose.high_fat",
                                  "low_dose.normal_fat", "low_dose.high_fat",
                                  "high_dose.normal_fat", "high_dose.high_fat"))

# import ensembl transcript -> gene id index file
tx2gene <- read.csv("data/liver/salmon_tx2gene.tsv", sep="\t")[, 1:2]
names(tx2gene) <- c("TXNAME", "GENEID") # rename

# create intermediate files for analysis
coldata <- metadata %>% filter(generation %in% "F0") # filter generations
files <- coldata %>% # vector of sample directories
  dplyr::select(path) %>%
  unlist 
names(files) <- coldata$sample_id # add sample names

# load transcript files
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# generate DESeq-object
dds <- DESeqDataSetFromTximport(txi = txi, colData = coldata, 
                                design = ~ exposure)
dds <- DESeq(dds)

# load ensembl gene data
ens <- read.delim("data/ensembl_annotations_mouse.txt")

# generate differential expression results
DE <- list()
#DE$low_dose <- make_resultstable(2) # exposure: low dose vs control animals
DE$high_dose <- make_resultstable(3) # exposure: high dose vs control animals


# ================================= figures ================================== #
deg <- DE$high_dose %>% 
  filter(padj <= 0.05) %>%
  dplyr::select(ensembl_gene_id, gene_name, sign)

deg$rownames <- deg$sign
deg$rownames <- sub("up", "Upregulated", deg$rownames)
deg$rownames <- sub("down", "Downregulated", deg$rownames)
deg$gene_name <- toupper(deg$gene_name)
deg$gene_name <- ifelse(is.na(deg$gene_name), "", deg$gene_name)

# DEG volcano
foo <- function(){
  library(ggrepel)
  ggplot() + 
    geom_point(data = subset(DE$high_dose, sign == "ns"), 
               aes(log2FoldChange, -log10(padj), color = sign), size = .5) +
    geom_hline(yintercept = -log10(.05), lty = 2, color = "navyblue", 
               alpha = .5) +
    geom_point(data = subset(DE$high_dose, sign != "ns"), 
               aes(log2FoldChange, -log10(padj), fill = sign), shape = 21, 
               size = 1.5) +
    scale_color_manual(values = c("ns"="grey40")) +
    scale_fill_manual(values = c("up"="darkseagreen3", "down"="deeppink2")) +
    theme_pubclean(base_size = 10) +
    theme(legend.position = "none", plot.title = element_text(hjust = .5)) + 
    labs(x = "log2 fold change", y = "-log10(adjusted pvalue)")
}
volcano_plot <- foo()

# DEG heatmap
foo <- function(){
  fontsize <- 8
  # color palette
  pal <- c("darkolivegreen2", "black", "firebrick2")
  ramp <- colorRamp2(c(-2, 0, 2), pal)
  
  genes <- deg$ensembl_gene_id
  
  #up; darkseagreen3, down; deeppink2
  
  treatment <- colData(dds)$exposure %>% as.vector
  samples <- treatment %in% c("zero_dose", "high_dose")
  
  exposure <- list(Exposure=c("Control" = "tomato1", "low_dose" = "darkorchid4", 
                              "High dose" = "darkslategray"))
  regulation <- list(Direction=c("Upregulated" = "darkseagreen4", 
                                 "Downregulated" = "deeppink4"))
  
  param <- list(labels_gp = gpar(fontsize = fontsize),
                title_gp = gpar(fontsize = fontsize, fontface = "bold"),
                border = "black",
                color_bar = "discrete",
                nrow = 1,
                title_position = "topcenter",
                #legend_direction = c("vertical"))
                legend_direction = "horizontal")
  
  dbp <- treatment[samples]
  dbp <- sub("zero_dose", "Control", dbp)
  dbp <- sub("high_dose", "High dose", dbp)
  
  # transform data
  p <- vst(dds)[genes, samples] %>% assay
  p <- t(scale(t(p), center=T, scale=T))
  rownames(p) <- deg$gene_name
  
  # annotation for the heatmap
  regulation_anno <- rowAnnotation(Direction = deg$rownames,
                                   col = regulation,
                                   show_annotation_name = FALSE, 
                                   simple_anno_size_adjust = TRUE, 
                                   width = unit(2, "mm"),
                                   border = TRUE,
                                   annotation_legend_param = param)
  
  exposure_anno <- HeatmapAnnotation(Exposure = dbp, 
                                     col = exposure, 
                                     show_annotation_name = FALSE, 
                                     simple_anno_size_adjust = TRUE, 
                                     height = unit(2, "mm"),
                                     border = TRUE,
                                     annotation_legend_param = param)
  
  hm <- Heatmap(p, 
                col = ramp, 
                name = "z-score", 
                cluster_columns = cluster_within_group(p, treatment[samples]), # manual clustering
                cluster_rows = TRUE,
                
                column_split = 2,
                column_title = NULL,
                
                bottom_annotation = exposure_anno,
                left_annotation = regulation_anno,
                
                row_title = NULL,
                row_split = 2,
                row_names_gp = gpar(fontsize = fontsize),
                column_names_gp = gpar(fontsize = fontsize),
                heatmap_legend_param = list(labels_gp = gpar(fontsize = fontsize),
                                            title_position = "topcenter",
                                            title_gp = gpar(fontsize = fontsize, 
                                                            fontface = "bold"),
                                            direction = "horizontal"),
                show_row_names = TRUE,
                show_column_dend = FALSE,
                show_row_dend = FALSE,
                border = T)
  
  grid.grabExpr(draw(hm, merge_legend = TRUE, heatmap_legend_side = "bottom"))
  #draw(hm, merge_legend = TRUE, heatmap_legend_side = "bottom")
}
deg_heatmap_plot <- foo()

# DEG gene ontology 
foo <- function(){
  out <- enrichGO(gene = deg$ensembl_gene_id,
                  OrgDb = org.Mm.eg.db,
                  keyType = "ENSEMBL",
                  ont = "ALL",
                  pvalueCutoff = .05,
                  pAdjustMethod = "BH",
                  universe = NULL,
                  readable = TRUE)
  
  
  # color = 'pvalue/p.adjust/qvalue'
  # x = 'GeneRatio/Count'
  # qvalue: proportion of false discovery rate in a population of significant 
  # p-values
  barplot(out, split = "ONTOLOGY", showCategory = 20, color = "qvalue", 
          x = "Count") +
    theme_pubr(base_size = 8) +
    theme(legend.position = "right") +
    facet_grid(ONTOLOGY ~., scales = "free", space = "free")
}
go_plot <- foo()

# gsea results
res <- DE$high_dose %>%
  dplyr::select(gene_name, stat) %>%
  na.omit() %>%
  filter(gene_name != "")

res <- res[!duplicated(res$gene_name), ]
res$symbol <- toupper(res$gene_name)

ranks <- res$stat
names(ranks) <- res$symbol
ranks <- sort(ranks, decreasing = TRUE)

gsea1 <- function(ranks, category, subcategory = NULL, pval = .05, to_rm = NULL,
                  y = "Normalised Enrichment score", title){
  set.seed(54321)
  #category <- "C5"
  #subcategory <- "GO:CC"
  species <- "Mus musculus"
  db <- msigdbr::msigdbr(species = species, category = category, 
                subcategory = subcategory)
  go_list <- split(x = toupper(db$gene_symbol), 
                   f = db$gs_name)

  res <- fgsea::fgsea(pathways = go_list, 
                      stats = ranks,
                      minSize = 15,
                      maxSize = 1000,
                      nperm = 10000)
  
  # collapse pathways, reduce redundancy
  collapsedPathways <- fgsea::collapsePathways(fgseaRes = res[res$padj < .05, ],
                                               pathways = go_list,
                                               stats = ranks)
  mainPathways <- res[pathway %in% collapsedPathways$mainPathways][order(-NES), 
                                                                   pathway]
  
  res <- res[pathway %in% mainPathways] %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval) %>% 
    arrange(desc(NES)) %>%
    mutate_if(is.numeric, round, digits = 4)
  
  if (!is.null(to_rm)) res$pathway <- stringr::str_replace(res$pathway, to_rm , 
                                                           "")
  res$Enrichment <- ifelse(res$NES > 0, "Up-regulated", "Down-regulated")
  
  print(paste0("number of pathways: ", nrow(res)))
  
  up <- subset(res, NES > 0)
  down <- subset(res, NES < 0)
  res <- rbind(head(up, n = 10), tail(down, n = 10))
  
  # trim the long names
  foo1 <- function(x, num){
    term <- strsplit(x, "")[[1]]
    underscores <- which(term %in% "_")
    n2 <- underscores[which.min(abs(underscores - num))]
    term[n2] <- "\n"
    paste(term, collapse="")
  }
  res$pathway <- unlist(lapply(res$pathway, function(x){
    num <- 40
    n <- nchar(x) > num
    if (n) foo1(x, num) else x
  }))
  
  # colors for the plot
  upcol <- colorRampPalette(colors = c("red4", "red1", 
                                       "lightpink"))(sum(res$Enrichment == 
                                                           "Up-regulated"))
  downcol <- colorRampPalette(colors = c("lightblue", "blue1", 
                                         "blue4"))(sum(res$Enrichment == 
                                                         "Down-regulated"))
  col <- c(upcol, downcol)
  names(col) <- 1:length(col)
  
  labelupcol <- c(rep("white", 5), rep("black", 5))[1:length(upcol)]
  labeldowncol <- rev(labelupcol[1:length(downcol)])
  labelcol <- c(labelupcol, labeldowncol)
  names(labelcol) <- 1:length(labelcol)
  
  labelcol <- c(labelupcol, labeldowncol)
  names(labelcol) <- 1:length(labelcol)
  
  res$index <- as.factor(1:nrow(res))
  res$edge_len <- unlist(lapply(res$leadingEdge, function(x) unlist(x) %>% length))
  
  ggplot(res, aes(reorder(pathway, NES), NES)) +
    geom_col(aes(fill = index)) + 
    geom_text(aes(label = edge_len, y = NES * .9), 
              color = labelcol, parse = TRUE, size = 3) + 
    geom_text(aes(label = edge_len, y = NES * .9), 
              color = labelcol, parse = TRUE, size = 3.05) + 
    geom_text(aes(label = edge_len, y = NES * .9),
              color = labelcol, parse = TRUE, size = 3.08) + 
    
    scale_fill_manual(values = col) + 
    coord_flip() + 
    labs(x= "", y = y, title = title) +
    #theme_minimal(base_size = 8) +
    theme_bw(base_size = 8) + 
    #theme_pubr(base_size = 8) +
    theme(legend.position = "none", 
          axis.text.y = element_text(color = "black"),
          axis.text.x = element_text(color = "black", vjust = 0.5))
}
p1 <- gsea1(ranks, "C5", "GO:BP", to_rm = "GOBP_", 
            title = "Biological Processes", y = NULL)
p2 <- gsea1(ranks, "C5", "GO:CC", to_rm = "GOCC_", 
            title = "Cellular Components", y = NULL)
p3 <- gsea1(ranks, "C5", "GO:MF", to_rm = "GOMF_", 
            title = "Molecular Functions")


# combine into one big plot
AB <- plot_grid(volcano_plot, go_plot, ncol = 2, rel_widths = c(.8, 1), 
                labels = c("A", "B"))
ht_plot <- plot_grid(NULL, deg_heatmap_plot, ncol = 2, rel_widths = c(.05, 1), 
                     labels = "C")
ABC <- plot_grid(AB, ht_plot, nrow = 2, rel_heights = c(.6, 1))
DEF <- plot_grid(p1, p2, p3, nrow = 3, align = "v", labels = c("D", "E", "F"))
ABCDEF <- plot_grid(ABC, DEF, ncol = 2, rel_widths = c(1, 1))

pdf("manuscript/figures/DEG5.0.pdf", width = 13, height = 9)
ABCDEF
dev.off()

# Figure legend
# (A) volcanoplot, x-axis: log2 fold change, y-axis: -log10 adjusted p-value. 
# DEG are highlighted in pink (downregultaed) and green (upregulated)
# (B) Gene ontology, DEG correlated with gene ontology terms showing enrichment.
# y-axis: GO terms, x-axis: number of DEGs found within the GO term, color of 
# bars indicate the q-value
# (C) Heatmap of DEGs, down- and upregulated genes are highlighted in pink 
# (downregultaed) and green (upregulated) on the left side. Experimental groups 
# are indicated in the bottom, orange (control) and darkgrey (high dose group). 
# The individual color for each DEG represents the z-score, ranging from -2 
# (green) to +2 (red). Z-score is a measure of the distance from the groups 
# mean. Gene names are presented on the right side.
# GSEA results using ontology gene sets (D) Biological processes, (E) Cellular 
# components and (F) Molecular functions downloaded from MSigdb: 
# http://www.gsea-msigdb.org/gsea/msigdb/index.jsp. Only the top and bottom 10 
# most significant resutls are shown. Numbers inside the bars represent the 
# number of genes in the leading edge



# ========================== include phenotypes? ============================= #

# load phenotypes
#phenotypes <- read.csv("data/dbp_project_phenotypes.csv")
