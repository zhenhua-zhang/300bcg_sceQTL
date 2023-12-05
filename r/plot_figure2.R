library(ggplot2)
library(readr)
library(dplyr)


gt.mtx2 = read.table("./gt_mtx.txt", row.names = 1)
snp38 = data.frame(read_tsv("./snp.info.txt"))

### boxplot: LCP1, IL6R ###

Snp = "rs2806897"
gt_sub = gt.mtx2[match(Snp, rownames(gt.mtx2)), ]
snp_info = snp38 %>% filter(snp == Snp)

features = "LCP1"

lcp1.exp = exp.ti[which(rownames(exp.ti) == "LCP1"), ]
gt.test = as.vector(unlist(gt_sub))[match(ids.shared, colnames(gt_sub))]
exps = as.vector(lcp1.exp)[match(ids.shared, names(lcp1.exp))]


bp = data.frame(gt = gt.test,
                exp = exps)

bp$gt = gsub("0", paste0(snp_info$ref,snp_info$ref), bp$gt)
bp$gt = gsub("1", paste0(snp_info$ref,snp_info$alt), bp$gt)
bp$gt = gsub("2", paste0(snp_info$alt,snp_info$alt), bp$gt)



pdf("./test.pdf")
ggboxplot(
  bp, x = "rs2806897_GT", y = "exp", fill = "rs2806897_GT", size = 0.5, 
  palette = c("#04458B", "#D92F1F", "#63B350"),
  add = "ggscatter", width = 0.6, xlab = "Genotypes", ylab = "FC of Expression: LCP1") 
dev.off()

paste(features, rownames(gt_sub), sep = "_")

### locuszoom ###

source("./LocusZooms-master/functions/locus_zoom.R")

eqtl <- read.delim("./lcp1-locus.txt", stringsAsFactors = FALSE, header = TRUE)
Unique.genes <- read.delim("./lcp-locus-genelist.txt", stringsAsFactors = FALSE, header = TRUE)
ld = read.table("./plink.ld", stringsAsFactors = FALSE, header = TRUE)
ld = ld %>% filter(SNP_A == "rs2806897")

locus.zoom(data = eqtl,                                    
           region = c(13, 46000000, 47000000),                             
           gene = "LCP1",
           snp = "rs2806897", 
           offset_bp = 0,                                                  
           ld.file = ld,                                           
           genes.data = Unique.genes,			                  
           plot.title = "LCP1-rs2806897",       
           file.name = "./lcp1_locus.jpg"
) 





