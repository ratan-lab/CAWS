#!/usr/bin/env Rscript 

library(tidyverse)
library(ggpubr)
library(corrplot)

reprod = c()
fragCount = NULL
for(hist in snakemake@input){
    filename = strsplit(hist, "/")[[1]][2]
    samplename = strsplit(filename, ".bin.bed")[[1]][1]
          
    if(is.null(fragCount)){
        fragCount = read.table(hist, header = FALSE) 
        colnames(fragCount) = c("chrom", "bin", samplename)
    }else{
        fragCountTmp = read.table(hist, header = FALSE)
        colnames(fragCountTmp) = c("chrom", "bin", samplename)
        fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
    }
}

M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs") 

constantfac <- dim(M)[1]/8
constantfac <- max(1, constantfac)
pdf(snakemake@output[[1]], width=7*constantfac, height=7*constantfac)
corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", addrect = 3, rect.col = "black", rect.lwd = 3,cl.pos = "b", tl.col = "indianred4", tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 2, number.cex = 1, col = colorRampPalette(c("midnightblue","white","darkred"))(100))
dev.off()
