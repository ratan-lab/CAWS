#!/usr/bin/env Rscript

library(tidyverse)
library(ggpubr)
library(viridis)

fragLen=c()

for(hist in snakemake@input){
    filename = strsplit(hist, "/")[[1]][2]
    samplename = strsplit(filename, ".fraglengths.txt")[[1]][1]

    fragLen=read.table(hist, header=FALSE) %>% 
            mutate(fragLen=V1 %>% as.numeric, 
                   fragCount=V2 %>% as.numeric, 
                   Weight=as.numeric(V2)/sum(as.numeric(V2)), 
                   sampleInfo=samplename) %>% 
            rbind(fragLen, .)
}
    
fragLen$sampleInfo=factor(fragLen$sampleInfo)

figA=fragLen %>% 
        ggplot(aes(x=sampleInfo, y=fragLen, weight=Weight, fill=sampleInfo)) +
        geom_violin(bw=5) +
        scale_y_continuous(breaks=seq(0, 800, 50)) +
        scale_fill_viridis(discrete=TRUE, begin=0.1, end=0.9, option="magma", alpha=0.8) +
        scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9) +
        theme_bw(base_size=20) +
        ggpubr::rotate_x_text(angle=20) +
        ylab("Fragment Length") +
        xlab("")

figB=fragLen %>% 
        ggplot(aes(x=fragLen, y=fragCount, color=sampleInfo, group=sampleInfo)) +
        geom_line(size=1) +
        scale_color_viridis(discrete=TRUE, begin=0.1, end=0.9, option="magma") +
        theme_bw(base_size=20) +
        xlab("Fragment Length") +
        ylab("Count") +
        coord_cartesian(xlim=c(0, 500))

constantfac <- length(unique(fragLen$sampleInfo))/12 
constantfac <- max(1, constantfac)
pdf(snakemake@output[[1]], width=14*constantfac, height=10)
ggarrange(figA, figB, ncol=2, common.legend = TRUE, legend="bottom")
dev.off()
