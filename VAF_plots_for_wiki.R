#!/usr/bin/env R

library(ggplot2)
library(data.table)
library(dplyr)
library(reshape)
library(tidyverse)

setwd("/mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq/ARCHER_processed/processed_variants")
files <- list.files(pattern = ".tsv")
temp <- lapply(files, fread, sep="\t", header = T)
data <- rbindlist( temp )
data <- as.data.frame(data)
data$key <- paste(data$PreferredSymbol, data$base_substitution, sep=" ")
colnames(data)[33] <- "AUAO" # Non unique column field (UAO) - taken from ARCHER sequencing information and probably not required, may remove in processing
dim(data)

setwd("/mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq/analysis/GrantApplication/")


output_filename <- "LBC_ARCHER.GrantApplication.VAF_Histogram.png"
png(output_filename, width = 1000, height = 1100, res = 200, pointsize =5)
p <- ggplot(data, aes(x=AF, colour = as.factor(wave))) + geom_histogram(bins = 100) + xlim(0, .6)
p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + ylab("Cumulative Count")
p <- p + xlab("Variant Allele Fraction")
p <- p + theme(axis.text.x = element_text(color="black", size=12),
               axis.text.y = element_text(color="black", size=12))
p <- p + theme(legend.title=element_text(), legend.position="bottom") + scale_colour_manual("Wave",values=c("red","orange", "green", "blue"))
p <- p + theme(axis.text=element_text(size=14),
               axis.title=element_text(size=14))
print(p)
dev.off()


output_filename <- "LBC_ARCHER.GrantApplication.VAF_Density.png"
png(output_filename, width = 1000, height = 1100, res = 200, pointsize =5)
p <- ggplot(data, aes(x=AF, colour = as.factor(wave), group = as.factor(wave))) + geom_density() + xlim(0, .45)
p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + ylab("Density")
p <- p + xlab("Variant Allele Fraction")
p <- p + theme(axis.text.x = element_text(color="black", size=12),
               axis.text.y = element_text(color="black", size=12))
p <- p + theme(legend.title=element_text(), legend.position="bottom") + scale_colour_manual("Wave",values=c("red","orange", "green", "blue"))
p <- p + theme(axis.text=element_text(size=14),
               axis.title=element_text(size=14))
print(p)
dev.off()

data <- subset(data, AF <= .45)
data.wave1 <- subset(data, wave == 1)
data.wave2 <- subset(data, wave == 2)
data.wave3 <- subset(data, wave == 3)
data.wave4 <- subset(data, wave == 4)

wave1_obs <- length(unique(data.wave1$participant_id))
wave2_obs <- length(unique(data.wave2$participant_id))
wave3_obs <- length(unique(data.wave3$participant_id))
wave4_obs <- length(unique(data.wave4$participant_id))

wave1_norm <- 34/wave1_obs
wave2_norm <- 34/wave2_obs
wave3_norm <- 34/wave3_obs
wave4_norm <- 34/wave4_obs

breaks <- c(0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45)
tags <- c("0.01","0.02","0.03","0.04","0.05","0.06","0.07","0.08","0.09","0.1","0.11","0.12","0.13","0.14","0.15","0.16","0.17","0.18","0.19","0.2","0.21","0.22","0.23","0.24","0.25","0.26","0.27","0.28","0.29","0.3","0.31","0.32","0.33","0.34","0.35","0.36","0.37","0.38","0.39","0.4","0.41","0.42","0.43","0.44","0.45")
group_tags1 <- cut(data.wave1$AF,
                  breaks=breaks, 
                  include.lowest=TRUE, 
                  right=FALSE, 
                  labels=tags)

group_tags2 <- cut(data.wave2$AF,
                   breaks=breaks, 
                   include.lowest=TRUE, 
                   right=FALSE, 
                   labels=tags)

group_tags3 <- cut(data.wave3$AF,
                   breaks=breaks, 
                   include.lowest=TRUE, 
                   right=FALSE, 
                   labels=tags)

group_tags4 <- cut(data.wave4$AF,
                   breaks=breaks, 
                   include.lowest=TRUE, 
                   right=FALSE, 
                   labels=tags)



summary(group_tags)
as_tibble(group_tags)

x <- rbind(table(group_tags1), table(group_tags2), table(group_tags3), table(group_tags4))
x <- as.data.frame(t(x))
colnames(x) <- c("wave1", "wave2", "wave3", "wave4")

x$wave1 <- x$wave1*wave1_norm
x$wave2 <- x$wave2*wave2_norm
x$wave3 <- x$wave3*wave3_norm
x$wave4 <- x$wave4*wave4_norm

x_labels <- c("","","","","","","","","","0.1","","","","","","","","","","0.2","","","","","","","","","","0.3","","","","","","","","","","0.4","","","","","")
y <- melt(setDT(x, keep.rownames = TRUE), "rn")

output_filename <- "LBC_ARCHER.GrantApplication.VAF_NormalisedHistogram.png"
png(output_filename, width = 1000, height = 1100, res = 200, pointsize =5)
p <- ggplot(y, aes(fill=variable, y=value, x=rn)) + 
  geom_bar(position="stack", stat="identity")
p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + ylab("Cumulative Count Normalised by Observations")
p <- p + xlab("Variant Allele Fraction")
p <- p + theme(axis.text.x = element_text(color="black", size=12),
               axis.text.y = element_text(color="black", size=12))
p <- p + theme(legend.title=element_text(), legend.position="bottom") + scale_fill_manual("",values=c("red","orange", "green", "blue"))
p <- p + theme(axis.text=element_text(size=14),
               axis.title=element_text(size=14))
p <- p + scale_x_discrete(labels = x_labels)
print(p)
dev.off()



