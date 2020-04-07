#!/usr/bin/env R

library(ggplot2)
library(data.table)
library(dplyr)
library(reshape)


setwd("/mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq/ARCHER_processed/processed_variants")
files <- list.files(pattern = ".tsv")
temp <- lapply(files, fread, sep="\t", header = T)
data <- rbindlist( temp )
data <- as.data.frame(data)
data$key <- paste(data$PreferredSymbol, data$base_substitution, sep=" ")
colnames(data)[33] <- "AUAO" # Non unique column field (UAO) - taken from ARCHER sequencing information and probably not required, may remove in processing
dim(data)

setwd("/mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq/analysis/CloneSizeDistributions")
variant <- "RUNX1 c.1265_1266delinsCC"
for (variant in (unique(data$key))) {
	print (variant)
	variant.data <- subset(data, key == variant)
	variant.data$key2 <- paste(variant.data$PreferredSymbol, variant.data$base_substitution, sep="_")

	output_filename <- paste("LBC_ARCHER.", variant.data$key2, ".CloneSizeDistributions.VAF_vs_Age.png", sep = "")
	png(output_filename, width = 1000, height = 1000, res = 200, pointsize =5)
	p <- ggplot(data=variant.data, aes(x=wave, y=AF, group = participant_id))
	p <- p + geom_line(data = variant.data, aes(x=wave, y=AF, group = participant_id), size = 1, col = "grey")
	p <- p + geom_point(data = variant.data, aes(x=wave, y=AF), colour = "grey")
	p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
									panel.background = element_blank(), axis.line = element_line(colour = "black"))
	p <- p + ylab("Variant Allele Fraction")
	p <- p + xlab("Wave")
	p <- p + ggtitle(variant)
	p <- p + scale_x_continuous(name="Wave",breaks=c(1,2,3,4), labels = c("1","2","3","4"))
	p <- p + scale_y_continuous(name="Variant Allele Fraction",limits=c(0,0.55),breaks=c(0,.1,.2,.3,.4,.5))
	p <- p + theme(axis.text.x = element_text(color="black", size=12),
					axis.text.y = element_text(color="black", size=12))
	p <- p + theme(axis.text=element_text(size=14),
					axis.title=element_text(size=14))
	print(p)
	dev.off()


	output_filename <-  paste("LBC_ARCHER.", variant.data$key2, ".CloneSizeDistributions.VAF_vs_Particpant.png", sep = "")
	png(output_filename, width = 1100, height = 1000, res = 200, pointsize =5)
	p <- ggplot(data=variant.data, aes(x=participant_id, y=AF, group = wave))
	p <- p + geom_line(data = variant.data, aes(x=participant_id, y=AF, group = wave, col = as.factor(wave)), size = 1)
	p <- p + geom_point(data = variant.data, aes(x=participant_id, y=AF, col = as.factor(wave)))
	p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
									panel.background = element_blank(), axis.line = element_line(colour = "black"))
	p <- p + ylab("Variant Allele Fraction")
	p <- p + xlab("Particicipant")
	p <- p + ggtitle(variant)
	p <- p + scale_y_continuous(name="Variant Allele Fraction",limits=c(0,0.55),breaks=c(0,.1,.2,.3,.4,.5))
	p <- p + theme(legend.title=element_blank(), legend.position="bottom")
	p <- p + theme(axis.text.x = element_text(color="black", size=12, angle = 45, hjust = 1),
					axis.text.y = element_text(color="black", size=12))
	p <- p + theme(axis.text=element_text(size=14),
					axis.title=element_text(size=14))
	print(p)
	dev.off()
}

