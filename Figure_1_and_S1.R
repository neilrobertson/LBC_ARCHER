library(ggplot2)
library(data.table)
library(dplyr)
library(dplyr)
library(reshape)
library(pheatmap)
library(maftools)
library(stringr)
library(ggsci)
library(Hmisc)

### Set Key Parameters/ Read in Data / Filter ###
setwd("/mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq/analysis/MS_Figures/Fig1_24May21")
input_data_prefix <- "LBC_ARCHER.2PCT_VAF.Mar21."
output_prefix <- "LBC_ARCHER.2PC_VAF.MS_Fig1_24May21."
data_origin_dir <- "/mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq/processed_data/Processed_0.5PCT_VAF_Mar2021/"
waves <- c(1,2,3,4,5)



study.info <- read.table("/mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq/ARCHER_PanelSeq_Metadata_SEH_deduped.All.tsv", header = T)
study.info <- subset(study.info, study.info$Original_Filename != "LBC0484V_07505000140_S46_R1_001.fastq.gz")
describe(study.info)
dim(study.info)


# https://stackoverflow.com/questions/9617348/reshape-three-column-data-frame-to-matrix-long-to-wide-format
data.non_syn <- read.table(paste(data_origin_dir, input_data_prefix, "non-synonymous.tsv", sep = ""), sep = "\t", header = T)
colnames(data.non_syn)[33] <- "AUAO" # Non unique column field (UAO) - taken from ARCHER sequencing information and probably not required, may remove in processing
data.non_syn$gene_key <- paste(data.non_syn$PreferredSymbol, data.non_syn$base_substitution, data.non_syn$participant_id, sep=" ")
data.non_syn$subject_wave <- paste(data.non_syn$participant_id, data.non_syn$wave, sep = "-")
data.non_syn$event_index <- paste(data.non_syn$participant_id, data.non_syn$base_substitution, sep = "_")
dim(data.non_syn)

complex.waves <- c("1\n~70 years\n~79 years","2\n~73 years\n~82 years","3\n~76 years\n~85 years","4\n~79 years\n~88 years","5\n~82 years\n~91 years")
names(complex.waves) <- waves
data.non_syn$complex_waves <- complex.waves[data.non_syn$wave]

#gg_color_hue <- function(n) {
#  hues = seq(15, 375, length = n + 1)
#  hcl(h = hues, l = 65, c = 100)[1:n]
#}
#n = 12
#cols = gg_color_hue(n)



table(data.non_syn$event_index)
data.unique_events <- data.non_syn[!duplicated(data.non_syn[,c('event_index')]), ]
table(data.unique_events$event_index)
dim(data.unique_events)

participant.events <- as.data.frame(table(data.unique_events$participant_id))
dim(participant.events)
mean(participant.events$Freq)
max(participant.events$Freq)

dim(subset(data.unique_events, data.unique_events$PreferredSymbol == "DNMT3A"))
length(unique(subset(data.unique_events, data.unique_events$PreferredSymbol == "DNMT3A")[["participant_id"]]))

dim(subset(data.unique_events, data.unique_events$PreferredSymbol == "TET2"))
length(unique(subset(data.unique_events, data.unique_events$PreferredSymbol == "TET2")[["participant_id"]]))

dim(subset(data.unique_events, data.unique_events$PreferredSymbol == "JAK2"))
length(unique(subset(data.unique_events, data.unique_events$PreferredSymbol == "JAK2")[["participant_id"]]))

dim(subset(data.unique_events, data.unique_events$PreferredSymbol == "NOTCH1"))
length(unique(subset(data.unique_events, data.unique_events$PreferredSymbol == "NOTCH1")[["participant_id"]]))

dim(subset(data.unique_events, data.unique_events$PreferredSymbol == "U2AF2"))
length(unique(subset(data.unique_events, data.unique_events$PreferredSymbol == "U2AF2")[["participant_id"]]))



### Variant Counts - All Unique Events ###
order <- as.data.frame(table(data.unique_events$clear_classification))
order <- order[ order(order$Freq, decreasing = F) , ]
p <- ggplot(data.unique_events) + geom_bar(aes(x = clear_classification, fill = simple_classification))
p <-  p + theme_classic() 
#+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + scale_x_discrete(limits=order$Var1)
p <- p + scale_y_continuous(expand = c(0,0),
                            limits=c(0,100),breaks=c(0,20,40,60,80,100))
p <- p + xlab("")
p <- p + ylab("Count")
p <- p + theme(legend.position = "none")
p <- p + theme(axis.text.x = element_text(color="black", size=10),
               axis.text.y = element_text(color="black", size=9))
p <- p + theme(axis.text=element_text(size=10),
               axis.title=element_text(size=10))
p <- p + coord_flip()
print(p)
p <- p + scale_fill_manual(values=c(Flank3="#00BBDA",Flank5="#E18A00",UTR5="#BE9C00",Frame_Shift_Del="#00BDC2",Frame_Shift_Ins="#24B700",In_Frame_Ins="#00BE70",Missense_Mutation="#00C1AB",Nonsense_Mutation="#F8766D",Nonstop_Mutation="#8B93FF",Splice_Region="#D575FE",Splice_Site="#F962DD"))
print(p)

output_filename <- paste("./png/", output_prefix, "VariantClassification_Count.AllWaves.UniqueEvents.BarChart.png", sep = "")
png(output_filename, width = 1000, height = 1100, res = 200, pointsize =5)
print(p)
dev.off()

output_filename <- paste("./svg/", output_prefix, "VariantClassification_Count.AllWaves.UniqueEvents.BarChart.svg", sep = "")
svg(output_filename, width = 3.2, height = 2, pointsize =8)
par(cex=1)
print(p)
dev.off()



### Per Gene Count ###
order <- as.data.frame(table(data.unique_events$PreferredSymbol))
order <- order[ order(order$Freq, decreasing = T) , ]
p <- ggplot(data.unique_events) + geom_bar(aes(x = PreferredSymbol, group = PreferredSymbol), fill = "grey")
p <- p + theme_classic() 
#+ theme(panel.grid.major = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + scale_x_discrete(limits=order$Var1)
p <- p + scale_y_continuous(expand = c(0,0),
                            limits = c(0,60), breaks=c(0,10,20,30,40,50,60))
p <- p + xlab("")
p <- p + ylab("Count")
p <- p + theme(legend.position = "none")
p <- p + theme(axis.text.x = element_text(color="black", angle = 45, size=6, hjust=0.95, vjust  = .8),
               axis.text.y = element_text(color="black", size=10))
p <- p + theme(axis.text=element_text(size=11),
               axis.title=element_text(size=11))
#p <- p + coord_flip()
print(p)

output_filename <- paste("./png/", output_prefix, "GeneSymbol_Count.AllWaves.UniqueEvents.BarChart.png", sep = "")
png(output_filename, width = 1800, height = 1100, res = 200, pointsize =5)
print(p)
dev.off()

output_filename <- paste("./svg/", output_prefix, "GeneSymbol_Count.AllWaves.UniqueEvents.BarChart.svg", sep = "")
svg(output_filename, width = 4.7, height = 2.2, pointsize =8)
par(cex=1)
print(p)
dev.off()




### Key Genes Largest Clone Bubble Heatmap ###
data.non_syn.key <- data.non_syn[data.non_syn$PreferredSymbol %in% c("DNMT3A", "TET2", "NOTCH1", "JAK2", "U2AF2", "JAK3", "ASXL1") , ]
data.largest_per_gene_AF <- data.non_syn.key %>%
  group_by(PreferredSymbol, subject_wave) %>% slice(which.max(AF))

data.largest_per_gene_AF <- as.data.frame(data.largest_per_gene_AF)

data.anno <- data.largest_per_gene_AF[ , colnames(data.largest_per_gene_AF) %in% c("subject_wave", "wave", "participant_id")]
data.anno <- unique(data.anno)
rownames(data.anno) <- data.anno$subject_wave
data.anno <- data.anno[ , !(colnames(data.anno) %in% c("subject_wave"))]

for_order <- reshape2::dcast(as.data.frame(data.largest_per_gene_AF), PreferredSymbol~subject_wave, value.var="AF")
for_order[is.na(for_order)] = 0
rownames(for_order) <- for_order$PreferredSymbol
for_order <- for_order[ , !(colnames(for_order) %in% c("PreferredSymbol")) ]
for_order <- for_order[apply(for_order, 1, function(x) any(x != 0 | is.na(x))), ]
cluster_order <- pheatmap(for_order, annotation_col = data.anno, show_rownames = F)

output_filename <- paste("./png/", output_prefix, "Largest_Clone.AllWaves.SelectGenes.pheatmap.png", sep = "")
png(output_filename, width = 2400, height = 800, res = 300, pointsize =4)
print(cluster_order)
dev.off()

output_filename <- paste("./svg/", output_prefix, "Largest_Clone.AllWaves.SelectGenes.pheatmap.svg", sep = "")
svg(output_filename, width = 12, height = 5, pointsize =8)
par(cex=1)
print(cluster_order)
dev.off()

p <- ggplot(data.largest_per_gene_AF, aes(subject_wave, PreferredSymbol)) +
  geom_point(aes(size = AF, group = colour_code, colour = simple_classification), alpha=0.3) 
p <- p + scale_size(range = c(0,15)) 
p <- p + xlab("") + ylab("")
p <- p + theme_minimal() + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                 legend.text=element_text(size=11))

p <- p + theme(axis.text.x = element_blank(),
               axis.text.y = element_text(color="black", size=9))
p <- p + theme(axis.text=element_text(size=16),
               axis.title=element_text(size=16, face = "bold"))
p <- p + scale_y_discrete(limits=rev(cluster_order$tree_row$labels[cluster_order$tree_row$order]))
p <- p + scale_x_discrete(limits=(cluster_order$tree_col$labels[cluster_order$tree_col$order]))
p <- p + guides(colour = guide_legend(override.aes = list(size=3), title = c("Variant Classification", "VAF")), labels = data.largest_per_gene_AF$Variant_Classification)
print(p)
p <- p + scale_colour_manual(values=c(Flank3="#00BBDA",Flank5="#E18A00",UTR5="#BE9C00",Frame_Shift_Del="#00BDC2",Frame_Shift_Ins="#24B700",In_Frame_Ins="#00BE70",Missense_Mutation="#00C1AB",Nonsense_Mutation="#F8766D",Nonstop_Mutation="#8B93FF",Splice_Region="#D575FE",Splice_Site="#F962DD"))
print(p)

output_filename <- paste("./png/", output_prefix, "Largest_Clone.AllWaves.SelectGenes.Bubble.png", sep = "")
png(output_filename, width = 2400, height = 1500, res = 250, pointsize =4)
print(p)
dev.off()

output_filename <- paste("./svg/", output_prefix, "Largest_Clone.AllWaves.SelectGenes.Bubble.svg", sep = "")
svg(output_filename, width = 10, height = 3, pointsize =8)
par(cex=1)
print(p)
dev.off()

output_filename <- paste("./svg/", output_prefix, "Largest_Clone.AllWaves.SelectGenes.Bubble.NoLegend.svg", sep = "")
svg(output_filename, width = 8, height = 2.5, pointsize =8)
par(cex=1)
p <- p + theme(legend.position = "none")
print(p)
dev.off()




### Largest Clone Bubble Heatmap ###
data.largest_per_gene_AF <- data.non_syn %>%
  group_by(PreferredSymbol, subject_wave) %>% slice(which.max(AF))

data.largest_per_gene_AF <- as.data.frame(data.largest_per_gene_AF)

data.anno <- data.largest_per_gene_AF[ , colnames(data.largest_per_gene_AF) %in% c("subject_wave", "wave", "participant_id")]
data.anno <- unique(data.anno)
rownames(data.anno) <- data.anno$subject_wave
data.anno <- data.anno[ , !(colnames(data.anno) %in% c("subject_wave"))]

for_order <- reshape2::dcast(as.data.frame(data.largest_per_gene_AF), PreferredSymbol~subject_wave, value.var="AF")
for_order[is.na(for_order)] = 0
rownames(for_order) <- for_order$PreferredSymbol
for_order <- for_order[ , !(colnames(for_order) %in% c("PreferredSymbol")) ]
for_order <- for_order[apply(for_order, 1, function(x) any(x != 0 | is.na(x))), ]
cluster_order <- pheatmap(for_order, annotation_col = data.anno, show_rownames = F)

output_filename <- paste("./png/", output_prefix, "Largest_Clone.AllWaves.pheatmap.png", sep = "")
png(output_filename, width = 2400, height = 1500, res = 300, pointsize =4)
print(cluster_order)
dev.off()

output_filename <- paste("./svg/", output_prefix, "Largest_Clone.AllWaves.pheatmap.svg", sep = "")
svg(output_filename, width = 12, height = 8, pointsize =8)
par(cex=1)
print(cluster_order)
dev.off()

p <- ggplot(data.largest_per_gene_AF, aes(subject_wave, PreferredSymbol)) +
  geom_point(aes(size = AF, group = simple_classification, colour = simple_classification), alpha=0.3) 
p <- p + scale_size(range = c(0,15)) 
p <- p + xlab("") + ylab("")
p <- p + theme_minimal() + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                 legend.text=element_text(size=11))

p <- p + theme(axis.text.x = element_blank(),
               axis.text.y = element_text(color="black", size=9))
p <- p + theme(axis.text=element_text(size=16),
               axis.title=element_text(size=16, face = "bold"))
p <- p + scale_y_discrete(limits=rev(cluster_order$tree_row$labels[cluster_order$tree_row$order]))
p <- p + scale_x_discrete(limits=(cluster_order$tree_col$labels[cluster_order$tree_col$order]))
p <- p + guides(colour = guide_legend(override.aes = list(size=3), title = c("Variant Classification", "VAF")), labels = data.largest_per_gene_AF$Variant_Classification)
print(p)
p <- p + scale_colour_manual(values=c(Flank3="#00BBDA",Flank5="#E18A00",UTR5="#BE9C00",Frame_Shift_Del="#00BDC2",Frame_Shift_Ins="#24B700",In_Frame_Ins="#00BE70",Missense_Mutation="#00C1AB",Nonsense_Mutation="#F8766D",Nonstop_Mutation="#8B93FF",Splice_Region="#D575FE",Splice_Site="#F962DD"))
print(p)


output_filename <- paste("./png/", output_prefix, "Largest_Clone.AllWaves.Bubble.png", sep = "")
png(output_filename, width = 2400, height = 1500, res = 250, pointsize =4)
print(p)
dev.off()

output_filename <- paste("./svg/", output_prefix, "Largest_Clone.AllWaves.Bubble.svg", sep = "")
svg(output_filename, width = 10, height = 5, pointsize =8)
par(cex=1)
print(p)
dev.off()

output_filename <- paste("./svg/", output_prefix, "Largest_Clone.AllWaves.Bubble.NoLegend.svg", sep = "")
svg(output_filename, width = 8, height = 5, pointsize =8)
par(cex=1)
p <- p + theme(legend.position = "none")
print(p)
dev.off()



### Interesting Gene Trajectories ###
genes_of_interest <- c("DNMT3A", "TET2", "JAK2", "ASXL1", "SF3B1")
for (id in genes_of_interest) {
  print (id)
  data.non_syn.id <- subset(data.non_syn, PreferredSymbol == id)
  
  print(dim(data.non_syn.id))
  
  p <- ggplot() + geom_line(data = data.non_syn.id, aes(x=wave, y=AF, colour = simple_classification, group = gene_key), size = 0.4, alpha = 1) + 
    geom_point(data = data.non_syn.id, aes(x=wave, y=AF, colour = simple_classification, group = gene_key))
  p <- p + theme_classic() + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                   axis.text=element_text(size=8),
                                   axis.title=element_text(size=8))
  #panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  p <- p + theme(axis.text.x = element_text(color="black", size=6.5),
                 axis.text.y = element_text(color="black", size=7))
  p <- p + ylab("VAF")
  p <- p + scale_x_continuous(name="",limits=c(1,5),breaks=c(1,2,3,4,5), labels = c("1\n~70 years\n~79 years","2\n~73 years\n~82 years","3\n~76 years\n~85 years","4\n~79 years\n~88 years","5\n~82 years\n~91 years"))
  #p <- p + ggtitle(id)
  #p <- p + theme(plot.title = element_text(vjust = - 4, size = 12))
  p <- p + theme(legend.position = "none")
  print(p)
  p <- p + scale_colour_manual(values=c(Flank3="#00BBDA",Flank5="#E18A00",UTR5="#BE9C00",Frame_Shift_Del="#00BDC2",Frame_Shift_Ins="#24B700",In_Frame_Ins="#00BE70",Missense_Mutation="#00C1AB",Nonsense_Mutation="#F8766D",Nonstop_Mutation="#8B93FF",Splice_Region="#D575FE",Splice_Site="#F962DD"))
  print(p)


  output_filename <- paste("./png/", output_prefix, id, ".GeneTimeCourse.OfInterest.png", sep = "")
  png(output_filename, width = 1000, height = 1100, res = 200, pointsize =5)
  print(p)
  dev.off()
  
  output_filename <- paste("./svg/", output_prefix, id, ".GeneTimeCourse.OfInterest.svg", sep = "")
  svg(output_filename, width = 2.8, height = 2.2, pointsize =8)
  par(cex=1)
  print(p)
  dev.off()
}






### Jitter VAF vs Gene ###
data.wave1 <- subset(data.non_syn, wave == 1)
newcsv <- data.wave1 %>%
  group_by(PreferredSymbol) %>%
  summarise(
    Total_AF = sum(AF)
  )
dim(data.non_syn)
newcsv <- as.data.frame(newcsv)
order <- newcsv[ order(-newcsv$Total_AF) , ]

p <- ggplot(data.wave1) + 
  geom_boxplot( aes(x = PreferredSymbol, y = AF), outlier.shape = NA) + 
  geom_jitter(aes(x = PreferredSymbol, y = AF),  position=position_jitter(0.2), fill="grey30", alpha=.4, shape = 19, size = 1.5)
p <- p + xlab("")
p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                 legend.position=c(.84,.75), legend.text=element_text(size=11))
p <- p + theme(axis.text.x = element_text(color="black", size=7, angle = 45, hjust = 1, face = "bold"),
               axis.text.y = element_text(color="black", size=7))
p <- p + theme(axis.text=element_text(size=11),
               axis.title=element_text(size=11, face = "bold"))
p <- p + scale_x_discrete(limits=order$PreferredSymbol)
p <- p + scale_y_continuous(name="VAF (Wave 1)",limits=c(0,.6),breaks=c(0,.1,.2,.3,.4,.5,.6))

output_filename <- paste("./png/", output_prefix, "VAF_ByGene_Jitter.png", sep = "")
png(output_filename, width = 2000, height = 1000, res = 150, pointsize =5)
print(p)
dev.off()

p <- ggplot(data.wave1 ) + 
  geom_boxplot( aes(x = PreferredSymbol, y = AF), outlier.shape = NA) + 
  geom_jitter(aes(x = PreferredSymbol, y = AF, fill = simple_classification),  position=position_jitter(0.2), alpha=.5, shape = 21, size = 1.2)
p <- p + xlab("")
p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x = element_text(color="black", size=6, angle = 45, hjust = 1),
               axis.text.y = element_text(color="black", size=8))
p <- p + theme(axis.text=element_text(size=10),
               axis.title=element_text(size=10))
p <- p + theme(legend.position = "none")
p <- p + scale_x_discrete(limits=order$PreferredSymbol)
p <- p + scale_y_continuous(name="VAF",limits=c(0,.5),breaks=c(0,.1,.2,.3,.4,.5))
print(p)
p <- p + scale_fill_manual(values=c(Flank3="#00BBDA",Flank5="#E18A00",UTR5="#BE9C00",Frame_Shift_Del="#00BDC2",Frame_Shift_Ins="#24B700",In_Frame_Ins="#00BE70",Missense_Mutation="#00C1AB",Nonsense_Mutation="#F8766D",Nonstop_Mutation="#8B93FF",Splice_Region="#D575FE",Splice_Site="#F962DD"))
print(p)

output_filename <- paste("./png/", output_prefix, "VAF_ByGene_Jitter.Coloured.png", sep = "")
png(output_filename, width = 2000, height = 1000, res = 150, pointsize =5)
print(p)
dev.off()

output_filename <- paste("./svg/", output_prefix, "VAF_ByGene_Jitter.Coloured.svg", sep = "")
svg(output_filename, width = 4.7, height = 3.5, pointsize =8)
par(cex=1)
print(p)
dev.off()



### Distribution of TimePoints Heatmap ###
study.info <- read.table("/mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq/ARCHER_PanelSeq_Metadata_SEH_deduped.All.tsv", header = T)
describe(study.info)
table(study.info$Cohort)
study.info <- subset(study.info, Cohort != "CONTROL")
study.info$Cohort <- ifelse(study.info$Cohort == "LBC21", "LBC1921", ifelse(study.info$Cohort == "LBC36", "LBC1936", ifelse(study.info$Cohort == "LBC1936", "LBC1936", ifelse(study.info$Cohort == "LBC1921", "LBC1921", "NA"))))
study.info <- subset(study.info, study.info$Original_Filename != "LBC0484V_07505000140_S46_R1_001.fastq.gz")
table(study.info$Cohort)


study.info.matrix <- study.info %>%
  group_by(Particiapant_ID, Wave) %>% select(Batch_Sample_Number)
for_order <- reshape2::dcast(as.data.frame(study.info.matrix), Particiapant_ID~Wave, value.var="Batch_Sample_Number")
rownames(for_order) <- for_order$Particiapant_ID
for_order <- for_order[ , !(colnames(for_order) %in% c("Particiapant_ID"))]

for_order[!is.na(for_order)] = 1
for_order[is.na(for_order)] = 0
for_order[] <- lapply(for_order, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
})

p <- pheatmap(as.matrix(for_order))
output_filename <- paste("./png/", output_prefix, "Participant-Wave.DistributionHeatmap.png", sep = "")
png(output_filename, width = 2000, height = 1000, res = 150, pointsize =5)
print(p)
dev.off()


mean(data.non_syn$AO)
median(data.non_syn$AO)
### Quality Metrics Scatter Plot ###
p <- ggplot(data = data.non_syn, aes(x = AO, y = UAO)) + geom_point(data = data.non_syn, aes(x=AO, y=UAO, fill = simple_classification, size = AF), alpha = 0.6, shape = 21) 
p <- p + geom_rug(data=data.non_syn, aes(col=simple_classification),alpha=0.6, size=0.6)
p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                 axis.text=element_text(size=10),
                                 axis.title=element_text(size=10))
p <- p + theme(axis.text.x = element_text(color="black", size=8),
               axis.text.y = element_text(color="black", size=8))
p <- p + ylab("UAO") + xlab("AU")
p <- p + scale_y_continuous(trans = "log10", limits = c(1,1000), breaks = c(1,10,100,500,1000)) 
p <- p + scale_x_continuous(trans = "log10", limits = c(1,5000), breaks = c(1,10,100,500,1000,5000))
p <- p + geom_hline(yintercept = 3, colour = "red", linetype = "dotted")
p <- p + geom_vline(xintercept = 5, colour = "red", linetype = "dotted")
print(p)
p <- p + scale_colour_manual(values=c(Flank3="#00BBDA",Flank5="#E18A00",UTR5="#BE9C00",Frame_Shift_Del="#00BDC2",Frame_Shift_Ins="#24B700",In_Frame_Ins="#00BE70",Missense_Mutation="#00C1AB",Nonsense_Mutation="#F8766D",Nonstop_Mutation="#8B93FF",Splice_Region="#D575FE",Splice_Site="#F962DD"))
p <- p + scale_fill_manual(values=c(Flank3="#00BBDA",Flank5="#E18A00",UTR5="#BE9C00",Frame_Shift_Del="#00BDC2",Frame_Shift_Ins="#24B700",In_Frame_Ins="#00BE70",Missense_Mutation="#00C1AB",Nonsense_Mutation="#F8766D",Nonstop_Mutation="#8B93FF",Splice_Region="#D575FE",Splice_Site="#F962DD"))
print(p)



output_filename <- paste("./png/", output_prefix, "Coverage_Quality.Scatter.png", sep = "")
png(output_filename, width = 800, height = 600, res = 250, pointsize =4)
print(p)
dev.off()

output_filename <- paste("./svg/", output_prefix, "Coverage_Quality.Scatter.svg", sep = "")
svg(output_filename, width = 4, height = 4, pointsize =8)
par(cex=1)
print(p)
dev.off()

output_filename <- paste("./svg/", output_prefix, "Coverage_Quality.Scatter.NoLegend.svg", sep = "")
svg(output_filename, width = 4, height = 4, pointsize =8)
par(cex=1)
p <- p + theme(legend.position = "none")
print(p)
dev.off()




# https://bioconductor.org/packages/release/bioc/vignettes/trackViewer/inst/doc/trackViewer.html#Introduction
library(Gviz)
library(rtracklayer)
library(trackViewer)
library(g3viz)
library(stringr)

lol.colour.keys <- data.non_syn[ , colnames(data.non_syn) %in% c("protein_substitution", "colour_code")]
lol.colour.keys <- distinct(lol.colour.keys)
lol.colour.keys <- as.data.table(lol.colour.keys[ lol.colour.keys$protein_substitution != "" , ])

all.p <- lol.colour.keys$protein_substitution
all.c <- lol.colour.keys$colour_code
names(all.c) <- all.p

### DNMT3A ###
#https://www.uniprot.org/uniprot/Q9Y6K1#function
data.non_syn$prot_subject_key <- paste(data.non_syn$participant_id, data.non_syn$protein_substitution, sep = "_")
dnmt3a <- data.non_syn[ data.non_syn$PreferredSymbol == "DNMT3A" , ]
order <- as.data.frame(table(unique(dnmt3a$prot_subject_key)))
order$p.key <- gsub("^.*_", "", order$Var1)
order <- order[ order$p.key != "" , ]
order <- as.data.frame(table(order$p.key))
order
order <- order[ order(order$Freq, decreasing = T) , ]
order
order$pos <- substr(order$Var1, 6, 8)
order$pos <- gsub("([0-9])([a-zA-Z])","\\1", order$pos)
order
order$cols <- all.c[order$Var1]

SNP <- as.numeric(as.character(order$pos))
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=order$Var1))
sample.gr$score <- order$Freq
sample.gr$color <- order$cols
sample.gr$SNPsideID <- sample(c("top", "bottom"), 
                              length(sample.gr),
                              replace=TRUE)
sample.gr$label.parameter.rot <- 60

features <- GRanges("chr1", IRanges(c(1, 292, 350, 482, 614, 634, 912), 
                                    width=c(292, 58, 132, 132, 20, 278, 0),
                                    names=c("","PWWP", "", "ADD", "", "SAM-Dependent MTase C5-type", "")))
features$height <- list(unit(1, "mm"), 
                        unit(3, "mm"), 
                        unit(1, "mm"),
                        unit(3, "mm"),
                        unit(1, "mm"),
                        unit(3, "mm"),
                        unit(1, "mm"))
features$fill <- c("#000000", "#FF8833", "#000000", "#51C6E6", "#000000", "#DFA32D", "#000000")

xaxis <- c(1, 200, 400, 600, 800, 912)
yaxis <- c(0, 1, 2, 3)
p <- lolliplot(sample.gr, features, xaxis = xaxis, yaxis = yaxis, type = "circle", cex=.6)
print(p)

output_filename <- paste("./svg/", output_prefix, "Lollipop.DNMT3A.svg", sep = "")
svg(output_filename, width = 4.5, height = 3.2, pointsize =8)
par(cex=1)
lolliplot(sample.gr, features, xaxis = xaxis, ylab = "", type = "circle", cex=.6)
dev.off()





### TET2 ###
#https://www.uniprot.org/uniprot/Q6N021#family_and_domains
data.non_syn$prot_subject_key <- paste(data.non_syn$participant_id, data.non_syn$protein_substitution, sep = "_")
tet2 <- data.non_syn[ data.non_syn$PreferredSymbol == "TET2" , ]
order <- as.data.frame(table(unique(tet2$prot_subject_key)))
order$p.key <- gsub("^.*_", "", order$Var1)
order <- order[ order$p.key != "" , ]
order <- as.data.frame(table(order$p.key))
order
order <- order[ order(order$Freq, decreasing = T) , ]
order
order$pos <- substr(order$Var1, 6, 8)
order$pos <- gsub("([0-9])([a-zA-Z])","\\1", order$pos)
order
order$cols <- all.c[order$Var1]



SNP <- as.numeric(as.character(order$pos))
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=order$Var1))
sample.gr$score <- order$Freq
sample.gr$color <- order$cols
sample.gr$SNPsideID <- sample(c("top", "bottom"), 
                              length(sample.gr),
                              replace=TRUE)
sample.gr$label.parameter.rot <- 60

features <- GRanges("chr1", IRanges(c(1, 1290, 1303, 1896, 1898, 1902, 1904), 
                                    width=c(1290, 13, 593, 2, 4, 2, 100),
                                    names=c("","Interaction w/ DNA", "", "2-Oxyglutarate Binding", "", "Substrate Binding", "")))
features$height <- list(unit(1, "mm"), 
                        unit(3, "mm"), 
                        unit(1, "mm"),
                        unit(3, "mm"),
                        unit(1, "mm"),
                        unit(3, "mm"),
                        unit(1, "mm"))
features$fill <- c("#000000", "#FF8833", "#000000", "#51C6E6", "#000000", "#DFA32D", "#000000")

xaxis <- c(1, 200, 400, 600, 800, 1000,1200,1400,1600,1800,2000)
yaxis <- c(0, 1, 2, 3)
p <- lolliplot(sample.gr, features, xaxis = xaxis, yaxis = yaxis, type = "circle", cex=.6)
print(p)

output_filename <- paste("./svg/", output_prefix, "Lollipop.TET2.svg", sep = "")
svg(output_filename, width = 4.5, height = 3.2, pointsize =8)
par(cex=1)
lolliplot(sample.gr, features, xaxis = xaxis, ylab = "", type = "circle", cex=.6)
dev.off()




### JAK2 ###
#https://www.uniprot.org/uniprot/O60674#family_and_domains
data.non_syn$prot_subject_key <- paste(data.non_syn$participant_id, data.non_syn$protein_substitution, sep = "_")
jak2 <- data.non_syn[ data.non_syn$PreferredSymbol == "JAK2" , ]
order <- as.data.frame(table(unique(jak2$prot_subject_key)))
order$p.key <- gsub("^.*_", "", order$Var1)
order <- order[ order$p.key != "" , ]
order <- as.data.frame(table(order$p.key))
order
order <- order[ order(order$Freq, decreasing = T) , ]
order
order$pos <- substr(order$Var1, 6, 8)
order$pos <- gsub("([0-9])([a-zA-Z])","\\1", order$pos)
order
order$cols <- all.c[order$Var1]

SNP <- as.numeric(as.character(order$pos))
sample.gr <- GRanges("chr1", IRanges(SNP, width=1, names=order$Var1))
sample.gr$score <- order$Freq
sample.gr$color <- order$cols
sample.gr$SNPsideID <- sample(c("top", "bottom"), 
                              length(sample.gr),
                              replace=TRUE)
sample.gr$label.parameter.rot <- 60

features <- GRanges("chr1", IRanges(c(1, 37, 380, 401, 482, 545, 809, 849, 1124), 
                                    width=c(37, 343, 21, 81, 63, 264, 40, 275, 8),
                                    names=c("","FERM", "", "SH2; atypical", "", "Protein Kinase 1", "", "Protein Kinase 2", "")))
features$height <- list(unit(1, "mm"), 
                        unit(3, "mm"), 
                        unit(1, "mm"),
                        unit(3, "mm"),
                        unit(1, "mm"),
                        unit(3, "mm"),
                        unit(1, "mm"),
                        unit(3, "mm"),
                        unit(1, "mm"))
features$fill <- c("#000000", "#FF8833", "#000000", "#51C6E6", "#000000", "#DFA32D", "#000000", "#DFA32D", "#000000")

xaxis <- c(1, 200, 400, 600, 800, 1000,1200)
p <- lolliplot(sample.gr, features, xaxis = xaxis, type = "circle", cex=.6)
print(p)

output_filename <- paste("./svg/", output_prefix, "Lollipop.JAK2.svg", sep = "")
svg(output_filename, width = 4.5, height = 3.2, pointsize =8)
par(cex=1)
lolliplot(sample.gr, features, xaxis = xaxis, ylab = "", type = "circle", cex=.6)
dev.off()





### Variant Counts ###
data.wave1 <- subset(data.non_syn, wave == 1)
order <- as.data.frame(table(data.wave1$clear_classification))
order <- order[ order(order$Freq, decreasing = F) , ]
p <- ggplot(data.wave1) + geom_bar(aes(x = clear_classification, fill = colour_code))
p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + scale_x_discrete(limits=order$Var1)
p <- p + scale_y_continuous(expand = c(0,0),
                            limits=c(0,90),breaks=c(0,15,30,45,60,75,90))
p <- p + xlab("")
p <- p + ylab("Count")
p <- p + theme(legend.position = "none")
p <- p + theme(axis.text.x = element_text(color="black", size=10),
               axis.text.y = element_text(color="black", size=9))
p <- p + theme(axis.text=element_text(size=10),
               axis.title=element_text(size=10))
p <- p + coord_flip()
print(p)

output_filename <- paste("./png/", output_prefix, "VariantClassification_Count.Wave1.BarChart.png", sep = "")
png(output_filename, width = 1000, height = 1100, res = 200, pointsize =5)
print(p)
dev.off()

output_filename <- paste("./svg/", output_prefix, "VariantClassification_Count.Wave1.BarChart.svg", sep = "")
svg(output_filename, width = 3.2, height = 2, pointsize =8)
par(cex=1)
print(p)
dev.off()



### Per Gene Count ###
data.wave1 <- subset(data.non_syn, wave == 1)
order <- as.data.frame(table(data.wave1$PreferredSymbol))
order <- order[ order(order$Freq, decreasing = T) , ]
order <- order[ !(order$Var1 %in% c("BRAF", "MYD88")) , ]
p <- ggplot(data.wave1) + geom_bar(aes(x = PreferredSymbol, group = PreferredSymbol), fill = "mediumturquoise")
p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + scale_x_discrete(limits=order$Var1)
p <- p + scale_y_continuous(expand = c(0,0),
                            limits = c(0,52))
p <- p + xlab("")
p <- p + ylab("Count")
p <- p + theme(legend.position = "none")
p <- p + theme(axis.text.x = element_text(color="black", size=6, angle = 45, hjust=1),
               axis.text.y = element_text(color="black", size=10))
p <- p + theme(axis.text=element_text(size=11),
               axis.title=element_text(size=11))
#p <- p + coord_flip()
print(p)

output_filename <- paste("./png/", output_prefix, "GeneSymbol_Count.Wave1.BarChart.png", sep = "")
png(output_filename, width = 1800, height = 1100, res = 200, pointsize =5)
print(p)
dev.off()

output_filename <- paste("./svg/", output_prefix, "GeneSymbol_Count.Wave1.BarChart.svg", sep = "")
svg(output_filename, width = 4.7, height = 2.2, pointsize =8)
par(cex=1)
print(p)
dev.off()




setwd("/mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq")
meta <- read.table("ARCHER_PanelSeq_Metadata_SEH_deduped.All.tsv", header = T, sep = "\t")
meta <- subset(meta, meta$Wave != "CONTROL")
dim(meta)

cohort <- readRDS("/mnt/tchandra-lab/Neil/projects/LBC_ARCH_meth_clock/analysis/Figures_Resub/CurrentBiology_ARCH_revision_data_20May2019.rds")
dim(cohort)
length(unique(cohort$ID_raw))
ids <- cohort$ID_raw
sex <- cohort$sex
names(sex) <- ids

meta.1 <- subset(meta, meta$Wave == 1)
table(sex[meta.1$Particiapant_ID])

meta.2 <- subset(meta, meta$Wave == 2)
table(sex[meta.2$Particiapant_ID])

meta.3 <- subset(meta, meta$Wave == 3)
table(sex[meta.3$Particiapant_ID])

meta.4 <- subset(meta, meta$Wave == 4)
table(sex[meta.4$Particiapant_ID])

meta.5 <- subset(meta, meta$Wave == 5)
table(sex[meta.5$Particiapant_ID])

meta$Cohort <- ifelse(meta$Cohort == "LBC21", "LBC21", ifelse( meta$Cohort == "LBC36", "LBC36", ifelse( meta$Cohort == "LBC1936", "LBC36", ifelse( meta$Cohort == "LBC1921", "LBC21", "NA"))))
table( meta$Cohort)
meta.21 <- subset(meta, meta$Cohort == "LBC21")
meta.36 <- subset(meta, meta$Cohort == "LBC36")
dim(meta.36)
dim(meta.21)
length(unique(meta.36$Particiapant_ID))
length(unique(meta.21$Particiapant_ID))

write.table(unique(meta.36$Particiapant_ID), "LBC36.UniqueIDs.tsv", sep = "\t", row.names = F, quote = F)
write.table(unique(meta.21$Particiapant_ID), "LBC21.UniqueIDs.tsv", sep = "\t", row.names = F, quote = F)

id.map <- read.table("Manuscript.LBC_IDs.ObfusticationMap.tsv", header = F)
id.map["LBC0046H"]
cur <- id.map$V2
names(cur) <- id.map$V1
id.map = cur
id.map

qc <- read.table("ARCHER_PanelSeq_Metadata.QC_info.tsv", header = T)
qc$new_ids <- id.map[qc$Particiapant_ID]

write.table(qc, "ARCHER_PanelSeq_Metadata.QC_info.AlteredIDs.tsv", row.names = F, sep = "\t")
