###
#cd /mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq
#python ./code/parse_ARCHER_syn.py --manifest ARCHER_PanelSeq_Metadata_SEH_deduped.All.tsv --output_dir ./processed_data/Processed_0.5PCT_VAF_Mar2021/


### Set Key Parameters/ Read in Data / Filter ###
setwd("/mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq/processed_data/Processed_0.5PCT_VAF_Mar2021/")
data_origin_dir <- "/mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq/processed_data/Processed_0.5PCT_VAF_Mar2021"
waves <- c(1,2,3,4,5)

library(ggplot2)
library(data.table)
library(dplyr)
library(reshape)
library(pheatmap)
library(maftools)
library(stringr)
library(ggsci)

CONSEQUENCE_FILTER <- c("5_prime_UTR_variant", "coding_sequence_variant", "feature_elongation", "feature_truncation", "frameshift_variant", "incomplete_terminal_codon_variant", "inframe_deletion", "inframe_insertion", "missense_variant", "protein_altering_variant", "start_lost", "stop_gained", "stop_lost", "transcript_ablation", "transcript_amplification", "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant")

ERROR_CALLS_BASESUB <- read.table("/mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq/resources/LBC_ARCHER.ErrorCalls.Mar21.txt", header = F, skip = 1, sep = "\t")
PROBABLE_GERMLINE <- read.table("/mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq/resources/LBC_ARCHER.GermlineCalls.Mar21.txt", header = F, skip = 1, sep = "\t")


vep.vc = c("splice_acceptor_variant", "splice_donor_variant", "transcript_ablation",
           "exon_loss_variant", "stop_gained", "stop_lost", "initiator_codon_variant",
           "start_lost", "missense_variant", "coding_sequence_variant",
           "conservative_missense_variant", "rare_amino_acid_variant", "transcript_amplification",
           "intron_variant", "INTRAGENIC", "intragenic_variant", "splice_region_variant",
           "incomplete_terminal_codon_variant", "synonymous_variant", "stop_retained_variant",
           "NMD_transcript_variant", "mature_miRNA_variant", "exon_variant",
           "non_coding_exon_variant", "non_coding_transcript_exon_variant",
           "non_coding_transcript_variant", "nc_transcript_variant", "5_prime_UTR_variant",
           "5_prime_UTR_premature_start_codon_gain_variant", "3_prime_UTR_variant",
           "TF_binding_site_variant", "regulatory_region_variant", "regulatory_region",
           "intergenic_variant", "intergenic_region", "upstream_gene_variant",
           "downstream_gene_variant", "disruptive_inframe_deletion", "inframe_deletion", "inframe_insertion", "disruptive_inframe_insertion")

#Corresponding MAF Variant Classifications
maf.vc = c("Splice_Site", "Splice_Site", "Splice_Site", "Splice_Site",
           "Nonsense_Mutation", "Nonstop_Mutation", "Translation_Start_Site",
           "Translation_Start_Site", "Missense_Mutation", "Missense_Mutation",
           "Missense_Mutation", "Missense_Mutation", "Intron", "Intron",
           "Intron", "Intron", "Splice_Region", "Silent", "Silent", "Silent",
           "Silent", "RNA", "RNA", "RNA", "RNA", "RNA", "RNA", "5'UTR",
           "5'UTR", "3'UTR", "IGR", "IGR", "IGR", "IGR", "IGR", "5'Flank",
           "3'Flank", "In_Frame_Del", "In_Frame_Del", "In_Frame_Ins", "In_Frame_Ins")
names(maf.vc) = vep.vc

### Set Unified Colour Settings ###
var_types <- c("3'Flank","5'Flank","5'UTR","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Region","Splice_Site")
col.values <- c("#00BBDA","#E18A00","#BE9C00","#00BDC2","#24B700","#00BE70","#00C1AB","#F8766D","#8B93FF","#D575FE","#F962DD")
#col.values <- c("lightgrey","grey","darkgrey","red","blue","darkred","darkblue","green","darkviolet","darkorchid3","seagreen2","wheat1")
names(col.values) <- var_types

clear.classifications <- c("3'Flank","5'Flank","5'UTR","Frame Shift Deletion","Frame Shift Insertion","In Frame Insertion","Missense Mutation","Nonsense Mutation","Nonstop Mutation","Splice Region","Splice Site")
names(clear.classifications) <- var_types

simple_classifications <- c("Flank3","Flank5","UTR5","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Region","Splice_Site")
names(simple_classifications) <- var_types

##############################
### 1% SYNONYMOUS VARIANTS ###
##############################
output_prefix <- "LBC_ARCHER.1PCT_VAF.Mar21."

files <- list.files(path = paste(data_origin_dir, "/processed_synonymous_variants", sep = ""), pattern = ".tsv", full.names = T)
temp <- lapply(files, fread, sep="\t", header = T)
data <- rbindlist( temp )
data <- as.data.frame(data)
data$key <- paste(data$PreferredSymbol, data$base_substitution, sep=" ")
data$p_key <- paste(data$PreferredSymbol, data$protein_substitution, sep=" ")
data$event_key <- paste(data$PreferredSymbol, data$base_substitution, data$participant_id, sep="_")
colnames(data)[33] <- "AUAO" # Non unique column field (UAO) - taken from ARCHER sequencing information and probably not required, may remove in processing
dim(data)

data.1pct <- subset(data, AF >= 0.01)
dim(data.1pct)
data <- data[ data$event_key %in% data.1pct$event_key, ]
dim(data)

data <- subset(data, !(key %in% ERROR_CALLS_BASESUB$V1))
dim(data)
data <- subset(data, !(key %in% PROBABLE_GERMLINE$V1))
dim(data)
table(data$consequence)
data <- data[!(data$consequence %in% CONSEQUENCE_FILTER) , ]
dim(data)
table(data$consequence)
data$Variant_Classification = maf.vc[as.character(data$consequence)]
data$Variant_Classification <- ifelse((data$consequence == "frameshift_variant" & data$type == "ins"), "Frame_Shift_Ins", ifelse((data$consequence == "frameshift_variant" & data$type == "del"), "Frame_Shift_Del", data$Variant_Classification))
data$TYPE <- ifelse(data$type == "mnp", "INS", toupper(data$type))
dim(data)
write.table(data, paste(output_prefix, "synonymous.tsv", sep = ""), sep = "\t", row.names = F)



##################################
### 1% NON-SYNONYMOUS VARIANTS ###
##################################
output_prefix <- "LBC_ARCHER.1PCT_VAF.Mar21."

files <- list.files(path = paste(data_origin_dir, "/processed_variants", sep = ""), pattern = ".tsv", full.names = T)
temp <- lapply(files, fread, sep="\t", header = T)
data <- rbindlist( temp )
data <- as.data.frame(data)
data$key <- paste(data$PreferredSymbol, data$base_substitution, sep=" ")
data$p_key <- paste(data$PreferredSymbol, data$protein_substitution, sep=" ")
data$event_key <- paste(data$PreferredSymbol, data$base_substitution, data$participant_id, sep="_")
colnames(data)[33] <- "AUAO" # Non unique column field (UAO) - taken from ARCHER sequencing information and probably not required, may remove in processing
dim(data)

data.1pct <- subset(data, AF >= 0.01)
dim(data.1pct)
data <- data[ data$event_key %in% data.1pct$event_key, ]
dim(data)

data <- subset(data, !(key %in% ERROR_CALLS_BASESUB$V1))
dim(data)
data <- subset(data, !(key %in% PROBABLE_GERMLINE$V1))
dim(data)

#https://github.com/PoisonAlien/maftools/blob/master/R/icgc_to_maf.R
data$Variant_Classification = maf.vc[as.character(data$consequence)]
data$Variant_Classification <- ifelse((data$consequence == "frameshift_variant" & data$type == "ins"), "Frame_Shift_Ins", ifelse((data$consequence == "frameshift_variant" & data$type == "del"), "Frame_Shift_Del", data$Variant_Classification))
data$TYPE <- ifelse(data$type == "mnp", "INS", toupper(data$type))
#data$colour_code <- col.values[data$Variant_Classification]
dim(data)
write.table(data, paste(output_prefix, "non-synonymous.tsv", sep = ""), sep = "\t", row.names = F)



##############################
### 2% SYNONYMOUS VARIANTS ###
##############################
output_prefix <- "LBC_ARCHER.2PCT_VAF.Mar21."

files <- list.files(path = paste(data_origin_dir, "/processed_synonymous_variants", sep = ""), pattern = ".tsv", full.names = T)
temp <- lapply(files, fread, sep="\t", header = T)
data <- rbindlist( temp )
data <- as.data.frame(data)
data$key <- paste(data$PreferredSymbol, data$base_substitution, sep=" ")
data$p_key <- paste(data$PreferredSymbol, data$protein_substitution, sep=" ")
data$event_key <- paste(data$PreferredSymbol, data$base_substitution, data$participant_id, sep="_")
colnames(data)[33] <- "AUAO" # Non unique column field (UAO) - taken from ARCHER sequencing information and probably not required, may remove in processing
dim(data)

data.2pct <- subset(data, AF >= 0.02)
dim(data.2pct)
data <- data[ data$event_key %in% data.2pct$event_key, ]
dim(data)

data <- subset(data, !(key %in% ERROR_CALLS_BASESUB$V1))
dim(data)
data <- subset(data, !(key %in% PROBABLE_GERMLINE$V1))
dim(data)
table(data$consequence)
data <- data[!(data$consequence %in% CONSEQUENCE_FILTER) , ]
dim(data)
table(data$consequence)
data$Variant_Classification = maf.vc[as.character(data$consequence)]
data$Variant_Classification <- ifelse((data$consequence == "frameshift_variant" & data$type == "ins"), "Frame_Shift_Ins", ifelse((data$consequence == "frameshift_variant" & data$type == "del"), "Frame_Shift_Del", data$Variant_Classification))
data$TYPE <- ifelse(data$type == "mnp", "INS", toupper(data$type))
dim(data)
write.table(data, paste(output_prefix, "synonymous.tsv", sep = ""), sep = "\t", row.names = F)




##################################
### 2% NON-SYNONYMOUS VARIANTS ###
##################################
output_prefix <- "LBC_ARCHER.2PCT_VAF.Mar21."

files <- list.files(path = paste(data_origin_dir, "/processed_variants", sep = ""), pattern = ".tsv", full.names = T)
temp <- lapply(files, fread, sep="\t", header = T)
data <- rbindlist( temp )
data <- as.data.frame(data)
data$key <- paste(data$PreferredSymbol, data$base_substitution, sep=" ")
data$p_key <- paste(data$PreferredSymbol, data$protein_substitution, sep=" ")
data$event_key <- paste(data$PreferredSymbol, data$base_substitution, data$participant_id, sep="_")
colnames(data)[33] <- "AUAO" # Non unique column field (UAO) - taken from ARCHER sequencing information and probably not required, may remove in processing
dim(data)

data.2pct <- subset(data, AF >= 0.02)
dim(data.2pct)
data <- data[ data$event_key %in% data.2pct$event_key, ]
dim(data)

data <- subset(data, !(key %in% ERROR_CALLS_BASESUB$V1))
dim(data)
data <- subset(data, !(key %in% PROBABLE_GERMLINE$V1))
dim(data)

#https://github.com/PoisonAlien/maftools/blob/master/R/icgc_to_maf.R
data$Variant_Classification = maf.vc[as.character(data$consequence)]
data$Variant_Classification <- ifelse((data$consequence == "frameshift_variant" & data$type == "ins"), "Frame_Shift_Ins", ifelse((data$consequence == "frameshift_variant" & data$type == "del"), "Frame_Shift_Del", data$Variant_Classification))
data$TYPE <- ifelse(data$type == "mnp", "INS", toupper(data$type))
data$colour_code <- col.values[data$Variant_Classification]
data$clear_classification <- clear.classifications[data$Variant_Classification]
data$simple_classification <- simple_classifications[data$Variant_Classification]
dim(data)
table(data$colour_code)
write.table(data, paste(output_prefix, "non-synonymous.tsv", sep = ""), sep = "\t", row.names = F)


x <- read.table("LBC_ARCHER.2PCT_VAF.Mar21.non-synonymous.tsv", header = T, sep = "\t")
y <- read.table("LBC_ARCHER.1PCT_VAF.Mar21.non-synonymous.tsv", header = T, sep = "\t")

y.2 <- subset(y, AF >= 0.02)
length(unique(y.2$key))
length(unique(x$key))
setdiff(unique(y.2$key), unique(x$key))
setdiff(unique(x$key), unique(y.2$key))


### Create MAF Format / Plot Onco ###
require(maftools)
output_prefix <- "LBC_ARCHER.2PCT_VAF.Mar21."
data <- read.table(paste(output_prefix, "non-synonymous.tsv", sep = ""), sep = "\t", header = T)
data$gene_key <- paste(data$PreferredSymbol, data$base_substitution, data$participant_id, sep=" ")

colnames(data)[33] <- "AUAO" # Non unique column field (UAO) - taken from ARCHER sequencing information and probably not required, may remove in processing
dim(data)
data <- subset(data, wave == 1)
dim(data)

maf.data <- data.table::data.table(Hugo_Symbol = data$PreferredSymbol, Entrez_Gene_Id = NA, Center = NA, NCBI_Build = "GRCh38", Chromosome = data$chromosome,
                                   Start_Position = data$position, End_Position = (data$position + nchar(as.character(data$mutation))), Strand = '+',
                                   Variant_Classification = data$Variant_Classification, Variant_Type = data$TYPE, Reference_Allele = data$reference,
                                   Tumor_Seq_Allele1 = data$reference, Tumor_Seq_Allele2 = data$mutation, dbSNP_RS = NA, dbSNP_Val_Status  = NA,
                                   Tumor_Sample_Barcode = data$participant_id, Verification_Status = "NA", Mutation_Status = "Somatic",
                                   Sequence_Source = "ArcherDX/Illumina", Validation_Method = "NA", HGVSc = data$base_substitution, HGVSp = data$protein_substitution, 
                                   Consequence = data$consequence, Symbol = data$PreferredSymbol, AF = data$AF)
write.table(maf.data, paste(output_prefix, "non-synonymous.wave1.2PCT_VAF.maf", sep = ""), sep = "\t", row.names = F)

