###
#cd /mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq
#python ./code/parse_ARCHER_syn.py --manifest ARCHER_PanelSeq_Metadata_SEH_deduped.All.tsv --output_dir ./processed_data/Processed_1PCT_VAF_Jan2021/


### Set Key Parameters/ Read in Data / Filter ###
setwd("/mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq/processed_data/Processed_1PCT_VAF_Jan2021/")
output_prefix <- "LBC_ARCHER.1PCT_VAF.Jan2021."
data_origin_dir <- "/mnt/tchandra-lab/Neil/data/DNA-seq/LBC_ARCHER_Panel-seq/processed_data/Processed_1PCT_VAF_Jan2021"
waves <- c(1,2,3,4,5)

library(ggplot2)
library(data.table)
library(dplyr)
library(dplyr)
library(reshape)
library(pheatmap)
library(maftools)
library(stringr)
library(ggsci)


CONSEQUENCE_FILTER <- c("5_prime_UTR_variant", "coding_sequence_variant", "feature_elongation", "feature_truncation", "frameshift_variant", "incomplete_terminal_codon_variant", "inframe_deletion", "inframe_insertion", "missense_variant", "protein_altering_variant", "start_lost", "stop_gained", "stop_lost", "transcript_ablation", "transcript_amplification", "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant")
AO_FILTER <- 5
UAO_FILTER <- 3
GNOMAD_AF_FILTER <- 0.05
AF_FILTER <- 0.01
SYNONYMOUS_AF_FILTER <- 0.01
AF_GERMLINE_LIMIT <- 0.45
SEQ_BIAS_FILTER <- 'Yes'

ERROR_CALLS <- c("RUNX1 p.Glu422Pro","RUNX1 p.Glu422Ala","RUNX1 p.Glu422Gly","RUNX1 p.Glu422Val","RUNX1 p.Arg423Gly","RUNX1 p.Arg423Pro","RUNX1 p.Ser424Pro","RUNX1 p.Ser424Ala")
ERROR_CALLS_BASESUB <- c("TP53 c.74+38_74+39delinsGC", "TP53 c.280T>G", "PHF6 c.588+14del", "PHF6 c.588+14dup", "ASXL1 c.1934del", "ASXL1 c.1934dup", "ATRX c.5698-7del", "ATRX c.5698-7dup", "ATRX c.2658_2659del", "ATRX c.3786_3788del", "ATRX c.3619del", "BCOR c.3746-7T>G", "BCOR c.3746-8G>C", "BCOR c.3104A>C", "BCOR c.3052-4C>A", "BCOR c.3746-1G>T", "BCORL1 c.224A>C", "BCORL1 c.1454T>G", "BCORL1 c.1736A>C", "BCORL1 c.1445C>A", "BRAF c.1208del", "CBL c.1380_1382del", "CBL c.718G>T", "CBL c.591-2A>C", "CDKN2A c.52A>C", "CDKN2A c.221A>C", "CDKN2A c.442G>A", "CDKN2A c.150+13T>G", "CDKN2A c.229A>C", "CEBPA c.805A>G", "CEBPA c.292A>C", "CEBPA c.767T>G", "CSFR3 c.2383A>C", "CSFR3 c.2347A>C", "CUX1 c.1256-24491_1256-24489del", "CUX1 c.1255+6474T>G", "CUX1 c.1256-24646T>G", "CUX1 c.440-6C>A", "CUX1 c.1256-24315A>C", "CUX1 c.1256-24577A>G", "CUX1 c.1680+2T>G", "CUX1 c.1256-24643A>C", "CUX1 c.440-6C>T", "CUX1 c.1256-24319T>G", "CUX1 c.1256-24317T>G", "CUX1 c.1256-24580A>C", "CUX1 c.1452A>C", "CUX1 c.1255+38503A>C", "CUX1 c.1255+6485A>C", "CUX1 c.1255+6492A>C", "CUX1 c.1255+6455A>C", "CUX1 c.1255+6464A>C", "CUX1 c.1256-24644_1256-24640inv", "CUX1 c.1256-24638C>A", "CUX1 c.1255+6504T>G", "CUX1 c.1255+6461A>C", "CUX1 c.1255+6497T>G", "CUX1 c.1255+6009A>C", "CUX1 c.1255+6063A>C", "CUX1 c.1255+1055A>C", "DDX41 c.138+15T>G", "DDX41 c.1223T>A", "DDX41 c.138+17C>G", "DDX41 c.138+15T>C", "DDX41 c.1103A>C", "DDX41 c.1230+2T>C", "DDX41 c.1230+2_1230+3delinsCT", "DDX41 c.138+10A>C", "DDX41 c.138+19G>A", "ETV6 c.427_428del", "EZH2 c.118-5_118-4del", "EZH2 c.118-4dup", "EZH2 c.118-3A>T", "EZH2 c.118-6_118-4del", "EZH2 c.118-2A>T", "EZH2 c.1048A>C", "FLT3 c.1419-4del", "FLT3 c.1419-4dup", "FLT3 c.1419-5_1419-4del", "GATA2 c.341A>C", "GATA2 c.-6_-4del", "JAK2 c.3291+16del", "JAK2 c.3291+16dup", "JAK3 c.2483T>G", "JAK3 c.2438A>C", "KDM6A c.-264T>G", "KDM6A c.-183T>G", "KDM6A c.2703-4C>T", "KDM6A c.-218A>C", "KDM6A c.-36del", "KDM6A c.-262T>G", "KDM6A c.-194T>G", "KMT2A c.269C>T", "KMT2A c.173del", "KMT2A c.397T>G", "LUC7L2 c.757_758del", "LUC7L2 c.757_758dup", "LUC7L2 c.17A>C", "LUC7L2 c.-14G>C", "LUC7L2 c.19A>C", "LUC7L2 c.1A>C", "LUC7L2 c.28A>C", "MYC c.154_156del", "NF1 c.61-3C>A", "NF1 c.41T>G", "NF1 c.38T>G", "NOTCH1 c.7244_7246del", "NOTCH1 c.7502A>C", "NPM c.847-4C>T", "PHF6 c.588+14del", "PHF6 c.588+14dup", "RAD21 c.815-6T>A", "RAD21 c.1162-6T>A", "RAD21 c.1162-5A>T", "RAD21 c.815-2A>G", "RAD21 c.815-2A>T", "RAD21 c.815-5T>A", "RAD21 c.815-1G>T", "RAD21 c.815-7T>A", "RAD21 c.815-1_815delinsTG", "RAD21 c.815T>G", "RAD21 c.815-4_815-3delinsTA", "RAD21 c.815-3C>G", "RUNX1 c.614-7A>C", "RUNX1 c.614-3T>A", "SH2B3 c.1234A>G", "SH2B3 c.1235A>G", "SRSF2 c.362+12T>G", "SRSF2 c.559T>G", "SRSF2 c.610T>G", "SRSF2 c.301A>C", "SRSF2 c.339T>G", "STAG2 c.3278-4C>A", "USAF2 c.259_261del",
                         "ABL1 c.620A>C","ASXL1 c.1281dup","ASXL1 c.1854dup","ASXL1 c.1933_1934dup","ASXL1 c.4127dup","ASXL1 c.25_27del","ATRX c.1074dup","ATRX c.2518dup","ATRX c.4744dup","ATRX c.6111-3dup","BCOR c.1005dup","BCOR c.1199dup","BCOR c.1722dup","BCOR c.2514dup","BCOR c.266T>G","BCOR c.2811dup","BCOR c.3746-7T>C","BCOR c.4157A>C","BCOR c.4268_4273del","BCOR c.4760dup","BCOR c.4834dup","BCORL1 c.1103C>A","BCORL1 c.1309T>C","BCORL1 c.1403C>T","BCORL1 c.1442T>C","BCORL1 c.1444_1445insAGGTAGGA","BCORL1 c.1447_1452inv","BCORL1 c.1447_1454del","BCORL1 c.172dup","BCORL1 c.175A>G","BCORL1 c.1929del","BCORL1 c.1942A>C","BCORL1 c.1947dup","BCORL1 c.1954A>C","BCORL1 c.220dup","BCORL1 c.2222A>C","BCORL1 c.249dup","BCORL1 c.2518T>G","BCORL1 c.2911dup","BCORL1 c.3001dup","BCORL1 c.331_332delinsCG","BCORL1 c.3496dup","BCORL1 c.5042dup","BCORL1 c.617dup","BCORL1 c.644A>G","BCORL1 c.823dup","BRAF c.1208dup","CBLC c.1303C>T","CBLC c.1332dup","CDKN2A c.407dup","CDKN2A c.97dup","CEBPA c.14A>C","CEBPA c.347G>A","CEBPA c.457C>G","CEBPA c.53A>C","CEBPA c.561_562delinsCT","CEBPA c.566C>A","CEBPA c.568T>C","CEBPA c.584_589dup","CEBPA c.61A>C","CEBPA c.68dup","CEBPA c.71A>C","CEBPA c.722T>G","CEBPA c.79A>C","CEBPA c.82A>C","CEBPA c.97T>G","CSF3R c.2076C>T","CSF3R c.2346dup","CSF3R c.2383A>C","CSF3R c.2461dup","CUX1 c.1255+1097dup","CUX1 c.1255+38593A>G","CUX1 c.1255+5923dup","CUX1 c.1255+6033A>C","CUX1 c.1255+6440A>C","CUX1 c.1255+6449A>C","CUX1 c.1255+6465A>C","CUX1 c.1256-24473A>C","CUX1 c.1256-24482dup","CUX1 c.1256-24491A>C","CUX1 c.1256-24505A>G","CUX1 c.1256-24527A>C","CUX1 c.1256-24529C>G","CUX1 c.1468A>C","CUX1 c.1967+12dup","DDX41 c.1394dup","DNMT3A c.1238dup","DNMT3A c.176dup","DNMT3A c.327dup","ETV6 c.613dup","EZH2 c.1184dup","EZH2 c.2187dup","FBXW7 c.1761dup","GATA2 c.302dup","GATA2 c.568dup","GATA2 c.716A>C","GATA2 c.742A>C","GATA2 c.818dup","IDH1 c.13dup","IDH1 c.96del","IDH2 c.435dup","IDH2 c.679-6dup","IKZF1 c.316A>C","IKZF1 c.56dup","JAK3 c.1915-5dup","JAK3 c.1978dup","JAK3 c.251dup","KDM6A c.-21T>C","KDM6A c.-36dup","KDM6A c.2080dup","KDM6A c.3712G>C","KDM6A c.3714_3715inv","KDM6A c.3717G>C","KMT2A c.11271dup","KMT2A c.1142dup","KMT2A c.1147dup","KMT2A c.1660dup","KMT2A c.173dup","KMT2A c.338A>C","KMT2A c.395T>G","KMT2A c.406A>C","KMT2A c.432+9dup","KMT2A c.766dup","KMT2A c.7887dup","KMT2A c.8445dup","LUC7L2 c.11A>C","LUC7L2 c.32T>G","LUC7L2 c.61+6T>G","MYC c.173dup","MYD88 c.683+22T>G","NF1 c.1527+14dup","NF1 c.2033dup","NF1 c.3197+9dup","NF1 c.3590C>T","NF1 c.4076dup","NF1 c.4110+3865dup","NF1 c.60+4A>C","NF1 c.7638dup","NF1 c.833del","NOTCH1 c.6392dup","NOTCH1 c.7465A>C","NOTCH1 c.7497C>G","NOTCH1 c.7499A>C","NOTCH1 c.7508A>C","NOTCH1 c.7520A>C","NOTCH1 c.7531A>C","NOTCH1 c.7541_7542del","PPM1D c.1636dup","PPM1D c.1589del","RAD21 c.1640dup","SH2B3 c.1038dup","SH2B3 c.1566dup","SH2B3 c.314A>C","SMC1A c.109+581dup","SMC1A c.1771dup","SMC1A c.2878dup","SMC1A c.3069dup","SMC1A c.3453dup","SMC3 c.2567dup","SRSF2 c.20dup","SRSF2 c.287dup","SRSF2 c.547T>G","SRSF2 c.630dup","STAG2 c.1400dup","TET2 c.1842dup","TET2 c.1985dup","TP53 c.-64T>G","TP53 c.1083dup","TP53 c.216dup","TP53 c.404G>T","TP53 c.548C>T","TP53 c.554G>A","U2AF1 c.-19C>T","U2AF1 c.1043dup","U2AF1 c.134A>C","U2AF1 c.231-2A>C","U2AF1 c.290dup","U2AF1 c.487-2A>C","U2AF1 c.531dup","U2AF1 c.742+11dup","XPO1 c.1664dup","ZRSR2 c.1072G>C","ZRSR2 c.1076A>G","ZRSR2 c.1084G>C","ZRSR2 c.1338_1343dup",
                         "BCORL1 c.644A>C","DDX41 c.1760dup","DDX41 c.701T>C","CUX1 c.1255+1097del","CUX1 c.1255+6070dup","DNMT3A c.2479-2A>T","DNMT3A c.2500dup","DNMT3A c.869_870del","DNMT3A c.899_900del","FLT3 c.2317dup","GATA2 c.599dup","KDM6A c.-183T>C","NOTCH1 c.7489A>C","NOTCH1 c.7517A>G","TP53 c.-72T>G","TP53 c.-78T>G","TP53 c.74+56dup","U2AF2 c.134A>C","U2AF2 c.231-2A>C","U2AF2 c.290dup","U2AF2 c.487-2A>C","U2AF2 c.531dup","U2AF2 c.742+11dup","BCORL1 c.4853+8T>G","BCORL1 c.872dup","DNMT3A c.1015-3dup","LUC7L2 c.-12T>C","LUC7L2 c.-6G>C","NOTCH1 c.7495A>C","SH2B3 c.1236+2_1236+4delinsGGG","TET2 c.267dup","U2AF c.-19C>T")

PROBABLE_GERMLINE <- c("ABL1 c.1546C>T","ABL1 c.1594A>T","ABL1 c.740A>G","ANKRD26 c.-140C>G","ASXL1 c.3306G>T","ASXL1 c.3306G>T","ATRX c.5787-8_5787-5del","ATRX c.1041T>G","ATRX c.3505A>C","CALR c.1226T>G","CBLC c.1364G>C","CCND2 c.785G>A","CDKN2A c.442G>A","CEBPA c.1044C>G","CSFR3 c.2422G>A","CSFR3 c.2041-75G>A","CUX1 c.1636A>C","DCK c.364C>T","ETV6 c.602T>C","FLT3 c.970G>A","FLT3 c.1249A>C","GATA2 c.481C>G","GATA2 c.30G>T","HRAS c.-10C>T","IDH2 c.94T>G","IZKF1 c.336C>G","JAK2 c.1776+5G>A","KDM6A c.1527+3G>A","KDM6A c.3098C>A","KDM6A c.-194T>G","KIT c.67+4G>A","KMT2A c.3634+4G>A","KMT2A c.5219C>T","KMT2A c.9659A>T","MPL c.1666G>T","MYC c.77A>G","BCORL1 c.3158A>G","KDM6A c.-167_-165dup","PTEN c.1001A>G","TET2 c.100C>T","TET2 c.1064G>A","TP53 c.108G>A","CBL c.1227+4C>T","CUX1 c.1256-24491_1256-24489dup","RUNX1 c.1354G>A","U2AF2 c.1411C>A","DDX41 c.1549+7G>T","MPL c.1565+5C>T","CUX1 c.1573C>G","IKZF1 c.161-8350A>C","RUNX1 c.167T>C","SH2B3 c.1697G>A","STAG2 c.1705G>A","CSF3R c.2041-35T>C","CSF3R c.2041-75G>A","SF3B1 c.2078-8T>A","ASXL1 c.2110G>A","CBL c.227C>T","CSF3R c.2422G>A","ATRX c.2875G>T","DDX41 c.291G>A","CUX1 c.295G>A","BCORL1 c.3158A>G","TET2 c.3251A>C","IKZF1 c.336C>G","JAK2 c.3371A>C","KMT2A c.3796C>T","ASXL1 c.4493C>T","ATRX c.5579A>G","SH2B3 c.557G>T","WT1 c.55C>T","NF1 c.6084+8C>G","NOTCH1 c.6853G>A","NOTCH1 c.7400C>T","SRSF2 c.75C>T","GATA2 c.832T>A","GATA1 c.85G>C","KMT2A c.89C>G","IDH1 c.94T>G","BCORL1 c.3765_3767del","NRAS c.553C>T","TET2 c.2599T>C","BCORL1 c.1606A>G","ATRX c.2595C>G","ATRX c.2595C>G","CUX1 c.1633A>G","JAK2 c.3323A>G","KMT2A c.7174T>C","NF1 c.2032C>G","NOTCH1 c.6745G>A","SH2B3 c.1523G>A","CSF3R c.2041-30C>T")

SYN_ERROR_CALLS <- c("RUNX1 p.Gly420=", "RUNX1 p.Gly421=", "RUNX1 p.Ser424=", "RUNX1 p.Pro425=")
SYN_ERROR_CALLS_BASESUB <- c("ABL1 c.1497A>G","ASXL1 c.12A>G","ASXL1 c.3759T>C","BCOR c.1260T>C","BCOR c.1692A>G","BROR c.2787C>T","BCORL1 c.1080A>C","BCORL1 c.1731A>C","BCROL1 c.1923A>C","CBLB c.1341A>C","CEBPA c.561G>C","CEBPA c.690G>T","CSF3R c.1254T>C","CSF3R c.1260T>C","CUX1 c.1137G>A","CUX1 c.1255+6283T>G","CUX1 c.1255+6430A>C","CUX1 c.1256-24516G>A","CUX1 c.1256-24516G>C","CUX1 c.1266A>G","CUX1 c.1458A>C","CUX1 c.1461A>C","CUX1 c.1464A>C","CUX1 c.1827T>G","CUX1 c.6G>A","CXCR4 c.783C>T","DCK c.300C>T","DDX41 c.1107A>C","DDX41 c.1200C>T","DDX41 c.897G>A","DNMT3A c.1266G>A","DNMT3A c.27C>T","DNMT3A c.318T>G","ETV6 c.258G>A","EZH2 c.1731G>A","EZH2 c.234C>T","FBXW7 c.1368G>A","FLT3 c.1683A>G","GATA2 c.42G>T","GATA2 c.768T>G","GATA2 c.819A>G","IDH1 c.315C>T","JAK2 c.2490G>A","KIT c.1638A>G","KIT c.2586G>C","KMT2A c.1710A>C","KMT2A c.4284A>C","KMT2A c.5670A>G","KMT2A c.7245C>T","NF1 c.5694G>A","NF1 c.702G>A","NOTCH1 c.4965T>G","NOTCH1 c.6555C>T","PDGFRA c.2472C>T","PTEN c.900C>T","RAD21 c.1440T>C","RBBP6 c.5019A>G","SF3B1 c.2631T>C","SH2B3 c.1257A>C","SH2B3 c.243C>T","SH2B3 c.318A>C","WT1 c.261T>G","WT1 c.330C>T","ZRSR2 c.864C>T")

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

col.vars <- c("3'Flank","5'Flank","5'UTR","Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Intron","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Region","Splice_Site")
col.values <- c("#00BBDA","#E18A00","#BE9C00","#8CAB00","#24B700","#00BE70","#00C1AB","#F8766D","#00ACFC","#8B93FF","#D575FE","#F962DD","#FF65AC")
names(col.values) <- col.vars

###########################
### SYNONYMOUS VARIANTS ###
###########################
files <- list.files(path = paste(data_origin_dir, "/processed_synonymous_variants", sep = ""), pattern = ".tsv", full.names = T)
temp <- lapply(files, fread, sep="\t", header = T)
data <- rbindlist( temp )
data <- as.data.frame(data)
data$key <- paste(data$PreferredSymbol, data$base_substitution, sep=" ")
data$p_key <- paste(data$PreferredSymbol, data$protein_substitution, sep=" ")
colnames(data)[33] <- "AUAO" # Non unique column field (UAO) - taken from ARCHER sequencing information and probably not required, may remove in processing

dim(data)
#data <- subset(data, AO >= AO_FILTER)
#data <- subset(data, UAO >= UAO_FILTER)
#data <- subset(data, AF >= AF_FILTER)
#data <- subset(data, AF <= AF_GERMLINE_LIMIT)
#data <- subset(data, HasSeqDirBias != SEQ_BIAS_FILTER)
data <- subset(data, !(p_key %in% SYN_ERROR_CALLS))
data <- subset(data, !(key %in% SYN_ERROR_CALLS_BASESUB))
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



###############################
### NON-SYNONYMOUS VARIANTS ###
###############################
files <- list.files(path = paste(data_origin_dir, "/processed_variants", sep = ""), pattern = ".tsv", full.names = T)
temp <- lapply(files, fread, sep="\t", header = T)
data <- rbindlist( temp )
data <- as.data.frame(data)
data$key <- paste(data$PreferredSymbol, data$base_substitution, sep=" ")
data$p_key <- paste(data$PreferredSymbol, data$protein_substitution, sep=" ")
colnames(data)[33] <- "AUAO" # Non unique column field (UAO) - taken from ARCHER sequencing information and probably not required, may remove in processing

dim(data)
#data <- subset(data, AO >= AO_FILTER)
#data <- subset(data, UAO >= UAO_FILTER)
#data <- subset(data, AF >= SYNONYMOUS_AF_FILTER)
#data <- subset(data, HasSeqDirBias != SEQ_BIAS_FILTER)
data <- subset(data, !(p_key %in% ERROR_CALLS))
data <- subset(data, !(key %in% ERROR_CALLS_BASESUB))
data <- subset(data, !(key %in% PROBABLE_GERMLINE))
dim(data)
#data.wave1 <- subset(data, wave == 1)
#variant_counts <- as.data.frame(table(data.wave1$key))
#variant_counts <- variant_counts[order(-variant_counts$Freq), ]
#head(variant_counts, n = 35)
#variant_counts <- subset(variant_counts, Freq > 9)
#dim(variant_counts)
#data <- data[ !(data$key %in% variant_counts$Var1),]
#https://github.com/PoisonAlien/maftools/blob/master/R/icgc_to_maf.R
data$Variant_Classification = maf.vc[as.character(data$consequence)]
data$Variant_Classification <- ifelse((data$consequence == "frameshift_variant" & data$type == "ins"), "Frame_Shift_Ins", ifelse((data$consequence == "frameshift_variant" & data$type == "del"), "Frame_Shift_Del", data$Variant_Classification))
data$TYPE <- ifelse(data$type == "mnp", "INS", toupper(data$type))
data$colour_code <- col.values[data$Variant_Classification]
dim(data)
write.table(data, paste(output_prefix, "non-synonymous.tsv", sep = ""), sep = "\t", row.names = F)



### Plot VAF Histograms ###

data.syn <- read.table(paste(output_prefix, "synonymous.tsv", sep = ""), sep = "\t", header = T)
data.non_syn <- read.table(paste(output_prefix, "non-synonymous.tsv", sep = ""), sep = "\t", header = T)
data.non_syn$gene_key <- paste(data$PreferredSymbol, data$base_substitution, data$participant_id, sep=" ")

colnames(data.non_syn)[33] <- "AUAO" # Non unique column field (UAO) - taken from ARCHER sequencing information and probably not required, may remove in processing
dim(data.non_syn)
data <- data.non_syn

dir.create("./plotting/")

output_filename <- paste("./plotting/", output_prefix, "VAF_Histogram.png", sep = "")
png(output_filename, width = 1000, height = 1100, res = 200, pointsize =5)
p <- ggplot(data, aes(x=AF, colour = as.factor(wave))) + geom_histogram(bins = 100) + xlim(0, .6)
p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + ylab("Cumulative Count")
p <- p + xlab("Variant Allele Fraction")
p <- p + theme(axis.text.x = element_text(color="black", size=12),
               axis.text.y = element_text(color="black", size=12))
p <- p + theme(legend.title=element_text(), legend.position="bottom") + scale_colour_manual("Wave",values=c("red","orange", "green", "blue", "pink"))
p <- p + theme(axis.text=element_text(size=14),
               axis.title=element_text(size=14))
print(p)
dev.off()


output_filename <- paste("./plotting/", output_prefix, "VAF_Density.png", sep = "")
png(output_filename, width = 1000, height = 1100, res = 200, pointsize =5)
p <- ggplot(data, aes(x=AF, colour = as.factor(wave), group = as.factor(wave))) + geom_density() + xlim(0, .45)
p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + ylab("Density")
p <- p + xlab("Variant Allele Fraction")
p <- p + theme(axis.text.x = element_text(color="black", size=12),
               axis.text.y = element_text(color="black", size=12))
p <- p + theme(legend.title=element_text(), legend.position="bottom") + scale_colour_manual("Wave",values=c("red","orange", "green", "blue", "pink"))
p <- p + theme(axis.text=element_text(size=14),
               axis.title=element_text(size=14))
print(p)
dev.off()

data <- subset(data, AF <= .45)
data.wave1 <- subset(data, wave == 1)
data.wave2 <- subset(data, wave == 2)
data.wave3 <- subset(data, wave == 3)
data.wave4 <- subset(data, wave == 4)
data.wave5 <- subset(data, wave == 5)

wave1_obs <- length(unique(data.wave1$participant_id))
wave2_obs <- length(unique(data.wave2$participant_id))
wave3_obs <- length(unique(data.wave3$participant_id))
wave4_obs <- length(unique(data.wave4$participant_id))
wave5_obs <- length(unique(data.wave5$participant_id))

minimum_obs <- min(c(wave1_obs,wave2_obs,wave3_obs,wave4_obs,wave5_obs))
wave1_norm <- minimum_obs/wave1_obs
wave2_norm <- minimum_obs/wave2_obs
wave3_norm <- minimum_obs/wave3_obs
wave4_norm <- minimum_obs/wave4_obs
wave5_norm <- minimum_obs/wave5_obs

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

group_tags5 <- cut(data.wave5$AF,
                   breaks=breaks, 
                   include.lowest=TRUE, 
                   right=FALSE, 
                   labels=tags)
#summary(group_tags)
#as_tibble(group_tags)
x <- rbind(table(group_tags1), table(group_tags2), table(group_tags3), table(group_tags4), table(group_tags5))
x <- as.data.frame(t(x))
colnames(x) <- c("wave1", "wave2", "wave3", "wave4", "wave5")

x$wave1 <- x$wave1*wave1_norm
x$wave2 <- x$wave2*wave2_norm
x$wave3 <- x$wave3*wave3_norm
x$wave4 <- x$wave4*wave4_norm
x$wave5 <- x$wave5*wave5_norm

x_labels <- c("","","","","","","","","","0.1","","","","","","","","","","0.2","","","","","","","","","","0.3","","","","","","","","","","0.4","","","","","")
y <- melt(setDT(x, keep.rownames = TRUE), "rn")

output_filename <- paste("./plotting/", output_prefix, "VAF_NormalisedHistogram.png", sep = "")
png(output_filename, width = 1000, height = 1100, res = 200, pointsize =5)
p <- ggplot(y, aes(fill=variable, y=value, x=rn)) + 
  geom_bar(position="stack", stat="identity")
p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + ylab("Cumulative Count Normalised by Observations")
p <- p + xlab("Variant Allele Fraction")
p <- p + theme(axis.text.x = element_text(color="black", size=12),
               axis.text.y = element_text(color="black", size=12))
p <- p + theme(legend.title=element_text(), legend.position="bottom") + scale_fill_manual("",values=c("red","orange", "green", "blue", "pink"))
p <- p + theme(axis.text=element_text(size=14),
               axis.title=element_text(size=14))
p <- p + scale_x_discrete(labels = x_labels)
print(p)
dev.off()




### Plot Gene/VAF Distribution ###
data.non_syn <- read.table(paste(output_prefix, "non-synonymous.tsv", sep = ""), sep = "\t", header = T)
data.non_syn$gene_key <- paste(data.non_syn$PreferredSymbol, data.non_syn$base_substitution, data.non_syn$participant_id, sep=" ")

data.non_syn <- subset(data.non_syn, wave == 1)
newcsv <- data.non_syn %>%
  group_by(PreferredSymbol) %>%
  summarise(
    Total_AF = sum(AF)
  )
dim(data.non_syn)
newcsv <- as.data.frame(newcsv)
order <- newcsv[ order(-newcsv$Total_AF) , ]

output_filename <- paste("./plotting/", output_prefix, "VAF_ByGene_BarChart.png", sep = "")
png(output_filename, width = 2000, height = 1000, res = 150, pointsize =5)

p <- ggplot(data.non_syn ) + 
  geom_boxplot( aes(x = PreferredSymbol, y = AF), outlier.shape = NA) + 
  geom_jitter(aes(x = PreferredSymbol, y = AF),  position=position_jitter(0.2), fill="grey30", alpha=.4, shape = 19, size = 0.7)
p <- p + xlab("")
p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                 legend.position=c(.84,.75), legend.text=element_text(size=11))
p <- p + theme(axis.text.x = element_text(color="black", size=8, angle = 90, hjust = 1, face = "bold"),
               axis.text.y = element_text(color="black", size=11))
p <- p + theme(axis.text=element_text(size=16),
               axis.title=element_text(size=16, face = "bold"))
p <- p + scale_x_discrete(limits=order$PreferredSymbol)
p <- p + scale_y_continuous(name="VAF (Wave 1)",limits=c(0,.6),breaks=c(0,.1,.2,.3,.4,.5,.6))
p
dev.off()

output_filename <- paste("./plotting/", output_prefix, "VAF_ByGene_BarChart.Coloured_Variant_Classification.png", sep = "")
png(output_filename, width = 2000, height = 1000, res = 150, pointsize =5)

p <- ggplot(data.non_syn ) + 
  geom_boxplot( aes(x = PreferredSymbol, y = AF), outlier.shape = NA) + 
  geom_jitter(aes(x = PreferredSymbol, y = AF, colour = Variant_Classification),  position=position_jitter(0.2), fill="grey30", alpha=.7, shape = 19, size = 0.7)
p <- p + xlab("")
p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                 legend.position=c(.84,.75), legend.text=element_text(size=11))
p <- p + theme(axis.text.x = element_text(color="black", size=8, angle = 90, hjust = 1, face = "bold"),
               axis.text.y = element_text(color="black", size=11))
p <- p + theme(axis.text=element_text(size=16),
               axis.title=element_text(size=16, face = "bold"))
p <- p + scale_x_discrete(limits=order$PreferredSymbol)
p <- p + scale_y_continuous(name="VAF (Wave 1)",limits=c(0,.6),breaks=c(0,.1,.2,.3,.4,.5,.6))
p <- p + guides(colour = guide_legend(override.aes = list(size=3), title = "Variant Classification"))
p
dev.off()






### Plot Cumulative VAF Heatmaps ###
data.non_syn <- read.table(paste(output_prefix, "non-synonymous.tsv", sep = ""), sep = "\t", header = T)
data.non_syn$gene_key <- paste(data.non_syn$PreferredSymbol, data.non_syn$base_substitution, data.non_syn$participant_id, sep=" ")

colnames(data.non_syn)[33] <- "AUAO" # Non unique column field (UAO) - taken from ARCHER sequencing information and probably not required, may remove in processing
dim(data.non_syn)
data <- data.non_syn
#data <- subset(data, AF >= .2)
for (cur_wave in waves) {
  # https://stackoverflow.com/questions/9617348/reshape-three-column-data-frame-to-matrix-long-to-wide-format
  data.wave1 <- subset(data, wave == cur_wave)
  newcsv <- data.wave1 %>%
    group_by(PreferredSymbol, participant_id) %>%
    summarise(
      Total_AF = sum(AF)
    )
  newcsv <- as.data.frame(newcsv)
  colnames(newcsv) <- c("PreferredSymbol", "participant_id", "Total_AF")
  newcsv <- reshape2::dcast(as.data.frame(newcsv), PreferredSymbol~participant_id, value.var="Total_AF")
  newcsv[is.na(newcsv)] = 0
  rownames(newcsv) <- newcsv$PreferredSymbol
  newcsv <- newcsv[ , !(colnames(newcsv) %in% c("PreferredSymbol")) ]
  p <- pheatmap(newcsv)
  output_filename <- paste("./plotting/", output_prefix, "Wave", cur_wave ,"_Cumulative_VAF_Heatmap.png", sep = "")
  png(output_filename, width = 2500, height = 2000, res = 200, pointsize =5)
  print(p)
  dev.off()
}



### Plot Largest Clone Variant Type Heatmaps ###
# https://stackoverflow.com/questions/9617348/reshape-three-column-data-frame-to-matrix-long-to-wide-format
data.non_syn <- read.table(paste(output_prefix, "non-synonymous.tsv", sep = ""), sep = "\t", header = T)
data.non_syn$gene_key <- paste(data$PreferredSymbol, data$base_substitution, data$participant_id, sep=" ")

colnames(data.non_syn)[33] <- "AUAO" # Non unique column field (UAO) - taken from ARCHER sequencing information and probably not required, may remove in processing
dim(data.non_syn)
data <- data.non_syn

for (cur_wave in waves) {
  data.wave1 <- subset(data, wave == cur_wave)
  newcsv <- data.wave1 %>%
    group_by(PreferredSymbol, participant_id) %>%
    summarise(
      largest_clone_AF = max(AF)
    )
  newcsv <- as.data.frame(newcsv)
  colnames(newcsv) <- c("PreferredSymbol", "participant_id", "largest_clone_AF")
  AF_map <- newcsv
  newcsv <- newcsv[ , colnames(newcsv) %in% c("PreferredSymbol", "participant_id", "largest_clone_AF")]
  newcsv <- reshape2::dcast(as.data.frame(newcsv), PreferredSymbol~participant_id, value.var="largest_clone_AF")
  newcsv[is.na(newcsv)] = 0
  rownames(newcsv) <- newcsv$PreferredSymbol
  newcsv <- newcsv[ , !(colnames(newcsv) %in% c("PreferredSymbol")) ]
  newcsv <- newcsv[apply(newcsv, 1, function(x) any(x != 0 | is.na(x))), ]
  p <- pheatmap(newcsv)
  
  output_filename <- paste("./plotting/", output_prefix, "Wave", cur_wave, "_Largest_Clone_VAF.png", sep = "")
  png(output_filename, width = 2500, height = 2000, res = 200, pointsize =5)
  print(p)
  dev.off()
}

### Plot Largest Clone Bubble Heatmaps ###
# https://stackoverflow.com/questions/9617348/reshape-three-column-data-frame-to-matrix-long-to-wide-format
data.non_syn <- read.table(paste(output_prefix, "non-synonymous.tsv", sep = ""), sep = "\t", header = T)
data.non_syn$gene_key <- paste(data.non_syn$PreferredSymbol, data.non_syn$base_substitution, data.non_syn$participant_id, sep=" ")
colnames(data.non_syn)[33] <- "AUAO" # Non unique column field (UAO) - taken from ARCHER sequencing information and probably not required, may remove in processing
dim(data.non_syn)
data <- data.non_syn

for (cur_wave in waves) {
  # https://stackoverflow.com/questions/9617348/reshape-three-column-data-frame-to-matrix-long-to-wide-format
  data.wave1 <- subset(data, wave == cur_wave)
  
  newcsv <- data.wave1 %>%
    group_by(PreferredSymbol, participant_id) %>% slice(which.max(AF))
  
  newcsv <- as.data.frame(newcsv)
  for_order <- reshape2::dcast(as.data.frame(newcsv), PreferredSymbol~participant_id, value.var="AF")
  for_order[is.na(for_order)] = 0
  rownames(for_order) <- for_order$PreferredSymbol
  for_order <- for_order[ , !(colnames(for_order) %in% c("PreferredSymbol")) ]
  for_order <- for_order[apply(for_order, 1, function(x) any(x != 0 | is.na(x))), ]
  cluster_order <- pheatmap(for_order)
  
  
  p <- ggplot(newcsv, aes(participant_id, PreferredSymbol)) +
    geom_point(aes(size = AF, colour = Variant_Classification), alpha=0.8) 
  p <- p + scale_size(range = c(0,15), limits = c(0,1)) 
  p <- p + xlab("") + ylab("")
  p <- p + theme_minimal() + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                   legend.text=element_text(size=11))
  
  p <- p + theme(axis.text.x = element_text(color="black", size=8, angle = 90, hjust = 1, face = "bold"),
                 axis.text.y = element_text(color="black", size=9))
  p <- p + theme(axis.text=element_text(size=16),
                 axis.title=element_text(size=16, face = "bold"))
  p <- p + scale_y_discrete(limits=rev(cluster_order$tree_row$labels[cluster_order$tree_row$order]))
  p <- p + scale_x_discrete(limits=(cluster_order$tree_col$labels[cluster_order$tree_col$order]))
  p <- p + guides(colour = guide_legend(override.aes = list(size=3), title = c("Variant Classification", "VAF")), labels = data$Variant_Classification)
  p
  
  output_filename <- paste("./plotting/", output_prefix, "Wave", cur_wave, "_Largest_Clone_VAF.Bubble.png", sep = "")
  png(output_filename, width = 2500, height = 2000, res = 200, pointsize =5)
  par(mfrow = c(1, 1), mar = c(4, 5, 6, 1))
  print(p)
  dev.off()
}





### Interesting Gene Trajectories ###
data.non_syn <- read.table(paste(output_prefix, "non-synonymous.tsv", sep = ""), sep = "\t", header = T)
data.non_syn$gene_key <- paste(data.non_syn$PreferredSymbol, data.non_syn$base_substitution, data.non_syn$participant_id, sep=" ")
genes <- unique(data.non_syn$PreferredSymbol)

genes_of_interest <- c("DNMT3A", "TET2", "JAK2", "ASXL1")
id <- "DNMT3A"
for (id in genes_of_interest) {
  
  print (id)
  data.non_syn.id <- subset(data.non_syn, PreferredSymbol == id)
  
  print(dim(data.non_syn.id))
  
  p <- ggplot() + geom_line(data = data.non_syn.id, aes(x=wave, y=AF, col = Variant_Classification, group = gene_key), size = 1.2) + 
    geom_point(data = data.non_syn.id, aes(x=wave, y=AF, col = Variant_Classification, group = gene_key))
  p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                   legend.text = element_text(size = 8),
                                   axis.text=element_text(size=16),
                                   axis.title=element_text(size=16))
  p <- p + ylab("Variant Allele Fraction")
  p <- p + scale_x_continuous(name="Wave",limits=c(1,5),breaks=c(1,2,3,4,5))
  p <- p + ggtitle(paste(id, "Non-Synonymous Variants", sep = " "))
  
  output_filename <- paste("./plotting/", output_prefix, id, ".GeneTimeCourse.OfInterest.png", sep = "")
  png(output_filename, width = 2400, height = 1000, res = 200, pointsize =5)
  print(p)
  dev.off()
}










### Plot Participant TimeCourses ###
data.non_syn <- read.table(paste(output_prefix, "non-synonymous.tsv", sep = ""), sep = "\t", header = T)
data.non_syn$gene_key <- paste(data.non_syn$PreferredSymbol, data.non_syn$base_substitution, data.non_syn$participant_id, sep=" ")

colnames(data.non_syn)[33] <- "AUAO" # Non unique column field (UAO) - taken from ARCHER sequencing information and probably not required, may remove in processing
dim(data.non_syn)

dir.create("./plotting/participant_timecourses/")
sample_ids <- unique(data.non_syn$participant_id)
#id <- "LBC360021"
for (id in sample_ids) {
  print (id)
  data.non_syn.id <- subset(data.non_syn, participant_id == id)
  print(dim(data.non_syn.id))
  
  p <- ggplot() + geom_line(data = data.non_syn.id, aes(x=wave, y=AF, col = PreferredSymbol, group = key), size = 1.2) + 
    geom_point(data = data.non_syn.id, aes(x=wave, y=AF, col = PreferredSymbol, group = key))
  p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                   legend.text = element_text(size = 8),
                                   axis.text=element_text(size=16),
                                   axis.title=element_text(size=16))
  p <- p + ylab("Variant Allele Fraction")
  p <- p + xlab("Wave")
  p <- p + ggtitle(id)
  
  output_filename <- paste("./plotting/participant_timecourses/", output_prefix, id, ".CrudeChronology.png", sep = "")
  png(output_filename, width = 1400, height = 1000, res = 200, pointsize =5)
  print(p)
  dev.off()
}


### Plot Gene TimeCourses ###
dir.create("./plotting/gene_timecourses/")
data.non_syn <- read.table(paste(output_prefix, "non-synonymous.tsv", sep = ""), sep = "\t", header = T)
data.non_syn$gene_key <- paste(data.non_syn$PreferredSymbol, data.non_syn$base_substitution, data.non_syn$participant_id, sep=" ")
genes <- unique(data.non_syn$PreferredSymbol)

for (id in genes) {
  
  print (id)
  data.non_syn.id <- subset(data.non_syn, PreferredSymbol == id)
  
  print(dim(data.non_syn.id))
  
  p <- ggplot() + geom_line(data = data.non_syn.id, aes(x=wave, y=AF, col = key, group = gene_key), size = 1.2) + 
    geom_point(data = data.non_syn.id, aes(x=wave, y=AF, col = key, group = gene_key))
  p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                   legend.text = element_text(size = 8),
                                   axis.text=element_text(size=16),
                                   axis.title=element_text(size=16))
  p <- p + ylab("Variant Allele Fraction")
  p <- p + xlab("Wave")
  p <- p + ggtitle(paste(id, "Non-Synonymous Variants", sep = " "))
  
  output_filename <- paste("./plotting/gene_timecourses/", output_prefix, id, ".GeneTimeCourse.png", sep = "")
  png(output_filename, width = 2400, height = 1000, res = 200, pointsize =5)
  print(p)
  dev.off()
}


### Plot Synonymous TimeCourses ###
dir.create("./plotting/gene_timecourses_synonymous/")
data.syn <- read.table(paste(output_prefix, "synonymous.tsv", sep = ""), sep = "\t", header = T)
data.syn$gene_key <- paste(data.syn$PreferredSymbol, data.syn$base_substitution, data.syn$participant_id, sep=" ")
genes <- unique(data.syn$PreferredSymbol)

for (id in genes) {
  print (id)
  data.syn.id <- subset(data.syn, PreferredSymbol == id)
  print(dim(data.syn.id))
  
  p <- ggplot() + geom_line(data = data.syn.id, aes(x=wave, y=AF, col = p_key, group = gene_key), size = 1.2) + 
    geom_point(data = data.syn.id, aes(x=wave, y=AF, col = p_key, group = gene_key))
  p <- p + theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                   legend.text = element_text(size = 8),
                                   axis.text=element_text(size=16),
                                   axis.title=element_text(size=16))
  p <- p + ylab("Variant Allele Fraction")
  p <- p + xlab("Wave")
  p <- p + ggtitle(paste(id, "Synonymous Variants", sep = " "))
  
  output_filename <- paste("./plotting/gene_timecourses_synonymous/", output_prefix, id, ".GeneTimeCourse.Synonymous.png", sep = "")
  png(output_filename, width = 1600, height = 1000, res = 200, pointsize =5)
  print(p)
  dev.off()
}


### Create MAF Format / Plot Onco ###
require(maftools)
data <- read.table(paste(output_prefix, "non-synonymous.tsv", sep = ""), sep = "\t", header = T)
data$gene_key <- paste(data$PreferredSymbol, data$base_substitution, data$participant_id, sep=" ")

colnames(data)[33] <- "AUAO" # Non unique column field (UAO) - taken from ARCHER sequencing information and probably not required, may remove in processing
dim(data)
data <- subset(data, wave == 1)
data <- subset(data, AF >= 0.0198)
dim(data)

maf.data <- data.table::data.table(Hugo_Symbol = data$PreferredSymbol, Entrez_Gene_Id = NA, Center = NA, NCBI_Build = "GRCh38", Chromosome = data$chromosome,
                                   Start_Position = data$position, End_Position = (data$position + nchar(as.character(data$mutation))), Strand = '+',
                                   Variant_Classification = data$Variant_Classification, Variant_Type = data$TYPE, Reference_Allele = data$reference,
                                   Tumor_Seq_Allele1 = data$reference, Tumor_Seq_Allele2 = data$mutation, dbSNP_RS = NA, dbSNP_Val_Status  = NA,
                                   Tumor_Sample_Barcode = data$participant_id, Verification_Status = "NA", Mutation_Status = "Somatic",
                                   Sequence_Source = "ArcherDX/Illumina", Validation_Method = "NA", HGVSc = data$base_substitution, HGVSp = data$protein_substitution, 
                                   Consequence = data$consequence, Symbol = data$PreferredSymbol, AF = data$AF)

write.table(maf.data, paste(output_prefix, "non-synonymous.wave1.2PCT_VAF.maf", sep = ""), sep = "\t", row.names = F)

dir.create("./plotting/maftools")

vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                 "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                 "In_Frame_Ins", "Missense_Mutation","Frameshift_Variant")
data = read.maf(maf.data)
getSampleSummary(data)
getGeneSummary(data)
getFields(data)

output_filename <- paste("./plotting/maftools/", output_prefix,"MaftoolsSummary.png", sep = "")
png(output_filename, width = 1400, height = 1000, res = 300, pointsize =5)
plotmafSummary(maf = data, rmOutlier = TRUE, addStat = 'median', textSize = 0.4, dashboard = TRUE, titvRaw = FALSE, showBarcodes = F)
dev.off()

output_filename <- paste("./plotting/maftools/", output_prefix,"Oncoplot.Top30.png", sep = "")
png(output_filename, width = 1400, height = 1000, res = 300, pointsize =5)
oncoplot(maf = data, top = 30, fontSize = 0, showTumorSampleBarcodes= F)
dev.off()

oncostrip(maf = data, genes = c('TP53','JAK2', 'DNMT3A', 'TET2', 'IDH2'), showTumorSampleBarcodes = F)

output_filename <- paste("./plotting/maftools/", output_prefix,"TiTv.png", sep = "")
png(output_filename, width = 1400, height = 1000, res = 300, pointsize =5)
data.titv = titv(maf = data, plot = FALSE, useSyn = TRUE)
plotTiTv(res = data.titv, axisTextSize = c(0.8, 1),showBarcodes = F)
dev.off()

output_filename <- paste("./plotting/maftools/", output_prefix,"Lollipop.JAK2.png", sep = "")
png(output_filename, width = 1400, height = 1000, res = 300, pointsize =5)
lollipopPlot(maf = data, gene = 'JAK2', AACol = 'HGVSp', showMutationRate = TRUE,collapsePosLabel = T, cBioPortal = TRUE, labelPos = c("617"))
dev.off()

output_filename <- paste("./plotting/maftools/", output_prefix,"Lollipop.DNMT3A.png", sep = "")
png(output_filename, width = 1400, height = 1000, res = 300, pointsize =5)
lollipopPlot(maf = data, gene = 'DNMT3A', AACol = 'HGVSp', showMutationRate = TRUE,collapsePosLabel = TRUE, cBioPortal = TRUE)
dev.off()

output_filename <- paste("./plotting/maftools/", output_prefix,"Lollipop.TET2.png", sep = "")
png(output_filename, width = 1400, height = 1000, res = 300, pointsize =5)
lollipopPlot(maf = data, gene = 'TET2', AACol = 'HGVSp', showMutationRate = TRUE,collapsePosLabel = TRUE, cBioPortal = TRUE)
dev.off()

output_filename <- paste("./plotting/maftools/", output_prefix,"Lollipop.ASXL1.png", sep = "")
png(output_filename, width = 1400, height = 1000, res = 300, pointsize =5)
lollipopPlot(maf = data, gene = 'ASXL1', AACol = 'HGVSp', showMutationRate = TRUE,collapsePosLabel = TRUE, cBioPortal = TRUE)
dev.off()

