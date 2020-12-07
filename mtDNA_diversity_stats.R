## calculate nucleotide diversity for mitochondrial haplotypes


library("ape")
library("pegas")
library("adegenet")

apat_align <- read.dna("/Users/idanaughton/Documents/Sanger_Sequences/aphaenogaster_final/Aphaenogaster_all_sequences_CO1_final.fasta", format="fasta")
apat_align

#set populations in alignment
SBA <- apat_align[1:5,]
SCL <- apat_align[6:9,]
GDL <- apat_align[10:11,]
SNI <- apat_align[13:18,]
SCA <- rbind(apat_align[12,],apat_align[19:21,])


#calculate nucleotide diversity per population
SBA_div <- nuc.div(SBA, variance = TRUE)
SBA_div
SCL_div <- nuc.div(SCL, variance = TRUE)
SCL_div
GDL_div <- nuc.div(GDL, variance = TRUE)
GDL_div
SNI_div <- nuc.div(SNI, variance = TRUE)
SNI_div
SCA_div <- nuc.div(SCA, variance = TRUE)
SCA_div

#tajima's D per population
SBA_taj <- tajima.test(SBA)
SBA_taj
SCL_taj <- tajima.test(SCL)
SCL_taj
GDL_taj <- tajima.test(GDL)
GDL_taj
SNI_taj <- tajima.test(SNI)
SNI_taj
SCA_taj <- tajima.test(SCA)
SCA_taj

#Haplotype diversity per population
SBA_hap <- hap.div(SBA, variance = TRUE)
SBA_hap
SCL_hap <- hap.div(SCL, variance = TRUE)
SCL_hap
GDL_hap <- hap.div(GDL, variance = TRUE)
GDL_hap
SNI_hap <- hap.div(SNI, variance = TRUE)
SNI_hap
SCA_hap <- hap.div(SCA, variance = TRUE)
SCA_hap

