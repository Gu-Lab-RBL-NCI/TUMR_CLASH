setwd("~/Desktop/P1 Pri-miRNA scaffold/Analysis cleavage Bioinformatic/")
miRBase <- read.table(file="high_confidence_miRBase21-master.tsv", header = TRUE, fill = TRUE, sep =  " ")
#miRBase <- miRBase[which(miRBase$high_conf=="TRUE"),]

setwd("~/Desktop/R_scripts/Public datasets/CLASH_Helwak/")
clash_seq <- read.table("1-s2.0-S009286741300439X-mmc_clean.txt", header = TRUE)
clash_seq$count <- 1

miRNA_list <- aggregate(clash_seq$count, by=list(clash_seq$miRNA_seq), FUN=sum, na.rm=FALSE)
colnames(miRNA_list) <- c("miRNA_seq","count")

names <- unique(miRBase[,c(1,6)])
colnames(names) <- c("miRNA", "miRNA_seq")
miRNA_list <- merge(names, miRNA_list,by="miRNA_seq", all=TRUE)

miRNA_list <- miRNA_list[order(-miRNA_list$count),]
miRNA_list <- miRNA_list[which(miRNA_list$count>=20),]

miRNA_list$p_seed7 <- NA
miRNA_list$p_last7 <- NA
miRNA_list$p_last7A <- NA 
miRNA_list$b_last7A <- NA 
miRNA_list$p_last7C <- NA 
miRNA_list$p_last7G <- NA 
miRNA_list$p_last7T <- NA 
library(seqinr)
require(stringdist)

i <- 1
for (i in 1:nrow(miRNA_list)) {
  miRNA_seq <- as.character(miRNA_list$miRNA_seq[i])
  data <- clash_seq[which(clash_seq$miRNA_seq==miRNA_seq),]
  background <- clash_seq[which(clash_seq$miRNA_seq!=miRNA_seq),]
  
  rev_compl <- toupper(c2s(rev(comp(s2c(miRNA_seq)))))
  seed7 <- substr(rev_compl,nchar(rev_compl)-7,nchar(rev_compl)-1)
  maxdist <- 1
  last7 <- substr(rev_compl,1,7)
  target_seedless <- gsub("A","[AG]",gsub("C","[CT]",last7))
  
  heptamers_seedless <- expand.grid(n1 = c("A","T","G","C"), n2 = c("A","T","G","C"), n3 = c("A","T","G","C"), n4 = c("A","T","G","C"), n5 = c("A","T","G","C"), n6 = c("A","T","G","C"), n7 = c("A","T","G","C"))
  heptamers_seedless$seedT <-paste0(heptamers_seedless$n1,heptamers_seedless$n2,heptamers_seedless$n3,heptamers_seedless$n4,heptamers_seedless$n5,heptamers_seedless$n6, heptamers_seedless$n7) 
  heptamers_seedless$match <- regexpr(pattern = target_seedless, text = heptamers_seedless$seedT, perl = TRUE)
  
  heptamers_seedless$dist <- stringdist(heptamers_seedless$seedT, last7)
  heptamers_seedless <- heptamers_seedless[which(heptamers_seedless$match!=-1),]
  heptamers_seedless <- heptamers_seedless[which(heptamers_seedless$dist<=maxdist),]
  heptamers_seedless$GC <- 1-(nchar(gsub("C", "", gsub("G", "", heptamers_seedless$seedT)))/nchar(heptamers_seedless$seedT))
  heptamers_seedless$seedTA <- paste0("A",heptamers_seedless$seedT)
  heptamers_seedless$seedTC <- paste0("C",heptamers_seedless$seedT)
  heptamers_seedless$seedTG <- paste0("G",heptamers_seedless$seedT)
  heptamers_seedless$seedTT <- paste0("T",heptamers_seedless$seedT)
  
  toMatch <- unique(heptamers_seedless$seedT)
  toMatchA <- unique(heptamers_seedless$seedTA)
  toMatchC <- unique(heptamers_seedless$seedTC)
  toMatchG <- unique(heptamers_seedless$seedTG)
  toMatchT <- unique(heptamers_seedless$seedTT)
  
  data$seed7 <- gregexpr(seed7, text = data$mRNA_seq_extended, perl = TRUE)
  data$last7 <- gregexpr(paste(toMatch,collapse="|"), text = data$mRNA_seq_extended, perl = TRUE)
  data$last7A <- gregexpr(paste(toMatchA,collapse="|"), text = data$mRNA_seq_extended, perl = TRUE)
  data$last7C <- gregexpr(paste(toMatchC,collapse="|"), text = data$mRNA_seq_extended, perl = TRUE)
  background$last7A <- gregexpr(paste(toMatchA,collapse="|"), text = background$mRNA_seq_extended, perl = TRUE)
  data$last7G <- gregexpr(paste(toMatchG,collapse="|"), text = data$mRNA_seq_extended, perl = TRUE)
  data$last7T <- gregexpr(paste(toMatchT,collapse="|"), text = data$mRNA_seq_extended, perl = TRUE)
  
  p_seed7 <- (nrow(data[which(data$seed7!="-1"),]))/(nrow(data))*100
  p_last7 <- (nrow(data[which(data$last7!="-1"),]))/(nrow(data))*100
  p_last7A <- (nrow(data[which(data$last7A!="-1"),]))/(nrow(data))*100
  b_last7A <- (nrow(background[which(background$last7A!="-1"),]))/(nrow(background))*100
  p_last7C <- (nrow(data[which(data$last7C!="-1"),]))/(nrow(data))*100
  p_last7G <- (nrow(data[which(data$last7G!="-1"),]))/(nrow(data))*100
  p_last7T <- (nrow(data[which(data$last7T!="-1"),]))/(nrow(data))*100

  miRNA_list$p_seed7[i] <- p_seed7
  miRNA_list$p_last7[i] <- p_last7
  miRNA_list$p_last7A[i] <- p_last7A 
  miRNA_list$b_last7A[i] <- b_last7A 
  miRNA_list$p_last7C[i] <- p_last7C 
  miRNA_list$p_last7G[i] <- p_last7G 
  miRNA_list$p_last7T[i] <- p_last7T 

  print(i)
}

over50_list <- miRNA_list[which(miRNA_list$count>=25),]
over50_list$frac <- over50_list$p_last7A/over50_list$p_last7
write.table(over50_list, "CLASH_Helwak_TUMR_extended.txt", sep="\t", append = FALSE, row.names = FALSE, col.names = TRUE)

sig2 <- mean(over50_list$b_last7A)+(2*sd(over50_list$b_last7A))
sig3 <- mean(over50_list$b_last7A)+(3*sd(over50_list$b_last7A))
sig4 <- mean(over50_list$b_last7A)+(4*sd(over50_list$b_last7A))
sig6 <- mean(over50_list$b_last7A)+(6*sd(over50_list$b_last7A))

oversig_list <- over50_list[which(over50_list$p_last7A>=sig3),]
oversig_list$notA <- oversig_list$p_last7-oversig_list$p_last7A

mean(oversig_list$p_last7A)
mean(oversig_list$p_last7C)
mean(oversig_list$p_last7G)
mean(oversig_list$p_last7T)
write.table(oversig_list, "CLASH_Helwak_TUMR_overnoise.txt", sep="\t", append = FALSE, row.names = FALSE, col.names = TRUE)

