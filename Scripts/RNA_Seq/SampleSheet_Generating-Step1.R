#!/usr/local/bin/Rscript
# CAD_mt709_0001: Monocyte, Macrophage, bone marrow monocyte
# RNASeq for mouse with diet treatment
#
# Analysis Performed by Xiaohui Zhao
# School of Medicine, University of Cambridge, UK
# Copyright Xiaohui Zhao (xz289@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

message("+-------           Basic settings                           --------+")

library("readxl")
## setup working directory
Basedir <- Basedir
setwd(Basedir) 

message("+-------------------------------------------------------------------+")
message("+   Generate the full sample sheet for second batch, Macrophage     +")
message("+-------------------------------------------------------------------+")

SamTab     <- read.csv("./Data/First_batch/SLX-17236.H533HBBXY.s_6.contents.csv", header=T)
SamTab     <- SamTab[order(SamTab$Barcode),]
Condition  <- rep(c("iWD", "cWD"), each=3)
Replicates <- rep(c(1:3), length=6)
Cond       <- paste0(Condition, "_REP", Replicates)
fastq_1    <- list.files(path="./Data/First_batch", pattern="*.r_1.fq.gz")
fastq_2    <- ""
seqID_1    <- unlist(lapply(fastq_1, function(x) strsplit(x, split="[.]")[[1]][2]))
seqID_1    <- gsub("_", "-", seqID_1)
f1s        <- cbind(fastq_1, seqID_1); 
fastq.files<- data.frame(fastq_1, seqID_1)
colnames(fastq.files) <- c("fastq_1", "Barcode")

SamTall    <- cbind(SamTab, Replicates, Condition, Cond)
SamTable   <- merge(SamTall, fastq.files, by = "Barcode")


write.csv(SamTable, 
          file = "./Data/First_batch/Macrophage_FirstBatch_SampleTable.csv", 
          row.names = F, quote = T)


message("+-----        Generate the sample sheet for Nextflow input     -----+")

sample     <- SamTable$Cond
fastq_1    <- SamTable$fastq_1
fastq_2    <- ""
strandedness <- "unstranded"

nfTable    <- data.frame(sample, fastq_1, fastq_2, strandedness)
write.csv(nfTable, 
          file = "./Data/First_batch/Macrophage_FirstBatch_Nextflow_SampleTable.csv", 
          row.names = F)

message("+-------------------------------------------------------------------+")
message("+   Generate the full sample sheet for second batch, Macrophage     +")
message("+-------------------------------------------------------------------+")


SamTab     <- read.csv("./Data/Second_batch/SLX-18681.HGVWJBBXY.s_5.contents.csv", header=T)
Condition  <- rep(c("cWD", "iWD"), each=3)
Replicates <- rep(c(1:3), length=6)
Cond       <- paste0(Condition, "_REP", Replicates)
fastq_1    <- list.files(path="./Data/Second_batch/", pattern="*.r_1.fq.gz")
fastq_2    <- ""
seqID_1    <- unlist(lapply(fastq_1, function(x) strsplit(x, split="[.]")[[1]][2]))
seqID_1    <- gsub("_", "-", seqID_1)
f1s        <- cbind(fastq_1, seqID_1); 
fastq.files<- data.frame(fastq_1, seqID_1)
colnames(fastq.files) <- c("fastq_1", "Barcode")

SamTall    <- cbind(SamTab, Replicates, Condition, Cond)
SamTable   <- merge(SamTall, fastq.files, by = "Barcode")


write.csv(SamTable, 
          file = "./Data/Second_batch/Macrophage_SecondBatch_SampleTable.csv", 
          row.names = F, quote = T)

message("+-----  Generate the sample sheet for Nextflow input           -----+")

sample     <- SamTable$Cond
fastq_1    <- SamTable$fastq_1
fastq_2    <- ""
strandedness <- "unstranded"

nfTable    <- data.frame(sample, fastq_1, fastq_2, strandedness)

write.csv(nfTable, 
          file = "./Data/Second_batch/Macrophage_SecondBatch_Nextflow_SampleTable.csv",
          row.names = F)

message("+-------------------------------------------------------------------+")
message("+        Generate the sample sheet for cWD Nrp1+/+ vs Nrp1 -/-      +")
message("+-------------------------------------------------------------------+")


cruk <- read.csv("./Data/SLX-22500/SLX-22500.HNLLVDRX2.s_2.contents.csv", header=T)
samT <- read_excel("./Data/SLX-22500/GTC280_TruSeq\ RNA\ UDI.xlsx", sheet = 2)
colnames(samT)[1] <- c("Sample.name")
samTable <- merge(samT, cruk, by = "Sample.name")
samTable <- samTable[order(samTable$Barcode), ]

fastq1   <- list.files(path = ".", pattern = "*.r1.fq.gz")
fastq2   <- list.files(path = ".", pattern = "*.r2.fq.gz")
samTable$fastq_1 <- fastq1
samTable$fastq_2 <- fastq2
samTable$sample  <- paste0(samTable$Group2, "_REP", samTable$Replicates)
samTable$strandedness <- samTable$Strand

write.csv(samTable, file = "./Data/SLX-22500/NRP1_SampleTable_cWD.csv", 
          row.names = F, quote = T)
write.csv(samTable[,c("sample", "fastq_1", "fastq_2", "strandedness")],
          file = "./Data/SLX-22500/NRP1_Nextflow_SampleTable_cWD.csv", 
          row.names = F)

##------------------------ Finish---------------------------------------------##