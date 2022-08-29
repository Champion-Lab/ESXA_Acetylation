# Document Created: 2022-05-20
# By: Simon D. Weaver (sweaver4@nd.edu)
# Data from Daniel Hu
# Lab: Matt Champion
# University of Notre Dame, Dept Chem & Biochem

# Title: ESXA_acetylation Data Analysis

# This is an example workflow for one biological replicate (November)
# Lines for other biological replicates have been commented out.

# Workflow analysis for evaluating the N-terminal protein acetylation from 3 
# sample types, 3 biological replicates with 3 triplicate analyses for each
# sample type/bio replicate combination. Each analysed by 2 different search
# methods in peaks: Database (standard) and LFQ


#libraries
library(stringr)
library(seqinr)
library(dplyr)
library(ggplot2)
library(tidyr)

#set working directory:
#setwd("~/04_Champion_Lab/02_N-terminal_Acetylation/Data Analysis/2022_05_20")
setwd("~/04_Champion_Lab/02_N-terminal_Acetylation/Data Analysis/For Publication Github 2022_05_24/ESXA_Acetylation")

# Data organization:
#   -3 biological replicates, denoted by their month: Sep, Nov, Nov.
#   -3 conditions: WT, Del, Complement
#   -3 technical replicates
#   -2 search techniques (Normal database and LFQ)

#source functions
source("functions/functions.R")


#IMPORTS
#import protein search results for normalization
#protein_Sep_df <- read.csv(file = "search_results/September lfq.proteins.csv")
#protein_Oct_df <- read.csv(file = "search_results/October lfq.proteins.csv")
protein_Nov_df <- read.csv(file = "search_results/November lfq.proteins.csv")

#import each database search result
#peptides_Sep_db <- read.csv(file = "search_results/Sep db.peptides.csv")
#peptides_Oct_db <- read.csv(file = "search_results/Oct db.peptides.csv")
peptides_Nov_db <- read.csv(file = "search_results/Nov db.peptides.csv")

#import database
db <- read.fasta(file = "database/10Feb2022_Mycobrowser_Mmarinum_v4.fasta",
                 seqtype = "AA",
                 as.string = T,
                 whole.header = F)

#create lists of each sample type
mos <- c("Sep", "Oct", "Nov")
samples <- c("WT", "Del", "Comp")
sample_names_27 <- list()
for (mo in mos){
  for (sample in samples) {
    for (rep in 1:3){
      sample_names_27[[length(sample_names_27)+1]] <- {
        paste0(mo, sample, "_", rep)
      }
    }
  }
}
#sample_names_Sep <- sample_names_27[1:9]
#sample_names_Oct <- sample_names_27[10:18]
sample_names_Nov <- sample_names_27[19:27]


#Pull out the protein area information for RpoA for each condition for norm.
rpoA_Sep <- protein_Sep_df[grepl(".*rpoA.*", protein_Sep_df$Accession),]
rpoA_Oct <- protein_Oct_df[grepl(".*rpoA.*", protein_Oct_df$Accession),]
rpoA_Nov <- protein_Nov_df[grepl(".*rpoA.*", protein_Nov_df$Accession),]

#save each rpoA area value in a list for easy access
# rpoA_areas <- list()
# rpoA_areas[["SepWT_1"]] <- rpoA_Sep$SepWT_1.Area[1]
# rpoA_areas[["SepWT_2"]] <- rpoA_Sep$SepWT_2.Area[1]
# rpoA_areas[["SepWT_3"]] <- rpoA_Sep$SepWT_3.Area[1]
# rpoA_areas[["SepDel_1"]] <- rpoA_Sep$SepDel_1.Area[1]
# rpoA_areas[["SepDel_2"]] <- rpoA_Sep$SepDel_2.Area[1]
# rpoA_areas[["SepDel_3"]] <- rpoA_Sep$SepDel_3.Area[1]
# rpoA_areas[["SepComp_1"]] <- rpoA_Sep$SepComp_1.Area[1]
# rpoA_areas[["SepComp_2"]] <- rpoA_Sep$SepComp_2.Area[1]
# rpoA_areas[["SepComp_3"]] <- rpoA_Sep$SepComp_3.Area[1]
# rpoA_areas[["OctWT_1"]] <- rpoA_Oct$OctWT_1.Area[1]
# rpoA_areas[["OctWT_2"]] <- rpoA_Oct$OctWT_2.Area[1]
# rpoA_areas[["OctWT_3"]] <- rpoA_Oct$OctWT_3.Area[1]
# rpoA_areas[["OctDel_1"]] <- rpoA_Oct$OctDel_1.Area[1]
# rpoA_areas[["OctDel_2"]] <- rpoA_Oct$OctDel_2.Area[1]
# rpoA_areas[["OctDel_3"]] <- rpoA_Oct$OctDel_3.Area[1]
# rpoA_areas[["OctComp_1"]] <- rpoA_Oct$OctComp_1.Area[1]
# rpoA_areas[["OctComp_2"]] <- rpoA_Oct$OctComp_2.Area[1]
# rpoA_areas[["OctComp_3"]] <- rpoA_Oct$OctComp_3.Area[1]
rpoA_areas[["NovWT_1"]] <- rpoA_Nov$NovWT_1.Area[1]
rpoA_areas[["NovWT_2"]] <- rpoA_Nov$NovWT_2.Area[1]
rpoA_areas[["NovWT_3"]] <- rpoA_Nov$NovWT_3.Area[1]
rpoA_areas[["NovDel_1"]] <- rpoA_Nov$NovDel_1.Area[1]
rpoA_areas[["NovDel_2"]] <- rpoA_Nov$NovDel_2.Area[1]
rpoA_areas[["NovDel_3"]] <- rpoA_Nov$NovDel_3.Area[1]
rpoA_areas[["NovComp_1"]] <- rpoA_Nov$NovComp_1.Area[1]
rpoA_areas[["NovComp_2"]] <- rpoA_Nov$NovComp_2.Area[1]
rpoA_areas[["NovComp_3"]] <- rpoA_Nov$NovComp_3.Area[1]
#turn the list into a dataframe
rpoA_areas_df <- data.frame(sample = attr(rpoA_areas, "names"),
                            rpoA_areas = unlist(rpoA_areas))


#Run the database analysis for Sep, Oct, and Nov
#source("September_replicate.R")
#source("October_replicate.R")
source("November_replicate.R")
