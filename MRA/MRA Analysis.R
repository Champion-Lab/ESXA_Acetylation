
#libraries
library(stringr)
library(seqinr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(ggrepel)

set.seed(1839)

#set working directory:
#setwd("~/04_Champion_Lab/02_N-terminal_Acetylation/Data Analysis/2022_05_20")
setwd("~/04_Champion_Lab/02_N-terminal_Acetylation/MRA/github/ESXA_Acetylation/")

#source functions
source("functions/functions.R")

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

sample_names_Nov <- sample_names_27[19:27]


maxMedian <- NA
minMin <- NA
for (samp in sample_names_Nov) {
  old_col <- paste0("Area.", samp)
  median <- median(peptides_Nov_db[[old_col]], na.rm = T)
  min <- min(peptides_Nov_db[[old_col]], na.rm = T)
  print(min)
  if (is.na(maxMedian)) {
    maxMedian <- median
  } else if (maxMedian < median){
    maxMedian <- median
  }
  if (is.na(minMin)) {
    minMin <- min
  } else if (minMin > min) {
    minMin <- min
  }
}


areas <- peptides_Nov_db %>% select(matches("Area\\.")) %>%
  pivot_longer(1:9) %>% group_by(name) %>%
  summarise(minimum = min(value, na.rm = T))

avgMin <- mean(areas$minimum)
sdMin <- sd(areas$minimum)


#custom median normalization
for (samp in sample_names_Nov) { #loops through all the sample names in Nov rep
  #create a string to create column 
  new_col <- paste0(samp, "_norm_area")
  #create a string to reference column
  old_col <- paste0("Area.", samp) 
  median <- median(peptides_Nov_db[[old_col]], na.rm = T)
  #set the value of the new column to the old column / rpoA area
  peptides_Nov_db[[new_col]] <- (peptides_Nov_db[[old_col]] - median) + maxMedian
}

#remove unused columns
columns_to_remove1 <- names(peptides_Nov_db[,grepl("Intensity",names(peptides_Nov_db))])
columns_to_remove2 <- names(peptides_Nov_db[,grepl("X.Spec",names(peptides_Nov_db))])
peptides_Nov_db <- peptides_Nov_db[ , -which(names(peptides_Nov_db) %in% c(columns_to_remove1, columns_to_remove2))]


#count the number of WT and Comp values that are not NA
peptides_Nov_db$WT_comp_count <- 0
temp_db <- peptides_Nov_db[,grepl("Nov[WC].*_norm_area",
                                  names(peptides_Nov_db))]
for (i in 1:nrow(peptides_Nov_db)) {
  peptides_Nov_db$WT_comp_count[i] <- rowSums(!is.na(temp_db[i,]))
}

#remove rows that do not contain at least 4 datapoints between Comp and WT
peptides_Nov_db <- peptides_Nov_db[peptides_Nov_db$WT_comp_count > 3,]


#create sample types for November
Nov_sample_avs <- c("NovWT", "NovDel", "NovComp")

#loop through each sample type and row and calculate an average
for (samp in Nov_sample_avs) {
  #create a new column name
  new_col <- paste0(samp, "_avg")
  new_col2 <- paste0(samp, "_stdev")
  new_col3 <- paste0(samp, "_rsd")
  #create the column and set it to NA
  peptides_Nov_db[[new_col]] <- NA
  peptides_Nov_db[[new_col2]] <- NA
  peptides_Nov_db[[new_col3]] <- NA
  #create a string to pull old column names
  old_cols <- paste0(samp, "_[1-3]_norm_area")
  #create a temporary db with only values from the 3 replicates of interests
  temp_db <- peptides_Nov_db[,grepl(old_cols,
                                    names(peptides_Nov_db))]
  #loop through each row
  for (i in 1:nrow(peptides_Nov_db)) {
    #create a vector of the replicate values from that line
    reps <- c(temp_db[i,1], temp_db[i,2], temp_db[i,3])
    #store the average of the reps in the new column, removing NA values
    average <- mean(reps, na.rm = T)
    peptides_Nov_db[[new_col]][i] <- average
    peptides_Nov_db[[new_col2]][i] <- sd(reps, na.rm = T)
    peptides_Nov_db[[new_col3]][i] <- average / (sd(reps, na.rm = T))
  }
}



#create a column that is just the amino acid sequence
peptides_Nov_db$sequence <- str_remove_all(peptides_Nov_db$Peptide, "[a-z1-9()+-:.]")

#create column for accession number
peptides_Nov_db$Accession_only <- ""
for (i in 1:nrow(peptides_Nov_db)) {
  peptides_Nov_db$Accession_only[i] <-  {
    str_split(peptides_Nov_db$Accession[i], "\\|")[[1]][1]
  }
}

#check for N-term Acetylation from the PTM column
peptides_Nov_db$Nterm_acet <- grepl(".*Acetylation \\(Protein N\\-term\\).*",
                                    peptides_Nov_db$PTM)




#return the start position of each peptide based on amino acid sequence in the database
peptides_Nov_db$start_position <- NA
for (i in 1:nrow(peptides_Nov_db)) {
  peptides_Nov_db$start_position[i] <- return_start_position(peptide = peptides_Nov_db$sequence[i],
                                                             accession = peptides_Nov_db$Accession_only[i],
                                                             database = db)
}

#remove peptides that start after amino acid 2
peptides_Nov_db <- peptides_Nov_db[peptides_Nov_db$start_position <= 2,]

#remove peptides that don't have a valid accession number
peptides_Nov_db <- peptides_Nov_db[!is.na(peptides_Nov_db$Accession_only),]

#set NA values to 0, only Del will have NA values
peptides_Nov_db$NovDel_avg[is.na(peptides_Nov_db$NovDel_avg)] <- 0

#count how many values have RSDs and print out the statistics
count_useable_nterms <- nrow(peptides_Nov_db)
WT_sd_nas <- sum(is.na(peptides_Nov_db$NovWT_stdev))
Del_sd_nas <- sum(is.na(peptides_Nov_db$NovDel_stdev))
Comp_sd_nas <- sum(is.na(peptides_Nov_db$NovComp_stdev))
print("Percent w/o std dev Nov:")
print(paste0("WT: ", round((WT_sd_nas/count_useable_nterms)*100, 2), "%"))
print(paste0("Del: ", round((Del_sd_nas/count_useable_nterms)*100, 2), "%"))
print(paste0("Comp: ", round((Comp_sd_nas/count_useable_nterms)*100, 2), "%"))

#calculate the average RSD for each sample type
av_rsd_WT <- mean(peptides_Nov_db$NovWT_rsd, na.rm = T)
av_rsd_Del <- mean(peptides_Nov_db$NovDel_rsd, na.rm = T)
av_rsd_Comp <- mean(peptides_Nov_db$NovComp_rsd, na.rm = T)

#set the RSD value to the average RSD value if it is NA
peptides_Nov_db$NovWT_rsd[is.na(peptides_Nov_db$NovWT_rsd)] <- av_rsd_WT
peptides_Nov_db$NovDel_rsd[is.na(peptides_Nov_db$NovDel_rsd)] <- av_rsd_Del
peptides_Nov_db$NovComp_rsd[is.na(peptides_Nov_db$NovComp_rsd)] <- av_rsd_Comp


#back calculate SD if SD is NA from the average RSD
for (i in 1:nrow(peptides_Nov_db)) {
  if (is.na(peptides_Nov_db$NovWT_stdev[i])) {
    stdev <- peptides_Nov_db$NovWT_avg[i] * (peptides_Nov_db$NovWT_rsd[i] / 100)
    peptides_Nov_db$NovWT_stdev[i] <- stdev
  }
  if (is.na(peptides_Nov_db$NovDel_stdev[i])) {
    stdev <- peptides_Nov_db$NovDel_avg[i] * (peptides_Nov_db$NovDel_rsd[i] / 100)
    peptides_Nov_db$NovDel_stdev[i] <- stdev
  }
  if (is.na(peptides_Nov_db$NovComp_stdev[i])) {
    stdev <- peptides_Nov_db$NovComp_avg[i] * (peptides_Nov_db$NovComp_rsd[i] / 100)
    peptides_Nov_db$NovComp_stdev[i] <- stdev
  }
}

#keep only the peptides with an acetylated N term
keeps_Nov <- peptides_Nov_db[peptides_Nov_db$Nterm_acet == T,]

#calculate a ratio of Del/WT and Del/Comp
keeps_Nov$Del_WT_ratio <- keeps_Nov$NovDel_avg / keeps_Nov$NovWT_avg
keeps_Nov$Del_Comp_ratio <- keeps_Nov$NovDel_avg / keeps_Nov$NovComp_avg

#Keep only those peptides where both ratios are less than 1
keeps_Nov <- keeps_Nov[keeps_Nov$Del_WT_ratio < 1 & keeps_Nov$Del_Comp_ratio < 1,]

#set a column to track whether the Del value is 0 or not
keeps_Nov$Del_zero <- NA
for (i in 1:nrow(keeps_Nov)) {
  if (keeps_Nov$NovDel_avg[i] == 0) {
    keeps_Nov$Del_zero[i] <- T
  } else {
    keeps_Nov$Del_zero[i] <- F
  }
}


#same analysis as above, but with the non acetylated versions
keeps_Nov_noA <- peptides_Nov_db[peptides_Nov_db$Nterm_acet == F,]

keeps_Nov_noA$Del_WT_ratio <- keeps_Nov_noA$NovDel_avg / keeps_Nov_noA$NovWT_avg
keeps_Nov_noA$Del_Comp_ratio <- keeps_Nov_noA$NovDel_avg / keeps_Nov_noA$NovComp_avg

#remove this filter step for non-acetylated
#keeps_Nov_noA <- keeps_Nov_noA[keeps_Nov_noA$Del_WT_ratio < 1 & keeps_Nov_noA$Del_Comp_ratio < 1,]

keeps_Nov_noA$Del_zero <- NA
for (i in 1:nrow(keeps_Nov_noA)) {
  if (keeps_Nov_noA$NovDel_avg[i] == 0) {
    keeps_Nov_noA$Del_zero[i] <- T
  } else {
    keeps_Nov_noA$Del_zero[i] <- F
  }
}

#create a new column on which to combine acet/nonacet. Take the peptide seq from
#non acet, and remove the nterm acet from the acet version
keeps_Nov_noA$Peptide_to_compare <- keeps_Nov_noA$Peptide
keeps_Nov$Peptide_to_compare <- remove_nacetyl_mod(keeps_Nov$Peptide)

#join together acet and nonacet dfs, keeping all rows
combined_Nov <- full_join(keeps_Nov, keeps_Nov_noA, by = "Peptide_to_compare",
                          suffix = c(".acet", ".non_A"))

#set NA values in sd and average to 0
combined_Nov$NovWT_avg.acet[is.na(combined_Nov$NovWT_avg.acet)] <- 0
combined_Nov$NovDel_avg.acet[is.na(combined_Nov$NovDel_avg.acet)] <- 0
combined_Nov$NovComp_avg.acet[is.na(combined_Nov$NovComp_avg.acet)] <- 0
combined_Nov$NovWT_avg.non_A[is.na(combined_Nov$NovWT_avg.non_A)] <- 0
combined_Nov$NovDel_avg.non_A[is.na(combined_Nov$NovDel_avg.non_A)] <- 0
combined_Nov$NovComp_avg.non_A[is.na(combined_Nov$NovComp_avg.non_A)] <- 0
combined_Nov$NovWT_stdev.acet[is.na(combined_Nov$NovWT_stdev.acet)] <- 0
combined_Nov$NovDel_stdev.acet[is.na(combined_Nov$NovDel_stdev.acet)] <- 0
combined_Nov$NovComp_stdev.acet[is.na(combined_Nov$NovComp_stdev.acet)] <- 0
combined_Nov$NovWT_stdev.non_A[is.na(combined_Nov$NovWT_stdev.non_A)] <- 0
combined_Nov$NovDel_stdev.non_A[is.na(combined_Nov$NovDel_stdev.non_A)] <- 0
combined_Nov$NovComp_stdev.non_A[is.na(combined_Nov$NovComp_stdev.non_A)] <- 0

#combine the areas of acet and non acet
combined_Nov$total_area_WT_Nov <- combined_Nov$NovWT_avg.acet + combined_Nov$NovWT_avg.non_A
combined_Nov$total_area_Del_Nov <- combined_Nov$NovDel_avg.acet + combined_Nov$NovDel_avg.non_A
combined_Nov$total_area_Comp_Nov <- combined_Nov$NovComp_avg.acet + combined_Nov$NovComp_avg.non_A

#propogate error
combined_Nov$total_area_WT_sd_Nov <- sqrt((combined_Nov$NovWT_stdev.acet^2) + (combined_Nov$NovWT_stdev.non_A^2))
combined_Nov$total_area_Del_sd_Nov <- sqrt((combined_Nov$NovDel_stdev.acet^2) + (combined_Nov$NovDel_stdev.non_A^2))
combined_Nov$total_area_Comp_sd_Nov <- sqrt((combined_Nov$NovComp_stdev.acet^2) + (combined_Nov$NovComp_stdev.non_A^2))

#group peptides with custom function, basically adds up all the values by accession ID
#see functions.R
final_Nov <- group_peptides(combined_Nov, "Nov")

#calculate all RSDs, fill NAs with 0
final_Nov$WT_rsd_acet <- (final_Nov$WT_sd_acet / final_Nov$WT_av_acet) * 100
final_Nov$WT_rsd_acet[is.nan(final_Nov$WT_rsd_acet)] <- 0
final_Nov$Del_rsd_acet <- (final_Nov$Del_sd_acet / final_Nov$Del_av_acet) * 100
final_Nov$Del_rsd_acet[is.nan(final_Nov$Del_rsd_acet)] <- 0
final_Nov$Comp_rsd_acet <- (final_Nov$Comp_sd_acet / final_Nov$Comp_av_acet) * 100
final_Nov$Comp_rsd_acet[is.nan(final_Nov$Comp_rsd_acet)] <- 0
final_Nov$WT_rsd_non_A <- (final_Nov$WT_sd_non_A / final_Nov$WT_av_non_A) * 100
final_Nov$WT_rsd_non_A[is.nan(final_Nov$WT_rsd_non_A)] <- 0
final_Nov$Del_rsd_non_A <- (final_Nov$Del_sd_non_A / final_Nov$Del_av_non_A) * 100
final_Nov$Del_rsd_non_A[is.nan(final_Nov$Del_rsd_non_A)] <- 0
final_Nov$Comp_rsd_non_A <- (final_Nov$Comp_sd_non_A / final_Nov$Comp_av_non_A) * 100
final_Nov$Comp_rsd_non_A[is.nan(final_Nov$Comp_rsd_non_A)] <- 0
final_Nov$Total_area_WT_rsd <- (final_Nov$total_area_WT_sd / final_Nov$total_area_WT) * 100
final_Nov$Total_area_WT_rsd[is.nan(final_Nov$Total_area_WT_rsd)] <- 0
final_Nov$Total_area_Del_rsd <- (final_Nov$total_area_Del_sd / final_Nov$total_area_Del) * 100
final_Nov$Total_area_Del_rsd[is.nan(final_Nov$Total_area_Del_rsd)] <- 0
final_Nov$Total_area_Comp_rsd <- (final_Nov$total_area_Comp_sd / final_Nov$total_area_Comp) * 100
final_Nov$Total_area_Comp_rsd[is.nan(final_Nov$Total_area_Comp_rsd)] <- 0

#calculate percent acetylated along with error propogation by sample
final_Nov$percent_acet_WT <- (final_Nov$WT_av_acet / final_Nov$total_area_WT) *100
final_Nov$percent_acet_WT_rsd <- NA
for (i in 1:nrow(final_Nov)) {
  final_Nov$percent_acet_WT_rsd[i] <- sqrt((final_Nov$WT_rsd_acet[i]) + (final_Nov$Total_area_WT_rsd[i])^2)
}
final_Nov$percent_acet_WT_sd <- (final_Nov$percent_acet_WT_rsd / 100) * final_Nov$percent_acet_WT

final_Nov$percent_acet_Del <- (final_Nov$Del_av_acet / final_Nov$total_area_Del) *100
final_Nov$percent_acet_Del_rsd <- NA
for (i in 1:nrow(final_Nov)) {
  final_Nov$percent_acet_Del_rsd[i] <- sqrt((final_Nov$Del_rsd_acet[i]) + (final_Nov$Total_area_Del_rsd[i])^2)
}
final_Nov$percent_acet_Del_sd <- (final_Nov$percent_acet_Del_rsd / 100) * final_Nov$percent_acet_Del

final_Nov$percent_acet_Comp <- (final_Nov$Comp_av_acet / final_Nov$total_area_Comp) *100
final_Nov$percent_acet_Comp_rsd <- NA
for (i in 1:nrow(final_Nov)) {
  final_Nov$percent_acet_Comp_rsd[i] <- sqrt((final_Nov$Comp_rsd_acet[i]) + (final_Nov$Total_area_Comp_rsd[i])^2)
}
final_Nov$percent_acet_Comp_sd <- (final_Nov$percent_acet_Comp_rsd / 100) * final_Nov$percent_acet_Comp

library(flexclust)
library(factoextra)


#clustering with ratios
Nov_cluster_acet <- data.frame(protein = final_Nov$accession,
                               WT = final_Nov$WT_av_acet,
                               Del = final_Nov$Del_av_acet,
                               Comp = final_Nov$Comp_av_acet)
row.names(Nov_cluster_acet) <- Nov_cluster_acet$protein
Nov_cluster_acet <- subset(Nov_cluster_acet, select = -(protein))
Nov_cluster_acet <- Nov_cluster_acet[rowSums(Nov_cluster_acet) > 0,]

Nov_cluster_acet_temp <- Nov_cluster_acet

Nov_cluster_acet$Del_Comp_rat <- Nov_cluster_acet$Del / Nov_cluster_acet$Comp
Nov_cluster_acet$Del_WT_rat <- Nov_cluster_acet$Del / Nov_cluster_acet$WT
Nov_cluster_acet$Comp_WT_rat <- Nov_cluster_acet$Comp / Nov_cluster_acet$WT

Nov_cluster_acet <- subset(Nov_cluster_acet, select = -c(WT, Del, Comp))
Nov_cluster_acet <- as.data.frame(scale(Nov_cluster_acet, center = F))

#creates the clustered object
#Nov_cluster_acet is the name of the dataframe with data, 2 is number of clusters
# 25 is the number of starting positions (higher nstart takes longer but produces more
# reliable results)
Nov_kmeans_acet <- kmeans(Nov_cluster_acet, 2, nstart = 25)
#print the clustered object
Nov_kmeans_acet


#add the clusters to the original dataframe, save as new object
Nov_cluster_acet_4 <- cbind(Nov_cluster_acet_temp, cluster = Nov_kmeans_acet$cluster)

#visualize the clusters (basically creates a PCA plot)
clusterplot <- fviz_cluster(Nov_kmeans_acet, Nov_cluster_acet, ggtheme = theme_bw(), geom = "point",
                            ellipse = T, main = F) +
  theme_bw(base_size = 20) +
  scale_color_manual(values = c("darkgoldenrod1", "magenta", "blue")) +
  scale_fill_manual(values = c("darkgoldenrod1", "magenta", "blue"))
#view the plot
clusterplot

#everything before the theme_bw line
clusterplot_to_save <- clusterplot <- fviz_cluster(Nov_kmeans_acet, Nov_cluster_acet, ggtheme = theme_bw(), geom = "point",
                                                   ellipse = T, main = F)

Nov_cluster_acet_4$protein <- row.names(Nov_cluster_acet_4)
Nov_cluster_acet_4_long <- pivot_longer(Nov_cluster_acet_4, c(WT, Del, Comp), names_to = "condition", values_to = "area")
Nov_cluster_acet_4_long$cluster <- as.factor(Nov_cluster_acet_4_long$cluster)

Nov_cluster_acet_4_long$area[Nov_cluster_acet_4_long$area == 0] <- 0.00001

write.csv(Nov_cluster_acet_4_long, "MRA/clusterexports.csv", row.names = F)








nov <- read.csv("MRA/clusterexports.csv")

nov$cluster <- as.character(nov$cluster)
nov$areaLog <- log(nov$area, base= 10)

for (i in 1:nrow(nov)) {
  if (nov$areaLog[i] < 0) {
    nov$areaLog[i] <- log(rnorm(1, mean = avgMin, sd = sdMin),base = 10)
  }
}

nov$condition <- factor(nov$condition, levels = c("WT", "Del", "Comp"))

labs <- c("EsxA Cluster", "Other")
names(labs) <- c("2", "1")

nov$cluster <- factor(nov$cluster, levels = c("2", "1"))

nov$category2[nov$cluster == "2"] <- "EsxA Cluster"
nov$category2[nov$cluster != "2"] <- "Other"

nov$label[nov$protein == "MMAR_5450"] <- "EsxA"

ggplot(nov, aes(x = condition, y = areaLog)) +
  geom_jitter((aes(color = cluster)),
              height = 0, width = 0.2) +
  theme_bw(base_size = 12) +
  geom_boxplot(aes(color = cluster, alpha = 0)) +
  labs(y = "N-terminal Acetylated Peptide\n Median Normalized Area\n(log10 transformed)",
       x = element_blank()) +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  facet_wrap("cluster", labeller = labeller(cluster = labs)) +
  # stat_compare_means(method = "anova") +
  scale_x_discrete(labels = c("Wild Type", "\u0394emp1", "\u0394emp1/Complement"))+
  stat_compare_means(method = "t.test",
                     comparisons = list(c("WT", "Comp"), c("WT", "Del"), c("Del", "Comp")),
                     label = "p.signif",tip.length = 0.01,step.increase = 0.1) +
  geom_text_repel(aes(label = label),nudge_y = 1,
                  nudge_x = 1, position = position_jitter(seed = 1)) +
  coord_cartesian(ylim = c(1.2,6.2))

ggsave("MRA/pairwise_ttest_analysis_of_EXSA_cluster.png")

ggplot(nov, aes(x = condition, y = areaLog)) +
  geom_jitter((aes(color = cluster)),
              height = 0, width = 0.2) +
  theme_bw(base_size = 12) +
  geom_boxplot(aes(color = cluster, alpha = 0)) +
  labs(y = "N-terminal Acetylated Peptide\n Median Normalized Area\n(log10 transformed)",
       x = element_blank()) +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  facet_wrap("cluster", labeller = labeller(cluster = labs)) +
  stat_compare_means(method = "anova") +
  scale_x_discrete(labels = c("Wild Type", "\u0394emp1", "\u0394emp1/Complement"))+
  # stat_compare_means(method = "t.test",
  #                    comparisons = list(c("WT", "Comp"), c("WT", "Del"), c("Del", "Comp")),
  #                    label = "p.signif",tip.length = 0.01,step.increase = 0.1) +
  coord_cartesian(ylim = c(1,5.5))

ggsave("MRA/anova_analysis_of_EXSA_cluster.png")





summary <- nov %>% group_by(cluster, condition) %>%
  summarise(n())
summary












ggplot(nov, aes(x = condition, y = areaLog)) +
  geom_point((aes(color = cluster)),
             position = position_jitter(width = 0.3, seed = 1)) +
  theme_bw(base_size = 12) +
  geom_boxplot(aes(color = cluster, alpha = 0)) +
  labs(y = "N-terminal Acetylated Peptide\n Median Normalized Area\n(log10 transformed)",
       x = element_blank()) +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  facet_wrap("cluster", labeller = labeller(cluster = labs)) +
  # stat_compare_means(method = "anova") +
  scale_x_discrete(labels = c("WT",
                              expression(paste("\u0394",italic("emp1"))),
                              expression(paste("\u0394",italic("emp1"), "/comp"))))+
  stat_compare_means(method = "t.test",
                     comparisons = list(c("WT", "Comp"), c("WT", "Del"), c("Del", "Comp")),
                     label = "p.signif",tip.length = 0.01,step.increase = 0.1) +
  geom_text_repel(aes(label = label),box.padding = 1.5,
                  position = position_jitter(width = 0.3,seed = 1)) +
  coord_cartesian(ylim = c(1.2,6.2))

ggsave("MRA/pairwise_ttest_analysis_of_EXSA_cluster.png")






ggplot(nov, aes(x = condition, y = areaLog)) +
  geom_point((aes(color = cluster)),
             position = position_jitter(width = 0.3, seed = 1)) +
  theme_bw(base_size = 12) +
  geom_boxplot(aes(color = cluster, alpha = 0)) +
  labs(y = "N-terminal Acetylated Peptide\n Median Normalized Area\n(log10 transformed)",
       x = element_blank()) +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  facet_wrap("cluster", labeller = labeller(cluster = labs)) +
  scale_x_discrete(labels = c("WT",
                              expression(paste("\u0394",italic("emp1"))),
                              expression(paste("\u0394",italic("emp1"), "/comp"))))+
  stat_compare_means(method = "anova") +
  stat_compare_means(method = "t.test",
                     comparisons = list(c("WT", "Comp"), c("WT", "Del"), c("Del", "Comp")),
                     label = "p.signif",tip.length = 0.01,step.increase = 0.1) +
  geom_text_repel(aes(label = label),box.padding = 1.5,
                  position = position_jitter(width = 0.3,seed = 1)) +
  coord_cartesian(ylim = c(1.2,6.2))



#tukey test following anova



stat.test.esxAcluster <- aov(areaLog ~ condition, data = nov[nov$category2 == "EsxA Cluster",]) %>%
  tukey_hsd()

stat.test.othercluster <- aov(areaLog ~ condition, data = nov[nov$category2 == "Other",]) %>%
  tukey_hsd()


summary(aov(areaLog ~ condition, data = nov[nov$category2 == "EsxA Cluster",]))
summary(aov(areaLog ~ condition, data = nov[nov$category2 == "Other",]))

stat.test.esxAcluster$cluster <- 2
stat.test.othercluster$cluster <- 1

stat.test <- bind_rows(stat.test.esxAcluster, stat.test.othercluster)

ggplot(nov, aes(x = condition, y = areaLog)) +
  geom_point((aes(color = cluster)),
             position = position_jitter(width = 0.3, seed = 1)) +
  theme_bw(base_size = 12) +
  geom_boxplot(aes(color = cluster, alpha = 0)) +
  labs(y = "N-terminal Acetylated Peptide\n Median Normalized Area\n(log10 transformed)",
       x = element_blank()) +
  theme(legend.position = "none",
        panel.grid = element_blank()) +
  facet_wrap("cluster", labeller = labeller(cluster = labs)) +
  scale_x_discrete(labels = c("WT",
                              expression(paste("\u0394",italic("emp1"))),
                              expression(paste("\u0394",italic("emp1"), "/comp"))))+
  geom_text_repel(aes(label = label),box.padding = 1.8,
                  position = position_jitter(width = 0.3,seed = 1)) +
  coord_cartesian(ylim = c(1.2,6.2)) +
  stat_pvalue_manual(
    stat.test, label = "p.adj.signif", 
    y.position = c(5.5, 5.75, 6),
    tip.length = 0.005) +
  ylab(expression(paste("N-terminal Acetylated Peptide Area (", log[10], " transformed)", sep = "")))

ggsave("MRA/tukey_analysis_of_EXSA_cluster.png")

