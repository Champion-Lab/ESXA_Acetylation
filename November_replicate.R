#November replicate
#loop through dfs and create normalized area columns
for (samp in sample_names_Nov) { #loops through all the sample names in Nov rep
  #create a string to create column 
  new_col <- paste0(samp, "_norm_area")
  #create a string to reference column
  old_col <- paste0("Area.", samp) 
  #pull rpoA area
  rpoA_area <- rpoA_areas_df[rpoA_areas_df$sample == samp,2]
  #set the value of the new column to the old column / rpoA area
  peptides_Nov_db[[new_col]] <- peptides_Nov_db[[old_col]] / rpoA_area
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

#keep only peptides that contain an nterm acetylated version of the peptide
final_Nov_acetylated <- final_Nov[final_Nov$peptides_acet_n > 0 & final_Nov$peptides_non_A_n > 0 ,]



#create a combined dataframe with relevant columns
Nov_WT <- data.frame(protein = final_Nov_acetylated$accession,
                     percent_acetylated = final_Nov_acetylated$percent_acet_WT,
                     percent_acetylated_sd = final_Nov_acetylated$percent_acet_WT_sd,
                     area_acet = final_Nov_acetylated$WT_av_acet,
                     sd_acet = final_Nov_acetylated$WT_sd_acet,
                     area_non_A = final_Nov_acetylated$WT_av_non_A,
                     sd_non_A = final_Nov_acetylated$WT_sd_non_A,
                     Condition = "WT")
Nov_Del <- data.frame(protein = final_Nov_acetylated$accession,
                      percent_acetylated = final_Nov_acetylated$percent_acet_Del,
                      percent_acetylated_sd = final_Nov_acetylated$percent_acet_Del_sd, 
                      area_acet = final_Nov_acetylated$Del_av_acet,
                      sd_acet = final_Nov_acetylated$Del_sd_acet,
                      area_non_A = final_Nov_acetylated$Del_av_non_A,
                      sd_non_A = final_Nov_acetylated$Del_sd_non_A,
                      Condition = "Del")
Nov_Comp <- data.frame(protein = final_Nov_acetylated$accession,
                       percent_acetylated = final_Nov_acetylated$percent_acet_Comp,
                       percent_acetylated_sd = final_Nov_acetylated$percent_acet_Comp_sd, 
                       area_acet = final_Nov_acetylated$Comp_av_acet,
                       sd_acet = final_Nov_acetylated$Comp_sd_acet,
                       area_non_A = final_Nov_acetylated$Comp_av_non_A,
                       sd_non_A = final_Nov_acetylated$Comp_sd_non_A,
                       Condition = "Comp")
plotdf_Nov <- bind_rows(Nov_WT, Nov_Del, Nov_Comp)

#plot % acetylated by sample type for each protein
Novplot <- ggplot(plotdf_Nov, aes(x = protein, ymin = percent_acetylated - percent_acetylated_sd,
                                  ymax = percent_acetylated + percent_acetylated_sd, fill = Condition)) +
  geom_bar(aes(y = percent_acetylated),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(position = position_dodge()) +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8)) +
  labs(y = "% acetylated", title = "November Replicate")
Novplot

#create dataframe with different line for each peptide acet version
plotdflong_Nov <- pivot_longer(plotdf_Nov, c(area_acet, area_non_A), names_to = "acetylated", values_to = "area")
plotdflong_Nov$sd <- NA
for (i in 1:nrow(plotdflong_Nov)) {
  if (plotdflong_Nov$acetylated[i] == "area_acet") {
    plotdflong_Nov$sd[i] <- plotdflong_Nov$sd_acet[i]
  } else {
    plotdflong_Nov$sd[i] <- plotdflong_Nov$sd_non_A[i]
  }
}

# create list of unique accession IDs, create a protein plot for each one
unique_accessions <- plotdflong_Nov$protein
Nov_protein_plots <- list()
for (accession in unique_accessions) {
  accession_df <- plotdflong_Nov[plotdflong_Nov$protein == accession,]
  Nov_protein_plots[[accession]] <- ggplot(accession_df,
                                                               aes(x = Condition,
                                                                   ymin = area - sd,
                                                                   ymax = area + sd,
                                                                   fill = acetylated)) +
    geom_bar(aes(y = area), stat = "identity", position = "dodge") +
    geom_errorbar(position = position_dodge(width = 0.9), width = 0.5) +
    theme_bw(base_size = 20) +
    labs(y = "Normalized Area", title = accession) +
    theme(panel.grid = element_blank()) +
    scale_fill_manual(labels = c("acetylated", "non-acetylated"), values = c("dodgerblue2", "green4"))
}



#save each protein plot as a png
for (protein in unique_accessions) {
  png(paste0("exports_2022_05_23/Individual_protein_plots/November/",protein, "_November_plot.png"), width = 1200, height = 700)
  show(Nov_protein_plots[[protein]])
  dev.off()
}


