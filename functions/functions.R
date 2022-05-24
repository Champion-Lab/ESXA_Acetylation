return_start_position <- function(peptide, accession, database){
  protein_names <- names(database)
  protein_index <- which(grepl(accession, protein_names))
  if(length(protein_index) > 1) {
    return(NA)
  } else if (length(protein_index) < 1) {
    return(NA)
  } else {
    protein <- db[[protein_index[1]]]
  }
  matches_df <- str_locate(protein, peptide)
  if (length(matches_df) < 1) {
    return(NA)
  } else {
    start_position <- matches_df[[1,1]]
  }
  return(start_position)
}

error_prop_addition <- function(sd_vector) {
  sd_vector <- sd_vector[!is.na(sd_vector)]
  values_n <- length(sd_vector)
  squared_sds <- rep(0, values_n)
  for (i in 1:values_n) {
    squared_sds[i] <- sd_vector[i]^2
  }
  sum_squares <- sum(squared_sds)
  return(sqrt(sum_squares))
}

remove_nacetyl_mod <- function(peptide) {
  return_string <- str_remove(peptide, "\\(\\+42\\.01\\)")
  return(return_string)
}

group_peptides <- function(combined_df, month = "Nov") {
  unique_proteins <- unique(combined_df$Accession_only.non_A)
  unique_proteins<-unique_proteins[!is.na(unique_proteins)]
  unique_proteins_n <- length(unique_proteins)
  colnames(combined_df) <- str_replace_all(string = names(combined_df),
                                         pattern = month,
                                         replacement = "")
  combined_df$Accession_only.non_A[is.na(combined_df$Accession_only.non_A)] <- "No Value"
  combined_df$Accession_only.acet[is.na(combined_df$Accession_only.acet)] <- "No Value"
  return_df <- data.frame(matrix(ncol = 24, nrow = unique_proteins_n))
  colnames(return_df) <- c("accession",
                           "peptides_acet",
                           "peptides_acet_n",
                           "peptides_non_A",
                           "peptides_non_A_n",
                           "total_peptide_variants",
                           "WT_av_acet",
                           "WT_sd_acet",
                           "Del_av_acet",
                           "Del_sd_acet",
                           "Comp_av_acet",
                           "Comp_sd_acet",
                           "WT_av_non_A",
                           "WT_sd_non_A",
                           "Del_av_non_A",
                           "Del_sd_non_A",
                           "Comp_av_non_A",
                           "Comp_sd_non_A",
                           "total_area_WT",
                           "total_area_WT_sd",
                           "total_area_Del",
                           "total_area_Del_sd",
                           "total_area_Comp",
                           "total_area_Comp_sd")
  
  for (i in 1:unique_proteins_n) {
    accession <- unique_proteins[i]
    temp_df <- combined_df[(combined_df$Accession_only.non_A == accession) | (combined_df$Accession_only.acet == accession),]
    peptides_acet <- paste0(temp_df$Peptide.acet[!is.na(temp_df$Peptide.acet)], collapse = ";")
    peptides_acet_n <- sum(!is.na(temp_df$Peptide.acet))
    peptides_non_A <- paste0(temp_df$Peptide.non_A[!is.na(temp_df$Peptide.non_A)], collapse = ";")
    peptides_non_A_n <- sum(!is.na(temp_df$Peptide.non_A))
    total_peptide_variants <- nrow(temp_df)
    WT_av_acet <- sum(temp_df$WT_avg.acet, na.rm = T)
    WT_sd_acet <- error_prop_addition(temp_df$WT_stdev.acet)
    Del_av_acet <- sum(temp_df$Del_avg.acet, na.rm = T)
    Del_sd_acet <- error_prop_addition(temp_df$Del_stdev.acet)
    Comp_av_acet <- sum(temp_df$Comp_avg.acet, na.rm = T)
    Comp_sd_acet <- error_prop_addition(temp_df$Comp_stdev.acet)
    WT_av_non_A <- sum(temp_df$WT_avg.non_A, na.rm = T)
    WT_sd_non_A <- error_prop_addition(temp_df$WT_stdev.non_A)
    Del_av_non_A <- sum(temp_df$Del_avg.non_A, na.rm = T)
    Del_sd_non_A <- error_prop_addition(temp_df$Del_stdev.non_A)
    Comp_av_non_A <- sum(temp_df$Comp_avg.non_A, na.rm = T)
    Comp_sd_non_A <- error_prop_addition(temp_df$Comp_stdev.non_A)
    total_area_WT <- sum(temp_df$total_area_WT_, na.rm = T)
    total_area_WT_sd <- error_prop_addition(temp_df$total_area_WT_sd_)
    total_area_Del <- sum(temp_df$total_area_Del_, na.rm = T)
    total_area_Del_sd <- error_prop_addition(temp_df$total_area_Del_sd_)
    total_area_Comp <- sum(temp_df$total_area_Comp_, na.rm = T)
    total_area_Comp_sd <- error_prop_addition(temp_df$total_area_Comp_sd_)
    
    protein_line <- data.frame(accession = accession,
                               peptides_acet = peptides_acet,
                               peptides_acet_n = peptides_acet_n,
                               peptides_non_A = peptides_non_A,
                               peptides_non_A_n = peptides_non_A_n,
                               total_peptide_variants = total_peptide_variants,
                               WT_av_acet = WT_av_acet,
                               WT_sd_acet = WT_sd_acet,
                               Del_av_acet = Del_av_acet,
                               Del_sd_acet = Del_sd_acet,
                               Comp_av_acet = Comp_av_acet,
                               Comp_sd_acet = Comp_sd_acet,
                               WT_av_non_A = WT_av_non_A,
                               WT_sd_non_A = WT_sd_non_A,
                               Del_av_non_A = Del_av_non_A,
                               Del_sd_non_A = Del_sd_non_A,
                               Comp_av_non_A = Comp_av_non_A,
                               Comp_sd_non_A = Comp_sd_non_A,
                               total_area_WT = total_area_WT,
                               total_area_WT_sd = total_area_WT_sd,
                               total_area_Del = total_area_Del,
                               total_area_Del_sd = total_area_Del_sd,
                               total_area_Comp = total_area_Comp,
                               total_area_Comp_sd = total_area_Comp_sd)
    return_df[i,] <- protein_line
  }
  
  return(return_df)
}

