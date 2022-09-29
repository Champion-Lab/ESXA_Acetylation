####Whole Proteome Analysis

#prepare a whole proteome class for marinum and a dag logo
#class of the 26 substrates
WholeMMAR_proteome <- prepareProteome(fasta = 'C:\\Users\\colla\\Documents\\Research Summer 2022\\R scripts and Excel Files\\MMARproteome.fasta', 
                                      species = 'Mycobacterium marinum')
#load in the List of Substrates
dat <- unlist(read.delim('C:\\Users\\colla\\Documents\\Research Summer 2022\\SubstrateFasta.csv',
                         header = FALSE, as.is = TRUE))
seq_whole <- formatSequence(seq = dat, proteome = WholeMMAR_proteome, upstreamOffset =8,
                            downstreamOffset = 8)
#Build Background Model of substrates vs. Whole Proteome
whole_BG <- bg<-buildBackgroundModel(seq_whole, proteome=WholeMMAR_proteome,
                                     numSubsamples=10, testType = "fisher")


#custom grouping scheme for base Logo
color = c(Threonine = "#FF0066", Serine = "#FF0066", 
          Cysteine = "#FF0066", Asparagine = "#FF0066",
          Glutamine = "#FF0066", Tyrosine = "#FF0066",
          Lysine = "#107F80", Arginine = "#107F80",
          Histidine = "#107F80", Glycine = '#40007F', 
          Alanine = '#40007F', Valine = '#40007F',
          Leucine = '#40007F', Methionine = '#40007F',
          Isoleucine = '#40007F', Phenylalanine = '#40007F',
          Tryptophan = '#40007F', Proline = '#40007F',
          Aspartate = "#AA66FF", Glutamate = "#AA66FF")
symbol = c(Threonine = "T", Serine = "S", Cysteine = "C", 
           Asparagine = "N", Glutamine = "Q", Tyrosine = "Y",
           Lysine = "K", Arginine = "R", Histidine = "H",
           Glycine = "G", Alanine = "A", Valine = "V", 
           Leucine = "L", Methionine = "M", Isoleucine = "I",
           Phenylalanine = "F", Tryptophan = "W", Proline = "P", 
           Aspartate = "D", Glutamate = "E")
group = list(
  Threonine = c("T"), Serine = c("S"), Cysteine = c("C"), 
  Asparagine = c("N"), Glutamine = c("Q"), Tyrosine = c("Y"),
  Lysine = c("K"), Arginine = c("R"), Histidine = c("H"),
  Glycine = c("G"), Alanine = c("A"), Valine = c("V"), 
  Leucine = c("L"), Methionine = c("M"), Isoleucine = c("I"),
  Phenylalanine = c("F"), Tryptophan = c("W"), Proline = c("P"), 
  Aspartate = c("D"), Glutamate = c("E"))
addScheme(color = color, symbol = symbol, group = group)



#Statistical Test of Background Model
t0.1 <- testDAU(seq_whole, whole_BG)
t0.2 <- testDAU(seq_whole, whole_BG, groupingScheme = "custom_group")
#Heat map of the basic set up...
dagHeatmap(testDAUresults = t0.1, type = "diff")
dagHeatmap(testDAUresults = t1.1, type = "diff")
#Statistical Test for Base Logo
t1.1 <- testDAU(dagPeptides = seq_whole, dagBackground = whole_BG,
                groupingScheme = "hydrophobicity_KD")
t2.1 <- testDAU(dagPeptides = seq_whole, dagBackground = whole_BG,
                groupingScheme = "charge_group")
t3.1 <- testDAU(dagPeptides = seq_whole, dagBackground = whole_BG,
                groupingScheme = "chemistry_property_Mahler")
t4.1 <- testDAU(dagPeptides = seq_whole, dagBackground = whole_BG,
                groupingScheme = "hydrophobicity_KD_group")
#Base Logo of hydrophobicity
dagLogo(t0.2, legend = FALSE, labels = c('2', '3', '4', '5', '6', '7', '8', '9', '10', '11'))

dagLogo(t1.1, groupingSymbol = getGroupingSymbol(t1@group), 
        legend = TRUE, labels = c('2', '3', '4', '5', '6', '7', '8', '9', '10', '11'))
dagLogo(t2.1, groupingSymbol = getGroupingSymbol(t2@group),legend = TRUE, 
        labels = c('2', '3', '4', '5', '6', '7', '8', '9', '10', '11'))
dagLogo(t3.1, groupingSymbol = getGroupingSymbol(t3@group), legend = TRUE, 
        labels = c('2', '3', '4', '5', '6', '7', '8', '9', '10', '11'))
dagLogo(t4.1, groupingSymbol = getGroupingSymbol(t4@group), legend = TRUE,
        labels = c('2', '3', '4', '5', '6', '7', '8', '9', '10', '11'))




