###########Making Nice Figure############
#Read in N-terminal Acetylome from the path on my device
Acetylome<-read.csv('C:\\Users\\colla\\OneDrive\\Documents\\Research Summer 2022\\N-terminal Acetylome.csv')
AcetylomeFasta<-AcetylomeFasta

####Creating Dag object classess needed for analysis
dat <- unlist(read.delim('C:\\Users\\colla\\OneDrive\\Documents\\Research Summer 2022\\SubstrateFasta.csv',
                         header = FALSE, as.is = TRUE))

#prepare an object of Proteome Class from a fasta file
proteome <- prepareProteome(fasta = 'C:\\Users\\colla\\OneDrive\\Documents\\Research Summer 2022\\R scripts and Excel Files\\AcetylomeFasta.csv',
                            species = "M.marinum")

#prepare an object of dagPeptides Class
seq <- formatSequence(seq = dat, proteome = proteome, upstreamOffset = 1,
                      downstreamOffset = 1)
#build acetylome background model
AcetylomeModel<-buildBackgroundModel(
  dagPeptides=seq,
  background = "inputSet",
  model = "any",
  targetPosition = "Nterminus",
  uniqueSeq = FALSE,
  numSubsamples = 300L,
  rand.seed = 1,
  replacement = FALSE,
  testType = "fisher",
  proteome
)


###Background scripts as described in paper
bg<-buildBackgroundModel(seq, proteome=proteome,
                         numSubsamples=10)
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
#Heat Map of Results
t0 <- testDAU(seq, bg, groupingScheme = "custom_group")
dagHeatmap(testDAUresults = t0, type = "diff")
dagLogo(t0, legend = TRUE, labels = c('2', '3', '4', '5', '6', '7', '8', '9', '10', '11'),
        fontsize = 12, )
              
