#Putative MMAR_1839 Target Analysis
#Generating Sequence Logos of 23 putative protein targets of MMAR_1839

#libraries
library(ggplot2)
library(ggseqlogo)
#Read in csv of putative target amino acid sequences
Substrates<-read.csv('C:\\Users\\colla\\Downloads\\1839_substrates_10aa.csv')

#Add color blind scheme from Prism, attributed to amino acid classification
cs1 = make_col_scheme(chars=c('G','P','A', 'V', 'L', 'I','M','C',
                              'F','Y','W','H','K','R','Q','N','E',
                              'D','S','T'), 
                      groups=c('Hydrophobic', 'Hydrophobic', 'Hydrophobic',
                               'Hydrophobic', 'Hydrophobic', 'Hydrophobic', 
                               'Hydrophobic', 'Polar', 'Hydrophobic', 'Polar',
                               'Hydrophobic', 'Basic', 'Basic', 'Basic', 
                               'Polar', 'Polar', 'Acidic', 'Acidic', 'Polar', 
                               'Polar'), 
                      cols=c('#40007F', '#40007F','#40007F','#40007F','#40007F',
                             '#40007F','#40007F','#FF0066', '#40007F', '#FF0066',
                             '#40007F', '#107F80', '#107F80', '#107F80','#FF0066',
                             '#FF0066', '#AA66FF', '#AA66FF','#FF0066','#FF0066'))
#Build Sequence Logo graph with R package ggseqlogog
ggseqlogo(Substrates$First.10.AA, col_scheme=cs1,)+
  scale_x_continuous(breaks=1:10, labels=2:11)+
  xlab("amino acid position")
                                            

                                                               