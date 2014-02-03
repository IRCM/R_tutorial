library(cummeRbund)

setwd("/stockage/Bio-Analyses/LEC/LEC-BIF-P6/Analysis/cuffdiff/reference_only/cyto_vs_total")

cuff <- readCufflinks()

cuff

disp <- dispersionPlot(genes(cuff))
disp

genes.scv<-fpkmSCVPlot(genes(cuff))
genes.scv

isoforms.scv<-fpkmSCVPlot(isoforms(cuff))#sthash.i28kJowP.dpuf
isoforms.scv

s<-csScatterMatrix(genes(cuff))
s

v<-csVolcanoMatrix(genes(cuff))
v

dend<-csDendro(genes(cuff))
d

b<-csBoxplot(genes(cuff))
b

myGeneIds <- c("FBgn0000100", "FBgn0002466", "FBgn0002590", "FBgn0002593", "FBgn0003274", "FBgn0003517", "FBgn0004404", "FBgn0004828", "FBgn0004867", "FBgn0010173", "FBgn0010198", "FBgn0010303", "FBgn0010741", "FBgn0013756", "FBgn0014028", "FBgn0015393", "FBgn0016119", "FBgn0016917", "FBgn0020497", "FBgn0020545", "FBgn0023211", "FBgn0024238", "FBgn0024558", "FBgn0026566", "FBgn0026872", "FBgn0026879", "FBgn0027561", "FBgn0029937")
 
 myGenes<-getGenes(cuff,myGeneIds)
 
h<-csHeatmap(myGenes,cluster='both')
h

s<-csScatter(myGenes,"Cyto","Total",smooth=T) 
s

myGeneId <- "FBgn0000100"
myGene<-getGene(cuff,myGeneId)
myGene

gl<-expressionPlot(myGene)
gl