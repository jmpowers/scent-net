rm(list=ls())

install.packages("devtools")
install.packages("vegan")
install.packages("vegan3d")
install.packages("MASS")
install.packages("tidyr")
install.packages("reshape2")
install.packages("plot3D")
install.packages("beepr")
install.packages("RColorBrewer")

library(devtools)
library(vegan)
library(vegan3d)
library(MASS)
library(tidyr)
library(reshape2)
library(plot3D)
library(RColorBrewer)
library(beepr)
library(dplyr)

beepr::beep(6)

read.shimadzu = function(filepath) {
  con = file(filepath, "r")
  mcpeaktable = list()
  simsearch = list()
  filelist = list()
  f = 0
  alreadyonheader = FALSE
  while ( TRUE ) {
    if(!alreadyonheader){
      line = readLines(con, n = 1)
      if ( length(line) == 0 ) { break } 
    } 
    if(line=="[Header]") {
      alreadyonheader=FALSE
      f = f + 1
      line = readLines(con, n = 1)
      filelist[[f]] <- basename(gsub("\\\\", "/", strsplit(line, "\t")[[1]][2]) )
    } 
    if(line=="[MC Peak Table]") {
      line = readLines(con, n = 1) # # of Peaks
        npeaks <- as.integer(strsplit(line[[1]], "\t", fixed=T)[[1]][[2]])
        line = readLines(con, n = 1) #Mass	TIC
        line = readLines(con, n = npeaks+1)#include the table header
        mcpeaktable[[f]] = read.delim(text=line)
        line = readLines(con, n = 1) #empty line after table
    }
    if(line == "[MS Similarity Search Results for Spectrum Process Table]") { 
      line = readLines(con, n=1) # # of spectra
      line = readLines(con, n=1) # table header
      l = 1
      simsearchlist = list()
      while (length(line) != 0 && line != "[Header]" ){
        simsearchlist[[l]] <- line
        l = l+1 
        line <- readLines(con, n=1)
      }
      simsearchstring <- paste(unlist(simsearchlist), collapse = "\n")
      simsearch[[f]] = read.delim(text=simsearchstring)
      if(length(line) != 0 && line=="[Header]") alreadyonheader= TRUE
      if ( length(line) == 0 ) { break }
    }
  }
  close(con)
  names(mcpeaktable) <- unlist(filelist)
  return(setNames(list(filelist, mcpeaktable, simsearch), c("filelist", "mcpeaktable", "simsearch")))
  
  } 
beep(4)

setwd("C:/Users/lucas/OneDrive - Christopher Newport University/Desktop/RMBL 2019/Data/R Volatiles/RVolatiles/Paul's Meadow")
myfiles <- read.shimadzu("Paul'sMeadowASCIIData.txt")

write.csv(myfiles$filelist, file = "filelist.csv")

myfiles$simsearch1 <- lapply(myfiles$simsearch, FUN = function(tab) tab[tab$Hit.. == 1,])
myfiles$godtable <- mapply(FUN = left_join, MoreArgs = list(by=c("Peak."="Spectrum.")), myfiles$mcpeaktable, myfiles$simsearch1, SIMPLIFY=FALSE)
goddesstable <- bind_rows(myfiles$godtable, .id = "Sample")

library(reshape2)
vol.all <- dcast(goddesstable, Sample~Name.x, sum, value.var="Area")
vol.all[,2] <- NULL
rownames(vol.all) <- vol.all[,1]
vol.all[,1] <- NULL
sapply(row.names(vol.all),FUN = function(old){paste(strsplit(old, "_")[[1]][1:3], collapse = ".")})
row.names(vol.all) <- sapply(row.names(vol.all),FUN = function(old){paste(strsplit(old, "_")[[1]][1:3], collapse = ".")})
row.names(vol.all)
type.vol <- factor(substr(row.names(vol.all),1,2))
type.vol

nmds.vol.all <- metaMDS(sqrt(vol.all), dist="bray", autotransform = FALSE)
ordiplot(nmds.vol.all, display =c("sites","species"), type = "p")
ordisample <- ordiplot(nmds.vol.all, display =c("sites","species"), type = "p")
species <- as.factor(substr(row.names(vol.all),1,2))
species
text(nmds.vol.all, display = "sites", labels = species, col = "black")

Compounds <- (data.frame((goddesstable [,c("Sample", "Ret.Time", "Area", "Name.x", "SI")])))
#better sorting
#Compounds$Date <- sapply(goddesstable[,"Sample"] ,FUN = function(dooty)(strsplit(dooty, "_")[[1]][3]))
#sapply(Compounds[,"Sample"] ,FUN = function(simple){paste(strsplit(simple, "_")[[1]][1:2], collapse = ".")})
#Compounds$Sample <- sapply(Compounds[,"Sample"] ,FUN = function(simple){paste(strsplit(simple, "_")[[1]][1:2], collapse = ".")})



CompoundsSorted <- subset(Compounds, Compounds[,"Name.x"] != "")
#type <- factor(substr(CompoundsSorted$Sample,1,2))
#type
#COMPOUNDS <- subset(CompoundsSorted, type != "Bl")
#str(CompoundsSorted)
#(colSums(CompoundsSorted["Area"]))

#filter for air and blank area from all other areas
#I want to sort by date and if the areas of air are greater than all else, I want the column to be null
#library(dplyr)
#group_by(CompoundsSorted,"Date", add = FALSE)
#ungroup(CompoundsSorted)
#CompoundsSorted$Date

#str(CompoundsSorted)
#setkey(CompoundsSorted, "Date")
#vol.filt <- !colMeans(CompoundsSorted[type == "Ai",])>10000






##############
####BOQUET####
#############
library(devtools)
install_github("jmpowers/bouquet")
library(bouquet)

library(magrittr)
# pipe operator %>% to chain functions
#contains sample GCMS_output and GCMS_metadata
GCMS_metadata<-read.csv ("Paul'sMeadowData2019.csv",header= TRUE,sep = ",", fileEncoding = "UTF-8-BOM")
str(GCMS_metadata)

GCMS_output<- CompoundsSorted
metadata <- load_metadata (GCMS_metadata, "SampleDate", "FileName", group = "Species", "Type", amount = "Biomass")
longdata <- load_longdata (GCMS_output, "Sample", "Ret.Time", "Name.x", "Area","SI", maxmatch=100)
sampletable <- make_sampletable (longdata)
str(longdata)

write.csv(row.names(sampletable), file = "filelist.csv")

row.names(sampletable) == metadata$sample


chemtable <-
  make_chemtable (longdata, metadata) %>%
  filter_RT(4, 20) %>%
  filter_match (0.80) %>%
  filter_freq (0.2, group = TRUE) %>%
  filter_contaminant (cont.list = c("Caprolactam")) %>%
  filter_area (min_maximum = 400000) %>%
  filter_ambient_ttest (sampletable, metadata, alpha = 0.0001)%>%
  combine_filters ()

dim(chemtable)

chemtable$filter_final = with(chemtable, freq.blank < 0.4 & filter_RT =="OK" & filter_match=="OK" & filter_contaminant=="OK" &filter_area=="OK" & filter_ambient_ttest=="OK" & filters_passed >= 6)
chemtable$filter_final
table(chemtable$filter_final)
finaltable <- prune_sampletable (sampletable, chemtable, metadata) 
names(finaltable)

dim(sampletable)
dim(metadata)

library(vegan)


vol <- finaltable
md <- metadata[metadata$sample %in% rownames(vol),]
vol <- vol / md$amount[row(vol)]

library(ggplot2)
ggplot(md, aes(y=amount,x=species)) + geom_point()  + scale_y_log10()

library(vegan3d)
library(rgl)

nmds.vol <- metaMDS(sqrt(vol), dist="bray", autotransform = FALSE)
ordiplot(nmds.vol, display =c("sites","species"))
text(nmds.vol, display = "species", cex=0.7, col="grey", labels=ifelse(colSums(vol)>1.5e8, colnames(vol),""))
species <- as.factor(substr(row.names(vol),1,2))
ordispider(nmds.vol, species, col= rainbow(nlevels(species)))
text(nmds.vol, labels=species, col = rainbow(nlevels(species))[as.integer(species)])
ordiellipse(nmds.vol, species, kind = "sd", conf = 0.80, draw = "polygon", col =(rainbow(nlevels(species))))

dist.vol <- vegdist(vol, method = "euclidean")
clust.vol <- hclust(dist.vol)
ordicluster(nmds.vol, clust.vol, col = "black")

nmds.vol.3d <- metaMDS(sqrt(vol), dist="bray", k = 3)
ordirgl(nmds.vol.3d, display = "species", col = "magenta", ax.col = "black", arr.col = "white", radius = 0.01)
orglspider(nmds.vol.3d, group = species, display = "sites", col = rainbow(nlevels(species)))
orglellipse(nmds.vol.3d, groups = species, display = "sites", kind ="sd", conf = 0.6, col = rainbow(nlevels(species)))
orgltext(nmds.vol.3d, text = species, cex = 1, col = rainbow(nlevels(species))[as.integer(species)])

writeASY(title = "ScentSpaceMap", outtype = "pdf")
writeASY(scene = ,
         title = "ScentSpaceMap", 
         outtype = "pdf")

ScentSpaceMap
#cluster
dist.vol <- vegdist(vol, method = "euclidean")
clust.vol <- hclust(dist.vol)
orglcluster(nmds.vol.3d, clust.vol, prune = 0, display = "sites",
            col = "black")

write.csv(finaltable, file = "FinalTableBigBoysDon'tCry.csv")
write.csv(chemtable, file = "ChemTableBigBoysDon'tCry.csv")

