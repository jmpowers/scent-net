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
library(bouquet)
library(magrittr)
library(scales)

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
myfiles <- read.shimadzu("DryMeadowTotal.txt")


myfiles$simsearch1 <- lapply(myfiles$simsearch, FUN = function(tab) tab[tab$Hit.. == 1,])
myfiles$godtable <- mapply(FUN = left_join, MoreArgs = list(by=c("Peak."="Spectrum.")), myfiles$mcpeaktable, myfiles$simsearch1, SIMPLIFY=FALSE)
goddesstable <- bind_rows(myfiles$godtable, .id = "Sample")

library(reshape2)
vol.all <- dcast(goddesstable, Sample~Name.x, sum, value.var="Area")
vol.all[,2] <- NULL
rownames(vol.all) <- vol.all[,1]
vol.all[,1] <- NULL
row.names(vol.all)
type.vol <- factor(substr(row.names(vol.all),1,2))
type.vol

#nmds.vol.all <- metaMDS(sqrt(vol.all), dist="bray", autotransform = FALSE)
#ordiplot(nmds.vol.all, display =c("sites","species"), type = "p")
#ordisample <- ordiplot(nmds.vol.all, display =c("sites","species"), type = "p")
#species <- as.factor(substr(row.names(vol.all),1,2))
#species
#text(nmds.vol.all, display = "sites", labels = species, col = "black")

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
####BOUQUET####
#############

install_github("jmpowers/bouquet")

# pipe operator %>% to chain functions
#contains sample GCMS_output and GCMS_metadata
GCMS_metadata<-read.csv ("BIGMETADATA.csv",header= TRUE,sep = ",", fileEncoding = "UTF-8-BOM")

GCMS_metadata$Year <- as.factor(GCMS_metadata$Year)
str(GCMS_metadata)

GCMS_output<- CompoundsSorted
metadata <- load_metadata (GCMS_metadata, "SampleDate", "FileName", group = "Specimen", "Type", amount = "Biomass")
longdata <- load_longdata (GCMS_output, "Sample", "Ret.Time", "Name.x", "Area","SI", maxmatch=100)
sampletable <- make_sampletable (longdata)
str(longdata)
metadata$type

chemtable <-
  make_chemtable (longdata, metadata) %>%
  filter_RT(4, 20) %>%
  filter_match (0.8) %>%
  filter_freq (0.2, group = TRUE) %>%
  filter_contaminant (cont.list = c("Caprolactam")) %>%
  filter_ambient_ratio(sampletable, metadata, ratio = 3)%>%
  filter_area (min_maximum = 400000) %>%
  filter_ambient_ttest (sampletable, metadata, alpha = 0.0001)%>%
  combine_filters ()

myfilestokeep <- 
chemtable$filter_keep <- chemtable$name %in% mylistocfcompoundstokeep
SixChemtable
dim(chemtable)
dim(sampletable)
dim(metadata)

row.names(sampletable) == metadata$sample

chemtable$filter_final = with(chemtable, freq.blank < 0.4 & filter_RT =="OK" & filter_match=="OK" & filter_contaminant=="OK" &filter_area=="OK" & filter_ambient_ttest=="OK" & filters_passed >= 7)
chemtable$filter_final
table(chemtable$filter_final)
finaltable <- prune_sampletable (sampletable, chemtable, metadata) 
names(finaltable)


vol <- finaltable
md <- metadata[metadata$sample %in% rownames(vol),]
vol <- vol / md$amount[row(vol)]

#PROPORTION#
PresentFile<- read.csv("PresentationFinalTableThing.csv", header = TRUE, sep = ",", fileEncoding = "UTF-8-BOM")
row.names(PresentFile) <- PresentFile$Sample
PresentFile$Sample = NULL
row.names(PresentFile)
row.names(PresentMetadata) <- PresentMetadata$Sample
PresentMetadata$Sample = NULL
row.names(PresentMetadata)

sampletable6 <- sampletable[rownames(sampletable) %in% rownames(PresentFile),]
alphaPinene <- sampletable6$`(1R)-2,6,6-Trimethylbicyclo[3.1.1]hept-2-ene`
alphaPinene
PresentFile <- cbind(PresentFile, alphaPinene)

write.csv(PresentFile, file = "FilteredCompoundTable.csv")

PresentMetadata <- read.csv("PresentationMetaData.csv", header = TRUE, sep = ",", fileEncoding = "UTF-8-BOM")

dim(PresentFile)
names(PresentFile)
Presentvol <- decostand(PresentFile, "total")
write.csv(Presentvol, "ProportionFilteredCompoundTableProp.csv")

pals <- c("Greys","Greens","Purples","RdPu","Blues","Oranges")
levels(SpeciesYear)
spyr.pal <- unlist(lapply(lapply(pals, brewer.pal, n=9), function(x) x[c(5,7,9)]))
par(bg="white")
Presentnmds.vol <- metaMDS(Presentvol, dist="bray", autotransform = FALSE)
ordiplot(Presentnmds.vol, type = "n", xlim=c(-1.0,1.0), ylim=c(-1.9,1.5))
points(Presentnmds.vol, display = "sites", col = spyr.pal[as.integer(SpeciesYear)], bg = alpha(spyr.pal, alpha=0.5)[as.integer(SpeciesYear)], pch = 21)

colSums(Presentvol)
text(Presentnmds.vol, display = "species", cex=0.7, col="black", labels=ifelse(colSums(Presentvol)>2.5, colnames(Presentvol),""))

species <- PresentMetadata$Species
YearMan <- as.factor(PresentMetadata$Year)
SpeciesYear <- factor(paste0(species,substr(as.character(YearMan), 4,4)))

ordispider(Presentnmds.vol, SpeciesYear, col= spyr.pal, lwd=2)
text(Presentnmds.vol, labels=SpeciesYear, col = spyr.pal[as.integer(SpeciesYear)])
ordiellipse(Presentnmds.vol, SpeciesYear, kind = "sd", conf = 0.75, draw = "polygon", col =spyr.pal)
legend("topleft", legend=as.character(levels(SpeciesYear)), fill=alpha(spyr.pal, alpha=0.5))
title(main = "Comparison of Volatiles Across Species & Years")

#LogTotalPeakArea
PresentFile$Total <- rowSums(PresentFile)
PresentFile$LogTotal <- log(PresentFile$Total)
names(PresentMetadata)
DryMeadowLog <- PresentFile$LogTotal
DryMeadowSpecies <- PresentMetadata$Species
DryMeadowYear <- PresentMetadata$Year
DryMeadowSpeciesYear <- factor(paste0(DryMeadowSpecies,substr(as.character(DryMeadowYear), 4,4)))

drymeadowpals <- c("Greys","Greens","Purples","RdPu","Blues","Oranges")
drymeadow.pal <- unlist(lapply(lapply(drymeadowpals, brewer.pal, n=9), function(x) x[c(4, 6, 8)]))
install.packages("grDevices")
library(grDevices)
par(mar=c(5,5,5,2), mgp=c(3,1,0))
plot(DryMeadowSpeciesYear, DryMeadowLog, xlab = "Species Across Years", ylab = "Log(TotalPeakArea)", col = alpha(c("green3", "red3", "blue3"), alpha = 0.5), xaxt = "n")
plot(DryMeadowSpeciesYear, DryMeadowLog, xlab = "Species Across Years", ylab = "Log(TotalPeakArea)", col = drymeadow.pal, xaxt = "n")
axis(1, at=c(2,5,8,11,14,17), labels=(levels(DryMeadowSpecies)))
axis(3, at=1:18, labels=c("2017","2018", "2019", "2017","2018", "2019","2017","2018", "2019","2017","2018", "2019","2017","2018", "2019","2017","2018","2019" ), outer = FALSE, las = 2)
title( main = "Total Volatiles by Species and Year", line = +2)
legend("bottomright", legend=as.character(levels(YearMan)), fill=alpha(c("green3", "red3", "blue3"), alpha=0.5))

cor(Nectar, NectarLog)
IpoNectarCorLM <- summary(lm(LogTotal~X24Hr.Nectar.Production, data = NectarProduction))
IpoNectarCorLM


#Stats#
adonis(Presentvol ~ species * YearMan, permu = 999)

Present.cap <- capscale(Presentvol ~ species * YearMan)
anova(Present.cap, by="term")
plot3d(Present.cap)

#YEAR#
YearMan
ypals <- c("Reds", "Greens","Blues")
yr.pal <- unlist(lapply(lapply(ypals, brewer.pal, n=9), function(x) x[c(7)]))
ordiplot(Presentnmds.vol, type = "n", xlim=c(-1,1), ylim=c(-1.3,1))
points(Presentnmds.vol, display = "sites", col = yr.pal[as.integer(YearMan)], bg = alpha(yr.pal, alpha=0.5)[as.integer(YearMan)], pch = 21)
ordispider(Presentnmds.vol, groups = species, col= sp.pal, lwd=2)
text(Presentnmds.vol, labels=SpeciesYear, col = sp.pal [as.integer(species)])
ordiellipse(Presentnmds.vol, YearMan, kind = "sd", conf = 0.80, draw = "polygon", col = yr.pal)
legend("topleft", legend=as.character(levels(YearMan)), fill=alpha(yr.pal, alpha=0.5))
title(main = "Comparison of Volatiles Across Years")

#Species#
species
ordiplot(Presentnmds.vol, type = "n", xlim=c(-1,1), ylim=c(-1.4,1.2))
points(Presentnmds.vol, display = "sites", col = sp.pal[as.integer(species)], bg = alpha(sp.pal, alpha=0.5)[as.integer(species)], pch = 21)
spals <- c("Greys","Greens","Purples","RdPu","Blues","Oranges")
sp.pal <- unlist(lapply(lapply(spals, brewer.pal, n=9), function(x) x[c(7)]))
ordispider(Presentnmds.vol, groups = species, col= sp.pal, lwd=2)
text(Presentnmds.vol, labels=SpeciesYear, col = sp.pal [as.integer(species)])
ordiellipse(Presentnmds.vol, species, kind = "sd", conf = 0.80, draw = "polygon", col = sp.pal)
legend("topleft", legend=as.character(levels(species)), fill=alpha(sp.pal, alpha=0.5))
title(main = "Comparison of Volatiles Across Species")

Wetness <- read.csv("envfitwetness.csv", header = TRUE, sep = ",", fileEncoding = "UTF-8-BOM")
Wetness$Year <- as.factor(Wetness$Year)
envfit(Presentnmds.vol, Wetness, permu = 999)

#AC#
row.names(PresentFile)
AeCoMetaData <- subset(PresentMetadata, PresentMetadata$Species == "AC")
row.names(AeCoMetaData)
AeCoVolData <- PresentFile[rownames(PresentFile) %in% rownames(AeCoMetaData),]
ACVol<- decostand(AeCoVolData, "total")
ACnmds.vol <- metaMDS(ACVol, dist="bray", autotransform = FALSE)
ACSpecies <- as.factor(AeCoMetaData$Species)
ACYear <- as.factor(AeCoMetaData$Year)
ACSpeciesYear <- factor(paste0(ACSpecies,substr(as.character(ACYear), 4,4)))
ACpals <- "Greys"
ACpals
levels(ACYear)
ACYear.pal <- unlist(lapply(lapply(ACpals, brewer.pal, n=9), function(x) x[c(4,7,9)]))
ordiplot(ACnmds.vol, type="n", xlim=c(-1,1), ylim=c(-1.5,1.5))
points(ACnmds.vol, display = "sites", col = ACYear.pal[as.integer(ACYear)], bg = alpha(ACYear.pal, alpha=0.5)[as.integer(ACYear)], pch = 21)
ordispider(ACnmds.vol, ACYear, col= ACYear.pal, lwd=2)
text(ACnmds.vol, labels=ACSpeciesYear, col = ACYear.pal[as.integer(ACSpeciesYear)])
ordiellipse(ACnmds.vol, ACSpeciesYear, kind = "sd", conf = 0.75, draw = "polygon", col =ACYear.pal)
text(ACnmds.vol, labels=ACSpeciesYear, col = ACYear.pal[as.integer(ACSpeciesYear)])
legend("topleft", legend=as.character(levels(ACYear)), fill=alpha(ACYear.pal, alpha=0.5))
title(main = substitute(paste(italic("Arenaria congesta"))))
adonis(ACVol ~ ACYear)
envfit
#AM#
row.names(PresentFile)
AcMiMetaData <- subset(PresentMetadata, PresentMetadata$Species == "AM")
row.names(AcMiMetaData)
AcMiVolData <- PresentFile[rownames(PresentFile) %in% rownames(AcMiMetaData),]
AMVol<- decostand(AcMiVolData, "total")
AMnmds.vol <- metaMDS(AMVol, dist="bray", autotransform = FALSE)
AMSpecies <- as.factor(AcMiMetaData$Species)
AMYear <- as.factor(AcMiMetaData$Year)
AMSpeciesYear <- factor(paste0(AMSpecies,substr(as.character(AMYear), 4,4)))
AMpals <- "Greens"
levels(AMYear)
AMYear.pal <- unlist(lapply(lapply(AMpals, brewer.pal, n=9), function(x) x[c(4,7,9)]))
ordiplot(AMnmds.vol, type="n", xlim=c(-1,1), ylim=c(-1.5,1.5))
points(AMnmds.vol, display = "sites", col = AMYear.pal[as.integer(AMYear)], bg = alpha(AMYear.pal, alpha=0.5)[as.integer(AMYear)], pch = 21)
ordispider(AMnmds.vol, AMYear, col= AMYear.pal, lwd=2)
text(AMnmds.vol, labels=AMSpeciesYear, col = AMYear.pal[as.integer(AMSpeciesYear)])
ordiellipse(AMnmds.vol, AMSpeciesYear, kind = "sd", conf = 0.75, draw = "polygon", col =AMYear.pal)
text(AMnmds.vol, labels=AMSpeciesYear, col = AMYear.pal[as.integer(AMSpeciesYear)])
legend("topleft", legend=as.character(levels(AMYear)), fill=alpha(AMYear.pal, alpha=0.5))
title(main = substitute(paste(italic("Achillea millefolium"))))
adonis(AMVol ~ AMYear)

#CR#
row.names(PresentFile)
CaRoMetaData <- subset(PresentMetadata, PresentMetadata$Species == "CR")
row.names(CaRoMetaData)
CaRoVolData <- PresentFile[rownames(PresentFile) %in% rownames(CaRoMetaData),]
CRVol<- decostand(CaRoVolData, "total")
CRnmds.vol <- metaMDS(CRVol, dist="bray", autotransform = FALSE)
CRSpecies <- as.factor(CaRoMetaData$Species)
CRYear <- as.factor(CaRoMetaData$Year)
CRSpeciesYear <- factor(paste0(CRSpecies,substr(as.character(CRYear), 4,4)))
CRpals <- "Purples"
levels(CRYear)
CRYear.pal <- unlist(lapply(lapply(CRpals, brewer.pal, n=9), function(x) x[c(4,6,9)]))
ordiplot(CRnmds.vol, type="n", xlim=c(-1,1), ylim=c(-1,1))
points(CRnmds.vol, display = "sites", col = CRYear.pal[as.integer(CRYear)], bg = alpha(CRYear.pal, alpha=0.5)[as.integer(CRYear)], pch = 21)
ordispider(CRnmds.vol, CRYear, col= CRYear.pal, lwd=2)
text(CRnmds.vol, labels=CRSpeciesYear, col = CRYear.pal[as.integer(CRSpeciesYear)])
ordiellipse(CRnmds.vol, CRSpeciesYear, kind = "sd", conf = 0.75, draw = "polygon", col =CRYear.pal)
legend("topleft", legend=as.character(levels(CRYear)), fill=alpha(CRYear.pal, alpha=0.5))
title(main = substitute(paste(italic("Campanula rotundifolia"))))
adonis(CRVol ~ CRYear)

#GB#
row.names(PresentFile)
GaBoMetaData <- subset(PresentMetadata, PresentMetadata$Species == "GB")
row.names(GaBoMetaData)
GaBoVolData <- PresentFile[rownames(PresentFile) %in% rownames(GaBoMetaData),]
GBVol<- decostand(GaBoVolData, "total")
GBnmds.vol <- metaMDS(GBVol, dist="bray", autotransform = FALSE)
GBSpecies <- as.factor(GaBoMetaData$Species)
GBYear <- as.factor(GaBoMetaData$Year)
GBSpeciesYear <- factor(paste0(GBSpecies,substr(as.character(GBYear), 4,4)))
GBpals <- "Blues"
levels(GBYear)
GBYear.pal <- unlist(lapply(lapply(GBpals, brewer.pal, n=9), function(x) x[c(5,7,9)]))
ordiplot(GBnmds.vol, type="n", xlim=c(-1,1), ylim=c(-1.5,1.5))
points(GBnmds.vol, display = "sites", col = GBYear.pal[as.integer(GBYear)], bg = alpha(GBYear.pal, alpha=0.5)[as.integer(GBYear)], pch = 21)
ordispider(GBnmds.vol, GBYear, col= GBYear.pal, lwd=2)
text(GBnmds.vol, labels=GBSpeciesYear, col = GBYear.pal[as.integer(GBSpeciesYear)])
ordiellipse(GBnmds.vol, GBSpeciesYear, kind = "sd", conf = 0.75, draw = "polygon", col =GBYear.pal)
legend("topleft", legend=as.character(levels(GBYear)), fill=alpha(GBYear.pal, alpha=0.5))
title(main = substitute(paste(italic("Galium boreale"))))
adonis(GBVol ~ GBYear)

#ES#
row.names(PresentFile)
ErSuMetaData <- subset(PresentMetadata, PresentMetadata$Species == "ES")
row.names(ErSuMetaData)
ErSuVolData <- PresentFile[rownames(PresentFile) %in% rownames(ErSuMetaData),]
ESVol<- decostand(ErSuVolData, "total")
ESnmds.vol <- metaMDS(ESVol, dist="bray", autotransform = FALSE)
ESSpecies <- as.factor(ErSuMetaData$Species)
ESYear <- as.factor(ErSuMetaData$Year)
ESSpeciesYear <- factor(paste0(ESSpecies,substr(as.character(ESYear), 4,4)))
ESpals <- "RdPu"
levels(ESYear)
ESYear.pal <- unlist(lapply(lapply(ESpals, brewer.pal, n=9), function(x) x[c(5,7,9)]))
ordiplot(ESnmds.vol, type="n", xlim=c(-1,1), ylim=c(-1.5,1.5))
points(ESnmds.vol, display = "sites", col = ESYear.pal[as.integer(ESYear)], bg = alpha(ESYear.pal, alpha=0.5)[as.integer(ESYear)], pch = 21)
ordispider(ESnmds.vol, ESYear, col= ESYear.pal, lwd=2)
text(ESnmds.vol, labels=ESSpeciesYear, col = ESYear.pal[as.integer(ESSpeciesYear)])
ordiellipse(ESnmds.vol, ESSpeciesYear, kind = "sd", conf = 0.75, draw = "polygon", col =ESYear.pal)
legend("topleft", legend=as.character(levels(ESYear)), fill=alpha(ESYear.pal, alpha=0.5))
title(main = substitute(paste(italic("Eriogonum subalpinum"))))
adonis(ESVol ~ ESYear)

#HV#
row.names(PresentFile)
HeViMetaData <- subset(PresentMetadata, PresentMetadata$Species == "HV")
row.names(HeViMetaData)
HeViVolData <- PresentFile[rownames(PresentFile) %in% rownames(HeViMetaData),]
HVVol<- decostand(HeViVolData, "total")
HVnmds.vol <- metaMDS(HVVol, dist="bray", autotransform = FALSE)
HVSpecies <- as.factor(HeViMetaData$Species)
HVYear <- as.factor(HeViMetaData$Year)
HVSpeciesYear <- factor(paste0(HVSpecies,substr(as.character(HVYear), 4,4)))
HVpals <- "Oranges"
levels(HVYear)
HVYear.pal <- unlist(lapply(lapply(HVpals, brewer.pal, n=9), function(x) x[c(4,7,9)]))
ordiplot(HVnmds.vol, type="n", xlim=c(-1,1), ylim=c(-1.5,1.5))
points(HVnmds.vol, display = "sites", col = HVYear.pal[as.integer(HVYear)], bg = alpha(HVYear.pal, alpha=0.5)[as.integer(HVYear)], pch = 21)
ordispider(HVnmds.vol, HVYear, col= HVYear.pal, lwd=2)
text(HVnmds.vol, labels=HVSpeciesYear, col = HVYear.pal[as.integer(HVSpeciesYear)])
ordiellipse(HVnmds.vol, HVSpeciesYear, kind = "sd", conf = 0.75, draw = "polygon", col =HVYear.pal)
legend("topleft", legend=as.character(levels(HVYear)), fill=alpha(HVYear.pal, alpha=0.5))
title(main = substitute(paste(italic("Heterotheca villosa"))))

adonis(HVVol ~ HVYear)

#TotalPeakArea#
TotalPeakArea <- PresentMetadata
TotalPeakArea
TotalPeakArea$Total <- rowSums(PresentFile)
TotalPeakArea$Log <- log(TotalPeakArea$Total)
names(TotalPeakArea)
LogVol <- TotalPeakArea$Log
LogSpecies <- as.factor(TotalPeakArea$Species)
LogYear <- as.factor(TotalPeakArea$Year)
plot(x = LogSpecies, y = LogVol, xlab = "Species", ylab = "Log(TotalPeakArea)", main  = "Total Peak Area by Species")
hist(LogVol, breaks = LogSpecies)
install.packages("car")
library(car)
Anova(lm(LogVol~LogSpecies*LogYear,contrasts=list(LogSpecies=contr.sum, LogYear=contr.sum)),type=3)

#3D#
library("rgl")
Presentvol.3d <- metaMDS(sqrt(Presentvol), dist="bray", k = 3)
ordirgl(Presentvol.3d, type = "n", col = "magenta", ax.col = "black", arr.col = "white", radius = 0.01)
play3d(spin3d(axis = c(0, 1, 0), rpm = 5), duration = 100)
orglspider(Presentvol.3d, group = SpeciesYear, display = "sites", col = spyr.pal)
orgltext(Presentvol.3d, text = SpeciesYear, cex = 1, col = spyr.pal[as.integer(SpeciesYear)])
orglellipse(Presentvol.3d, groups = SpeciesYear, display = "sites", kind ="sd", conf = 0.7, col = spyr.pal)

orglspider(Presentvol.3d, group = species, col= rainbow(nlevels(species)), lwd=2)
orgltext(Presentvol.3d, text = species, cex = 1, col = rainbow(nlevels(species))[as.integer(species)])
orglellipse(Presentvol.3d, groups = species, display = "sites", kind ="sd", conf = 0.7, col = rainbow(nlevels(species)))

install.packages("purrr")
install.packages("magick")
install.packages("readr")
library(dplyr) # use for fixing up data
library(purrr) # for mapping over a function
library(readr)
library(magick)

Angle1 <- 5
gif.delay <- 3

plot3d(tree3D, aspect="iso", size=1, axes=F, col="grey21")
Angle <- rep(Angle1 * pi / 180, 360/Angle1)

dir.create("animation2")
Animation2.dir <- paste(getwd(), "/animation2/", sep="")

for (i in seq(Angle)) {
  view3d(userMatrix = rotate3d(par3d("userMatrix"),
                               Angle[i], 0, 1, 0))
  
  rgl.snapshot(filename=paste(paste(Animation2.dir, "frame-", sep=""),
                              sprintf("%03d", i), ".png", sep=""))
}

move.to.animation2.dir <- paste("cd", Animation2.dir, "&&")
ImageMagick.code <- paste("convert -delay ", gif.delay,
                          " -loop 0 frame*.png animated2.gif", sep="")
run.into.console <- paste(move.to.animation2.dir, ImageMagick.code)

system(run.into.console)

list.files(path = "/animation/", pattern = "*.png", full.names = T) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=2) %>% # animates, can opt for number of loops
  image_write("nmdsScentSpace.gif")

#Nectar#
IpomopsisTotal <-read.shimadzu("IpomopsisTotal.txt")
IpomopsisTotal$simsearch1 <- lapply(IpomopsisTotal$simsearch, FUN = function(tab) tab[tab$Hit.. == 1,])
IpomopsisTotal$godtable <- mapply(FUN = left_join, MoreArgs = list(by=c("Peak."="Spectrum.")), IpomopsisTotal$mcpeaktable, IpomopsisTotal$simsearch1, SIMPLIFY=FALSE)
IpoGoddesstable <- bind_rows(IpomopsisTotal$godtable, .id = "Sample")

library(reshape2)
IpoVol.all <- dcast(IpoGoddesstable, Sample~Name.x, sum, value.var="Area")
IpoVol.all[,2] <- NULL
rownames(IpoVol.all) <- IpoVol.all[,1]
IpoVol.all[,1] <- NULL
row.names(IpoVol.all)

IpoCompounds <- (data.frame((IpoGoddesstable [,c("Sample", "Ret.Time", "Area", "Name.x", "SI")])))
IpoCompoundsSorted <- subset(IpoCompounds, IpoCompounds[,"Name.x"] != "")
row.names(IpoVol.all)

IpoGCMS_metadata<-read.csv ("IpoMetaData.csv",header= TRUE,sep = ",", fileEncoding = "UTF-8-BOM")

IpoGCMS_output<- IpoCompoundsSorted
names(IpoGCMS_metadata)
Ipometadata <- load_metadata (IpoGCMS_metadata, "FileName","Date", group = "Site", "Type", amount = "Biomass")
Ipolongdata <- load_longdata (IpoGCMS_output, "Sample", "Ret.Time", "Name.x", "Area","SI", maxmatch=100)
Iposampletable <- make_sampletable (Ipolongdata)
str(Ipolongdata)
Ipometadata$type

Ipochemtable <- 
  make_chemtable (Ipolongdata, Ipometadata)%>%
  filter_RT(4, 20) %>%
  filter_match (0.8) %>%
  filter_freq (0.2, group = TRUE) %>%
  filter_contaminant (cont.list = c("Caprolactam")) %>%
  filter_ambient_ratio(Iposampletable, Ipometadata, ratio = 3)%>%
  filter_ambient_ttest (Iposampletable, Ipometadata, alpha = 0.0001)%>%
  combine_filters ()
names(Ipochemtable)
Ipochemtable$filter
Ipochemtable$filter_final = with(Ipochemtable, freq.blank < 0.4 & filter_RT =="OK" & filter_match=="OK" & filter_contaminant=="OK"  & filter_ambient_ttest=="OK" & filters_passed >= 5)
Ipochemtable$filter_final
table(Ipochemtable$filter_final)

names(Ipochemtable)
Ipofinaltable <- prune_sampletable (Iposampletable, Ipochemtable, Ipometadata) 
names(Ipofinaltable)

write.csv(Ipofinaltable, file = "IpoFinalTable.csv")
write.csv(Iposampletable, file = "IpoSampleTable.csv")

IpoFiltered <- read.csv("IpoManualFilteredTable.csv", header = TRUE, sep = ",", fileEncoding = "UTF-8-BOM")
row.names(IpoFiltered) <- IpoFiltered$X
IpoFiltered$X = NULL
IpoFiltered$Total <- rowSums(IpoFiltered)
IpoFiltered$LogTotal <- log(IpoFiltered$Total)
write.csv(IpoFiltered, file = "IpoFilterLog.csv")

#Nectar Production#
IpoNectar<- read.csv("IpoDataSheet.csv", header = TRUE, sep = ",", fileEncoding = "UTF-8-BOM")
names(IpoNectar)
LogTotalPeakArea<- IpoNectar$LogTotal
IpoBiomass <- IpoNectar$Biomass
IpoVWC <- IpoNectar$VWC
IpoSites <- IpoNectar$Site
dim(IpoNectar)
names(IpoNectar)

NectarProduction <- IpoNectar
NectarLog <- NectarProduction$LogTotal
NectarSites <- NectarProduction$Site
Nectar <- NectarProduction$X24Hr.Nectar.Production
plot(Nectar, NectarLog, xlab = "Nectar Production", ylab = "Log(TotalPeakArea)")
points(x = Nectar, y = NectarLog, col = as.integer(NectarSites), pch = 21, bg = as.integer(NectarSites))
legend("bottomright", legend=levels(IpoSites), fill= c("black", "red"))
abline(lm(NectarLog~Nectar), col="blue") # regression line (y~x) 
title( main = "Correlation of Nectar Production to Total Volatiles")
cor(Nectar, NectarLog)
IpoNectarCorLM <- summary(lm(LogTotal~X24Hr.Nectar.Production, data = NectarProduction))
IpoNectarCorLM
#Biomass#
plot(IpoBiomass, LogTotalPeakArea, xlab = "Biomass (g)", ylab = "Log(TotalPeakArea)")
points(x = IpoBiomass, y = LogTotalPeakArea, col = as.integer(IpoSites), pch = 21, bg = as.integer(IpoSites))
legend("bottomright", legend=levels(IpoSites), fill= c("black", "red"))
abline(lm(LogTotalPeakArea~IpoBiomass), col="blue") # regression line (y~x) 
title( main = "Correlation of Biomass to Total Volatiles")

IpoBiomassCorLM <- summary(lm(LogTotal~IpoBiomass, data = IpoNectar))
IpoBiomassCorLM

#VWC#
plot(IpoVWC, LogTotalPeakArea, xlab = "Volumetric Water Content (%)", ylab = "Log(TotalPeakArea)")
points(x = IpoVWC, y = LogTotalPeakArea, col = as.integer(IpoSites), pch = 21, bg = as.integer(IpoSites))
legend("bottomright", legend=levels(IpoSites), fill= c("black", "red"))
abline(lm(LogTotalPeakArea~IpoVWC), col="blue") 
title( main = "Correlation of VWC to Total Volatiles")

IpoVWCCorLM <- summary(lm(LogTotal~IpoVWC, data = IpoNectar))
IpoVWCCorLM




library(ggplot2)
install.packages("ggpubr")
library(ggpubr)
IpoNectarCor <- ggplot(IpoNectar, aes(X24Hr.Nectar.Production, LogTotal))+
   geom_point()+
   labs(title = "Nectar Production Correlations to Total Volatiles", x = "Nectar Production", y = "Log(TotalPeakArea)")+
  stat_smooth(method=lm)
IpoNectarCor

IpoBiomassCor <- ggplot(IpoNectar, aes(Biomass, LogTotal))+
  geom_point()+
  labs(title = "Biomass Compared to Total Volatiles", x = "Biomass", y = "Log(TotalPeakArea)")+
  stat_smooth(method=lm)
IpoBiomassCor
  
IpoNectarCorLM <- summary(lm(X24Hr.Nectar.Production~LogTotal, data = IpoNectar))
IpoNectarCorLM
IpoBiomassCorLM <- summary(lm(Biomass~LogTotal, data = IpoNectar))
IpoBiomassCorLM

