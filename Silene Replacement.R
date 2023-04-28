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
library(grDevices)
library(crayon)

#No significant difference in sent composition
ordiplot(ACnmds.vol, type="n", xlim=c(-1,1), ylim=c(-1.5,1.5))
points(ACnmds.vol, display = "sites", col = ACYear.pal[as.integer(ACYear)], bg = alpha(ACYear.pal, alpha=0.5)[as.integer(ACYear)], pch = 21)
ordispider(ACnmds.vol, ACYear, col= ACYear.pal, lwd=2)
text(ACnmds.vol, labels=ACSpeciesYear, col = ACYear.pal[as.integer(ACSpeciesYear)])
ordiellipse(ACnmds.vol, ACSpeciesYear, kind = "sd", conf = 0.75, draw = "polygon", col =ACYear.pal)
text(ACnmds.vol, labels=ACSpeciesYear, col = ACYear.pal[as.integer(ACSpeciesYear)])
legend("topleft", legend=as.character(c("S. latifolia", "S. dioica", "S. vulgaris")), fill=alpha(ACYear.pal, alpha=0.5))
title(main = substitute(paste(italic("Arenaria congesta"))))

ACSpeciesYear
NewNames <- ACSpeciesYear
levels(NewNames) <- c(levels(NewNames), "SL", "SD", "SV")
NewNames[NewNames == 'AC9'] <- 'SD'
NewNames[NewNames == 'AC8'] <- 'SL'
NewNames[NewNames == 'AC7'] <- 'SV'
NewNames
ordiplot(ACnmds.vol, type="n", xlim=c(-1,1), ylim=c(-1.5,1.5))
ordispider(ACnmds.vol, ACYear, col= ACYear.pal, lwd=2)
legend("topleft", legend=as.character(c("S. vulgaris", "S. latifolia", "S. dioica" )), text.font = 3, fill=alpha(ACYear.pal, alpha=0.5))
ordiellipse(ACnmds.vol, ACSpeciesYear, kind = "sd", conf = 0.75, draw = "polygon", col =ACYear.pal)
text(ACnmds.vol, labels=NewNames, col = "black", cex = 1.2)
text(ACnmds.vol, labels=NewNames, col = ACYear.pal[as.integer(ACSpeciesYear)])
title(main = expression("Scent Composition of Three"~italic("Silene")~"Species"))

#Significant depiction of pollinator selection
ESSpeciesYear
NewNames2 <- ESSpeciesYear
levels(NewNames2) <- c(levels(NewNames2), "SL", "SD", "SV")
NewNames2[NewNames2 == 'ES9'] <- 'SD'
NewNames2[NewNames2 == 'ES8'] <- 'SL'
NewNames2[NewNames2 == 'ES7'] <- 'SV'
NewNames2

ordiplot(ESnmds.vol, type="n", xlim=c(-1,1), ylim=c(-1.5,1.5))
ordispider(ESnmds.vol, ESYear, col= c("blue", "green", "red"), lwd=2)
ordiellipse(ESnmds.vol, ESSpeciesYear, kind = "sd", conf = 0., draw = "polygon", col =alpha(c("blue", "green", "red"), alpha=0.5))
legend("topleft", legend=as.character(c("Male", "Hermaphrodite", "Female" )), fill=alpha(c("blue", "green", "red"), alpha=0.5))

text(ESnmds.vol, labels=NewNames2, col = "black")
title(main = expression("VOC Emissions in Three Sexes of"~italic("Silene")))

#Attempt 2
GBSpeciesYear
NewNames3 <- GBSpeciesYear
levels(NewNames3) <- c(levels(NewNames3), "SL", "SD", "SV")
NewNames3[NewNames3 == 'GB9'] <- 'SD'
NewNames3[NewNames3 == 'GB8'] <- 'SL'
NewNames3[NewNames3 == 'GB7'] <- 'SV'
NewNames3
ordiplot(GBnmds.vol, type="n", xlim=c(-1,1), ylim=c(-1,1))
ordispider(GBnmds.vol, GBYear, col= c("yellow", "pink", "magenta"), lwd=2)
ordiellipse(GBnmds.vol, GBSpeciesYear, kind = "sd", conf = 0.75, draw = "polygon", col =alpha(c("yellow", "pink", "magenta"), alpha = 0.75))
legend("topleft", legend=as.character(c("S. vulgaris", "S. latifolia", "S. dioica" )), fill=alpha(c("yellow", "pink", "magenta"), alpha = 0.75))
text(GBnmds.vol, labels=NewNames3, col = "black")
title(main = expression("Scent Composition of Three"~italic("Silene")~"Species"))
adonis(GBVol ~ GBYear)

#Fake Physical Characteristics

morpho <- read.table("FakeMorphoSilene.csv", header = TRUE, sep = ",", fileEncoding = "UTF-8-BOM")
morpho$Colors <- as.character(morpho$Colors)
plot(morpho$CorollaWidth, morpho$TotalEmissions, xlab = "Corolla Width (cm)", ylab = "Total Emissions (ng)", cex = 1.5)
points(x = morpho$CorollaWidth, y = morpho$TotalEmissions, col = "black", pch = 21, bg = alpha(morpho$Colors, alpha = .75), cex = 1.5)
legend("bottomright", legend=as.character(c("S. latifolia", "S. dioica", "S. vulgaris" )), fill= alpha(c("pink", "magenta", "yellow"), alpha = 0.75), text.font = 3)

abline(lm(morpho$TotalEmissions~morpho$CorollaWidth), col="blue") # regression line (y~x) 
title(main = "Correlation of Corolla Width to Total Emissions")

fakeboy <- summary(lm(TotalEmissions~CorollaWidth, data = morpho))
fakeboy
