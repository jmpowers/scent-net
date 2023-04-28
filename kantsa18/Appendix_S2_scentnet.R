##########################################################################################
##########################################################################################
######  Appendix S2: "paco: implementing Procrustean Approach to Cophylogeny in R"  ######
######      Code to undertake example paco analysis on the Mendoza pollination      ######   
######                network described in the application note.                    ######
##########################################################################################
##########################################################################################

# first we read in the data and load required packages
library(paco)
library(ape)
library(ggplot2)

# start the paco procedure
poll_phy <- polls.classtree$phylo
pla_phy <- as.phylo(hclust(vegdist(scent.net), "mcquitty"))
int <- net.lump
H <- cophenetic(polls.classtree$phylo)
P <- cophenetic(as.phylo(hclust(vegdist(scent.net), "mcquitty")))
D <- prepare_paco_data(H, P, net.lump)
D <- add_pcoord(D, correction='cailliez')
# now we are ready for cophylogenetic analysis
D <- PACo(D, nperm=1000, seed=13, method='quasiswap', symmetric=TRUE)
# and to investigate the contribution of individual links
res <- residuals_paco(D$proc)

# to visualise this we use the ape function cophyloplot weighted by interaction contribution
# first we must make a list out of the interaction matrix
assoc <- data.frame(pol=rownames(net.lump)[which(net.lump>0, arr.ind=TRUE)[,'row']], pla=colnames(net.lump)[which(net.lump>0, arr.ind=TRUE)[,'col']])
# to weight the interactions we use the cophylogenetic contribution transformed to best show
# the differences graphically
weight <- (res^-2)/50

mnet.lump <- melt(net.lump)
mnet.lump <- mnet.lump[mnet.lump$value>0,]
#cophylo with interaction strength
cophyloplot(pla_phy, poll_phy, mnet.lump[,1:2], lwd=sqrt(mnet.lump[,3])/4, show.tip.label=F,  col='steelblue', length.line=0, gap=-20, space=60)
#cophylo with residuals
cophyloplot(pla_phy, poll_phy, mnet.lump[,1:2], lwd=weight/2, show.tip.label=F, col='steelblue', length.line=0, gap=-20, space=60)

cophyloplot(poll_phy, pla_phy, assoc, show.tip.label=FALSE, use.edge.length=FALSE,
            lwd=weight, col='steelblue', length.line=0, gap=-20, space=60)


# to analyse the links further we will split the interactions based on the plot
cophy_int <- c(grep('oxalis', names(res)), grep('viola', names(res)))
# remove the one interaction of Oxalis compacta that is not with the Lepidopterans
cophy_int <- cophy_int[-grep('Anthidium', names(res[cophy_int]))]
noncophy <- res[-cophy_int]
cophy <- res[cophy_int]

# Visualise residuals of  interactions against the rest
noncophy_dat <- data.frame(res=noncophy)
cophy_dat <- data.frame(res=cophy)
ggplot(noncophy_dat, aes(x=res))+
geom_density(fill='grey70')+
theme_bw()+
geom_vline(data=cophy_dat, aes(xintercept=res), col='darkorange1')+
theme(panel.grid=element_blank())+
xlab('Procrustes residuals')+
ylab('Frequency')+
scale_x_continuous(limits=c(0.05, 0.2), expand=c(0.,0), breaks=c(0.05, 0.10, 0.15, 0.20))+
scale_y_continuous(expand=c(0,0), limits=c(0.00, 25))

# Welch's t-test to test the difference between cophylogenetic signal of interactions
ttest <- t.test(cophy, noncophy)

# visualise the difference with a box and whisker plot
dat <- rbind(data.frame(cophy=cophy, level='high'), data.frame(cophy=noncophy, level='rest'))
ggplot(dat, aes(x=level, y=cophy, fill=level))+
geom_boxplot(alpha=0.85)+
scale_fill_brewer(palette='Paired')+
scale_x_discrete(labels=c('L-M/O','Rest'))+
ylab('Procrustes residual')+
theme_bw()+
theme(
	axis.title.x=element_blank(),
	panel.grid=element_blank(),
	legend.position='none',
	axis.text.x=element_text(size=14),
	axis.title.y=element_text(size=14)
	)

# test the influence of degree on cophylogenetic signal
special <- data.frame(cophy.sig=res, pla_deg=NA, row.names=names(res))
f <- function(x) length(grep(strsplit(x, '-')[[1]][2], names(res)))
special$pla_deg <- sapply(rownames(special), f)

plant_degree <- summary(lm(cophy.sig ~ pla_deg, data=special))
