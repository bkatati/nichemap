######################
# Manuscript Title:
# Niche Partitioning Association of Fungal Genera Correlated with
# Lower Fusarium and Fumonisin-B1 levels in Maize
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Section - [A]

#Retrieval of Internal Mycobiome Amplicon Sequence Variants (ASVs) from phyloseq object:

require(plyr)
ps <- readRDS("https://github.com/bkatati/nichemap/blob/main/phyloseq.rds")
# NB: if file path error occurs, download rds file "phyloseq" from site ""https://github.com/bkatati/nichemap" 
# On your PC, create appropriate local drive path for the file and change above file path.

psg0 <- tax_glom(ps, "Genus")
psg1 <- transform_sample_counts(psg0, function(OTU)100* OTU / sum(OTU))
otu_table(psg1) <- t(otu_table(psg1))
OTUg <- otu_table(psg1)
TAXg <- tax_table(psg1)[, "Genus"]
mycobiome <- merge(TAXg, OTUg, by=0, all = TRUE)
mycobiome$Row.names = NULL
mycobiome$Mean=rowMeans(mycobiome[,-c(1)], na.rm=TRUE)
mycobiome <- mycobiome[order(desc(mycobiome$Mean)),]
head(mycobiome)

# Save the ASVs as 'csv' file:
# write.csv(mycobiome, "C:/~your_path/InternoBiome.csv", row.names = F)

# The csv file contains both external and internal mycobiome ASVs of 2018/2019 sampling season (S1).
# Transpose the "csv" and delete all external mycobiome fields (letter 'E' in sample code).
# Last two digits are field number: e.g. "KASI01" is field number-1 internal mycobiome, whereas
# KASE01 is field number-1 external mycobiome.


# INTERNAL MYCOBIOME CENSUS - AMPLICON SEQUENCE VARIANTS:

# datai <- read.csv("https://github.com/bkatati/nichemap/blob/main/InternoBiome.csv")
# NB: if file path error occurs, download csv file "InternoBiome.csv" from site ""https://github.com/bkatati/nichemap" 
# On your PC, create appropriate local drive path for the file and change above file path.
##########

# Section - [B]

 ######################################################
 # SUPPLEMENTAL TABLE S2 - Correlation Coefficients ###
 ######################################################
 
# Fungal General Relative Abundance Correlations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

head(datai)
dim(datai)

# Note1: to create the the Spearman Correlation Matrix (Figure 2), skip 'databi' below" and go to # CORRELATION MATRIX."
# Note2: to generate the Spearman Correlation Coefficients per district (Supplemental Table S2) use the 'Subset Section'
# in the string before proceeding to # CORRELATION MATRIX" e.g:
# "District = Livingstone," will generate coefficients for Livingstone as a sampled district:

# Subset Section:
#***************************************

# databi <- subset(datai, variable == "wet")
# head(databi)
# dim(databi)
# datai <- databi
#***************************************

# CORRELATION MATRIX

# Create Correlation Coefficients (top 30 genera):

cormat <- cor(datai[,5:34], method = c("spearman"))
head(cormat)
dim(cormat)

require(reshape2)
melted_cormat <- melt(cormat)
head(melted_cormat)

# To generate only half (lower triangle of correlation matrix) coefficients and remove upper redundancy:
get_lower_tri <- function(cormat) {
  cormat[lower.tri(cormat)]<-
    NA
  return(cormat)
}
lower_tri <- get_lower_tri(cormat)
lower_tri

# Melt the correlation matrix to drop 'NA'
melted_cormat <- melt(lower_tri, na.rm = TRUE)

# Write the correlation coefficients for Supplemental Table S2:

# write.csv(melted_cormat, file="C:/~your_path/rho-spearman.csv", row.names = F)

# Check Significance of the Coefficients (using all genera):

cormatall <- cor(datai[,5:86], method = c("spearman"))
dim(cormatall)

require(Hmisc)
coefp <- rcorr(as.matrix(cormatall), type=c("spearman"))
myPval <- coefp$P
# View(myPval)

# You may write the p-values:
# write.csv(myPval, file="C:/~your_path/p-Spearman.csv", row.names = T)
# alternative:
# write.csv(coefp[["P"]], file="C:/~your_path/p-Spearman.csv", row.names = T)


##############################################################
# FIGURE 2  - Spearman Correlation Matrix of Fungal Genera ###
##############################################################


# You may reorder the correlation to visualise patterns better.
reorder_cormat <- function(cormat) {
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <- cormat[hc$order, hc$order]
}

# Reorder
cormat2 <- reorder_cormat(cormat)
lower_tri <- get_lower_tri(cormat2)

# melt the correlation matrix:
melted_cormat <- melt(lower_tri, na.rm = TRUE)

# Create heatmap with ggheatmap:
require(ggplot2)
q2 <- ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))
pp <- q2 + geom_tile(colour = "white") + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1), 
                       space = "Lab",
                       name =   "Spearman's Correlation Matrix
(Internal Mycobiome - Overall)")

qp <- pp + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                                              size = 9, hjust = 1)) + 
  theme(axis.text.y = element_text(angle = 0, vjust = 1, size = 9, hjust = 1)) +
  theme(text = element_text(size = 11)) +
  coord_fixed()

# Matrix Coefficient Plot
#********************
gh <- qp + geom_text(aes(Var2, Var1, label = round(value, digits = 1)), colour = "black", size = 2.5)
#********************

# Note3: To remove numbers (coefficients) from the plot you're generating, use the below "gh <- qp":
# gh <- qp

hm <- gh + theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
                 panel.grid.major = element_blank(), panel.border = element_blank(),
                 panel.background = element_blank(), axis.ticks = element_blank(), 
                 legend.position = c(0.38, .9),
                 legend.direction = "horizontal") + guides(fill = guide_colorbar(barwidth = 15, barheight = 1,
                                                                                 title.position = "top",
                                                                                 title.hjust = 0.5))
hm + theme(axis.text = element_text(face = "italic"))



#############################################################
# FIGURE 3 - Multi-Dimensional Scaling, Agronomic Factors  ##
#############################################################

# Read the dataframe of agronomic factors
mycobiome<-read.csv("https://github.com/bkatati/nichemap/blob/main/dataframe_FusB1.csv")
# NB: if file path error occurs, download csv file "dataframe_FusB1.csv" from site "https://github.com/bkatati/nichemap" 
# On your PC, create appropriate local drive path for the file and change above file path.

head(mycobiome)
dim(mycobiome)

# Generate matrix from dataframe
mycobiome.matrix<-as.matrix(mycobiome[,7:12])
head(mycobiome.matrix)

# Transform matrix to square roots to minimize influence of most abundant genera:

mycobiome.mat<-sqrt(mycobiome.matrix)
head(mycobiome.mat)

require(vegan)
set.seed(79) #for reproducible results

mycobiomeMDS<-metaMDS(mycobiome.mat, distance="bray", k=2, trymax=35, autotransform=TRUE) #k = number of dimensions
mycobiomeMDS

# Note4: Stress value < 0.2, is acceptable

library(ggplot2); packageVersion("ggplot2")
MDS1 <- mycobiomeMDS$points[,1]
MDS2 <- mycobiomeMDS$points[,2]
mycobiome.plot<-cbind(mycobiome, MDS1, MDS2)
head(mycobiome.plot)

require(dplyr) 
fit<-envfit(mycobiomeMDS, mycobiome.mat)
arrow<-data.frame(fit$vectors$arrows,R = fit$vectors$r, P = fit$vectors$pvals)
arrow$Agronomic_Factors <- rownames(arrow)

# Collect only significant (p < 0.05 predictors, irrespective of level of contribution R):
arrow.p<-filter(arrow, P < 0.05, R > 0)
arrow.p

# Note5: you may write arrow values (P-values and R-sq):
# write.csv(arrow.p, file = "C:/~your_path/arrows.csv", row.names = F)

set.seed(07) # run each time for reproducible coordinates
p <- ggplot(data=mycobiome.plot, aes(MDS1,MDS2)) + theme_classic() +
  geom_point(data=mycobiome.plot, aes(MDS1, MDS2, color=FB1_level), position=position_jitter(1)) +##separates overlapping points
  stat_ellipse(aes(fill=FB1_level), alpha=0.25,type='t',size = 1, segments = 360, level = 0.67, geom="polygon", show.legend = NA) +##changes shading on ellipses
  theme_classic() +
  geom_segment(data=arrow.p, aes(x=0, y=0, xend=NMDS1, yend=NMDS2, lty = Agronomic_Factors), arrow=arrow(length=unit(.99, "cm")*arrow.p$R), colour = "blue", size = 0.3, alpha = 1) ##add arrows (scaled by R-squared value)

# MDS for Fumonisin-B1:
set.seed(79)
mdsFB1 <- p + theme(legend.position = "right", text = element_text(size = 12), axis.title = element_text(size = 13),
                    axis.text = element_text(size = 10)) + annotate("text", x=(1.1), y=(-0.5), label=paste('Stress =',round(mycobiomeMDS$stress,3))) +
  annotate("text", x=(-1.1), y=(0.35), label=paste('early')) + annotate("text", x=(1.2), y=(0.5), label=paste('late')) +
  annotate("text", x=(-0.2), y=(-1.1), label=paste('pest')) + annotate("text", x=(0.3), y=(-1.1), label=paste('med')) +
  labs(title = expression(paste("b) MDS of Agronomic Factors that may Influence FB1 levels in Maize")))+
  theme(title =element_text(size=10))
mdsFB1


set.seed(07) # run each time for reproducible coordinates
p <- ggplot(data=mycobiome.plot, aes(MDS1,MDS2)) + theme_classic() +
  geom_point(data=mycobiome.plot, aes(MDS1, MDS2, color=Fus_Level), position=position_jitter(1)) +##separates overlapping points
  stat_ellipse(aes(fill=Fus_Level), alpha=0.25,type='t',size = 1, segments = 360, level = 0.67, geom="polygon", show.legend = NA) +##changes shading on ellipses
  theme_classic() +
  geom_segment(data=arrow.p, aes(x=0, y=0, xend=NMDS1, yend=NMDS2, lty = Agronomic_Factors), arrow=arrow(length=unit(.99, "cm")*arrow.p$R), colour = "blue", size = 0.3, alpha = 1) ##add arrows (scaled by R-squared value)


# MDS for Fusarium:
set.seed(79)
mdsFus <- p + theme(legend.position = "right", text = element_text(size = 12), axis.title = element_text(size = 13),
                    axis.text = element_text(size = 10)) + annotate("text", x=(1.2), y=(-0.5), label=paste('Stress =',round(mycobiomeMDS$stress,3))) +
  annotate("text", x=(-1.1), y=(0.35), label=paste('early')) + annotate("text", x=(1.2), y=(0.5), label=paste('late')) +
  annotate("text", x=(-0.2), y=(-1.1), label=paste('pest')) + annotate("text", x=(0.3), y=(-1.1), label=paste('med')) +
  labs(title = expression(paste("a) MDS of Agronomic Factors that may Influence", italic(" Fusarium"), " Abundance in Maize"))) +
  theme(title =element_text(size=10))
mdsFus

# Note6: You may change R values in "arrow.p" above to filter off certain agronomic factors as needed.

# Two-in-one plot:

require(patchwork)
set.seed(79)
mdsFus / mdsFB1 +
  plot_layout(heights = c(2,2))
 ######################################################
 # Influence of fungal genus abundance on FB1 levels ##
 ######################################################

# [a] Spearman rank correlation rho of genera relative abundances against FB1:

Corr <- read.csv("https://github.com/bkatati/nichemap/blob/main/Spearman.csv")
# NB: if file path error occurs, download csv file "Spearman.csv" from site "https://github.com/bkatati/nichemap" 
# On your PC, create appropriate local drive path for the file and change above file path.

head(Corr)

# i) Sarocladium-FB1
Spear1 <- cor.test(x=Corr$Sar, y=Corr$FB1, method = "spearman")
Spear1

# ii) Stenocarpella-FB1
Spear2 <- cor.test(x=Corr$Sten, y=Corr$FB1, method = "spearman")
Spear2

# ii) Fusarium-FB1
Spear3 <- cor.test(x=Corr$FusGib, y=Corr$FB1, method = "spearman")
Spear3

#######################################################################
## SUPPLEMENTAL DATA S3 - Fusarium v Sarocladium, External Mycobiome ##
#######################################################################

Corr2 <- read.csv("https://github.com/bkatati/nichemap/blob/main/S3_ExternoBiome.csv")

# NB: if file path error occurs, download csv file "ExternoBiome.csv" from site "https://github.com/bkatati/nichemap" 
# On your PC, create appropriate local drive path for the file and change above file path.

head(Corr2)
Spear4 <- cor.test(x=Corr2$sar, y=Corr2$fus, method = "spearman")
Spear4

###############################END####################################################
