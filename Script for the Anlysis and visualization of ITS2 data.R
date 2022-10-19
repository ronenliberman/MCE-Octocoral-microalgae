#### Script for "Symbiodiniaceae *ITS2* analyses"
author: "Ronen Liberman -- ronenliberman@gmail.com"

## This scripts is using 
# About this document


  #### All analyses preformed with R version `r getRversion()`.
##  This is the code accompanies the publication MS ## titled : "Octocorals in the Gulf of Aqaba exhibit high photosymbiont fidelity"
#If you download my entire accompanying github directory you should be able to re-run these analyses by
##following along with the code chunks in R Studio. If you download the code separately or you are 
##using this pipeline on your own data, you may need to change the working directory to where the associated files are housed 
  #(i.e. ```setwd("~/path/to/directory/with/data")```).   

#The data used for this analysis are calculated with *SymPortal* . 
#The raw Symbiodiniaceae *ITS2* sequences obtained from octocorals samples can be found in ### 
#Hopefully you are able to follow along with this file and find it useful to use with your own data!


## Loading required packages
#For the following analyses we will require the use of a number of different R packages.
if (!require("pacman")) install.packages("pacman")
pacman::p_load("ggplot2", "dplyr","plyr","ggsci", "reshape2", "MCMC.OTU", "vegan", "edgeR", "tidyverse",
               "gt", "vegan", "pairwiseAdonis", "RColorBrewer", "Redmonder")


#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
#library("phyloseq")

########################
if (!require("edgeR")){BiocManager::install("edgeR", update = FALSE) 
  library(edgeR)}

setwd("C:/Users/Ronen/Desktop/Symportal_Analysis_Liberman") #put here the folder with the supp. data

## Loading *ITS2* sequence data into R
its2RAW = read.delim("177_20210916_all_03_DBV_20210917T005252.seqs.absolute.abund_only.txt", header = TRUE, check.names = FALSE)
head(its2Seq)
View(its2RAW)

MetaData = read.csv("metadata.csv")


#merging metadata and absolute reads 
its2Seq= merge(MetaData,its2RAW,  by="sample_uid", all = T) #combine both df's

head(its2Seq)
#View(its2Seq)


its2Seq$Depthgroup= as.factor(its2Seq$Depthgroup)

# Working only on Caldocopium 

## Purging outlier sequences and normalizing reads

# Here we remove low abundance (< 0.001%) sequences and normalize - default settings 
#sequence counts with weighted trimmed mean of M-values (TMM; [Robinson and Oshlack 2010](https://doi.org/10.1186/gb-2010-11-3-r25){target="_blank"})
#. This helps to account for disparity in sequencing depth across libraries.<br>

cladoits <- its2Seq %>%  ##here I choose only host species with cladocopium 
  filter(Genus!="Litophyton")


goods = purgeOutliers(cladoits, count.columns = 7:length(its2Seq),otu.cut = 0.001)

#View(goods)
#head(its2Seq)
its2SeqTransposed = t(goods[, 7:length(goods[1, ])])
its2SeqList = DGEList(counts = its2SeqTransposed,remove.zeros = T)
head(its2SeqList$samples)
#Now we can use TMM normalization in *edgeR*

its2SeqNorm =  calcNormFactors(its2SeqList, method = "TMM")
head(its2SeqNorm$samples)
its2TMM = t(cpm(its2SeqNorm, normalized.lib.sizes = TRUE))
#bind with only cladocopium metadata 
Meta_clado = read.csv("metadata_cladocopium.csv") #using the metadata of cladocopium species only. This df is found in the supp folder

# View(Meta_clado)
its2SeqNorm = cbind(Meta_clado[,c(1:6)], its2TMM)
head(its2SeqNorm)


colOrder = order(colSums(its2SeqNorm[7:length(its2SeqNorm[1,])]), decreasing = FALSE) + 6

its2SeqPerc = cbind(its2SeqNorm[,c(1:6)], its2SeqNorm[,c(colOrder)])

its2SeqPerc$sum = apply(its2SeqPerc[, c(7:length(its2SeqPerc[1,]))], 1, function(x) {
  sum(x, na.rm = T)
})


#NORMALIZATION
its2SeqPerc = cbind(its2SeqPerc[, c(1:6)], (its2SeqPerc[, c(7:(ncol(its2SeqPerc)-1))] 
                                            / its2SeqPerc$sum))

#Now a quick sanity check. If this worked the sum of each row should = 100% (i.e. "1").

apply(its2SeqPerc[, c(7:(ncol(its2SeqPerc)))], 1, function(x) {
  sum(x, na.rm = T) 
})   

#I added an additional column to sort better for the stacked barplot. This was just a work around to get the ```facet_grid()``` function to play nice with our data. 
#I added a coulumn "barPlotOrder" and for each population I filled in a series 1:n foreach sample at each Site:Depth combo,
#so now there's no large blank expanses on the plot.


g = otuStack(its2SeqPerc, count.columns = c(7:length(its2SeqPerc[1, ]))
             ,
             condition.columns = c(1:6))[1:19521,]# remove summ rows
#View(g) 

## remove summ rows
Seq <- g %>% 
  filter(otu!="summ") %>% 
  mutate(Depthgroup = fct_relevel(Depthgroup, "Shallow", "MCE")) # level shallow first then MCE
Seq$uniq_id <- paste(Seq$sample,Seq$Depthgroup)

View(Seq)

CladoNorm <- its2SeqNorm ## this I will need later in order to do statistics on Host-species with Cladocopium 

#checking my results by looking at one genus - this is possible if you want to re-check data 

# Ovabunda <-  Seq %>% 
#  filter(Genus=="Ovabunda")
#View(Ovabunda)

###Construct ITS plot ### 

colorCount = length(c(7:length(its2SeqPerc[1,])))
getPalette = colorRampPalette(redmonder.pal(8, "qMSOGnYl"), bias = 1.7 )

#Ploting the ITS-type results for host species with Cladocopium 
its2SeqPlotA = ggplot(Seq, aes(x = sample, y = count*100, fill = as.factor(otu))) +
  geom_bar(position = "stack", stat = "identity", color = "black",
           size = 0.25) + 
  ylab("Relative abundance (%)") +
  # scale_colour_manual(values = c(mypal(12), mypal2(12))) +
  scale_fill_manual(values=rev(getPalette(colorCount)))+ 
  facet_wrap(Species~Depthgroup, scales = "free_x")  + 
  labs(fill = expression(paste(italic("ITS2"), " sequence"))) +
  guides(fill = guide_legend(ncol = 9, reverse = TRUE))+
  theme_bw()

its2SeqPlot = its2SeqPlotA +                                       ### adding unique themes
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "bottom",
        legend.justification = "left",
        legend.direction = "horizontal",
        legend.title = element_text(color = "black", size = 12, hjust = 0.5),
        legend.text = element_text(color = "black", size = 10),
        legend.key = element_blank(),
        legend.key.size = unit(0.4,"line"),
        legend.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        strip.text.x = element_text(size = 12, face = "italic"),
        strip.text.y = element_text(size = 12),
        strip.background = element_rect(fill = "white", size = 0.9)
  )
its2SeqPlot

Cladocopium_itsSeq_plot =  its2SeqPlot #save this plot as cladocopium plot 


################################################
###  preparing the Symbiodinium data and plot , Similar to the Cladocopium data. Here I subset for species with symbiodinium which is Litophyton only and the 
## purge for a cutoff of 0.001. Afterwards I normalize the counts. 

Symbio <- its2Seq %>% 
  filter(Genus=="Litophyton")
View(Symbio)

goods = purgeOutliers(Symbio, count.columns = 7:length(its2Seq),otu.cut = 0.001)
head(goods)
#View(goods)

its2SeqTransposed = t(goods[, 7:length(goods[1, ])])
its2SeqList = DGEList(counts = its2SeqTransposed,remove.zeros = T)
head(its2SeqList$samples)
#Now we can use TMM normalization in *edgeR*

its2SeqNorm =  calcNormFactors(its2SeqList, method = "TMM")
head(its2SeqNorm$samples)
its2TMM = t(cpm(its2SeqNorm, normalized.lib.sizes = TRUE))

#bind with only cladocopium metadata 
Meta_symbio = read.csv("metadata_symbiodinium.csv")
View(Meta_symbio)

its2SeqNorm = cbind(Meta_symbio[,c(1:6)], its2TMM)
head(its2SeqNorm)
#View(its2SeqNorm)

colOrder = order(colSums(its2SeqNorm[7:length(its2SeqNorm[1,])]), decreasing = FALSE) + 6

its2SeqPerc = cbind(its2SeqNorm[,c(1:6)], its2SeqNorm[,c(colOrder)])

its2SeqPerc$sum = apply(its2SeqPerc[, c(7:length(its2SeqPerc[1,]))], 1, function(x) {
  sum(x, na.rm = T)
})
#View(its2SeqPerc)

its2SeqPerc = cbind(its2SeqPerc[, c(1:6)], (its2SeqPerc[, c(7:(ncol(its2SeqPerc)-1))] 
                                            / its2SeqPerc$sum))

#Now a quick sanity check. If this worked the sum of each row should = 100% (i.e. "1").

apply(its2SeqPerc[, c(7:(ncol(its2SeqPerc)))], 1, function(x) {
  sum(x, na.rm = T) 
})   

#I added an additional column to sort better for the stacked barplot. This was just a work around to get the ```facet_grid()``` function to play nice with our data. 
#I added a coulumn "barPlotOrder" and for each population I filled in a series 1:n foreach sample at each Site:Depth combo,
#so now there's no large blank expanses on the plot.


g = otuStack(its2SeqPerc, count.columns = c(7:length(its2SeqPerc[1, ]))
             ,
             condition.columns = c(1:6))[1:19521,]# remove summ rows
#View(g) 

## remove summ rows
Seq <- g %>% 
  filter(otu!="summ") %>% 
  mutate(Depthgroup = fct_relevel(Depthgroup, "Shallow", "MCE")) # level shallow first then MCE

Seq$uniq_id <- paste(Seq$sample,Seq$Depthgroup)
levels(Seq$otu)

#Check data 

#View(Seq)

#View(its2SeqPerc)


write.csv(Seq, "Symbiodinium_filtered_dat.csv")   ##If you wish to export 

#Seq <- read.csv("Symbiodinium_dat.csv")

#levels(its2SeqPerc$Species)


###Construct ITS plot ### 

colorCount = length(c(7:length(its2SeqPerc[1,])))
getPalette = colorRampPalette(redmonder.pal(8, "qMSOOrRd"), bias = 1.7 )


levels(Seq$Depthgroup) = c ("Shallow", "MCE")  #relevel the depth with shallow first 


Seq <- Seq %>% 
  arrange(desc(grouping))
#levels(Seq$otu) = c("C1", "C1dw" , "C1b","C107c", "C3", "C3cq", "C3dx","C3it", "C3hx", "C3hy","C64a", "C65","C65a","C65d" ,"26673_C")
levels(Seq$otu)


its2SeqPlotA = ggplot(Seq, aes(x = uniq_id, y = count*100, fill = factor(otu))) +
  geom_bar(position = "stack", stat = "identity", color = "black",
           size = 0.25) + 
  ylab("Relative abundance (%)") +
  scale_fill_manual(values=rev(getPalette(colorCount)))+ facet_wrap(Species~Depthgroup, scales = "free_x")  + 
  labs(fill = expression(paste(italic("ITS2"), " type"))) +
  guides(fill = guide_legend(ncol = 10, reverse = TRUE))+
  theme_bw()

its2SeqPlot = its2SeqPlotA +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "bottom",
        legend.justification = "left",
        legend.direction = "horizontal",
        legend.title = element_text(color = "black", size = 10, hjust = 0.5),
        legend.text = element_text(color = "black", size = 7),
        legend.key = element_blank(),
        legend.key.size = unit(0.3,"line"),
        legend.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        strip.text.x = element_text(size = 12, face = "italic"),
        strip.text.y = element_text(size = 12),
        strip.background = element_rect(fill = "white", size = 0.9)
  )

its2SeqPlot

#for nephtheids 
Symbio_its_barplot=  its2SeqPlot

#export plots 
ggsave(Symbio_its_barplot, file="Symbio_ITS_barplot_0.1Perc.tiff", width = 12, height = 20, units = "cm",dpi=600)
ggsave(Cladocopium_itsSeq_plot, file="Cladocopium_ITS_type_0.01perc_barplot color_try.tiff", width = 20, height = 20,dpi=600, units = "cm")

################################################################################################
# ITS2 type profiles
# We can now look at the ITS2 type profiles predicted by SymPortal.


profRAW = read.delim("177_20210916_all_03_DBV_20210917T005252.profiles.absolute.abund_only.txt", header = TRUE, check.names = FALSE)
head(its2Seq)
View(profRAW)
MetaData = read.csv("metadata.csv")

p.Seq= merge(MetaData,profRAW,  by="sample_uid", all = T)
#edit the single Ovabunda profi
#profSeq=edit(p.Seq)

head(p.Seq)

its2ProfsTransposed = t(p.Seq[, 7:length(p.Seq[1, ])])

its2ProfsList = DGEList(counts = its2ProfsTransposed)
head(its2ProfsList$samples)

its2ProfsNorm =  calcNormFactors(its2ProfsList, method = "TMM")
head(its2ProfsNorm$samples)
its2TMM = t(cpm(its2ProfsNorm, normalized.lib.sizes = TRUE))
its2ProfsNorm = cbind(p.Seq[,c(1:6)], its2TMM)

head(its2ProfsNorm)
SymbioNorm <- its2ProfsNorm #save this for doing statistics later on

##Preparing ITS2 type profiles for plotting

colOrder2 = order(colSums(its2ProfsNorm[7:length(its2ProfsNorm[1,])]), decreasing = TRUE) + 6

its2ProfsPerc = cbind(its2ProfsNorm[,c(1:6)],its2ProfsNorm[,c(colOrder2)])
its2ProfsPerc$sum = apply(its2ProfsPerc[, c(7:length(its2ProfsPerc[1,]))], 1, function(x) {
  sum(x, na.rm = T)
})

its2ProfsPerc = cbind(its2ProfsPerc[, c(1:6)], (its2ProfsPerc[, c(7:(ncol(its2ProfsPerc)-1))] 
                                                / its2ProfsPerc$sum))
head(its2ProfsPerc)
##verifying that sum OTU of each sample = 1 (100%)
apply(its2ProfsPerc[, c(7:(ncol(its2ProfsPerc)))], 1, function(x) {
  sum(x, na.rm = T)
})

gssProf = otuStack(its2ProfsPerc, count.columns = c(7:length(its2ProfsPerc[1, ])),
                   condition.columns = c(1:6))[1:3133,] # remove summ rows
View(gssProf)
levels(gssProf$Species)
## remove summ rows
prf <- gssProf %>% 
  filter(otu!="summ") %>% 
  mutate(Depthgroup = fct_relevel(Depthgroup, "Shallow", "MCE")) %>%  # level shallow first then MCE
  mutate(Species = fct_relevel(Species, "Ovabunda sp.", "R. fulvum", "S. eilatensis","S. mesophotica","S. leptoclados", "S. vrijmoethi",  " L. arboreum  "," L. savigney")) # level shallow first then MCE

#prf$Species

View(prf)
colorCount2 = length(c(7:length(its2ProfsPerc[1,]))) +1
getPalette2 = colorRampPalette(redmonder.pal(8, "qPBI"))

#change the lables of the legend 
lab= read.csv("real_type profile.csv") # here I load the real names of the profiles. This is important because in the original datasheet, 
#every profile is numbered, but these numbers have names. they represent real biological entities, i.e profiles
head(lab)

#head(prf)
colnames(lab)[1]<- "otu"

mergd_for_lab= prfp.Seq= merge(lab,prf,  by="otu", all = T)
levels(mergd_for_lab$Depthgroup) = c ("Shallow", "MCE")
levels(mergd_for_lab$Depthgroup)
view(mergd_for_lab)
colnames(mergd_for_lab)[5]<- "sample"
head(mergd_for_lab)
###plot
library(ggsci)

its2ProfsPlotA = ggplot(mergd_for_lab, aes(x = sample, y = count*100, fill = factor(ITS2.type.profile))) +
  geom_bar(position = "stack", stat = "identity", color = "black", size = 0.25) + 
  ylab("Relative abundance (%)") +
  # scale_fill_jco()+
  scale_fill_manual(values = getPalette(colorCount2))+
  # scale_fill_manual(values=getPalette(colorCount))+
  #scale_fill_manual(values = c(getPalette2(colorCount2)[2:12], 
  #                           carto_pal(n = 7, "Tropic")[c(6,7)]))+
  labs(fill = expression(paste(italic("ITS2"), " type profile"))) +
  guides(fill = guide_legend(ncol = 1, reverse = FALSE)) +
  #scale_fill_manual(values=rev(getPalette(colorCount)))+
  facet_wrap(Species~Depthgroup, scales = "free_x",nrow= 2) +
  theme_bw()
its2ProfsPlotA

its2PRFPlot = its2ProfsPlotA +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12),
        legend.position = "right",
        legend.title = element_text(color = "black", size = 10, hjust = 0.5),
        legend.text = element_text(color = "black", size = 8),
        legend.key = element_blank(),
        legend.key.size = unit(0.3,"line"),
        legend.background = element_blank(),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_blank(),
        strip.text.x = element_text(size = 12, face = "italic"),
        strip.text.y = element_text(size = 12),
        strip.background = element_rect(fill = "white", size = 0.9)
  )

its2PRFPlot
#export plot
ggsave(its2PRFPlot, file="ITS_Profile.barplot_legend_right.tiff", width = 20, height = 20,dpi=600, units = "cm")
 
 ###################################################
 
#Statistical analyses on ITS2 communities  
 


 #Im using the same database as for the ITS-type - the normalized one - check above 'its2SeqNorm'
 
 #First I'm using PERMANOVA to test if there is a difference at the 'community'all species' level
 # across depth or between species - and for this im using a nested-Permnova 
 
 #Here I Work only on Cladocopium host species - I do not mix Clades! 
 #Here im using the normalized its2seq generated above 
 
 #S. mesophotica is a depth-specialist, found only in the MCE and is unbalancing - i want to filter it out
 
 itsNomesophotica=  its2SeqNorm %>% 
   #filter(Genus!="Litophyton") %>% 
   filter(Species!="S. mesophotica")
 #check
 itsNomesophotica$Species
 
 set.seed(694) #setting seed allows repetition of randomized processes
 
 #So, is there a difference across depth? 
 #perm <- how(nperm = 9999)
 setBlocks(perm) <- with(itsNomesophotica, Species)
 
 its2Adonis_nest_depth = adonis2(itsNomesophotica[, c(7:ncol(itsNomesophotica))] ~ Depthgroup ,
                                 data = itsNomesophotica,  method = "bray")
 its2Adonis_nest_depth
#No  p> 0.1 - this is part of TABLE 2
 
 #perm <- how(nperm = 9999)
 setBlocks(perm) <- with(itsNomesophotica, Depthgroup)
 
 its2Adonis_nest_species = adonis2(itsNomesophotica[, c(7:ncol(itsNomesophotica))] ~ Species ,
                                   data = itsNomesophotica, method = "bray")
 its2Adonis_nest_species
 
 #Yes! p= 0.001  this is part of TABLE 2 

#OK so overall there is no difference across depth but between species! 
 
 ## but, looking at my figures, few species seem to have differences across depth , therefore i must check for each species.
 
 #Working on each species Cladocopium host species separately to check diff. across depth 
 
 set.seed(694)
 
 its_ovabunda=  its2SeqNorm %>% 
   filter(Species=="Ovabunda sp.")   

 ova_Adonis = adonis(its_ovabunda[, c(7:ncol(its_ovabunda))] ~ Depthgroup , 
                     data = its_ovabunda, method = "bray")
 
 ova_Adonis
 
 #Pairwise PERMANOVA for multiple comparisons
 
 ova_its2Pairadonis = data.frame(pairwise.adonis(its_ovabunda[, c(7:ncol(its_ovabunda))],
                                                 factors = its_ovabunda$Depthgroup,
                                                 sim.method = "bray", p.adjust.m = "BH"))
 
 ova_its2Pairadonis_df <- ova_its2Pairadonis %>% 
   mutate(species= its_ovabunda$Species[1],n= NROW(its_ovabunda))#add species name and n replicates = METADATA
 head(ova_its2Pairadonis_df)
 
 ########## LEPTOCLADOS#####################
 its_lepto=  its2SeqNorm %>% 
   filter(Species=="S. leptoclados")   
 
lepto_its2Adonis = adonis(its_lepto[, c(7:ncol(its_lepto))] ~ Depthgroup , 
                           data = its_lepto, permutations = 9999, method = "bray")
lepto_its2Adonis

#Pairwise PERMANOVA for multiple comparisons

lepto_its2Pairadonis = data.frame(pairwise.adonis(its_lepto[, c(7:ncol(its_lepto))],
                                                  factors = its_lepto$Depthgroup,
                                                  sim.method = "bray", p.adjust.m = "BH"))

lepto_its2Pairadonis_df <- lepto_its2Pairadonis %>% 
  mutate(species= its_lepto$Species[1],n=NROW(its_lepto))
head(lepto_its2Pairadonis_df)
 #######################
##RHYTISMA FULVUM

its_rhy=  its2SeqNorm %>% 
  filter(Species=="R. fulvum")   

rhy_its2Adonis = adonis(its_rhy[, c(7:ncol(its_rhy))] ~ Depthgroup , 
                          data = its_rhy,  method = "bray",p.adjust.m = "BH")

rhy_its2Adonis
#Pairwise PERMANOVA for multiple comparisons

rhy_its2Pairadonis = data.frame(pairwise.adonis(its_rhy[, c(7:ncol(its_rhy))],
                                                factors = its_rhy$Depthgroup,
                                                sim.method = "bray", p.adjust.m = "BH"))

rhy_its2Pairadonis_df <- rhy_its2Pairadonis %>% 
  mutate(species= its_rhy$Species[1], n=NROW(its_rhy))
head(rhy_its2Pairadonis_df)

################
#SINLARIA VRIJMOETHI

its_vrij=  its2SeqNorm %>% 
  filter(Species=="S. vrijmoethi")   

vrij_its2Adonis = adonis(its_vrij[, c(7:ncol(its_vrij))] ~ Depthgroup , 
                         data = its_vrij, method = "bray")

vrij_its2Adonis
#Pairwise PERMANOVA for multiple comparisons

vrij_its2Pairadonis = data.frame(pairwise.adonis(its_vrij[, c(7:ncol(its_vrij))],
                                                 factors = its_vrij$Depthgroup,
                                                 sim.method = "bray", p.adjust.m = "BH"))

vrij_its2Pairadonis_df <- vrij_its2Pairadonis %>% 
  mutate(species= its_vrij$Species[1],n = NROW(its_vrij))
head(vrij_its2Pairadonis_df)

#Eilatensis
 its_eilat=  its2SeqNorm %>% 
   filter(Species=="S. eilatensis")   
 
 eilat_its2Adonis = adonis(its_eilat[, c(7:ncol(its_eilat))] ~ Depthgroup , 
                           data = its_eilat, method = "bray")
 
 eilat_its2Adonis

 eilat_its2Pairadonis = data.frame(pairwise.adonis(its_eilat[, c(7:ncol(its_eilat))],
                                                   factors = its_eilat$Depthgroup,
                                                   sim.method = "bray", p.adjust.m = "BH"))
 
 eilat_its2Pairadonis_df <- eilat_its2Pairadonis %>% 
   mutate(species= its_eilat$Species[1], n = NROW(its_eilat))
 head(eilat_its2Pairadonis_df)
 
 
 #Pairwise PERMANOVA for multiple comparisons shows that only iN
 #S. eilatensis and litophyton sp. there is a sig. diff across depth 
 
 #combine all the pair.adonis results to one df
 head(rhy_its2Pairadonis_df)
 
Table_2 <- rbind(ova_its2Pairadonis_df, rhy_its2Pairadonis_df,lepto_its2Pairadonis_df, vrij_its2Pairadonis_df, eilat_its2Pairadonis_df,Symbio_its2Pairadonis_df)
 View(Table_2) 
 
 write.csv(file="pairwise.adonis.csv", Table_2)## export and add the all species (Nested PERMANOVA) calculations
 
 all.adonis <- read.csv("all.adonis.csv")# from here i took the data for each individual species in table 2 
 
 head(all.adonis)
 
 # Export table of pairwise permanova comparisons as a picture file, possible, not part of the MS.
 table1 <- all.adonis %>% gt %>% 
   tab_header(
     title = md("")) %>% 
   fmt_number( # A column (numeric data)
     columns = SumsOfSqs, 
     decimals = 3 
   ) %>%
   fmt_number( 
     columns =  F.Model, 
     decimals = 3 
   ) %>%  
   fmt_number( 
     columns =  R2, 
     decimals = 3 
   ) %>% 
   fmt_number( 
     columns =   p.value, 
     decimals = 3 
   ) %>% 
   cols_label( SumsOfSqs =" Sum of Squeres") %>% 
   
   cols_hide( # Then hide variables
     columns = c(p.adjusted))  %>% 
   tab_style(
     style = list(
       cell_text(style = "italic")
     ),
     # different location
     cells_body(
       # entire column
       columns = Compared.taxa)) %>% 
   cols_align(
     align = "center",
     columns = c(F.Model,SumsOfSqs, p.value,R2) #at the moment KEEP (Statisticla.test , Compared.taxa ,Compared.variable) = left 
   )
 table1 %>% 
   gtsave(
     "Results of pairwise permanova comparisons across depth.png") #export to folder
 table1
 
 # Now lets see from where the significant between species comes from - for that i need to do a pairwise comparison between
 # all Cladocopium hosting species , and i can use the dataset including S. mesophotica becuase depth is not a factor here.
 #Well, depth could be but nested pairwise PERMANOVA was not successful for me. 
 
 
 set.seed(694)
 
 Table_3 = pairwise.adonis(its2SeqNorm[, c(7:ncol(its2SeqNorm))],
                                       factors = c(its2SeqNorm$Species),
                                       sim.method = "bray", p.adjust.m = "BH", perm = 999)
 
 Table_3 
 #So there are significant differences almost between all species 
 
 #Structure of Table 3 ## I can export it , but its not part of the MS because of table editing restrictions 
 Table_3_export <- itsClado_Pairadonis %>% gt %>% 
   tab_header(
     title = md("")) %>% 
   fmt_number(columns = c(SumsOfSqs, F.Model, p.adjusted), decimals = 3) %>% 
   cols_label( SumsOfSqs =" Sum of Squeres") %>% 
   tab_style(
     style = list(
       cell_text(style = "italic")
     ),
     # different location
     cells_body(
       # entire column
       columns = pairs)) %>% 
   cols_align(align = "center",
              columns = c(F.Model,SumsOfSqs, p.value,R2)) %>% 
   cols_hide( # Then hide variables
     columns = c(sig))%>%
   cols_label(pairs= md("Pairs"), Df= md("df")) %>% 
   cols_move(columns = Df,
             after = R2)%>% 
   gtsave(
     "PairwiseAdonis_only_Cldocopium_with mesophotica.png")
 
 #Can i use similar comparison on the 2 Litophyton species?  
 
 ############ Symbiodinium species ########### 
 #Working on the db Symbio that was determined before - in the Graphics part 
 
 View(Symbio)
 #Here there is no nested because depthgroup =  2 species  
 
 Symbio_its2Adonis = adonis(Symbio[, c(7:ncol(Symbio))] ~ Depthgroup , 
                            data = Symbio, method = "bray")
 
 Symbio_its2Adonis
 
 set.seed(694)
 Symbio_its2Pairadonis = data.frame(pairwise.adonis(Symbio[, c(7:ncol(Symbio))],
                                                    factors = Symbio$Depthgroup,
                                                    sim.method = "bray", p.adjust.m = "BH"))
 
 Symbio_its2Pairadonis_df <- Symbio_its2Pairadonis %>% 
   mutate(species=  Symbio$Genus[1], n=NROW(Symbio))
 head(Symbio_its2Pairadonis_df)
 
 #p = 0.0273 It is also significant ' but here the difference is either between species or depth - 1 factor actually
 ####
 
 ##### So there are differences between many species - But which ITS Types contribute the most? 
 #This can be found using Similarity percentage test (SIMPER) 
 
 SimperSpecies = simper(sqrt(its2SeqNorm[, c(7:ncol(its2SeqNorm))]), its2SeqNorm$Species,permutations = how())
 
 # creating simper output tables
 # remove $overall (only 1 element, not of the same length as other vectors in list)
 OTU.simper.df <- lapply(
   SimperSpecies,
   function(x) {
     do.call(data.frame,
             x[-3])
   }
 )
 
 #arrange top contributors in a table 
 fs_by_factor <- ldply(OTU.simper.df)#combining a list of df
 
 #change col names to something that make sense
 colnames(fs_by_factor)[1] = "Species_compared"
 colnames(fs_by_factor)[2] = "ITS2 type"
 
 #filter the top contributes
 Supp_Table_1 <- fs_by_factor %>% 
   filter(p<0.05&cusum<0.8&average>0.02) %>%  
   arrange(Species_compared, desc(cusum))
 str(Supp_Table_1)
 View(Supp_Table_1) #to view all the SIMPER resutls
 write.csv(file="pariwise_simper top contributer.csv", top_ave) #export as an excel doc.
 
 # Export table 
table_S1 <- Supp_Table_1 %>% gt %>% 
   tab_header(
     title = md("Top ITS2 types significantly contributing to pairwise differences")) %>% 
  fmt_percent( # A column (numeric data)
    columns = average,
    decimals = 2 
   ) %>%
   fmt_number( # A column (numeric data)
     columns = cusum, 
     decimals = 3 
   ) %>% 
   cols_label(average= "Average contribution", cusum = "Cumulative contribution",
              ava= "Average abundance sp.1", avb="Average abundance sp.1") %>% 
   cols_hide( # Then hide variables
     columns = c(sd, ratio,ord,p)
   )%>% 
  cols_align(align = "center",
             columns = c(average,cusum, ava,avb))%>% 
   gtsave(
     "Table s1_Simper_cladocopium pairwise comparison.png") #export to folder
table_S1

### Simper on the Litophyton pairwise comparison 

symbio_Simper =simper(sqrt(Symbio[, c(7:ncol(Symbio))]), Symbio$Depthgroup,permutations = how())
View(eilat_Simper)

# creating simper output tables
# remove $overall (only 1 element, not of the same length as other vectors in list)
symbio.OTU.simper.df <- lapply(
  symbio_Simper,
  function(x) {
    do.call(data.frame,
            x[-3])
  }
)
symbio.OTU.simper.df

#arrange top contributors in a table 

library(plyr)
simper_by_factor <- ldply(symbio.OTU.simper.df)#combining a list of df

#change col names to something that make sense
colnames(simper_by_factor)[1] = "Species_compared"
colnames(simper_by_factor)[2] = "ITS2 type"

#filter the top contributers
Supp_Table_2 <- simper_by_factor %>% 
  filter(p<0.05&cusum<0.8&average>0.02) %>%  
  arrange(Species_compared, desc(cusum))
#str(top_ave)
Supp_Table_2

# Export table 
table_S2 <- Supp_Table_2 %>% gt %>% 
  tab_header(
    title = md("Top ITS2 type significantly contributing for Litophyton diff across depth")) %>% 
  fmt_percent( # A column (numeric data)
    columns = average,
    decimals = 2 
  ) %>%
  fmt_number( # A column (numeric data)
    columns = cusum, 
    decimals = 3 
  ) %>% 
  cols_label(average= "Average contribution", cusum = "Cumulative contribution",
             ava= "Average abundance in L. arboreum", avb="Average abundance in L. savingey") %>% 
  cols_hide( # Then hide variables
    columns = c(Species_compared, sd, ratio,ord,p)) %>% 
  cols_align(align = "center",
             columns = c(average,cusum, ava,avb)) %>% 
  cols_align(align = "center",
             columns = c(average,cusum, ava,avb)) %>% 
  gtsave(
    "Table S2_Simper_litophyton.png") #export to folder

##########

#I can also calculate this way what are the ITS2 that affect the difference in S. eilatensis along depth  
View(its_eilat)

 eilat_Simper =simper(sqrt(its_eilat[, c(7:ncol(its_eilat))]), its_eilat$Depthgroup,permutations = how())
 View(eilat_Simper)

 # creating simper output tables
 
 # remove $overall (only 1 element, not of the same length as other vectors in list)
 OTU.simper.df <- lapply(
   eilat_Simper,
   function(x) {
     do.call(data.frame,
             x[-3])
   }
 )
 OTU.simper.df
 
 #arrange top contributors in a table 
 
 library(plyr)
simper_by_factor <- ldply(OTU.simper.df)#combining a list of df
 
 #change col names to something that make sense
 colnames(simper_by_factor)[1] = "Species_compared"
 colnames(simper_by_factor)[2] = "ITS2 type"
 
 #filter the top contributers
 Supp_Table_3 <- simper_by_factor %>% 
   filter(p<0.05&cusum<0.8&average>0.02) %>%  
   arrange(Species_compared, desc(cusum))
 str(top_ave)
 Supp_Table_3
 
 # Export table 
table_S3<- Supp_Table_3 %>% gt %>% 
   tab_header(
     title = md("ITS2 types Significantly contributing for S. eilatensis diff. across depth")) %>% 
  fmt_percent( # A column (numeric data)
    columns = average,
    decimals = 2 
   ) %>%
   fmt_number( # A column (numeric data)
     columns = cusum, 
     decimals = 3 
   ) %>% 
   cols_label(average= "Average contribution", cusum = "Cumulative contribution",
              ava= "Average abundance in shallow", avb="Average abundance in MCE") %>% 
   cols_hide( # Then hide variables
     columns = c(Species_compared, sd, ratio,ord,p)) %>% 
  gtsave(
    "table_S3_Simper_S. eilatensis.png") #export to folder

table_S3


