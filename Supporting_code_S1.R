#Supporting code S1

#'Neighbours, not consumers, drive local intraspecific phytochemical changes in two grassland species'

library(tidyverse); library(ggpubr); library(vegan); library(car); library(dendextend)

rm(list=ls())

##############################################################################

#Enemy damage

damage <- read.csv("enemydamage.csv", header=T, stringsAsFactors = T)
damage <- damage %>%
  mutate(logdamageI = log10(damageI + 1)) %>%
  mutate(logdamageF = log10(damageF + 1)) %>%
  mutate(totaldamage = damageI + damageF) %>%
  mutate(logtotaldamage = log10(totaldamage + 1))

mod1 <- lm(logdamageI ~ Treatment * Diversity, data=damage)
summary(mod1)

mod2 <- lm(logdamageF ~ Treatment * Diversity, data=damage)
summary(mod2)

mod3 <- lm(logtotaldamage ~ Treatment * Diversity, data=damage)
summary(mod3)

TukeyHSD(aov(mod1))
TukeyHSD(aov(mod2))
TukeyHSD(aov(mod3))

damage %>%
  group_by(Treatment) %>%
  summarise(mean(damageI))
damage %>%
  group_by(Treatment) %>%
  summarise(mean(damageF))
damage %>%
  group_by(Treatment) %>%
  summarise(mean(totaldamage))
damage %>%
  group_by(Diversity) %>%
  summarise(mean(damageI))
damage %>%
  group_by(Diversity) %>%
  summarise(mean(damageF))
damage %>%
  group_by(Treatment, Diversity) %>%
  summarise(mean(damageI))
damage %>%
  group_by(Treatment, Diversity) %>%
  summarise(mean(damageF))
damage %>%
  group_by(Treatment, Diversity) %>%
  summarise(mean(totaldamage))
damage %>%
  group_by(Species, Treatment, Diversity) %>%
  summarise(mean(damageI))
damage %>%
  group_by(Species, Treatment, Diversity) %>%
  summarise(mean(damageF))

#Plot the damage

totaldamage <- ggplot(damage, aes(x=Diversity, y=totaldamage, fill=Treatment)) +
  geom_boxplot(position=position_dodge(0.8)) +
  ylab("Total consumer damage (%)") + xlab("Competition Treatment") +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme_bw() + theme(legend.position="top") + 
  theme(axis.title.x=element_blank()) + 
  ylim(0, 12)
totaldamage

overallinsectdamage <- ggplot(damage, aes(x=Diversity, y=damageI, fill=Treatment)) +
  geom_boxplot(position=position_dodge(0.8)) +
  ylab("Insect damage (% leaf removed)") + xlab("Competition Treatment") +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme_bw() + theme(legend.position="none")  + 
  ylim(0, 12)
overallinsectdamage

overallfungaldamage <- ggplot(damage, aes(x=Diversity, y=damageF, fill=Treatment)) +
  geom_boxplot(position=position_dodge(0.8)) +
  ylab("Fungal damage (% leaf diseased)") + xlab("Competition Treatment") +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme_bw() + theme(legend.position="none")  + 
  ylim(0, 12)
overallfungaldamage

#####################################################################################

#NDVI analysis

ndvi <- read.csv("NDVI.csv", header=T, stringsAsFactors = T) %>%
  unite('Combo', Diversity:Treatment, remove=F)

mod1 <- lm(NDVI ~ Treatment*Diversity, data=ndvi)
anova(mod1)
##Response: NDVI
##                    Df  Sum Sq Mean Sq F value    Pr(>F)    
##Diversity            1 0.38428 0.38428 75.1038 2.209e-11 ***
##Treatment            1 0.07393 0.07393 14.4485 0.0004063 ***
##Diversity:Treatment  1 0.00025 0.00025  0.0496 0.8247799    
##Residuals           48 0.24560 0.00512 

ndvimeans <- ndvi %>%
  group_by(Treatment, Diversity) %>%
  dplyr::summarise(mean(NDVI), 1.96*sd(NDVI)/sqrt(length(NDVI))) %>%
  dplyr::rename(mean = `mean(NDVI)`, se = `1.96 * sd(NDVI)/sqrt(length(NDVI))`) %>%
  unite('Combo', Diversity:Treatment, remove=F)

#Pesticide treatment NDVI 1.17 times higher than control (1.14 for Lc, 1.13 for Ag) for monocultures
#Polyculture treatment NDVI 1.40 times higher than mono (for controls), 1.37 times higher (for pesticide)

points <- data.frame(Species=c('L.capitata', 'L.capitata', 'A.gerardii', 'A.gerardii'),
                     Treatment=c('Control', 'Pesticide', 'Control', 'Pesticide'),
                     NDVI=c(0.47, 0.539, 0.452, 0.509))

ndviplot <- ggplot() +
  geom_point(data=ndvi, aes(x=Treatment, y=NDVI, shape=Combo), 
             size=1.8, alpha=0.15, stroke=1.0, position = position_dodge(0.9), 
             colour="darkgreen", fill="darkgreen") +
  geom_pointrange(data=ndvimeans, aes(x=Treatment, y=mean, shape=Combo,
                                      ymin=mean-se, ymax=mean+se),
                  size=0.6, position=position_dodge(0.9), colour="darkgreen", fill="darkgreen")+
  theme_bw() +
  theme(legend.position = "top") +
  scale_shape_manual(values=c(1, 2, 21, 24)) +
  ylab("NDVI") + 
  theme(legend.position = "top") +
  geom_point(data=points, aes(x=Treatment, y=NDVI), size=2.4, 
             colour=c("#3d3d3dff", "#3d3d3dff", "#ff8b00ff","#ff8b00ff"), 
             shape=c(1, 2, 1, 2), stroke=1.0)
ndviplot

#Aboveground biomass

biomass <- read.csv("biomass.csv", header=T, stringsAsFactors = T) %>%
  unite('Combo', Diversity:Treatment, remove=F)

mod2 <- lm(Biomass ~ Treatment*Diversity, data=biomass)
anova(mod2)
##                    Df Sum Sq Mean Sq F value    Pr(>F)    
##Treatment            1  14326   14326  1.6707    0.2024    
##Diversity            1 296046  296046 34.5232 3.881e-07 ***
##Treatment:Diversity  1     49      49  0.0058    0.9399    
##Residuals           48 411613    8575 

biomassmeans <- biomass %>%
  group_by(Treatment, Diversity) %>%
  dplyr::summarise(mean(Biomass), 1.96*sd(Biomass)/sqrt(length(Biomass))) %>%
  dplyr::rename(mean = `mean(Biomass)`, se = `1.96 * sd(Biomass)/sqrt(length(Biomass))`) %>%
  unite('Combo', Diversity:Treatment, remove=F)

#Pesticide treatment biomass 1.62 times higher than control for monocultures, 1.14 times higher for poly
#Polyculture treatment biomass 3.88 times higher than mono (for controls), 2.72 times higher (for pesticide)

sixteenspbiomass <- read.csv("16spbiomass.csv", header=T, stringsAsFactors = T) %>%
  unite('Combo', Diversity:Treatment, remove=F)

sixteenspbiomassmeans <- sixteenspbiomass %>%
  group_by(Species, Treatment, Diversity) %>%
  dplyr::summarise(mean(Biomass), 1.96*sd(Biomass)/sqrt(length(Biomass))) %>%
  dplyr::rename(mean = `mean(Biomass)`, se = `1.96 * sd(Biomass)/sqrt(length(Biomass))`) %>%
  unite('Combo', Diversity:Treatment, remove=F)
sixteenspbiomassmeans[is.na(sixteenspbiomassmeans)] <- 0

biomassplot <- ggplot() +
  geom_point(data=biomass, aes(x=Treatment, y=Biomass, shape=Combo), 
             size=1.8, alpha=0.15, stroke=1.0, position = position_dodge(0.9), 
             colour="darkgreen", fill="darkgreen") +
  geom_pointrange(data=biomassmeans, aes(x=Treatment, y=mean, shape=Combo,
                                         ymin=mean-se, ymax=mean+se),
                  size=0.6, position=position_dodge(0.9), colour="darkgreen", fill="darkgreen")+
  geom_point(data=sixteenspbiomass, aes(x=Treatment, y=Biomass, shape=Combo, 
                                        colour=Species, fill=Species), 
             size=1.8, alpha=0.15, stroke=1.0, 
             position = position_dodge(0.4)) +
  geom_pointrange(data=sixteenspbiomassmeans, aes(x=Treatment, y=mean, shape=Combo, 
                                                  colour=Species, fill=Species,        
                                                  ymin=mean-se, ymax=mean+se),
                  size=0.6, position=position_dodge(0.4))+           
  theme_bw() +
  theme(legend.position = "top") +
  scale_shape_manual(values=c(1, 2, 21, 24)) +
  scale_color_manual(values=c("#ff8b00ff", "#3d3d3dff")) + 
  scale_fill_manual(values=c("#ff8b00ff", "#3d3d3dff")) +
  ylab("Aboveground biomass") + 
  theme(legend.position = "top") 
biomassplot

abovegroundv2 <- ggarrange(ndviplot, biomassplot, ncol=2, labels=c("a", "b"))
ggsave("abovegroundv2.svg", width=180, height=100, units="mm")

#####################################################################################

##OVERALL PCAs ACROSS BOTH SPECIES FOR PHYTOCHEMICAL COMPOSITION

#Aqueous extract

aq <- read.csv("aqueousdatanormalised2.csv", header=T)
aq$Sample <- as.factor(aq$Sample)
aq$Sample <- as.character.factor(aq$Sample)

variables <- aq %>%
  select(Sample, Species, Diversity, Treatment)

aq.pca <- prcomp(aq[ ,c(5:174)])
summary(aq.pca)
aq.pca.axes <- aq.pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var="Sample")
aq.pca.data <- full_join(variables, aq.pca.axes)

#Plot

aq.pca.data.mod <- aq.pca.data %>%
  unite('Combo', Diversity:Treatment, remove=F) %>%
  mutate(SpDiv = paste0(Species, "_", Diversity))

aq.fullplot.mod <- ggplot(aq.pca.data.mod, aes(x=PC1, y=PC2, colour=Species, shape=Combo, 
                                               fill=Species)) +
  geom_point(size=2.5, alpha=0.7, stroke=1.0) + 
  scale_colour_manual(values=c("#ff8b00ff", "#3d3d3dff")) + 
  scale_fill_manual(values=c("#ff8b00ff", "#3d3d3dff")) +
  scale_shape_manual(values=c(1, 2, 21, 24)) +
  xlab("Aqueous PC1 (65.3%)") + ylab("Aqueous PC2 (17.8%)") + 
  theme_classic() +
  theme(legend.position = "top") +
  xlim(-3000, 5000) + ylim(-3000, 2000)
aq.fullplot.mod

#######################################

#Plot these PCA coordinates but for each species individually

aq.pca.data.andro <- aq.pca.data %>% 
  filter(Species=="Andropogon gerardii") %>%
  unite('Combo', Diversity:Treatment, remove=F)

aq.fullplot.andro <- ggplot(aq.pca.data.andro, aes(x=PC1, y=PC2, colour=Combo, shape=Diversity)) +
  geom_point(size=2.5, alpha=0.7) + 
  scale_colour_manual(values=c("#ffb100ff", "#ff8b00ff", "#838383ff", "#3d3d3dff")) + 
  xlab("Aqueous PC1") + ylab("Aqueous PC2") + 
  theme_classic() +
  theme(legend.position = "top") 
aq.fullplot.andro

aq.pca.data.lespe <- aq.pca.data %>% 
  filter(Species=="Lespedeza capitata") %>%
  unite('Combo', Diversity:Treatment, remove=F)

aq.fullplot.lespe <- ggplot(aq.pca.data.lespe, aes(x=PC1, y=PC2, colour=Combo, shape=Diversity)) +
  geom_point(size=2.5, alpha=0.7) + 
  scale_colour_manual(values=c("#ffb100ff", "#ff8b00ff", "#838383ff", "#3d3d3dff")) + 
  xlab("Aqueous PC1") + ylab("Aqueous PC2") + 
  theme_classic() +
  theme(legend.position = "top") 
aq.fullplot.lespe

#######################################

#Lipid extract  
  
li <- read.csv("lipiddatanormalised2.csv", header=T)
li$Sample <- as.factor(li$Sample)
li$Sample <- as.character.factor(li$Sample)

li.pca <- prcomp(li[ ,c(5:163)])
summary(li.pca)
li.pca.axes <- li.pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var="Sample")
li.pca.data <- full_join(variables, li.pca.axes)

li.pca.data.mod <- li.pca.data %>%
  unite('Combo', Diversity:Treatment, remove=F) %>%
  mutate(SpDiv = paste0(Species, "_", Diversity))

li.fullplot.mod <- ggplot(li.pca.data.mod, aes(x=PC1, y=PC2, colour=Species, shape=Combo, 
                                               fill=Species)) +
  geom_point(size=2.5, alpha=0.7, stroke=1.0) + 
  scale_colour_manual(values=c("#ff8b00ff", "#3d3d3dff")) + 
  scale_fill_manual(values=c("#ff8b00ff", "#3d3d3dff")) +
  scale_shape_manual(values=c(1, 2, 21, 24)) +
  xlab("Lipid PC1 (67.8%)") + ylab("Lipid PC2 (13.7%)") + 
  theme_classic() +
  theme(legend.position = "top") 
li.fullplot.mod

fullplotmod <- ggarrange(aq.fullplot.mod, li.fullplot.mod, ncol=2, labels=c("a", "b"))
ggsave("fullplotmod.svg", width=180, height=100, units="mm")

#Plot these PCA coordinates but for each species individually (for lipid)

li.pca.data.andro <- li.pca.data %>% 
  filter(Species=="Andropogon gerardii") %>%
  unite('Combo', Diversity:Treatment, remove=F)

li.fullplot.andro <- ggplot(li.pca.data.andro, aes(x=PC1, y=PC2, colour=Combo, shape=Diversity)) +
  geom_point(size=2.5, alpha=0.7) + 
  scale_colour_manual(values=c("#ffb100ff", "#ff8b00ff", "#838383ff", "#3d3d3dff")) + 
  xlab("Lipid PC1") + ylab("Lipid PC2") + 
  theme_classic() +
  theme(legend.position = "top") 
li.fullplot.andro

li.pca.data.lespe <- li.pca.data %>% 
  filter(Species=="Lespedeza capitata") %>%
  unite('Combo', Diversity:Treatment, remove=F)

li.fullplot.lespe <- ggplot(li.pca.data.lespe, aes(x=PC1, y=PC2, colour=Combo, shape=Diversity)) +
  geom_point(size=2.5, alpha=0.7) + 
  scale_colour_manual(values=c("#ffb100ff", "#ff8b00ff", "#838383ff", "#3d3d3dff")) + 
  xlab("Lipid PC1") + ylab("Lipid PC2") + 
  theme_classic() +
  theme(legend.position = "top") 
li.fullplot.lespe

#NOW SPLIT SPECIES BY SPECIES

#Andropogon gerardi

#Aqueous extract

androaq <- aq %>% 
  filter(Species=="Andropogon gerardii")

variablesandro <- androaq %>%
  select(Species, Diversity, Treatment) %>%
  rownames_to_column(var="Sample")

androaq.pca <- prcomp(androaq[ ,c(5:174)])
summary(androaq.pca)
androaq.pca.axes <- androaq.pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var="Sample")
androaq.pca.data <- full_join(variablesandro, androaq.pca.axes) %>%
   unite('Combo', Diversity:Treatment, remove=F)

androaq.plot <- ggplot(androaq.pca.data, aes(x=PC1, y=PC2, color=Species, shape=Combo, 
                                             fill=Species)) +
  geom_point(size=2.5, alpha=0.7, stroke=1.0) +
  xlab("Aqueous PC1 (73.3%)") + ylab("Aqueous PC2 (8.4%)") +
  scale_colour_manual(values="#ff8b00ff") + 
  scale_fill_manual(values="#ff8b00ff") +
  scale_shape_manual(values=c(1, 2, 21, 24)) +
  theme_classic() +
  theme(legend.position = "none") + 
  xlim(-3500, 4500) + ylim(-2500, 1000)
androaq.plot

#Lipid extract

androli <- li %>% 
  filter(Species=="Andropogon gerardii")

androli.pca <- prcomp(androli[ ,c(5:163)])
summary(androli.pca)
androli.pca.axes <- androli.pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var="Sample")
androli.pca.data <- full_join(variablesandro, androli.pca.axes) %>%
  unite('Combo', Diversity:Treatment, remove=F)

androli.plot <- ggplot(androli.pca.data, aes(x=PC1, y=PC2, color=Species, shape=Combo, 
                                             fill=Species)) +
  geom_point(size=2.5, alpha=0.7, stroke=1.0) +
  xlab("Lipid PC1 (52.5%)") + ylab("Lipid PC2 (28.4%)") +
  scale_colour_manual(values="#ff8b00ff") + 
  scale_fill_manual(values="#ff8b00ff") +
  scale_shape_manual(values=c(1, 2, 21, 24)) +
  theme_classic() +
  theme(legend.position = "none") +
  xlim(-450, 350) + ylim(-250, 250)
androli.plot

fullplotandro <- ggarrange(androaq.plot, androli.plot, ncol=2, labels=c("a", "b"))
ggsave("fullplotandro.svg", width=180, height=100, units="mm")

#Lespedeza capitata

#Aqueous extract

lespeaq <- aq %>% 
  filter(Species=="Lespedeza capitata")

variableslespe <- lespeaq %>%
  select(Species, Diversity, Treatment) %>%
  rownames_to_column(var="Sample")

lespeaq.pca <- prcomp(lespeaq[ ,c(5:174)])
summary(lespeaq.pca)

lespeaq.pca.axes <- lespeaq.pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var="Sample")
lespeaq.pca.data <- full_join(variableslespe, lespeaq.pca.axes) %>%
  unite('Combo', Diversity:Treatment, remove=F)

lespeaq.plot <- ggplot(lespeaq.pca.data, aes(x=PC1, y=PC2, color=Species, shape=Combo, 
                                             fill=Species)) +
  geom_point(size=2.5, alpha=0.7, stroke=1.0) +
  xlab("Aqueous PC1 (48.5%)") + ylab("Aqueous PC2 (19.1%)") +
  scale_colour_manual(values="#3d3d3dff") + 
  scale_fill_manual(values="#3d3d3dff") +
  scale_shape_manual(values=c(1, 2, 21, 24)) +
  theme_classic() +
  theme(legend.position = "none") +
  xlim(-1800, 2000) + ylim(-1200, 800)
lespeaq.plot

#Lipid extract

lespeli <- li %>% 
  filter(Species=="Lespedeza capitata")

lespeli.pca <- prcomp(lespeli[ ,c(5:163)])
summary(lespeli.pca)
lespeli.pca.axes <- lespeli.pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var="Sample")
lespeli.pca.data <- full_join(variableslespe, lespeli.pca.axes) %>%
  unite('Combo', Diversity:Treatment, remove=F)

lespeli.plot <- ggplot(lespeli.pca.data, aes(x=PC1, y=PC2, color=Species, shape=Combo, 
                                             fill=Species)) +
  geom_point(size=2.5, alpha=0.7, stroke=1.0) +
  xlab("Lipid PC1 (45.4%)") + ylab("Lipid PC2 (21.8%)") +
  scale_colour_manual(values="#3d3d3dff") + 
  scale_fill_manual(values="#3d3d3dff") +
  scale_shape_manual(values=c(1, 2, 21, 24)) +
  theme_classic() +
  theme(legend.position = "none") +
  xlim(-250, 300) + ylim(-250, 200)
lespeli.plot

fullplotlespe <- ggarrange(lespeaq.plot, lespeli.plot, ncol=2, labels=c("a", "b"))
ggsave("fullplotlespe.svg", width=180, height=100, units="mm")

#All four within-species PCAs together

allpca <- ggarrange(lespeaq.plot, lespeli.plot, androaq.plot, androli.plot, 
                    ncol=2, nrow=2, labels=c("a", "b", "c", "d"))
ggsave("allpca.svg", width=180, height=180, units="mm")

###########################################################################

#Try do some statistical tests using PERMANOVA

#Make new variable with combo of both treatments

androaq <- androaq %>%
  unite('Combo', Diversity:Treatment, remove=F)
androli <- androli %>%
  unite('Combo', Diversity:Treatment, remove=F)
lespeaq <- lespeaq %>%
  unite('Combo', Diversity:Treatment, remove=F)
lespeli <- lespeli %>%
  unite('Combo', Diversity:Treatment, remove=F)

#Compare species

adonis2(aq[ ,c(5:174)] ~ aq$Species, 
        permutations=10000, method="euclidean", by="margin")

adonis2(li[ ,c(5:163)] ~  li$Species, 
        permutations=10000, method="euclidean", by="margin")

#For the below, interactions were tested (i.e. Diversity*Treatment), and found
#to be non-significant

#Andropogon gerardii aqueous (androaq) 

adonis2(androaq[ ,c(6:175)] ~ androaq$Diversity + androaq$Treatment, 
       permutations=10000, method="euclidean", by="margin")

#Andropogon gerardii lipid (androli)

adonis2(androli[ ,c(6:164)] ~  androli$Diversity + androli$Treatment, 
       permutations=10000, method="euclidean", by="margin")

#Lespedeza capitata aqueous (lespeaq)

adonis2(lespeaq[ ,c(6:175)] ~ lespeaq$Diversity + lespeaq$Treatment, 
       permutations=10000,  method="euclidean", by="margin")

#Lespedeza capitata lipid (lespeli)

adonis2(lespeli[ ,c(6:164)] ~ lespeli$Diversity + lespeli$Treatment, 
       permutations=10000,  method="euclidean", by="margin")

#Try testing each combo as the grouping

#Andropogon gerardii aqueous (androaq)

adonis2(androaq[ ,c(6:175)] ~ androaq$Combo, 
        permutations=10000, method="euclidean")

#Andropogon gerardii lipid (androli)

adonis2(androli[ ,c(6:164)] ~  androli$Combo, 
        permutations=10000, method="euclidean")

#Lespedeza capitata aqueous (lespeaq)

adonis2(lespeaq[ ,c(6:175)] ~ lespeaq$Combo, 
        permutations=10000,  method="euclidean")

#Lespedeza capitata lipid (lespeli)

adonis2(lespeli[ ,c(6:164)] ~ lespeli$Combo, 
        permutations=10000,  method="euclidean")

#Try now doing ANOSIM (using rank order of dissimilarities)

anosim(aq[ ,c(5:174)], aq$Species, 
        permutations=10000, distance="euclidean")

anosim(li[ ,c(5:163)], li$Species, 
        permutations=10000, distance="euclidean")

#Andropogon gerardii aqueous (androaq)

anosim(androaq[ ,c(6:175)], androaq$Diversity, 
        permutations=10000, distance="euclidean")

anosim(androaq[ ,c(6:175)], androaq$Treatment, 
       permutations=10000, distance="euclidean")

#Andropogon gerardii lipid (androli)

anosim(androli[ ,c(6:164)], androli$Diversity, 
        permutations=10000, distance="euclidean")

anosim(androli[ ,c(6:164)], androli$Treatment, 
       permutations=10000, distance="euclidean")

#Lespedeza capitata aqueous (lespeaq)

anosim(lespeaq[ ,c(6:175)], lespeaq$Diversity, 
        permutations=10000,  distance="euclidean")

anosim(lespeaq[ ,c(6:175)], lespeaq$Treatment, 
       permutations=10000,  distance="euclidean")

#Lespedeza capitata lipid (lespeli)

anosim(lespeli[ ,c(6:164)], lespeli$Diversity, 
        permutations=10000,  distance="euclidean")

anosim(lespeli[ ,c(6:164)], lespeli$Treatment, 
       permutations=10000,  distance="euclidean")

#Try testing each combo as the grouping
#Nothing of significance

#Andropogon gerardii aqueous (androaq)

anosim(androaq[ ,c(6:175)], androaq$Combo, 
        permutations=10000, distance="euclidean")

#Andropogon gerardii lipid (androli)

anosim(androli[ ,c(6:164)], androli$Combo, 
        permutations=10000, distance="euclidean")

#Lespedeza capitata aqueous (lespeaq)

anosim(lespeaq[ ,c(6:175)], lespeaq$Combo, 
        permutations=10000,  distance="euclidean")

#Lespedeza capitata lipid (lespeli)

anosim(lespeli[ ,c(6:164)], lespeli$Combo, 
        permutations=10000,  distance="euclidean")

#See if there are differences in dispersion (multivariate Levene's)
#Nothing significant re: dispersion

#Andropogon aqueous

anova(betadisper(vegdist(androaq[ ,c(6:175)], method="euclidean"), androaq$Diversity))
anova(betadisper(vegdist(androaq[ ,c(6:175)], method="euclidean"), androaq$Treatment))
anova(betadisper(vegdist(androaq[ ,c(6:175)], method="euclidean"), androaq$Combo))

#Andropogon lipid

anova(betadisper(vegdist(androli[ ,c(6:164)], method="euclidean"), androli$Diversity))
anova(betadisper(vegdist(androli[ ,c(6:164)], method="euclidean"), androli$Treatment))
anova(betadisper(vegdist(androli[ ,c(6:164)], method="euclidean"), androli$Combo))

#Lespedeza aqueous

anova(betadisper(vegdist(lespeaq[ ,c(6:175)], method="euclidean"), lespeaq$Diversity))
anova(betadisper(vegdist(lespeaq[ ,c(6:175)], method="euclidean"), lespeaq$Treatment))
anova(betadisper(vegdist(lespeaq[ ,c(6:175)], method="euclidean"), lespeaq$Combo))

#Lespedeza lipid

anova(betadisper(vegdist(lespeli[ ,c(6:164)], method="euclidean"), lespeli$Diversity))
anova(betadisper(vegdist(lespeli[ ,c(6:164)], method="euclidean"), lespeli$Treatment))
anova(betadisper(vegdist(lespeli[ ,c(6:164)], method="euclidean"), lespeli$Combo))

anova(betadisper(vegdist(aq[ ,c(5:174)], method="euclidean"), aq$Species))

anova(betadisper(vegdist(li[ ,c(5:163)], method="euclidean"), li$Species))

#Additional statistical tests: diversity and treatment against individual PCs

#Andropogon aqueous

mod1 <- lm(PC1 ~ Diversity, data=androaq.pca.data)
anova(mod1)

mod2 <- lm(PC1 ~ Treatment, data=androaq.pca.data)
anova(mod2)

mod2a <- lm(PC1 ~ Treatment*Diversity, data=androaq.pca.data)
anova(mod2a)

mod3 <- lm(PC2 ~ Diversity, data=androaq.pca.data)
anova(mod3)

mod4 <- lm(PC2 ~ Treatment, data=androaq.pca.data)
anova(mod4)

mod4a <- lm(PC2 ~ Treatment*Diversity, data=androaq.pca.data)
anova(mod4a)

mod5 <- lm(PC3 ~ Diversity, data=androaq.pca.data)
anova(mod5)

mod6 <- lm(PC3 ~ Treatment, data=androaq.pca.data)
anova(mod6)

mod6a <- lm(PC3 ~ Treatment*Diversity, data=androaq.pca.data)
anova(mod6a)

#Andropogon lipid

mod7 <- lm(PC1 ~ Diversity, data=androli.pca.data)
anova(mod7)

mod8 <- lm(PC1 ~ Treatment, data=androli.pca.data)
anova(mod8)

mod8a <- lm(PC1 ~ Treatment*Diversity, data=androli.pca.data)
anova(mod8a)

mod9 <- lm(PC2 ~ Diversity, data=androli.pca.data)
anova(mod9)

mod10 <- lm(PC2 ~ Treatment, data=androli.pca.data)
anova(mod10)

mod10a <- lm(PC2 ~ Treatment*Diversity, data=androli.pca.data)
anova(mod10a)

mod11 <- lm(PC3 ~ Diversity, data=androli.pca.data)
anova(mod11)

mod12 <- lm(PC3 ~ Treatment, data=androli.pca.data)
anova(mod12)

mod12a <- lm(PC3 ~ Treatment*Diversity, data=androli.pca.data)
anova(mod12a)

#Lespedeza aqueous

mod13 <- lm(PC1 ~ Diversity, data=lespeaq.pca.data)
anova(mod13)

mod14 <- lm(PC1 ~ Treatment, data=lespeaq.pca.data)
anova(mod14)

mod14a <- lm(PC1 ~ Treatment*Diversity, data=lespeaq.pca.data)
anova(mod14a)

mod15 <- lm(PC2 ~ Diversity, data=lespeaq.pca.data)
anova(mod15)

mod16 <- lm(PC2 ~ Treatment, data=lespeaq.pca.data)
anova(mod16)

mod16a <- lm(PC2 ~ Treatment*Diversity, data=lespeaq.pca.data)
anova(mod16a)

mod17 <- lm(PC3 ~ Diversity, data=lespeaq.pca.data)
anova(mod17)

mod18 <- lm(PC3 ~ Treatment, data=lespeaq.pca.data)
anova(mod18)

mod18a <- lm(PC3 ~ Treatment*Diversity, data=lespeaq.pca.data)
anova(mod18a)

#Lespedeza lipid

mod19 <- lm(PC1 ~ Diversity, data=lespeli.pca.data)
anova(mod19)

mod20 <- lm(PC1 ~ Treatment, data=lespeli.pca.data)
anova(mod20)

mod20a <- lm(PC1 ~ Treatment*Diversity, data=lespeli.pca.data)
anova(mod20a)

mod21 <- lm(PC2 ~ Diversity, data=lespeli.pca.data)
anova(mod21)

mod22 <- lm(PC2 ~ Treatment, data=lespeli.pca.data)
anova(mod22)

mod22a <- lm(PC2 ~ Treatment*Diversity, data=lespeli.pca.data)
anova(mod22a)

mod23 <- lm(PC3 ~ Diversity, data=lespeli.pca.data)
anova(mod23)

mod24 <- lm(PC3 ~ Treatment, data=lespeli.pca.data)
anova(mod24)

mod24a <- lm(PC3 ~ Treatment*Diversity, data=lespeli.pca.data)
anova(mod24a)
#interaction here, but not investigated further

#Test for correlations among PCs
#Use the overall data (because for the within-species data, 
#PC1 and PC2 will be totally orthogonal so have 0 correlation)

lespeaqcorr <- aq.pca.data %>% 
  filter(Species=="Lespedeza capitata")

androaqcorr <- aq.pca.data %>% 
  filter(Species=="Andropogon gerardii")

lespelicorr <- li.pca.data %>% 
  filter(Species=="Lespedeza capitata")

androlicorr <- li.pca.data %>% 
  filter(Species=="Andropogon gerardii")

mod25 <- lm(PC1 ~ PC2, data=lespeaqcorr)
anova(mod25)

mod26 <- lm(PC1 ~ PC2, data=androaqcorr)
anova(mod26)

mod27 <- lm(PC1 ~ PC2, data=lespelicorr)
anova(mod27)

mod28 <- lm(PC1 ~ PC2, data=androlicorr)
anova(mod28)

#Test whether the species are different from each other in the overall PCA

mod29 <- lm(PC1 ~ Species, data=aq.pca.data)
anova(mod29)

mod30 <- lm(PC2 ~ Species, data=aq.pca.data)
anova(mod30)

mod31 <- lm(PC1 ~ Species, data=li.pca.data)
anova(mod31)

mod32 <- lm(PC2 ~ Species, data=li.pca.data)
anova(mod32)

###########################################################################

#Try to compare clusterings from the aqueous and lipid data for each species
#(do samples cluster in same way across the two extracts?)

#Need to compare: lespeaq with lespeli; androaq with androli
#lespeaq[ ,c(6:175)]
#lespeli[ ,c(6:164)]
#androaq[ ,c(6:175)]
#androli[ ,c(6:164)]

lespeaqdend <- column_to_rownames(lespeaq, var="Sample")
lespeaqdend <- lespeaqdend[ ,c(5:174)] %>%
  vegdist(method = "euclidean") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram()
plot(lespeaqdend)

lespelidend <- column_to_rownames(lespeli, var="Sample")
lespelidend <- lespelidend[ ,c(5:163)] %>%
  vegdist(method = "euclidean") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram()
plot(lespelidend)

lespecor <- cor_bakers_gamma(lespeaqdend, lespelidend)
lespecor

androaqdend <- column_to_rownames(androaq, var="Sample")
androaqdend <- androaqdend[ ,c(5:174)] %>%
  vegdist(method = "euclidean") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram()
plot(androaqdend)

androlidend <- column_to_rownames(androli, var="Sample")
androlidend <- androlidend[ ,c(5:163)] %>%
  vegdist(method = "euclidean") %>%
  hclust(method = "ward.D2") %>%
  as.dendrogram()
plot(androlidend)

androcor <- cor_bakers_gamma(androaqdend, androlidend)
androcor

#Now do permutation test

#For lespe

R <- 1000
cor_bakers_gamma_results_lespe <- numeric(R)
dend_mixed <- lespeaqdend
for(i in 1:R) {
  dend_mixed <- sample.dendrogram(dend_mixed, replace = FALSE)
  cor_bakers_gamma_results_lespe[i] <- cor_bakers_gamma(lespeaqdend, dend_mixed)
}
plot(density(cor_bakers_gamma_results_lespe),
     main = "Baker's gamma distribution under H0 (Lespe)",
     xlim = c(-1,1))
abline(v = 0, lty = 2)
abline(v = lespecor, lty = 2, col = 2)

sum(lespecor < cor_bakers_gamma_results_lespe)/ R

#For andro

cor_bakers_gamma_results_andro <- numeric(R)
dend_mixed <- androaqdend
for(i in 1:R) {
  dend_mixed <- sample.dendrogram(dend_mixed, replace = FALSE)
  cor_bakers_gamma_results_andro[i] <- cor_bakers_gamma(androaqdend, dend_mixed)
}
plot(density(cor_bakers_gamma_results_andro),
     main = "Baker's gamma distribution under H0 (Andro)",
     xlim = c(-1,1))
abline(v = 0, lty = 2)
abline(v = androcor, lty = 2, col = 2)

sum(androcor < cor_bakers_gamma_results_andro)/ R


###########################################################################

#Check whether phytochemical diversity (a la Richards et al. 2015) relates to anything at all

diversity <- read.csv("diversity.csv", header=T, stringsAsFactors = T)

#Plot aqueous diversity 

aqdivplot <- ggplot(diversity, aes(x=Species, y=diversityaq, 
                                   colour=Species, shape=Combo, fill=Species)) +
  geom_point(size=2.5, alpha=0.7, stroke=1.0, position = position_dodge(0.5)) +
  theme_bw() +
  theme(legend.position = "top") +
  scale_colour_manual(values=c("#ff8b00ff", "#3d3d3dff")) + 
  scale_fill_manual(values=c("#ff8b00ff", "#3d3d3dff")) +
  scale_shape_manual(values=c(1, 2, 21, 24)) +
  ylab("Aqueous phytochemical diversity") + 
  theme(legend.position = "none")
aqdivplot

#Plot lipid diversity

lidivplot <- ggplot(diversity, aes(x=Species, y=diversityli, 
                                   colour=Species, shape=Combo, fill=Species)) +
  geom_point(size=2.5, alpha=0.7, stroke=1.0, position = position_dodge(0.5)) +
  theme_bw() +
  theme(legend.position = "top") +
  scale_colour_manual(values=c("#ff8b00ff", "#3d3d3dff")) + 
  scale_fill_manual(values=c("#ff8b00ff", "#3d3d3dff")) +
  scale_shape_manual(values=c(1, 2, 21, 24)) +
  ylab("Lipid phytochemical diversity") + 
  theme(legend.position = "none")
lidivplot

fulldivplot <- ggarrange(aqdivplot, lidivplot, ncol=2, labels=c("a", "b"))
ggsave("fulldivplot.svg", width=180, height=100, units="mm")

#Test for differences among species

mod1 <- lm(diversityaq ~ Species, data=diversity)
summary(mod1)
Anova(mod1)

mod2 <- lm(diversityli ~ Species, data=diversity)
summary(mod2)
Anova(mod2)

#Test for differences among treatments within species

diversitylespe <- diversity %>%
  filter(Species=="Lespedeza capitata")
diversityandro <- diversity %>%
  filter(Species=="Andropogon gerardii")

#Lespedeza capitata

mod3 <- lm(diversityaq ~ Diversity + Treatment, data=diversitylespe)
summary(mod3)
Anova(mod3)

mod3a <- lm(diversityaq ~ Diversity*Treatment, data=diversitylespe)
summary(mod3a)
Anova(mod3a)

mod4 <- lm(diversityli ~ Diversity + Treatment, data=diversitylespe)
summary(mod4)
Anova(mod4)

mod4a <- lm(diversityli ~ Diversity*Treatment, data=diversitylespe)
summary(mod4a)
Anova(mod4a)

#Andropogon gerardii

mod5 <- lm(diversityaq ~ Diversity + Treatment, data=diversityandro)
summary(mod5)
Anova(mod5)

mod5a <- lm(diversityaq ~ Diversity*Treatment, data=diversityandro)
summary(mod5a)
Anova(mod5a)

mod6 <- lm(diversityli ~ Diversity + Treatment, data=diversityandro)
summary(mod6)
Anova(mod6)

mod6a <- lm(diversityli ~ Diversity*Treatment, data=diversityandro)
summary(mod6a)
Anova(mod6a)

###########################################################################

#Peak module analysis

#based on Richards et al. 2018

library(WGCNA)

#for complete data

mixpip <- read.csv("aqueousdata.csv", header=T, row.names=1)

powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(mixpip, powerVector = powers )

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
#based on this, power should be 10

#calculate modules

netpip= blockwiseModules(mixpip, power =10,
                         TOMType = "unsigned", minModuleSize =3,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "mixoutTOM.csv")
#list of modules and number of nodes in them:
table(netpip$colors)
## 0  1  2  3  4  5  6 
## 4 84 48 17  7  6  4

#plot dendrogram of modules
mergedColors = labels2colors(netpip$colors)
plotDendroAndColors(netpip$dendrograms[[1]], mergedColors[netpip$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels = netpip$colors
moduleColors = labels2colors(netpip$colors)
MEs = netpip$MEs;
geneTree = netpip$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "networkConstruction-auto.RData")
lnames = load(file = "networkConstruction-auto.RData");
MEs0 = moduleEigengenes(mixpip, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
lnames
nGenes = ncol(mixpip);
nSamples = nrow(mixpip);

sizeGrWindow(6,6)
par(cex=1.0)
plotEigengeneNetworks(MEs, "Module dendrogram", marDendro = c(0,4,1,0), plotHeatmaps=FALSE)
plotEigengeneNetworks(MEs, "Module adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle=90)

write.csv(MEs,"MEall.csv") # gives a file of modules and eigenvalues across samples, samples are just numbered
write.csv(moduleColors,"modulesall.csv")# gives a file of list the chemical shift (node) for each module

#just for Lespedeza

mixpip <- read.csv("aqueousdatalespe.csv", header=T, row.names=1)

powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(mixpip, powerVector = powers )

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
#based on this, power should be 6

#calculate modules

netpip= blockwiseModules(mixpip, power =6,
                         TOMType = "unsigned", minModuleSize =3,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "mixoutTOM.csv")
#list of modules and number of nodes in them:
table(netpip$colors)
## 1  2  3  4  5  6  7  8  9 10 11 12 
## 67 27 16 13  8  8  6  6  6  5  5  3

#plot dendrogram of modules
mergedColors = labels2colors(netpip$colors)
plotDendroAndColors(netpip$dendrograms[[1]], mergedColors[netpip$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels = netpip$colors
moduleColors = labels2colors(netpip$colors)
MEs = netpip$MEs;
geneTree = netpip$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "networkConstruction-auto.RData")
lnames = load(file = "networkConstruction-auto.RData");
MEs0 = moduleEigengenes(mixpip, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
lnames
nGenes = ncol(mixpip);
nSamples = nrow(mixpip);

sizeGrWindow(6,6)
par(cex=1.0)
plotEigengeneNetworks(MEs, "Module dendrogram", marDendro = c(0,4,1,0), plotHeatmaps=FALSE)
plotEigengeneNetworks(MEs, "Module adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle=90)

write.csv(MEs,"MElespe.csv") # gives a file of modules and eigenvalues across samples, samples are just numbered
write.csv(moduleColors,"moduleslespe.csv")# gives a file of list the chemical shift (node) for each module

#just for Andropogon

mixpip <- read.csv("aqueousdataandro.csv", header=T, row.names=1)

powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(mixpip, powerVector = powers )

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, col="red")
#based on this, power should be 16

#calculate modules

netpip= blockwiseModules(mixpip, power =16,
                         TOMType = "unsigned", minModuleSize =3,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "mixoutTOM.csv")
#list of modules and number of nodes in them:
table(netpip$colors)
## 0   1   2   3   4   5   6   7   8   9 
##13 101  22   9   6   5   4   4   3   3

#plot dendrogram of modules
mergedColors = labels2colors(netpip$colors)
plotDendroAndColors(netpip$dendrograms[[1]], mergedColors[netpip$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels = netpip$colors
moduleColors = labels2colors(netpip$colors)
MEs = netpip$MEs;
geneTree = netpip$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "networkConstruction-auto.RData")
lnames = load(file = "networkConstruction-auto.RData");
MEs0 = moduleEigengenes(mixpip, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
lnames
nGenes = ncol(mixpip);
nSamples = nrow(mixpip);

sizeGrWindow(6,6)
par(cex=1.0)
plotEigengeneNetworks(MEs, "Module dendrogram", marDendro = c(0,4,1,0), plotHeatmaps=FALSE)
plotEigengeneNetworks(MEs, "Module adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle=90)

write.csv(MEs,"MEandro.csv") # gives a file of modules and eigenvalues across samples, samples are just numbered
write.csv(moduleColors,"modulesandro.csv")# gives a file of list the chemical shift (node) for each module

#Plot PCA loadings against chemical shifts

#Andro aqueous

androaqloadings <- androaq.pca[["rotation"]] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Chemical_shift") %>%
  mutate(Chemical_shift = as.numeric(gsub("X", "", Chemical_shift)))

modulesandro <- read.csv("modulesandro_joshedits.csv", header=T) #a cleaned-up version of above file
#join to loadings data 

androaqloadings <- left_join(androaqloadings, modulesandro, by="Chemical_shift")

androaqloadingspc1 <- ggplot(androaqloadings, aes(x=Chemical_shift, y=PC1, 
                                                  fill=module)) +
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("\u03b4 (ppm)") + ylab(expression(paste(italic("A. gerardii") , ' PC1 loading (aq)'))) +
  scale_fill_manual(values=c("black", "blue", "brown", "green","grey", "magenta",
                             "pink", "red", "turquoise", "yellow")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(legend.position="none")
androaqloadingspc1

androaqloadingspc2 <- ggplot(androaqloadings, aes(x=Chemical_shift, y=PC2, 
                                                  fill=module)) +
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("\u03b4 (ppm)") + ylab(expression(paste(italic("A. gerardii") , ' PC2 loading (aq)'))) +
  scale_fill_manual(values=c("black", "blue", "brown", "green","grey", "magenta",
                             "pink", "red", "turquoise", "yellow")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(legend.position="none")
androaqloadingspc2

androaqloadingsplot <- ggarrange(androaqloadingspc1, androaqloadingspc2,
                                 nrow=2, labels=c("a", "b"))
ggsave("androaqloadingsplot.svg", width=180, height=160, units="mm")

#Andro lipid

androliloadings <- androli.pca[["rotation"]] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Chemical_shift") %>%
  mutate(Chemical_shift = as.numeric(gsub("X", "", Chemical_shift)))

androliloadingspc1 <- ggplot(androliloadings, aes(x=Chemical_shift, y=PC1)) +
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("\u03b4 (ppm)") + ylab(expression(paste(italic("A. gerardii") , ' PC1 loading (li)'))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
androliloadingspc1

androliloadingspc2 <- ggplot(androliloadings, aes(x=Chemical_shift, y=PC2)) +
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("\u03b4 (ppm)") + ylab(expression(paste(italic("A. gerardii") , ' PC2 loading (li)'))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
androliloadingspc2

androliloadingsplot <- ggarrange(androliloadingspc1, androliloadingspc2,
                                 nrow=2, labels=c("a", "b"))
ggsave("androliloadingsplot.svg", width=180, height=160, units="mm")

#Lespe aqueous

lespeaqloadings <- lespeaq.pca[["rotation"]] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Chemical_shift") %>%
  mutate(Chemical_shift = as.numeric(gsub("X", "", Chemical_shift)))

moduleslespe <- read.csv("moduleslespe_joshedits.csv", header=T) #a cleaned-up version of above file
#join to loadings data 

lespeaqloadings <- left_join(lespeaqloadings, moduleslespe, by="Chemical_shift")

lespeaqloadingspc1 <- ggplot(lespeaqloadings, aes(x=Chemical_shift, y=PC1,
                                                  fill=module)) +
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("\u03b4 (ppm)") + ylab(expression(paste(italic("L. capitata") , ' PC1 loading (aq)'))) +
  scale_fill_manual(values=c("black", "blue", "brown", "green","greenyellow", "magenta",
                             "pink", "purple", "red","tan", "turquoise", "yellow")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(legend.position = "none")
lespeaqloadingspc1

lespeaqloadingspc2 <- ggplot(lespeaqloadings, aes(x=Chemical_shift, y=PC2,
                                                  fill=module)) +
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("\u03b4 (ppm)") + ylab(expression(paste(italic("L. capitata") , ' PC2 loading (aq)'))) +
  scale_fill_manual(values=c("black", "blue", "brown", "green","greenyellow", "magenta",
                             "pink", "purple", "red","tan", "turquoise", "yellow")) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(legend.position = "none")
lespeaqloadingspc2

lespeaqloadingsplot <- ggarrange(lespeaqloadingspc1, lespeaqloadingspc2,
                                 nrow=2, labels=c("a", "b"))
ggsave("lespeaqloadingsplot.svg", width=180, height=160, units="mm")

#Lespe lipid

lespeliloadings <- lespeli.pca[["rotation"]] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Chemical_shift") %>%
  mutate(Chemical_shift = as.numeric(gsub("X", "", Chemical_shift)))

lespeliloadingspc1 <- ggplot(lespeliloadings, aes(x=Chemical_shift, y=PC1)) +
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("\u03b4 (ppm)") + ylab(expression(paste(italic("L. capitata") , ' PC1 loading (li)'))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
lespeliloadingspc1

lespeliloadingspc2 <- ggplot(lespeliloadings, aes(x=Chemical_shift, y=PC2)) +
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("\u03b4 (ppm)") + ylab(expression(paste(italic("L. capitata") , ' PC2 loading (li)'))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
lespeliloadingspc2

lespeliloadingsplot <- ggarrange(lespeliloadingspc1, lespeliloadingspc2,
                                 nrow=2, labels=c("a", "b"))
ggsave("lespeliloadingsplot.svg", width=180, height=160, units="mm")

#Andro vs. Lespe aqueous

fullaqloadings <- aq.pca[["rotation"]] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Chemical_shift") %>%
  mutate(Chemical_shift = as.numeric(gsub("X", "", Chemical_shift)))

#read in module information

modules <- read.csv("modulesall_joshedits.csv", header=T)
#join to loadings data 

fullaqloadings <- left_join(fullaqloadings, modules, by="Chemical_shift")

aqloadingspc1 <- ggplot(fullaqloadings, aes(x=Chemical_shift, y=PC1, fill=module)) +
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("\u03b4 (ppm)") + ylab("PC1 loading (aq)") +
  scale_fill_manual(values=c("blue", "brown", "green", "grey", "red", 
                             "turquoise", "yellow")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(legend.position = "none")
aqloadingspc1

aqloadingspc2 <- ggplot(fullaqloadings, aes(x=Chemical_shift, y=PC2, fill=module)) +
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("\u03b4 (ppm)") + ylab("PC2 loading (aq)") +
  scale_fill_manual(values=c("blue", "brown", "green", "grey", "red", 
                             "turquoise", "yellow")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(legend.position = "none")
aqloadingspc2

aqloadingsplot <- ggarrange(aqloadingspc1, aqloadingspc2,
                                 nrow=2, labels=c("a", "b"))
ggsave("aqloadingsplot.svg", width=180, height=160, units="mm")

#Andro vs. Lespe lipid

fullliloadings <- li.pca[["rotation"]] %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Chemical_shift") %>%
  mutate(Chemical_shift = as.numeric(gsub("X", "", Chemical_shift)))

liloadingspc1 <- ggplot(fullliloadings, aes(x=Chemical_shift, y=PC1)) +
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("\u03b4 (ppm)") + ylab("PC1 loading (li)") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
liloadingspc1

liloadingspc2 <- ggplot(fullliloadings, aes(x=Chemical_shift, y=PC2)) +
  geom_bar(stat="identity") +
  theme_bw() +
  xlab("\u03b4 (ppm)") + ylab("PC2 loading (li)") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
liloadingspc2

liloadingsplot <- ggarrange(liloadingspc1, liloadingspc2,
                                 nrow=2, labels=c("a", "b"))
ggsave("liloadingsplot.svg", width=180, height=160, units="mm")

#Try and plot the modules/see how they correlate to treatments
#use geom_tile

lespemodules <- read.csv("MElespe_joshedits.csv", header=T, stringsAsFactors = T)
lespemodules$Module <- as.factor(lespemodules$Module)
andromodules <- read.csv("MEandro_joshedits.csv", header=T, stringsAsFactors = T)
andromodules$Module <- as.factor(andromodules$Module)

#Lespedeza capitata

lcmodules <- ggplot(lespemodules, aes(x=Module, y=Treatment)) +
  geom_raster(aes(fill = Eigenvalue)) +
  scale_fill_gradient2(low="red", mid="white", high="blue")
lcmodules
#modules 1, 4, 5, 6, 9 all re: diversity

#Andropogon gerardii

agmodules <- ggplot(andromodules, aes(x=Module, y=Treatment)) +
  geom_raster(aes(fill = Eigenvalue)) +
  scale_fill_gradient2(low="red", mid="white", high="blue")
agmodules
#maybe module 5 re: pesticides

modulesplot <- ggarrange(lcmodules, agmodules,
                            nrow=2, labels=c("a Lespedeza capitata",
                                             "b Andropogon gerardii"))
ggsave("modulesplot.svg", width=180, height=160, units="mm")

#Test for differences in modules with treatment - Lespe

lespemodule1 <- filter(lespemodules, Module==1) %>%
  select(Eigenvalue) %>%
  rename(Module1=Eigenvalue)
lespemodule2 <- filter(lespemodules, Module==2) %>%
  select(Eigenvalue) %>%
  rename(Module2=Eigenvalue)
lespemodule3 <- filter(lespemodules, Module==3) %>%
  select(Eigenvalue) %>%
  rename(Module3=Eigenvalue)
lespemodule4 <- filter(lespemodules, Module==4) %>%
  select(Eigenvalue) %>%
  rename(Module4=Eigenvalue)
lespemodule5 <- filter(lespemodules, Module==5) %>%
  select(Eigenvalue) %>%
  rename(Module5=Eigenvalue)
lespemodule6 <- filter(lespemodules, Module==6) %>%
  select(Eigenvalue) %>%
  rename(Module6=Eigenvalue)
lespemodule7 <- filter(lespemodules, Module==7) %>%
  select(Eigenvalue) %>%
  rename(Module7=Eigenvalue)
lespemodule8 <- filter(lespemodules, Module==8) %>%
  select(Eigenvalue) %>%
  rename(Module8=Eigenvalue)
lespemodule9 <- filter(lespemodules, Module==9) %>%
  select(Eigenvalue) %>%
  rename(Module9=Eigenvalue)
lespemodule10 <- filter(lespemodules, Module==10) %>%
  select(Eigenvalue) %>%
  rename(Module10=Eigenvalue)
lespemodule11 <- filter(lespemodules, Module==11) %>%
  select(Eigenvalue) %>%
  rename(Module11=Eigenvalue)
lespemodule12 <- filter(lespemodules, Module==12) %>%
  select(Eigenvalue) %>%
  rename(Module12=Eigenvalue)

dependent_vars <- as.data.frame(cbind(lespemodule1, lespemodule2, lespemodule3, lespemodule4, lespemodule5, lespemodule6,
                        lespemodule7, lespemodule8, lespemodule9, lespemodule10, lespemodule11, lespemodule12))

modulegeneral <- filter(lespemodules, Module==1)
independent_var <- modulegeneral$Treatment

lespedataframe <- cbind(independent_var, dependent_vars)

diverse <- c("Mono", "Mono", "Mono", "Poly", "Poly", "Poly",  "Mono", "Mono", "Mono", "Poly", "Poly", "Poly")
enemies <- c("Control", "Control", "Control", "Control", "Control", "Control", "Pest", "Pest", "Pest", "Pest", "Pest", "Pest")
lespedataframe <- cbind(lespedataframe, enemies)
lespedataframe <- cbind(lespedataframe, diverse)

summary(lm(Module1 ~ enemies + diverse, data=lespedataframe))
summary(lm(Module2 ~ enemies + diverse, data=lespedataframe))
summary(lm(Module3 ~ enemies + diverse, data=lespedataframe))
summary(lm(Module4 ~ enemies + diverse, data=lespedataframe))
summary(lm(Module5 ~ enemies + diverse, data=lespedataframe))
summary(lm(Module6 ~ enemies + diverse, data=lespedataframe))
summary(lm(Module7 ~ enemies + diverse, data=lespedataframe))
summary(lm(Module8 ~ enemies*diverse, data=lespedataframe))
summary(lm(Module9 ~ enemies + diverse, data=lespedataframe))
summary(lm(Module10 ~ enemies + diverse, data=lespedataframe))
summary(lm(Module11 ~ enemies + diverse, data=lespedataframe))
summary(lm(Module12 ~ enemies + diverse, data=lespedataframe))

#Test for differences in modules with treatment - Andro

andromodule1 <- filter(andromodules, Module==1) %>%
  select(Eigenvalue) %>%
  rename(Module1=Eigenvalue)
andromodule2 <- filter(andromodules, Module==2) %>%
  select(Eigenvalue) %>%
  rename(Module2=Eigenvalue)
andromodule3 <- filter(andromodules, Module==3) %>%
  select(Eigenvalue) %>%
  rename(Module3=Eigenvalue)
andromodule4 <- filter(andromodules, Module==4) %>%
  select(Eigenvalue) %>%
  rename(Module4=Eigenvalue)
andromodule5 <- filter(andromodules, Module==5) %>%
  select(Eigenvalue) %>%
  rename(Module5=Eigenvalue)
andromodule6 <- filter(andromodules, Module==6) %>%
  select(Eigenvalue) %>%
  rename(Module6=Eigenvalue)
andromodule7 <- filter(andromodules, Module==7) %>%
  select(Eigenvalue) %>%
  rename(Module7=Eigenvalue)
andromodule8 <- filter(andromodules, Module==8) %>%
  select(Eigenvalue) %>%
  rename(Module8=Eigenvalue)
andromodule9 <- filter(andromodules, Module==9) %>%
  select(Eigenvalue) %>%
  rename(Module9=Eigenvalue)
andromodule10 <- filter(andromodules, Module==10) %>%
  select(Eigenvalue) %>%
  rename(Module10=Eigenvalue)

dependent_vars <- as.data.frame(cbind(andromodule1, andromodule2, andromodule3, andromodule4, andromodule5, andromodule6,
                                      andromodule7, andromodule8, andromodule9, andromodule10))

androdataframe <- cbind(independent_var, dependent_vars)

diverse <- c("Mono", "Mono", "Mono", "Poly", "Poly", "Poly",  "Mono", "Mono", "Mono", "Poly", "Poly", "Poly")
enemies <- c("Control", "Control", "Control", "Control", "Control", "Control", "Pest", "Pest", "Pest", "Pest", "Pest", "Pest")
androdataframe <- cbind(androdataframe, enemies)
androdataframe <- cbind(androdataframe, diverse)

summary(lm(Module1 ~ enemies + diverse, data=androdataframe))
summary(lm(Module2 ~ enemies + diverse, data=androdataframe))
summary(lm(Module3 ~ enemies + diverse, data=androdataframe))
summary(lm(Module4 ~ enemies + diverse, data=androdataframe))
summary(lm(Module5 ~ enemies + diverse, data=androdataframe))
summary(lm(Module6 ~ enemies + diverse, data=androdataframe))
summary(lm(Module7 ~ enemies + diverse, data=androdataframe))
summary(lm(Module8 ~ enemies + diverse, data=androdataframe))
summary(lm(Module9 ~ enemies + diverse, data=androdataframe))
summary(lm(Module10 ~ enemies + diverse, data=androdataframe))

############################################################

#Redo analyses using only 50% of the characteristics (top 50%)

#Extract loadings, redo with only informative characteristics (just manipulate in excel
#and load improved versions)

aqloadings <- aq.pca[["rotation"]]
write.csv(aqloadings, "rawaqloadings.csv")
liloadings <- li.pca[["rotation"]]
write.csv(liloadings, "rawliloadings.csv")

#Reload in with half the characteristics gone and repeat whole procedure

##OVERALL PCAs ACROSS BOTH SPECIES

#Aqueous extract

aq50 <- read.csv("aqueousdatanormalised250%.csv", header=T)
aq50$Sample <- as.factor(aq50$Sample)
aq50$Sample <- as.character.factor(aq50$Sample)

variables50 <- aq50 %>%
  select(Sample, Species, Diversity, Treatment)

aq50.pca <- prcomp(aq50[ ,c(5:89)])
summary(aq50.pca)
aq50.pca.axes <- aq50.pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var="Sample")
aq50.pca.data <- full_join(variables50, aq50.pca.axes)

aq50.fullplot <- ggplot(aq50.pca.data, aes(x=PC1, y=PC2, colour=Species)) +
  geom_point(size=2.5, alpha=0.7) + 
  stat_ellipse(type="norm") +
  scale_colour_manual(values=c('#009E73', '#AA4499')) + 
  xlab("Aqueous PC1 (67.1%)") + ylab("Aqueous PC2 (18.7%)") + 
  theme_classic() +
  theme(legend.position = "top") 
aq50.fullplot

#Lipid extract  

li50 <- read.csv("lipiddatanormalised250%.csv", header=T)
li50$Sample <- as.factor(li50$Sample)
li50$Sample <- as.character.factor(li50$Sample)

li50.pca <- prcomp(li50[ ,c(5:84)])
summary(li50.pca)
li50.pca.axes <- li50.pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var="Sample")
li50.pca.data <- full_join(variables, li50.pca.axes)

li50.fullplot <- ggplot(li50.pca.data, aes(x=PC1, y=PC2, colour=Species)) +
  geom_point(size=2.5, alpha=0.7) + 
  stat_ellipse(type="norm") +
  scale_colour_manual(values=c('#009E73', '#AA4499')) + 
  xlab("Lipid PC1 (70.4%)") + ylab("Lipid PC2 (14.1%)") +
  theme_classic()  +
  theme(legend.position = "top") 
li50.fullplot

fullplot50 <- ggarrange(aq50.fullplot, li50.fullplot, ncol=2, labels=c("a", "b"))
ggsave("fullplot.svg", width=180, height=100, units="mm")

#NOW SPLIT SPECIES BY SPECIES

#Andropogon gerardii

#Aqueous extract

androaq50 <- aq50 %>% 
  filter(Species=="Andropogon gerardii")

variablesandro50 <- androaq50 %>%
  select(Species, Diversity, Treatment) %>%
  rownames_to_column(var="Sample")

androaq50.pca <- prcomp(androaq50[ ,c(5:89)])
summary(androaq50.pca)
androaq50.pca.axes <- androaq50.pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var="Sample")
androaq50.pca.data <- full_join(variablesandro50, androaq50.pca.axes)

androaq.plot <- ggplot(androaq50.pca.data, aes(x=PC1, y=PC2, color=Treatment, shape=Diversity)) +
  geom_point(size=3, alpha=0.7) +
  xlab("Aqueous PC1 (76.9%)") + ylab("Aqueous PC2 (7.2%)") +
  scale_colour_manual(values=c('#009E73', '#0de319ff')) + 
  theme_classic() +
  theme(legend.position = "top") 
androaq.plot

#Lipid extract

androli50 <- li50 %>% 
  filter(Species=="Andropogon gerardii")

variablesandroli <- androli50 %>%
  select(Species, Diversity, Treatment) %>%
  rownames_to_column(var="Sample")

androli50.pca <- prcomp(androli50[ ,c(5:84)])
summary(androli50.pca)
androli50.pca.axes <- androli50.pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var="Sample")
androli50.pca.data <- full_join(variablesandroli, androli50.pca.axes)

androli50.plot <- ggplot(androli50.pca.data, aes(x=PC1, y=PC2, color=Treatment, shape=Diversity)) +
  geom_point(size=3, alpha=0.7) +
  xlab("Lipid PC1 (54.8%)") + ylab("Lipid PC2 (28.3%)") +
  scale_colour_manual(values=c('#009E73', '#0de319ff')) + 
  theme_classic() +
  theme(legend.position = "top") 
androli50.plot

fullplotandro50 <- ggarrange(androaq50.plot, androli50.plot, ncol=2, labels=c("a", "b"))
ggsave("fullplotandro.svg", width=180, height=100, units="mm")

#Lespedeza capitata

#Aqueous extract

lespeaq50 <- aq50 %>% 
  filter(Species=="Lespedeza capitata")

variableslespe50 <- lespeaq50 %>%
  select(Species, Diversity, Treatment) %>%
  rownames_to_column(var="Sample")

lespeaq50.pca <- prcomp(lespeaq50[ ,c(5:89)])
summary(lespeaq50.pca)
lespeaq50.pca.axes <- lespeaq50.pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var="Sample")
lespeaq50.pca.data <- full_join(variableslespe50, lespeaq50.pca.axes)

lespeaq50.plot <- ggplot(lespeaq50.pca.data, aes(x=PC1, y=PC2, color=Treatment, shape=Diversity)) +
  geom_point(size=3, alpha=0.7) +
  xlab("Aqueous PC1 (51.5%)") + ylab("Aqueous PC2 (21.0%)") +
  scale_colour_manual(values=c('#AA4499', '#6e4bf9ff')) + 
  theme_classic() +
  theme(legend.position = "top") 
lespeaq50.plot

#Lipid extract

lespeli50 <- li50 %>% 
  filter(Species=="Lespedeza capitata")

lespeli50.pca <- prcomp(lespeli50[ ,c(5:84)])
summary(lespeli50.pca)
lespeli50.pca.axes <- lespeli50.pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var="Sample")
lespeli50.pca.data <- full_join(variableslespe50, lespeli50.pca.axes)

lespeli50.plot <- ggplot(lespeli50.pca.data, aes(x=PC1, y=PC2, color=Treatment, shape=Diversity)) +
  geom_point(size=3, alpha=0.7) +
  xlab("Lipid PC1 (51.4%)") + ylab("Lipid PC2 (20.7%)") +
  scale_colour_manual(values=c('#AA4499', '#6e4bf9ff')) + 
  theme_classic() +
  theme(legend.position = "top") 
lespeli50.plot

###########################

#Try do some statistical tests using PERMANOVA

#NO effects (even in models when tested seperately)

#Andropogon gerardii aqueous (androaq50)

adonis2(androaq50[ ,c(5:89)] ~ androaq50$Diversity + androaq50$Treatment, 
        permutations=10000, method="euclidean", by="margin")

#Andropogon gerardii lipid (androli50)

adonis2(androli50[ ,c(5:84)] ~  androli50$Diversity + androli50$Treatment, 
        permutations=10000, method="euclidean", by="margin")

#Lespedeza capitata aqueous (lespeaq50)

adonis2(lespeaq50[ ,c(5:89)] ~ lespeaq50$Diversity + lespeaq50$Treatment, 
        permutations=10000,  method="euclidean", by="margin")

#Lespedeza capitata lipid (lespeli50)

adonis2(lespeli50[ ,c(5:84)] ~ lespeli50$Diversity + lespeli50$Treatment, 
        permutations=10000,  method="euclidean", by="margin")

#Re-do with only phenolics section (>5ppm)

#Exclude one outlier

diversity5ppm <- diversity[-18, ]

#Plot aqueous diversity >5ppm

aqdivplot5ppm <- ggplot(diversity5ppm, aes(x=Species, y=diversityaq5ppm, fill=Combo)) +
  geom_dotplot(binaxis='y', stackdir='center', position = position_dodge(0.5)) +
  theme_bw() +
  theme(legend.position = "top") +
  scale_fill_manual(values=c("#ffb100ff", "#ff8b00ff", "#838383ff", "#3d3d3dff"))+
  ylab("Aqueous phytochemical diversity") 
aqdivplot5ppm

#Test for differences among species

mod7 <- lm(diversityaq5ppm ~ Species, data=diversity5ppm)
summary(mod7)
Anova(mod7)

#Test for differences among treatments within species

diversitylespe5ppm <- diversity5ppm %>%
  filter(Species=="Lespedeza capitata")
diversityandro5ppm <- diversity5ppm %>%
  filter(Species=="Andropogon gerardii")

#Lespedeza capitata

mod8 <- lm(diversityaq5ppm ~ Diversity + Treatment, data=diversitylespe5ppm)
summary(mod8)
Anova(mod8)

#Andropogon gerardii

mod9 <- lm(diversityaq5ppm ~ Diversity + Treatment, data=diversityandro5ppm)
summary(mod9)
Anova(mod9)

###########################################################################

#Test whether any different composition signals using >5ppm from aqueous extract (phenolics)

aq5 <- read.csv("aqueousdatanormalised2_5ppm.csv", header=T)
aq5$Sample <- as.factor(aq5$Sample)
aq5$Sample <- as.character.factor(aq5$Sample)

variables5 <- aq5 %>%
  select(Sample, Species, Diversity, Treatment)

aq5.pca <- prcomp(aq5[ ,c(5:94)])
summary(aq5.pca)
aq5.pca.axes <- aq5.pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var="Sample")
aq5.pca.data <- full_join(variables5, aq5.pca.axes)

aq5.fullplot <- ggplot(aq5.pca.data, aes(x=PC1, y=PC2, colour=Species)) +
  geom_point(size=2.5, alpha=0.7) + 
  stat_ellipse(type="norm") +
  scale_colour_manual(values=c('#009E73', '#AA4499')) + 
  xlab("Aqueous PC1 (58.2%)") + ylab("Aqueous PC2 (20.6%)") + 
  theme_classic() +
  theme(legend.position = "top") 
aq5.fullplot

#Andro only 5ppm

androaq5 <- aq5 %>% 
  filter(Species=="Andropogon gerardii")

variablesandro5 <- androaq5 %>%
  select(Species, Diversity, Treatment) %>%
  rownames_to_column(var="Sample")

androaq5.pca <- prcomp(androaq5[ ,c(5:94)])
summary(androaq5.pca)
androaq5.pca.axes <- androaq5.pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var="Sample")
androaq5.pca.data <- full_join(variablesandro5, androaq5.pca.axes)

androaq5.plot <- ggplot(androaq5.pca.data, aes(x=PC1, y=PC2, color=Treatment, shape=Diversity)) +
  geom_point(size=3, alpha=0.7) +
  xlab("Aqueous PC1 (73.3%)") + ylab("Aqueous PC2 (8.4%)") +
  scale_colour_manual(values=c('#009E73', '#0de319ff')) + 
  theme_classic() +
  theme(legend.position = "top") 
androaq5.plot

#Lespe only 5ppm

lespeaq5 <- aq5 %>% 
  filter(Species=="Lespedeza capitata")

variableslespe5 <- lespeaq5 %>%
  select(Species, Diversity, Treatment) %>%
  rownames_to_column(var="Sample")

lespeaq5.pca <- prcomp(lespeaq5[ ,c(5:94)])
summary(lespeaq5.pca)
lespeaq5.pca.axes <- lespeaq5.pca$x %>%
  as.data.frame() %>% 
  rownames_to_column(var="Sample")
lespeaq5.pca.data <- full_join(variableslespe5, lespeaq5.pca.axes)

lespeaq5.plot <- ggplot(lespeaq5.pca.data, aes(x=PC1, y=PC2, color=Treatment, shape=Diversity)) +
  geom_point(size=3, alpha=0.7) +
  xlab("Aqueous PC1 (48.5%)") + ylab("Aqueous PC2 (19.1%)") +
  scale_colour_manual(values=c('#AA4499', '#6e4bf9ff')) + 
  theme_classic() +
  theme(legend.position = "top") 
lespeaq5.plot

#Testing these
#Nothing doing

adonis2(androaq5[ ,c(5:94)] ~ androaq5$Diversity + androaq5$Treatment, 
        permutations=10000, method="euclidean", by="margin")

adonis2(lespeaq5[ ,c(5:94)] ~ lespeaq5$Diversity + lespeaq5$Treatment, 
        permutations=10000,  method="euclidean", by="margin")







