
######### Code - Can endocranial volume and locomotion predict petrosal lobule size? ############

#Libraries (possibly needs more libraries)
library(phytools) #open tree
library(ggplot2) #plots
library(nlme) # GLS analysis
library(AICcmodavg) # AIC 
library(RRPP) #pairwise comparaisons
library(geiger) # pANCOVA
library(evomap) # pANCOVA
require(devtools)
install_github("JeroenSmaers/evomap")

#directory
setwd("~/Desktop/Squirrel_June_8_2020/Code")

#Import squirrel data
squirrel.data1<-read.csv("squirrels_PEQ_res.csv", header=T)

#Delete a specific row and taxon in tree
squirrel.data<-squirrel.data1[-c(23),]

#Import tree
tree_squirrel1<-read.newick("Calibrated_tree_meng")

#delete taxa
tree_squirrel<-drop.tip(tree_squirrel1, c("Petinomys_setosus"))

#Transform data to log10
squirrel.data$Petrosal_lobule_volume_mm3<-log10(squirrel.data$Petrosal_lobule_volume_mm3)
names(squirrel.data)[names(squirrel.data) == "Petrosal_lobule_volume_mm3"] <- "PL"

squirrel.data$Brain_volume_mm3<-log10(squirrel.data$Brain_volume_mm3)
names(squirrel.data)[names(squirrel.data) == "Brain_volume_mm3"] <- "Brain.mg"

#Select other variables
Locomotion<-squirrel.data$Locomotion
abbreviation<-squirrel.data$abbreviation

##########################  Analyses -- OLS ###############################

# Look at the correlation among data point by Family
ggplot(squirrel.data, aes(Brain.mg, PL, color = Family)) +
  theme_light() + theme(legend.position = "top") + 
  geom_point(data = dplyr::filter(squirrel.data, Family == "Sciuridae"), shape = 16, size = 3,
             aes(color = "#4DBBD5FF")) +
  geom_point(data = dplyr::filter(squirrel.data, Family == "Aplodontidae"), shape = 16, size = 3,
             aes(color = "#E64B35FF")) + 
  geom_point(data = dplyr::filter(squirrel.data, Family == "Ischyromyidae"), shape = 16, size = 3,
             aes(color = "#3C5488FF")) +
  scale_color_manual(name = "", values = c("#3C5488FF","#4DBBD5FF","#E64B35FF"),labels = c("Ischyromyidae","Sciuridae", "Aplodontidae")) + 
  labs(x = "log(Body mass)", y = "log(Petrosal lobule volume)") +
  geom_text(data = dplyr::filter(squirrel.data, Family == "Sciuridae"), color = "#4DBBD5FF",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1) +
  geom_text(data = dplyr::filter(squirrel.data, Family == "Ischyromyidae"), color = "#3C5488FF",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1) +
  geom_text(data = dplyr::filter(squirrel.data, Family == "Aplodontidae"), color = "#E64B35FF",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1) 

#Model using OLS (no phylogeny)
m.ols<-gls(PL ~ Brain.mg + Locomotion, data=squirrel.data, method="ML")
m.ols.body<-gls(PL ~ Brain.mg + Locomotion*Brain.mg, data=squirrel.data, method="ML")
m.ols.loco<-gls(PL ~ Brain.mg*Locomotion, data=squirrel.data, method="ML") #same as above

#Chosing the model to use
anova(m.ols, m.ols.body, m.ols.loco) # m.ols with lowest AIC
summary(m.ols)

#Residuals vs fitted plot
plot(fitted(m.ols), residuals(m.ols)) 
abline(0,0) #pattern visible so not good to use OLS

# Plot residuals ordered “by phylogeny” (Goldbogen et al., 2019; Supp. p.58)
is_tip <- tree_squirrel$edge[,2]<=length(tree_squirrel$tip.label)
ordered_tips <- tree_squirrel$edge[is_tip,2] # extract the order of tree tips 
oj <- residuals(m.ols)
tl <- tree_squirrel$tip.label[ordered_tips] # check order in tree to put them in same order in csv. file

#Plot residuals against phylogeny 
Spe<-squirrel.data$Sp
ggplot(squirrel.data, aes(x=Spe, y = oj, color = Family))+
  theme_light() + theme(legend.position = "top") +
  labs(x = "Species index", y = "OLS residuals") +
  geom_point() +
  geom_abline(intercept = 0, slope = 0) ## pattern, i.e., green group lower compared to blue group

##########################  Analyses -- PGLS ###############################

### Select PGLS model

# Model with locomotion only
Lambda<-gls(PL ~ Brain.mg + Locomotion, data=squirrel.data, 
            correlation=corPagel(value=1,phy=tree_squirrel), method="ML")
Brownian<-gls(PL ~ Brain.mg + Locomotion, data=squirrel.data, 
              correlation=corBrownian(1,phy=tree_squirrel), method="ML")
OU<-gls(PL ~ Brain.mg + Locomotion, data=squirrel.data, 
        correlation=corMartins(1,phy=tree_squirrel,fixed = TRUE), method="ML")
Blomberg<-gls(PL ~ Brain.mg + Locomotion, data=squirrel.data, 
              correlation=corBlomberg(1.5,phy=tree_squirrel,fixed = TRUE), method="ML") #if 1, this is a Brownian model

# Model with interaction betwenn Locomotion and Body mass
Lambda1<-gls(PL ~ Brain.mg + Locomotion*Brain.mg, data=squirrel.data, 
             correlation=corPagel(value=1,phy=tree_squirrel), method="ML")
Brownian1<-gls(PL ~ Brain.mg + Locomotion*Brain.mg, data=squirrel.data, 
               correlation=corBrownian(1,phy=tree_squirrel), method="ML")
OU1<-gls(PL ~ Brain.mg + Locomotion*Brain.mg, data=squirrel.data, 
         correlation=corMartins(1,phy=tree_squirrel, fixed = TRUE), method="ML")
Blomberg1<-gls(PL ~ Brain.mg + Locomotion*Brain.mg, data=squirrel.data, 
               correlation=corBlomberg(1.5,phy=tree_squirrel,fixed = TRUE), method="ML") #if 1, this is a Brownian model

#Chosing the best model using AIC
anova(Lambda,Brownian,OU, Blomberg) # Lambda is the best 
anova(Lambda1,Brownian1,OU1, Blomberg1)
anova(Lambda,Lambda1) #Final model
summary(Lambda)

#Lambda = 0.54 (PL relative size shows some phylogenetic signal)
#Coefficients: p-value interpretation:
#Endocranial volume and locomotion alone can predict brain size (p-value = 0 and 0.002)

################# post ad-hoc test on locomotion ################

#tutorial RRPP: https://cran.r-project.org/web/packages/RRPP/vignettes/Using.RRPP.html

#Import squirrel data - to run the following
squirrel.data1<-read.csv("squirrels_PEQ_res.csv", header=T)
squirrel.data<-squirrel.data1[-c(23),]
tree_squirrel1<-read.newick("Calibrated_tree_meng")
tree_squirrel<-drop.tip(tree_squirrel1, c("Petinomys_setosus"))

PL<-log10(squirrel.data$Petrosal_lobule_volume_mm3)
Brain<-log10(squirrel.data$Brain_volume_mm3)
Ecology<-squirrel.data$Ecology

#Pairwise comparisons - non-phylogenetic! (Weisbecker et al., 2019; line 552 in Analyses)
interaction_frame <- rrpp.data.frame(pl=PL,brain=Brain,loco=Ecology)
fit<-lm.rrpp(pl~brain+loco,SS.type = c("I"),data=interaction_frame)
summary(fit, formula = FALSE)
anova(fit)

Interactions <- pairwise(fit,covariate=interaction_frame$body,
                         groups=interaction_frame$loco)

summary(Interactions, test.type="dist") # shows which locomotor modes are significantly different

# Based on the fit model, how does brain size predicted to vary for each locomotor category? (what I understand this does)
sizeDF <- data.frame(loco = c("Arboreal", "Scansorial","Glider","Terrestrial","Fossorial"))
rownames(sizeDF) <- c("Arboreal", "Scansorial","Glider","Terrestrial","Fossorial")
sizePreds <- predict(fit,sizeDF)
plot(sizePreds) # The brain size of Arboreal and Glider are more similar than with the other locomotor categories

################ Make graph with PGLS regression lines for each locomotion ###########

#Import data
squirrel.data1<-read.csv("squirrels_PEQ_res.csv", header=T)
squirrel.data<-squirrel.data1[-c(23),]

#Import tree
tree_squirrel1<-read.newick("Calibrated_tree_meng")
tree_squirrel<-drop.tip(tree_squirrel1, c("Petinomys_setosus"))

#Transform data to log10
squirrel.data$Petrosal_lobule_volume_mm3<-log10(squirrel.data$Petrosal_lobule_volume_mm3)
names(squirrel.data)[names(squirrel.data) == "Petrosal_lobule_volume_mm3"] <- "PL"
squirrel.data$Brain_volume_mm3<-log10(squirrel.data$Brain_volume_mm3)
names(squirrel.data)[names(squirrel.data) == "Brain_volume_mm3"] <- "Brain.mg"

#subset dataset per locomotion to use in model
Arb<- squirrel.data[ which(squirrel.data$Ecology=='Arboreal'), ]
Fos<- squirrel.data[ which(squirrel.data$Ecology=='Fossorial'), ]
Gli<- squirrel.data[ which(squirrel.data$Ecology=='Glider'), ]
Sca<- squirrel.data[ which(squirrel.data$Ecology=='Scansorial'), ]
Ter<- squirrel.data[ which(squirrel.data$Ecology=='Terrestrial'), ]

#Find taxa with specific locomotion
Arb_sp<-subset(squirrel.data, Ecology=="Arboreal", select=Species)
Fos_sp<-subset(squirrel.data, Ecology=="Fossorial", select=Species)
Gli_sp<-subset(squirrel.data, Ecology=="Glider", select=Species)
Sca_sp<-subset(squirrel.data, Ecology=="Scansorial", select=Species)
Ter_sp<-subset(squirrel.data, Ecology=="Terrestrial", select=Species)

#Select taxa per locomotion for the tree
Arb_species<-c("Protoxerus_stangeri","Heliosciurus_rufobrachium","Callosciurus_sp",
               "Sciurus_carolinensis","Sciurus_granatensis","Tamiasciurus_hudsonicus",
               "Protosciurus_rachelae","Ratufa_affinis","Cedromus_wilsoni",
               "Prosciurus_relictus")

Fos_species<-c("Aplodontia_rufa","Mesogaulus_paniensis","Ischyromys_typus",
               "Pseudotomus_oweni","Pseudotomus_petersoni","Pseudotomus_hians",
               "Pseudotomus_horribilis")

Gli_species<-c("Pteromyscus_pulverulentus","Aeromys_tephromelas","Pteromys_volans",
               "Petaurista_petaurista","Hylopetes_spadiceus",
               "Glaucomys_volans")

Sca_species<-c("Paraxerus_cepapi","Funisciurus_pyrropus","Tamias_minimus",
               "Dremomys_rufigenis","Rapamys_atramontis","Reithroparamys_delicatissimus",
               "Paramys_copei","Paramys_delicatus")

Ter_species<-c("Xerus_rutilus","Cynomys_ludovicianus","Urocitellus_richardsonii",
               "Marmota_marmota","Lariscus_insignis","Rhinosciurus_laticaudatus")

#subset tree per locomotion
Arb_Tree<-drop.tip(tree_squirrel,tree_squirrel$tip.label[-match(Arb_species, tree_squirrel$tip.label)])
Fos_Tree<-drop.tip(tree_squirrel,tree_squirrel$tip.label[-match(Fos_species, tree_squirrel$tip.label)])
Gli_Tree<-drop.tip(tree_squirrel,tree_squirrel$tip.label[-match(Gli_species, tree_squirrel$tip.label)])
Sca_Tree<-drop.tip(tree_squirrel,tree_squirrel$tip.label[-match(Sca_species, tree_squirrel$tip.label)])
Ter_Tree<-drop.tip(tree_squirrel,tree_squirrel$tip.label[-match(Ter_species, tree_squirrel$tip.label)])

#Create model PGLS regression line for each locomotion
Arbline_Br_B <-gls(PL ~ Brain.mg, correlation=corPagel (1,phy=Arb_Tree,fixed=T), data=Arb)
Fosline_Br_B <-gls(PL ~ Brain.mg, correlation=corPagel (1,phy=Fos_Tree), data=Fos)
Gliline_Br_B <-gls(PL ~ Brain.mg, correlation=corPagel (1,phy=Gli_Tree), data=Gli)
Scaline_Br_B <-gls(PL ~ Brain.mg, correlation=corPagel (1,phy=Sca_Tree), data=Sca)
Terline_Br_B <-gls(PL ~ Brain.mg, correlation=corPagel (1,phy=Ter_Tree), data=Ter)

#Prepare PGLS for each locomotor mode
pgls.fit.Arb <- predict(Arbline_Br_B)
predframe.Arb <- with(Arb, data.frame(Species, Ecology, Brain.mg, PL = pgls.fit.Arb))

pgls.fit.Fos <- predict(Fosline_Br_B)
predframe.Fos <- with(Fos, data.frame(Species, Ecology, Brain.mg, PL = pgls.fit.Fos))

pgls.fit.Gli <- predict(Gliline_Br_B)
predframe.Gli <- with(Gli, data.frame(Species, Ecology, Brain.mg, PL = pgls.fit.Gli))

pgls.fit.Sca <- predict(Scaline_Br_B)
predframe.Sca <- with(Sca, data.frame(Species, Ecology, Brain.mg, PL = pgls.fit.Sca))

pgls.fit.Ter <- predict(Terline_Br_B)
predframe.Ter <- with(Ter, data.frame(Species, Ecology, Brain.mg, PL = pgls.fit.Ter))

#Make graph with PGLS corrected regressions
ggplot(squirrel.data, aes(Brain.mg, PL, color = Ecology)) +
  geom_point(data = dplyr::filter(squirrel.data, Ecology == "Arboreal"),
             size = 2, aes(color = "#0EAF28")) +
  geom_point(data = dplyr::filter(squirrel.data, Ecology == "Fossorial"),
             size = 2, aes(color = "#975822")) +
  geom_point(data = dplyr::filter(squirrel.data, Ecology == "Glider"),
             size = 2, aes(color = "#73DAF3")) +
  geom_point(data = dplyr::filter(squirrel.data, Ecology == "Scansorial"), 
             size = 2, aes(color = "#F31616")) + 
  geom_point(data = dplyr::filter(squirrel.data, Ecology == "Terrestrial"), 
             size = 2, aes(color = "#F3B116")) +
  geom_line(data = dplyr::filter(predframe.Arb, Ecology == "Arboreal"), color = "#0EAF28",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.Fos, Ecology == "Fossorial"), color = "#975822",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.Gli, Ecology == "Glider"), color = "#73DAF3",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.Sca, Ecology == "Scansorial"), color = "#F31616",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.Ter, Ecology == "Terrestrial"), color = "#F3B116",
            linetype = 1.5) +
  theme_minimal() + 
  #theme(legend.position = "top") +
  scale_color_manual(name = "", values = c("#0EAF28","#73DAF3","#975822", "#F31616", "#F3B116"),labels = c("Arboreal",
            "Glider","Fossorial", "Scansorial", "Terrestrial")) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12,
              face = "bold")) +
  labs(x = "log(Endocranial volume)", y = "log(Petrosal lobule volume)") +
  geom_text(data = dplyr::filter(squirrel.data, Ecology == "Arboreal"), color = "#0EAF28",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1)+ 
  geom_text(data = dplyr::filter(squirrel.data, Ecology == "Fossorial"), color = "#975822",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1) +           
  geom_text(data = dplyr::filter(squirrel.data, Ecology == "Glider"), color = "#73DAF3",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1) +
  geom_text(data = dplyr::filter(squirrel.data, Ecology == "Scansorial"), color = "#F31616",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1) +
  geom_text(data = dplyr::filter(squirrel.data, Ecology == "Terrestrial"), color = "#F3B116",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1)

##### END!
