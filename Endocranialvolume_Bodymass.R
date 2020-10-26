
######### Code - Can body mass and locomotion predict endocranial size? ############

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
squirrel.data<-read.csv("squirrels_PEQ_res.csv", header=T)

#Import tree
tree_squirrel<-read.newick("Calibrated_tree_meng")

#Transform data to log10
squirrel.data$Brain_volume_cm3<-log10(squirrel.data$Brain_volume_cm3)
names(squirrel.data)[names(squirrel.data) == "Brain_volume_cm3"] <- "Brain"

squirrel.data$Body_mass_g<-log10(squirrel.data$Body_mass_g)
names(squirrel.data)[names(squirrel.data) == "Body_mass_g"] <- "Body"

#Select other variables
Locomotion<-squirrel.data$Locomotion
abbreviation<-squirrel.data$abbreviation

##########################  Analyses -- OLS ###############################

# Look at the correlation among data point by Family
ggplot(squirrel.data, aes(Body, Brain, color = Family)) +
  theme_light() + theme(legend.position = "top") + 
  geom_point(data = dplyr::filter(squirrel.data, Family == "Sciuridae"), shape = 16, size = 3,
             aes(color = "#4DBBD5FF")) +
  geom_point(data = dplyr::filter(squirrel.data, Family == "Aplodontidae"), shape = 16, size = 3,
             aes(color = "#E64B35FF")) + 
  geom_point(data = dplyr::filter(squirrel.data, Family == "Ischyromyidae"), shape = 16, size = 3,
             aes(color = "#3C5488FF")) +
  scale_color_manual(name = "", values = c("#3C5488FF","#4DBBD5FF","#E64B35FF"),labels = c("Ischyromyidae","Sciuridae", "Aplodontidae")) + 
  labs(x = "log(Body mass)", y = "log(Endocranial volume)") +
  geom_text(data = dplyr::filter(squirrel.data, Family == "Sciuridae"), color = "#4DBBD5FF",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1) +
  geom_text(data = dplyr::filter(squirrel.data, Family == "Ischyromyidae"), color = "#3C5488FF",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1) +
  geom_text(data = dplyr::filter(squirrel.data, Family == "Aplodontidae"), color = "#E64B35FF",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1) 

#Model using OLS (no phylogeny)
m.ols<-gls(Brain ~ Body + Locomotion, data=squirrel.data, method="ML")
m.ols.body<-gls(Brain ~ Body + Locomotion*Body, data=squirrel.data, method="ML")
m.ols.loco<-gls(Brain ~ Body*Locomotion, data=squirrel.data, method="ML") #same as above

#Chosing the model to use
anova(m.ols, m.ols.body, m.ols.loco) # m.ols.body (and m.ols.loco) with lowest AIC
summary(m.ols.body) # Body mass alone + Interaction between Body mass and locomotion have a siginificant effect on brain size

#Residuals vs fitted plot
plot(fitted(m.ols.body), residuals(m.ols.body)) 
abline(0,0) #pattern visible so not good to use OLS

# Plot residuals ordered “by phylogeny” (Goldbogen et al., 2019; Supp. p.58)
is_tip <- tree_squirrel$edge[,2]<=length(tree_squirrel$tip.label)
ordered_tips <- tree_squirrel$edge[is_tip,2] # extract the order of tree tips 
oj <- residuals(m.ols.body)
tl <- tree_squirrel$tip.label[ordered_tips] # check order in tree to put them in same order in csv. file

#Plot residuals against phylogeny 
Spe<-squirrel.data$Sp
ggplot(squirrel.data, aes(x=Spe, y = oj, color = Family))+
  theme_light() + theme(legend.position = "top") +
  labs(x = "Species index", y = "OLS residuals") +
  geom_point() +
  geom_abline(intercept = 0, slope = 0) ## pattern, i.e., red ones are very low so need PGLS

##########################  Analyses -- PGLS ###############################

### Select PGLS model
# Model with locomotion only
Lambda<-gls(Brain ~ Body + Locomotion, data=squirrel.data, 
            correlation=corPagel(value=1,phy=tree_squirrel), method="ML")
Brownian<-gls(Brain ~ Body + Locomotion, data=squirrel.data, 
            correlation=corBrownian(1,phy=tree_squirrel), method="ML")
OU<-gls(Brain ~ Body + Locomotion, data=squirrel.data, 
            correlation=corMartins(1,phy=tree_squirrel), method="ML")
Blomberg<-gls(Brain ~ Body + Locomotion, data=squirrel.data, 
            correlation=corBlomberg(1.5,phy=tree_squirrel,fixed = TRUE), method="ML") #if 1, this is a Brownian model

# Model with interaction betwenn Locomotion and Body mass
Lambda1<-gls(Brain ~ Body + Locomotion*Body, data=squirrel.data, 
             correlation=corPagel(value=1,phy=tree_squirrel), method="ML")
Brownian1<-gls(Brain ~ Body + Locomotion*Body, data=squirrel.data, 
             correlation=corBrownian(1,phy=tree_squirrel), method="ML")
OU1<-gls(Brain ~ Body + Locomotion*Body, data=squirrel.data, 
             correlation=corMartins(1,phy=tree_squirrel, fixed = TRUE), method="ML")
Blomberg1<-gls(Brain ~ Body + Locomotion*Body, data=squirrel.data, 
             correlation=corBlomberg(1.5,phy=tree_squirrel,fixed = TRUE), method="ML") #if 1, this is a Brownian model

#Chosing the best model using AIC
anova(Lambda,Brownian,OU, Blomberg) # Lambda is the best 
anova(Lambda1,Brownian1,OU1, Blomberg1) # Lambda1 is the best
anova(Lambda,Lambda1) # Lambda1 is the best model overall
summary(Lambda1) #Final model

#Lambda = 0.72 (brain size shows some phylogenetic signal)
#Coefficients: p-value interpretation:
#Locomotion alone does not predict brain size (p-value = 0.12)
#Body mass alone and the interaction between Body mass and locomotion 
#can predict brain size (p-value = 0 and 0.02)

################# post ad-hoc test on locomotion ################

#tutorial RRPP: https://cran.r-project.org/web/packages/RRPP/vignettes/Using.RRPP.html

#Import squirrel data - to run pairwise test
squirrel.data<-read.csv("squirrels_PEQ_res.csv", header=T)
Brain<-log10(squirrel.data$Brain_volume_cm3)
Body<-log10(squirrel.data$Body_mass_g)
Ecology<-squirrel.data$Ecology

#Pairwise comparisons - non-phylogenetic! (Weisbecker et al., 2019; line 552 in Analyses)
interaction_frame <- rrpp.data.frame(brain=Brain,body=Body,loco=Ecology)
fit<-lm.rrpp(brain~body+loco*body,SS.type = c("I"),data=interaction_frame)
summary(fit, formula = FALSE)
anova(fit) # When phylogeny is NOT taken into account, all predict brain size

Interactions <- pairwise(fit,covariate=interaction_frame$body,
                            groups=interaction_frame$loco)

summary(Interactions, test.type="dist") # shows which locomotor modes are significantly different

#Ex: Arboreal rodents have significantly different brain size compared to Fossorial
#Scansorial, and Terrestrial when accounting for body mass but not phylogeny

# Based on the fit model, how does brain size predicted to vary for each locomotor category? (what I understand this does)
sizeDF <- data.frame(loco = c("Arboreal", "Scansorial","Glider","Terrestrial","Fossorial"))
rownames(sizeDF) <- c("Arboreal", "Scansorial","Glider","Terrestrial","Fossorial")
sizePreds <- predict(fit,sizeDF)
plot(sizePreds) # The brain size of Arboreal and Glider are more similar than with the other locomotor categories

################ Make graph with PGLS regression lines for each locomotion ###########

#Import squirrel data
squirrel.data<-read.csv("squirrels_PEQ_res.csv", header=T)
squirrel.data$Brain_volume_cm3<-log10(squirrel.data$Brain_volume_cm3)
names(squirrel.data)[names(squirrel.data) == "Brain_volume_cm3"] <- "Brain"
squirrel.data$Body_mass_g<-log10(squirrel.data$Body_mass_g)
names(squirrel.data)[names(squirrel.data) == "Body_mass_g"] <- "Body"

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
               "Petaurista_petaurista","Hylopetes_spadiceus","Petinomys_setosus",
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

#Create model PGLS regression line for each locomotion -- same warning
Arbline_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=Arb_Tree), data=Arb)
Fosline_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=Fos_Tree), data=Fos)
Gliline_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=Gli_Tree), data=Gli)
Scaline_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=Sca_Tree), data=Sca)
Terline_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=Ter_Tree), data=Ter)

#Prepare PGLS for each locomotor mode
pgls.fit.Arb <- predict(Arbline_Br_B) #predict values for brain size
predframe.Arb <- with(Arb, data.frame(Species, Ecology, Body, Brain = pgls.fit.Arb))

pgls.fit.Fos <- predict(Fosline_Br_B) #predict values for brain size
predframe.Fos <- with(Fos, data.frame(Species, Ecology, Body, Brain = pgls.fit.Fos))

pgls.fit.Gli <- predict(Gliline_Br_B) #predict values for brain size
predframe.Gli <- with(Gli, data.frame(Species, Ecology, Body, Brain = pgls.fit.Gli))

pgls.fit.Sca <- predict(Scaline_Br_B) #predict values for brain size
predframe.Sca <- with(Sca, data.frame(Species, Ecology, Body, Brain = pgls.fit.Sca))

pgls.fit.Ter <- predict(Terline_Br_B) #predict values for brain size
predframe.Ter <- with(Ter, data.frame(Species, Ecology, Body, Brain = pgls.fit.Ter))

#Make graph with PGLS corrected regressions
ggplot(squirrel.data, aes(Body, Brain, color = Ecology)) +
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
  labs(x = "log(Body mass)", y = "log(Endocranial volume)") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12,
            face = "bold")) +
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

########### Plot Arboreal rodents against other locomotor behaviours #############

#subset dataset per locomotion to use in model
Arb<- squirrel.data[ which(squirrel.data$EcoA=='Arboreal'), ]
NArb<- squirrel.data[ which(squirrel.data$EcoA=='Other'), ]

#Select taxa per locomotion for the tree
Arb_species<-c("Protoxerus_stangeri","Heliosciurus_rufobrachium","Callosciurus_sp",
               "Sciurus_carolinensis","Sciurus_granatensis","Tamiasciurus_hudsonicus",
               "Protosciurus_rachelae","Ratufa_affinis","Cedromus_wilsoni",
               "Prosciurus_relictus")

NArb_species<-c("Aplodontia_rufa","Mesogaulus_paniensis","Ischyromys_typus",
                "Pseudotomus_oweni","Pseudotomus_petersoni","Pseudotomus_hians",
                "Pseudotomus_horribilis","Pteromyscus_pulverulentus","Aeromys_tephromelas","Pteromys_volans",
                "Petaurista_petaurista","Hylopetes_spadiceus","Petinomys_setosus",
                "Glaucomys_volans","Paraxerus_cepapi","Funisciurus_pyrropus","Tamias_minimus",
                "Dremomys_rufigenis","Rapamys_atramontis","Reithroparamys_delicatissimus",
                "Paramys_copei","Paramys_delicatus","Xerus_rutilus","Cynomys_ludovicianus","Urocitellus_richardsonii",
                "Marmota_marmota","Lariscus_insignis","Rhinosciurus_laticaudatus")

#subset tree per locomotion
Arb_Tree<-drop.tip(tree_squirrel,tree_squirrel$tip.label[-match(Arb_species, tree_squirrel$tip.label)])
NArb_Tree<-drop.tip(tree_squirrel,tree_squirrel$tip.label[-match(NArb_species, tree_squirrel$tip.label)])

#Create model PGLS regression line for each locomotion -- same warning
Arbline_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=Arb_Tree), data=Arb)
NArbline_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=NArb_Tree), data=NArb)

#Prepare PGLS for each locomotor mode
pgls.fit.Arb <- predict(Arbline_Br_B) #predict values for brain size
predframe.Arb <- with(Arb, data.frame(Species, EcoA, Body, Brain = pgls.fit.Arb))

pgls.fit.NArb <- predict(NArbline_Br_B) #predict values for brain size
predframe.NArb <- with(NArb, data.frame(Species, EcoA, Body, Brain = pgls.fit.NArb))

#Make graph for arboreal vs. others
ggplot(squirrel.data, aes(Body, Brain, color = EcoA)) +
  geom_point(data = dplyr::filter(squirrel.data, EcoA == "Arboreal"),
             size = 2, aes(color = "#0EAF28")) +
  geom_point(data = dplyr::filter(squirrel.data, EcoA == "Other"),
             size = 2, aes(color = "#BDBAB8")) +
  geom_line(data = dplyr::filter(predframe.Arb, EcoA == "Arboreal"), color = "#0EAF28",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.NArb, EcoA == "Other"), color = "#BDBAB8",
            linetype = 1.5) +
  scale_color_manual(name = "", values = c("#0EAF28","#BDBAB8"),
            labels = c("Arboreal","Other"))+
  labs(x = "log(Body mass)", y = "log(Endocranial volume)") +
  theme_minimal() +
  #theme(legend.position = "top") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12,
            face = "bold")) +
  geom_text(data = dplyr::filter(squirrel.data, EcoA == "Arboreal"), color = "#0EAF28",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1)+ 
  geom_text(data = dplyr::filter(squirrel.data, EcoA == "Other"), color = "#BDBAB8",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1) 

########### Plot Fossorial rodents against other locomotor behaviours #############

#subset dataset per locomotion to use in model
Fos<- squirrel.data[ which(squirrel.data$EcoF=='Fossorial'), ]
NFos<- squirrel.data[ which(squirrel.data$EcoF=='Other'), ]

#Select taxa per locomotion for the tree
Fos_species<-c("Aplodontia_rufa","Mesogaulus_paniensis","Ischyromys_typus",
               "Pseudotomus_oweni","Pseudotomus_petersoni","Pseudotomus_hians",
               "Pseudotomus_horribilis")

NFos_species<-c("Protoxerus_stangeri","Heliosciurus_rufobrachium","Callosciurus_sp",
                "Sciurus_carolinensis","Sciurus_granatensis","Tamiasciurus_hudsonicus",
                "Protosciurus_rachelae","Ratufa_affinis","Cedromus_wilsoni",
                "Prosciurus_relictus","Pteromyscus_pulverulentus","Aeromys_tephromelas","Pteromys_volans",
                "Petaurista_petaurista","Hylopetes_spadiceus","Petinomys_setosus",
                "Glaucomys_volans","Paraxerus_cepapi","Funisciurus_pyrropus","Tamias_minimus",
                "Dremomys_rufigenis","Rapamys_atramontis","Reithroparamys_delicatissimus",
                "Paramys_copei","Paramys_delicatus","Xerus_rutilus","Cynomys_ludovicianus","Urocitellus_richardsonii",
                "Marmota_marmota","Lariscus_insignis","Rhinosciurus_laticaudatus")

#subset tree per locomotion
Fos_Tree<-drop.tip(tree_squirrel,tree_squirrel$tip.label[-match(Fos_species, tree_squirrel$tip.label)])
NFos_Tree<-drop.tip(tree_squirrel,tree_squirrel$tip.label[-match(NFos_species, tree_squirrel$tip.label)])

#Create model PGLS regression line for each locomotion -- same warning
Fosline_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=Fos_Tree), data=Fos)
NFosline_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=NFos_Tree), data=NFos)

#Prepare PGLS for each locomotor mode
pgls.fit.Fos <- predict(Fosline_Br_B) #predict values for brain size
predframe.Fos <- with(Fos, data.frame(Species, EcoF, Body, Brain = pgls.fit.Fos))

pgls.fit.NFos <- predict(NFosline_Br_B) #predict values for brain size
predframe.NFos <- with(NFos, data.frame(Species, EcoF, Body, Brain = pgls.fit.NFos))

#Make graph for fossorial vs. others
ggplot(squirrel.data, aes(Body, Brain, color = EcoF)) +
  geom_point(data = dplyr::filter(squirrel.data, EcoF == "Fossorial"),
             size = 2, aes(color = "#975822")) +
  geom_point(data = dplyr::filter(squirrel.data, EcoF == "Other"),
             size = 2, aes(color = "#BDBAB8")) +
  geom_line(data = dplyr::filter(predframe.Fos, EcoF == "Fossorial"), color = "#975822",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.NFos, EcoF == "Other"), color = "#BDBAB8",
            linetype = 1.5) +
  scale_color_manual(name = "", values = c("#975822","#BDBAB8"),
            labels = c("Fossorial","Other"))+
  labs(x = "log(Body mass)", y = "log(Endocranial volume)") +
  theme_minimal() + 
  #theme(legend.position = "top") +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12,
           face = "bold")) +
geom_text(data = dplyr::filter(squirrel.data, EcoF == "Fossorial"), color = "#975822",
          aes(label = abbreviation), hjust = -0.3, vjust = 1.1)+ 
  geom_text(data = dplyr::filter(squirrel.data, EcoF == "Other"), color = "#BDBAB8",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1) 

########### Plot Glider rodents against other locomotor behaviours #############

#subset dataset per locomotion to use in model
Gli<- squirrel.data[ which(squirrel.data$EcoG=='Glider'), ]
NGli<- squirrel.data[ which(squirrel.data$EcoG=='Other'), ]

#Select taxa per locomotion for the tree
Gli_species<-c("Pteromyscus_pulverulentus","Aeromys_tephromelas","Pteromys_volans",
                "Petaurista_petaurista","Hylopetes_spadiceus","Petinomys_setosus",
                "Glaucomys_volans")

NGli_species<-c("Aplodontia_rufa","Mesogaulus_paniensis","Ischyromys_typus",
               "Pseudotomus_oweni","Pseudotomus_petersoni","Pseudotomus_hians",
               "Pseudotomus_horribilis","Protoxerus_stangeri","Heliosciurus_rufobrachium","Callosciurus_sp",
                "Sciurus_carolinensis","Sciurus_granatensis","Tamiasciurus_hudsonicus",
                "Protosciurus_rachelae","Ratufa_affinis","Cedromus_wilsoni",
                "Prosciurus_relictus","Paraxerus_cepapi","Funisciurus_pyrropus","Tamias_minimus",
                "Dremomys_rufigenis","Rapamys_atramontis","Reithroparamys_delicatissimus",
                "Paramys_copei","Paramys_delicatus","Xerus_rutilus","Cynomys_ludovicianus","Urocitellus_richardsonii",
                "Marmota_marmota","Lariscus_insignis","Rhinosciurus_laticaudatus")

#subset tree per locomotion
Gli_Tree<-drop.tip(tree_squirrel,tree_squirrel$tip.label[-match(Gli_species, tree_squirrel$tip.label)])
NGli_Tree<-drop.tip(tree_squirrel,tree_squirrel$tip.label[-match(NGli_species, tree_squirrel$tip.label)])

#Create model PGLS regression line for each locomotion -- same warning
Gliline_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=Gli_Tree), data=Gli)
NGliline_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=NGli_Tree), data=NGli)

#Prepare PGLS for each locomotor mode
pgls.fit.Gli <- predict(Gliline_Br_B) #predict values for brain size
predframe.Gli <- with(Gli, data.frame(Species, EcoG, Body, Brain = pgls.fit.Gli))

pgls.fit.NGli <- predict(NGliline_Br_B) #predict values for brain size
predframe.NGli <- with(NGli, data.frame(Species, EcoG, Body, Brain = pgls.fit.NGli))

#Make graph for glider vs. others
ggplot(squirrel.data, aes(Body, Brain, color = EcoG)) +
  geom_point(data = dplyr::filter(squirrel.data, EcoG == "Glider"),
             size = 2, aes(color = "#73DAF3")) +
  geom_point(data = dplyr::filter(squirrel.data, EcoG == "Other"),
             size = 2, aes(color = "#BDBAB8")) +
  geom_line(data = dplyr::filter(predframe.Gli, EcoG == "Glider"), color = "#73DAF3",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.NGli, EcoG == "Other"), color = "#BDBAB8",
            linetype = 1.5) +
  scale_color_manual(name = "", values = c("#73DAF3","#BDBAB8"),
            labels = c("Glider","Other locomotor behaviours"))+
  labs(x = "log(Body mass)", y = "log(Endocranial volume)") +
  theme_light() + 
  theme(legend.position = "top") +
  #theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16,
            #face = "bold")) +
  geom_text(data = dplyr::filter(squirrel.data, EcoG == "Glider"), color = "#73DAF3",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1)+ 
  geom_text(data = dplyr::filter(squirrel.data, EcoG == "Other"), color = "#BDBAB8",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1) 

########### Plot Scansorial rodents against other locomotor behaviours #############

#subset dataset per locomotion to use in model
Sca<- squirrel.data[ which(squirrel.data$EcoS=='Scansorial'), ]
NSca<- squirrel.data[ which(squirrel.data$EcoS=='Other'), ]

#Select taxa per locomotion for the tree
Sca_species<-c("Paraxerus_cepapi","Funisciurus_pyrropus","Tamias_minimus",
               "Dremomys_rufigenis","Rapamys_atramontis","Reithroparamys_delicatissimus",
               "Paramys_copei","Paramys_delicatus")

NSca_species<-c("Pteromyscus_pulverulentus","Aeromys_tephromelas","Pteromys_volans",
                "Petaurista_petaurista","Hylopetes_spadiceus","Petinomys_setosus",
                "Glaucomys_volans","Aplodontia_rufa","Mesogaulus_paniensis","Ischyromys_typus",
                "Pseudotomus_oweni","Pseudotomus_petersoni","Pseudotomus_hians",
                "Pseudotomus_horribilis","Protoxerus_stangeri","Heliosciurus_rufobrachium","Callosciurus_sp",
                "Sciurus_carolinensis","Sciurus_granatensis","Tamiasciurus_hudsonicus",
                "Protosciurus_rachelae","Ratufa_affinis","Cedromus_wilsoni",
                "Prosciurus_relictus","Xerus_rutilus","Cynomys_ludovicianus","Urocitellus_richardsonii",
                "Marmota_marmota","Lariscus_insignis","Rhinosciurus_laticaudatus")

#subset tree per locomotion
Sca_Tree<-drop.tip(tree_squirrel,tree_squirrel$tip.label[-match(Sca_species, tree_squirrel$tip.label)])
NSca_Tree<-drop.tip(tree_squirrel,tree_squirrel$tip.label[-match(NSca_species, tree_squirrel$tip.label)])

#Create model PGLS regression line for each locomotion -- same warning ## error
Scaline_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=Sca_Tree), data=Sca)
NScaline_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=NSca_Tree), data=NSca)

#Prepare PGLS for each locomotor mode
pgls.fit.Sca <- predict(Scaline_Br_B) #predict values for brain size
predframe.Sca <- with(Sca, data.frame(Species, EcoS, Body, Brain = pgls.fit.Sca))

pgls.fit.NSca <- predict(NScaline_Br_B) #predict values for brain size
predframe.NSca <- with(NSca, data.frame(Species, EcoS, Body, Brain = pgls.fit.NSca))

#Make graph for scansorial vs. others
ggplot(squirrel.data, aes(Body, Brain, color = EcoS)) +
  geom_point(data = dplyr::filter(squirrel.data, EcoS == "Scansorial"),
             size = 2, aes(color = "#F31616")) +
  geom_point(data = dplyr::filter(squirrel.data, EcoS == "Other"),
             size = 2, aes(color = "#BDBAB8")) +
  geom_line(data = dplyr::filter(predframe.Sca, EcoS == "Scansorial"), color = "#F31616",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.NSca, EcoS == "Other"), color = "#BDBAB8",
            linetype = 1.5) +
  scale_color_manual(name = "", values = c("#BDBAB8", "#F31616"),
            labels = c("Other locomotor behaviours", "Scansorial"))+
  labs(x = "log(Body mass)", y = "log(Endocranial volume)") +
  theme_light() + 
  theme(legend.position = "top") +
  #theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16,
            #face = "bold")) +
  geom_text(data = dplyr::filter(squirrel.data, EcoS == "Scansorial"), color = "#F31616",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1)+ 
  geom_text(data = dplyr::filter(squirrel.data, EcoS == "Other"), color = "#BDBAB8",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1) 

########### Plot Terrestrial rodents against other locomotor behaviours #############
#subset dataset per locomotion to use in model
Ter<- squirrel.data[ which(squirrel.data$EcoT=='Terrestrial'), ]
NTer<- squirrel.data[ which(squirrel.data$EcoT=='Other'), ]

#Select taxa per locomotion for the tree
Ter_species<-c("Xerus_rutilus","Cynomys_ludovicianus","Urocitellus_richardsonii",
               "Marmota_marmota","Lariscus_insignis","Rhinosciurus_laticaudatus")

NTer_species<-c("Pteromyscus_pulverulentus","Aeromys_tephromelas","Pteromys_volans",
                "Petaurista_petaurista","Hylopetes_spadiceus","Petinomys_setosus",
                "Glaucomys_volans","Aplodontia_rufa","Mesogaulus_paniensis","Ischyromys_typus",
                "Pseudotomus_oweni","Pseudotomus_petersoni","Pseudotomus_hians",
                "Pseudotomus_horribilis","Protoxerus_stangeri","Heliosciurus_rufobrachium","Callosciurus_sp",
                "Sciurus_carolinensis","Sciurus_granatensis","Tamiasciurus_hudsonicus",
                "Protosciurus_rachelae","Ratufa_affinis","Cedromus_wilsoni",
                "Prosciurus_relictus","Paraxerus_cepapi","Funisciurus_pyrropus","Tamias_minimus",
                "Dremomys_rufigenis","Rapamys_atramontis","Reithroparamys_delicatissimus",
                "Paramys_copei","Paramys_delicatus")

#subset tree per locomotion
Ter_Tree<-drop.tip(tree_squirrel,tree_squirrel$tip.label[-match(Ter_species, tree_squirrel$tip.label)])
NTer_Tree<-drop.tip(tree_squirrel,tree_squirrel$tip.label[-match(NTer_species, tree_squirrel$tip.label)])

#Create model PGLS regression line for each locomotion -- same warning
Terline_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=Ter_Tree), data=Ter)
NTerline_Br_B <-gls(Brain ~ Body, correlation=corPagel (1,phy=NTer_Tree), data=NTer)

#Prepare PGLS for each locomotor mode
pgls.fit.Ter <- predict(Terline_Br_B) #predict values for brain size
predframe.Ter <- with(Ter, data.frame(Species, EcoT, Body, Brain = pgls.fit.Ter))

pgls.fit.NTer <- predict(NTerline_Br_B) #predict values for brain size
predframe.NTer <- with(NTer, data.frame(Species, EcoT, Body, Brain = pgls.fit.NTer))

#Make graph for terrestrial vs. others
ggplot(squirrel.data, aes(Body, Brain, color = EcoT)) +
  geom_point(data = dplyr::filter(squirrel.data, EcoT == "Terrestrial"),
             size = 2, aes(color = "#F3B116")) +
  geom_point(data = dplyr::filter(squirrel.data, EcoT == "Other"),
             size = 2, aes(color = "#BDBAB8")) +
  geom_line(data = dplyr::filter(predframe.Ter, EcoT == "Terrestrial"), color = "#F3B116",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.NTer, EcoT == "Other"), color = "#BDBAB8",
            linetype = 1.5) +
  scale_color_manual(name = "", values = c("#BDBAB8","#F3B116"),
            labels = c("Other locomotor behaviours", "Terrestrial"))+
  labs(x = "log(Body mass)", y = "log(Endocranial volume)") +
  theme_light() + 
  theme(legend.position = "top") +
  #theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16,
            #face = "bold")) +
  geom_text(data = dplyr::filter(squirrel.data, EcoT == "Terrestrial"), color = "#F3B116",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1)+ 
  geom_text(data = dplyr::filter(squirrel.data, EcoT == "Other"), color = "#BDBAB8",
            aes(label = abbreviation), hjust = -0.3, vjust = 1.1)

####################### pANCOVA - Different different allometric of one locomotor mode
###### againstother locomotions

#tutorial at https://smaerslab.com/gls-ancova/

#Import tree
tree_squirrel<-read.newick("Calibrated_tree_meng")
squirrel.data<-read.csv("squirrels_PEQ_res.csv", header=T, row.names = 1)

#Select data
data<-squirrel.data[,8:9]
Y<-"Brain_volume_cm3"
X<-"Body_mass_g"
data<-data[,c(which(colnames(data)==Y),which(colnames(data)==X)),drop=F]
data<-log(data)

tree<-treedata(tree_squirrel,data,sort=T,warnings=T)$phy   #match the tree to the data
data<-as.data.frame(treedata(tree,data,sort=T,warnings=T)$data)   #match the data to the tree
colnames(data)<-c("Dependent","Independent")

#Select taxa and get row number
Arboreal<-which(squirrel.data$Ecology=="Arboreal")
Fossorial<-which(squirrel.data$Ecology=="Fossorial")
Glider<-which(squirrel.data$Ecology=="Glider")
Scansorial<-which(squirrel.data$Ecology=="Scansorial")
Terrestrial<-which(squirrel.data$Ecology=="Terrestrial")

####### Arboreal 
#For differences in slope: 
grpSA<-rep("A",length(rownames(data))) 
grpSA[Arboreal]<-"B" 
grpSA<-as.factor(grpSA) 
names(grpSA)<-rownames(data)

#For differences in intercept:
grpIA<-rep("A",length(rownames(data))) 
grpIA[Arboreal]<-"B" 
grpIA<-as.factor(grpIA) 
names(grpIA)<-rownames(data)

ModelA<-model.matrix(as.formula(Dependent~Independent),data)

#Three types of models can be conceived that have more parameters than the 
#baseline model (those that vary slope but not intercept; vary intercept but 
#not slope; vary both):
Model_SA<-model.matrix(as.formula(Dependent~grpSA:Independent),data) # Differences in slopes, holding intercept constant
Model_IA<-model.matrix(as.formula(Dependent~grpIA + Independent),data) #Differences in intercept, holding slopes constant
Model_SIA<-model.matrix(as.formula(Dependent~grpIA + grpSA:Independent),data) # Differences in slopes and differences in intercept

#The following tests each of the above three models against the baseline model 
#of one slope and one intercept:
gls.ancova(Dependent~Independent,vcv(tree),ModelA,Model_SA) # Differences in slopes, holding intercept constant
gls.ancova(Dependent~Independent,vcv(tree),ModelA,Model_IA) # Differences in intercept, holding slopes constant
gls.ancova(Dependent~Independent,vcv(tree),ModelA,Model_SIA) # Differences in slopes and differences in intercept

####### Fossorial
#For differences in slope: 
grpSF<-rep("A",length(rownames(data))) 
grpSF[Fossorial]<-"B" 
grpSF<-as.factor(grpSF) 
names(grpSF)<-rownames(data)

#For differences in intercept:
grpIF<-rep("A",length(rownames(data))) 
grpIF[Fossorial]<-"B" 
grpIF<-as.factor(grpIF) 
names(grpIF)<-rownames(data)

ModelF<-model.matrix(as.formula(Dependent~Independent),data)

Model_SF<-model.matrix(as.formula(Dependent~grpSF:Independent),data) # Differences in slopes, holding intercept constant
Model_IF<-model.matrix(as.formula(Dependent~grpIF + Independent),data) #Differences in intercept, holding slopes constant
Model_SIF<-model.matrix(as.formula(Dependent~grpIF + grpSF:Independent),data) # Differences in slopes and differences in intercept

gls.ancova(Dependent~Independent,vcv(tree),ModelF,Model_SF) # Differences in slopes, holding intercept constant
gls.ancova(Dependent~Independent,vcv(tree),ModelF,Model_IF) # Differences in intercept, holding slopes constant
gls.ancova(Dependent~Independent,vcv(tree),ModelF,Model_SIF) # Differences in slopes and differences in intercept

####### Glider
#For differences in slope: 
grpSG<-rep("A",length(rownames(data))) 
grpSG[Glider]<-"B" 
grpSG<-as.factor(grpSG) 
names(grpSG)<-rownames(data)

#For differences in intercept:
grpIG<-rep("A",length(rownames(data))) 
grpIG[Glider]<-"B" 
grpIG<-as.factor(grpIG) 
names(grpIG)<-rownames(data)

ModelG<-model.matrix(as.formula(Dependent~Independent),data)

Model_SG<-model.matrix(as.formula(Dependent~grpSG:Independent),data) # Differences in slopes, holding intercept constant
Model_IG<-model.matrix(as.formula(Dependent~grpIG + Independent),data) #Differences in intercept, holding slopes constant
Model_SIG<-model.matrix(as.formula(Dependent~grpIG + grpSG:Independent),data) # Differences in slopes and differences in intercept

gls.ancova(Dependent~Independent,vcv(tree),ModelG,Model_SG) # Differences in slopes, holding intercept constant
gls.ancova(Dependent~Independent,vcv(tree),ModelG,Model_IG) # Differences in intercept, holding slopes constant
gls.ancova(Dependent~Independent,vcv(tree),ModelG,Model_SIG) # Differences in slopes and differences in intercept

####### Scansorial
#For differences in slope: 
grpSS<-rep("A",length(rownames(data))) 
grpSS[Scansorial]<-"B" 
grpSS<-as.factor(grpSS) 
names(grpSS)<-rownames(data)

#For differences in intercept:
grpIS<-rep("A",length(rownames(data))) 
grpIS[Scansorial]<-"B" 
grpIS<-as.factor(grpIS) 
names(grpIS)<-rownames(data)

ModelS<-model.matrix(as.formula(Dependent~Independent),data)

Model_SS<-model.matrix(as.formula(Dependent~grpSS:Independent),data) # Differences in slopes, holding intercept constant
Model_IS<-model.matrix(as.formula(Dependent~grpIS + Independent),data) #Differences in intercept, holding slopes constant
Model_SIS<-model.matrix(as.formula(Dependent~grpIS + grpSS:Independent),data) # Differences in slopes and differences in intercept

gls.ancova(Dependent~Independent,vcv(tree),ModelS,Model_SS) # Differences in slopes, holding intercept constant
gls.ancova(Dependent~Independent,vcv(tree),ModelS,Model_IS) # Differences in intercept, holding slopes constant
gls.ancova(Dependent~Independent,vcv(tree),ModelS,Model_SIS) # Differences in slopes and differences in intercept

####### Terrestrial
#For differences in slope: 
grpST<-rep("A",length(rownames(data))) 
grpST[Terrestrial]<-"B" 
grpST<-as.factor(grpST) 
names(grpST)<-rownames(data)

#For differences in intercept:
grpIT<-rep("A",length(rownames(data))) 
grpIT[Terrestrial]<-"B" 
grpIT<-as.factor(grpIT) 
names(grpIT)<-rownames(data)

ModelT<-model.matrix(as.formula(Dependent~Independent),data)

Model_ST<-model.matrix(as.formula(Dependent~grpST:Independent),data) # Differences in slopes, holding intercept constant
Model_IT<-model.matrix(as.formula(Dependent~grpIT + Independent),data) #Differences in intercept, holding slopes constant
Model_SIT<-model.matrix(as.formula(Dependent~grpIT + grpST:Independent),data) # Differences in slopes and differences in intercept

gls.ancova(Dependent~Independent,vcv(tree),ModelT,Model_ST) # Differences in slopes, holding intercept constant
gls.ancova(Dependent~Independent,vcv(tree),ModelT,Model_IT) # Differences in intercept, holding slopes constant
gls.ancova(Dependent~Independent,vcv(tree),ModelT,Model_SIT) # Differences in slopes and differences in intercept

##### END
  
  
  
  