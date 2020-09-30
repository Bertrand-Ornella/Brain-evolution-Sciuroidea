
######### Code - Can body mass and locomotion predict endocranial size? ############

#Libraries (possibly needs more libraries)
library(phytools) #open tree
library(ggplot2) #plots
library(nlme) # GLS analysis
library(AICcmodavg) # AIC 
library(RRPP) #pairwise comparaisons

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
ggplot(squirrel.data, aes(Brain, y = value, color = Family)) +
  theme_light() + theme(legend.position = "top") + 
  geom_point(data = dplyr::filter(squirrel.data, Family == "Sciuridae"), shape = 16, size = 3,
             aes(y = Body, color = "#4DBBD5FF"))+ 
  geom_point(data = dplyr::filter(squirrel.data, Family == "Aplodontidae"), shape = 16, size = 3,
             aes(y = Body, color = "#E64B35FF")) + 
  geom_point(data = dplyr::filter(squirrel.data, Family == "Ischyromyidae"), shape = 16, size = 3,
             aes(y = Body, color = "#3C5488FF")) +
  scale_color_manual(name = "", values = c("#3C5488FF","#4DBBD5FF","#E64B35FF"),labels = c("Ischyromyidae","Sciuridae", "Aplodontidae")) + 
  labs(x = "log(Body mass)", y = "log(Endocranial volume") 
  #add "+" at the end of the line above to get species names
  geom_text(data = dplyr::filter(squirrel.data, Family == "Sciuridae"), color = "#4DBBD5FF",
            aes(y = Body, label = abbreviation), hjust = -0.3, vjust = 1.1) +
  geom_text(data = dplyr::filter(squirrel.data, Family == "Ischyromyidae"), color = "#3C5488FF",
            aes(y = Body, label = abbreviation), hjust = -0.3, vjust = 1.1) +
  geom_text(data = dplyr::filter(squirrel.data, Family == "Aplodontidae"), color = "#E64B35FF",
            aes(y = Body, label = abbreviation), hjust = -0.3, vjust = 1.1) 

#Model using OLS (no phylogeny)
m.ols<-gls(Brain ~ Body + Locomotion, data=squirrel.data, method="ML")
m.ols.body<-gls(Brain ~ Body + Locomotion*Body, data=squirrel.data, method="ML")
m.ols.loco<-gls(Brain ~ Body*Locomotion, data=squirrel.data, method="ML") #same as above

#Chosing the model to use
anova(m.ols, m.ols.body, m.ols.loco) # m.ols.body (and m.ols.loco) with lowest AIC

#Residuals vs fitted plot
plot(fitted(m.ols.body), residuals(m.ols.body)) 
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
  geom_abline(intercept = 0, slope = 0) ## pattern, i.e., red ones are very low so need PGLS

##########################  Analyses -- PGLS ###############################

### Select PGLS model ### NOTE: I get a Warning message for all of the PGLS analyses

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

#Other method = Calculate AICs -- close to values from the anova
LL <- c(Brownian$logLik, Lambda$logLik, OU$logLik, Blomberg$logLik)
Ks <- c(length(Brownian$coefficients), length(Lambda$coefficients), 
        length(OU$coefficients), length(Blomberg$coefficients))
Modnames <- c("Brownian", "Lambda", "OU", "Blomberg")
aictabCustom(logL = LL, K = Ks, modnames = Modnames, nobs = Brownian$dims$N,sort = TRUE)

################# post ad-hoc test on locomotion ################

#tutorial RRPP: https://cran.r-project.org/web/packages/RRPP/vignettes/Using.RRPP.html

#Import squirrel data - to run the following
squirrel.data<-read.csv("squirrels_PEQ_res.csv", header=T)
Brain<-log10(squirrel.data$Brain_volume_cm3)
Body<-log10(squirrel.data$Body_mass_g)
Ecology<-squirrel.data$Ecology

#Pairwise comparisons - non-phylogenetic! (Weisbecker et al., 2019; line 552 in Analyses)
interaction_frame <- rrpp.data.frame(brain=Brain,body=Body,loco=Ecology)
fit<-lm.rrpp(brain~body+loco*body,SS.type = c("I"),data=interaction_frame)
summary(fit, formula = FALSE)
anova(fit) # When phylogeny is NOT taken into account, locomotion alone predicts brain size

Interactions <- pairwise(fit,covariate=interaction_frame$body,
                            groups=interaction_frame$loco)

summary(Interactions, test.type="dist") # shows which locomotor modes are significantly different

#Ex: Arboreal rodents have significantly different brain size compared to Fossorial
#Scansorial, and Terrestrial when accounting for body mass but not phylogeny

# Based on the fit model, how does brain size predicted to vary for each locomotor category? (what I understand this does)
sizeDF <- data.frame(loco = c("Arboreal", "Scansorial","Glider","Terrestrial","Fossorial"))
rownames(sizeDF) <- c("Arboreal", "Scansorial","Glider","Terrestrial","Fossorial")
sizePreds <- predict(fit,sizeDF)
plot(sizePreds) # The brain size of Arboreal and Gliding taxa are more similar than to the other locomotor categories

################ Make graph with PGLS regression lines for each locomotion ###########

#Both work 
#Import squirrel data - to run the following
Brain<-log10(squirrel.data$Brain_volume_cm3)
Body<-log10(squirrel.data$Body_mass_g)
#Or this one
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

#Create model PGLS regression line for each locomotion -- same warning messages as above
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
  geom_point(data = dplyr::filter(squirrel.data, Ecology == "Arboreal"), color = "#0EAF28",
             size = 2) +
  geom_point(data = dplyr::filter(squirrel.data, Ecology == "Fossorial"), color = "#975822",
             size = 2) +
  geom_point(data = dplyr::filter(squirrel.data, Ecology == "Glider"), color = "#73DAF3",
             size = 2) +
  geom_point(data = dplyr::filter(squirrel.data, Ecology == "Scansorial"), color = "#F31616",
             size = 2) + 
  geom_point(data = dplyr::filter(squirrel.data, Ecology == "Terrestrial"), color = "#F3B116",
             size = 2) +
  geom_line(data = dplyr::filter(predframe.Arb, Ecology == "Arboreal"), color = "#0EAF28",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.Fos, Ecology == "Fossorial"), color = "#975822",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.Gli, Ecology == "Glider"), color = "#73DAF3",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.Sca, Ecology == "Scansorial"), color = "#F31616",
            linetype = 1.5) +
  geom_line(data = dplyr::filter(predframe.Ter, Ecology == "Terrestrial"), color = "#F3B116",
            linetype = 1) +
  guides(size = FALSE, color = FALSE) +
  labs(x = "log(Body mass)", y = "log(Endocranial volume)") +
  theme_light() + 
  theme(legend.position = "top") +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 16,
            face = "bold")) 
#add "+" at the end of the line above to get species names
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

##### END! :)
  
  
  
  
