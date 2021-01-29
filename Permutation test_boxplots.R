#Permutation test

library(coin) #one_way test
library(rcompanion) #pairwisePermutationTest
library(ggpubr) #ggboxplot
library(car) #Levene's test

#directory
setwd("~/Desktop/Squirrel_June_8_2020/Code")

#Import squirrel data
squirrel.data1<-read.csv("squirrels_PEQ_res.csv", header=T)
squirrel.data<-squirrel.data1[, c("PEQ_Meng", "Ecology")]
squirrel.data2<-squirrel.data1[, c("OB_percentage", "Ecology")]
squirrel.data3<-squirrel.data1[, c("PL_percentage", "Ecology")]
squirrel.data4<-squirrel.data1[, c("Neocortex_surface_percentage", "Ecology")]

comp<- list(c("Arboreal", "Fossorial"),c("Arboreal", "Scansorial"), 
            c("Arboreal", "Terrestrial"),c("Fossorial", "Glider"),c("Fossorial", "Scansorial"))

###############PEQ permutation test

##ggplot - boxplot - PEQ
ggboxplot(squirrel.data1,x="Ecology", y="PEQ_Meng", fill="Ecology", palette=c("chartreuse3", "chocolate4", "steelblue1","firebrick2","yellow"))+
  geom_point() + 
  geom_text(aes(label = abbreviation ), hjust = 0, nudge_x = 0.05)+
  labs(x='Locomotion', y='Phylogenetic encephalization quotient')

#Fisher-Pitman permutation test
oneway_test(PEQ_Meng~Ecology,data=squirrel.data)

squirrel.data$Ecology = factor(squirrel.data$Ecology, levels = c("Arboreal","Fossorial","Glider",
                      "Scansorial","Terrestrial")) 

PT_PEQ<-pairwisePermutationTest(PEQ_Meng~Ecology,data=squirrel.data,method="fdr")
PT_PEQ

#Test if data normally distributed (Normally distributed)
shapiro.test(squirrel.data$PEQ_Meng)

# Bartlett test when data are normally distrubuted (YES there is homogeneity of variances)
bartlett.test(PEQ_Meng ~ Ecology, data = squirrel.data)

############Olfactory bulb permutation test

##ggplot - boxplot - OB
ggboxplot(squirrel.data1,x="Ecology", y="OB_percentage", fill="Ecology", palette=c("chartreuse3", "chocolate4", "steelblue1","firebrick2","yellow"))+
  geom_point() + 
  geom_text(aes(label = abbreviation ), hjust = 0, nudge_x = 0.05)+
  labs(x='Locomotion', y='Olfactory bulb volume percentage')

#Fisher-Pitman permutation test
oneway_test(OB_percentage~Ecology,data=squirrel.data2)

squirrel.data2$Ecology = factor(squirrel.data2$Ecology, levels = c("Arboreal","Fossorial","Glider",
                                                                 "Scansorial","Terrestrial")) 

PT_OB<-pairwisePermutationTest(OB_percentage~Ecology,data=squirrel.data2,method="fdr")
PT_OB

#Test if data normally distributed (not normally distributed)
shapiro.test(squirrel.data2$OB_percentage)

# Levene's test when data are not normally distrubuted (YES there is homogeneity of variances)
leveneTest(OB_percentage ~ Ecology, data = squirrel.data2)

#############Petrosal lobule permutation test

##ggplot - boxplot - Petrosal lobules
ggboxplot(squirrel.data1,x="Ecology", y="PL_percentage", fill="Ecology", palette=c("chartreuse3", "chocolate4", "steelblue1","firebrick2","yellow"))+
  geom_point() + 
  geom_text(aes(label = abbreviation ), hjust = 0, nudge_x = 0.05)+
  labs(x='Locomotion', y='Petrosal lobule volume percentage')

#Fisher-Pitman permutation test
oneway_test(PL_percentage~Ecology,data=squirrel.data3)

squirrel.data3$Ecology = factor(squirrel.data3$Ecology, levels = c("Arboreal","Fossorial","Glider",
                                                                 "Scansorial","Terrestrial")) 

PT_PL<-pairwisePermutationTest(PL_percentage~Ecology,data=squirrel.data3,method="fdr")
PT_PL

#Test if data normally distributed (Normally distributed)
shapiro.test(squirrel.data3$PL_percentage)

# Bartlett test when data are normally distrubuted (YES there is homogeneity of variances)
bartlett.test(PL_percentage ~ Ecology, data = squirrel.data3)

############ Neocortex permutation test

##ggplot - boxplot - Neocortex
ggboxplot(squirrel.data1,x="Ecology", y="Neocortex_surface_percentage", fill="Ecology", palette=c("chartreuse3", "chocolate4", "steelblue1","firebrick2","yellow"))+
  geom_point() + 
  geom_text(aes(label = abbreviation ), hjust = 0, nudge_x = 0.05)+
  labs(x='Locomotion', y='Neocortical surface area percentage')

#Fisher-Pitman permutation test
oneway_test(Neocortex_surface_percentage~Ecology,data=squirrel.data4)
#independence_test(Neocortex_surface_percentage~Ecology,data=squirrel.data)

squirrel.data4$Ecology = factor(squirrel.data4$Ecology, levels = c("Arboreal","Fossorial","Glider",
                                                                 "Scansorial","Terrestrial")) 

PT_Neo<-pairwisePermutationTest(Neocortex_surface_percentage~Ecology,data=squirrel.data4,method="fdr")
PT_Neo

#Test if data normally distributed (Not normally distributed)
shapiro.test(squirrel.data4$Neocortex_surface_percentage)

# Levene's test when data are not normally distrubuted (NO there is not homogeneity of variances)
leveneTest(Neocortex_surface_percentage ~ Ecology, data = squirrel.data4)

### END


