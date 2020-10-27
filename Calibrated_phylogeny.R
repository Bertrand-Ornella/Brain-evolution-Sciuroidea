#### Make a calibrated tree 

library("paleotree")
library("phytools")

#directory
setwd("~/Desktop/Squirrel_June_8_2020/Code")

#Import tree
tree_squirrel<-read.nexus("squirrel-tree_Meng.nex")

#### bin_cal3TimePaleoPhy method 

#taxa range for age
taxon.times<-read.csv("taxon_times3.csv",row.names=1)
int.times<-read.csv("int_times3.csv",row.names=1)

squirrel.range<-list(int.times,taxon.times)

#bin_cal3TimePaleoPhy method
likFun<-make_durationFreqDisc(squirrel.range)

spRes<-optim(parInit(likFun),likFun,lower=parLower(likFun),upper=parUpper(likFun),
             method="L-BFGS-B",control=list(maxit=1000000))

#sampling PROBABILITY per bin
sProb <- spRes[[1]][2]

#calculate meanInt. We need to use an average int.length (intervals not the same duration)
intLength<--apply(squirrel.range[[1]],1,diff)
hist(intLength)
meanInt <--mean(apply(squirrel.range[[1]], 1, diff)) #close to 1.8 Million years

sRate <- sProb2sRate(sProb,int.length = meanInt)
# we also need extinction rate and branching rate (see above)
# need to divide by int.length...

divRate <- spRes[[1]][1]/meanInt

#calibrated tree
tree_squirrel1 <- bin_cal3TimePaleoPhy(tree_squirrel, squirrel.range,brRate = divRate, extRate = divRate, 
                                       sampRate = sRate, ntrees = 100, plot = FALSE)
multiDiv(tree_squirrel1)

#Average tree - takes about ~7 mins
tree_squirrel3<-averageTree(tree_squirrel1, method="quadratic.path.difference")

#Meng tree
tree_squirrel2<-root(tree_squirrel3,outgroup=c("Paramys_copei","Paramys_delicatus",
                                               "Pseudotomus_horribilis","Pseudotomus_oweni",
                                               "Pseudotomus_petersoni","Pseudotomus_hians"),
                                                resolve.root=TRUE)

plot(tree_squirrel2)

############ Save trees

write.tree(tree_squirrel2, file = "Calibrated_tree_meng")
tree_squirrel2<-read.newick("Calibrated_tree_meng")
