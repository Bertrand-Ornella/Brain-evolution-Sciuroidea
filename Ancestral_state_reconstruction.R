library("ape")
library("strap") #geological scale
library("phytools") #ASR

#directory
setwd("~/Desktop/Squirrel_June_8_2020/Code")

#Import squirrel data 
squirrel.data2<-read.csv("squirrels_PEQ_res.csv", header=T, row.names = 1)

#Import tree
tree_squirrel2<-read.newick("Calibrated_tree_meng")

######## Node numbering ############

### Node number on non-calibrated tree
tree_squirrel<-read.nexus("squirrel-tree_Meng.nex")
plotTree(tree_squirrel,pts=F,node.numbers=T,fsize=0.6)

### Node number on non-calibrated tree - without Petinomys
tree_squirrelEX<-drop.tip(tree_squirrel, c("Petinomys_setosus"))
plotTree(tree_squirrelEX,pts=F,node.numbers=T,fsize=0.6)

####Ancestral state reconstructions - Continous characters

squirrel.data<-squirrel.data2[,5:20] # Meng

for (i in 1:ncol(squirrel.data)){
  
  dat1<-as.matrix(squirrel.data)[,i]
  dat2<-na.exclude(dat1)#remove Nas
  dat3<-as.numeric(dat2)
  dat3<-setNames(dat3, names(dat2))
  
  #Drop the missing taxa 
  if (anyNA(dat1, recursive = FALSE)==TRUE){
    phy<-drop.tip(tree_squirrel2, names(which(is.na(dat1))))
    anc.eq<-fastAnc(phy,dat3,vars=TRUE,CI=TRUE) # estimate states  with variances & 95% confidence intervals for each node
    
    num.anc<-anc.eq$ace
    num.anc2<-round(num.anc,digits =2)
    write.table(num.anc2, paste(colnames(squirrel.data)[i],"_ancestral.csv"))
    
    #plot ancestral states in the tree
    
    #Represent ancestral states as colors 
    plt<-contMap(phy,dat3, plot=FALSE)
    plt<-setMap(plt, colors=c("black","darkmagenta","violetred2", "orange", "lightgoldenrod1"))
  
    # Open a pdf file
    pdf(paste(colnames(squirrel.data[i]),"plot.pdf"),width = 18.27, height = 11.69, paper= "a4r")
    
    # 2. Create a plot
    plot(plt,legend=FALSE,ylim=c(1-0.09*(Ntip(plt$tree)-1),Ntip(plt$tree)))
    
    #Bar for ancestral state range
    add.color.bar(leg=0.3*max(nodeHeights(phy)),plt$cols,colnames(squirrel.data[i]),
                  lims=plt$lims,digits=3,prompt=FALSE,x=20,
                  y=1-0.08*(Ntip(plt$tree)-60),lwd=7,fsize=0.8,subtitle="")
    
    #Bar time scale
    phy$root.time <- 53
    
    time.tree<- geoscalePhylo(tree=ladderize(phy,right=FALSE), units=c("Period", "Epoch", "Age"), boxes="Epoch",
                              cex.tip=0.8, cex.age=0.8, cex.ts=1.0, label.offset=0, x.lim=c(0,53), lwd=3, width=2)
    
    #Label for nodes tree ASR
    nodelabels(num.anc2, cex=0.8, frame="none", adj = c(1.2, -0.3) )
    
    # Close the pdf file
    dev.off() 
    
    #use a phenogram to represent the phenotypes
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    
    # Open a pdf file
    pdf(paste(colnames(squirrel.data[i]),"plot2.pdf"),width = 18.27, height = 11.69, paper= "a4r")
    
    # 2. Create a plot
    phenogram(phy,dat3,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))
    
    # 2. Create the phenogram plot
    plotTree(phy,pts=F,node.numbers=T,fsize=0.6) #node labels

    # Close the pdf file
    dev.off() 
    
  }else {
    
    anc.eq<-fastAnc(tree_squirrel2,dat3,vars=TRUE,CI=TRUE) # estimate states  with variances & 95% confidence intervals for each node
    num.anc<-anc.eq$ace
    num.anc2<-round(num.anc,digits =2)
    
    write.table(num.anc2, paste(colnames(squirrel.data)[i],"_ancestral.csv"))
    
    #plot ancestral states in the tree
    
    #Represent ancestral states as colors 
    plt<-contMap(tree_squirrel2,dat3, plot=FALSE)
    plt<-setMap(plt, colors=c("black","darkmagenta","violetred2", "orange", "lightgoldenrod1"))

    # Open a pdf file
    pdf(paste(colnames(squirrel.data[i]),"plot.pdf"), width = 18.27, height = 11.69, paper= "a4r")
    
    # 2. Create a plot
    plot(plt,legend=FALSE,ylim=c(1-0.09*(Ntip(plt$tree)-1),Ntip(plt$tree)))
    
    #Bar for ancestral state range
    add.color.bar(leg=0.3*max(nodeHeights(tree_squirrel2)),plt$cols,colnames(squirrel.data[i]),
                  lims=plt$lims,digits=3,prompt=FALSE,x=20,
                  y=1-0.08*(Ntip(plt$tree)-60),lwd=7,fsize=0.8,subtitle="")
    
    #Bar time scale
    tree_squirrel2$root.time <- 53
    
    time.tree<- geoscalePhylo(tree=ladderize(tree_squirrel2,right=FALSE), units=c("Period", "Epoch", "Age"), boxes="Epoch",
                              cex.tip=0.8, cex.age=0.8, cex.ts=1.0, label.offset=0, x.lim=c(0,53), lwd=3, width=2)
    
    #Label for nodes tree ASR
    nodelabels(num.anc2, cex=0.8, frame="none", adj = c(1.2, -0.3) )

    # Close the pdf file
    dev.off()
    
    #use a phenogram to represent the phenotypes
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    
    # Open a pdf file
    pdf(paste(colnames(squirrel.data[i]),"plot2.pdf"),width = 18.27, height = 11.69,  paper= "a4r")
    
    # 2. Create a plot
    phenogram(tree_squirrel2,dat3,fsize=0.6,ftype="i",spread.labels=TRUE, spread.cost=c(1,0))
    
    # 2. Create the phenogram plot
    plotTree(tree_squirrel2,pts=F,node.numbers=T,fsize=0.6) #node labels
   
    # Close the pdf file
    dev.off()
  }
}

############################################################################################################
#Discrete characters: Locomotion 

#Import squirrel data 
dataloc<-read.csv("squirrels_PEQ_res.csv", header=T, row.names = 1)
loco<-setNames(dataloc[,24],rownames(dataloc))

plotTree(tree_squirrel2,type="phylogram",fsize=0.7,ftype="i",lwd=1)
cols<-setNames(c("chartreuse3", "chocolate4", "steelblue1","firebrick2","yellow"),levels(loco))
tiplabels(pie=to.matrix(loco[tree_squirrel2$tip.label],
                        levels(loco)),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=-2*par()$usr[1],
                  y=-40*par()$usr[3],fsize=0.8)

###### Reconstruction ##########

#MCMC approach -- sample character histories from their posterior probability distribution. 
mtrees<-make.simmap(tree_squirrel2,loco,nsim=1000, model="ER")

#summarize a set of stochastic maps
pd<-summary(mtrees)
pd
plot(pd,fsize=0.6,ftype="i",colors=cols,lwd=1,ylim=c(-2,Ntip(tree_squirrel2)))

add.simmap.legend(colors=cols,prompt=FALSE,x=5,
                  y=8,fsize=0.6)

#### END