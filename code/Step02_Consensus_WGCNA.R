### Converting PS19Array Data to Ensembl IDs; Using mm9, NCBIM37, dec 2011 and Affy from biomart

load('../output/Regressed_data_Heart.rda')
datExpr.Heart=as.data.frame(t(normExpr.reg))
targets.Heart=targets
rm(normExpr.reg,targets)

load('../output/Regressed_data_DRG.rda')
datExpr.DRG=as.data.frame(t(normExpr.reg))
targets.DRG=targets
rm(normExpr.reg,targets)

load('../output/Regressed_data_CBL.rda')
datExpr.CBL=as.data.frame(t(normExpr.reg))
targets.CBL=targets
rm(normExpr.reg,targets)



gnS=intersect(colnames(datExpr.Heart),intersect(colnames(datExpr.DRG),colnames(datExpr.CBL)))
length(gnS) ##16920

##Subsetting data
datExpr.Heart= datExpr.Heart[,match(gnS,colnames(datExpr.Heart))]
datExpr.DRG= datExpr.DRG[,match(gnS,colnames(datExpr.DRG))]
datExpr.CBL= datExpr.CBL[,match(gnS,colnames(datExpr.CBL))]

save(list=ls(),file="../output/ConsensusData.rda")






#######################Consensus WGCNA #############################################
nSets=3
setLabels=c("Heart","DRG","CBL")
shortLabels=setLabels

multiExpr=vector(mode="list",length=nSets)

multiExpr[[1]] = list(data=as.data.frame(datExpr.Heart)) # Heart
names(multiExpr[[1]]$data)=colnames(datExpr.Heart)
rownames(multiExpr[[1]]$data)=rownames(datExpr.Heart)

multiExpr[[2]] = list(data=as.data.frame(datExpr.DRG)) #DRG
names(multiExpr[[2]]$data)=colnames(datExpr.DRG)
rownames(multiExpr[[2]]$data)=rownames(datExpr.DRG)


multiExpr[[3]] = list(data=as.data.frame(datExpr.CBL)) #CBL
names(multiExpr[[3]]$data)=colnames(datExpr.CBL)
rownames(multiExpr[[3]]$data)=rownames(datExpr.CBL)


checkSets(multiExpr) # check data size


multiMeta=list(Heart =list(data=targets.Heart),DRG=list(data=targets.DRG),CBL=list(data=targets.CBL))

save(list=ls(),file="../output/ConsensusData.rda")


#######################Consensus  #############################################

## Network Construction

# Choose a set of soft-thresholding powers
powers = c(seq(1,10,by=1), seq(12,40, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
verbose = 5,networkType="signed",corFnc="bicor")[[2]]);


# Plot the results:
pdf("../output/1_Power.pdf", height=10, width=18)

colors = c("blue","red","green")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
"Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
for (col in 1:length(plotCols))
{
ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
}
}
# Plot the quantities in the chosen columns vs. the soft thresholding power

par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
if (set==1)
{
plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
main = colNames[col]);
addGrid();
}
if (col==1)
{
text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
labels=powers,cex=cex1,col=colors[set]);
} else
text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
labels=powers,cex=cex1,col=colors[set]);
if (col==1)
{
legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
} else
legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off()


save(list=ls(),file="../output/ConsensusData.rda")


###Actual Network
softPower=10


net=blockwiseConsensusModules(multiExpr, blocks = NULL,
                                         maxBlockSize = 30000, ## This should be set to a smaller size if the user has limited RAM
                                         randomSeed = 12345,
                                         corType = "pearson", ## no use for bicor
                                         power = softPower,
                                         consensusQuantile = 0.2,
                                         networkType = "signed",
                                         TOMType = "unsigned",
                                         TOMDenom = "min",
                                         scaleTOMs = TRUE, scaleQuantile = 0.8,
                                         sampleForScaling = TRUE, sampleForScalingFactor = 1000,
                                         useDiskCache = TRUE, chunkSize = NULL,
                                         deepSplit = 2,
                                         detectCutHeight = 0.995, minModuleSize = 100,
                                         mergeCutHeight = 0.2,
                                         saveTOMs = TRUE,
                                         consensusTOMFileNames = "ConsensusTOM-block.%b.rda")


save(list=ls(),file="../output/consensus.rda")
consMEs = net$multiMEs;
moduleColors = net$colors;
consTree = net$dendrograms[[1]];

pdf("../output/SignedDendro_Consensus.pdf",height=10, width=15)
plotDendroAndColors(consTree, moduleColors, "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
main = "Consensus gene dendrogram and module colors")
dev.off()

save(list=ls(),file="../output/consensus.rda")

load("ConsensusTOM-block.1.rda") # consensus TOM

# Various Tree Cutting Params
    consTree= hclust(1-consTomDS,method="average");

mColorh <- mLabelh <- colorLabels <- NULL
  for (minModSize in c(40,100)) {
    for (dthresh in c(0.2,0.25)) {
      for (ds in c(2,4)) {
          print("Trying parameters:")
          print(c(minModSize,dthresh,ds))
          tree = cutreeHybrid(dendro = consTree, pamStage=FALSE,
            minClusterSize = minModSize, cutHeight = 0.995,
            deepSplit = ds, distM = as.matrix(1-consTomDS))

          merged <- mergeCloseModules(exprData = multiExpr,colors = tree$labels,
                                      cutHeight = dthresh)
          mColorh <- cbind(mColorh,labels2colors(merged$colors))
          mLabelh <- c(mLabelh,paste("DS=",ds," mms=\n",minModSize," dcor=",dthresh))
        }
      }
    }
    save(list=ls(),file="../output/consensus.rda")

    pdf("../output/SignedDendro_Consensus.pdf",height=10, width=15)
    plotDendroAndColors(consTree,mColorh,groupLabels=mLabelh,addGuide=TRUE,dendroLabels=FALSE,main="Dendrogram With Different Module Cutting Parameters")
    dev.off()








### All Networks
tmpMulti =vector(mode="list",length=nSets)

		  for (x in 1:3){ ## Looping over Heart, DRG and CBL
		  	trgets <- multiMeta[[x]]$data
		  thisExpr <- multiExpr[[x]]$data

		   datTraits<- trgets[,c(12,6,11,14)]

		tmpMulti[[x]]$traitmat = as.data.frame(cbind(as.numeric(factor(datTraits[,1],c("ND","Dox","DoxR")))-1,as.factor(datTraits[,2]),as.factor(datTraits[,3]),as.numeric(datTraits[,4])))
		rownames(tmpMulti[[x]]$traitmat )=rownames(datTraits)
		colnames(tmpMulti[[x]]$traitmat )=c("Treatment","Group","Sex","RIN")


		geneSigs=matrix(NA,nrow=4,ncol=ncol(thisExpr)) # create a vector to hold the data

		for(i in 1:ncol(geneSigs)) {

      exprvec=as.numeric(thisExpr[,i]) # get the expression vector for ith gene
    	treatmentr=bicor(exprvec, tmpMulti[[x]]$traitmat[,1],use="pairwise.complete.obs")
    	groupr=sqrt(max(summary(lm(exprvec~as.factor(tmpMulti[[x]]$traitmat[,2])))$adj.r.squared,0))
    	sexr=sqrt(max(summary(lm(exprvec~as.factor(tmpMulti[[x]]$traitmat[,3])))$adj.r.squared,0)) # calculate adjusted R^2s square-root for categorical variables
    	rinr=bicor(exprvec, tmpMulti[[x]]$traitmat[,4],use="pairwise.complete.obs")

    	geneSigs[,i]=c(treatmentr, groupr, sexr,rinr)


			#cat('Done for gene...',i,'\n')
		}


    geneSigs[1,] =numbers2colors(as.numeric(geneSigs[1,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))
    geneSigs[2,] =numbers2colors(as.numeric(geneSigs[2,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1)) # For categorical variables
    geneSigs[3,] =numbers2colors(as.numeric(geneSigs[3,]),signed=FALSE,centered=FALSE,blueWhiteRed(100)[51:100],lim=c(0,1)) # For categorical variables like strain or wt_tg we do not want values, thus lim=c(0,1), and signed and centered=F
    geneSigs[4,] =numbers2colors(as.numeric(geneSigs[4,]),signed=TRUE,centered=TRUE,blueWhiteRed(100),lim=c(-1,1))

    rownames(geneSigs)=c("Treatment","Group","Sex","RIN")



		 tmpMulti[[x]]$genecols <- geneSigs

		  tmpMulti[[x]]$netData$netName <- c(paste("Signed bicor consensus network quantile of 0.2 and power of 10"))
		  tmpMulti[[x]]$netData$TOMdendrogram <- consTree
		  tmpMulti[[x]]$netData$moduleColors <- mColorh
		  tmpMulti[[x]]$netData$cutParameters <- mLabelh


		}



	###### Combine data	and plot final dendrogram




		  mColorh1 <- cbind(mColorh,t(tmpMulti[[1]]$genecols),t(tmpMulti[[2]]$genecols),t(tmpMulti[[3]]$genecols))
		  mLabelh1 <- c(mLabelh,rownames(tmpMulti[[1]]$genecols),rownames(tmpMulti[[2]]$genecols),rownames(tmpMulti[[3]]$genecols))

		  pdf("../output/Dendrogram_with_GeneSigs.pdf",height=30,width=25)
		  plotDendroAndColors(consTree, mColorh1, groupLabels = mLabelh1,addGuide=TRUE,dendroLabels=FALSE,main="FXN Consensus Mouse")
		  multiData.Heart <- tmpMulti[[1]]
		  multiData.DRG <- tmpMulti[[2]]
		  multiData.CBL <- tmpMulti[[3]]

		  save(multiData.Heart, multiData.DRG, multiData.CBL ,file="Individual_Consensus_data.rda")
		  dev.off()



save(list=ls(),file="../output/consensus.rda")




### Final redraw

load("../output/consensus.rda")
load("../output/Individual_Consensus_data.rda")
mms=100
ds=4
dthresh=0.2


tree = cutreeHybrid(dendro = consTree, pamStage=FALSE,
              minClusterSize = mms, cutHeight = 0.995,
              deepSplit = ds, distM = as.matrix(1-consTomDS))

         merged <- mergeCloseModules(exprData = multiExpr,colors = tree$labels,
                                        cutHeight = dthresh)
        mColorh <- cbind(labels2colors(merged$colors),t(multiData.Heart$genecols), t(multiData.DRG$genecols),t(multiData.CBL$genecols))
        mLabelh <- c("Merged Colors",rownames(multiData.Heart$genecols),rownames(multiData.DRG$genecols),rownames(multiData.CBL$genecols))

  pdf("../output/ConsensusTOM_FinalDendro.pdf",height=15,width=25)
  plotDendroAndColors(consTree, mColorh, groupLabels = mLabelh,addGuide=TRUE,dendroLabels=FALSE,main= paste("Signed bicor network with power =10, mms=",mms,"ds=",ds,"dthresh=",dthresh));
  dev.off()





moduleColors_Cons = labels2colors(merged$colors);

cols<-labels2colors(unique(merged$colors))
labels<-as.data.frame(cbind(unique(merged$colors),cols))
colnames(labels)<-c('Labels','colors')

#MEs_Heart=merged$newMEs[[1]]$data
#MEs_DRG=merged$newMEs[[2]]$data
#MEs_CBL=merged$newMEs[[3]]$data

consKME1=consensusKME(multiExpr=multiExpr, moduleColors_Cons,
                  multiEigengenes = NULL,
                  consensusQuantile = 0.2,
                  signed = TRUE)

 consensus.KMEs=consKME1[,regexpr('consensus.kME',names(consKME1))>0]

 MEs_DRG=moduleEigengenes(datExpr.DRG, colors = moduleColors_Cons,softPower= softPower, nPC=1)$eigengenes
 MEs_DRG=orderMEs(MEs_DRG)

 MEs_Heart=moduleEigengenes(datExpr.Heart, colors = moduleColors_Cons,softPower= softPower, nPC=1)$eigengenes
 MEs_Heart=orderMEs(MEs_Heart)

 MEs_CBL=moduleEigengenes(datExpr.CBL, colors = moduleColors_Cons,softPower= softPower, nPC=1)$eigengenes
 MEs_CBL=orderMEs(MEs_CBL)

 rm(consTomDS)
  rownames(consensus.KMEs)=colnames(multiExpr[[1]]$data)

ensembl=read.csv('../data/ensembl_GeneSymbol_Mouse.csv')
ensembl=ensembl[na.omit(match(rownames(consensus.KMEs),ensembl$Ensembl.Gene.ID)),]

geneInfo.cons=cbind(ensembl$Ensembl.Gene.ID,ensembl$Associated.Gene.Name,moduleColors_Cons,consensus.KMEs)
colnames(geneInfo.cons)[c(1:2)]<-c('Ensembl.Gene.ID','GeneSymbol')
save(list=ls(),file="../output/ConsensusKME_cquant0.2.rda")




##

group=factor(multiMeta[[1]]$data[,'Drug.x'])
group=factor(group,c('T0_WT_Dox','T0_TG_ND','T0_TG_Dox','T1_WT_Dox','T1_TG_ND','T1_TG_Dox','T3_WT_Dox','T3_TG_ND','T3_TG_Dox','T4_WT_Dox','T4_TG_ND','T4_TG_Dox','T5_WT_Dox','T5_TG_ND','T5_TG_Dox','T4_TG_DoxR','T5_TG_DoxR'))

uniquemodcolors=as.character(unique(geneInfo.cons$moduleColors_Cons))
pdf("ME_Heart.pdf",height=4,width=24)
for (j in 1:length(uniquemodcolors)){
	   thismod= uniquemodcolors[j]
	thisME <- MEs_Heart[,paste("ME",thismod,sep="")]
	boxplot(thisME~group,main=paste(thismod, "Module ME by Condition"), ylab="Module Eigengene Value")
}
dev.off()


###

group=factor(multiMeta[[2]]$data[,'Drug.x'])
group=factor(group,c('T0_WT_Dox','T0_TG_ND','T0_TG_Dox','T1_WT_Dox','T1_TG_ND','T1_TG_Dox','T3_WT_Dox','T3_TG_ND','T3_TG_Dox','T4_WT_Dox','T4_TG_ND','T4_TG_Dox','T5_WT_Dox','T5_TG_ND','T5_TG_Dox','T4_TG_DoxR','T5_TG_DoxR'))

uniquemodcolors=as.character(unique(geneInfo.cons$moduleColors_Cons))
pdf("ME_DRG.pdf",height=4,width=24)
for (j in 1:length(uniquemodcolors)){
	   thismod= uniquemodcolors[j]
	thisME <- MEs_DRG[,paste("ME",thismod,sep="")]
	boxplot(thisME~group,main=paste(thismod, "Module ME by Condition"), ylab="Module Eigengene Value")
}
dev.off()

###

group=factor(multiMeta[[3]]$data[,'Drug.x'])
group=factor(group,c('T0_WT_Dox','T0_TG_ND','T0_TG_Dox','T1_WT_Dox','T1_TG_ND','T1_TG_Dox','T3_WT_Dox','T3_TG_ND','T3_TG_Dox','T4_WT_Dox','T4_TG_ND','T4_TG_Dox','T5_WT_Dox','T5_TG_ND','T5_TG_Dox','T4_TG_DoxR','T5_TG_DoxR'))

uniquemodcolors=as.character(unique(geneInfo.cons$moduleColors_Cons))
pdf("ME_CBL.pdf",height=4,width=24)
for (j in 1:length(uniquemodcolors)){
	   thismod= uniquemodcolors[j]
	thisME <- MEs_CBL[,paste("ME",thismod,sep="")]
	boxplot(thisME~group,main=paste(thismod, "Module ME by Condition"), ylab="Module Eigengene Value")
}
dev.off()









#####GO Elite for Modules
dir.create("./geneInfo")
dir.create("./geneInfo/background/")
dir.create("./geneInfo/input/")
dir.create("./geneInfo/output/")


geneInfo.cons$SystemCode =rep("En",length=nrow(geneInfo.cons))
background=geneInfo.cons[,"Ensembl.Gene.ID"]
background=as.data.frame(background)

## Output files for GO elite
background <- cbind(background,rep("En",length=length(background)))
colnames(background) <- c("Source Identifier","SystemCode")
write.table(background,"./geneInfo/background/denominator.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")

uniquemodcolors=unique(moduleColors_Cons)
uniquemodcolors=uniquemodcolors[uniquemodcolors!='grey']
for(i in 1:length(uniquemodcolors)){
	thismod= uniquemodcolors[i]
	ind=which(colnames(geneInfo.cons)==paste("consensus.kME",thismod,sep=""))
	thisInfo=geneInfo.cons[geneInfo.cons$moduleColors_Cons==thismod,c(1,23,ind)] ##1=Ensembl.ID, 21="SystemCode",ind=kME value
	colnames(thisInfo) <- c("Source Identifier","SystemCode","kME")
	write.table(thisInfo,file=paste("./geneInfo/input/",thismod,"_Module.txt",sep=""),row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")


}



#Module Plots
geneInfo.vis=geneInfo.cons
nams=colnames(geneInfo.vis)[-c(1:3)]
colnames(geneInfo.vis)[-c(1:3)]<- substring(nams,14,30)

colnames(geneInfo.vis)[1]= "Ensembl.Gene.ID"
colnames(geneInfo.vis)[2]= "GeneSymbol"
colnames(geneInfo.vis)[3]= "Initially.Assigned.Module.Color"


library(igraph);
library(RColorBrewer);

  load("ConsensusTOM-block.1.rda")

TOM.matrix = as.matrix(consTomDS);



pdf("../output/ModuleNetworks.pdf",height=9,width=10);
# for (i in 1:length(sigmodcolors))  {
  # mod=sigmodcolors[i];
  # numgenesingraph = 50;
  # numconnections2keep = 1500;
for (mod in uniquemodcolors)  {
	# mod="turquoise"
  numgenesingraph = 100;
  numconnections2keep = 1500;
  cat('module:',mod,'\n');
  geneInfo.vis=geneInfo.vis[geneInfo.vis$GeneSymbol!="NA",]
  colind = which(colnames(geneInfo.vis)==mod);
  rowind = which(geneInfo.vis[,3]==mod);
  cat(' ',length(rowind),'probes in module\n');
  submatrix = geneInfo.vis[rowind,];
  orderind = order(submatrix[,colind],decreasing=TRUE);
  if (length(rowind) < numgenesingraph) {
    numgenesingraph = length(rowind);
    numconnections2keep = numgenesingraph * (numgenesingraph - 1);
  }
  cat('Making network graphs, using top',numgenesingraph,'probes and',numconnections2keep,'connections of TOM\n');
  submatrix = submatrix[orderind[1:numgenesingraph],];
  #Identify the columns in the TOM that correspond to these hub probes
  matchind = match(submatrix$Ensembl.Gene.ID,colnames(datExpr.Heart));
  reducedTOM = TOM.matrix[matchind,matchind];

  orderind = order(reducedTOM,decreasing=TRUE);
  connections2keep = orderind[1:numconnections2keep];
  reducedTOM = matrix(0,nrow(reducedTOM),ncol(reducedTOM));
  reducedTOM[connections2keep] = 1;

  g0 <- graph.adjacency(as.matrix(reducedTOM[1:10,1:10]),mode="undirected",weighted=TRUE,diag=FALSE)
    layoutMata <- layout.circle(g0)

    g0 <- graph.adjacency(as.matrix(reducedTOM[11:50,11:50]),mode="undirected",weighted=TRUE,diag=FALSE)
    layoutMatb <- layout.circle(g0)

   g0 <- graph.adjacency(as.matrix(reducedTOM[51:ncol(reducedTOM),51:ncol(reducedTOM)]),mode="undirected",weighted=TRUE,diag=FALSE)
   layoutMatc <- layout.circle(g0)
    g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
    layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.8, layoutMatc)

  plot(g1,edge.color="grey",vertex.color=mod,vertex.label=as.character(submatrix$GeneSymbol),vertex.label.cex=0.7,vertex.label.dist=0.45,vertex.label.degree=-pi/4,vertex.label.color="black",layout= layoutMat,vertex.size=submatrix[,colind]^2*8,main=paste(mod,"module"))

}
dev.off();
