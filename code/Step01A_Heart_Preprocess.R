
## Load normalized Illumina expression data for Heart
load('../data/Normalized_data_Heart.rda')

##CollapseRows


ensembl=read.csv("../data/illumina_ENSG_Mouse.csv")
ensembl=subset(ensembl,!duplicated(ensembl$Illumina.Mouse.Ref.8.V2.probe))
ensembl=ensembl[-1,]

ind1=intersect(ensembl$Illumina.Mouse.Ref.8.V2.probe,rownames(normExpr))
normExpr_ID= normExpr[na.omit(match(ind1,rownames(normExpr))),]
ensembl1 = ensembl[na.omit(match(ind1,ensembl$Illumina.Mouse.Ref.8.V2.probe)),]

normExpr_GeneID=collapseRows(normExpr_ID, as.character(ensembl1$Ensembl.Gene.ID), ensembl1$Illumina.Mouse.Ref.8.V2.probe)$datETcollapsed  ### changing to GeneSymbols using collapseRows

normExpr<-normExpr_GeneID

save(normExpr,targets,file="../output/ENSG_Normalized_data_Heart.rda")


#########PCA


pdf("../output/PCA_Heart.pdf",height=20,width=24)

	thisdat.HTSC <- t(scale(t(normExpr),scale=F))
	PC.HTSC <- prcomp(thisdat.HTSC,center=F);
	topPC1 <- PC.HTSC$rotation[,1:5];
	varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
	topvar <- varexp[1:5]
	colnames(topPC1) <- paste("Normalized.Expression\n",colnames(topPC1)," (",signif(100*topvar[1:5],2),"%)",sep="")

        pairsdat <- data.frame(treatment=as.numeric(factor(targets$Treatment)),Weight=as.numeric(targets$Body.weight),Sex=as.numeric(factor(targets$Sex)),TimePoint=as.numeric(factor(targets$TimePoint)),RIN=as.numeric(targets$RNA_RIN),Batch=as.numeric(factor(targets$Batch)),Extraction.Batch=as.numeric(factor(targets$Extraction.Batch)),Dissection.Date=as.numeric(factor(targets$Diss.date)))


		cond=labels2colors(targets$Treatment)  ## colors



		panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
		  usr <- par("usr"); on.exit(par(usr))
		  par(usr = c(0, 1, 0, 1))
		  r <- abs(cor(x, y,use="pairwise.complete.obs",method="pearson"))
		  txt <- format(c(r, 0.123456789), digits = digits)[1]
		  txt <- paste0(prefix, txt)
		  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
		  text(0.5, 0.5, txt, cex = cex.cor * r)
		}


	panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
		  usr <- par("usr"); on.exit(par(usr))
		  par(usr = c(0, 1, 0, 1))
		  if (class(x) == "numeric" & class(y) == "numeric") {
		    r <- abs(cor(x, y,use="pairwise.complete.obs",method="pearson"))
		  } else {
		    lmout <- lm(y~x)
		    r <- sqrt(summary(lmout)$adj.r.squared)
		  }
		  txt <- format(c(r, 0.123456789), digits = digits)[1]
		  txt <- paste0(prefix, txt)
		  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
		  text(0.5, 0.5, txt, cex = cex.cor * r)
	}

		pairs(cbind(topPC1,pairsdat),col= cond,pch=19,upper.panel = panel.cor,main="Covariates and Expression Comparison -- |Spearman's rho| correlation values")

dev.off()



###
sdout <- 4
  ## Remove outliers
  ##Calculate signed, weighted biweight midcorrelation
  normadj <- (0.5+0.5*bicor(normExpr)^2)

  ## Calculate connectivity
  netsummary <- fundamentalNetworkConcepts(normadj)
  ku <- netsummary$Connectivity
  z.ku <- ku-(mean(ku))/sqrt(var(ku))
  ## Declare as outliers those samples which are more than sdout sd above the mean connectivity based on the chosen measure
  outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
  print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep=""))
  print(colnames(normExpr)[outliers])
  print(table(outliers))

  normExpr <- normExpr[,!outliers]
  targets <- targets[!outliers,]


##Combat Date of dissection



library(sva);

par.prior=TRUE; ## TRUE means parametric adjustment will be used, false is non-parametric adjustment
prior.plots=FALSE;  ##Check with  parametric adjustment, if red and black lines are overlapping use parametric (fast), else non-parametric (slow)



	batch.all = targets
	batch = as.numeric(factor(targets$Diss.date));
	mod0 = model.matrix(~1,data=batch.all);

	datExpr.Combat = ComBat(dat=normExpr,batch=batch,mod=mod0,par.prior= par.prior,prior.plots= prior.plots);


	batch = as.numeric(factor(targets$Extraction.Batch));
	mod0 = model.matrix(~1,data=batch.all);

	datExpr.Combat = ComBat(dat=datExpr.Combat,batch=batch,mod=mod0,par.prior= par.prior,prior.plots= prior.plots);



pdf("../output/PCA_Heart_Combat.pdf",height=20,width=24)

	thisdat.HTSC <- t(scale(t(datExpr.Combat ),scale=F))
	PC.HTSC <- prcomp(thisdat.HTSC,center=F);
	topPC1 <- PC.HTSC$rotation[,1:5];
	varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
	topvar <- varexp[1:5]
	colnames(topPC1) <- paste("Normalized.Expression\n",colnames(topPC1)," (",signif(100*topvar[1:5],2),"%)",sep="")

        pairsdat <- data.frame(treatment=as.numeric(factor(targets$Treatment)),Weight=as.numeric(targets$Body.weight),Sex=as.numeric(factor(targets$Sex)),TimePoint=as.numeric(factor(targets$TimePoint)),RIN=as.numeric(targets$RNA_RIN),Batch=as.numeric(factor(targets$Batch)),Extraction.Batch=as.numeric(factor(targets$Extraction.Batch)),Dissection.Date=as.numeric(factor(targets$Diss.date)))


		cond=labels2colors(targets$Treatment)  ## colors



		panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
		  usr <- par("usr"); on.exit(par(usr))
		  par(usr = c(0, 1, 0, 1))
		  r <- abs(cor(x, y,use="pairwise.complete.obs",method="pearson"))
		  txt <- format(c(r, 0.123456789), digits = digits)[1]
		  txt <- paste0(prefix, txt)
		  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
		  text(0.5, 0.5, txt, cex = cex.cor * r)
		}


	panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
		  usr <- par("usr"); on.exit(par(usr))
		  par(usr = c(0, 1, 0, 1))
		  if (class(x) == "numeric" & class(y) == "numeric") {
		    r <- abs(cor(x, y,use="pairwise.complete.obs",method="pearson"))
		  } else {
		    lmout <- lm(y~x)
		    r <- sqrt(summary(lmout)$adj.r.squared)
		  }
		  txt <- format(c(r, 0.123456789), digits = digits)[1]
		  txt <- paste0(prefix, txt)
		  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
		  text(0.5, 0.5, txt, cex = cex.cor * r)
	}

		pairs(cbind(topPC1,pairsdat),col= cond,pch=19,upper.panel = panel.cor,main="Covariates and Expression Comparison -- |Spearman's rho| correlation values")

dev.off()



save(datExpr.Combat,targets,file="../output/CombatCorrected_Heart.rda")



########
####################Regression

options(stringsAsFactors=FALSE)
library(boot)
boot <- TRUE
numboot <- 1000
bs <- function(formula, data, indices) {
  d <- data[indices,] # allows bootstrap function to select samples
  fit <- lm(formula, data=d)
  return(coef(fit))
}

	normExpr=datExpr.Combat

	 Treatment<-as.numeric(factor(targets$Treatment,c("ND","Dox","DoxR")))-1
	Group<- as.numeric(factor(targets$Drug.x))-1
	Sex<-as.numeric(factor(targets$Sex))-1
	RIN<-as.numeric(targets$RNA_RIN)


	  varnames <- c("Treatment","Group","Sex","RIN")

regvars <- as.data.frame(cbind(Treatment,Group,Sex,RIN))


## Run the regression
normExpr.reg <- matrix(NA,nrow=nrow(normExpr),ncol=ncol(normExpr))
rownames(normExpr.reg) <- rownames(normExpr)
colnames(normExpr.reg) <- colnames(normExpr)
coefmat <- matrix(NA,nrow=nrow(normExpr),ncol=ncol(regvars)+1)## change it to ncol(regvars)+1 when condition has 2 levels

  set.seed(8675309)
  for (i in 1:nrow(normExpr)) {
    if (i%%1000 == 0) {print(i)}
    thisexp <- as.numeric(normExpr[i,])

   bs.results <- boot(data=data.frame(thisexp,regvars), statistic=bs,
                     R=numboot, formula=thisexp~Treatment+Group+Sex+RIN)

    ## get the median - we can sometimes get NA values here... so let's exclude these - old code #bs.stats <- apply(bs.results$t,2,median)
    bs.stats <- rep(NA,ncol(bs.results$t)) ##ncol is 3 here (thisexp, construct and extracted)
    for (n in 1:ncol(bs.results$t)) {
      bs.stats[n] <- median(na.omit(bs.results$t[,n]))
    }
    coefmat[i,] <- bs.stats
    normExpr.reg[i,] <- thisexp - bs.stats[5]*regvars[,"RIN"]
  }

quantile(normExpr[,1],c(0,0.025,0.25,0.5,0.75,0.975,1))
quantile(normExpr.reg[,1],c(0,0.025,0.25,0.5,0.75,0.975,1))


save(normExpr.reg,targets,file="../output/Regressed_data_Heart.rda")


library(WGCNA)
pdf("../output/PCA_Heart_RegressionRIN.pdf",height=20,width=24)

	thisdat.HTSC <- t(scale(t(normExpr.reg),scale=F))
	PC.HTSC <- prcomp(thisdat.HTSC,center=F);
	topPC1 <- PC.HTSC$rotation[,1:5];
	varexp <- (PC.HTSC$sdev)^2 / sum(PC.HTSC$sdev^2)
	topvar <- varexp[1:5]
	colnames(topPC1) <- paste("Normalized.Expression\n",colnames(topPC1)," (",signif(100*topvar[1:5],2),"%)",sep="")

        pairsdat <- data.frame(treatment=as.numeric(factor(targets$Treatment)),Weight=as.numeric(targets$Body.weight),Sex=as.numeric(factor(targets$Sex)),TimePoint=as.numeric(factor(targets$TimePoint)),RIN=as.numeric(targets$RNA_RIN),Batch=as.numeric(factor(targets$Batch)),Extraction.Batch=as.numeric(factor(targets$Extraction.Batch)),Dissection.Date=as.numeric(factor(targets$Diss.date)))


		cond=labels2colors(targets$Treatment)  ## colors



		panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
		  usr <- par("usr"); on.exit(par(usr))
		  par(usr = c(0, 1, 0, 1))
		  r <- abs(cor(x, y,use="pairwise.complete.obs",method="pearson"))
		  txt <- format(c(r, 0.123456789), digits = digits)[1]
		  txt <- paste0(prefix, txt)
		  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
		  text(0.5, 0.5, txt, cex = cex.cor * r)
		}


	panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { ## Useful function for comparing multivariate data
		  usr <- par("usr"); on.exit(par(usr))
		  par(usr = c(0, 1, 0, 1))
		  if (class(x) == "numeric" & class(y) == "numeric") {
		    r <- abs(cor(x, y,use="pairwise.complete.obs",method="pearson"))
		  } else {
		    lmout <- lm(y~x)
		    r <- sqrt(summary(lmout)$adj.r.squared)
		  }
		  txt <- format(c(r, 0.123456789), digits = digits)[1]
		  txt <- paste0(prefix, txt)
		  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
		  text(0.5, 0.5, txt, cex = cex.cor * r)
	}

		pairs(cbind(topPC1,pairsdat),col= cond,pch=19,upper.panel = panel.cor,main="Covariates and Expression Comparison -- |Spearman's rho| correlation values")

dev.off()
