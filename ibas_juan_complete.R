returnMinSec <- function(x) {
  result <- paste(floor(x/60)," Min and ", round(x%%60)," Sec",sep="")
  return (result)
}

computeProbDist <- function(IMMgrados,g,bootstrapSet,cvSet){
  
  #calculate the prob for each of the 2001 (assuming # of cross validation genes in config file=2000) cross validation genes
  hyperGeometricSample <- c()
  
  for(cvi in cvSet){
    
    #get the interaction partners of the cross validation node i
    neighbours <- IMMgrados[[cvi]]
    n <- length(neighbours)
    
    #create the reduced bootstrapSet set with the current bootstrap gene removed
    bootstrapSetMinusOne <- bootstrapSet[!bootstrapSet %in% g]
    
    kSet <- neighbours %in% bootstrapSetMinusOne
    
    #sum the number of TRUE elements which mean they were also in the boostrap
    k <- sum(kSet)
    
    #the size of the boostrap (subtract one so as to not count the held out gene)
    K <- length(bootstrapSet) - 1
    
    #compute the probability
    hyperGeometricProb <- dhyper(k, #white balls drawn from urn
                                 K, #white balls in urn
                                 ReferenceNetworkSize - K, #black balls in urn
                                 n) #total balls drawn from urn
    
    #add hypergeometric probability to the set
    hyperGeometricSample <- c(hyperGeometricSample,hyperGeometricProb )
  }
  
  return (hyperGeometricSample)
}

lookupLikelihoodRatio <- function(prob, roc){
  
  #so here is our candidate list of ROC points
  potentialROCPpoints <- roc[roc$tpr>0 & roc$fpr>0,]
  
  #find the differences and min distance between each point and the input probability 
  potentialROCPpoints$difference <- abs(potentialROCPpoints$rocSpace - prob)
  minDifference <- min(potentialROCPpoints$difference)
  
  #use the min distance to extract the corresponding likelihood ratio
  predictedROCPoint <- potentialROCPpoints[potentialROCPpoints$difference == minDifference, ]
  
  if(nrow(predictedROCPoint)>1){
    predictedROCPoint <- predictedROCPoint[1,]
  }
  return (predictedROCPoint$tpr/predictedROCPoint$fpr)
  
}

# Score the entire reference networks using the provided gene and its optimal ROC curve

scoreReferenceNetwork <-function(gene,cell,optimalROCCurve, bootstrapGenes){
  scores <- c()
  for (genesInRefNet in V(inweb.g)$name){
    # get interaction partners based on the cell
    interactionPartners <- grados.IMM[[cell]][[genesInRefNet]]
    n <- length(interactionPartners)
    
    #create the reduced bootstrapGenes set with the current gene removed   
    bootstrapGenesMinusOne <- bootstrapGenes[!bootstrapGenes %in% gene]
    
    
    kSet <- interactionPartners %in% bootstrapGenesMinusOne
    
    #sum the number of TRUE elements which mean they were also in the boostrap
    k <- sum(kSet)
    
    #the size of the boostrap (subtract one so as to not count the held out gene)
    K <- length(bootstrapGenes) - 1
    
    hyperGeometricProb <- dhyper(k, #white balls drawn from urn
                                 K, #white balls in urn
                                 ReferenceNetworkSize - K, #black balls in urn
                                 n) #total balls drawn from urn
    
    likelihoodRatio <- lookupLikelihoodRatio(hyperGeometricProb,optimalROCCurve)
    
    scores <- c(scores, likelihoodRatio)
    
  }
  return (scores)
}

# ROC curves
# Description: This function takes in a list of probabilities from the cross validation set (similar genes plus the holdout appended to the end)
#              It iterates through the ROCCutoffs specified in the config file and creates the true positive rates and false positive rates by cutoff value
# Parameters:
#        probSet: list of enrichment probabilities 
#        measure: an indicator for saving the roc curve 
#        outFile: the location to save the roc curve

createROCCurve <- function(probSet,measure,outFile){
  
  tpRates  <- c()
  fpRates  <- c()
  numCVGenes <- length(probSet)
  
  savetp  <- c()
  savep   <- c()
  savefp  <- c()
  saven   <- c()
  fpRates <- c()
  tpRates <- c()
  
  useThisROC <- c()
  if(useOriginalIbasROCCurve){
    useThisROC <- c(1,0.01)
    lowerBound <- min(probSet)
    insertNewValue <- 0.01/5
    while(insertNewValue>lowerBound){
      useThisROC <- c(useThisROC,insertNewValue)
      insertNewValue <- insertNewValue/5
    }
    useThisROC <- sort(useThisROC,decreasing = FALSE)
  }else{
    useThisROC <- customROCCutoffs
  }
  
  for(cutoff in useThisROC){
    
    p = 1
    n = numCVGenes - 1
    
    tp = 0
    fp = 0
    
    position <- 1
    
    #iterate through each of the enrichment probabilities
    for (prob in probSet){
      
      #I always append the holdout gene to the last element so I only need to check the 
      #last element for the true positive
      if(position == numCVGenes){
        if( prob <= cutoff ){
          tp <- tp + 1
        }
      }else{
        #these are the real cross-validation genes
        if( prob <= cutoff ){
          fp <- fp + 1
        }
      }
      position <- position + 1
    }
    
    tpRate  <- tp/p
    fpRate  <- fp/n
    tpRates <- c(tpRates,tpRate)  
    fpRates <- c(fpRates,fpRate)  
    
    savetp  <- c(savetp,tp) 
    savep   <- c(savep,p)
    savefp  <- c(savefp,fp)
    saven   <- c(saven,n)  
  }
  
  rocCurve <- data.frame(tpr = tpRates,fpr = fpRates, rocSpace = useThisROC)
  
  #save the results if requested
  if(measure){
    df <- data.frame( cutoff=ROCCutoffs, truePositives=savetp, positives=savep, falsePositives=savefp,negatives=saven,falsePositiveRate=fpRates,truePositiveRates=tpRates )
    write.table(df, file = outFile, sep = ",",row.names = FALSE) 
  }
  
  return (rocCurve)
}

# Name: outputROCsAndAUCs
# Description: This function takes in a label and cross validation sample (of probabilities)
#              and outputs the ROC curves, ROC curve data, AUC distribution and AUC data
# Parameters:
#             label:  descriptive name of the test network
#          CVSample:  a list of lists containing the cross validation probabilities for each bootstrap gene

outputROCsAndAUCs <- function(label,CVSample){
  
  index <- 1
  outputAUCData <- c() #save all of the AUC's in this list
  rocOutput <- paste(outputLocation,label,"/ROCs/",sep="")
  aucOutput <- paste(outputLocation,label,"/AUCs/",sep="")
  
  for(l in CVSample){
    
    #get the probability valuesfrom this list of lists
    probValues <- unlist(l)
    
    #location to put the actual roc values
    rocOutDataLocation <- paste(outputLocation,label,"/ROCs/data/roc",index,".csv",sep="")
    
    #create and return the ROC curve as a data.frame
    roc <- createROCCurve(probValues,TRUE,rocOutDataLocation)
    
    #plot the ROC curve and output the result
    png(paste(rocOutput,label,"_CV_ROC", as.character(index), ".png",sep=""))
    plot(roc[,2], roc[,1],main = paste("ROC Curve for Gene",index," in Bootstrap of ", label,sep=""),xlab="False Positive %",ylab="True Positive %")
    lines(roc[,2], roc[,1])
    polygon( c(min(roc[,2]), roc[,2], max(roc[,2])), c( min(roc[,1]), roc[,1], min(roc[,1])), density=20 ) 
    dev.off()
    
    #now compute the AUC and save it to the AUC list
    auc <- AUC(roc)
    outputAUCData <- c(outputAUCData, auc)
    
    index <- index + 1
  }
  
  #create and output the AUC distribution (AUC's for all of the genes in this bootstrap)
  df <- data.frame(probs=outputAUCData)
  names(df)[1] <- "aucData"
  m <- ggplot(df, aes(x=aucData))
  m + geom_histogram(aes(fill = ..count..), binwidth=0.025,origin=0) +
    scale_fill_gradient("Count", low = "lightgreen", high = "red") + 
    ggtitle(paste("Distribution of AUC's for Genes in the Bootstrap of ",label,sep="")) + 
    xlab("Area Under The Curve") +
    ylab("Count of Bootstrap Genes") 
  ggsave(paste(aucOutput,"AUCDist",".png",sep=""))
  
  write.table(df, file = paste(aucOutput,"AUCDist.csv",sep=""), sep = ",",row.names = FALSE)
  
}

outputROCsAndAUCs <- function(label,CVSample){
  
  index <- 1
  outputAUCData <- c() #save all of the AUC's in this list
  rocOutput <- paste(outputLocation,label,"/ROCs/",sep="")
  aucOutput <- paste(outputLocation,label,"/AUCs/",sep="")
  
  for(l in CVSample){
    
    #get the probability valuesfrom this list of lists
    probValues <- unlist(l)
    
    #location to put the actual roc values
    rocOutDataLocation <- paste(outputLocation,label,"/ROCs/data/roc",index,".csv",sep="")
    
    #create and return the ROC curve as a data.frame
    roc <- createROCCurve(probValues,TRUE,rocOutDataLocation)
    
    #plot the ROC curve and output the result
    png(paste(rocOutput,label,"_CV_ROC", as.character(index), ".png",sep=""))
    plot(roc[,2], roc[,1],main = paste("ROC Curve for Gene",index," in Bootstrap of ", label,sep=""),xlab="False Positive %",ylab="True Positive %")
    lines(roc[,2], roc[,1])
    polygon( c(min(roc[,2]), roc[,2], max(roc[,2])), c( min(roc[,1]), roc[,1], min(roc[,1])), density=20 ) 
    dev.off()
    
    #now compute the AUC and save it to the AUC list
    auc <- AUC(roc)
    outputAUCData <- c(outputAUCData, auc)
    
    index <- index + 1
  }
  
  #create and output the AUC distribution (AUC's for all of the genes in this bootstrap)
  df <- data.frame(probs=outputAUCData)
  names(df)[1] <- "aucData"
  m <- ggplot(df, aes(x=aucData))
  m + geom_histogram(aes(fill = ..count..), binwidth=0.025,origin=0) +
    scale_fill_gradient("Count", low = "lightgreen", high = "red") + 
    ggtitle(paste("Distribution of AUC's for Genes in the Bootstrap of ",label,sep="")) + 
    xlab("Area Under The Curve") +
    ylab("Count of Bootstrap Genes") 
  ggsave(paste(aucOutput,"AUCDist",".png",sep=""))
  
  write.table(df, file = paste(aucOutput,"AUCDist.csv",sep=""), sep = ",",row.names = FALSE)
  
}


# Name : AUC
# Description: This function takes in a rocCurve data.frame and compute the area under the curve 
#
# Parameters:
#       rocCurve: a dataframe with the true positive rate and false positive rate
#       

AUC <- function(rocCurve){
  
  ptpr   <- 0
  pfpr   <- 0
  aucSum <- 0
  
  for(row in 1:nrow(rocCurve)){
    
    if(row==1){
      ptpr <- 0
      pfpr <- 0
    }else{
      ptpr <- rocCurve[row-1,1]
      pfpr <- rocCurve[row-1,2]
    }
    
    tpr <- rocCurve[row,1]
    fpr <- rocCurve[row,2]
    
    #article: from kasper's email https://ccrma.stanford.edu/workshops/mir2009/references/ROCintro.pdf
    aucSum <- aucSum + (tpr*(fpr-pfpr) - 0.5 * (tpr - ptpr) * (fpr - pfpr))
    
  }
  
  return (aucSum)
}

# Name : outputAUCs
# Description: This function takes in a list of bootstrap gene's cross validation probability distributions and creates a ROC curve and then AUC for each
#              It then outputs all AUC's to a file and plots the distribution of AUC's
# Parameters:
#       rocCurve: a dataframe with the true positive rate and false positive rate
#      

outputAUCs <- function(label,CVSample){
  
  index <- 1
  output <- paste(outputLocation,label,"/AUCs/",sep="")
  outputAUCData <- c()
  
  for(l in CVSample){
    probValues <- unlist(l)
    roc <- createROCCurve(probValues,FALSE, "")
    auc <- AUC(roc)
    outputAUCData <- c(outputAUCData, auc)
    index <- index + 1
  }
  
  df <- data.frame(probs=outputAUCData)
  names(df)[1] <- "aucData"
  m <- ggplot(df, aes(x=aucData))
  m + geom_histogram(aes(fill = ..count..), binwidth=0.025,origin=0) +
    scale_fill_gradient("Count", low = "lightgreen", high = "red") + 
    ggtitle(paste("Distribution of AUC's for Genes in the Bootstrap of ",label,sep="")) + 
    xlab("Area Under The Curve") +
    ylab("Count of Bootstrap Genes") 
  ggsave(paste(output,"AUCDist",".png",sep=""))
  
  write.table(df, file = paste(output,"AUCDist.csv",sep=""), sep = ",",row.names = FALSE)
}


################################################################################
#############################    IBAS      #####################################
################################################################################

### Reading inweb3 data
read.table("/Users/jfertaj/INVESTIGACION/WELLCOME/Data_Integration/IBAS_Juan/InWeb3/InWeb3_score.eda", skip=1) -> inweb3
inweb3 <- inweb3[,c(1,3,5)]
inweb.g <- graph.data.frame(inweb3, directed=F)

#######################
######  Specs   #######
#######################

debug <- FALSE  #you can turn print statements for debugging on/off
saveData <- TRUE # you can save data
notesToThisRun <-  "Testing to see about IBAS running in r"
numCVgenes <- 400  #number of cross validation genes
numBootstraps <- 10  #number of boostraps (sampled networks) we wish to use
useOriginalIbasROCCurve <- TRUE
ReferenceNetworkSize <- vcount(inweb.g)

#######################
######  Outputs #######
#######################
#this list will contain the average max AUC from the bootstrap genes for each of the bootstraps
aveMaxAUCsForAllBootstraps <- c()
#if save is on we will save the auc by cell to see if there is any systematic behavior 
aucByCell <- data.frame(ThisCell= character(0), AUCResult= numeric(0))
# this data.frame will contain columns that are the average posterior from the optimal cells in the bootstrap
scoredNetworkForAllBootstraps <-  data.frame(init=replicate(ReferenceNetworkSize,0))
# this data.frame will contain colums that are the max LR from each bootstrap 
maxLRForAllBootstraps <- data.frame(init=replicate(ReferenceNetworkSize,0))

numberOfROCPoints <- c()
allPercentilesOfHeldoutGeneProbs <-c()

##########################
######  Algorithm  #######
##########################

setwd("~/INVESTIGACION/WELLCOME/Data_Integration/IBAS_Juan")

WNTNodeNumbers<-c('11961','1236','13075','13308','13971','16377','22173','22730','329','3687','7590')
TGFNodeNumbers <- c('11388', '12424','16681','17783','22488','22780','3803','4554', '6216', '8964', '9049')
RR3KNodeNumbers <- c('10090', '1064','11032', '12222', '1635', '4695', '4792', '5746', '7622', '8226','857')

#init the testNetwork to WNT for now
testNetwork <- WNTNodeNumbers

networkTime <- proc.time()

ISs <- c(0, 0.05, 0.1, 0.154, 0.2, 0.25, 0.3, 0.35)
IMM.subg <- lapply(ISs, function(x){subg <- subgraph.edges(inweb.g, which(E(inweb.g)$V5>=x))})

#IMM is a list of pairs containing the c(network order, interaction threshold)
IMM <- c()
Orders <- seq(1,2,1) #network order
maxOrder <- max(Orders) 

for(o in Orders){
  for(is in ISs){
    IMM <- c(IMM,list(c(o,is)))
  }
}

immCellToKey <- function(cell){
  return(paste(as.character(unlist(cell)[1]),as.character(unlist(cell)[2]),sep="-"))
}

grados.IMM <- list()
for(thres in 1:length(IMM.subg)){
  #while(thres < length(IMM.subg){
  if(thres==1){
    vecinos.1 <- neighborhood(IMM.subg[[thres]], 1, V(IMM.subg[[thres]]))
    vecinos.1 <- lapply(vecinos.1, function(x) V(IMM.subg[[thres]])$name[x])
    names(vecinos.1) <- V(IMM.subg[[thres]])$name ### cambiando esto
    vecinos.2 <- neighborhood(IMM.subg[[thres]], 2, V(IMM.subg[[thres]]))
    vecinos.2 <- lapply(vecinos.2, function(x) V(IMM.subg[[thres]])$name[x])
    names(vecinos.2) <- V(IMM.subg[[thres]])$name
    for (h in 1:length(vecinos.2)){
      vecinos.2[[h]] <- vecinos.2[[h]][!vecinos.2[[h]] %in% vecinos.1[[h]]]
    }
    grados.IMM[[1]] <- vecinos.1
    grados.IMM[[2]] <- vecinos.2
  }else{
    j <- thres*2
    vecinos.1 <- neighborhood(IMM.subg[[thres]], 1, V(IMM.subg[[thres]]))
    vecinos.1 <- lapply(vecinos.1, function(x) V(IMM.subg[[thres]])$name[x])
    names(vecinos.1) <- V(IMM.subg[[thres]])$name
    vecinos.2 <- neighborhood(IMM.subg[[thres]], 2, V(IMM.subg[[thres]]))
    vecinos.2 <- lapply(vecinos.2, function(x) V(IMM.subg[[thres]])$name[x])
    names(vecinos.2) <- V(IMM.subg[[thres]])$name
    for (h in 1:length(vecinos.2)){
      vecinos.2[[h]] <- vecinos.2[[h]][!vecinos.2[[h]] %in% vecinos.1[[h]]]
    }
    grados.IMM[[j-1]] <- vecinos.1
    grados.IMM[[j]] <- vecinos.2
  }
}
names(grados.IMM) <- unlist(lapply(IMM, immCellToKey))

grados.inweb <- data.frame(grado=degree(inweb.g), genes=names(degree(inweb.g)))
inweb.grados <- split(grados.inweb, grados.inweb[,1])
inweb.grados <- lapply(inweb.grados, '[[', 2)
inweb.grados <- lapply(inweb.grados, as.character)

l3 <- list()
ll1 <- length(inweb.grados)
length(l3) <- ll1
for (i in 1:ll1){
  vec1 <- inweb.grados[[i]]
  jl <- 1;jr<-1;
  while (length(vec1) < 2000){
    if(i==1 || i-jl==0) {
      vec1 <- c(vec1, inweb.grados[[i+jr]])
      jr <- jr+1
    } else if (i==ll1 || jr+i==ll1 ){
      vec1 <- c(vec1, inweb.grados[[i-jl]])
      jl <- jl+1
    }else {
      vec1 <- c(vec1, inweb.grados[[i-jl]], inweb.grados[[i+jr]])
      jl <- jl+1
      jr <- jr+1
    } 
  } 
  #l3[[i]] <- sample2(vec1, 10) 
  l3[[i]] <- vec1  
}

names(l3) <- names(table(degree(inweb.g)))

networkTimeString <- returnMinSec((proc.time() - networkTime)[3])
print(paste("Network preparation took: ", networkTimeString, sep=""))

getBootstrapSampleGenes <- function(nodes){
  sample(names(which(degree(inweb.g) >= (1-0.15)*nodes & degree(inweb.g) <= (1+0.15)*nodes)),1)
} # we need to replace 0.15 for a variable call thres

ibasStartTime <- proc.time()
for(i in 1:numBootstraps){
  degreesForTestNetwork <- degree(inweb.g, testNetwork)
  if(i==1){
    bootstrapGenes <- testNetwork
  }else{
    bootstrapGenes <-sapply(degreesForTestNetwork, getBootstrapSampleGenes)
  }	
  
  if (debug) print(paste("STARTING BOOTSTRAP...", i, sep=" "))
  #bootstrapGenes <- testNetwork
  #degreesForBootstrap <- degree(inweb.g, bootstrapGenes) ### move to inside saveData
  
  #bootstrapGenes <-sapply(degreesForBootstrap, getBootstrapSampleGenes)
  
  #if save is on then save the degrees for each of the bootstrap genes
  if (saveData){
    degreesForBootstrap <- degree(inweb.g, bootstrapGenes) 
    bootstrapColumn <- replicate(length(bootstrapGenes),i)
    df <- data.frame(degree=degreesForBootstrap, bootstrap=bootstrapColumn)
    write.table(df, file = paste(outputLocation,"SavedData/","degreeByBootstraps.txt",sep=""), 
                append = TRUE, sep = "\t", eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names=FALSE)
  } 
  
  if (debug) print("Successfully sampled bootstrapGenes...")
  #print("Successfully sampled bootstrapGenes...") # change to debug
  
  #initialize allScoredNetwork to hold all of the scored reference networks 
  allScoredNetwork <- data.frame(init=replicate(ReferenceNetworkSize,0))
  averageMaxAUC <- 0
  for(gene in bootstrapGenes){
    
    if (debug) print(paste("For bootstrap gene ", gene, sep=""))
    #print(paste("For bootstrap gene ", gene, sep="")) # change to debug
    grado <- degree(inweb.g, gene)
    CVset <- l3[[grado]]
    CVset <- CVset[!CVset %in% gene]
    CVset <- c(sample(l3[[grado]], numCVgenes, replace=FALSE), gene)
    
    if (debug) print("Successfully sampled CVSet...")
    #print("Successfully sampled CVSet...") # change to debug
    
    #for each cell in the IMM
    maxAUC <- 0
    optimalCell <- ""
    optimalROCCurve <- data.frame()
    for (cell in names(grados.IMM)){
      
      #prepare the cell
      #keyForIMMHash <- immCellToKey(cell)
      if (debug) print(paste("For cell ", cell, "...", sep=""))
      #print(paste("For cell ", cell, "...", sep="")) # change to debug
      probValues <- computeProbDist(grados.IMM, gene,bootstrapGenes, CVset)
      if (debug) print(paste("probValues completed successfully ...", sep=""))
      #print(paste("probValues completed successfully ...", sep="")) # change to debug
      
      probOfHeldoutGene <- probValues[length(probValues)]
      sortedProbs <- sort(probValues, decreasing=TRUE)
      #percentage of bootstrap p-values larger than that of test-network
      percentileOfHeldoutGeneProb <- (match(probOfHeldoutGene, sortedProbs)-1)/(length(sortedProbs)-1)
      #allPercentilesOfHeldoutGeneProbs <- c(allPercentilesOfHeldoutGeneProbs,percentileOfHeldoutGeneProb)
      
      #if save is on then save the degrees for each of the bootstrap genes
      if (saveData){
        ITMColumn <- replicate(length(probValues),cell)
        df <- data.frame(cell=ITMColumn, probs=probValues)
        write.table(df, file = paste(outputLocation,"SavedData/","probabilitiesByITMCells.csv",sep=""), 
                    append = TRUE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names=FALSE)
      }
      
      roc <- createROCCurve(probValues,FALSE,"")
      numROCPoints <- nrow(roc)
      numberOfROCPoints <- c()
      numberOfROCPoints <- c(numberOfROCPoints,nrow(roc))
      
      if (saveData){
        ITMColumn <- replicate(1,cell)
        df <- data.frame(bootstrap=i,cell=ITMColumn, nrowROC=numROCPoints,percOfHeldoutGeneProb=percentileOfHeldoutGeneProb)
        write.table(df, file = paste(outputLocation,"SavedData/","nrowROCAndProbPercByITMCells.csv",sep=""), 
                    append = TRUE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names=FALSE)
        
        #Save the Viz of the ROC curve for this bootstrap, gene, cell
        saveAndVizROC(roc,i, gene, keyForIMMHash,numROCPoints,percentileOfHeldoutGeneProb)
      }
      
      #now compute the AUC and save it to the AUC list
      auc <- AUC(roc)
      
      if (saveData){
        aucByCell <- rbind(aucByCell, data.frame(ThisCell = keyForIMMHash, AUCResult = auc))
      }
      
      #update the max AUC/cell/optimal ROC curve if necessary
      if(auc > maxAUC){
        maxAUC <- auc
        optimalCell <- cell
        optimalROCCurve <- roc
      }
      if (debug) print("Successfully computed ROC/AUC")
    }
    
    averageMaxAUC <- averageMaxAUC + maxAUC
    
    #we have the optimal cell and max AUC for this gene so we want to save it to a master file
    #score the entire inweb network with the optimal cell for this gene
    scoredNetwork <- scoreReferenceNetwork(gene,optimalCell,optimalROCCurve, bootstrapGenes)
    if (debug) print(paste("For bootstrap gene ", gene, " also successfully scored reference network", sep=""))
    allScoredNetwork <- cbind(allScoredNetwork, scoredNetwork)     
  }
  
  #take the average of the max auc's and save it for the test network p-value calculation
  averageMaxAUC <- averageMaxAUC / length(bootstrapGenes)
  aveMaxAUCsForAllBootstraps <- c(aveMaxAUCsForAllBootstraps,averageMaxAUC)  
  
  #drop the initial column from the 
  allScoredNetwork <- data.frame(allScoredNetwork[,!(names(allScoredNetwork) %in% c("init"))])
  
  #average the optimal and average etc posterior and AUC across bootstrap 
  #genes and save these because we will use them in the p-value calculations
  #transform this list of lists into a dataframe or just make a dataframe in the beginning
  
  allScoredNetwork$aveLikelihood = apply(allScoredNetwork,1,mean,na.rm=TRUE)
  
  #get the max likelihood ratio for the bootstrap
  maxLR <- apply(allScoredNetwork,1,max,na.rm=TRUE)
  maxLRForAllBootstraps <- cbind(maxLRForAllBootstraps,maxLR)
  
  # now add allScoredNetwork$aveLikelihood to some global data.frame scoredNetworkForAllBootstraps
  scoredNetworkForAllBootstraps <- cbind(scoredNetworkForAllBootstraps,allScoredNetwork$aveLikelihood)
  print(paste("BOOTSTRAP: ", i, " COMPLETED.", sep=""))
  
}
runTimeString <- returnMinSec((proc.time() - ibasStartTime)[3])
print(paste("num of BOOTSTRAPS = ", numBootstraps, " with ", numCVgenes,  " took: ", runTimeString, sep=""))

if (saveData){
  
  #########################################################################
  # First create and save the degree distribution histograms (test vs Other)
  #########################################################################
  source(paste(codeLocation, "visualize/","vizDegreeDistByBootstrap.r", sep=""))
  
  #########################################################################
  # Save the numberOfROCPoints compiled during the algorithm
  #########################################################################
  write.table(numberOfROCPoints, 
              file = paste(outputLocation,"SavedData/","numberOfROCPoints.csv",sep=""), 
              append = FALSE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = FALSE)
  
  #########################################################################
  # Save the aveMaxAUCsForAllBootstraps and aucByCell compiled
  # during the algorithm
  #########################################################################
  # save aveMaxAUCsForAllBootstraps to a file for viz 
  write.table(aveMaxAUCsForAllBootstraps, 
              file = paste(outputLocation,"SavedData/","aveMaxAUCsForAllBootstraps.csv",sep=""), 
              append = FALSE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = FALSE)
  
  # save aucByCell to a file for viz 
  write.table(aucByCell, 
              file = paste(outputLocation,"SavedData/","aucByCell.csv",sep=""), 
              append = FALSE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names = FALSE)
  
  ############################################################################
  # Save the scoredNetworkForAllBootstraps which has the entire reference
  # network scored (average likelihoods from optimal cells) for each bootstrap
  ############################################################################
  
  # first remove the init column (superfluous)
  scoredNetworkForAllBootstraps <- data.frame(scoredNetworkForAllBootstraps[,!(names(scoredNetworkForAllBootstraps) %in% c("init"))])
  
  #assign the column names the bootstrap iteration "Bootstrap i"
  colnames(scoredNetworkForAllBootstraps) <- makeVarNamesList("Bootstrap ",numBootstraps)
  
  write.table(scoredNetworkForAllBootstraps, 
              file = paste(outputLocation,"SavedData/","scoredNetworkForAllBootstraps.csv",sep=""), 
              append = FALSE, sep = ",", eol = "\n", na = "NA", dec = ".", row.names = FALSE,col.names=TRUE)
  
  #compare the average likelihood across the bootstraps to the likelihood from the test network (ave of the test network)
  saveAndVizReferenceNetworkScores(scoredNetworkForAllBootstraps,"Compare Ave Likelihood: Test vs Ave of Other Bootstraps","Ave Likelihood","Gene Count","CompareAveLikelihoodScores-TestVSOtherBSs")
  
  ############################################################################
  # Gene Specific P-Value: Compare the Average Likelihood from the Test
  # Network to the Average Likelihoods from the other bootstraps
  ############################################################################
  #calculate gene specific p-value, returned as an array
  geneSpecificPV <- c()
  for (i in 1:ReferenceNetworkSize){
    sortedScores <- sort(scoredNetworkForAllBootstraps[i,], decreasing=TRUE)
    PV <- (match(scoredNetworkForAllBootstraps[i,][[1]],sortedScores)-1)/length(sortedScores)
    geneSpecificPV <- c(geneSpecificPV, PV)
  }
  
  #save and visualize gene specific p-value
  saveAndVizHistogram(geneSpecificPV,
                      title="Gene Specific P-Value Distribution",
                      xAxisLabel="P-Values",
                      yAxisLabel="Count of Genes",
                      "GeneSpecificPValues")
  
  ############################################################################
  # Gene Based P-Value: Compare the Max Likelihood from the Test
  # Network to the Max Likelihoods from the other bootstraps
  ############################################################################
  #first remove the init column (superfluous)
  maxLRForAllBootstraps <- data.frame(maxLRForAllBootstraps[,!(names(maxLRForAllBootstraps) %in% c("init"))])
  #calculate gene-based p-value, returned as an array
  
  #compare the average MAX likelihood across the bootstraps to the max likelihood from the test network
  saveAndVizReferenceNetworkScores(maxLRForAllBootstraps,"Compare MAX Likelihood: Test vs Ave of Other Bootstraps","Max Likelihood","Gene Count","CompareMAXLikelihoodScores-TestVSOtherBSs")
  
  geneBasedPV <- c()
  for (i in 1:ReferenceNetworkSize){
    sortedScores <- sort(maxLRForAllBootstraps[i,], decreasing=TRUE)
    PV <- (match(maxLRForAllBootstraps[i,][[1]],sortedScores)-1)/length(sortedScores)
    geneBasedPV <- c(geneBasedPV, PV) 
  }
  
  #save and visualize gene based p-value
  saveAndVizHistogram(geneBasedPV,
                      title="Gene Based P-Value Distribution",
                      xAxisLabel="P-Values",
                      yAxisLabel="Count of Genes",
                      "GeneBasedPValues")
  
  ############################################################################
  # Visualize the distribution of enrichment probabilities by IMM cell
  ############################################################################
  #during the alg we saved all probs to probabilitiesByITMCells.csv, this script reads and makes visual for the probs
  source(paste(vizCodesLocation, "vizProbabilitiesByCell.r", sep=""))
  
  ############################################################################
  # Calculate and Visualize the Test Network P-Value 
  ############################################################################
  #calculate test-network p-value
  testNetworkMaxAUC <- aveMaxAUCsForAllBootstraps[1]
  sortedAveMaxAUC <- sort(aveMaxAUCsForAllBootstraps, decreasing=TRUE)
  #percentage of bootstrap p-values larger than that of test-network
  testNetworkPV <- (match(testNetworkMaxAUC, sortedAveMaxAUC)-1)/length(aveMaxAUCsForAllBootstraps)
  #save the distribution of MAx UACs that determine the network P-Value
  vizAndSaveBootstrapAUCDist(sortedAveMaxAUC, testNetworkMaxAUC, testNetworkPV)
  
  ############################################################################
  # Finally, save pertinent statistics and characteristics of this run
  ############################################################################
  saveRunStats(ibasStartTime, proc.time(),paste(outputLocation,"SavedRunTimes/runTimes.log",sep=""),notesToThisRun)
  
}



