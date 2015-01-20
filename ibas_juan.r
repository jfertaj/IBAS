library('igraph')
library('plyr')

WNTNodeNumbers<-c('11961','1236','13075','13308','13971','16377','22173','22730','329','3687','7590')
TGFNodeNumbers <- c('11388', '12424','16681','17783','22488','22780','3803','4554', '6216', '8964', '9049')
RR3KNodeNumbers <- c('10090', '1064','11032', '12222', '1635', '4695', '4792', '5746', '7622', '8226','857')

#init the testNetwork to WNT for now
testNetwork <- WNTNodeNumbers

networkTime <- proc.time()
read.table("/Users/jfertaj/INVESTIGACION/WELLCOME/Data_Integration/IBAS_Juan/InWeb3/InWeb3_score.eda", skip=1) -> inweb3
inweb3 <- inweb3[,c(1,3,5)]
inweb.g <- graph.data.frame(inweb3, directed=F)

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
runTimeString <- returnMinSec((proc.time() - ibasstartTime)[3])
print(paste("num of BOOTSTRAPS = ", numBootstraps, " took: ", runTimeString, sep=""))














    getGenesWithSameDegree <- function(gene){
    	grado <- degree(inweb.g, gene)
		genes <- names(which(degree(inweb.g)==grado))
		genes <- genes[-which(genes==)]
      
     	numFound <- length(genes)
     	if (numCVGenes > numFound){




     	}else{sample(genes, numCVgenes)}
      

length(names(which(inweb.g, '11947')+1)))


      #first get the cross validation set of genes for this bootstrap gene
      CVSet <- getCVSet(gene)

 prueba <- data.frame(grado=degree(inweb.g), genes=names(degree(inweb.g)))

 prueba.split <- split(prueba, prueba[,1])
 prueba.split <- lapply(prueba.split, '[[', 2)
 prueba.split <- lapply(prueba.split, as.character)

l3 <- list()
ll1 <- length(prueba.split)
length(l3) <- ll1
for (i in 1:ll1){
   vec1 <- prueba.split[[i]]
   jl <- 1;jr<-1;
    while (length(vec1) < 2000){
       if(i==1 || i-jl==0) {
          vec1 <- c(vec1, prueba.split[[i+jr]])
          jr <- jr+1
        } else if (i==ll1 || jr+i==ll1 ){
           vec1 <- c(vec1, prueba.split[[i-jl]])
           jl <- jl+1
 }else {
            vec1 <- c(vec1, prueba.split[[i-jl]], prueba.split[[i+jr]])
        jl <- jl+1
        jr <- jr+1
          } 
  } 
    #l3[[i]] <- sample2(vec1, 10) 
    l3[[i]] <- vec1  
}



