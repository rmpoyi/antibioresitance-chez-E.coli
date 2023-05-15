library(arules)
library(dplyr)
library(plyr)

#1er essai-----------------------------------------

log18.trans <- as(log18, "transactions")

log18.sets <- apriori(log18.trans, parameter = list(support = 1/length(log18.trans), maxlen = 22, minlen = 2, target = "frequent itemsets"))
log18.sets@quality$CrossSupRatio <- interestMeasure(log18.sets, "CrossSupportRatio", log18.trans, reuse = TRUE)

calculs <- inspect(log18.sets)

summary(calculs)

simul <- data.frame()

for (bd in 1:100) {
  x <- randlog18[[bd]]
  x <- mutate(x, db = rep(bd, nrow(x)))
  simul <- rbind(simul, x)
  
}
simul.trans <- as(simul, "transactions")

simul.sets <- apriori(simul.trans, parameter = list(support = 1/length(simul.trans), maxlen = 22, minlen = 2, target = "frequent itemsets"))
simul.sets@quality$CrossSuperRatio <- interestMeasure(simul.sets, "CrossSupportRatio", simul.trans, reuse = TRUE)


data2018<-read.csv(file="EC-2018.csv",sep=";")

data2018<-data2018[data2018$type.prelevement=="URINES",]

atb2018_res <- data2018[,11:49]

atb2018_res <- subset(atb2018_res, select = -c(TCC, CN, AZT, CS, TGC, C3G, FQ))


atb2018_res[atb2018_res == ''] <- NA
atb2018_res[atb2018_res == 'S'] <- FALSE
atb2018_res[atb2018_res == 'I'] <- FALSE
atb2018_res[atb2018_res == 'R'] <- TRUE


atb2018_res_log <- atb2018_res %>% mutate_all(as.logical)

data2019<-read.csv(file="EC-2019.csv",sep=";")

data2019<-data2019[data2019$type.prelevement=="URINES",]

atb2019_res <- data2019[,11:49]

atb2019_res <- subset(atb2019_res, select = -c(TCC, CN, AZT, CS, TGC, C3G, FQ))


atb2019_res[atb2019_res == ''] <- NA
atb2019_res[atb2019_res == 'S'] <- FALSE
atb2019_res[atb2019_res == 'I'] <- FALSE
atb2019_res[atb2019_res == 'R'] <- TRUE

atb2019_res_log <- atb2019_res %>% mutate_all(as.logical)


# EQM fonction------- 

expectedQM <- function(setName, transName, oTransName, dbName, setList){
  if(any(dim(transName)!=dim(oTransName))){ #transaction and optimistic transaction databases must have same columns
    warning("Different dimensions of transactions and optimistic transactions. STOP")
  }
  
  if(length(setName)!=length(setList)){ #setList must correspond to the sets
    warning("Different number of sets in setName and setList. STOP")
  }
  
  #calculate pSup, oSup, eSup for each itemset
  pSup <- interestMeasure(setName, "support", transName, reuse=TRUE)
  oSup <- interestMeasure(setName, "support", oTransName, reuse=FALSE)
  eSup <- pSup/(1- (oSup - pSup))
  
  #check that expected conditions are met
  if (all((pSup - eSup) < 1e-10)==FALSE){ #pSup should be <= eSup although a rounding error due to division in calculating eSup may result in small positive value
    warning("Expected Support is less than Pessimistic Support.")
  }
  
  if (all((eSup - oSup) < 1e-10)==FALSE){ #eSup should be <= oSup although a rounding error due to division in calculating eSup may result in small positive value
    warning("Optimistic Support is less than Expected Support.")
  }
  
  #calculate pSup, oSup, eSup for each individual AM
  pItemSupport <- itemFrequency(transName)
  oItemSupport <- itemFrequency(oTransName)
  eItemSupport <- pItemSupport/(1-(oItemSupport-pItemSupport))
  
  #eLift
  eLift <- eSup / sapply(setList, function(i) prod(eItemSupport[i]))
  
  #eCSR
  eCSR <- sapply(setList, function(i) min(eItemSupport[i])) / sapply(setList, function(i) max(eItemSupport[i]))
  
  #cLift
  AMsInSet <- lapply(setList, function(i) itemLabels(transName)[i]) #for each itemset, a character vector of the AMs in the set
  OnlyAMsInSet <- lapply(AMsInSet, function(i) subset(dbName, select=i) %>% na.omit()) #in the binary database, select only AM columns in the set and isolates that had been tested against all of them (no NA). creates a list of binary databases, one for each itemset
  setCount <- interestMeasure(setName, "support", transName, reuse=TRUE) * length(transName) #for each itemset, the number of isolates with the set (all AM in the set must have been tested in these isolates for them to be counted in the support)
  numerator <- setCount/ldply(OnlyAMsInSet, nrow) #P(set | all AM in set are tested)
  denominator <- ldply(OnlyAMsInSet, function(i) prod(sapply(i, mean))) #find the support of each individual AM in the set (mean applied to logical vector) and multiply them together
  cLift <- numerator/denominator
  
  #save together as dataframe
  eQM <- as.data.frame(cbind(pSup, oSup, eSup, eLift, eCSR, cLift))
  colnames(eQM) <- c("pSup", "oSup", "eSup", "eLift", "eCSR", "cLift")
  eQM
}


#transactions databases------


dbNames <- c("atb2018_res_log") #, "atb2019_res_log") # "atb2020_res_log", "atb2021_res_log"
varNames <- c("atb2018") # , "atb2019") # "atb2020", "atb2021"

for (i in 1:length(varNames)){ #for each year
  x <- as(get(dbNames[i]), "transactions") #make binary data into transactions
  label <- paste('atb',as.character(paste(2017+i)), sep="")
  assign(label, x) #save
  rm(x)
  rm(label)
}


#atb2018_res_log <- cbind(atb2018_res_log, Year = rep(2018, nrow(atb2018_res_log)))

#atb2019_res_log <- cbind(atb2019_res_log, Year = rep(2019, nrow(atb2019_res_log)))

#atb2020_res_log <- cbind(atb2020_res_log, Year = rep(2020, nrow(atb2020_res_log)))

#atb2021_res_log <- cbind(atb2021_res_log, Year = rep(2021, nrow(atb2021_res_log)))

#atb.db <- rbind(atb2018_res_log, atb2019_res_log, atb2020_res_log, atb2021_res_log)


oTransNames <- c("atb2018_o") #, "atb2019_o") # "atb2020_o", "atb2021_o")
for (i in 1:length(oTransNames)){ #for each year
  x <- atb2018_res_log[itemLabels(get(varNames[i]))] #select only AM columns that are in the actual transaction database (e.g have at least one isolate resistant), otherwise columns with all NA will be included and then turned into all resistant
  x[is.na(x)] <- as.logical("TRUE") #replace NA with resistant 
  x <- as(x, "transactions") #make binary data into transactions
  label <- paste('atb',as.character(paste(2017+i)), "_o", sep="")
  assign(label, x) #save
  rm(x)
  rm(label)
}


setNames=c("atb2018_sets") #, "atb2019_sets", "atb2020_sets", "atb2021_sets")

atb2018_sets <- apriori(atb2018, parameter=list(support=1/length(atb2018)))

for (i in seq_along(varNames)){ #for each year
  #get relevant transaction, optimistic transaction, and binary databases
  data <- get(varNames[i]) 
  odata <- get(oTransNames[i])
  db <- get(dbNames[i])
  minsup <- 1/length(data) #minimum support = 1 isolate resistant
  #mine sets
  sets <- apriori(data, parameter=list(support=minsup,  maxlen=length(data2018), minlen=2, target="frequent itemsets"))
  
  #get quality measures
  itemset_list <- LIST(items(sets), decode = FALSE) #required for expectedQM
  CrossSupRatio <- interestMeasure(sets, "crossSupportRatio", data, reuse=TRUE)
  lift <- interestMeasure(sets, "lift", data, reuse=TRUE)
  sets@quality$CrossSupRatio <- CrossSupRatio #save to set database
  sets@quality$lift <- lift
  eQM <- expectedQM(sets, data, odata, db, itemset_list) #calculate expected quality measures (QM)
  sets@quality$eSup <- eQM$eSup
  sets@quality$eCSR <- eQM$eCSR
  sets@quality$eLift <- eQM$eLift
  sets@quality$cLift <- eQM$cLift
  
  #save
  assign(setNames[i], sets)
  
  rm(data, odata, db, minsup, sets, itemset_list, CrossSupRatio, lift, eQM)
}

 #eSupport function calculates expected Support for each itemset in a set database (currently implemented for one set database at a time so that eSupport can be implemented within apriori mining step)
#requires:
#transName: a character label for the transaction database used to mine set database
#setName: a character label for the set database
#oTransName: a character label for the optimistic transaction database with missing values replaced by TRUE. must have same columns as transaction database
#setList: itemset_list <- LIST(items(setName), decode = FALSE). Faster to create this variable outside the function since it is used in eSupport, eCSR, and eLift
expectedSupport <- function(setName, transName, oTransName){
  if(any(dim(transName)!=dim(oTransName))){ #transaction and optimistic transaction databases must have same columns
    warning("Different dimensions of transactions and optimistic transactions. STOP")
  }
  
  pSup <- interestMeasure(setName, "support", transName, reuse=TRUE) #pSup is standard support calculated by arules
  oSup <- interestMeasure(setName, "support", oTransName, reuse=FALSE) #oSup is the standard support calculated when missing values are set to TRUE/resistant
  eSup <- pSup/(1- (oSup - pSup)) #expected support
  
  if (all((pSup - eSup) < 1e-10)==FALSE){ #pSup should be <= eSup although a rounding error due to division in calculating eSup may result in small positive value
    warning("Expected Support is less than Pessimistic Support.") #if pSup>eSup, there may be something wrong with optimistic transaction database
  }
  
  if (all((eSup - oSup) < 1e-10)==FALSE){ #eSup should be <= oSup although a rounding error due to division in calculating eSup may result in small positive value
    warning("Optimistic Support is less than Expected Support.") #if eSup>oSup, there may be something wrong with optimistic transaction database
  }
  
  eSup
}

esupp <- expectedSupport(atb2018_sets, atb2018, atb2018_o)
