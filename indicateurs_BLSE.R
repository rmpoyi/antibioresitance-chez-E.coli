library(magrittr)
library(plyr)
library(dplyr)
library(arules)
library(data.table)

#Préparation

data2018<-read.csv(file="EC-2018.csv",sep=";")


BLSE2018 <- data2018
unique(BLSE2018$phenotype) #on choisit les phénotypes dit BLSE, BLSE+IMPER (4) et BLSE+CARBA (2)

BLSE2018 <- filter(BLSE2018, type.prelevement == "URINES" & phenotype == "BLSE" | phenotype == "BLSE+CARBA" | phenotype == "BLSE+IMPER")

unique(BLSE2018$phenotype)
unique(BLSE2018$type.prelevement)

BLSE2018<- BLSE2018[-which(BLSE2018$type.prelevement == "PUS-SUPERFICIEL"),]


BLSE2018 <- subset(BLSE2018, select = -c(TCC, CN, AZT, CS, TGC, C3G, FQ, AMX, AMP)) #supprimer les colonnnes pas pertinentes 

BLSE2018_res <- BLSE2018[,11:40]

BLSE2018_res[BLSE2018_res == ''] <- NA
BLSE2018_res[BLSE2018_res == 'S'] <- FALSE
BLSE2018_res[BLSE2018_res == 'I'] <- FALSE
BLSE2018_res[BLSE2018_res == 'R'] <- TRUE


BLSE2018_res_log <- BLSE2018_res %>% mutate_all(as.logical) 
BLSE2018_res_log <- data.table(BLSE2018_res_log)

BLSE18trans <- as(BLSE2018_res_log, "transactions")

#atb18trans_ql = atb2018_res_log %>%
#distinct()

otrans18 <- BLSE2018_res_log[,itemLabels((BLSE18trans))]
otrans18[is.na(otrans18)] <- as.logical("TRUE")
otrans18 <- as(otrans18, "transactions")

set18 <- apriori(BLSE18trans, parameter = list(support=1/length(BLSE18trans), maxlen=22, minlen=2, maxtime= 25, target="frequent itemsets"))
set18@quality$CrossSupRatio <- interestMeasure(set18, "crossSupportRatio", BLSE18trans, reuse = TRUE)
set18@quality$lift <- interestMeasure(set18, "lift", BLSE18trans, reuse=TRUE)
inspect(head(set18))

pSup <- interestMeasure(set18, "support", BLSE18trans, reuse=TRUE)
oSup <- interestMeasure(set18, "support", otrans18, reuse=FALSE)
eSup <- pSup/(1- (oSup - pSup))

set18@quality$eSupport <- eSup

inspect(head(set18))

setlist <- LIST(items(set18), decode = FALSE)

ecsr <- expectedCSR(BLSE18trans, otrans18, setlist)

set18@quality$eCSR <- ecsr

inspect(head(set18))

elift <- expectedLift(BLSE18trans, otrans18, setlist, eSup)

set18@quality$eLift <- elift
inspect(head(set18))

assign("SetName", set18)

rm(set18, ecsr, elift, eSup)

save(BLSE2018_res_log, file = "dataBLSE18.Rdata")
save(BLSE18trans, file = "BLSE18trans.Rdata")
save(otrans18, file = "BLSE18otrans.Rdata")
save(setlist, file = "BLSE18setlist.Rdata")
save(SetName, file = "BLSE18ruleset.Rdata")



testBLSE <- DATAFRAME(SetName)

save(testBLSE, file = "BLSE18.Rdata")

write.csv(testBLSE, "test_data.csv", row.names = T)
fwrite(testBLSE, "test_data.csv")
testBLSE <- data.table(testBLSE)

BLSE2018_res_log <- data.frame(BLSE2018_res_log)

"-------------------------------------------"
  
rand_dbNames="rand_BLSE18"

rand_setNames="rand_BLSE18sets"

dbNames <- "BLSE2018_res_log"
transNames <- "BLSE18trans"
setNames <- "SetName"

maxlen = 22
minlen = 2

load("dataBLSE18.Rdata")
load("BLSE18trans.Rdata")
load("BLSE18ruleset.Rdata")
load("BLSE18otrans.Rdata")


Dist.out <- PValue.sets(dbNames, varNames, setNames, rand_dbNames, rand_setNames, 2, length(BLSE2018_res_log))

for (j in seq_along(rand_dbNames)){
  set.seed(500) #set seed so that function is repeatable
  db <- rep(list(get(dbNames[j])),100) #copy actual data to maintain NA positions
  odb <- list() #storage for optimistic databases without missing data
  n_obs <- nrow(get(transNames[j])) #number of observations
  n_nonNA <- n_obs - colSums(is.na(get(dbNames[j]))) #vector of number of non-NA positions for each antimicrobial
  items <- sapply(get(dbNames[j]), sum, na.rm=TRUE)/n_nonNA #frequency of each item, excluding NA from denominator because binomial generator will also ingore NA positions
  
  
  
  for (i in 1:100){ #for each of the 100 random databases
    for (k in 1:length(items)){ #for each item
      bin_dat <- rbinom(n=n_nonNA[k], size=1, prob=items[k]) #generate random binary data (as 0,1)
      rand_dat <- db[[i]][,k] #get column for AM k
      rand_dat[!is.na(rand_dat)] <- bin_dat #replace non-NA values with random binary data
      db[[i]][,k] <- as.logical(rand_dat) #save
    }
    

    
        
    #make optimistic random db by replacing NA with TRUE
    #need to remove AM that were not tested at all since those will be excluded when db is made into transactions. Also need to remove AM that had only FALSE and NA since those will also be exluded from db when it is made into transactions
    odb[[i]] <- db[[i]][colSums(db[[i]], na.rm=TRUE)>0] #keep only AM with at least some TRUE values
    odb[[i]][is.na(odb[[i]])]  <-  TRUE #make optimisitic database where missing data are replaced by resistance
  }
  assign(rand_dbNames[j],list("db"=db, "odb"=odb)) #save to rand_db
}

rm(bin_dat, i, j, k, n_nonNA, n_obs, rand_dat, db, odb)


for (j in seq_along(rand_setNames)){ #for each database
  rand_sets <- vector("list", length=100) #storage for sets
  db <- get(rand_dbNames[j])[["db"]] #get the appropriate random_db and optimistic random_db
  odb <- get(rand_dbNames[j])[["odb"]] 
  
  for (i in 1:100){ #mine each of the 100 rand_db
    #first make transaction sets
    dat <- db[[i]]
    odat <- odb[[i]]
    trans <- as(dat, "transactions") #make transaction database
    otrans <- as(odat, "transactions") #optimistic transaction database
    
    #warn if dimensions are not equal and item frequencies not appropriate--could indicate issue with dropping AM with no resistance from odb
    if(any(dim(trans)!=dim(otrans))){
      warning("Different dimensions in transaction and optmistic_transactions. STOP")
    }
    if(any(itemFrequency(trans)>itemFrequency(otrans))){
      warning("issue with item freq in trans databases")
    }
    
    #next mine sets from trans
    minsup <- 1/length(trans) #min support is default 1 isolate
    sets <- apriori(trans, 
                    parameter=list(support=minsup, 
                                   maxlen=maxlen, 
                                   minlen=minlen, 
                                   target="frequent itemsets"),
                    control=list(verbose=FALSE)) #mine sets
    
    #get quality measures
    itemset_list <- LIST(items(sets), decode = FALSE) #need for eSup, eLift, eCSR
    #qualmeasures <- interestMeasure(sets, c("crossSupportRatio", "lift"), trans, reuse=TRUE)
    esupmeasures <- expectedSupport(sets, trans, otrans)
    ecsrmeasures <- expectedCSR(trans, otrans, itemset_list)
    eliftmeasures <- expectedLift(trans, otrans, itemset_list, esupmeasures)
    
    #store quality measures
    #sets@quality$csr <- qualmeasures$crossSupportRatio
    #sets@quality$lift <- qualmeasures$lift
    sets@quality$eSup <- esupmeasures
    sets@quality$eLift <- eliftmeasures
    sets@quality$eCSR <- ecsrmeasures
    #sets@quality$cLift <- e_qualmeasures$cLift
    rand_sets[[i]] <- sets #store sets
    rm(esupmeasures, ecsrmeasures, eliftmeasures, sets, itemset_list, minsup, trans, otrans)
  }
  assign(rand_setNames[j],rand_sets) #store sets from all 100 random_db
}
