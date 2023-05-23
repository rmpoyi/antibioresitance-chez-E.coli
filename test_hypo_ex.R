library(arules)
library(purrr)
library(dplyr)
library(plyr)

example2 <- data.frame("Ampicillin" = c("TRUE", "FALSE", "FALSE", "FALSE", "TRUE", "TRUE", "FALSE", "FALSE", "TRUE", "TRUE"),
                       "Ciprofloxacin" = c("TRUE", "NA", "FALSE", "FALSE", "TRUE", "TRUE", "NA", "FALSE", "FALSE", "FALSE"),
                       "Azithromycin" = c("TRUE", "FALSE", "NA", "NA", "TRUE", "TRUE", "FALSE", "TRUE", "NA", "FALSE"),
                       "Tetracycline" = c("FALSE", "FALSE", "FALSE", "TRUE", "FALSE", "FALSE", "FALSE", "FALSE", "TRUE", "TRUE"))

example2 <- as.data.frame(sapply(example2, as.logical))
example2.trans <- as(example2, "transactions")
summary(example2.trans)
example2 
summary(example2) 

example2.sets <- apriori(example2.trans, parameter=list(support=1/length(example2.trans),  maxlen=4, minlen=2, target="frequent itemsets"))
example2.sets@quality$CrossSupRatio <- interestMeasure(example2.sets, "crossSupportRatio", example2.trans, reuse=TRUE)
example2.sets@quality$lift <- interestMeasure(example2.sets, "lift", example2.trans, reuse=TRUE)
inspect(example2.sets) 

example2.otrans <- example2
example2.otrans[is.na(example2.otrans)] <- as.logical("TRUE")
example2.otrans <- as(example2.otrans, "transactions")

varNames <- "example2.trans"
oTransNames <- "example2.otrans"
dbNames <- "example2"
setNames <- "example2.sets"


for (i in seq_along(varNames)){ #for each year
  #get relevant transaction, optimistic transaction, and binary databases
  data <- get(varNames[i]) 
  odata <- get(oTransNames[i])
  db <- get(dbNames[i])
  minsup <- 1/length(data) #minimum support = 1 isolate resistant
  #mine sets
  sets <- apriori(data, parameter=list(support=minsup,  maxlen=length(example2), minlen=2, target="frequent itemsets"))
  
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

df_example2 <- DATAFRAME(example2.sets)

rand_dbNames = "rand_example2"
rand_setNames = "rand_example2.sets"

Dist.out <- PValue.sets(dbNames, varNames, setNames, rand_dbNames, rand_setNames, 2, length(example2))

for (i in seq_along(setNames)){ #for each db of sets
  sets <- get(setNames[i]) #get sets
  
  #find percentile of each actual QM
  perc.ecsr <- as.numeric(row.names(Dist.out[["set.ecsr_pct_avg"]])
                          [adply(sets@quality$eCSR, .margins=1, #margins=1 splits by rows--so for each set eCSR value x
                                 function(x) detect_index(Dist.out[["set.ecsr_pct_avg"]][,i], #search within the null distribution for the QM and the relevant db i
                                                          function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% #find where this function is true--the avg null QM value is <= set eCSR value
                              mapvalues(0,1, warn_missing=FALSE)]) 
  
  perc.csr <- as.numeric(row.names(Dist.out[["set.csr_pct_avg"]])
                         [adply(sets@quality$CrossSupRatio, .margins=1, 
                                function(x) detect_index(Dist.out[["set.csr_pct_avg"]][,i], 
                                                         function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% 
                             mapvalues(0,1, warn_missing=FALSE)]) 
  
  perc.eLift <- as.numeric(row.names(Dist.out[["set.elift_pct_avg"]])
                           [adply(sets@quality$eLift, .margins=1, 
                                  function(x) detect_index(Dist.out[["set.elift_pct_avg"]][,i], 
                                                           function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% 
                               mapvalues(0,1, warn_missing=FALSE)]) 
  
  perc.lift <- as.numeric(row.names(Dist.out[["set.lift_pct_avg"]])
                          [adply(sets@quality$lift, .margins=1, 
                                 function(x) detect_index(Dist.out[["set.lift_pct_avg"]][,i], 
                                                          function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% 
                              mapvalues(0,1, warn_missing=FALSE)]) 
  
  perc.cLift <- as.numeric(row.names(Dist.out[["set.clift_pct_avg"]])
                           [adply(sets@quality$cLift, .margins=1, 
                                  function(x) detect_index(Dist.out[["set.clift_pct_avg"]][,i], 
                                                           function(y) y<=x, .dir="backward"), .id=NULL)[,1] %>% 
                               mapvalues(0,1, warn_missing=FALSE)]) 
  
  #save as P-values
  sets@quality$Pval.eCSR <- 1-perc.ecsr
  sets@quality$Pval.CSR <- 1-perc.csr
  sets@quality$Pval.eLift <- 1-perc.eLift
  sets@quality$Pval.lift <- 1-perc.lift
  sets@quality$Pval.cLift <- 1-perc.cLift
  
  assign(setNames[i], sets)
  rm(sets, perc.ecsr, perc.csr, perc.eLift, perc.lift, perc.cLift)
}

df_final_example2 <- DATAFRAME(example2.sets)
