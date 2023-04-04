   #(Simulations d') Analyses du projet d'antibior?sistance d'E.coli 
#script cr?e le 24/03/2023 par Rolly Mpoyi

##data prep--------------------------------------------------------------------------------
#matrice vide
matrix <- matrix(, nrow = 10, ncol = 6)

#noms des colonnes 
colnames(matrix) <- c("Isolate ID","Ampicillin (AMP)", "Ciprofloxacin (CIP)", "Azithromycin (AZI)", "Tetracycline(TET)", "Year")

#transf en data frame 
matrix.df <- as.data.frame(matrix)

# ajouter les ann?es
x <- 2018:2022

#remplissage des rows  (R == 1, S == 0) 
matrix.df[1,] <- c(1, 1, 1, 1, 0, sample(x, 1, replace = TRUE))
matrix.df[2,] <- c(2, 0, NA, 0, 0, sample(x, 1, replace = TRUE))
matrix.df[3,] <- c(3, 0, 0, NA, 0, sample(x, 1, replace = TRUE))
matrix.df[4,] <- c(4, 0, 0, NA, 1, sample(x, 1, replace = TRUE))
matrix.df[5,] <- c(5, 1, 1, 1, 0, sample(x, 1, replace = TRUE))
matrix.df[6,] <- c(6, 1, 1, 1, 0, sample(x, 1, replace = TRUE))
matrix.df[7,] <- c(7, 0, NA, 0, 0, sample(x, 1, replace = TRUE))
matrix.df[8,] <- c(8, 0, 0, 1, 0, sample(x, 1, replace = TRUE))
matrix.df[9,] <- c(9, 1, 0, NA, 1, sample(x, 1, replace = TRUE))
matrix.df[10,] <- c(10, 1, 0, 0, 1, sample(x, 1, replace = TRUE))

matrix.df

# data values en num?rique
matrix_num <- sapply(matrix.df, as.numeric)
matrix_num <- as.data.frame(matrix_num)
# data values as logical TRUE = resist, FALSE = sensible
matrix_logic <- matrix.df
matrix_logic <- sapply(matrix_logic[,2:5], as.logical)
matrix_logic <- as.data.frame(matrix_logic)
matrix_logic <- cbind(matrix.df[1],matrix_logic, matrix.df[6])

matrix_logic


#fonction ? modifier pour changer les valeurs "R/S/I" en forme logique vrai/faux plus tard)
if (is.na(bp.index)==TRUE){ #if there bp.index for the i MIC column is NA, then there is no bp for that AM
  data[,i]<-NA #replace all MIC values with NA
  warning("missing breakpoint for ", colnames(data)[i], "\n") #warn that there is no bp for a column

  } else{ #if the bp.index is not NA, then there is a bp for that AM
  data[,i]<-discretize(data[,i], method="fixed", 
                       breaks=c(-Inf, bp$NSbp[bp.index], Inf), 
                       right=TRUE, labels=c("FALSE", "TRUE")) #discretize the MIC values at the bp. Intervals are right closed: (-Inf,bp]; (bp, Inf). Hence, MIC values >bp are given "True" (NS) and MIC values <=bp are given "False" (S)
  data[,i]=as.logical(data[,i]) #make logical
}

#### m?thode de Cazer et.al
require(arules)
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





##1?re m?thode: analyses descriptives-------------------------------------------------------------------------------
require(dplyr)
require(tidyr)

par(mfrow = c(2,2))
hist(matrix.df$`Ampicillin (AMP)`, main = "Fr?q de res d'AMP", xlab = "Resistance (1 or 0)")
hist(matrix.df$`Ciprofloxacin (CIP)`, main = "Fr?q de res de CIP", xlab = "Resistance (1 or 0)")
hist(matrix.df$`Azithromycin (AZI)`, main = "Fr?q de res d'AZI", xlab = "Resistance (1 or 0)")
hist(matrix.df$`Tetracycline(TET)`, main = "Fr?q de res de TET", xlab = "Resistance (1 or 0)")


#nombre d'isolats test?s par ann?es 

matrix_logic %>%
  count(Year) -> year_count

year_count

#group? par l'ann?e, ordon?e par ann?e descendante
matrix_logic %>%
  group_by(Year) %>%
  arrange(desc(Year)) -> year_group

#test antibiogramme
anti_bio <- function(data, index) {  #Propre fonction cr??e (sans groupement)
  #data = le data frame r?unissant toutes les donn?es 
  #l'indice des colonnes d'antibiotiques ? analyser
  
  #calcul pr?valence de chaque ATB arrondi ? 2 chiffres
  eSupport <- round(sapply(data[,index], mean, na.rm = TRUE),2)
  
  #nombre d'isolats r?sistants ? chaque ATB 
  Nombre_Isolat <- sapply(data[,index], sum, na.rm = TRUE)
  
  #total_isolat test?es 
  matrix_noNa <- !is.na(data[,index])
  matrix_noNa <- as.data.frame(matrix_noNa)
  Total_Isolat_test <- sapply(matrix_noNa, sum)
  
  #cr?ation de l'antibiogramme
  abg <- rbind(eSupport, Nombre_Isolat, Total_Isolat_test)
  abg
} 
antibiogramme <- anti_bio(matrix_logic, 2:5)

antibiogramme

antibiogram <- function(data, index, group) { #fonction reprise group?e

  require(dplyr)
  require(tidyr)
  
  abgm <- group_by(data[,c(group, index)] , data[,group]) %>%
  summarise_all(list(Prevalence=~mean(.,na.rm=TRUE),Number_Isolates_Tested=~sum(!is.na(.))))
  
  overall_prev <- summarise_all(data[,c(group,index)], list(Prevalence=~mean(.,na.rm=TRUE),Number_Isolates_Tested=~sum(!is.na(.))))
  overall_prev <- cbind("data[, group]"="Overall", overall_prev)
  abgm <- rbind(abgm, overall_prev)
  

  abgm <- abgm[,order(names(abgm))]
  abgm 
  
}
Year_abg <- antibiogram(matrix_logic,2:5, 6)

Year_abg

#test atb profile ie. les noms des antibiotiques test?es pour chaque isolat
atb_profile <- function(data, index){
  ATB_NOMS <- names(data)[index]  #prendre les noms des ATBS qui se trouvent dans chaque isolat
  atb_panels <- rep("", dim(data)[1]) #initialise un vecteur de caract?res blanc aussi long que le nombre d'isolats
  
  for(atb in ATB_NOMS) {	#pour chaque atb
    x <- with(data, get(atb)) #x est la col des valeurs pour les ATB 
    x <- ifelse(is.na(x), '', paste(atb, '')) #si NA, colle un espace blanc
    atb_panels <- paste(atb_panels, x, sep = '') #colle la col x. pour chaque isolat, atb_panels va accumuler les noms des ATB qui ont ?t? test?s pour cet isolat
  } 
  atb_panels #Retourne le profil test? de la sensibilit? atbiotique de chaque isolat
}

ATB_atb <- atb_profile(matrix.df,2:5)

ATB_atb

#test mdr_profile
#Pour chaque isolat, identifie tout les antibiotiques resist?s par cet isolat, donne le nombre d'atbs qu'il r?siste
#data: data des ATBS test?s 
#index: indices des colonnes d'ATBS

mdr_profile <- function (data, index){
  #packages requis
  require(stringr)
  require(mgsub)
  
  z <- data[,index] #colonne des ATBs
  
  mdr_profile <- data.frame(Profile=rep("", dim(z)[1]), NombreRes=numeric(dim(z)[1])) #cr?ation dune df pour stocker les profils
  profile <- rep("", dim(z)[1]) #Initialise un caract?re blanc aussi long que le nombre d'isolat dans le df
  
  
  mdr_profile[,2] <- rowSums(z, na.rm=TRUE) #nombre de resistances aux ATBs pour chaque isolat
  
  for(i in 1:ncol(z)) { #pour chaque ATB
    z[which(z[,i]==TRUE), i] <- as.character(colnames(z)[i]) #si resistant, colle le nom de l'ATB dans la cellule
    z[which(z[,i]==FALSE), i] <- as.character("") #si sensible, colle un espace blanc dans la cellule
    z[is.na(z[,i]), i] <- as.character("") #ceux qui n'ont pas ?t? test?es, colle un espace blanc
    profile <- paste(profile, z[,i], sep = " ") #colle tous les cols d'ATB ensemble pour avoir une chaine de caracs de tous les ATB pour lesquels chaque isolat (chaque ligne) est resistant
  }
  profile <- trimws(profile, which="both") #suppr les espaces blancs
  mdr_profile[,1] <- profile #stocker les profils de resistance
  
  mdr_profile #sortie
}
profil_res <- mdr_profile(matrix_logic,2:5)

profil_res
unique(profil_res) #retourne les diff?rents profils de r?sistances 

#calculs ph?notype uniques, fr?q de ph?notypes 
#pour chaque nombre de r?sistances aux ATBs, calcule le nombre d'isolats et de ph?notypes uniques
multi_res <- as.data.frame(table(profil_res$NombreRes)) #nombre de pattern de resistance et d'isolats 
names(multi_res) <- c("NombreResATB", "NombreIsolats")

multi_res

#nombre de pattern de res uniques parmi le nombre de resistances retenu
for (i in 0:max(profil_res$NombreRes)){ #pour chauqe pattern de resistance 
  multi_res[i+1,3] <- nrow(unique(profil_res[which(profil_res$NombreRes==i),])) #nombre de res uniques pour  chaque profil de multi-resistances
}
names(multi_res)[3]="NombreUniqPatterns"

multi_res

##2?me m?thode: test statistique d'ind?pendance--------------------------------------------------------------------------------
#disribution pr?dite des multiresistances et pattern de res selon l'hypoth?se nulle de l'ind?pendance des traits de r?sistances 
#simule 100 jeu de donn?es 

set.seed(500) #


rand_matrix_logic <- rep(list(matrix_logic[,2:5]),100) #copie les donn?es des ATBs
n_obs <- nrow(matrix.df) #nombre d'observations
n_nonNA <- n_obs - colSums(is.na(matrix_logic[,2:5])) #vecteur du nombre de non-NA pour chaque atb
items <- sapply(matrix_logic[,2:5], sum, na.rm=TRUE)/n_nonNA #calcul de la prevalence de chaque resistance

rand_matrix_logic

#pour chacune des 100 jeu de donnees
for (i in 1:100){ #pour chacune des 100 jeu de donnees 
  for (k in 1:length(items)){ #pour chaque antibiotique
    bin_dat <- rbinom(n=n_nonNA[k], size=1, prob=items[k]) #valeurs binomiales ind?pendantes g?ner?es, taille = nombre de non-NA, probabilit? =% de res
    rand_dat <- rand_matrix_logic[[i]][,k] #donne la col pour les ATB de la base avec les m?mes NA
    rand_dat[!is.na(rand_dat)] <- bin_dat #replace les valeurs non-NA avec des valeurs binomiales ind?pendantes
    rand_matrix_logic[[i]][,k] <- as.logical(rand_dat) #sauvegarde dans la base de donnees
    
  }
} 

#exemples des jeux de donn?es simul?s 
rand_matrix_logic[[2]]
rand_matrix_logic[[3]]
rand_matrix_logic[[4]]

rand_MDRProfiles <- data.frame()

for (i in 1:100){ #pour chacunes des 100 data frame
  x <- mdr_profile(rand_matrix_logic[[i]], seq(1:ncol(rand_matrix_logic[[i]])))
  x <- mutate(x, db=rep(i, nrow(x))) #ajoute l'identification de la data frame
  rand_MDRProfiles <- rbind(rand_MDRProfiles, x) #joint tout les profils de multi-resistance des 100 df ensembler
}

rand_MDRProfiles #profils de multi-resistance pour chaque isolat dans chaque jeu de donnees simul?s

#creer une table des multi-r?sistances
#calculer le nombre de resistances pour chaque isolat
rand_MDRProfiles$NombreRes <- as.factor(rand_MDRProfiles$NombreRes)
rand_MDRProfiles$db <- as.factor(rand_MDRProfiles$db)

#à faire pour i de 1 à 100
unique(mdr_profile(rand_matrix_logic[[i]],1:4))
