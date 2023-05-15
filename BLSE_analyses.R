library(magrittr)
library(plyr)
library(dplyr)
library(ggplot2)

#test 2018-------------------
data2018<-read.csv(file="EC-2018.csv",sep=";")


BLSE2018 <- data2018
unique(BLSE2018$phenotype) #on choisit les phénotypes dit BLSE, BLSE+IMPER (4) et BLSE+CARBA (2)

BLSE2018 <- filter(BLSE2018, type.prelevement == "URINES" & phenotype == "BLSE" | phenotype == "BLSE+CARBA" | phenotype == "BLSE+IMPER")

unique(BLSE2018$phenotype)
unique(BLSE2018$type.prelevement)

BLSE2018<- BLSE2018[-which(BLSE2018$type.prelevement == "PUS-SUPERFICIEL"),]


BLSE2018 <- subset(BLSE2018, select = -c(TCC, CN, AZT, CS, TGC, C3G, FQ, AMX, AMP)) #supprimer les colonnnes pas pertinentes 

BLSE2018 <- BLSE2018[,11:40]

#test antibiotique 
BLSE2018test <- BLSE2018

BLSE2018test[BLSE2018test == ''] <- 0
BLSE2018test[BLSE2018test == 'S'] <- 1
BLSE2018test[BLSE2018test == 'I'] <- 1
BLSE2018test[BLSE2018test == 'R'] <- 1

SeqBLSE2018test<-unique(BLSE2018test)
dim(SeqATB2018)[1]
#Ajout d'une colonne indiquant la fréquence de chaque jeu d'atb testé
FreqBLSE2018test<-BLSE2018test %>% group_by_all() %>% count() %>% ungroup()
#Dix jeux d'antibiotiques testés le plus fréquemment
TenTopBLSE2018test<-FreqBLSE2018test %>%
  arrange(desc(n)) %>%
  slice_head(n=10)

#Proportion des isolats testés concernés par ces dix jeux d'antibiotiques
sum(TenTopBLSE2018test[,33])/dim(BLSE2018test)[1]
#Dataframe décrivant ces dix jeux, dans l'ordre de fréquence décroissante
TenTopBLSE2018test <- TenTopBLSE2018test[order(TenTopBLSE2018test$n,decreasing=TRUE),]
#Génération d'une table au format txt à partir de cette dataframe
write.table(TenTopBLSE2018test,sep = ",", quote = FALSE, row.names = F)
#Possible de lire cette table dans Word ou Excel
#Calcul de la fréquence de test de chaque atb
BLSE2018testN <- BLSE2018test %>% mutate_all(as.numeric)
PropBLSE2018<-colSums(BLSE2018testN)/dim(BLSE2018test)[1]

write.table(PropBLSE2018,sep = ",", quote = FALSE, row.names = T)


#test 2019---------------

data2019<-read.csv(file="EC-2019.csv",sep=";")

BLSE2019 <- data2019
unique(BLSE2019$phenotype) #on choisit les phénotypes dit BLSE, BLSE+IMPER (4) et BLSE+CARBA (2)

BLSE2019 <- BLSE2019[BLSE2019$phenotype == "BLSE",]

unique(BLSE2019$phenotype)

BLSE2019 <- subset(BLSE2019, select = -c(AMC,AMC_urine, AMC03, TCC, TZP, C3G, FQ)) #supprimer les colonnnes pas pertinentes 

BLSE2019 <- BLSE2019[,11:42]

#test antibiotique 
BLSE2019test <- BLSE2019

BLSE2019test[BLSE2019test == ''] <- 0
BLSE2019test[BLSE2019test == 'S'] <- 1
BLSE2019test[BLSE2019test == 'I'] <- 1
BLSE2019test[BLSE2019test == 'R'] <- 1

SeqBLSE2019test<-unique(BLSE2019test)
dim(SeqBLSE2019test)[1]
#Ajout d'une colonne indiquant la fréquence de chaque jeu d'atb testé
FreqBLSE2019test<-BLSE2019test %>% group_by_all() %>% count() %>% ungroup()
#Dix jeux d'antibiotiques testés le plus fréquemment
TenTopBLSE2019test<-FreqBLSE2019test %>%
  arrange(desc(n)) %>%
  slice_head(n=10)

#Proportion des isolats testés concernés par ces dix jeux d'antibiotiques
sum(TenTopBLSE2019test[,33])/dim(BLSE2019test)[1]
#Dataframe décrivant ces dix jeux, dans l'ordre de fréquence décroissante
TenTopBLSE2019test <- TenTopBLSE2019test[order(TenTopBLSE2019test$n,decreasing=TRUE),]
#Génération d'une table au format txt à partir de cette dataframe
write.table(TenTopBLSE2019test,sep = ",", quote = FALSE, row.names = F)
#Possible de lire cette table dans Word ou Excel
#Calcul de la fréquence de test de chaque atb
BLSE2019testN <- BLSE2019test %>% mutate_all(as.numeric)
PropBLSE2019<-colSums(BLSE2019testN)/dim(BLSE2019test)[1]

write.table(PropBLSE2019,sep = ",", quote = FALSE, row.names = T)


#test 2020-------------------------------

data2020<-read.csv(file="EC-2020.csv",sep=";")

BLSE2020 <- data2020
unique(BLSE2020$phenotype) #on choisit les phénotypes dit BLSE, BLSE+IMPER (4) et BLSE+CARBA (2)

BLSE2020 <- BLSE2020[BLSE2020$phenotype == "BLSE",]

unique(BLSE2020$phenotype)

BLSE2020 <- subset(BLSE2020, select = -c(AMC,AMC_urine, AMC03, TCC, TZP, C3G, FQ)) #supprimer les colonnnes pas pertinentes 

BLSE2020 <- BLSE2020[,11:42]

#test antibiotique 
BLSE2020test <- BLSE2020

BLSE2020test[BLSE2020test == ''] <- 0
BLSE2020test[BLSE2020test == 'S'] <- 1
BLSE2020test[BLSE2020test == 'I'] <- 1
BLSE2020test[BLSE2020test == 'R'] <- 1

SeqBLSE2020test<-unique(BLSE2020test)
dim(SeqBLSE2020test)[1]
#Ajout d'une colonne indiquant la fréquence de chaque jeu d'atb testé
FreqBLSE2020test<-BLSE2020test %>% group_by_all() %>% count() %>% ungroup()
#Dix jeux d'antibiotiques testés le plus fréquemment
TenTopBLSE2020test<-FreqBLSE2020test %>%
  arrange(desc(n)) %>%
  slice_head(n=10)

#Proportion des isolats testés concernés par ces dix jeux d'antibiotiques
sum(TenTopBLSE2020test[,33])/dim(BLSE2020test)[1]
#Dataframe décrivant ces dix jeux, dans l'ordre de fréquence décroissante
TenTopBLSE2020test <- TenTopBLSE2020test[order(TenTopBLSE2020test$n,decreasing=TRUE),]
#Génération d'une table au format txt à partir de cette dataframe
write.table(TenTopBLSE2020test,sep = ",", quote = FALSE, row.names = F)
#Possible de lire cette table dans Word ou Excel
#Calcul de la fréquence de test de chaque atb
BLSE2020testN <- BLSE2020test %>% mutate_all(as.numeric)
PropBLSE2020<-colSums(BLSE2020testN)/dim(BLSE2020test)[1]

write.table(PropBLSE2020,sep = ",", quote = FALSE, row.names = T)

#test 2021--------
data2021<-read.csv(file="EC-2021.csv",sep=";")

BLSE2021 <- data2021
unique(BLSE2021$phenotype) #on choisit les phénotypes dit BLSE, BLSE+IMPER (4) et BLSE+CARBA (2)

BLSE2021 <- BLSE2021[BLSE2021$phenotype == "BLSE",]

unique(BLSE2021$phenotype)

BLSE2021 <- subset(BLSE2021, select = -c(AMC,AMC_urine, AMC03, TCC, TZP, C3G, FQ)) #supprimer les colonnnes pas pertinentes 

BLSE2021 <- BLSE2021[,11:42]

#test antibiotique 
BLSE2021test <- BLSE2021

BLSE2021test[BLSE2021test == ''] <- 0
BLSE2021test[BLSE2021test == 'S'] <- 1
BLSE2021test[BLSE2021test == 'I'] <- 1
BLSE2021test[BLSE2021test == 'R'] <- 1

SeqBLSE2021test<-unique(BLSE2021test)
dim(SeqBLSE2021test)[1]
#Ajout d'une colonne indiquant la fréquence de chaque jeu d'atb testé
FreqBLSE2021test<-BLSE2021test %>% group_by_all() %>% count() %>% ungroup()
#Dix jeux d'antibiotiques testés le plus fréquemment
TenTopBLSE2021test<-FreqBLSE2021test %>%
  arrange(desc(n)) %>%
  slice_head(n=10)

#Proportion des isolats testés concernés par ces dix jeux d'antibiotiques
sum(TenTopBLSE2021test[,33])/dim(BLSE2021test)[1]
#Dataframe décrivant ces dix jeux, dans l'ordre de fréquence décroissante
TenTopBLSE2021test <- TenTopBLSE2021test[order(TenTopBLSE2021test$n,decreasing=TRUE),]
#Génération d'une table au format txt à partir de cette dataframe
write.table(TenTopBLSE2021test,sep = ",", quote = FALSE, row.names = F)
#Possible de lire cette table dans Word ou Excel
#Calcul de la fréquence de test de chaque atb
BLSE2021testN <- BLSE2021test %>% mutate_all(as.numeric)
PropBLSE2021<-colSums(BLSE2021testN)/dim(BLSE2021test)[1]

write.table(PropBLSE2021,sep = ",", quote = FALSE, row.names = T)


#res 2018-------
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


BLSE2018_res_log <- BLSE2018_res %>% mutate_all(as.logical)  #forme logique pour pouvoir executer les fonctions
str(BLSE2018_res_log)
antibio_BLSE2018 <- anti_bio(BLSE2018_res_log,1:30)
antibio_BLSE2018 <- as.data.frame(antibio_BLSE2018)

#calcul du nombre de res antibiotiques de chaque isolat

BLSE2018_res_num <- BLSE2018

BLSE2018_res_num[BLSE2018_res_num == ""] <- NA
BLSE2018_res_num[BLSE2018_res_num == "S"] <- 0
BLSE2018_res_num[BLSE2018_res_num == "I"] <- 0
BLSE2018_res_num[BLSE2018_res_num == "R"] <- 1


FreqBLSE2018 <-BLSE2018_res_num %>% group_by_all() %>% count() %>% ungroup()

uniqfreqBLSE2018 <- unique(FreqBLSE2018)

TenTopBLSE2018<-uniqfreqBLSE2018 %>%
  arrange(desc(n)) %>%
  slice_head(n=10)

TenTopBLSE2018 <- TenTopBLSE2018[order(TenTopBLSE2018$n,decreasing=TRUE),]


#montrer les atbs ayant le plus de souches testes resistantes 

#changer les colonnes et les lignes 
final_antibio_BLSE2018 <- as.data.frame(t(antibio_BLSE2018))
top_mean_BLSE2018 <- top_n(final_antibio_BLSE2018, 10, eSupport) 

top_mean_BLSE2018 %>%  #graphe top 10 resistances
  ggplot(., aes(x=row.names(.), y= eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Les 10 antibiotiques ayant la plus grande prevalence de souches testees 'BLSE' resistantes en 2018") +
  guides(x = guide_axis(angle = 30))

#distribution des resistances d'atb par isolats (a modifier)
final_antibio_BLSE2018 %>% 
  ggplot(., aes(reorder(row.names(.), -eSupport, sum), eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Fréquence des résistances antibiotiques des BLSE en 2018") +
  guides(x = guide_axis(angle = 60)) +
  theme(axis.text = element_text(face = "bold"))
ggsave("freq_res_BLSE18.png")

#profils de resistances
tot_profils_BLSE2018 <- res_profile(BLSE2018_res_log, 1:30)
uniq_profils_BLSE2018 <-unique(tot_profils_BLSE2018)


freq_profils_BLSE2018 <- tot_profils_BLSE2018 %>% group_by_all() %>% count() %>% ungroup()



# changement de la colonne freq pour avoir un % du nb total isolats 

sum(freq_profils_BLSE2018$n)
freq_profils_BLSE2018_perc <- freq_profils_BLSE2018
freq_profils_BLSE2018_perc[,3] <- round(freq_profils_BLSE2018[,3]/sum(freq_profils_BLSE2018$n),5)*100

top10_freq_profils_BLSE2018 <- top_n(freq_profils_BLSE2018_perc, 10, n)
sum(top10_freq_profils_BLSE2018[,3]) #correspond à % des resistances 
top10_freq_profils_BLSE2018 %>%  #graphe profil res
  ggplot(., aes(reorder(Profile, -n, sum),n)) + 
  geom_col() +
  labs(x = "Profils de resistance", y = "Frequence(%)", title = "Les 10 profils de resistances des BLSE les plus frequents en 2018") +
  guides(x = guide_axis(angle = 30)) +
  theme(axis.text = element_text(face = "bold"))
ggsave("top10res_BLSE2018.png")

#distribution du nb d'atb presents dans les multiresitances

freq_profils_BLSE2018_perc %>%
  ggplot(., aes(x=NombreRes, y= n)) +
  geom_col() +
  labs(x = "Nombre de résistances", y = "Prevalence (%)", title = "Distribution du nombre de résistances antibiotiques chez les BLSE en 2018")

write.table(top10_freq_profils_BLSE2018, sep = ",", quote = FALSE, row.names = F)

write.table(final_antibio_BLSE2018, sep = ",", quote = FALSE, row.names = T)

#res 2019---------
data2019<-read.csv(file="EC-2019.csv",sep=";")


BLSE2019 <- data2019
unique(BLSE2019$phenotype) #on choisit les phénotypes dit BLSE, BLSE+IMPER (4) et BLSE+CARBA (2)

BLSE2019 <- filter(BLSE2019, type.prelevement == "URINES" & phenotype == "BLSE" | phenotype == "BLSE+CARBAPENEMASE (OXA-48)" | phenotype == "BLSE+IMPER" | phenotype == "BlSE" | phenotype == "bLSE" | phenotype == "blse" | phenotype == "BLSE+CASEHN" | phenotype == "BLSE (PAS DE CARBA)" | phenotype == "BLSE+imper" )

unique(BLSE2019$phenotype)
unique(BLSE2019$type.prelevement)


BLSE2019 <- subset(BLSE2019, select = -c(TCC, CN, AZT, CS, TGC, C3G, FQ, AMX, AMP)) #supprimer les colonnnes pas pertinentes 

BLSE2019_res <- BLSE2019[,11:40]

BLSE2019_res[BLSE2019_res == ''] <- NA
BLSE2019_res[BLSE2019_res == 'S'] <- FALSE
BLSE2019_res[BLSE2019_res == 'I'] <- FALSE
BLSE2019_res[BLSE2019_res == 'R'] <- TRUE


BLSE2019_res_log <- BLSE2019_res %>% mutate_all(as.logical)  #forme logique pour pouvoir executer les fonctions
str(BLSE2019_res_log)
antibio_BLSE2019 <- anti_bio(BLSE2019_res_log,1:30)
antibio_BLSE2019 <- as.data.frame(antibio_BLSE2019)

#calcul du nombre de res antibiotiques de chaque isolat

BLSE2019_res_num <- BLSE2019[,11:40]

BLSE2019_res_num[BLSE2019_res_num == ""] <- NA
BLSE2019_res_num[BLSE2019_res_num == "S"] <- 0
BLSE2019_res_num[BLSE2019_res_num == "I"] <- 0
BLSE2019_res_num[BLSE2019_res_num == "R"] <- 1


FreqBLSE2019 <-BLSE2019_res_num %>% group_by_all() %>% count() %>% ungroup()

uniqfreqBLSE2019 <- unique(FreqBLSE2019)

TenTopBLSE2019<-uniqfreqBLSE2019 %>%
  arrange(desc(n)) %>%
  slice_head(n=10)

TenTopBLSE2019 <- TenTopBLSE2019[order(TenTopBLSE2019$n,decreasing=TRUE),]


#montrer les atbs ayant le plus de souches testes resistantes 

#changer les colonnes et les lignes 
final_antibio_BLSE2019 <- as.data.frame(t(antibio_BLSE2019))
top_mean_BLSE2019 <- top_n(final_antibio_BLSE2019, 10, eSupport) 

top_mean_BLSE2019 %>%  #graphe top 10 resistances
  ggplot(., aes(x=row.names(.), y= eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Les 10 antibiotiques ayant la plus grande prevalence de souches testees 'BLSE' resistantes en 2019") +
  guides(x = guide_axis(angle = 30))

#distribution des resistances d'atb par isolats (a modifier)
final_antibio_BLSE2019 %>% 
  ggplot(., aes(reorder(row.names(.), -eSupport,sum), eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Fréquence des résistances antibiotiques des BLSE en 2019") +
  guides(x = guide_axis(angle = 60))+
  theme(axis.text = element_text(face = "bold"))
ggsave("freq_res_BLSE19.png")

#profils de resistances
tot_profils_BLSE2019 <- res_profile(BLSE2019_res_log, 1:30)
uniq_profils_BLSE2019 <-unique(tot_profils_BLSE2019)


freq_profils_BLSE2019 <- tot_profils_BLSE2019 %>% group_by_all() %>% count() %>% ungroup()



# changement de la colonne freq pour avoir un % du nb total isolats 

sum(freq_profils_BLSE2019$n)
freq_profils_BLSE2019_perc <- freq_profils_BLSE2019
freq_profils_BLSE2019_perc[,3] <- round(freq_profils_BLSE2019[,3]/sum(freq_profils_BLSE2019$n),5)*100

top10_freq_profils_BLSE2019 <- top_n(freq_profils_BLSE2019_perc, 10, n)
sum(top10_freq_profils_BLSE2019[,3]) #correspond à % des resistances 
top10_freq_profils_BLSE2019 %>%  #graphe profil res
  ggplot(., aes(reorder(Profile, -n,sum), n)) + 
  geom_col() +
  labs(x = "Profils de resistance", y = "Frequence(%)", title = "Les 10 profils de resistances des BLSE les plus frequents en 2019") +
  guides(x = guide_axis(angle = 30))+
  theme(axis.text = element_text(face = "bold"))
ggsave("top10res_BLSE19.png")

#distribution du nb d'atb presents dans les multiresitances

freq_profils_BLSE2019_perc %>%
  ggplot(., aes(x=NombreRes, y= n)) +
  geom_col() +
  labs(x = "Nombre de résistances", y = "Prevalence (%)", title = "Distribution du nombre de résistances antibiotiques chez les BLSE en 2019")

write.table(top10_freq_profils_BLSE2019, sep = ",", quote = FALSE, row.names = F)

write.table(final_antibio_BLSE2019, sep = ",", quote = FALSE, row.names = T)


#res 2020------

data2020<-read.csv(file="EC-2020.csv",sep=";")


BLSE2020 <- data2020
unique(BLSE2020$phenotype) #on choisit les phénotypes dit BLSE, BLSE+IMPER (4) et BLSE+CARBA (2)

BLSE2020 <- filter(BLSE2020, type.prelevement == "URINES" & phenotype == "BLSE" | phenotype == "BLSE+CARBAPENEMASE (OXA-48)" | phenotype == "BLSE+IMPER" | phenotype == "BLSE+CASEHN" | phenotype == "bLSE" | phenotype == "blse" | phenotype == "BLSE+CASE+CARBAPENEMASE (OXA-48)" | phenotype == "BLSE+carba" | phenotype == "BLSE+NDM" | phenotype ==  "BLSE+CARBA"| phenotype == "CASE+BLSE" | phenotype == "BLSE+CASE" | phenotype == "BLSE+CARBAPENEMASE (OXA48)")

unique(BLSE2020$phenotype)
unique(BLSE2020$type.prelevement)

BLSE2020<- BLSE2020[-which(BLSE2020$type.prelevement == "PUS-SUPERFICIEL (AUTRES)"),]
BLSE2020<- BLSE2020[-which(BLSE2020$type.prelevement == "HEMOCULTURE"),]

BLSE2020 <- subset(BLSE2020, select = -c(TCC, CN, AZT, CS, TGC, C3G, FQ, AMX, AMP)) #supprimer les colonnnes pas pertinentes 

BLSE2020_res <- BLSE2020[,11:40]

BLSE2020_res[BLSE2020_res == ''] <- NA
BLSE2020_res[BLSE2020_res == 'S'] <- FALSE
BLSE2020_res[BLSE2020_res == 'I'] <- FALSE
BLSE2020_res[BLSE2020_res == 'R'] <- TRUE


BLSE2020_res_log <- BLSE2020_res %>% mutate_all(as.logical)  #forme logique pour pouvoir executer les fonctions
str(BLSE2020_res_log)
antibio_BLSE2020 <- anti_bio(BLSE2020_res_log,1:30)
antibio_BLSE2020 <- as.data.frame(antibio_BLSE2020)

#calcul du nombre de res antibiotiques de chaque isolat

BLSE2020_res_num <- BLSE2020[,11:40]

BLSE2020_res_num[BLSE2020_res_num == ""] <- NA
BLSE2020_res_num[BLSE2020_res_num == "S"] <- 0
BLSE2020_res_num[BLSE2020_res_num == "I"] <- 0
BLSE2020_res_num[BLSE2020_res_num == "R"] <- 1


FreqBLSE2020 <-BLSE2020_res_num %>% group_by_all() %>% count() %>% ungroup()

uniqfreqBLSE2020 <- unique(FreqBLSE2020)

TenTopBLSE2020<-uniqfreqBLSE2020 %>%
  arrange(desc(n)) %>%
  slice_head(n=10)

TenTopBLSE2020 <- TenTopBLSE2020[order(TenTopBLSE2020$n,decreasing=TRUE),]


#montrer les atbs ayant le plus de souches testes resistantes 

#changer les colonnes et les lignes 
final_antibio_BLSE2020 <- as.data.frame(t(antibio_BLSE2020))
top_mean_BLSE2020 <- top_n(final_antibio_BLSE2020, 10, eSupport) 

top_mean_BLSE2020 %>%  #graphe top 10 resistances
  ggplot(., aes(x=row.names(.), y= eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Les 10 antibiotiques ayant la plus grande prevalence de souches testees 'BLSE' resistantes en 2020") +
  guides(x = guide_axis(angle = 30))

#distribution des resistances d'atb par isolats (a modifier)
final_antibio_BLSE2020 %>% 
  ggplot(., aes(reorder(row.names(.), -eSupport, sum), eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Fréquence des résistances antibiotiques des BLSE en 2020") +
  guides(x = guide_axis(angle = 60))+
  theme(axis.text = element_text(face = "bold"))
ggsave("freq_res_BLSE20.png")

#profils de resistances
tot_profils_BLSE2020 <- res_profile(BLSE2020_res_log, 1:30)
uniq_profils_BLSE2020 <-unique(tot_profils_BLSE2020)


freq_profils_BLSE2020 <- tot_profils_BLSE2020 %>% group_by_all() %>% count() %>% ungroup()



# changement de la colonne freq pour avoir un % du nb total isolats 

sum(freq_profils_BLSE2020$n)
freq_profils_BLSE2020_perc <- freq_profils_BLSE2020
freq_profils_BLSE2020_perc[,3] <- round(freq_profils_BLSE2020[,3]/sum(freq_profils_BLSE2020$n),5)*100

top10_freq_profils_BLSE2020 <- top_n(freq_profils_BLSE2020_perc, 10, n)
sum(top10_freq_profils_BLSE2020[,3]) #correspond à % des resistances 
top10_freq_profils_BLSE2020 %>%  #graphe profil res
  ggplot(., aes(reorder(Profile, -n,sum), n)) + 
  geom_col() +
  labs(x = "Profils de resistance", y = "Frequence(%)", title = "Les 10 profils de resistances des BLSE les plus frequents en 2020") +
  guides(x = guide_axis(angle = 30))+
  theme(axis.text = element_text(face = "bold"))
ggsave("top10res_BLSE20.png")

#distribution du nb d'atb presents dans les multiresitances

freq_profils_BLSE2020_perc %>%
  ggplot(., aes(x=NombreRes, y= n)) +
  geom_col() +
  labs(x = "Nombre de résistances", y = "Prevalence (%)", title = "Distribution du nombre de résistances antibiotiques chez les BLSE en 2020")

write.table(top10_freq_profils_BLSE2020, sep = ",", quote = FALSE, row.names = F)

write.table(final_antibio_BLSE2020, sep = ",", quote = FALSE, row.names = T)

#res 2021--------

data2021<-read.csv(file="EC-2021.csv",sep=";")


BLSE2021 <- data2021
unique(BLSE2021$phenotype) #on choisit les phénotypes dit BLSE, BLSE+IMPER (4) et BLSE+CARBA (2)

BLSE2021 <- filter(BLSE2021, type.prelevement == "URINES" & phenotype == "BLSE+CARBAPENEMASE (NDM)" | phenotype == "BLSE" | phenotype == "BLSE+CASE" | phenotype == "BLSE+CASE+IMPER" | phenotype == "BLSE+IMPER" | phenotype == "BLSE+CARBAPENEMASE (OXA-48)" | phenotype == "blse" | phenotype == "BLSE+CASEHN")
BLSE2021 <- filter(BLSE2021, type.prelevement == "URINES")

unique(BLSE2021$phenotype)
unique(BLSE2021$type.prelevement)

BLSE2021 <- subset(BLSE2021, select = -c(TCC, CN, AZT, CS, TGC, C3G, FQ, AMX, AMP)) #supprimer les colonnnes pas pertinentes 

BLSE2021_res <- BLSE2021[,11:40]

BLSE2021_res[BLSE2021_res == ''] <- NA
BLSE2021_res[BLSE2021_res == 'S'] <- FALSE
BLSE2021_res[BLSE2021_res == 'I'] <- FALSE
BLSE2021_res[BLSE2021_res == 'R'] <- TRUE


BLSE2021_res_log <- BLSE2021_res %>% mutate_all(as.logical)  #forme logique pour pouvoir executer les fonctions
str(BLSE2021_res_log)
antibio_BLSE2021 <- anti_bio(BLSE2021_res_log,1:30)
antibio_BLSE2021 <- as.data.frame(antibio_BLSE2021)

#calcul du nombre de res antibiotiques de chaque isolat

BLSE2021_res_num <- BLSE2021[,11:40]

BLSE2021_res_num[BLSE2021_res_num == ""] <- NA
BLSE2021_res_num[BLSE2021_res_num == "S"] <- 0
BLSE2021_res_num[BLSE2021_res_num == "I"] <- 0
BLSE2021_res_num[BLSE2021_res_num == "R"] <- 1


FreqBLSE2021 <-BLSE2021_res_num %>% group_by_all() %>% count() %>% ungroup()

uniqfreqBLSE2021 <- unique(FreqBLSE2021)

TenTopBLSE2021<-uniqfreqBLSE2021 %>%
  arrange(desc(n)) %>%
  slice_head(n=10)

TenTopBLSE2021 <- TenTopBLSE2021[order(TenTopBLSE2021$n,decreasing=TRUE),]


#montrer les atbs ayant le plus de souches testes resistantes 

#changer les colonnes et les lignes 
final_antibio_BLSE2021 <- as.data.frame(t(antibio_BLSE2021))
top_mean_BLSE2021 <- top_n(final_antibio_BLSE2021, 10, eSupport) 

top_mean_BLSE2021 %>%  #graphe top 10 resistances
  ggplot(., aes(x=row.names(.), y= eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Les 10 antibiotiques ayant la plus grande prevalence de souches testees 'BLSE' resistantes en 2021") +
  guides(x = guide_axis(angle = 30))

#distribution des resistances d'atb par isolats (a modifier)
final_antibio_BLSE2021 %>% 
  ggplot(., aes(reorder(row.names(.), -eSupport, sum),eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Fréquence des résistances antibiotiques des BLSE en 2021") +
  guides(x = guide_axis(angle = 60))+
  theme(axis.text = element_text(face = "bold"))
ggsave("freq_res_BLSE21.png")

#profils de resistances
tot_profils_BLSE2021 <- res_profile(BLSE2021_res_log, 1:30)
uniq_profils_BLSE2021 <-unique(tot_profils_BLSE2021)


freq_profils_BLSE2021 <- tot_profils_BLSE2021 %>% group_by_all() %>% count() %>% ungroup()



# changement de la colonne freq pour avoir un % du nb total isolats 

sum(freq_profils_BLSE2021$n)
freq_profils_BLSE2021_perc <- freq_profils_BLSE2021
freq_profils_BLSE2021_perc[,3] <- round(freq_profils_BLSE2021[,3]/sum(freq_profils_BLSE2021$n),5)*100

top10_freq_profils_BLSE2021 <- top_n(freq_profils_BLSE2021_perc, 10, n)
sum(top10_freq_profils_BLSE2021[,3]) #correspond à % des resistances 
top10_freq_profils_BLSE2021 %>%  #graphe profil res
  ggplot(., aes(reorder(Profile, -n, sum),n)) + 
  geom_col() +
  labs(x = "Profils de resistance", y = "Frequence(%)", title = "Les 10 profils de resistances des BLSE les plus frequents en 2021") +
  guides(x = guide_axis(angle = 30)) +
  theme(axis.text = element_text(face = "bold"))
ggsave("top10res_BLSE21.png", scale = 1, width = 24, height = 12, units = "in")

#distribution du nb d'atb presents dans les multiresitances

freq_profils_BLSE2021_perc %>%
  ggplot(., aes(x=NombreRes, y= n)) +
  geom_col() +
  labs(x = "Nombre de résistances", y = "Prevalence (%)", title = "Distribution du nombre de résistances antibiotiques chez les BLSE en 2021")

write.table(top10_freq_profils_BLSE2021, sep = ",", quote = FALSE, row.names = F)

write.table(final_antibio_BLSE2021, sep = ",", quote = FALSE, row.names = T)


#sans BLSE2018------------------

data2018<-read.csv(file="EC-2018.csv",sep=";")


sBLSE2018 <- data2018
unique(sBLSE2018$phenotype) #on choisit les phénotypes dit BLSE, BLSE+IMPER (4) et BLSE+CARBA (2)


sBLSE2018 <- sBLSE2018[-which(sBLSE2018$phenotype == "BLSE" | sBLSE2018$phenotype == "BLSE+CARBA" | sBLSE2018$phenotype == "BLSE+IMPER"),]

sBLSE2018 <- filter(sBLSE2018, type.prelevement == "URINES")

unique(sBLSE2018$phenotype)
unique(sBLSE2018$type.prelevement)



sBLSE2018 <- subset(sBLSE2018, select = -c(TCC, CN, AZT, CS, TGC, C3G, FQ)) #supprimer les colonnnes pas pertinentes 

sBLSE2018_res <- sBLSE2018[,11:42]

sBLSE2018_res[sBLSE2018_res == ''] <- NA
sBLSE2018_res[sBLSE2018_res == 'S'] <- FALSE
sBLSE2018_res[sBLSE2018_res == 'I'] <- FALSE
sBLSE2018_res[sBLSE2018_res == 'R'] <- TRUE


sBLSE2018_res_log <- sBLSE2018_res %>% mutate_all(as.logical)  #forme logique pour pouvoir executer les fonctions
str(sBLSE2018_res_log)
antibio_sBLSE2018 <- anti_bio(sBLSE2018_res_log,1:32)
antibio_sBLSE2018 <- as.data.frame(antibio_sBLSE2018)

#calcul du nombre de res antibiotiques de chaque isolat

sBLSE2018_res_num <- sBLSE2018

sBLSE2018_res_num[sBLSE2018_res_num == ""] <- NA
sBLSE2018_res_num[sBLSE2018_res_num == "S"] <- 0
sBLSE2018_res_num[sBLSE2018_res_num == "I"] <- 0
sBLSE2018_res_num[sBLSE2018_res_num == "R"] <- 1


FreqsBLSE2018 <-sBLSE2018_res_num %>% group_by_all() %>% count() %>% ungroup()

uniqfreqsBLSE2018 <- unique(FreqsBLSE2018)

TenTopsBLSE2018<-uniqfreqsBLSE2018 %>%
  arrange(desc(n)) %>%
  slice_head(n=10)

TenTopsBLSE2018 <- TenTopsBLSE2018[order(TenTopsBLSE2018$n,decreasing=TRUE),]


#montrer les atbs ayant le plus de souches testes resistantes 

#changer les colonnes et les lignes 
final_antibio_sBLSE2018 <- as.data.frame(t(antibio_sBLSE2018))
top_mean_sBLSE2018 <- top_n(final_antibio_sBLSE2018, 10, eSupport) 

top_mean_sBLSE2018 %>%  #graphe top 10 resistances
  ggplot(., aes(x=row.names(.), y= eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Les 10 antibiotiques ayant la plus grande prevalence de souches testees hors BLSE resistantes en 2018") +
  guides(x = guide_axis(angle = 30))

#distribution des resistances d'atb par isolats (a modifier)
final_antibio_sBLSE2018 %>% 
  ggplot(., aes(reorder(row.names(.), -eSupport,sum),eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Fréquence des résistances antibiotiques des souches hors BLSE en 2018") +
  guides(x = guide_axis(angle = 60))+
  theme(axis.text = element_text(face = "bold"))
ggsave("freq_res_sBLSE18.png")

#profils de resistances
tot_profils_sBLSE2018 <- res_profile(sBLSE2018_res_log, 1:32)
uniq_profils_sBLSE2018 <-unique(tot_profils_sBLSE2018)


freq_profils_sBLSE2018 <- tot_profils_sBLSE2018 %>% group_by_all() %>% count() %>% ungroup()



# changement de la colonne freq pour avoir un % du nb total isolats 

sum(freq_profils_sBLSE2018$n)
freq_profils_sBLSE2018_perc <- freq_profils_sBLSE2018
freq_profils_sBLSE2018_perc[,3] <- round(freq_profils_sBLSE2018[,3]/sum(freq_profils_sBLSE2018$n),5)*100

top10_freq_profils_sBLSE2018 <- top_n(freq_profils_sBLSE2018_perc, 10, n)
sum(top10_freq_profils_sBLSE2018[,3]) #correspond à % des resistances 
top10_freq_profils_sBLSE2018 %>%  #graphe profil res
  ggplot(., aes(reorder(Profile, -n, sum), n)) + 
  geom_col() +
  labs(x = "Profils de resistance", y = "Frequence(%)", title = "Les 10 profils de resistances des souches hors BLSE les plus frequents en 2018") +
  guides(x = guide_axis(angle = 30))+
  theme(axis.text = element_text(face = "bold")) 
ggsave("top10res_sBLSE18.png")

#distribution du nb d'atb presents dans les multiresitances

freq_profils_sBLSE2018_perc %>%
  ggplot(., aes(x=NombreRes, y= n)) +
  geom_col() +
  labs(x = "Nombre de résistances", y = "Prevalence (%)", title = "Distribution du nombre de résistances antibiotiques des souches hors BLSE en 2018")

write.table(top10_freq_profils_sBLSE2018, sep = ",", quote = FALSE, row.names = F)

write.table(final_antibio_sBLSE2018, sep = ",", quote = FALSE, row.names = T)

#sans BLSE19--------------

data2019<-read.csv(file="EC-2019.csv",sep=";")


sBLSE2019 <- data2019
unique(sBLSE2019$phenotype) #on choisit les phénotypes dit BLSE, BLSE+IMPER (4) et BLSE+CARBA (2)

sBLSE2019 <- filter(sBLSE2019, type.prelevement == "URINES")
sBLSE2019 <- sBLSE2019[-which(sBLSE2019$phenotype == "BLSE" | sBLSE2019$phenotype == "BLSE+CARBAPENEMASE (OXA-48)" | sBLSE2019$phenotype == "BLSE+IMPER" | sBLSE2019$phenotype == "BlSE" | sBLSE2019$phenotype == "bLSE" | sBLSE2019$phenotype == "blse" | sBLSE2019$phenotype == "BLSE+CASEHN" | sBLSE2019$phenotype == "BLSE (PAS DE CARBA)" | sBLSE2019$phenotype == "BLSE+imper"),]

unique(sBLSE2019$phenotype)
unique(BLSE2019$type.prelevement)

sBLSE2019 <- subset(sBLSE2019, select = -c(TCC, CN, AZT, CS, TGC, C3G, FQ)) #supprimer les colonnnes pas pertinentes 

sBLSE2019_res <- sBLSE2019[,11:42]

sBLSE2019_res[sBLSE2019_res == ''] <- NA
sBLSE2019_res[sBLSE2019_res == 'S'] <- FALSE
sBLSE2019_res[sBLSE2019_res == 'I'] <- FALSE
sBLSE2019_res[sBLSE2019_res == 'R'] <- TRUE


sBLSE2019_res_log <- sBLSE2019_res %>% mutate_all(as.logical)  #forme logique pour pouvoir executer les fonctions
str(sBLSE2019_res_log)
antibio_sBLSE2019 <- anti_bio(sBLSE2019_res_log,1:32)
antibio_sBLSE2019 <- as.data.frame(antibio_sBLSE2019)

#calcul du nombre de res antibiotiques de chaque isolat

sBLSE2019_res_num <- sBLSE2019

sBLSE2019_res_num[sBLSE2019_res_num == ""] <- NA
sBLSE2019_res_num[sBLSE2019_res_num == "S"] <- 0
sBLSE2019_res_num[sBLSE2019_res_num == "I"] <- 0
sBLSE2019_res_num[sBLSE2019_res_num == "R"] <- 1


FreqsBLSE2019 <-sBLSE2019_res_num %>% group_by_all() %>% count() %>% ungroup()

uniqfreqsBLSE2019 <- unique(FreqsBLSE2019)

TenTopsBLSE2019<-uniqfreqsBLSE2019 %>%
  arrange(desc(n)) %>%
  slice_head(n=10)

TenTopsBLSE2019 <- TenTopsBLSE2019[order(TenTopsBLSE2019$n,decreasing=TRUE),]


#montrer les atbs ayant le plus de souches testes resistantes 

#changer les colonnes et les lignes 
final_antibio_sBLSE2019 <- as.data.frame(t(antibio_sBLSE2019))
top_mean_sBLSE2019 <- top_n(final_antibio_sBLSE2019, 10, eSupport) 

top_mean_sBLSE2019 %>%  #graphe top 10 resistances
  ggplot(., aes(x=row.names(.), y= eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Les 10 antibiotiques ayant la plus grande prevalence de souches testees hors BLSE resistantes en 2019") +
  guides(x = guide_axis(angle = 30))

#distribution des resistances d'atb par isolats (a modifier)
final_antibio_sBLSE2019 %>% 
  ggplot(., aes(reorder(row.names(.), -eSupport,sum), eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Fréquence des résistances antibiotiques des souches hors BLSE en 2019") +
  guides(x = guide_axis(angle = 60))+
  theme(axis.text = element_text(face = "bold"))
ggsave("freq_res_sBLSE19.png")

#profils de resistances
tot_profils_sBLSE2019 <- res_profile(sBLSE2019_res_log, 1:32)
uniq_profils_sBLSE2019 <-unique(tot_profils_sBLSE2019)


freq_profils_sBLSE2019 <- tot_profils_sBLSE2019 %>% group_by_all() %>% count() %>% ungroup()



# changement de la colonne freq pour avoir un % du nb total isolats 

sum(freq_profils_sBLSE2019$n)
freq_profils_sBLSE2019_perc <- freq_profils_sBLSE2019
freq_profils_sBLSE2019_perc[,3] <- round(freq_profils_sBLSE2019[,3]/sum(freq_profils_sBLSE2019$n),5)*100

top10_freq_profils_sBLSE2019 <- top_n(freq_profils_sBLSE2019_perc, 10, n)
sum(top10_freq_profils_sBLSE2019[,3]) #correspond à % des resistances 
top10_freq_profils_sBLSE2019 %>%  #graphe profil res
  ggplot(., aes(reorder(Profile, -n, sum), n)) + 
  geom_col() +
  labs(x = "Profils de resistance", y = "Frequence(%)", title = "Les 10 profils de resistances des souches hors BLSE les plus frequents en 2019") +
  guides(x = guide_axis(angle = 30))+
  theme(axis.text = element_text(face = "bold"))
ggsave("top10res_sBLSE19.png")

#distribution du nb d'atb presents dans les multiresitances

freq_profils_sBLSE2019_perc %>%
  ggplot(., aes(x=NombreRes, y= n)) +
  geom_col() +
  labs(x = "Nombre de résistances", y = "Prevalence (%)", title = "Distribution du nombre de résistances antibiotiques des souches hors BLSE en 2019")

write.table(top10_freq_profils_sBLSE2019, sep = ",", quote = FALSE, row.names = F)

write.table(final_antibio_sBLSE2019, sep = ",", quote = FALSE, row.names = T)

#sans BLSE20------------------

data2020<-read.csv(file="EC-2020.csv",sep=";")


sBLSE2020 <- data2020
unique(sBLSE2020$phenotype) #on choisit les phénotypes dit BLSE, BLSE+IMPER (4) et BLSE+CARBA (2)

sBLSE2020 <- filter(sBLSE2020, type.prelevement == "URINES")
sBLSE2020 <- sBLSE2020[-which(sBLSE2020$phenotype == "BLSE" | sBLSE2020$phenotype == "BLSE+CARBAPENEMASE (OXA-48)" | sBLSE2020$phenotype == "BLSE+IMPER" | sBLSE2020$phenotype == "BLSE+CASEHN" | sBLSE2020$phenotype == "bLSE" | sBLSE2020$phenotype == "blse" | sBLSE2020$phenotype == "BLSE+CASE+CARBAPENEMASE (OXA-48)" | sBLSE2020$phenotype == "BLSE+carba" | sBLSE2020$phenotype == "BLSE+NDM" | sBLSE2020$phenotype ==  "BLSE+CARBA"| sBLSE2020$phenotype == "CASE+BLSE" | sBLSE2020$phenotype == "BLSE+CASE" | sBLSE2020$phenotype == "BLSE+CARBAPENEMASE (OXA48)"),]

unique(sBLSE2020$phenotype)
unique(sBLSE2020$type.prelevement)

sBLSE2020 <- subset(sBLSE2020, select = -c(TCC, CN, AZT, CS, TGC, C3G, FQ)) #supprimer les colonnnes pas pertinentes 

sBLSE2020_res <- sBLSE2020[,11:42]

sBLSE2020_res[sBLSE2020_res == ''] <- NA
sBLSE2020_res[sBLSE2020_res == 'S'] <- FALSE
sBLSE2020_res[sBLSE2020_res == 'I'] <- FALSE
sBLSE2020_res[sBLSE2020_res == 'R'] <- TRUE


sBLSE2020_res_log <- sBLSE2020_res %>% mutate_all(as.logical)  #forme logique pour pouvoir executer les fonctions
str(sBLSE2020_res_log)
antibio_sBLSE2020 <- anti_bio(sBLSE2020_res_log,1:32)
antibio_sBLSE2020 <- as.data.frame(antibio_sBLSE2020)

#calcul du nombre de res antibiotiques de chaque isolat

sBLSE2020_res_num <- sBLSE2020

sBLSE2020_res_num[sBLSE2020_res_num == ""] <- NA
sBLSE2020_res_num[sBLSE2020_res_num == "S"] <- 0
sBLSE2020_res_num[sBLSE2020_res_num == "I"] <- 0
sBLSE2020_res_num[sBLSE2020_res_num == "R"] <- 1


FreqsBLSE2020 <-sBLSE2020_res_num %>% group_by_all() %>% count() %>% ungroup()

uniqfreqsBLSE2020 <- unique(FreqsBLSE2020)

TenTopsBLSE2020<-uniqfreqsBLSE2020 %>%
  arrange(desc(n)) %>%
  slice_head(n=10)

TenTopsBLSE2020 <- TenTopsBLSE2020[order(TenTopsBLSE2020$n,decreasing=TRUE),]


#montrer les atbs ayant le plus de souches testes resistantes 

#changer les colonnes et les lignes 
final_antibio_sBLSE2020 <- as.data.frame(t(antibio_sBLSE2020))
top_mean_sBLSE2020 <- top_n(final_antibio_sBLSE2020, 10, eSupport) 

top_mean_sBLSE2020 %>%  #graphe top 10 resistances
  ggplot(., aes(x=row.names(.), y= eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Les 10 antibiotiques ayant la plus grande prevalence de souches testees hors BLSE resistantes en 2020") +
  guides(x = guide_axis(angle = 30))

#distribution des resistances d'atb par isolats (a modifier)
final_antibio_sBLSE2020 %>% 
  ggplot(., aes(reorder(row.names(.), -eSupport,sum), eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Fréquence des résistances antibiotiques des souches hors BLSE en 2020") +
  guides(x = guide_axis(angle = 60))+
  theme(axis.text = element_text(face = "bold"))
ggsave("freq_res_sBLSE20.png")

#profils de resistances
tot_profils_sBLSE2020 <- res_profile(sBLSE2020_res_log, 1:32)
uniq_profils_sBLSE2020 <-unique(tot_profils_sBLSE2020)


freq_profils_sBLSE2020 <- tot_profils_sBLSE2020 %>% group_by_all() %>% count() %>% ungroup()



# changement de la colonne freq pour avoir un % du nb total isolats 

sum(freq_profils_sBLSE2020$n)
freq_profils_sBLSE2020_perc <- freq_profils_sBLSE2020
freq_profils_sBLSE2020_perc[,3] <- round(freq_profils_sBLSE2020[,3]/sum(freq_profils_sBLSE2020$n),5)*100

top10_freq_profils_sBLSE2020 <- top_n(freq_profils_sBLSE2020_perc, 10, n)
sum(top10_freq_profils_sBLSE2020[,3]) #correspond à % des resistances 
top10_freq_profils_sBLSE2020 %>%  #graphe profil res
  ggplot(., aes(reorder(Profile, -n, sum), n)) + 
  geom_col() +
  labs(x = "Profils de resistance", y = "Frequence(%)", title = "Les 10 profils de resistances des souches hors BLSE les plus frequents en 2020") +
  guides(x = guide_axis(angle = 30))+
  theme(axis.text = element_text(face = "bold"))
ggsave("top10res_sBLSE20.png")

#distribution du nb d'atb presents dans les multiresitances

freq_profils_sBLSE2020_perc %>%
  ggplot(., aes(x=NombreRes, y= n)) +
  geom_col() +
  labs(x = "Nombre de résistances", y = "Prevalence (%)", title = "Distribution du nombre de résistances antibiotiques des souches hors BLSE en 2020")

write.table(top10_freq_profils_sBLSE2020, sep = ",", quote = FALSE, row.names = F)

write.table(final_antibio_sBLSE2020, sep = ",", quote = FALSE, row.names = T)

#sans BLSE 21------- 

data2021<-read.csv(file="EC-2021.csv",sep=";")


sBLSE2021 <- data2021
unique(sBLSE2021$phenotype) #on choisit les phénotypes dit BLSE, BLSE+IMPER (4) et BLSE+CARBA (2)

sBLSE2021 <- filter(sBLSE2021, type.prelevement == "URINES")

sBLSE2021 <- sBLSE2021[-which(sBLSE2021$phenotype == "BLSE+CARBAPENEMASE (NDM)" | sBLSE2021$phenotype == "BLSE" | sBLSE2021$phenotype == "BLSE+CASE" | sBLSE2021$phenotype == "BLSE+CASE+IMPER" | sBLSE2021$phenotype == "BLSE+IMPER" | sBLSE2021$phenotype == "BLSE+CARBAPENEMASE (OXA-48)" | sBLSE2021$phenotype == "blse" | sBLSE2021$phenotype == "BLSE+CASEHN"),]

unique(sBLSE2021$phenotype)
unique(sBLSE2021$type.prelevement)

sBLSE2021 <- subset(sBLSE2021, select = -c(TCC, CN, AZT, CS, TGC, C3G, FQ)) #supprimer les colonnnes pas pertinentes 

sBLSE2021_res <- sBLSE2021[,11:42]

sBLSE2021_res[sBLSE2021_res == ''] <- NA
sBLSE2021_res[sBLSE2021_res == 'S'] <- FALSE
sBLSE2021_res[sBLSE2021_res == 'I'] <- FALSE
sBLSE2021_res[sBLSE2021_res == 'R'] <- TRUE


sBLSE2021_res_log <- sBLSE2021_res %>% mutate_all(as.logical)  #forme logique pour pouvoir executer les fonctions
str(sBLSE2021_res_log)
antibio_sBLSE2021 <- anti_bio(sBLSE2021_res_log,1:32)
antibio_sBLSE2021 <- as.data.frame(antibio_sBLSE2021)

#calcul du nombre de res antibiotiques de chaque isolat

sBLSE2021_res_num <- sBLSE2021

sBLSE2021_res_num[sBLSE2021_res_num == ""] <- NA
sBLSE2021_res_num[sBLSE2021_res_num == "S"] <- 0
sBLSE2021_res_num[sBLSE2021_res_num == "I"] <- 0
sBLSE2021_res_num[sBLSE2021_res_num == "R"] <- 1


FreqsBLSE2021 <-sBLSE2021_res_num %>% group_by_all() %>% count() %>% ungroup()

uniqfreqsBLSE2021 <- unique(FreqsBLSE2021)

TenTopsBLSE2021<-uniqfreqsBLSE2021 %>%
  arrange(desc(n)) %>%
  slice_head(n=10)

TenTopsBLSE2021 <- TenTopsBLSE2021[order(TenTopsBLSE2021$n,decreasing=TRUE),]


#montrer les atbs ayant le plus de souches testes resistantes 

#changer les colonnes et les lignes 
final_antibio_sBLSE2021 <- as.data.frame(t(antibio_sBLSE2021))
top_mean_sBLSE2021 <- top_n(final_antibio_sBLSE2021, 10, eSupport) 

top_mean_sBLSE2021 %>%  #graphe top 10 resistances
  ggplot(., aes(x=row.names(.), y= eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Les 10 antibiotiques ayant la plus grande prevalence de souches testees hors BLSE resistantes en 2021") +
  guides(x = guide_axis(angle = 30))

#distribution des resistances d'atb par isolats (a modifier)
final_antibio_sBLSE2021 %>% 
  ggplot(., aes(reorder(row.names(.), -eSupport, sum), eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Fréquence des résistances antibiotiques des souches hors BLSE en 2021") +
  guides(x = guide_axis(angle = 60))+
  theme(axis.text = element_text(face = "bold"))
ggsave("freq_res_sBLSE21.png")

#profils de resistances
tot_profils_sBLSE2021 <- res_profile(sBLSE2021_res_log, 1:32)
uniq_profils_sBLSE2021 <-unique(tot_profils_sBLSE2021)


freq_profils_sBLSE2021 <- tot_profils_sBLSE2021 %>% group_by_all() %>% count() %>% ungroup()



# changement de la colonne freq pour avoir un % du nb total isolats 

sum(freq_profils_sBLSE2021$n)
freq_profils_sBLSE2021_perc <- freq_profils_sBLSE2021
freq_profils_sBLSE2021_perc[,3] <- round(freq_profils_sBLSE2021[,3]/sum(freq_profils_sBLSE2021$n),5)*100

top10_freq_profils_sBLSE2021 <- top_n(freq_profils_sBLSE2021_perc, 10, n)
sum(top10_freq_profils_sBLSE2021[,3]) #correspond à % des resistances 
top10_freq_profils_sBLSE2021 %>%  #graphe profil res
  ggplot(., aes(reorder(Profile, -n, sum), n)) + 
  geom_col() +
  labs(x = "Profils de resistance", y = "Frequence(%)", title = "Les 10 profils de resistances des souches hors BLSE les plus frequents en 2021") +
  guides(x = guide_axis(angle = 30))+
  theme(axis.text = element_text(face = "bold"))
ggsave("top10res_sBLSE21.png")

#distribution du nb d'atb presents dans les multiresitances

freq_profils_sBLSE2021_perc %>%
  ggplot(., aes(x=NombreRes, y= n)) +
  geom_col() +
  labs(x = "Nombre de résistances", y = "Prevalence (%)", title = "Distribution du nombre de résistances antibiotiques des souches hors BLSE en 2021")

write.table(top10_freq_profils_sBLSE2021, sep = ",", quote = FALSE, row.names = F)

write.table(final_antibio_sBLSE2021, sep = ",", quote = FALSE, row.names = T)
