library(magrittr)
library(plyr)
library(dplyr)

#Lecture du fichier de données pour 2018
data2018<-read.csv(file="EC-2018.csv",sep=";")
#Filtre pour ne garder que les prélèvements urinaires et enlever les colonnes ajoutées à la main par PRIMO
data2018<-data2018[data2018$type.prelevement=="URINES",]
data2018<-subset(data2018, select = -c(C3G,FQ) )
#Dataframe ne comptant que les données d'antibiogramme, réordonnée dans l'ordre voulu
atb2018<-data2018[,11:47]
atb2018<-atb2018[,c("AKN","GEN","TOB","FOS","TEM","TIC","AMP","AMX","MEC","PIP","AMC","AMC_urine","AMC03","TCC","TZP","ERT","IMP","MEM","CN","CXM","FOX","CAZ","CFM","CRO","CTX","FEP","FUR","CIP","LEV","NOR","OFL","NAL","AZT","CS","TGC","TSU","TMP")]
#Transformation pour n'avoir qu'une indicatrice de test d'antibiotique
atb2018[atb2018 == ''] <- 0
atb2018[atb2018 == 'S'] <- 1
atb2018[atb2018 == 'I'] <- 1
atb2018[atb2018 == 'R'] <- 1
#Nombre de jeux d'antibiotiques uniques testés sur les isolats
SeqATB2018<-unique(atb2018)
dim(SeqATB2018)[1]
#Ajout d'une colonne indiquant la fréquence de chaque jeu d'atb testé
FreqATB2018<-atb2018 %>% group_by_all() %>% count() %>% ungroup()
#Dix jeux d'antibiotiques testés le plus fréquemment
TenTopATB2018<-FreqATB2018 %>%
  arrange(desc(n)) %>%
  slice_head(n=10)

#Proportion des isolats testés concernés par ces dix jeux d'antibiotiques
sum(TenTopATB2018[,38])/dim(atb2018)[1]
#Dataframe décrivant ces dix jeux, dans l'ordre de fréquence décroissante
TenTopATB2018 <- TenTopATB2018[order(TenTopATB2018$n,decreasing=TRUE),]
#Génération d'une table au format txt à partir de cette dataframe
write.table(TenTopATB2018,sep = ",", quote = FALSE, row.names = F)
#Possible de lire cette table dans Word ou Excel
#Calcul de la fréquence de test de chaque atb
atb2018N <- atb2018 %>% mutate_all(as.numeric)
PropATB2018<-colSums(atb2018N)/dim(atb2018)[1]

write.table(PropATB2018,sep = ",", quote = FALSE, row.names = F)


#Lecture du fichier de données pour 2019
data2019<-read.csv(file="EC-2019.csv",sep=";")
#Filtre pour ne garder que les prélèvements urinaires et enlever les colonnes ajoutées à la main par PRIMO
data2019<-data2019[data2019$type.prelevement=="URINES",]
data2019<-subset(data2019, select = -c(C3G,FQ) )
#Dataframe ne comptant que les données d'antibiogramme, réordonnée dans l'ordre voulu
atb2019<-data2019[,11:47]
atb2019<-atb2019[,c("AKN","GEN","TOB","FOS","TEM","TIC","AMP","AMX","MEC","PIP","AMC","AMC_urine","AMC03","TCC","TZP","ERT","IMP","MEM","CN","CXM","FOX","CAZ","CFM","CRO","CTX","FEP","FUR","CIP","LEV","NOR","OFL","NAL","AZT","CS","TGC","TSU","TMP")]
#Transformation pour n'avoir qu'une indicatrice de test d'antibiotique
atb2019[atb2019 == ''] <- 0
atb2019[atb2019 == 'S'] <- 1
atb2019[atb2019 == 'I'] <- 1
atb2019[atb2019 == 'R'] <- 1
#Nombre de jeux d'antibiotiques uniques testés sur les isolats
SeqATB2019<-unique(atb2019)
dim(SeqATB2019)[1]
#Ajout d'une colonne indiquant la fréquence de chaque jeu d'atb testé
FreqATB2019<-atb2019 %>% group_by_all() %>% count() %>% ungroup()
#Dix jeux d'antibiotiques testés le plus fréquemment
TenTopATB2019<-top_n(FreqATB2019,10,n)
#Proportion des isolats testés concernés par ces dix jeux d'antibiotiques
sum(TenTopATB2019[,38])/dim(atb2019)[1]
#Dataframe décrivant ces dix jeux, dans l'ordre de fréquence décroissante
TenTopATB2019 <- TenTopATB2019[order(TenTopATB2019$n,decreasing=TRUE),]
#Génération d'une table au format txt à partir de cette dataframe
write.table(TenTopATB2019,sep = ",", quote = FALSE, row.names = F)
#Possible de lire cette table dans Word ou Excel
#Calcul de la fréquence de test de chaque atb
atb2019N <- atb2019 %>% mutate_all(as.numeric)
PropATB2019<-colSums(atb2019N)/dim(atb2019)[1]
write.table(PropATB2019,sep = ",", quote = FALSE, row.names = F)

#Lecture du fichier de données pour 2020
data2020<-read.csv(file="EC-2020.csv",sep=";")
#Filtre pour ne garder que les prélèvements urinaires et enlever les colonnes ajoutées à la main par PRIMO
data2020<-data2020[data2020$type.prelevement=="URINES",]
data2020<-subset(data2020, select = -c(C3G,FQ) )
#Dataframe ne comptant que les données d'antibiogramme, réordonnée dans l'ordre voulu
atb2020<-data2020[,11:47]
atb2020<-atb2020[,c("AKN","GEN","TOB","FOS","TEM","TIC","AMP","AMX","MEC","PIP","AMC","AMC_urine","AMC03","TCC","TZP","ERT","IMP","MEM","CN","CXM","FOX","CAZ","CFM","CRO","CTX","FEP","FUR","CIP","LEV","NOR","OFL","NAL","AZT","CS","TGC","TSU","TMP")]
#Transformation pour n'avoir qu'une indicatrice de test d'antibiotique
atb2020[atb2020 == ''] <- 0
atb2020[atb2020 == 'S'] <- 1
atb2020[atb2020 == 'I'] <- 1
atb2020[atb2020 == 'R'] <- 1
#Nombre de jeux d'antibiotiques uniques testés sur les isolats
SeqATB2020<-unique(atb2020)
dim(SeqATB2020)[1]
#Ajout d'une colonne indiquant la fréquence de chaque jeu d'atb testé
FreqATB2020<-atb2020 %>% group_by_all() %>% count() %>% ungroup()
#Dix jeux d'antibiotiques testés le plus fréquemment
TenTopATB2020<-top_n(FreqATB2020,10,n)
#Proportion des isolats testés concernés par ces dix jeux d'antibiotiques
sum(TenTopATB2020[,38])/dim(atb2020)[1]
#Dataframe décrivant ces dix jeux, dans l'ordre de fréquence décroissante
TenTopATB2020 <- TenTopATB2020[order(TenTopATB2020$n,decreasing=TRUE),]
#Génération d'une table au format txt à partir de cette dataframe
write.table(TenTopATB2020,sep = ",", quote = FALSE, row.names = F)
#Possible de lire cette table dans Word ou Excel
#Calcul de la fréquence de test de chaque atb
atb2020N <- atb2020 %>% mutate_all(as.numeric)
PropATB2020<-colSums(atb2020N)/dim(atb2020)[1]
write.table(PropATB2020,sep = ",", quote = FALSE, row.names = F)

#Lecture du fichier de données pour 2021
data2021<-read.csv(file="EC-2021.csv",sep=";")
#Filtre pour ne garder que les prélèvements urinaires et enlever les colonnes ajoutées à la main par PRIMO
data2021<-data2021[data2021$type.prelevement=="URINES",]
data2021<-subset(data2021, select = -c(C3G,FQ) )
#Dataframe ne comptant que les données d'antibiogramme, réordonnée dans l'ordre voulu
atb2021<-data2021[,11:47]
atb2021<-atb2021[,c("AKN","GEN","TOB","FOS","TEM","TIC","AMP","AMX","MEC","PIP","AMC","AMC_urine","AMC03","TCC","TZP","ERT","IMP","MEM","CN","CXM","FOX","CAZ","CFM","CRO","CTX","FEP","FUR","CIP","LEV","NOR","OFL","NAL","AZT","CS","TGC","TSU","TMP")]
#Transformation pour n'avoir qu'une indicatrice de test d'antibiotique
atb2021[atb2021 == ''] <- 0
atb2021[atb2021 == 'S'] <- 1
atb2021[atb2021 == 'I'] <- 1
atb2021[atb2021 == 'R'] <- 1
#Nombre de jeux d'antibiotiques uniques testés sur les isolats
SeqATB2021<-unique(atb2021)
dim(SeqATB2021)[1]
#Ajout d'une colonne indiquant la fréquence de chaque jeu d'atb testé
FreqATB2021<-atb2021 %>% group_by_all() %>% count() %>% ungroup()
#Dix jeux d'antibiotiques testés le plus fréquemment
TenTopATB2021<-top_n(FreqATB2021,10,n)
#Proportion des isolats testés concernés par ces dix jeux d'antibiotiques
sum(TenTopATB2021[,38])/dim(atb2021)[1]
#Dataframe décrivant ces dix jeux, dans l'ordre de fréquence décroissante
TenTopATB2021 <- TenTopATB2021[order(TenTopATB2021$n,decreasing=TRUE),]
#Génération d'une table au format txt à partir de cette dataframe
write.table(TenTopATB2021,sep = ",", quote = FALSE, row.names = F)
#Possible de lire cette table dans Word ou Excel
#Calcul de la fréquence de test de chaque atb
atb2021N <- atb2021 %>% mutate_all(as.numeric)
PropATB2021<-colSums(atb2021N)/dim(atb2021)[1]
write.table(PropATB2021,sep = ",", quote = FALSE, row.names = F)

#Evolution temporelle des antibiotiques testés
PropATB<-rbind(PropATB2018,PropATB2019,PropATB2020,PropATB2021)
Year<-c(2018,2019,2020,2021)
trendTest<-rep(0, times = 37) #questions sur ça
for (i in 1:37) {
  trendTest[i]<-(cor.test(Year,PropATB[,i],method="spearman"))$p.value
}
trendTest

#analyses resistances 

atb2018_res <- data2018[,11:47]

atb2018_res[atb2018_res == ''] <- NA
atb2018_res[atb2018_res == 'S'] <- FALSE
atb2018_res[atb2018_res == 'I'] <- FALSE
atb2018_res[atb2018_res == 'R'] <- TRUE


atb2018_res_log <- atb2018_res %>% mutate_all(as.logical)  #forme logique pour pouvoir executer les fonctions
str(atb2018_res_log)
antibio_res <- anti_bio(atb2018_res_log,1:37)
antibio_res <- as.data.frame(antibio_res)

#calcul du nombre de res antibiotiques de chaque isolat

atb2018_res_num <- data2018[,11:47]
atb2018_res_num[atb2018_res_num == ""] <- NA
atb2018_res_num[atb2018_res_num == "S"] <- 0
atb2018_res_num[atb2018_res_num == "I"] <- 0
atb2018_res_num[atb2018_res_num == "R"] <- 1

FreqATBres2018<-atb2018_res_num %>% group_by_all() %>% count() %>% ungroup()
TenTopATB2018res<-FreqATBres2018 %>%
  arrange(desc(n)) %>%
  slice_head(n=10)

TenTopATB2018res <- TenTopATB2018res[order(TenTopATB2018res$n,decreasing=TRUE),]

TenTopATB2018res[TenTopATB2018res == "TRUE"] <- 1
TenTopATB2018res[TenTopATB2018res == "FALSE"] <- 0


#montrer les atbs ayant le plus de souches testes resistantes 

#changer les colonnes et les lignes 
final_antibio_res_2018 <- as.data.frame(t(antibio_res))
top_mean_res_2018 <- top_n(final_antibio_res_2018, 10, eSupport) 

top_mean_res_2018 %>%  #graphe top 10 resistances
  ggplot(., aes(x=row.names(.), y= eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Les 10 antibiotiques ayant la plus grande prevalence de souches testees resistantes en 2018") +
  guides(x = guide_axis(angle = 30))

#distribution des resistances d'atb par isolats (a modifier)
final_antibio_res_2018 %>% 
  ggplot(., aes(x=row.names(.), y= eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Distribution de la prevalence des antibiotiques ayant des souches testees resistantes en 2018") +
  guides(x = guide_axis(angle = 60))

#profils de resistances
tot_profils_res <- res_profile(atb2018_res_log, 1:37)
uniq_profils_res <-unique(tot_profils_res)

freq_profils_res <- tot_profils_res %>% group_by_all() %>% count() %>% ungroup()

# changement de la colonne freq pour avoir un % du nb total isolats 

sum(freq_profils_res$n)
freq_profils_res_perc <- freq_profils_res
freq_profils_res_perc[,3] <- round(freq_profils_res[,3]/sum(freq_profils_res$n),5)*100

top10_freq_profils_res <- top_n(freq_profils_res_perc, 10, n)
sum(top10_freq_profils_res[,3]) #correspond à 70% des resistances 
top10_freq_profils_res %>%  #graphe profil res
  ggplot(., aes(x=Profile, y= n)) + 
  geom_col() +
  labs(x = "Profils de resistance", y = "Frequence(%)", title = "Les 10 profils de resistances les plus frequents en 2018") +
  guides(x = guide_axis(angle = 30))

#distribution du nb d'atb presents dans les multiresitances

freq_profils_res_perc %>%
  ggplot(., aes(x=NombreRes, y= n)) +
  geom_col() +
  labs(x = "Nombre d'antibiotiques resistants", y = "Prevalence (%)", title = "Distribution du nombre d'antibiotiques present dans les multiresistances en 2018") 

# faire la même chose pour 2019, 2020, 2021 
#2019--------------------------- 
atb2019_res <- data2019[,11:47]

atb2019_res[atb2019_res == ''] <- NA
atb2019_res[atb2019_res == 'S'] <- FALSE
atb2019_res[atb2019_res == 'I'] <- FALSE
atb2019_res[atb2019_res == 'R'] <- TRUE


atb2019_res_log <- atb2019_res %>% mutate_all(as.logical)  #forme logique pour pouvoir executer les fonctions
str(atb2019_res_log)
antibio_res_2019 <- anti_bio(atb2019_res_log,1:37)
antibio_res_2019 <- as.data.frame(antibio_res_2019)

#calcul du nombre de res antibiotiques de chaque isolat

FreqATBres2019<-atb2019_res %>% group_by_all() %>% count() %>% ungroup()

#montrer les atbs ayant le plus de souches testes resistantes 

#changer les colonnes et les lignes 
final_antibio_res_2019 <- as.data.frame(t(antibio_res_2019))
top_mean_res_2019 <- top_n(final_antibio_res_2019, 10, eSupport) 

top_mean_res_2019 %>%  #graphe top10 res
  ggplot(., aes(x=row.names(.), y= eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Les 10 antibiotiques ayant la plus grande prevalence de souches testees resistantes en 2019") +
  guides(x = guide_axis(angle = 30))

#distribution des resistances d'atb par isolats (a modifier)
final_antibio_res_2019 %>% 
  ggplot(., aes(x=row.names(.), y= eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Distribution de la prevalence des antibiotiques ayant des souches testees resistantes en 2019") +
  guides(x = guide_axis(angle = 60))

#profils de resistances
tot_profils_res_2019 <- res_profile(atb2019_res_log, 1:37)
uniq_profils_res_2019 <-unique(tot_profils_res_2019)

freq_profils_res_2019 <- tot_profils_res_2019 %>% group_by_all() %>% count() %>% ungroup()

sum(freq_profils_res_2019$n)
freq_profils_res_perc_2019 <- freq_profils_res_2019
freq_profils_res_perc_2019[,3] <- round(freq_profils_res_2019[,3]/sum(freq_profils_res_2019$n),5)*100

top10_freq_profils_res_perc_2019 <- top_n(freq_profils_res_perc_2019, 10, n)
sum(top10_freq_profils_res_perc_2019[,3]) #correspond à ~62% des resistances 
top10_freq_profils_res_perc_2019 %>%  #graphe profils de res les + frequents 
  ggplot(., aes(x=Profile, y= n)) + 
  geom_col() + 
  labs(x = "Profils de resistance", y = "Prevalence (%)", title = "Les 10 profils de resistances les plus frequents en 2019") +
  guides(x = guide_axis(angle = 30))

#distribution du nb d'atb presents dans les multiresitances

freq_profils_res_perc_2019 %>%
  ggplot(., aes(x=NombreRes, y= n)) +
  geom_col() +
  labs(x = "Nombre d'antibiotiques resistants", y = "Prevalence (%)", title = "Distribution du nombre d'antibiotiques present dans les multiresistances en 2019") 

#2020-------------------------------------------------------

atb2020_res <- data2020[,11:47]

atb2020_res[atb2020_res == ''] <- NA
atb2020_res[atb2020_res == 'S'] <- FALSE
atb2020_res[atb2020_res == 'I'] <- FALSE
atb2020_res[atb2020_res == 'R'] <- TRUE


atb2020_res_log <- atb2020_res %>% mutate_all(as.logical)  #forme logique pour pouvoir executer les fonctions
str(atb2020_res_log)
antibio_res_2020 <- anti_bio(atb2020_res_log,1:37)
antibio_res_2020 <- as.data.frame(antibio_res_2020)

#montrer les atbs ayant le plus de souches testes resistantes 

#changer les colonnes et les lignes 
final_antibio_res_2020 <- as.data.frame(t(antibio_res_2020))
top_mean_res_2020 <- top_n(final_antibio_res_2020, 10, eSupport) 

top_mean_res_2020 %>%  #graphe top 10 res
  ggplot(., aes(x=row.names(.), y= eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Les 10 antibiotiques ayant la plus grande prevalence de souches testees resistantes en 2020") +
  guides(x = guide_axis(angle = 30))

#distribution des resistances d'atb par isolats (a modifier)
final_antibio_res_2020 %>% 
  ggplot(., aes(x=row.names(.), y= eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Distribution de la prevalence des antibiotiques ayant des souches testees resistantes en 2020") +
  guides(x = guide_axis(angle = 60))


tot_profils_res_2020 <- res_profile(atb2020_res_log, 1:37)
uniq_profils_res_2020 <-unique(tot_profils_res_2020)

freq_profils_res_2020 <- tot_profils_res_2020 %>% group_by_all() %>% count() %>% ungroup()

freq_profils_res_perc_2020 <- freq_profils_res_2020
freq_profils_res_perc_2020[,3] <- round(freq_profils_res_2020[,3]/sum(freq_profils_res_2020$n),5)*100

top10_freq_profils_res_perc_2020 <- top_n(freq_profils_res_perc_2020, 10, n)
sum(top10_freq_profils_res_perc_2020[,3]) #correspond à 62% des resistances 
top10_freq_profils_res_perc_2020 %>%  #graphe top10 profils res
  ggplot(., aes(x=Profile, y= n)) + 
  geom_col() + 
  labs(x = "Profils de resistance", y = "Prevalence (%)", title = "Les 10 profils de resistances les plus frequents en 2020") +
  guides(x = guide_axis(angle = 30))

#distribution du nb d'atb presents dans les multiresitances

freq_profils_res_perc_2020 %>%
  ggplot(., aes(x=NombreRes, y= n)) +
  geom_col() +
  labs(x = "Nombre d'antibiotiques resistants", y = "Prevalence (%)", title = "Distribution du nombre d'antibiotiques present dans les multiresistances en 2020")

#2021---------------------------

atb2021_res <- data2021[,11:47]

atb2021_res[atb2021_res == ''] <- NA
atb2021_res[atb2021_res == 'S'] <- FALSE
atb2021_res[atb2021_res == 'I'] <- FALSE
atb2021_res[atb2021_res == 'R'] <- TRUE


atb2021_res_log <- atb2021_res %>% mutate_all(as.logical)  #forme logique pour pouvoir executer les fonctions
str(atb2021_res_log)
antibio_res_2021 <- anti_bio(atb2021_res_log,1:37)
antibio_res_2021 <- as.data.frame(antibio_res_2021)

#montrer les atbs ayant le plus de souches testes resistantes 

#changer les colonnes et les lignes 
final_antibio_res_2021 <- as.data.frame(t(antibio_res_2021))
top_mean_res_2021 <- top_n(final_antibio_res_2021, 10, eSupport) 

top_mean_res_2021 %>%  #top 10 res
  ggplot(., aes(x=row.names(.), y= eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Les 10 antibiotiques ayant la plus grande prevalence de souches testees resistantes en 2021") +
  guides(x = guide_axis(angle = 30))

#distribution des resistances d'atb par isolats (a modifier)
final_antibio_res_2021 %>% 
  ggplot(., aes(x=row.names(.), y= eSupport)) + 
  geom_col() +
  labs(x = "Code de l'antibiotique", y = "Prevalence (%)", title = "Distribution de la prevalence des antibiotiques ayant des souches testees resistantes en 2021") +
  guides(x = guide_axis(angle = 60))

tot_profils_res_2021 <- res_profile(atb2021_res_log, 1:37)
uniq_profils_res_2021 <-unique(tot_profils_res_2021)

freq_profils_res_2021 <- tot_profils_res_2021 %>% group_by_all() %>% count() %>% ungroup()

freq_profils_res_perc_2021 <- freq_profils_res_2021
freq_profils_res_perc_2021[,3] <- round(freq_profils_res_2021[,3]/sum(freq_profils_res_2021$n),5)*100

top10_freq_profils_res_perc_2021 <- top_n(freq_profils_res_perc_2021, 10, n)
sum(top10_freq_profils_res_perc_2021[,3]) #correspond à 62% des resistances 
top10_freq_profils_res_perc_2021 %>%  #top 10 profils res
  ggplot(., aes(x=Profile, y= n)) + 
  geom_col()  + 
  labs(x = "Profils de resistance", y = "Prevalence (%)", title = "Les 10 profils de resistances les plus frequents en 2021") +
  guides(x = guide_axis(angle = 30))

#distribution du nb d'atb presents dans les multiresitances

freq_profils_res_perc_2021 %>%
  ggplot(., aes(x=NombreRes, y= n)) +
  geom_col() +
  labs(x = "Nombre d'antibiotiques resistants", y = "Prevalence (%)", title = "Distribution du nombre d'antibiotiques present dans les multiresistances en 2021")

#ecrire tout les data frames sur tableau excel----------------

write.table(TenTopATB2018,sep = ",", quote = FALSE, row.names = F)

write.table(TenTopATB2018res, sep = ",", quote = FALSE, row.names = F)

atb2018N <- atb2018 %>% mutate_all(as.numeric)
PropATB2018<-colSums(atb2018N)/dim(atb2018)[1]

write.table(PropATB2018,sep = ",", quote = FALSE, row.names = F)

#refaire les memes analyses en supprimant les ATB avec le moins de 10% de souches testees 

atb_moins_10 <- c("TCC", "CN", "AZT", "CS", "TGC")

#refaire les memes analyses pour les BLSE et une sans les BLSE

BLSE <- c("AMC","AMC_urine", "AMC03", "TCC", "TZP")

