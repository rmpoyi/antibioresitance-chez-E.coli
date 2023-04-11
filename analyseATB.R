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
TenTopATB2018 <- TenTopATB2018[order(TenTopATB2018$freq,decreasing=TRUE),]
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
FreqATB2019<-atb2019 %>% group_by_all() %>% count()
#Dix jeux d'antibiotiques testés le plus fréquemment
TenTopATB2019<-top_n(FreqATB2019,10,freq)
#Proportion des isolats testés concernés par ces dix jeux d'antibiotiques
sum(TenTopATB2019[,38])/dim(atb2019)[1]
#Dataframe décrivant ces dix jeux, dans l'ordre de fréquence décroissante
TenTopATB2019 <- TenTopATB2019[order(TenTopATB2019$freq,decreasing=TRUE),]
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
FreqATB2020<-atb2020 %>% group_by_all() %>% count()
#Dix jeux d'antibiotiques testés le plus fréquemment
TenTopATB2020<-top_n(FreqATB2020,10,freq)
#Proportion des isolats testés concernés par ces dix jeux d'antibiotiques
sum(TenTopATB2020[,38])/dim(atb2020)[1]
#Dataframe décrivant ces dix jeux, dans l'ordre de fréquence décroissante
TenTopATB2020 <- TenTopATB2020[order(TenTopATB2020$freq,decreasing=TRUE),]
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
FreqATB2021<-atb2021 %>% group_by_all() %>% count()
#Dix jeux d'antibiotiques testés le plus fréquemment
TenTopATB2021<-top_n(FreqATB2021,10,freq)
#Proportion des isolats testés concernés par ces dix jeux d'antibiotiques
sum(TenTopATB2021[,38])/dim(atb2021)[1]
#Dataframe décrivant ces dix jeux, dans l'ordre de fréquence décroissante
TenTopATB2021 <- TenTopATB2021[order(TenTopATB2021$freq,decreasing=TRUE),]
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
trendTest<-rep(0, times = 37)
for (i in 1:37) {
trendTest[i]<-(cor.test(Year,PropATB[,i],method="spearman"))$p.value
}
trendTest
