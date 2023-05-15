#disribution predite des multiresistances et pattern de res selon l'hypothese nulle de l'independance des traits de resistances 
#simule 100 jeu de donnees 

library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)

set.seed(22)

# sur le jeu de base sans TCC, CN, AZT, CS et TGC 

#2018------------- 

data18<-read.csv(file="EC-2018.csv",sep=";")

data18<-data18[data18$type.prelevement=="URINES",]

res18 <- data18[,11:49]

res18 <- subset(res18, select = -c(TCC, CN, AZT, CS, TGC, C3G, FQ))

res18[res18 == ''] <- NA
res18[res18 == 'S'] <- FALSE
res18[res18 == 'I'] <- FALSE
res18[res18 == 'R'] <- TRUE

log18 <- res18 %>% mutate_all(as.logical)

randlog18 <- rep(list(log18[,1:32]),100)
n_obs <- nrow(log18)
n_pasNA <- n_obs - colSums(is.na(log18[,1:32]))
prev <- sapply(log18[,1:32], sum, na.rm = TRUE)/n_pasNA

for (bd in 1:100) {
  for (atb in 1:length(prev)) {
    binom <- rbinom(n= n_pasNA[atb], size = 1, prob = prev[atb])
    rand <- randlog18[[bd]][,atb]
    rand[!is.na(rand)] <- binom 
    randlog18[[bd]][,atb] <- as.logical(rand)
    
  }
  
}

multires18 <- data.frame()

rm(binom, n_obs, prev, rand)

for (bd in 1:100) {
  x <- mdr_profile(randlog18[[bd]], seq(1:ncol(randlog18[[bd]])))
  x <- mutate(x, db= rep(bd, nrow(x)))
  multires18 <- rbind(multires18, x)
  
}

multires18$NombreRes <- as.factor(multires18$NombreRes)


multires18 %>% 
  group_by(db, NombreRes) %>%
  summarise(profiles = n_distinct(Profile)) %>% ungroup() %>%
  complete(db, NombreRes, fill= list(profiles=0)) -> unique_profil

multires18 %>%
  group_by(db, NombreRes) %>% 
  summarise(souches= n()) %>% ungroup() %>%
  complete(db, NombreRes, fill = list(souches=0)) -> compte

table_multires18 <- inner_join(compte, unique_profil, by=c("db","NombreRes"))
#table_multires18 <- table_multires18[,-4]
rm(unique_profil, compte)
names(table_multires18) <- c("Jeu", "NombreRes","NombreSouches","NombreProfils")


mdr18 <- mdr_profile(log18,1:32)

mdr18$NombreRes <- as.factor(mdr18$NombreRes)

mdr18 %>% 
  group_by(NombreRes) %>% 
  summarise(profiles = n_distinct(Profile)) %>% ungroup() %>%
  complete(NombreRes, fill= list(n=0)) -> mdr_profil

mdr18 %>%
  group_by(NombreRes) %>% 
  summarise(souches= n()) %>% ungroup() %>%
  complete(NombreRes, fill = list(souches=0)) -> mdr_compte

table_mdr18 <- inner_join(mdr_compte, mdr_profil, by=c("NombreRes"))
table_mdr18 <- table_mdr18[,-3]
rm(mdr_profil, mdr_compte)
names(table_mdr18) <- c("NombreRes","NombreSouches","NombreProfils")
table_mdr18 <- cbind("Jeu" = rep(101, nrow(table_mdr18)), table_mdr18)


plot_table_multires18 <-data.frame()
plot_table_multires18 <- bind_rows(table_multires18, table_mdr18)

plot_table_multires18 <- group_by(plot_table_multires18, Jeu) %>% mutate(StdNombreProfils = NombreProfils/sum(NombreProfils))

colors <- c(rep("black", 100), "red") #actual data is red; simulated data is black
names(colors) <- seq(1,101,1)
size <- c(rep(0.1,100), 1) #actual data is thicker line than simulated data
names(size) <- seq(1,101,1)
alpha <- c(rep(0.25,100), 1) #simulated data is partially transparent
names(alpha) <- seq(1,101,1)


pourc_souches_mdr_plot <- ggplot(plot_table_multires18, 
                                 aes(x=NombreRes, group = Jeu, color=as.character(Jeu), size=as.character(Jeu), alpha=as.character(Jeu))) +
  geom_line(aes(y=NombreSouches/nrow(log18)))+ #all db have same number of isolates (nrow(SA.db))
  scale_color_manual(values=colors)+
  scale_size_manual(values=size)+
  scale_alpha_manual(values=alpha)+
  scale_y_continuous(labels=scales::percent_format(accuracy=1), expand=c(0,0), limits=c(0,0.55))+
  scale_x_discrete(expand=c(0,0))+
  labs(y="Pourcentage de souches",
       x="Nombre de résistances",
       title = "Distribution du nombre de résistances antibiotiques en 2018")+
  theme_light()+
  theme(legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title=element_text(size=22), axis.text=element_text(size=20))
pourc_souches_mdr_plot

ggsave("test_indp_distrib_2018.png", width = 24, height = 12, units = "in")

freq_profils_mdr_plot <- ggplot(plot_table_multires18, 
                                aes(x=NombreRes, group = Jeu, color=as.character(Jeu), size=as.character(Jeu), alpha=as.character(Jeu))) +
  geom_line(aes(y=StdNombreProfils))+
  scale_color_manual(values=colors)+
  scale_size_manual(values=size)+
  scale_alpha_manual(values=alpha)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand=c(0,0), limits=c(0,0.55))+
  scale_x_discrete(expand=c(0,0))+
  #geom_line(data = table_multires18[table_multires18$Jeu!=1,], aes())
  labs(y="Fréquence\n des profils de resistances distincts",
       x="Nombre de résistances",
       title = "Fréquence des profils de multiresistances distincs en 2018")+
  theme_light()+
  theme(legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title=element_text(size=22), axis.text=element_text(size=20))


freq_profils_mdr_plot

ggsave("test_indp_freq_2018.png", width = 24, height = 12, units = "in")

#2019-----------------

data19<-read.csv(file="EC-2019.csv",sep=";")

data19<-data19[data19$type.prelevement=="URINES",]

res19 <- data19[,11:49]

res19 <- subset(res19, select = -c(TCC, CN, AZT, CS, TGC, C3G, FQ))

res19[res19 == ''] <- NA
res19[res19 == 'S'] <- FALSE
res19[res19 == 'I'] <- FALSE
res19[res19 == 'R'] <- TRUE

log19 <- res19 %>% mutate_all(as.logical)

randlog19 <- rep(list(log19[,1:32]),100)
n_obs <- nrow(log19)
n_pasNA <- n_obs - colSums(is.na(log19[,1:32]))
prev <- sapply(log19[,1:32], sum, na.rm = TRUE)/n_pasNA

for (bd in 1:100) {
  for (atb in 1:length(prev)) {
    binom <- rbinom(n= n_pasNA[atb], size = 1, prob = prev[atb])
    rand <- randlog19[[bd]][,atb]
    rand[!is.na(rand)] <- binom 
    randlog19[[bd]][,atb] <- as.logical(rand)
    
  }
  
}

multires19 <- data.frame()

for (bd in 1:100) {
  x <- mdr_profile(randlog19[[bd]], seq(1:ncol(randlog19[[bd]])))
  x <- mutate(x, db= rep(bd, nrow(x)))
  multires19 <- rbind(multires19, x)
  
}

multires19$NombreRes <- as.factor(multires19$NombreRes)


multires19 %>% 
  group_by(db, NombreRes) %>%
  summarise(profiles = n_distinct(Profile)) %>% ungroup() %>%
  complete(db, NombreRes, fill= list(profiles=0)) -> unique_profil

multires19 %>%
  group_by(db, NombreRes) %>% 
  summarise(souches= n()) %>% ungroup() %>%
  complete(db, NombreRes, fill = list(souches=0)) -> compte

table_multires19 <- inner_join(compte, unique_profil, by=c("db","NombreRes"))
#table_multires19 <- table_multires19[,-4]
rm(unique_profil, compte)
names(table_multires19) <- c("Jeu", "NombreRes","NombreSouches","NombreProfils")

mdr19 <- mdr_profile(log19,1:32)

mdr19$NombreRes <- as.factor(mdr19$NombreRes)

mdr19 %>% 
  group_by(NombreRes) %>% 
  summarise(profiles = n_distinct(Profile)) %>% ungroup() %>%
  complete(NombreRes, fill= list(n=0)) -> mdr_profil

mdr19 %>%
  group_by(NombreRes) %>% 
  summarise(souches= n()) %>% ungroup() %>%
  complete(NombreRes, fill = list(souches=0)) -> mdr_compte

table_mdr19 <- inner_join(mdr_compte, mdr_profil, by=c("NombreRes"))
#table_mdr19 <- table_mdr19[,-3]
rm(mdr_profil, mdr_compte)
names(table_mdr19) <- c("NombreRes","NombreSouches","NombreProfils")
table_mdr19 <- cbind("Jeu" = rep(101, nrow(table_mdr19)), table_mdr19)


plot_table_multires19 <-data.frame()
plot_table_multires19 <- bind_rows(table_multires19, table_mdr19)

plot_table_multires19 <- group_by(plot_table_multires19, Jeu) %>% mutate(StdNombreProfils = NombreProfils/sum(NombreProfils))

colors <- c(rep("black", 100), "red") #actual data is red; simulated data is black
names(colors) <- seq(1,101,1)
size <- c(rep(0.1,100), 1) #actual data is thicker line than simulated data
names(size) <- seq(1,101,1)
alpha <- c(rep(0.25,100), 1) #simulated data is partially transparent
names(alpha) <- seq(1,101,1)


pourc_souches_mdr_plot19 <- ggplot(plot_table_multires19, 
                                 aes(x=NombreRes, group = Jeu, color=as.character(Jeu), size=as.character(Jeu), alpha=as.character(Jeu))) +
  geom_line(aes(y=NombreSouches/nrow(log19)))+ #all db have same number of isolates (nrow(SA.db))
  scale_color_manual(values=colors)+
  scale_size_manual(values=size)+
  scale_alpha_manual(values=alpha)+
  scale_y_continuous(labels=scales::percent_format(accuracy=1), expand=c(0,0), limits=c(0,0.55))+
  scale_x_discrete(expand=c(0,0))+
  labs(y="Pourcentage de souches",
       x="Nombre de résistances",
       title = "Distribution du nombre de résistances antibiotiques en 2019")+
  theme_light()+
  theme(legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title=element_text(size=22), axis.text=element_text(size=20))
pourc_souches_mdr_plot19

ggsave("test_indp_distrib_2019.png", width = 24, height = 12, units = "in")

freq_profils_mdr_plot19 <- ggplot(plot_table_multires19, 
                                aes(x=NombreRes, group = Jeu, color=as.character(Jeu), size=as.character(Jeu), alpha=as.character(Jeu))) +
  geom_line(aes(y=StdNombreProfils))+
  scale_color_manual(values=colors)+
  scale_size_manual(values=size)+
  scale_alpha_manual(values=alpha)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand=c(0,0), limits=c(0,0.55))+
  scale_x_discrete(expand=c(0,0))+
  #geom_line(data = table_multires19[table_multires19$Jeu!=1,], aes())
  labs(y="Fréquence\n des profils de resistances distincts",
       x="Nombre de résistances",
       title = "Fréquence des profils de multiresistances distincts en 2019")+
  theme_light()+
  theme(legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title=element_text(size=22), axis.text=element_text(size=20))


freq_profils_mdr_plot19

ggsave("test_indp_freq_2019.png", width = 24, height = 12, units = "in")


#2020-----------------------

rm(multires18, multires19, data18, data19, log18, log19, mdr18, mdr19, randlog18, randlog19)
data20<-read.csv(file="EC-2020.csv",sep=";")

data20<-data20[data20$type.prelevement=="URINES",]

res20 <- data20[,11:49]

res20 <- subset(res20, select = -c(TCC, CN, AZT, CS, TGC, C3G, FQ))

res20[res20 == ''] <- NA
res20[res20 == 'S'] <- FALSE
res20[res20 == 'I'] <- FALSE
res20[res20 == 'R'] <- TRUE

log20 <- res20 %>% mutate_all(as.logical)

randlog20 <- rep(list(log20[,1:32]),100)
n_obs <- nrow(log20)
n_pasNA <- n_obs - colSums(is.na(log20[,1:32]))
prev <- sapply(log20[,1:32], sum, na.rm = TRUE)/n_pasNA


for (bd in 1:100) {
  for (atb in 1:length(prev)) {
    binom <- rbinom(n= n_pasNA[atb], size = 1, prob = prev[atb])
    rand <- randlog20[[bd]][,atb]
    rand[!is.na(rand)] <- binom 
    randlog20[[bd]][,atb] <- as.logical(rand)
    
  }
  
}

multires20 <- data.frame()

for (bd in 1:100) {
  x <- mdr_profile(randlog20[[bd]], seq(1:ncol(randlog20[[bd]])))
  x <- mutate(x, db= rep(bd, nrow(x)))
  multires20 <- rbind(multires20, x)
  
}

multires20$NombreRes <- as.factor(multires20$NombreRes)


multires20 %>% 
  group_by(db, NombreRes) %>%
  summarise(profiles = n_distinct(Profile)) %>% ungroup() %>%
  complete(db, NombreRes, fill= list(profiles=0)) -> unique_profil

multires20 %>%
  group_by(db, NombreRes) %>% 
  summarise(souches= n()) %>% ungroup() %>%
  complete(db, NombreRes, fill = list(souches=0)) -> compte

table_multires20 <- inner_join(compte, unique_profil, by=c("db","NombreRes"))
#table_multires20 <- table_multires20[,-4]
rm(unique_profil, compte)
names(table_multires20) <- c("Jeu", "NombreRes","NombreSouches","NombreProfils")

mdr20 <- mdr_profile(log20,1:32)

mdr20$NombreRes <- as.factor(mdr20$NombreRes)

mdr20 %>% 
  group_by(NombreRes) %>% 
  summarise(profiles = n_distinct(Profile)) %>% ungroup() %>%
  complete(NombreRes, fill= list(n=0)) -> mdr_profil

mdr20 %>%
  group_by(NombreRes) %>% 
  summarise(souches= n()) %>% ungroup() %>%
  complete(NombreRes, fill = list(souches=0)) -> mdr_compte

table_mdr20 <- inner_join(mdr_compte, mdr_profil, by=c("NombreRes"))
#table_mdr20 <- table_mdr20[,-3]
rm(mdr_profil, mdr_compte)
names(table_mdr20) <- c("NombreRes","NombreSouches","NombreProfils")
table_mdr20 <- cbind("Jeu" = rep(101, nrow(table_mdr20)), table_mdr20)


plot_table_multires20 <-data.frame()
plot_table_multires20 <- bind_rows(table_multires20, table_mdr20)

plot_table_multires20 <- group_by(plot_table_multires20, Jeu) %>% mutate(StdNombreProfils = NombreProfils/sum(NombreProfils))

colors <- c(rep("black", 100), "red") #actual data is red; simulated data is black
names(colors) <- seq(1,101,1)
size <- c(rep(0.1,100), 1) #actual data is thicker line than simulated data
names(size) <- seq(1,101,1)
alpha <- c(rep(0.25,100), 1) #simulated data is partially transparent
names(alpha) <- seq(1,101,1)


pourc_souches_mdr_plot20 <- ggplot(plot_table_multires20, 
                                   aes(x=NombreRes, group = Jeu, color=as.character(Jeu), size=as.character(Jeu), alpha=as.character(Jeu))) +
  geom_line(aes(y=NombreSouches/nrow(log20)))+ #all db have same number of isolates (nrow(SA.db))
  scale_color_manual(values=colors)+
  scale_size_manual(values=size)+
  scale_alpha_manual(values=alpha)+
  scale_y_continuous(labels=scales::percent_format(accuracy=1), expand=c(0,0), limits=c(0,0.55))+
  scale_x_discrete(expand=c(0,0))+
  labs(y="Pourcentage de souches",
       x="Nombre de résistances",
       title = "Distribution du nombre de résistances antibiotiques en 2020")+
  theme_light()+
  theme(legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title=element_text(size=22), axis.text=element_text(size=20))
pourc_souches_mdr_plot20

ggsave("test_indp_distrib_2020.png", width = 24, height = 12, units = "in")

freq_profils_mdr_plot20 <- ggplot(plot_table_multires20, 
                                  aes(x=NombreRes, group = Jeu, color=as.character(Jeu), size=as.character(Jeu), alpha=as.character(Jeu))) +
  geom_line(aes(y=StdNombreProfils))+
  scale_color_manual(values=colors)+
  scale_size_manual(values=size)+
  scale_alpha_manual(values=alpha)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand=c(0,0), limits=c(0,0.55))+
  scale_x_discrete(expand=c(0,0))+
  #geom_line(data = table_multires20[table_multires20$Jeu!=1,], aes())
  labs(y="Fréquence\n des profils de resistances distincts",
       x="Nombre de résistances",
       title = "Fréquence des profils de multiresistances distincts en 2020")+
  theme_light()+
  theme(legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title=element_text(size=22), axis.text=element_text(size=20))


freq_profils_mdr_plot20

ggsave("test_indp_freq_2020.png", width = 24, height = 12, units = "in")

#BLSE18-----------------

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

randBLSE18 <- rep(list(BLSE2018_res_log[,1:30]),100)
n_obs <- nrow(BLSE2018_res_log)
n_pasNA <- n_obs - colSums(is.na(BLSE2018_res_log[,1:30]))
prev <- sapply(BLSE2018_res_log[,1:30], sum, na.rm = TRUE)/n_pasNA

for (bd in 1:100) {
  for (atb in 1:length(prev)) {
    binom <- rbinom(n= n_pasNA[atb], size = 1, prob = prev[atb])
    rand <- randBLSE18[[bd]][,atb]
    rand[!is.na(rand)] <- binom 
    randBLSE18[[bd]][,atb] <- as.logical(rand)
    
  }
  
}

multires_BLSE18 <- data.frame()

for (bd in 1:100) {
  x <- mdr_profile(randBLSE18[[bd]], seq(1:ncol(randBLSE18[[bd]])))
  x <- mutate(x, db= rep(bd, nrow(x)))
  multires_BLSE18 <- rbind(multires_BLSE18, x)
  
}

multires_BLSE18$NombreRes <- as.factor(multires_BLSE18$NombreRes)


multires_BLSE18 %>% 
  group_by(db, NombreRes) %>%
  summarise(profiles = n_distinct(Profile)) %>% ungroup() %>%
  complete(db, NombreRes, fill= list(profiles=0)) -> unique_profil

multires_BLSE18 %>%
  group_by(db, NombreRes) %>% 
  summarise(souches= n()) %>% ungroup() %>%
  complete(db, NombreRes, fill = list(souches=0)) -> compte

table_multires_BLSE18 <- inner_join(compte, unique_profil, by=c("db","NombreRes"))
#table_multires18 <- table_multires18[,-4]
rm(unique_profil, compte)
names(table_multires_BLSE18) <- c("Jeu", "NombreRes","NombreSouches","NombreProfils")

mdr_BLSE18 <- mdr_profile(BLSE2018_res_log,1:30)

mdr_BLSE18$NombreRes <- as.factor(mdr_BLSE18$NombreRes)

mdr_BLSE18 %>% 
  group_by(NombreRes) %>% 
  summarise(profiles = n_distinct(Profile)) %>% ungroup() %>%
  complete(NombreRes, fill= list(n=0)) -> mdr_profil

mdr_BLSE18 %>%
  group_by(NombreRes) %>% 
  summarise(souches= n()) %>% ungroup() %>%
  complete(NombreRes, fill = list(souches=0)) -> mdr_compte

table_mdr_BLSE18 <- inner_join(mdr_compte, mdr_profil, by=c("NombreRes"))
table_mdr_BLSE18 <- table_mdr_BLSE18[,-3]
rm(mdr_profil, mdr_compte)
names(table_mdr_BLSE18) <- c("NombreRes","NombreSouches","NombreProfils")
table_mdr_BLSE18 <- cbind("Jeu" = rep(101, nrow(table_mdr_BLSE18)), table_mdr_BLSE18)


plot_table_multires_BLSE18 <-data.frame()
plot_table_multires_BLSE18 <- bind_rows(table_multires_BLSE18, table_mdr_BLSE18)

plot_table_multires_BLSE18 <- group_by(plot_table_multires_BLSE18, Jeu) %>% mutate(StdNombreProfils = NombreProfils/sum(NombreProfils))

colors <- c(rep("black", 100), "red") #actual data is red; simulated data is black
names(colors) <- seq(1,101,1)
size <- c(rep(0.1,100), 1) #actual data is thicker line than simulated data
names(size) <- seq(1,101,1)
alpha <- c(rep(0.25,100), 1) #simulated data is partially transparent
names(alpha) <- seq(1,101,1)


pourc_souches_mdr_plot <- ggplot(plot_table_multires_BLSE18, 
                                 aes(x=NombreRes, group = Jeu, color=as.character(Jeu), size=as.character(Jeu), alpha=as.character(Jeu))) +
  geom_line(aes(y=NombreSouches/nrow(BLSE2018_res_log)))+ #all db have same number of isolates (nrow(SA.db))
  scale_color_manual(values=colors)+
  scale_size_manual(values=size)+
  scale_alpha_manual(values=alpha)+
  scale_y_continuous(labels=scales::percent_format(accuracy=1), expand=c(0,0), limits=c(0,0.55))+
  scale_x_discrete(expand=c(0,0))+
  labs(y="Pourcentage de souches BLSE",
       x="Nombre de résistances",
       title = "Distribution du nombre de résistances antibiotiques des BLSE en 2018")+
  theme_light()+
  theme(legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title=element_text(size=22), axis.text=element_text(size=20))
pourc_souches_mdr_plot

ggsave("test_indp_distrib_BLSE2018.png", width = 24, height = 12, units = "in")

freq_profils_mdr_plot <- ggplot(plot_table_multires_BLSE18, 
                                aes(x=NombreRes, group = Jeu, color=as.character(Jeu), size=as.character(Jeu), alpha=as.character(Jeu))) +
  geom_line(aes(y=StdNombreProfils))+
  scale_color_manual(values=colors)+
  scale_size_manual(values=size)+
  scale_alpha_manual(values=alpha)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand=c(0,0), limits=c(0,0.55))+
  scale_x_discrete(expand=c(0,0))+
  #geom_line(data = table_multires18[table_multires18$Jeu!=1,], aes())
  labs(y="Fréquence\n des profils de resistances distincts des BLSE",
       x="Nombre de résistances",
       title = "Fréquence des profils de multiresistances distincs de BLSE en 2018")+
  theme_light()+
  theme(legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title=element_text(size=22), axis.text=element_text(size=20))


freq_profils_mdr_plot

ggsave("test_indp_freq_BLSE2018.png", width = 24, height = 12, units = "in")


#sans BLSE18----------

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


sBLSE2018_res_log <- sBLSE2018_res %>% mutate_all(as.logical)

randsBLSE18 <- rep(list(sBLSE2018_res_log[,1:32]),100)
n_obs <- nrow(sBLSE2018_res_log)
n_pasNA <- n_obs - colSums(is.na(sBLSE2018_res_log[,1:32]))
prev <- sapply(sBLSE2018_res_log[,1:32], sum, na.rm = TRUE)/n_pasNA

for (bd in 1:100) {
  for (atb in 1:length(prev)) {
    binom <- rbinom(n= n_pasNA[atb], size = 1, prob = prev[atb])
    rand <- randsBLSE18[[bd]][,atb]
    rand[!is.na(rand)] <- binom 
    randsBLSE18[[bd]][,atb] <- as.logical(rand)
    
  }
  
}

multires_sBLSE18 <- data.frame()

for (bd in 1:100) {
  x <- mdr_profile(randsBLSE18[[bd]], seq(1:ncol(randsBLSE18[[bd]])))
  x <- mutate(x, db= rep(bd, nrow(x)))
  multires_sBLSE18 <- rbind(multires_sBLSE18, x)
  
}

multires_sBLSE18$NombreRes <- as.factor(multires_sBLSE18$NombreRes)


multires_sBLSE18 %>% 
  group_by(db, NombreRes) %>%
  summarise(profiles = n_distinct(Profile)) %>% ungroup() %>%
  complete(db, NombreRes, fill= list(profiles=0)) -> unique_profil

multires_sBLSE18 %>%
  group_by(db, NombreRes) %>% 
  summarise(souches= n()) %>% ungroup() %>%
  complete(db, NombreRes, fill = list(souches=0)) -> compte

table_multires_sBLSE18 <- inner_join(compte, unique_profil, by=c("db","NombreRes"))
#table_multires18 <- table_multires18[,-4]
rm(unique_profil, compte)
names(table_multires_sBLSE18) <- c("Jeu", "NombreRes","NombreSouches","NombreProfils")

mdr_sBLSE18 <- mdr_profile(sBLSE2018_res_log,1:32)

mdr_sBLSE18$NombreRes <- as.factor(mdr_sBLSE18$NombreRes)

mdr_sBLSE18 %>% 
  group_by(NombreRes) %>% 
  summarise(profiles = n_distinct(Profile)) %>% ungroup() %>%
  complete(NombreRes, fill= list(n=0)) -> mdr_profil

mdr_sBLSE18 %>%
  group_by(NombreRes) %>% 
  summarise(souches= n()) %>% ungroup() %>%
  complete(NombreRes, fill = list(souches=0)) -> mdr_compte

table_mdr_sBLSE18 <- inner_join(mdr_compte, mdr_profil, by=c("NombreRes"))
table_mdr_sBLSE18 <- table_mdr_sBLSE18[,-3]
rm(mdr_profil, mdr_compte)
names(table_mdr_sBLSE18) <- c("NombreRes","NombreSouches","NombreProfils")
table_mdr_sBLSE18 <- cbind("Jeu" = rep(101, nrow(table_mdr_sBLSE18)), table_mdr_sBLSE18)


plot_table_multires_sBLSE18 <-data.frame()
plot_table_multires_sBLSE18 <- bind_rows(table_multires_sBLSE18, table_mdr_sBLSE18)

plot_table_multires_sBLSE18 <- group_by(plot_table_multires_sBLSE18, Jeu) %>% mutate(StdNombreProfils = NombreProfils/sum(NombreProfils))

colors <- c(rep("black", 100), "red") #actual data is red; simulated data is black
names(colors) <- seq(1,101,1)
size <- c(rep(0.1,100), 1) #actual data is thicker line than simulated data
names(size) <- seq(1,101,1)
alpha <- c(rep(0.25,100), 1) #simulated data is partially transparent
names(alpha) <- seq(1,101,1)


pourc_souches_mdr_plot <- ggplot(plot_table_multires_sBLSE18, 
                                 aes(x=NombreRes, group = Jeu, color=as.character(Jeu), size=as.character(Jeu), alpha=as.character(Jeu))) +
  geom_line(aes(y=NombreSouches/nrow(sBLSE2018_res_log)))+ #all db have same number of isolates (nrow(SA.db))
  scale_color_manual(values=colors)+
  scale_size_manual(values=size)+
  scale_alpha_manual(values=alpha)+
  scale_y_continuous(labels=scales::percent_format(accuracy=1), expand=c(0,0), limits=c(0,0.55))+
  scale_x_discrete(expand=c(0,0))+
  labs(y="Pourcentage de souches hors BLSE",
       x="Nombre de résistances",
       title = "Distribution du nombre de résistances antibiotiques hors BLSE en 2018")+
  theme_light()+
  theme(legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title=element_text(size=22), axis.text=element_text(size=20))
pourc_souches_mdr_plot

ggsave("test_indp_distrib_sBLSE2018.png", width = 24, height = 12, units = "in")

freq_profils_mdr_plot <- ggplot(plot_table_multires_sBLSE18, 
                                aes(x=NombreRes, group = Jeu, color=as.character(Jeu), size=as.character(Jeu), alpha=as.character(Jeu))) +
  geom_line(aes(y=StdNombreProfils))+
  scale_color_manual(values=colors)+
  scale_size_manual(values=size)+
  scale_alpha_manual(values=alpha)+
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand=c(0,0), limits=c(0,0.55))+
  scale_x_discrete(expand=c(0,0))+
  #geom_line(data = table_multires18[table_multires18$Jeu!=1,], aes())
  labs(y="Fréquence\n des profils de resistances distincts hors BLSE",
       x="Nombre de résistances",
       title = "Fréquence des profils de multiresistances distincs hors BLSE en 2018")+
  theme_light()+
  theme(legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title=element_text(size=22), axis.text=element_text(size=20))


freq_profils_mdr_plot

ggsave("test_indp_freq_sBLSE2018.png", width = 24, height = 12, units = "in")


#essai (pas bon)------------------
uniq_multires18 <- unique(multires18)

multires18 <- multires18 %>% group_by_all() %>% count() %>% ungroup()

sum18 <- sum(nrow(log18))

multires18[,4] <- round(multires18[,4]/sum18,5)*100

multires18 <- group_by(multires18, db)

#set up: 

ggplot(multires18, aes(x=NombreRes, group = db, color=as.character(db), size=as.character(db), alpha=as.character(db)))+
  geom_line(aes(y = n))
   
#essai bloucle for 
for (i in 1:nrow(multires18[,3])) {
  ggplot(multires18, aes(multires18[,2], multires18[,4]))+
  geom_line()+
  labs(x = "Nombre de résistances", y = "Prevalence (%)", title = "Distribution du nombre de résistances antibiotiques des souches testées en 2018 pour chaque jeu générée")
i <- i+1
}
    geom_line() +
    labs(x = "Nombre de résistances", y = "Prevalence (%)", title = "Distribution du nombre de résistances antibiotiques des souches testées en 2018")
  
  
plot_table_multires18 %>%
  select(Jeu, NombreRes, NombreSouches) %>%
  distinct() %>%
  ggplot(., 
         aes(x=as.numeric(NombreRes), group = Jeu))+ #, size=as.character(Jeu), alpha=as.character(Jeu))) +
  geom_line(aes(y=NombreSouches), col = "grey80") 

plot_table_multires18 %>%
  filter(Jeu != 101, !is.na(StdNombreProfils)) %>%
  mutate(NombreRes = as.character(NombreRes)) %>%
  group_by(Jeu, NombreRes) %>%
  summarise(freq = sum(NombreProfils)/unique(NombreSouches), .groups = "drop") %>%
  ggplot(., 
         aes(x=as.numeric(NombreRes), group = Jeu))+ #, size=as.character(Jeu), alpha=as.character(Jeu))) +
  geom_line(aes(y=freq), col = "grey80") #+ 
#geom_line(data = Jeu = 101 , aes())

#faire une boucle for et un graphe en lignes du nombre d'atbs resitants pour chaque data frame donné