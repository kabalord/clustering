#installation des packages avec une fonction pour la gestion des packages

list_packages_names <- c("fpc", "clusteval", "cluster")
for (packname in list_packages_names){
  if(packname %in% rownames(installed.packages()) == FALSE){
    install.packages(packname)
  }
}

library(fpc)
library(clusteval)
library(cluster)

#1.1 Données

n1<-20
n2<-30
n3<-15
n4<-5
set.seed(123) # important : fixer la graine du générateur
echantillon<-data.frame(rbind(matrix(rnorm(n1*2),ncol=2),
                              matrix(rnorm(n2*2,mean = 4),ncol=2),
                              matrix(c(rnorm(n3,mean = 4),rnorm(n3,mean=-4.5)),ncol=2),
                              matrix(c(rnorm(n4,mean = -2),rnorm(n4,mean=-2.5)),ncol=2)))
names(echantillon)<-c("Variable1","Variable2")


#1.2 Partition

#maxence fonction
random_clust <- function (data, x){
  sample(1:x, nrow(data),
         replace=T)
}

res <- random_clust(echantillon,3) # pour 4 classes 

print(res)

echantillon$Classe <- res
print(echantillon)


quad_clust <- function(df){
  
  med1 <- median(df$Variable1)
  
  med2 <- median(df$Variable2)
  
  df[df$Variable1<med1 & df$Variable2<med2,"Classe"] <- 1
  
  df[df$Variable1>med1 & df$Variable2>med2,"Classe"] <- 2
  
  df[df$Variable1>med1 & df$Variable2<med2,"Classe"] <- 3
  
  df[df$Variable1<med1 & df$Variable2>med2,"Classe"] <- 4
  
  return(df)
  
}


quad_clust(echantillon)


#1.3 Visualisation des classes obtenues


plot(x    = echantillon$Variable1 ,
     y    = echantillon$Variable2 ,
     pch  = 21 ,
     bg   = echantillon$Classe ,
     xlab = "Variable 1" ,
     ylab = "Variable 2")


#1.4 Description des classes obtenues


synthese<-do.call(data.frame,aggregate(echantillon, 
                                       by=list(echantillon$Classe),
                                       FUN = function(x){
                                         return(c(Effectif   = length(x),
                                                  Moyenne    = mean(x),
                                                  Mediane    = median(x),
                                                  Minimum    = min(x),
                                                  Maximum    = max(x),
                                                  Ecart_type = sd(x)))
                                       }))
# knitr::kable(synthese)
print(synthese)


#1.5 Évaluation comparaison


list_packages_names <- c("fpc", "clusteval")
for (packname in list_packages_names){
  if(packname %in% rownames(installed.packages()) == FALSE){
    install.packages(packname)
  }
}


#1.5.1 Visualisation avancée


x1 <- quad_clust(echantillon)

clusplot(pam(x1, 4))

plotcluster(x1, x1$Classe)



x2 <- random_clust(echantillon,3)

clusplot(pam(x2, 3))

plotcluster(x2, x2)


#1.6 Analyse des coefficients de silhouette.


si<-silhouette(echantillon$Classe,dist(echantillon[,1:2]))
plot(si)
print(si)

#Faire de même avec les résultats de quad_clust
resqc<-silhouette(x1$Classe,dist(x1[,1:2]))
plot(resqc)
print(resqc)


#1.7 Calcul d’inertie

# Fonction de calcul d'inertie
# df : un data frame contenant le nuage de points
# p : un vecteur de pondération
inertie<-function(df,p=NULL){
  # Si le vecteur de pondération n'est pas fourni
  # tous les poids valebt 1
  if (is.null(p)){
    p <- rep(1,nrow(df))
  }
  # Calcul du centre de gravité : moyenne de chaque colonne, pondérée par p
  g <- ( p %*% as.matrix(df) ) / sum(p)
  # calcul de l'inertie
  iner <- 0
  for (i in seq(nrow(df))){ # pour chaque point
    iner <- iner + sum((df[i,] - g)^2) * p[i] # ajouter à iner la distance au carré de ce point à g, pondérée par le poids du point
  }
  return(iner)
}

# Fonction de calcul de l'inertie intra-classe
# df : toujours le data frame
# cl : le vecteur des labels d'appartenance aux clusters
inertie_intra<-function(df,cl){
  res<-rep(0,length(unique(cl)))
  for (k in unique(cl)){
    res[k] <- inertie(df[which(cl==k),])
  }
  return(sum(res))
}



# Fonction de clacul de l'inertie interclasse
# df : toujours le data frame
# clu : idem
inertie_inter<-function(df,clu){
  M<-matrix(rep(0,ncol(df)*length(unique(clu))),ncol = ncol(df))
  for (k in unique(clu)){
    M[k,]<-unname(colMeans(df[which(clu==k),]))
  }
  return(inertie(data.frame(M),unname(table(clu))))
}

# -------
#  MAIN
# ------


# Chargement des données
data(cars)

# Partitionnement des données (ici par CAH)
hca <- hclust(dist(cars))
clu <- cutree(hca,h=40)

# Affichage
print(hca)
print(clu)

# Affichage des inerties calculées
print(inertie_intra(cars,clu))
print(inertie_inter(cars,clu))
print(inertie(cars))


print(inertie_intra(echantillon,random_clust(echantillon,3)))
print(inertie_inter(echantillon,random_clust(echantillon,3)))
print(inertie2(echantillon))


#1.8 Similarité

#2 Clustering sur des données « réelles » avec les algorithmes des k-means et CAH

#2.1 Installation du package
install.packages("cluster.datasets")

#2.2 Chargement des données
library(cluster.datasets)
data("acidosis.patients")

#2.3 Première description
summary(acidosis.patients)
str(acidosis.patients)
head(acidosis.patients,n = 4)


#2.4 CAH

#Représenter les données à l’aide d’une matrice de scatterplots
pairs(acidosis.patients)

#Centrer et réduire les données
acidosis.patients.scale <- scale(acidosis.patients)

#Calculer la matrice des distances euclidiennes entre individus
acidosis.patients.dist <-dist(acidosis.patients)

#Effectuer une CAH à partir de cette matrice
acidosis.patients.hclust <- hclust(acidosis.patients.dist)
plot(acidosis.patients.hclust)

#Effectuer le partitionnement correspondant. Afficher la matrice de scatterplots en coloriant les points en fonction de leur classe d’appartenance
dendrogram <- as.dendrogram(acidosis.patients.hclust)
plot(dendrogram)
rect.hclust(acidosis.patients.hclust, k=4)
pairs(acidosis.patients.hclust[,4], col = c("red", "cornflowerblue", "purple", "yellow"))


#2.5  k-means
kmeansmet <- kmeans(acidosis.patients.scale,k=4)
plot(kmeans)



#3 Nombre optimal de classes
install.packages(c("factoextra","NbClust"))
library(factoextra)
library(NbClust)
# Elbow method
fviz_nbclust(acidosis.patients, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")
# Silhouette method
fviz_nbclust(acidosis.patients, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")
# Gap statistic
# nboot = 50 to keep the function speedy.
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(acidosis.patients, kmeans, nstart = 25, method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")


library("NbClust")
nb <- NbClust(acidosis.patients, distance = "euclidean", min.nc = 2,
              max.nc = 10, method = "kmeans")
library("factoextra")
fviz_nbclust(nb)


#tests
setwd("/Users/walterroaserrano/Desktop/M22021/apprentissageNONsupervise")
getwd()
df <- read.csv("varied.csv")
# Elbow method
fviz_nbclust(df, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")
# Silhouette method
fviz_nbclust(df, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")
# Gap statistic
# nboot = 50 to keep the function speedy.
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
fviz_nbclust(df, kmeans, nstart = 25, method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")



# 7 Aides, indications et compléments

# 7.1 Analyses descriptives et multivariées

install.packages("ade4")
library(ade4)

varied <- read.csv(file = 'varied.csv',row.names=1, header=FALSE)

#ACP

dudi.pca(varied)



