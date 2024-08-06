library(tidyverse)
library(phytools)
library(geiger)
# matriz filogenetica
data <- read.csv ("dispersalsyndromes.csv", header = T )
## unifico tipo de cerdas
for(i in 1:150){ 
  if (data$tipo_cerda[i] == 2 | data$tipo_cerda[i] == 3){ 
    data$tipo_cerda[i] <- 1 
  }
}

#cambio de numero a string la clasificacion de la variable
data$tipo_cerda <- factor(data$tipo_cerda, levels = c(1,5,4,8,0), labels = c("Setose","Paleaceous","Aristate","Coroniform","Epappose" ))
colnames(data)[7] <- 'pappus_type'
spp <- data$spp
spp <- gsub('_','-', spp) ## cambio el guion bajo por medio
rownames(data) <- spp

# arbol ultrametrico
tree <- read.newick('ultratree.nwk')
tips <- tree$tip.label
## borro las comillas
for (i in 1:length(tips)){
  tips[i] <- substr(tips[i],2,nchar(tips[i])-1)
}
tree$tip.label <- tips

#podo arbol y data
nas <- 'Nastanthus-patagonicus'
data <- data[!rownames(data)%in%nas,] ## saco nasthantus
name.check(tree,data)
# antes tengo que seguir corrigiendo los nombres
right <- c('Calea-candolleana', "Chrysothamnus-nauseosus", "Eutrochium fistulosum",
           "Gladdiopappus-vernonoidea", "Helianthus-niveus","Lactuca-florida", "Pallenis-maritima",
           "Pteronia-camphorata ")
fake <- c('Calea-candollena', "Chrysothanmnus-nauseosus", "Eutrochium-fistulosum", 
          "Gladdiopappus-vernonioidea", "Helianthus-niveus-spp-tephrodes", "Lactuca-floridana", 
          "Pallenis-maritimus", "Pteronia-camphorata")
subdata <- data[fake,]
rownames(subdata) <- right
data <- rbind.data.frame(data,subdata) ## tengo llas filas mal escritas y las bien escritas en el mismo df.
## luego se cortan con el arbol
chk <- name.check(tree, data)
tree2 <- drop.tip(tree, chk$tree_not_data)
data <- data[tree2$tip.label,] ## corto y ordeno como esta en el arbol
name.check(tree2, data)
which(tree2$tip.label == "Lactuca-plumieri")
tree2$tip.label[87] <- "Lactuca-plumieri.1"
name.check(tree2, data)
# calculo variables

data$area<-data$long_cerda_max^2*pi
data$vol<-pi*(data$anch_aquenio/2)^2*(data$long_aquenio/2)
data$PL<-data$vol/(data$area+0.0000001)
data$log_PL<-log(data$vol/(data$area+0.0000001))
data$loglog_PL<-log(data$vol)/log(data$area+0.0000001)
data$sqrt_PL <- sqrt(data$PL)
name.check(tree2, data)
write.csv(data, file = 'final-df.csv')
write.tree(tree2, file = 'final-tree.nwk')
