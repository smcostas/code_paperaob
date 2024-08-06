library(tidyverse)
library(phytools)
library(geiger)
library(viridis)
# Data, tree and creating some new features ######
funcional <- read.csv("funcional.csv", header = T)
rownames(funcional) <- funcional$spp
funcional <- filter(funcional, S_long_cerda_max > 0)
summary(funcional)
funcional$S_vol_mm3 <- pi*(funcional$S_ancho_aq/2)^2*(funcional$S_long_aq/2)
funcional$S_area_mm2 <- funcional$S_long_cerda_max^2*pi
data <- read.csv('final-df.csv', header = T, row.names = 1, stringsAsFactors = T)
tree2 <- read.newick('final-tree.nwk')
name.check(tree2, data)
## to use later ####
data$pappus_type <- factor(data$pappus_type, levels = c('Setose', 'Paleaceous', 'Aristate', 'Coroniform', 'Epappose'))


# Phylomorphospace Fig. 4 #######
head(data)
traits <- data[,c('vol','area')]
traits$area <- traits$area+0.0000001
log_traits <- log(traits)
colnames(log_traits) <- c('log(vol mm3)', 'log(area mm2)')
colr_pappus <- setNames(magma(5, end = 0.94, direction = -1) , levels(as.factor(data$pappus_type)))
pappus <- setNames(data$pappus_type, rownames(data))

## plotting ########
#pdf(file = 'phylomorphospace_paper.pdf')
phylomorphospace(tree2, log_traits, label = "off", node.size = c(0,1))
tiplabels(pie = to.matrix(pappus[tree2$tip.label],levels(pappus)),piecol=colr_pappus,cex=0.35, offset= 1)
add.simmap.legend(colors=colr_pappus,prompt=FALSE,x= 4,
                  y= -5,fsize=0.5, shape = "circle")
points(log(funcional$S_vol_mm3),log(funcional$S_area_mm2+0.0000001), col = "#35B77980", pch=18, cex = 1.5)
setosa = filter(data, pappus_type == 'Setose') %>% select(c(pappus_type,vol,area))
abline(lm(log(area) ~ log(vol), data =   setosa), col = "#FDE2A2FF", lty = 'dashed',  lwd=2)
paleacea = filter(data, pappus_type == 'Paleaceous') %>% select(c(pappus_type,vol,area))
abline(lm(log(area) ~ log(vol), data =   paleacea), col = "#F7725CFF", lty = 'dashed',  lwd=2)
Aristada = filter(data, pappus_type == 'Aristate') %>% select(c(pappus_type,vol,area))
abline(lm(log(area) ~ log(vol), data =   Aristada), col = "#AA337DFF", lty = 'dashed',  lwd=2)
Coroniforme = filter(data, pappus_type == 'Coroniform') %>% select(c(pappus_type,vol,area))
abline(lm(log(area) ~ log(vol), data =   Coroniforme), col = "#4A1079FF", lty = 'dashed',  lwd=2)
Epaposa = filter(data, pappus_type == 'Epappose') %>% select(c(pappus_type,vol,area))
abline(lm(log(area+0.0000001) ~ log(vol), data =   Epaposa), col = "#000004FF", lty = 'dashed',  lwd=2)
abline(lm(log(funcional$S_area_mm2+0.0000001) ~ log(funcional$S_vol_mm3)), col = "#35B77980", lty = 'dashed',  lwd=2)
legend(x = 4, y = -10, legend = 'Functional data set', bty = 'n', cex = 0.5, col ="#35B77980", pch = 18)
dev.off()
