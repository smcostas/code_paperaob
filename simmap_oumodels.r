library(tidyverse);library(phytools);library(corHMM);library(geiger)
library(beepr);library(viridis);library(OUwie);library(patchwork);library(nlme)
library(reshape2)
# DATA #######
data <- read.csv('final-df.csv', header = T, row.names = 1, stringsAsFactors = T)
paleac <- data %>% filter(pappus_type == 'Paleaceous') %>% select(c('spp', 'Tribe'))
setoso <- data %>% filter(pappus_type == 'Setose') %>% select(c('spp', 'Tribe'))
setoso
# TREE ####
tree2 <- read.newick('final-tree.nwk')
name.check(tree2,data)
data$pappus_type
head(data)
summary(data)

t <- table(data$Tribe, data$pappus_type)
t

##  TO USE LATER ########
data$pappus_type <- factor(data$pappus_type, levels = c('Setose', 'Paleaceous', 'Aristate', 'Coroniform', 'Epappose'))
colr_pappus <- setNames(magma(5, end = 0.94, direction = -1) , levels(as.factor(data$pappus_type))) ### named vector to use in plots 
my_theme <-   theme(
  axis.ticks = element_blank(),
  axis.line = element_line(colour = "grey50"),
  panel.grid = element_line(color = "#b4aea9"),
  panel.grid.minor = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_line(linetype = "dashed"),
  panel.background = element_rect(fill = "white", color = "white"),
  plot.background = element_rect(fill = "white", color = "white")
)
## SEEE THE TREE ####
plot(tree2, ftype = 'off')
nodelabels(text=time,node=1:tree2$Nnode+Ntip(tree2), frame="none",adj=c(1.1,-0.4))
time = as.vector(branching.times(tree2))

# PAPPUS TYPE SIMMAP using hidden states model : three alternative transition models ####
smdata <- data.frame(
  spp = rownames(data),
  pappus_type = data$pappus_type
)
## ALL RATES DIFFERENT ####
ARD_pt <- corHMM(tree2, smdata, rate.cat = 1, model = 'ARD', node.states = 'marginal',
                 root.p = 'yang', n.cores = 9, get.tip.states = T) ## if you run this in windows you can't add n.cores > 1
ARD_pt

## SYMETRIC RATES #####
SYM_pt <- corHMM(tree2, smdata, rate.cat = 1, model = 'SYM', node.states = 'marginal',
                 root.p = 'yang', n.cores = 9, get.tip.states = T)
SYM_pt

## EQUAL RATES ######
ER_pt <- corHMM(tree2, smdata, rate.cat = 1, model = 'ER', node.states = 'marginal',
                root.p = 'yang', n.cores = 9, get.tip.states = T)
ER_pt

## AICc COMPARISON ######
aicc_cerda<-setNames(
  c(ARD_pt$AICc, SYM_pt$AICc, ER_pt$AICc),
  c("ARD","SYM", 'ER'))
aicc_cerda

AICC.W_cerda <- aic.w(aicc_cerda) ## here you can see which model explain better the data

## normalizing the number of stochastic maps according to each model weight #####
nsim<-1000

Nsim<-round(nsim*AICC.W_cerda)
d<-if(sum(Nsim)>nsim) -1 else 1
nsim<-Nsim+d*sample(c(rep(1,abs(nsim-sum(Nsim))),
                      rep(0,length(Nsim)-abs(nsim-sum(Nsim))))) ## to get 1000 trees
nsim ## named vector with number of trees per model 

## Fitting the model to get the simmaps ########
smp_ARD <- makeSimmap(tree = tree2, data = smdata,
                      model =  ARD_pt$solution,
                      rate.cat = 1, nSim = nsim['ARD'],
                      nCores = 9) ## same nCores > 1

smp_SYM <- makeSimmap(tree = tree2, data = smdata,
                      model =  SYM_pt$solution,
                      rate.cat = 1, nSim = nsim['SYM'],
                      nCores = 9)

smp_ER <- makeSimmap(tree = tree2, data = smdata,
                     model =  ER_pt$solution,
                     rate.cat = 1, nSim = nsim['ER'],
                     nCores = 9)
beep('fanfare')

smtrees_final <- c(smp_ARD, smp_SYM, smp_ER)

### Summarizing the Results of the 10000 SIMMAPS #####
sum_sm <- describe.simmap(smtrees_final);beep('fanfare') 
sum_sm
### root probabilities ####
rootP <- as.data.frame(sum_sm$ace[1,])
colnames(rootP) <- 'root_probabilities'
rootP
### to get transitions #####
tr_pappus <- matrix(0,5,5) 
rownames(tr_pappus)<-colnames(tr_pappus)<-levels(as.factor(smdata$pappus_type))
tr_pappus
for (i in 1:1000){
  obj <- describe.simmap(smtrees_final[[i]])
  tr_pappus <- tr_pappus + obj$Tr
  
}


beep()

cambios <- sum(tr_pappus)/1000
cambios_p <-  tr_pappus/1000
prop_cambios <- cambios_p/sum(cambios_p)
sum(prop_cambios[,1]) # aris
sum(prop_cambios[,2]) # coron
sum(prop_cambios[,3]) # epappose
sum(prop_cambios[,4]) # palea
sum(prop_cambios[,5]) # seta


### histogram for transitions from/to setose #####


hist_df <- data.frame(
  setose.other = rep(0,1000),
  other.setose = rep(0,1000),
  setose.epappose = rep(0,1000),
  epappose.setose = rep(0,1000),
  other.epappose = rep(0,1000),
  epappose.other = rep(0,1000)
) 

for (i in 1:1000){
  obj <- describe.simmap(smtrees_final[[i]])
  hist_df$setose.other[i] <- sum(obj$Tr[5,c(1,2,4)])
  hist_df$other.setose[i] <- sum(obj$Tr[c(1,2,4),5])
  hist_df$setose.epappose[i] <- sum(obj$Tr[5,3])
  hist_df$epappose.setose[i] <- sum(obj$Tr[3,5])
  hist_df$other.epappose[i] <- sum(obj$Tr[c(1,2,4), 3])
  hist_df$epappose.other[i] <- sum(obj$Tr[3, c(1,2,4)])
}
hist_df$tree <- seq(1,1000, by = 1)
hist_reshape <- melt(hist_df, id.vars = 'tree', value.name = 'total_change', variable.name='transition')
hist_sum <- hist_reshape %>% select(total_change, transition) %>% group_by(transition) %>% 
  summarise_all(lst(mean))
head(hist_sum)

ggplot(hist_reshape, aes(x = total_change, color = transition, fill = transition)) +  
  geom_histogram(position = 'identity', alpha = 0.5) + 
  geom_vline(data = hist_sum, aes(xintercept=mean, color = transition),
             linetype="dashed") +  scale_y_continuous(expand = c(0,0), breaks = c(0,250,500,750,1000), limits = c(0,1050)) +
  scale_fill_manual(values = c(viridis(4, alpha = 0.5, begin = 0.5, direction = -1, option = "A"), 
                               viridis(1, begin = 0.05), 'black')) + 
  scale_color_manual(values = c(viridis(4, alpha = 0.5, begin = 0.5, direction = -1, option = "A"), 
                                viridis(1, begin = 0.05), 'black')) +
  labs(x= 'total number of transitions', y = 'frequency across 1000 simulations') + 
  my_theme

## Plotting simmaps #####

### some changes to tribes ####
data$Tribe <- as.character(data$Tribe) ## to add tribes without problems...
data["Spilanthes-sp", "Tribe"] <- 'Heliantheae'
data["Guizotia-abissynica", "Tribe"] <- 'Milerieae'
data["Calea-candolleana", "Tribe"] <- 'Neurolaeneae'
data$Tribe <- as.factor(data$Tribe)

nodes <- levels(as.factor(data$Tribe))

### plotting with time branches and tribes #####

time = formatC(as.vector(branching.times(tree2)), format = 'f', flag='0', digits = 2)
#pdf('sm_pappustype_tribes_time.pdf', width = 12, height = 9)
par(mar = c(1.2,1.2,1.2,1.2))
plot(sum_sm, colors = colr_pappus,lwd = 2, offset= 1,ftype = 'off', cex = c(0.35,0.35), type="fan")
add.simmap.legend(leg=levels(data$pappus_type), colors=colr_pappus,prompt=FALSE,x=-120,
                  y=-40,fsize=0.8)
for (i in 1:length(nodes)){
  node <- filter(data, Tribe == nodes[i])
  node <- rownames(node)
  if (length(node)>1) {
    node <- getMRCA(tree2, node)
    arc.cladelabels(tree2,nodes[i],node,1.03,1.06,
                    mark.node=F)
  }
}
nodelabels(text=time,node=1:tree2$Nnode+Ntip(tree2), frame="none", adj=c(1.1,-0.4),cex=0.5)
dev.off()
### not so beautiful to see in R, some modifications were made in inkscape





## SIMMAP pappus  presence absence #### this was needed for evolutionary models  but wasn't discussed in the paper #########
p.a <- factor(levels = c('present', 'absent'))
for (i in 1:length(data$pappus_type)){
  if 
  (data$pappus_type[i] == 'Epappose') {p.a[i] <- 'absent'} 
  else
  {p.a[i] <- 'present'}  
}
p.a

smdata_pa <- data.frame(
  spp = rownames(data),
  p.a = p.a
)
smdata_pa
### models  #########
ARD_pa <- corHMM(tree2, smdata_pa, rate.cat = 1, model = 'ARD', node.states = 'marginal',
                 root.p = 'yang', n.cores = 9, get.tip.states = T)
ARD_pa

SYM_pa <- corHMM(tree2, smdata_pa, rate.cat = 1, model = 'SYM', node.states = 'marginal',
                 root.p = 'yang', n.cores = 9, get.tip.states = T)
SYM_pa

ER_pa <- corHMM(tree2, smdata_pa, rate.cat = 1, model = 'ER', node.states = 'marginal',
                root.p = 'yang', n.cores = 9, get.tip.states = T)
ER_pa

aicc_pa<-setNames(
  c(ARD_pa$AICc, SYM_pa$AICc, ER_pa$AICc),
  c("ARD","SYM", 'ER'))
aicc_pa

AICC.W_pa <- aic.w(aicc_pa)


nsim2<-1000

Nsim2<-round(nsim2*AICC.W_pa)
d<-if(sum(Nsim2)>nsim2) -1 else 1
nsim2<-Nsim2+d*sample(c(rep(1,abs(nsim2-sum(Nsim2))),
                        rep(0,length(Nsim2)-abs(nsim2-sum(Nsim2)))))
nsim2
### run the simmaps ######
sm_pa_ARD <- makeSimmap(tree = tree2, data = smdata_pa,
                        model =  ARD_pa$solution,
                        rate.cat = 1, nSim = nsim2['ARD'],
                        nCores = 9)

sm_pa_SYM <- makeSimmap(tree = tree2, data = smdata_pa,
                        model =  SYM_pa$solution,
                        rate.cat = 1, nSim = nsim2['SYM'], 
                        n.cores = 9)

sm_pa_ER <- makeSimmap(tree = tree2, data = smdata_pa,
                       model =  ER_pa$solution,
                       rate.cat = 1, nSim = nsim2['ER'],
                       nCores = 9)
beep()

smtrees_pa_final <- c(sm_pa_ARD, sm_pa_SYM, sm_pa_ER)

sum_sm_pa <- describe.simmap(smtrees_pa_final);beep('fanfare')
sum_sm_pa
### root probabilities ######
rootP_pa <- as.data.frame(sum_sm_pa$ace[1,])
colnames(rootP_pa) <- 'root_probabilities'
rootP_pa


### Plot simmaps ####
colr_pa<- setNames(viridis(2, end = 1, direction = -1) , levels(as.factor(smdata_pa$p.a)))
plot(sum_sm_pa, fsize=0.4, ftype="off", colors=colr_pa,
     ylim = c(-2, Ntip(tree2)), cex = c(0.4,0.4))
#pdf('sm_presence-absence.pdf', width = 12, height = 8)
plot(sum_sm_pa, colors = colr_pa,lwd = 2, offset= 1,fsize=0.6, ftype="off", cex = c(0.35,0.35), type="fan")
add.simmap.legend(leg=levels(as.factor(smdata_pa$p.a)), colors=colr_pa,prompt=FALSE,x=-120,
                  y=-40,fsize=0.8)
dev.off()

# Modelling evolution log(PL) with BM, OU and MULTI OU using OUwie #####
## preparing data to fit models #####
### Pappus types #######
ouwie_data_pt <- data.frame(
  spp = rownames(data),
  pappus_type = data$pappus_type,
  log_PL = data$log_PL
)
ouwie_data_pt
levels(ouwie_data_pt$pappus_type)

### Presence absence #######
ouwie_data_pa <- data.frame(
  spp = rownames(data),
  p.a = p.a,
  log_PL = data$log_PL
)
ouwie_data_pa

## Brownian motion model ####
fitBM <- OUwie(smtrees_final[[2]], ouwie_data, model="BM1", simmap.tree= T, root.station=FALSE)  
## OU model ########
fitOU <- OUwie(smtrees_final[[2]], ouwie_data, model="OU1", simmap.tree= T, root.station=FALSE)

## MODEL COMPARISONS
### df empty to fill with comparisons #######
comp <- data.frame(
  aiccBM = rep(fitBM$AICc,1000),
  aiccOU = rep(fitOU$AICc,1000),
  aiccMOU_pt = rep(0,1000), # valores obtenidos de aicc en cada ciclo (cada mapeo est) sera cargado en este vector
  aiccMOU_pa = rep(0,1000),
  model = factor(rep('OU',1000), levels = c('MOU_pt', 'MOU_pa', 'OU', 'BM')) ## este vector contendrá el modelo seleccionado 
)
levels(comp$model)

### comparisons #######
#### this could be very slow, take care
for (i in 1:1000){
  print(i)
  fitOUM_pt <- OUwie(smtrees_final[[i]], ouwie_data_pt, model="OUM", simmap.tree= T, root.station=FALSE)  ### note that each i run different tree so oyu get different aicc depending on the simmap tree regimes
  fitOUM_pa <- OUwie(smtrees_pa_final[[i]], ouwie_data_pa, model="OUM", simmap.tree= T, root.station=FALSE)
  comp$aiccMOU_pt[i] <- fitOUM_pt$AICc
  comp$aiccMOU_pa[i] <- fitOUM_pa$AICc
  if ((comp$aiccMOU_pt[i] <= comp$aiccMOU_pa[i]) & (comp$aiccMOU_pt[i] <= comp$aiccOU[i])) {
    comp$model[i] <- 'MOU_pt'
  }
  else
    if ((comp$aiccMOU_pa[i] < comp$aiccMOU_pt[i]) & (comp$aiccMOU_pa[i] < comp$aiccOU[i])) {
      comp$model[i] <- 'MOU_pa'
    }
  print(comp$model[i])
}
beep('fanfare')
### boxplot with aicc comparisons #####
AICc <- c(comp$aiccMOU_pt, comp$aiccMOU_pa, comp$aiccOU, comp$aiccBM)
aicc_comp <- data.frame(
  model = Model,
  aicc = AICc
)
boxplot(AICc ~ Model, data = aicc_comp)

### complete comparison between MULTI OU MODELS #####
fmodel <- factor(rep('MOU_pt', 1000*1000), levels = c('MOU_pt', 'MOU_pa'))
i = 1
for (c in 1:1000) {
  for (d in 1:1000){
    if (comp$aiccMOU_pa[c] < comp$aiccMOU_pt[d]) {fmodel[i] <-  'MOU_pa'}
    i = i+1
    print(i)
  }
}
beep()
table(fmodel)

### plotting the tree which get the aicc min model with node probabilities #######
colnames(comp) <- c('BM','OU1', 'OU5', 'OU2', 'model selected')
index_aiccmin <- which.min(comp$OU5)
pappus <- setNames(data$pappus_type, rownames(data))
#pdf('simmap_minaic.pdf', width = 12, height = 8)
plot(smtrees_final[[index_aiccmin]], colors = colr_pappus,lwd = 2, offset= 1,fsize=0.6, ftype="off", type="fan")
tiplabels(pie = to.matrix(pappus[tree2$tip.label],levels(pappus)),piecol=colr_pappus, cex=0.35)
colr_ace<- setNames(c("#AA337DFF", "#4A1079FF", "#000004FF", "#F7725CFF", "#FDE2A2FF") , colnames(sum_sm$ace))# lo tengo que hacer manual
nodelabels(pie=sum_sm$ace,piecol= colr_ace, cex=0.35)
add.simmap.legend(leg=levels(data$pappus_type), colors=colr_pappus,prompt=FALSE,x=-120,
                  y=-40,fsize=0.8)
dev.off()

### Confidence interval for the best model ####
final_model <- OUwie(smtrees_final[[index_aiccmin]], ouwie_data_pt, model="OUM", 
                     simmap.tree= T, root.station=FALSE, get.root.theta=T)


## parametros del modelo
theta <- final_model$theta[2:6,1]
theta0 <- final_model$theta[1,1]
alpha <- as.vector(final_model$solution[1,])
sigma <- as.vector(final_model$solution[2,])
ci_boot <- OUwie.boot(smtrees_final[[index_aiccmin]], ouwie_data_pt, 
                      model = 'OUM', nboot = 1000, alpha = alpha, 
                      sigma.sq = sigma, theta = theta, theta0 = theta0, 
                      simmap.tree = T) 
head(ci_boot)
ci_inf_alpha <- min(ci_boot[,1])
ci_sup_alpha <- max(ci_boot[,1])
ci_inf_sigma.sq <- min(ci_boot[,6])
ci_sup_sigma.sq <- max(ci_boot[,6])
ci_inf_theta_aris <- min(ci_boot[,'theta_Aristate'])
ci_sup_theta_aris <- max(ci_boot[,'theta_Aristate'])
ci_inf_theta_coron <- min(ci_boot[,'theta_Coroniform'])
ci_sup_theta_coron <- max(ci_boot[,'theta_Coroniform'])
ci_inf_theta_epap <- min(ci_boot[,'theta_Epappose'])
ci_sup_theta_epap <- max(ci_boot[,'theta_Epappose'])
ci_inf_theta_palea <- min(ci_boot[,'theta_Paleaceous'])
ci_sup_theta_palea <- max(ci_boot[,'theta_Paleaceous'])
ci_inf_theta_set <- min(ci_boot[,'theta_Setose'])
ci_sup_theta_set <- max(ci_boot[,'theta_Setose'])

table_ci <- data.frame(
  inf_alpha = ci_inf_alpha,
  sup_alpha = ci_sup_alpha,
  inf_sigama.sq = ci_inf_sigma.sq,
  sup_sigma.sq = ci_sup_sigma.sq,
  inf_theta_set = ci_inf_theta_set,
  sup_theta_set = ci_sup_theta_set,
  inf_theta_palea = ci_inf_theta_palea,
  sup_theta_palea = ci_sup_theta_palea,
  inf_theta_aris = ci_inf_theta_aris,
  sup_theta_aris = ci_sup_theta_aris,
  inf_theta_coron = ci_inf_theta_coron,
  sup_theta_coron = ci_sup_theta_coron,
  inf_theta_epap = ci_inf_theta_epap,
  sup_theta_epep = ci_sup_theta_epap
)
table_ci

## trees of Fig 2 ######
#pdf('evolutive-model-scenarios.pdf', height =10, width = 15)
par(mfcol=c(1,4))
plotTree(tree2,lwd = 1.5, offset= 0.5,fsize=0.4, ftype="off")
plotTree(tree2,lwd = 1.5, offset= 0.5,fsize=0.4, ftype="off")
plot(smtrees_pa_final[[306]], colors = colr_pa,lwd = 2, offset= 0.5,fsize=0.4, ftype="off")
plot(smtrees_final[[index_aiccmin]], colors = colr_pappus,lwd = 2, offset= 0.5,fsize=0.4, ftype = 'off')
dev.off()

## figure 5 phenogram with evolutionary optima ########
### preparing data #################
funcional <- read.csv("funcional.csv", header = T) ## needed to plot the mean
rownames(funcional) <- funcional$spp
funcional <- filter(funcional, S_long_cerda_max > 0)
p_type <- setNames(data$pappus_type, rownames(data))
p_type <- as.factor(p_type)
p_type
log_PL <- setNames(data$log_PL, rownames(data))
labels <-setNames(seq(-15,25,by= 40/(Ntip(tree2)-1)),names(sort(log_PL)))
log_pl_funcional = funcional$S_vol_mm3/(funcional$S_area_mm2+0.0000001)
log_pl_funcional = log(log_pl_funcional)

#pdf(file = 'phenogram-optima_fm.pdf', width = 12, height = 8 )
par(mai = c(0.3,0.3,0.3,0.3), lwd = 0.3, cex = 0.5 )
#obj <- phenogram(smtrees_final[[index_aiccmin]],log_PL,lwd=3, ylim = c(-15,25),colors=
#                   colr_pappus,fsize=0.4, label.pos=sample(labels), ftype="off", spread.range = c(-15,25), link = 8)
obj <- phenogram(tree2,log_PL,lwd=3, ylim = c(-15,25), colors = 'grey6',lwd = 0.1 ,fsize=0.4, label.pos=sample(labels), ftype="off", spread.range = c(-15,25), link = 8)
# lineas para pintar los tiplabels
#label.colors<-setNames(colr_pappus[p_type],names(p_type))
#lastPP<-get("last_plot.phylo",envir=.PlotPhyloEnv) 
#shadowtext(x=obj,labels=rownames(obj),
#           col=label.colors[rownames(obj)],bg="white",pos=4,
#           offset=lastPP$label.offset, cex = 0.4)
nodelabels(pie=sum_sm$ace,piecol=setNames(c("#AA337DFF", "#4A1079FF","#000004FF","#F7725CFF","#FDE2A2FF"),colnames(sum_sm$ace)), cex=0.25)
tiplabels(pie = to.matrix(pappus[tree2$tip.label],levels(pappus)),piecol=colr_pappus, cex=0.2, adj = c(1.2,0.55))
add.simmap.legend(x = 84, y= -10, fsize = 0.4, colors=colr_cerda,prompt=FALSE)
arrows(x0 = -3.5, x1 = -2, y0 = -4.91804428, y1 = -4.91804428, length = 0.1, col = colr_pappus[1],lwd = 3) ## setose optimum
arrows(x0 = -3.5, x1 = -2, y0 = -3.66978451, y1 = -3.66978451, length = 0.1, col = colr_pappus[2],lwd = 3) ## paleaceous optimum
arrows(x0 = -3.5, x1 = -2, y0 = -0.12777672, y1 = -0.12777672 , length = 0.1, col = colr_pappus[3],lwd = 3) ## aristate optimum
arrows(x0 = -3.5, x1 = -2, y0 = 0.04571855, y1 = 0.04571855 , length = 0.1, col = colr_pappus[4],lwd = 3) ## coroniform optimum
arrows(x0 = -3.5, x1 = -2, y0 = 17.19908256, y1 = 17.19908256 , length = 0.1, col = colr_pappus[5],lwd = 3) ## epappose optimum
arrows(x0 = -3.5, x1 = -2, y0 = mean(log_pl_funcional), y1 = mean(log_pl_funcional) , length = 0.1, col = "#35B77980",lwd = 3) ## mean funcional data
dev.off()


