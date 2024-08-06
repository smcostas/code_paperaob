library(tidyverse);library(phytools);library(geiger);library(viridis); library(nlme); library(AICcmodavg)

# Read Tree and Data #######
data <- read.csv('final-df.csv', header = T, row.names = 1, stringsAsFactors = T)
data$pappus_type <- factor(data$pappus_type, levels = c('Setose', 'Paleaceous', 'Aristate', 'Coroniform', 'Epappose'))
subdata <- data[,c('spp','vol','area', 'pappus_type')]
tree2 <- read.newick('final-tree.nwk')
name.check(tree2,data)

subdata$pappus_type <- as.factor(subdata$pappus_type)
subdata$logvol <- log(subdata$vol)
subdata$logarea <- log(subdata$area+0.0000001)

## Without Epappose species #####
subdata2 <- filter(subdata, pappus_type != 'Epappose')
subdata2$pappus_type <- factor(subdata2$pappus_type, 
                               levels = c("Setose","Paleaceous","Aristate","Coroniform"))
tree3 <- drop.tip(tree2, rownames(filter(subdata, pappus_type == 'Epappose')))
length(tree3$tip.label)
pagel <- corPagel(1,tree3) ## vcv matrix using pagel
name.check(tree3,subdata2)

## reorder and create new variable log_pl ######
subdata2 <- subdata2[tree3$tip.label,]
subdata2$log_pl <- log(subdata2$vol/subdata2$area)

# Fitting Pagel's lamda  PGLS #####
fitgls <- gls (logarea ~ logvol*pappus_type , data = subdata2, correlation = pagel)
summary(fitgls)
anova(fitgls, type="marginal")
intervals(fitgls)

## Diagnostic plots #####
plot(fitgls, ylim = c(-1.5,1.5))
qqnorm(fitgls, abline=c(0,1))
plot(fitgls$residuals ~ subdata2$logvol)
plot(fitgls$residuals ~ subdata2$pappus_type)

## Plotting pgls Fig 6 ####
colr_pappus <- setNames(magma(5, end = 0.94, direction = -1) , levels(as.factor(data$pappus_type)))
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

ypred <- predictSE.gls(fitgls, subdata2[,c('logvol','pappus_type')] , se.fit = T)
subdata2$ypred <- ypred$fit
subdata2$se <- ypred$se.fit
p <- ggplot(subdata2, aes(x=logvol, y =  logarea, col = pappus_type)) 
p <- p + geom_ribbon( aes(ymin = ypred-se, ymax = ypred+se, 
                          color = NULL, fill = pappus_type), 
                      alpha = .15) + geom_line(aes(y = ypred), 
                                               linewidth = 1) + geom_point()
p + labs(x= 'log(Vol mm3)', y = 'log(area mm2)', col = 'Pappus type', 
         fill = 'Pappus type') + 
  scale_y_continuous(breaks = seq(-5,5, by = 2.5)) +
  scale_color_manual(values = colr_pappus[1:4], labels = c('Setose', 'Palaceous', 'Aristate', 'Coroniform')) + 
  scale_fill_manual(values = colr_pappus[1:4] , labels = c('Setose', 'Palaceous', 'Aristate', 'Coroniform')) + my_theme +  
  theme(legend.key = element_blank()) 
