library(tidyverse)
library(phytools)
library(geiger)
library(viridis)
library(ggrepel)
## data #####
data <- read.csv('final-df.csv', header = T, row.names = 1)
tree2 <- read.newick('final-tree.nwk')
name.check(tree2,data)

## pANOVA ###########
log_PL <- setNames(data$log_PL, rownames(data))
pappus_type <- setNames(data$pappus_type, rownames(data))
pappus_type <- as.factor(pappus_type)
panova <- phylANOVA(tree2, pappus_type, log_PL, nsim = 10000)
panova
      
## Violin-plot fig S3 ###########

data_b <- data %>% 
  group_by(pappus_type) %>% mutate(n = n()) %>% 
  mutate(label = paste0(pappus_type,'\n(n = ',n,')')) ## change the name to add n for graphics
data_b$label <- as.factor(data_b$label)
levels(data_b$label)
data_b$label <- factor(data_b$label, levels = c("Epappose\n(n = 11)", ## change the order for graphics
                                                "Coroniform\n(n = 17)", 
                                                "Aristate\n(n = 15)", 
                                                "Paleaceous\n(n = 17)", 
                                                "Setose\n(n = 90)"))
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

violin <- ggplot(data = data_b, aes(x = label, y = log_PL))
violin <- violin + geom_violin(trim = F, fill = 'NA' ) + geom_boxplot(
  color = 'grey72', fill = 'NA', width=0.3, alpha = 0.7 ) + 
  geom_jitter(aes(color = pappus_type),
              size = 2,alpha = 0.7,width = 0.15)
violin <- violin + scale_color_manual(values = c("#AA337DFF", "#4A1079FF", "#000004FF", "#F7725CFF", "#FDE2A2FF"))
violin <- violin + #+ coord_flip()
  ylab('log(PL)') + xlab('Pappus type') + my_theme
violin

df_letters <- data.frame( # data frame to print optima
  x = c("Epappose\n(n = 11)", "Coroniform\n(n = 17)", "Aristate\n(n = 15)", 
         "Paleaceous\n(n = 17)", "Setose\n(n = 90)"),
  y = c(23, 7, 7, 2, 2 ),
  label = c('c', 'b', 'b', 'a', 'a')
)


violin + geom_text(data = df_letters, aes(x = x, y = y, label = label))

### plot optima same plot #####

df_regimen <- data.frame( # data frame to print optima
  pt = c("Epappose\n(n = 11)", "Coroniform\n(n = 17)", "Aristate\n(n = 15)", 
         "Paleaceous\n(n = 17)", "Setose\n(n = 90)"),
  optima = c(17.199, 0.046, -0.13, -3.67, -4.918 ) ## model optima
)

violin + geom_point(data = df_regimen, aes(x = pt, y = optima), color = '#238443',
                    size = 5)  + geom_label_repel(data = df_regimen, 
                                                  aes(x = pt, y = optima, label = 'theta'),
                                                  nudge_x = 0.5, nudge_y = 3,
                                                  segment.color = 'black', parse = T,
                                                  alpha = 0.5, segment.linetype = 2)

#####################################
