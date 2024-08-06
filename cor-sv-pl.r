library(tidyverse)
library(ggpubr)

funcional <- read.csv("funcional.csv", header = T)
rownames(funcional) <- funcional$spp
funcional <- filter(funcional, S_long_cerda_max > 0)
summary(funcional)
# Getting new features ####
## Volume ####
funcional$S_vol_mm3 <- pi*(funcional$S_ancho_aq/2)^2*(funcional$S_long_aq/2)

##  Pappus area ####
funcional$S_area_mm2 <- funcional$S_long_cerda_max^2*pi

## Plume loading ######
funcional$S_PL <- funcional$S_vol_mm3/funcional$S_area_mm2


funcional <- funcional[!rownames(funcional)%in%"Cynara_cardubnculus",] # cleaning this spp.

## Pearson correlation test #######
cor.test(funcional$B_SV,funcional$S_PL)
 
# Supplementary plot fig S2. #######
p1 <- ggplot(funcional, aes(x = S_PL, y = B_SV)) + 
  geom_smooth(method = "lm", formula = y~x, fill = "grey", color = "black") + 
  geom_point(color = "black")
p2 <- p1 + xlab ("Plume loading ") + ylab ("Settling Velocity (m/s)") + 
  scale_x_continuous(expand = c(0,0.001)) + scale_y_continuous(limits = c(0.18,2.65), expand = c(0,0)) +
  theme(panel.background = element_rect(fill = NA, colour = "white"), 
        axis.line = element_line (color = "black"),
        axis.text = element_text(color = "black")
  )
p2 + stat_cor(label.x = 0.05, label.y = 2.3)
dev.off()

########################################################################





