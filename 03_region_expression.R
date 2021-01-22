suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(tools)
  library(magrittr)
  library(tibble)
  library(plyr)
  library(pheatmap)
  library(readxl)
  library(gdata)
  library(biomaRt)
  library(data.table)
  library(pander)
  library(tidyr)
  library(viridis)
  library(RColorBrewer)
  library(gridExtra)})




path.data <- path.expand("BrainSpan/Kang")

ZMYND8_data <- read.csv("ZMYND8_metadata.csv", header = TRUE, row.names = 1)



### Region expression in prenatal period

pren_R<- as.data.frame(table(ZMYND8_data$structure_acronym[grep("pcw", ZMYND8_data$age)]))
pren_R_structures <- pren_R %>% 
  dplyr::filter(Freq >=10)

PreN_fig <- ZMYND8_data %>%  
  dplyr::filter(period == "prenatal") %>%
  dplyr::filter(structure_acronym %in% as.character(pren_R_structures$Var1)) %>%
  ggplot(aes(structure_acronym, ZMYND8_expression))  + 
  geom_boxplot(aes(fill=structure_acronym), alpha=0.65)  + 
  geom_jitter(size=2.5) +
  scale_color_viridis(discrete = TRUE, option = "magma")+ theme_bw() +
  scale_fill_viridis(discrete = TRUE, direction = 1)  +
  xlab("Region") + ylab("RPKM") + ylim(0, 30) +
  theme(legend.position = "none", 
        axis.text=element_text(size=11, face="bold", color="black"), 
        axis.title = element_text(siz=12, face="bold", color="black"), 
        axis.text.x = element_text(angle=90)) + ggtitle("Prenatal Expression")


### Region expression in postnatal period 

postn_R <- as.data.frame(table(ZMYND8_data$structure_acronym[-grep("pcw", ZMYND8_data$age)]))

postn_R_structures <- postn_R %>% 
  dplyr::filter(Freq >= 10)

PN_fig <- ZMYND8_data %>%  
  dplyr::filter(period == "postnatal") %>%
  dplyr::filter(structure_acronym %in% as.character(postn_R_structures$Var1)) %>%
  ggplot(aes(structure_acronym, ZMYND8_expression))  + 
  geom_boxplot(aes(fill=structure_acronym), alpha=0.65)  + 
  geom_jitter(size=2.5) +
  scale_color_viridis(discrete = TRUE, option = "magma")+ theme_bw() +
  scale_fill_viridis(discrete = TRUE, direction = 1)  +
  xlab("Region") + ylab("RPKM") + ylim(0, 30) +
  theme(legend.position = "none", 
        axis.text=element_text(size=11, face="bold", color="black"), 
        axis.title = element_text(siz=12, face="bold", color="black"), 
        axis.text.x = element_text(angle=90)) + ggtitle("Postnatal Expression")




pdf("ZMYND8_region_expression_brainSpan.pdf",width =14, height =12)
grid.arrange(PreN_fig, PN_fig, ncol=2, nrow=2)
dev.off()



### Linear model 


##### prenatal 
prenatal_structures <- ZMYND8_data %>%  
  dplyr::filter(period == "prenatal") %>%
  dplyr::filter(structure_acronym %in% as.character(pren_R_structures$Var1)) 


model_pren=lm(prenatal_structures$ZMYND8_expression ~ prenatal_structures$structure_acronym); aov(model_pren)
summary(aov(model_pren))


##### postnatal

postnatal_structures <- ZMYND8_data %>%  
  dplyr::filter(period == "postnatal") %>%
  dplyr::filter(structure_acronym %in% as.character(postn_R_structures$Var1))

model_postn=lm(postnatal_structures$value ~ postnatal_structures$structure_acronym)
summary(aov(model_postn))

