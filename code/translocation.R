Packages <- c(
  "ape", "lme4", "ggtern", "nlme", "lmerTest", "openxlsx", 
  "edgeR", "ggplot2", "vegan", "tidyverse", "ggalt", 
  "patchwork", "NST", "ggpmisc", "reshape2", "ggprism", 
  "ggpubr", "ggalluvial", "dplyr", "ggh4x", "car", "NST",
  "sjPlot", "performance", "MuMIn", "mgcv", "emmeans", 
  "multcomp", "stringr", "gamm4", "Matrix", "ggthemes", 
  "tidyr", "glmm.hp", "partR2", "corrplot", "gratia", 
  "ggstar", "ggmagnify", "ggbreak", "picante", "RColorBrewer", 
  "ggforce", "grid", "gridExtra", "linkET", "cowplot"
)
#install.packages(setdiff(Packages, rownames(installed.packages())))
lapply(Packages, library, character.only = TRUE)

setwd("D:/desktop/ASV/")

###########common functions################
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = sd(x[[col]]/sqrt(length(x[[col]]))))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

####translocation
diversity <- read.xlsx("field_translocation.xlsx", sheet = "group",colNames = T, rowNames = T)
diversity$Site <- as.factor(diversity$Site)
diversity$Popu <- as.factor(diversity$Popu)
diversity <- subset(diversity,Popu != "L" & Time != "Jun")

###alpha diversity in translocation
model <- lmer(log(ASVs.of.total) ~ Compartment * Site * Species * Popu + (1|pair), 
              na.action= na.omit,data = diversity)
car::Anova(model, type = 3)
ranova(model)

model <- lmer(ASVs ~ Compartment * Site * Species * Popu + (1|pair), 
              na.action= na.omit,data = diversity)
car::Anova(model, type = 3)
ranova(model)

model <- lmer(read.of.shared ~ Compartment * Site * Species * Popu + (1|pair), 
              na.action= na.omit,data = diversity)
car::Anova(model, type = 3)
ranova(model)

###beta diversity in translocation
###permanova analysis
rm(list=ls())
total <- readRDS("trans_bray.rds")
total <- data.frame(as.matrix(total))
group <- read.xlsx("field_translocation.xlsx", sheet = "group",colNames = T, rowNames = T)
group$Site <- as.factor(group$Site)
group$Popu <- as.factor(group$Popu)
group1 <- subset(group,Popu != "1590" & Time != "Jun")
data <- total[rownames(total) %in% rownames(group1),]
wunifrac <- data[,colnames(data) %in% rownames(group1),]
with(group1, adonis2(wunifrac ~ pair + Compartment * Site * Species * Popu, 
                     data = group1, permutations = 999, by = "term"))

total <- readRDS("trans_shared_bray.rds")
total <- data.frame(as.matrix(total))
group <- read.xlsx("field_translocation.xlsx", sheet = "group",colNames = T, rowNames = T)
group$Site <- as.factor(group$Site)
group$Popu <- as.factor(group$Popu)
group1 <- subset(group,Popu != "1590" & Time != "Jun")
data <- total[rownames(total) %in% rownames(group1),]
wunifrac <- data[,colnames(data) %in% rownames(group1),]
with(group1, adonis2(wunifrac ~ pair + Compartment * Site * Species * Popu, 
                     data = group1, permutations = 999, by = "term"))

total <- readRDS("trans_unique_bray.rds")
total <- data.frame(as.matrix(total))
group <- read.xlsx("field_translocation.xlsx", sheet = "group",colNames = T, rowNames = T)
group$Site <- as.factor(group$Site)
group$Popu <- as.factor(group$Popu)
group1 <- subset(group,Popu != "1590" & Time != "Jun")
data <- total[rownames(total) %in% rownames(group1),]
wunifrac <- data[,colnames(data) %in% rownames(group1),]
with(group1, adonis2(wunifrac ~ pair + Compartment * Site * Species * Popu, 
                     data = group1, permutations = 999, by = "term"))

###plant performance in tranlocation
plant <- read.xlsx("field_translocation.xlsx", sheet = "plant",colNames = T, rowNames = T)
group1 <- subset(plant, Elevational != "0")
plant$Elevational <- factor(plant$Elevational)
model <- lmer(lnHome_CK ~ Elevational * Species + (1|Popu), 
              na.action= na.omit,data = group1)
car::Anova(model, type = 3)
ranova(model)

############################################################################
#################################graph######################################
############################################################################
##############Fig3#######################
group <- read.xlsx("field_translocation.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)
group$ASVs.of.total <- log10(group$ASVs.of.total)

# Define general drawing functions
create_plot <- function(data, y_var, y_label, fill_values, color_values, y_breaks = NULL) {
  p <- ggplot(data, aes(x = Compartment, group = Type)) +
    geom_errorbar(aes_string(ymin = paste(y_var, "- se"), ymax = paste(y_var, "+ se"), color = "Type"),
                  position = position_dodge(width = 0.6), width = 0, size = 0.5) +
    geom_point(aes_string(y = y_var, color = "Type"), 
               position = position_dodge(width = 0.6), size = 2) +
    labs(x = "Groups", y = y_label) +
    scale_fill_manual(values = fill_values) +
    scale_color_manual(values = color_values) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid = element_blank(), 
          legend.position = "none",
          axis.title = element_text(color = 'black', size = 14),
          axis.ticks = element_line(color = 'black'),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = 'black', size = 11))
  
  if (!is.null(y_breaks)) {
    p <- p + scale_y_break(y_breaks[1:2], space = 0.001, scales = 0.5)
  }
  
  return(p)
}

# OTUs.of.total for high group
high <- subset(group, Type %in% c("HH", "HM", "HL"))
df1 <- data_summary(high, varname = "ASVs.of.total", groupnames = c("Compartment", "Type"))
df1$Type <- factor(df1$Type, levels = c("HH", "HM", "HL"))

p1 <- create_plot(df1, "ASVs.of.total", "No. of total ASVs (log)", 
                  c("#B7B6B6", "#70D3ED", "#B79F00"), c("#B7B6B6", "#70D3ED", "#B79F00"), c(2.0, 2.7))

p1

# OTUs.of.total for mid group
mid <- subset(group, Type %in% c("MM", "ML"))
df2 <- data_summary(mid, varname = "ASVs.of.total", groupnames = c("Compartment", "Type"))
df2$Type <- factor(df2$Type, levels = c("MM", "ML"))
p2 <- create_plot(df2, "ASVs.of.total", "No. of total ASVs (log)", 
                  c("#E09513", "#9962EF"), c("#E09513", "#9962EF"), c(2.0, 2.7))
p2

df1_shared <- data_summary(high, varname = "ASVs", groupnames = c("Compartment", "Type"))
df1_shared$Type <- factor(df1_shared$Type, levels = c("HH", "HM", "HL"))
p3 <- create_plot(df1_shared, "ASVs", "% ASVs shared", 
                  c("#B7B6B6", "#70D3ED", "#B79F00"), c("#B7B6B6", "#70D3ED", "#B79F00"), c(2.5, 5))
p3
# % OTUs shared for mid group
df2_shared <- data_summary(mid, varname = "ASVs", groupnames = c("Compartment", "Type"))
df2_shared$Type <- factor(df2_shared$Type, levels = c("MM", "ML"))
p4 <- create_plot(df2_shared, "ASVs", "% ASVs shared", 
                  c("#E09513", "#9962EF"), c("#E09513", "#9962EF"), c(2.4, 7))
p4
# Reads for high group
df1_reads <- data_summary(high, varname = "read.of.shared", groupnames = c("Compartment", "Type"))
df1_reads$Type <- factor(df1_reads$Type, levels = c("HH", "HM", "HL"))
p5 <- create_plot(df1_reads, "read.of.shared", "RA of shared ASVs (%)", 
                  c("#B7B6B6", "#70D3ED", "#B79F00"), c("#B7B6B6", "#70D3ED", "#B79F00"))
p5
# Reads for mid group
df2_reads <- data_summary(mid, varname = "read.of.shared", groupnames = c("Compartment", "Type"))
df2_reads$Type <- factor(df2_reads$Type, levels = c("MM", "ML"))
p6 <- create_plot(df2_reads, "read.of.shared", "RA of shared ASVs (%)", 
                  c("#E09513", "#9962EF"), c("#E09513", "#9962EF"))
p6


#############Fig4#############
###entire
ASV <- readRDS("trans_Jun_bray.rds")
group <- read.xlsx("field_translocation.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)
ASV1 <- data.frame(as.matrix(ASV))

perform_pcoa_CK <- function(group_data, title_suffix = "") {
  data <- ASV1[rownames(ASV1) %in% rownames(group_data),]
  rhizos1 <- data[, colnames(data) %in% rownames(group_data)]
  
  pcoa <- cmdscale(rhizos1, k = 2, eig = TRUE)
  poi <- as.data.frame(pcoa$points)
  colnames(poi) <- c("PCoA1", "PCoA2")
  poi$SampleID <- rownames(poi)
  
  sample_site <- cbind(poi, group_data[, c(3, 4, 5, 6, 8)])
  sample_site$Site <- factor(sample_site$Site)
  
  df1 <- data_summary(sample_site, varname = "PCoA1", groupnames = c("Type"))
  colnames(df1) <- c("Type", "mean_PCoA1", "sd_PCoA1", "se_PCoA1")
  
  df2 <- data_summary(sample_site, varname = "PCoA2", groupnames = c("Type"))
  colnames(df2) <- c("Type", "mean_PCoA2", "sd_PCoA2", "se_PCoA2")
  
  summary_df <- cbind(df1, df2)[, -5]
  
  summary_df$Type <- factor(summary_df$Type, levels = c("HH", "HM", "HL", "MM", "ML","CKH","CKM"))
  
  p <- ggplot() + 
    geom_point(sample_site, mapping = aes(x = PCoA1, y = PCoA2, fill = Type, shape = Species), size = 1.5, alpha = 0.8,stroken = 0.5) +
    geom_errorbarh(summary_df, mapping = aes(x = mean_PCoA1, y = mean_PCoA2, 
                                             xmin = mean_PCoA1 + se_PCoA1, 
                                             xmax = mean_PCoA1 - se_PCoA1, color = Type), height = 0, size = 0.5) +
    geom_errorbar(summary_df, mapping = aes(x = mean_PCoA1, y = mean_PCoA2, 
                                            ymin = mean_PCoA2 + se_PCoA2, 
                                            ymax = mean_PCoA2 - se_PCoA2, color = Type), width = 0, size = 0.5) +
    geom_point(summary_df, mapping = aes(x = mean_PCoA1, y = mean_PCoA2, color = Type), size = 4) +
    scale_shape_manual(values = c(21, 22, 24)) +
    scale_fill_manual(values = c(HH = "#717171", HL = "#B79F00", HM = "#4BCAE8", LL = "#1D4899", ML = "#9962EF", MM = "#DF9413", CKH = "black", CKM = "red")) +
    scale_color_manual(values = c(HH = "#717171", HL = "#B79F00", HM = "#4BCAE8", LL = "#1D4899", ML = "#9962EF", MM = "#DF9413", CKH = "black", CKM = "red")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed", size = 0.3) +
    geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.3) +
    scale_x_continuous(labels = scales::label_comma(accuracy = 0.01)) +
    scale_y_continuous(labels = scales::label_comma(accuracy = 0.01)) +
    labs(x = paste("PCoA1 (", format(100 * pcoa$eig[1] / sum(pcoa$eig), digits = 3), "%)", sep = ""), 
         y = paste("PCoA2 (", format(100 * pcoa$eig[2] / sum(pcoa$eig), digits = 3), "%)", sep = ""),
         title = paste("PCoA Plot", title_suffix)) +
    geom_magnify(from = c(-0.25, -0.05, -0.1, 0.1),
                 to = c(-0.1, -0.2, 0.1, -0.1))+
    theme(panel.background = element_rect(fill='white', colour='black'),
          panel.grid=element_blank(), 
          legend.position = "right",
          axis.title = element_text(color='black',size=14),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.text=element_text(colour='black',size=11))
  
  return(p)
}

# phyllosphere
group1 <- subset(group, Compartment == "L" & (Type %in% c("HH", "HM", "HL", "MM", "ML","CKH","CKM")))
p4a_1 <- perform_pcoa_CK(group1, "Phyllosphere")

perform_pcoa_CK <- function(group_data, title_suffix = "") {
  data <- ASV1[rownames(ASV1) %in% rownames(group_data),]
  rhizos1 <- data[, colnames(data) %in% rownames(group_data)]
  
  pcoa <- cmdscale(rhizos1, k = 2, eig = TRUE)
  poi <- as.data.frame(pcoa$points)
  colnames(poi) <- c("PCoA1", "PCoA2")
  poi$SampleID <- rownames(poi)
  
  sample_site <- cbind(poi, group_data[, c(3, 4, 5, 6, 8)])
  sample_site$Site <- factor(sample_site$Site)
  
  df1 <- data_summary(sample_site, varname = "PCoA1", groupnames = c("Type"))
  colnames(df1) <- c("Type", "mean_PCoA1", "sd_PCoA1", "se_PCoA1")
  
  df2 <- data_summary(sample_site, varname = "PCoA2", groupnames = c("Type"))
  colnames(df2) <- c("Type", "mean_PCoA2", "sd_PCoA2", "se_PCoA2")
  
  summary_df <- cbind(df1, df2)[, -5]
  
  summary_df$Type <- factor(summary_df$Type, levels = c("HH", "HM", "HL", "MM", "ML","CKM","CKH"))
  
  p <- ggplot() + 
    geom_point(sample_site, mapping = aes(x = PCoA1, y = PCoA2, fill = Type, shape = Species), size = 1.5, alpha = 0.8,stroken = 0.5) +
    geom_errorbarh(summary_df, mapping = aes(x = mean_PCoA1, y = mean_PCoA2, 
                                             xmin = mean_PCoA1 + se_PCoA1, 
                                             xmax = mean_PCoA1 - se_PCoA1, color = Type), height = 0, size = 0.5) +
    geom_errorbar(summary_df, mapping = aes(x = mean_PCoA1, y = mean_PCoA2, 
                                            ymin = mean_PCoA2 + se_PCoA2, 
                                            ymax = mean_PCoA2 - se_PCoA2, color = Type), width = 0, size = 0.5) +
    geom_point(summary_df, mapping = aes(x = mean_PCoA1, y = mean_PCoA2, color = Type), size = 4) +
    scale_shape_manual(values = c(21, 22, 24)) +
    scale_fill_manual(values = c(HH = "#717171", HL = "#B79F00", HM = "#4BCAE8", LL = "#1D4899", ML = "#9962EF", MM = "#DF9413", CKH = "black", CKM = "red")) +
    scale_color_manual(values = c(HH = "#717171", HL = "#B79F00", HM = "#4BCAE8", LL = "#1D4899", ML = "#9962EF", MM = "#DF9413", CKH = "black", CKM = "red")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed", size = 0.3) +
    geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.3) +
    scale_x_continuous(labels = scales::label_comma(accuracy = 0.001)) +
    scale_y_continuous(labels = scales::label_comma(accuracy = 0.001)) +
    labs(x = paste("PCoA1 (", format(100 * pcoa$eig[1] / sum(pcoa$eig), digits = 3), "%)", sep = ""), 
         y = paste("PCoA2 (", format(100 * pcoa$eig[2] / sum(pcoa$eig), digits = 3), "%)", sep = ""),
         title = paste("PCoA Plot", title_suffix)) +
    geom_magnify(from = c(-0.25, -0.05, -0.1, 0.1),
                 to = c(-0.1, -0.2, 0.1, -0.1))+
    scale_x_break(c(-0.62,0.00),scales = 3,expand = T)+
    theme(panel.background = element_rect(fill='white', colour='black'),
          panel.grid=element_blank(), 
          legend.position = "right",
          axis.title = element_text(color='black',size=14),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.text=element_text(colour='black',size=11))
  
  return(p)
}

# rhizosphere
group2 <- subset(group, Compartment == "R" & (Type %in% c("HH", "HM", "HL", "MM", "ML","CKH","CKM")))
p4a_2 <- perform_pcoa_CK(group2, "Rhizosphere")

print(p4a_1)
print(p4a_2)
p4a_1+p4a_2

####Fig4b
####soil
soil <- read.xlsx("field_translocation.xlsx", sheet = "source_sink_soil", colNames = TRUE, rowNames = FALSE)
soil$type <- paste0(soil$Soil,soil$Site)
soil1 <- subset(soil, Species != "F")
df1 <- data_summary(soil1, varname="Proportion",
                    groupnames=c("type", "Group"))

types <- c("HM","HL","ML")
plots <- list()
for (i in types) {
  df2 <- subset(df1, type == i)
  p <- ggplot() +
    geom_arc_bar(data = df2,
                 stat = "pie",
                 aes(x0 = 0, y0 = 0, r0 = 1.1, r = 3,
                     amount = Proportion, fill = Group)) +
    scale_fill_manual(values = c(local_leaf="#BDB985", local_soil="#F8830C", origin_soil="#C68267", Unknown="#D3D1D1")) +
    theme_void() +
    theme(legend.position = "none") +
    annotation_custom(grob = textGrob(i, 
                                      gp = gpar(fontsize = 20)), 
                      xmin = 0, xmax = 0, ymin = 0, ymax = 0) 
  
  plots[[i]] <- p
}
grid.arrange(grobs = plots, ncol = 3)

###leaf
leaf <- read.xlsx("field_translocation.xlsx", sheet = "source_sink_leaf", colNames = TRUE, rowNames = FALSE)
leaf$type <- paste0(leaf$Soil,leaf$Site)
leaf1 <- subset(leaf, Species != "F")
df1 <- data_summary(leaf1, varname="Proportion",
                    groupnames=c("type", "Group"))

types <- c("HM","HL","ML")
plots <- list()
for (i in types) {
  df2 <- subset(df1, type == i)
  p <- ggplot() +
    geom_arc_bar(data = df2,
                 stat = "pie",
                 aes(x0 = 0, y0 = 0, r0 = 1.1, r = 3,
                     amount = Proportion, fill = Group)) +
    scale_fill_manual(values = c(local_leaf="#BDB985", local_soil="#F8830C", origin_soil="#C68267", Unknown="#D3D1D1")) +
    theme_void() +
    theme(legend.position = "none") +
    annotation_custom(grob = textGrob(i, 
                                      gp = gpar(fontsize = 20)), 
                      xmin = 0, xmax = 0, ymin = 0, ymax = 0)
  
  plots[[i]] <- p
}
grid.arrange(grobs = plots, ncol = 3)

####Fig4c
group <- read.xlsx("field_translocation.xlsx", sheet = "group",colNames = T, rowNames = T)
bacteria_unifrac1 <- readRDS("trans_bray.rds")
bacteria_unifrac2 <- data.frame(as.matrix(bacteria_unifrac1))

leaf <- subset(group, Compartment == "L")
soil <- subset(group, Compartment == "R")
data <- bacteria_unifrac2[rownames(bacteria_unifrac2) %in% rownames(leaf),]
rhizos1 <- data[, colnames(data) %in% rownames(soil)]
diagonal_values <- diag(as.matrix(rhizos1))  
row_names_diag <- rownames(rhizos1)  
col_names_diag <- colnames(rhizos1)  
result <- data.frame(ID = row_names_diag, id = col_names_diag, Value = diagonal_values)
results <- left_join(result, leaf[c(1,3,5,6,8)],by = "ID")

####entire
results_entire <- results %>%
  mutate(elevation = case_when(
    Type == "HL" ~ 1051,
    Type == "HM" ~ 695,
    Type == "ML" ~ 356,
    TRUE ~ 0))
#write.table(results_entire, file ="results_entire.csv",sep =",", quote =FALSE)

fit<-gamm(formula(Value~s(elevation,bs='cr', k = 3)), 
          random=list(Species=~1),data=results_entire,
          family=gaussian,method = 'REML')
draw(fit$gam)+theme_bw()
summary(fit$gam)

p4c <- ggplot() +
  geom_point(results_entire, mapping = aes(x = elevation, y = Value, color = Species),size = 3,shape = 21,alpha = 0.5) +
  stat_smooth(data = results_entire,
              aes(x = elevation, y = Value),
              method = "gam", formula = y ~ s(x, k = 3),se = TRUE, linewidth = 1, color = "black") +
  theme_bw()+
  scale_color_manual(values = c("#6D6D67","#F6B16B","#B84A23")) +
  labs(x = "Transplantation elevational distance (m)",y = "Community dissimilarity", title = NULL) +
  theme_classic()+
  theme(axis.title = element_text(color='black',size=12),
        axis.title.x=element_text(colour='black', size=12,vjust = 0),
        axis.title.y=element_text(colour='black', size=12,vjust = 0),
        axis.text=element_text(colour='black',size=10),
        legend.title=element_text(size = 12,face = "bold"),
        legend.text=element_text(size= 10),
        legend.key=element_blank(),
        legend.position = "right",
        legend.background = element_rect(colour = "black"));p4c

#########################Fig5#################
group <- read.xlsx("field_translocation.xlsx", sheet = "plant",colNames = T, rowNames = T)
group1 <- subset(group, Elevational != "0")
fit<-gamm(formula(lnHome_CK~s(Elevational,bs='cr', k = 3)), 
          random=list(Species=~1),data=group,
          family=gaussian,method = 'REML')
draw(fit$gam)+theme_bw()
summary(fit$gam)

df <- data_summary(group1, varname = "lnHome_CK", groupnames = c("Elevational", "Species"))
df$Elevational <- factor(df$Elevational,levels = c("356","695","1051"))
df$CI <- 1.96*df$se

ggplot(data = df) +
  geom_errorbar(mapping = aes(x = Elevational, y = lnHome_CK, ymin = lnHome_CK - se, ymax = lnHome_CK + se, 
                              color = Species, group = Species),
                position = position_dodge(width = 0.8), width = 0, size = 0.5, alpha = 1) +
  geom_bar(mapping = aes(x = Elevational, y = lnHome_CK, fill = Species, group = Species, color = Species),
           position = position_dodge(width = 0.8), stat = "identity", size = 0.7, width = 0.7, alpha = 0.8) +
  scale_color_manual(values = c("#6D6D67","#F6B16B","#B84A23")) + 
  scale_fill_manual(values = c("#6D6D67","#F6B16B","#B84A23")) + 
  labs(y = "Effect size", title = NULL) +
  theme_bw() +
  geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.5) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.line = element_line(linetype = 1, color = "black", size = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5, lineend = "round"),
        axis.text = element_text(colour = 'black', size = 9),
        axis.title = element_text(colour = 'black', size = 10)) -> Fig1a;Fig1a

### Compare with 0
group1 <- subset(group, Elevational != "0")
elevation <- unique(group1$Elevational)
spe <- unique(group1$Species)
result <- data.frame(elevation = numeric(), species = character(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through the altitudes and species
for (i in elevation) {
  for (j in spe) {
    
    data <- subset(group1, Elevational == i & Species == j)
    
    # Perform a Wilcoxon test. Note that the test here should be based on the lnHome_CK column in data.
    if (nrow(data) > 0) {  
      result1 <- wilcox.test(data$lnHome_CK, mu = 0, alternative = "two.sided")
      
      
      result <- rbind(result, data.frame(elevation = i, species = j, p_value = result1$p.value))
    } else {
      result <- rbind(result, data.frame(elevation = i, species = j, p_value = NA))
    }
  }
}
result

###corrlation with lnHome_CK & distance
fit <- lmer(lnHome_CK~distance+(1|Species)+(1|Popu),group1)
summary(fit)
r.squaredGLMM(fit)

get.data = get_model_data(fit, type = "pred", terms = "distance [all]")
get.data = as.data.frame(get.data); colnames(get.data)[1] = "distance"
##### rotation

ggplot(get.data, aes(x = distance, y = predicted)) +
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high),fill="#6D6D67", alpha = I(0.2)) +
  geom_line(size=1, linetype = 1, color = "#6D6D67") +
  geom_point(data = group1, 
             aes(x = distance, y = lnHome_CK, color = Species), width = 0.2,
             height = 0, size = 2.5, shape = 21, 
             stroke = 0.5, alpha = 0.9,
             show.legend = FALSE) +
  theme_bw()+
  scale_color_manual(values = c("#6D6D67","#F6B16B","#B84A23"))+
  labs(x = "Community dissimilarity\n(Pair phyllosphere and rhizosphere)",y = "ln(PGR in transplanted site/PGR in origin site)") +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.position = "none",
        axis.line=element_line(linetype=1,color="black",size=0.5),
        axis.ticks=element_line(color="black",size=0.5,lineend = 2),
        axis.text=element_text(colour='black',size=9),
        axis.title=element_text(colour='black', size=10)) -> p1; p1

################FigS2#################
df1 <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)
df1$month <- factor(df1$month,levels = c("Sep","Oct","Dec","Jun","Aug"))

p1 <- ggplot(df1, aes(x = month, y = AT)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(color = "grey", size = 1, alpha = 0.5) +
  labs(x = "Month", y = "Air temperature (°C)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))

p2 <- ggplot(df1, aes(x = month, y = AM)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") + 
  geom_jitter(color = "grey", size = 1, alpha = 0.5) +  
  labs(x = "Month", y = "Air moisture (%)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))

p3 <- ggplot(df1, aes(x = month, y = AT)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
  geom_jitter(color = "grey", size = 1, alpha = 0.5) + 
  labs(x = "Month", y = "Specific leaf area (cm2/g)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))

p4 <- ggplot(df1, aes(x = month, y = AT)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") + 
  geom_jitter(color = "grey", size = 1, alpha = 0.5) + 
  labs(x = "Month", y = "Chlorophyll content (%)") +
  theme_bw() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5))

df <- read.xlsx("field_survey.xlsx", sheet = "height_june", colNames = TRUE, rowNames = TRUE)
p5 <- ggplot()+
  geom_point(df,mapping = aes(elevation, log(height),color = species), shape= 21)+
  labs(x = "Plant height (cm, log)", y = "Elevation (m)")+
  scale_color_manual(values = c("#974FA3", "#4CAB3B", "#3C72AF"))+
  stat_smooth(df,mapping = aes(elevation, log(height),color = species),
              method = "lm",
              formula = y ~ x,
              se = F,
              alpha = 0.5)+
  stat_poly_eq(df,mapping = aes(elevation, log(height),color = species, label =  paste(..rr.label..,..p.value.label..,  sep = "~~~~")),
               formula = y~x, parse = TRUE)+
  theme(panel.background = element_rect(fill='white', colour='black'),
        panel.grid=element_blank(), 
        axis.title.x=element_text(colour='black', size=12),
        axis.title.y=element_text(colour='black', size=12),
        axis.text=element_text(colour='black',size=10),
        legend.title=element_text(size = 10,face = "bold"),
        legend.text=element_text(size= 10),
        legend.position = "top");p5

combined_plot <- (p1 | p2) / (p3 | p4)

combined_plot <- combined_plot + 
  draw_plot_label(c("a", "b", "c", "d"), 
                  c(0.05, 0.05, 0.05, 0.05), 
                  c(0.95, 0.95, 0.45, 0.45), 
                  size = 5)

###############FigS5#############
group <- read.xlsx("field_translocation.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)
group$ASVs.of.total <- log10(group$ASVs.of.total)

create_plot <- function(data, y_var, y_label, fill_values, color_values, y_breaks = NULL) {
  p <- ggplot(data, aes(x = Compartment, group = Type)) +
    geom_errorbar(aes_string(ymin = paste(y_var, "- se"), ymax = paste(y_var, "+ se"), color = "Type"),
                  position = position_dodge(width = 0.8), width = 0.2, size = 0.5) +
    geom_bar(aes_string(y = y_var, fill = "Type"), 
             position = position_dodge(width = 0.8), stat = "identity", width = 0.7, size = 0.5) +
    labs(x = "Groups", y = y_label) +
    facet_wrap(~ Species, scales = "free")+
    scale_fill_manual(values = fill_values) +
    scale_color_manual(values = color_values) +
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid = element_blank(), 
          legend.position = "right",
          axis.title = element_text(color = 'black', size = 14),
          axis.ticks = element_line(color = 'black'),
          axis.line = element_line(colour = "black"),
          axis.text = element_text(colour = 'black', size = 11))
  
  if (!is.null(y_breaks)) {
    p <- p + scale_y_break(y_breaks[1:2], space = 0.001, scales = 0.5)
  }
  
  return(p)
}

# OTUs.of.total for high group
high <- subset(group, Type %in% c("HH", "HM", "HL"))
df1 <- data_summary(high, varname = "ASVs.of.total", groupnames = c("Compartment", "Type","Species"))
df1$Type <- factor(df1$Type, levels = c("HH", "HM", "HL"))

p1 <- create_plot(df1, "ASVs.of.total", "No. of total ASVs (log)", 
                  c("#6D9DFF", "#70D3ED", "#B79F00"), c("#6D9DFF", "#70D3ED", "#B79F00"))
p1

# OTUs.of.total for mid group
mid <- subset(group, Type %in% c("MM", "ML"))
df2 <- data_summary(mid, varname = "ASVs.of.total", groupnames = c("Compartment", "Type","Species"))
df2$Type <- factor(df2$Type, levels = c("MM", "ML"))
p2 <- create_plot(df2, "ASVs.of.total", "No. of total ASVs (log)", 
                  c("#E09513", "#9962EF"), c("#E09513", "#9962EF"))
p2

df1_shared <- data_summary(high, varname = "ASVs", groupnames = c("Compartment", "Type","Species"))
df1_shared$Type <- factor(df1_shared$Type, levels = c("HH", "HM", "HL"))
p3 <- create_plot(df1_shared, "ASVs", "% ASVs shared", 
                  c("#6D9DFF", "#70D3ED", "#B79F00"), c("#6D9DFF", "#70D3ED", "#B79F00"))
p3
# % OTUs shared for mid group
df2_shared <- data_summary(mid, varname = "ASVs", groupnames = c("Compartment", "Type","Species"))
df2_shared$Type <- factor(df2_shared$Type, levels = c("MM", "ML"))
p4 <- create_plot(df2_shared, "ASVs", "% ASVs shared", 
                  c("#E09513", "#9962EF"), c("#E09513", "#9962EF"))
p4
# Reads for high group
df1_reads <- data_summary(high, varname = "read.of.shared", groupnames = c("Compartment", "Type","Species"))
df1_reads$Type <- factor(df1_reads$Type, levels = c("HH", "HM", "HL"))
p5 <- create_plot(df1_reads, "read.of.shared", "RA of shared ASVs (%)", 
                  c("#6D9DFF", "#70D3ED", "#B79F00"), c("#6D9DFF", "#70D3ED", "#B79F00"))
p5
# Reads for mid group
df2_reads <- data_summary(mid, varname = "read.of.shared", groupnames = c("Compartment", "Type","Species"))
df2_reads$Type <- factor(df2_reads$Type, levels = c("MM", "ML"))
p6 <- create_plot(df2_reads, "read.of.shared", "RA of shared ASVs (%)", 
                  c("#E09513", "#9962EF"), c("#E09513", "#9962EF"))
p6
p1/p3/p5
p2/p4/p6

#####################FigS6##################
OTU <- readRDS("trans_unique_bray.rds")
group <- read.xlsx("field_translocation.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)
OTU1 <- data.frame(as.matrix(OTU))

perform_pcoa_plot <- function(group_data, title_suffix = "") {
  data <- OTU1[rownames(OTU1) %in% rownames(group_data),]
  rhizos1 <- data[, colnames(data) %in% rownames(group_data)]
  
  pcoa <- cmdscale(rhizos1, k = 2, eig = TRUE)
  poi <- as.data.frame(pcoa$points)
  colnames(poi) <- c("PCoA1", "PCoA2")
  poi$SampleID <- rownames(poi)
  
  sample_site <- cbind(poi, group_data[, c(3, 4, 5, 6, 8)])
  sample_site$Site <- factor(sample_site$Site)
  
  # Calculate mean and standard errors
  df1 <- data_summary(sample_site, varname = "PCoA1", groupnames = c("Type"))
  colnames(df1) <- c("Type", "mean_PCoA1", "sd_PCoA1", "se_PCoA1")
  
  df2 <- data_summary(sample_site, varname = "PCoA2", groupnames = c("Type"))
  colnames(df2) <- c("Type", "mean_PCoA2", "sd_PCoA2", "se_PCoA2")
  
  summary_df <- cbind(df1, df2)[, -5]
  
  summary_df$Type <- factor(summary_df$Type, levels = c("HH", "HM", "HL", "MM", "ML"))
  
  p <- ggplot() + 
    geom_point(sample_site, mapping = aes(x = PCoA1, y = PCoA2, fill = Type, shape = Species), size = 1.5, alpha = 0.5) +
    geom_errorbarh(summary_df, mapping = aes(x = mean_PCoA1, y = mean_PCoA2, 
                                             xmin = mean_PCoA1 + se_PCoA1, 
                                             xmax = mean_PCoA1 - se_PCoA1, color = Type), height = 0, size = 0.5) +
    geom_errorbar(summary_df, mapping = aes(x = mean_PCoA1, y = mean_PCoA2, 
                                            ymin = mean_PCoA2 + se_PCoA2, 
                                            ymax = mean_PCoA2 - se_PCoA2, color = Type), width = 0, size = 0.5) +
    geom_point(summary_df, mapping = aes(x = mean_PCoA1, y = mean_PCoA2, color = Type), size = 4) +
    scale_shape_manual(values = c(21, 22, 24)) +
    scale_fill_manual(values = c(HH = "#717171", HL = "#B79F00", HM = "#4BCAE8", LL = "#1D4899", ML = "#9962EF", MM = "#DF9413")) +
    scale_color_manual(values = c(HH = "#717171", HL = "#B79F00", HM = "#4BCAE8", LL = "#1D4899", ML = "#9962EF", MM = "#DF9413")) +
    geom_vline(aes(xintercept = 0), linetype = "dashed", size = 0.3) +
    geom_hline(aes(yintercept = 0), linetype = "dashed", size = 0.3) +
    scale_x_continuous(labels = scales::label_comma(accuracy = 0.01)) +
    scale_y_continuous(labels = scales::label_comma(accuracy = 0.01)) +
    labs(x = paste("PCoA1 (", format(100 * pcoa$eig[1] / sum(pcoa$eig), digits = 3), "%)", sep = ""), 
         y = paste("PCoA2 (", format(100 * pcoa$eig[2] / sum(pcoa$eig), digits = 3), "%)", sep = ""),
         title = paste("PCoA Plot", title_suffix)) +
    geom_magnify(from = c(-0.25, -0.05, -0.1, 0.1),
                 to = c(-0.1, -0.2, 0.1, -0.1))+
    theme(panel.background = element_rect(fill='white', colour='black'),
          panel.grid=element_blank(), 
          legend.position = "right",
          axis.title = element_text(color='black',size=14),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.text=element_text(colour='black',size=11))
  
  return(p)
}

# phyllosphere
group1 <- subset(group, Compartment == "L" & (Type %in% c("HH", "HM", "HL", "MM", "ML")))
p4c <- perform_pcoa_plot(group1, "Phyllosphere")

# rhizosphere
group2 <- subset(group, Compartment == "R" & (Type %in% c("HH", "HM", "HL", "MM", "ML")))
p4d <- perform_pcoa_plot(group2, "Rhizosphere")

print(p4a)
print(p4b)
print(p4c)
print(p4d)
(p4a+p4b)/(p4c+p4d)

