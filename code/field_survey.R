Packages <- c(
  "ape", "lme4", "ggtern", "nlme", "lmerTest", "openxlsx", 
  "edgeR", "ggplot2", "vegan", "tidyverse", "ggalt", 
  "patchwork", "NST", "ggpmisc", "reshape2", "ggprism", 
  "ggpubr", "ggalluvial", "dplyr", "ggh4x", "car", "NST",
  "sjPlot", "performance", "MuMIn", "mgcv", "emmeans", 
  "multcomp", "stringr", "gamm4", "Matrix", "ggthemes", 
  "tidyr", "glmm.hp", "partR2", "corrplot", "gratia", 
  "ggstar", "ggmagnify", "ggbreak", "picante", "RColorBrewer", 
  "ggforce", "grid", "gridExtra", "linkET", "cowplot","ggalluvial","kableExtra"
)
#install.packages(setdiff(Packages, rownames(installed.packages())))
lapply(Packages, library, character.only = TRUE)

setwd("D:/desktop/ASV/")

###########Common functions################
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
##########################Tables######################
###alpha diversity in field survey
set.seed(12345)
diversity <- read.xlsx("field_survey.xlsx", sheet = "group",colNames = T, rowNames = T)
#diversity$plot <- paste0(diversity$species,"_",diversity$elevation)
diversity$month <- factor(diversity$month,levels = c("Sep","Oct","Dec","Jun","Aug"))

####model of time ////elevation
model <- lmer(log(ASVs.of.total) ~ species * compartment * month + (1|pair)+(1|elevation), 
              na.action= na.omit, data = diversity)
car::Anova(model, type = 3)
ranova(model)

model <- lmer(ASVs ~ species * compartment * month + (1|pair)+(1|elevation), 
              na.action= na.omit, data = diversity)
car::Anova(model, type = 3)
ranova(model)

model <- lmer(read.of.shared ~ species * compartment * month + (1|pair)+(1|elevation), 
              na.action= na.omit, data = diversity)
car::Anova(model, type = 3)
ranova(model)

###beta diversity in field survey
###PERMANOVA ANALYSIS
total <- readRDS("field_bray.rds")
total <- data.frame(as.matrix(total))
group <- read.xlsx("field_survey.xlsx", sheet = "group",colNames = T, rowNames = T)

with(group, adonis2(total ~ pair + month * species * compartment, data = group, permutations = 999, strata = elevation, by = "term"))

with(group, adonis2(total ~ pair+species+compartment+pH+AT+AM+STN+LTN+SLA+Chl+
                      compartment:pH+compartment:AT+compartment:AM+compartment:STN+compartment:LTN+
                      compartment:SLA+compartment:Chl,
                    data = group, permutations = 999, strata = month, by = "term"))
##leaf & environment
group1 <- subset(group, compartment == "L")
data <- total[rownames(total) %in% rownames(group1),]
bray <- data[,colnames(data) %in% rownames(group1),]
with(group1, adonis2(bray ~ species+pH+AT+AM+STN+LTN+SLA+Chl+
                       species:pH+species:AT+species:AM+species:STN+species:LTN+species:SLA+
                       species:Chl,data = group1, permutations = 999,strata = elevation, by = "term"))

##root & environment
group1 <- subset(group, compartment == "R")
data <- total[rownames(total) %in% rownames(group1),]
bray <- data[,colnames(data) %in% rownames(group1),]
with(group1, adonis2(bray ~ species+pH+AT+AM+STN+LTN+SLA+Chl+
                       species:pH+species:AT+species:AM+species:STN+species:LTN+species:SLA+
                       species:Chl,data = group1, permutations = 999,strata = elevation, by = "term"))

###shared & unique analysis
total <- readRDS("field_shared_bray.rds")
total <- data.frame(as.matrix(total))
group <- read.xlsx("field_survey.xlsx", sheet = "group",colNames = T, rowNames = T)
#group$plot <- paste0(group$species,"_",group$elevation)

with(group, adonis2(total ~ pair+species+compartment+pH+AT+AM+STN+LTN+SLA+Chl+
                      compartment:pH+compartment:AT+compartment:AM+compartment:STN+compartment:LTN+
                      compartment:SLA+compartment:Chl,
                    data = group, permutations = 999, strata = elevation, by = "term"))

total <- readRDS("field_unique_bray.rds")
total <- data.frame(as.matrix(total))
group <- read.xlsx("field_survey.xlsx", sheet = "group",colNames = T, rowNames = T)
#group$plot <- paste0(group$species,"_",group$elevation)

with(group, adonis2(total ~ pair+species+compartment+pH+AT+AM+STN+LTN+SLA+Chl+
                      compartment:pH+compartment:AT+compartment:AM+compartment:STN+compartment:LTN+
                      compartment:SLA+compartment:Chl,
                    data = group, permutations = 999, strata = elevation, by = "term"))

############################################################################
#################################graph######################################
############################################################################
#######################Fig1####################
diversity <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)

diversity$logASVs.of.total <- log(diversity$ASVs.of.total)
df <- data_summary(diversity, varname = "logASVs.of.total", groupnames = c("compartment"))

create_plot <- function(data, response_var, y_label) {
  df <- data_summary(data, varname = response_var, groupnames = c("compartment", "month", "species"))
  df$month <- factor(df$month, levels = c("Sep", "Oct", "Dec", "Jun", "Aug"))
  
  ggplot() +
    geom_line(data = df, mapping = aes_string(x = "month", y = response_var, group = "compartment", color = "compartment"), size = 0.3) +
    geom_jitter(data = data, mapping = aes_string(x = "month", y = response_var, shape = "compartment", color = "compartment"), size = 0.7, width = 0.03, alpha = 0.3) +
    geom_point(data = df, mapping = aes_string(x = "month", y = response_var, color = "compartment", shape = "compartment"), size = 0.5) +
    geom_errorbar(data = df, mapping = aes_string(x = "month", y = response_var, ymin = paste0(response_var, " - se"), ymax = paste0(response_var, " + se")), color = "black", width = 0.1, size = 0.5) +
    scale_shape_manual(values = c("L" = 16, "R" = 16)) +
    scale_color_manual(values = c("L" = "#108b96", "R" = "#7A8234")) +
    labs(y = y_label) +
    theme_classic() +
    theme(
      axis.title = element_text(color = 'black', size = 12),
      axis.title.x = element_text(colour = 'black', size = 12, vjust = 0),
      axis.title.y = element_text(colour = 'black', size = 12, vjust = 0),
      axis.text = element_text(colour = 'black', size = 10),
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 10),
      legend.key = element_blank(),
      legend.position = "none",
      legend.background = element_rect(colour = "black")
    ) +
    facet_wrap(~species) +
    theme(strip.text = element_blank())
}

p1a <- create_plot(diversity, "logASVs.of.total", "No. of total ASVs (log)")
p1b <- create_plot(diversity, "ASVs", "% ASVs shared")
p1c <- create_plot(diversity, "read.of.shared", "RA of shared ASVs (%)")

layout <- "
AAAAAA
BBBBBB
CCCCCC
"
p1a + p1b + p1c + 
  plot_layout(design = layout)

###Multiple comparisons (based on mixed models)
set.seed(12345)
group <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)
group$ASVs.of.total <- log10(group$ASVs.of.total)
group$month <- factor(group$month,levels = c("Sep","Oct","Dec","Jun","Aug"))
group$plot <- paste0(group$species,"_",group$elevation)
##ASVs.of.total    |    ASVs     |    read.of.shared
###emmeans packages
model <- lmer(log(ASVs.of.total) ~ species * compartment * month + (1|pair)+(1|plot), 
              na.action= na.omit, data = group)
car::Anova(model, type = 3)
ranova(model)

emm1 = emmeans(model, specs = pairwise ~ compartment * month|species, type = 'response', adjust = 'tukey')
emm1_multi = multcomp::cld(emm1,alpha=0.05,Letters=letters,adjust="none",decreasing = T)
emm1_multi$.group <- trimws(emm1_multi$.group)
emm1_multi

################Fig2##############
###Fig2a
entire_dis <- readRDS("field_bray.rds")
group <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)
entire_dis <- as.matrix(entire_dis)

perform_pcoa_and_plot <- function(data, group_data, compartment1) {
  compartment_data <- subset(group_data, compartment == compartment1)
  filtered_data <- data[rownames(data) %in% rownames(compartment_data), ]
  leaf_data <- filtered_data[, colnames(filtered_data) %in% rownames(compartment_data)]
  
  pcoa_result <- cmdscale(leaf_data, k = 2, eig = TRUE)
  points_df <- as.data.frame(pcoa_result$points)
  names(points_df)[1:2] <- c("PCoA1", "PCoA2")
  points_df$SampleID <- rownames(points_df)
  
  group_data$SampleID <- rownames(group_data)
  sample_site <- merge(points_df, group_data)
  
  df <- aggregate(cbind(PCoA1, PCoA2) ~ month, sample_site, FUN = mean)
  
  colnames(df)[2:3] <- c("PCoA1_mean", "PCoA2_mean")
  
  sample_site$type <- paste0(sample_site$month)
  df$type <- paste0(df$month)
  sample_site2 <- sample_site %>% left_join(df[, c(2:4)], by = "type")
  
  df_line <- aggregate(cbind(PCoA1, PCoA2) ~ month, sample_site, FUN = mean)
  colnames(df_line)[2:3] <- c("PCoA1_line", "PCoA2_line")
  
  sample_site2$month <- factor(sample_site2$month, levels = c("Sep", "Oct", "Dec", "Jun", "Aug"))
  df$month <- factor(df$month, levels = c("Sep", "Oct", "Dec", "Jun", "Aug"))
  df_line$month <- factor(df_line$month, levels = c("Sep", "Oct", "Dec", "Jun", "Aug"))
  df_line <- df_line %>% arrange(month)

  pcoa_eig <- round(pcoa_result$eig[1:2] / sum(pcoa_result$eig), 3)
  
  ggplot() +
    geom_point(sample_site2, mapping = aes(PCoA1, PCoA2, color = month), shape = 21, size = 1.5, alpha = 0.8) +
    geom_point(df, mapping = aes(PCoA1_mean, PCoA2_mean, fill = month), size = 2.5, alpha = 1, shape = 21) +
    geom_point(df_line, mapping = aes(PCoA1_line, PCoA2_line), color = NA, size = 2.5) +
    scale_color_manual(values = c("#EEE90A", "#32B40E", "#7E549A", "#3B528B", "#21908C")) +
    scale_fill_manual(values = c("#EEE90A", "#32B40E", "#7E549A", "#3B528B", "#21908C")) +
    geom_encircle(data = sample_site2, mapping = aes(PCoA1, PCoA2, group = month, fill = month), expand = 0, spread = 0.5,
                  s_shape = 1, size = 1, linetype = 1, alpha = 0.1) +
    geom_curve(data = df_line, mapping = aes(x = lag(PCoA1_line), y = lag(PCoA2_line), 
                                             xend = PCoA1_line, yend = PCoA2_line),
               curvature = 0.3, arrow = arrow(length = unit(0.3, "cm")), size = 1) +
    labs(x = paste("PCoA1 (", format(100 * pcoa_eig[1], digits = 3), "%)", sep = ""),
         y = paste("PCoA2 (", format(100 * pcoa_eig[2], digits = 3), "%)", sep = "")) +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 12),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid = element_blank(), legend.position = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
    scale_x_continuous(labels = scales::label_comma(accuracy = 0.001)) +
    scale_y_continuous(labels = scales::label_comma(accuracy = 0.001))
}
p2a <- perform_pcoa_and_plot(entire_dis, group, "L")
p2b <- perform_pcoa_and_plot(entire_dis, group, "R")
p2a + p2b

####p2c
###F. nilgerrensis // P.asiatica // T. repens
entire_dis <- readRDS("field_bray.rds")
data <- as.matrix(entire_dis)
merge_and_filter <- function(data, group, columns) {
  group1 <- group[, columns]
  group1$group <- rownames(group)
  
  data_3col <- dist.3col(data) %>%
    setNames(c("group", "group2", "distance")) %>%
    left_join(group1, by = "group")
  
  group2 <- group[, columns]
  group2$group2 <- rownames(group)
  
  data_joined <- left_join(data_3col, group2, by = "group2")
  data_joined$species <- ifelse(data_joined$species.x == data_joined$species.y, data_joined$species.x, 0)
  data_filtered <- subset(data_joined, species != 0)
  
  return(data_filtered)
}

columns <- c(2, 3, 5, 6, 7)
data_filtered <- merge_and_filter(data, group, columns)

data_filtered$elevation <- ifelse(data_filtered$elevation.x == data_filtered$elevation.y, data_filtered$elevation.x, 0)
data_filtered <- subset(data_filtered, elevation != 0)

data_filtered$compartment <- ifelse(data_filtered$compartment.x == data_filtered$compartment.y, data_filtered$compartment.x, 0)
data_filtered <- subset(data_filtered, compartment != 0)

data_filtered$time <- abs(data_filtered$time.x - data_filtered$time.y)

leaf <- subset(data_filtered, compartment == "L" & species == "F")##//P or T
root <- subset(data_filtered, compartment == "R" & species == "F")##//P or T
###gamm model
fit_leaf <- gamm(distance ~ s(time, bs = 'cr', k = 3), random = list(elevation.x = ~1), data = leaf, family = gaussian, method = "REML")
fit_root <- gamm(distance ~ s(time, bs = 'cr', k = 3), random = list(elevation.x = ~1), data = root, family = gaussian, method = "REML")

draw(fit_leaf$gam) + theme_bw()
summary(fit_leaf$gam)

draw(fit_root$gam) + theme_bw()
summary(fit_root$gam)

###p1// p2 // p3
p1 <- ggplot() +
  geom_jitter(data = root, 
             aes(x = time, y = distance), width = 0.2,
            height = 0, size = 1, shape = 21, 
           stroke = 0.3, color = "#7A8234", alpha = 0.7,
          show.legend = FALSE) +
  geom_jitter(data = leaf, 
             aes(x = time, y = distance), width = 0.2,
            height = 0, size = 1, shape = 21, 
           stroke = 0.3, color = "#108b96", alpha = 0.7,
          show.legend = FALSE) +
  stat_smooth(data = leaf,
              aes(x = time, y = distance),
              method = "gam", formula = y ~ s(x, k = 3),se = TRUE, linewidth = 0.5, color = "#108b96", fill = "grey70") +
  stat_smooth(data = root,
              aes(x = time, y = distance),
              method = "gam", formula = y ~ s(x, k = 3),se = TRUE, linewidth = 0.5, color = "#7A8234", fill = "grey70") +
  labs(x = "Time interval (month)", y = "Community dissimilarity \n (weighted Unifrac distance)") +
  theme_bw() +
  scale_y_continuous(labels = scales::label_comma(accuracy = 0.001)) +
  theme(panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 12));p1

###Fig2c
entire_dis <- readRDS("field_bray.rds")
entire_dis <- data.frame(as.matrix(entire_dis))
group <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)
leaf <- subset(group, compartment == "L")
soil <- subset(group, compartment == "R")
data <- entire_dis[rownames(entire_dis) %in% rownames(leaf),]
rhizos1 <- data[, colnames(data) %in% rownames(soil)]
diagonal_values <- diag(as.matrix(rhizos1))  
row_names_diag <- rownames(rhizos1) 
col_names_diag <- colnames(rhizos1)  
result <- data.frame(ID = row_names_diag, id = col_names_diag, Value = diagonal_values)
result1 <- left_join(result, leaf[c(1,3,5,6)],by = "ID")
Fdata <- subset(result1, species == "T") ##F/P/T
fit<-gamm(formula(Value~s(time,bs='cr', k = 5)), 
          random=list(elevation=~1),data=Fdata,
          family=gaussian)
summary(fit$gam)
draw(fit$gam)+theme_bw()
##p1/p2/p3
Fig2c <- ggplot() +
  geom_jitter(data = Fdata, 
              aes(x = time, y = Value), width = 0.2,
              height = 0, size = 1, shape = 21, 
              stroke = 0.3, fill = "black", alpha = 0.3,
              show.legend = FALSE) +
  stat_smooth(data = Fdata,
              aes(x = time, y = Value),
              method = "gam", formula = y ~ s(x, k = 5),se = TRUE, linewidth = 0.5, color = "black", fill = "grey70") +
  theme_bw() +
  labs(y = "Community dissimilarity\n(Paired phyllosphere and rhizosphere)") +
  scale_y_continuous(labels = scales::label_comma(accuracy = 0.01)) +
  theme(panel.grid.major = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(color = "black", size = 10),
        panel.grid.minor = element_blank(),
        axis.title.y = element_text(color = "black", size = 12),
        axis.title.x = element_text(color = "black", size = 12));Fig2c

###Fig2d
##permanova analysis of three species
total <- readRDS("field_bray.rds")
total <- data.frame(as.matrix(total))
group <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)

compartments <- unique(group$compartment) 
species_list <- unique(group$species) 

results <- data.frame()

for (Compartment in compartments) {
  for (Species in species_list) {
    group_subset <- subset(group, compartment == Compartment & species == Species)
    if (nrow(group_subset) > 1) {
      data_subset <- total[rownames(total) %in% rownames(group_subset), ]
      bray <- data_subset[, colnames(data_subset) %in% rownames(group_subset)]
      if (ncol(bray) > 1) {
        adonis_result <- data.frame(with(group_subset, adonis2(bray ~ pH + AT + AM + STN + LTN + SLA + Chl,
                                                               data = group_subset, permutations = 999, 
                                                               strata = elevation, by = "term")))
        
        adonis_result$species = Species
        adonis_result$compartment = Compartment
        results <- rbind(results, adonis_result)
      }
    }
  }
}
results

###GRAPH
group <- read.xlsx("field_survey.xlsx", sheet = "PERMANOVA",colNames = T, rowNames = T)
summarised_results <- group %>%
  group_by(species, compartment, factor) %>%
  summarise(total_R2 = sum(R2, na.rm = TRUE), .groups = 'drop')
ggplot(group, aes(x = compartment, y = R2 *100, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +  
  labs(x = "", y = "Explained variance (%)", title = "Bacteria") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

group$alpha <- ifelse(group$p_value < 0.05, 0.6, 0.3)

group$variable <- factor(group$variable, levels=(c("AT","AM","pH","STN","LTN","SLA","Chl")))
Fig2d <- ggplot(group, aes(x = compartment, y = R2 *100,fill = compartment, alpha = alpha)) +
  geom_bar(group, mapping = aes(x = variable, y = R2 *100,fill = compartment),stat = "identity", position = "dodge") +
  labs(x = "Compartment", y = "Explained variance (%)", fill = "Group") +
  scale_fill_manual(values = c("L" = "#108b96", "R" = "#7A8234")) +
  theme_classic() +
  scale_alpha(range = c(0.2, 1)) +  
  #scale_fill_manual(values = c("#108b96",  "#9D302E")) +
  facet_wrap(~ species, scales = "free")+
  theme(axis.title = element_text(color = 'black', size = 12),
        axis.title.x = element_text(colour = 'black', size = 12, vjust = 0),
        axis.title.y = element_text(colour = 'black', size = 12, vjust = 0),
        axis.text.x = element_text(colour = 'black', size = 10,angle = 30),
        axis.text.y = element_text(colour = 'black', size = 10),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10),
        legend.key = element_blank(),
        legend.position = "none",
        legend.background = element_rect(colour = "black"));Fig2d

####################FigS3#################
field_taxonomy <- read.xlsx("field_survey.xlsx", sheet = "tax", colNames = T, rowNames = T)
field_taxonomy$otu = rownames(field_taxonomy)
colnames(field_taxonomy)

field_otu <- read.xlsx("field_survey.xlsx", sheet = "ASV1", colNames = T, rowNames = T)
#field_otu = field_otu[,-c(1,2)]
field_otu$otu = rownames(field_otu)
dim(field_otu)
field_otu = merge(field_otu, field_taxonomy, by = "otu", all.x = T) #%>% filter(!is.na(kindom))
dim(field_otu)
rownames(field_otu) = field_otu$otu
unique(field_otu$Kingdom)
field_otu$otu = NULL
field_otu[1:6,1:6]
dim(field_otu)
field_otu = field_otu[,c(1:488)]

field_group <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = T, rowNames = T)
field_group$sample = rownames(field_group)
field_group$position = factor(field_group$compartment, levels = c("L","R"))
head(field_group)
dim(field_group)

#Adding sampling week property
field_group$sampling_week[field_group$time == '9']  <- 1 
field_group$sampling_week[field_group$time == '10'] <- 2
field_group$sampling_week[field_group$time == '12'] <- 3
field_group$sampling_week[field_group$time == '18'] <- 4
field_group$sampling_week[field_group$time == '20'] <- 5
####Subsetting the data to the SPECIES and MONTHS
dim(field_otu)
field_otu[1:6,1:6]
field_otu=field_otu[rowSums(field_otu)>0,]
dim(field_otu)

rel_field_otu = decostand(field_otu,method= "total", MARGIN = 2)
colSums(rel_field_otu)

rel_field_otu[1:6,1:6]
re_otus_abun <- data.frame(otu = row.names(rel_field_otu), rel_field_otu) %>%
  gather(sample, abun, -otu)
head(re_otus_abun)
dim(re_otus_abun)
#library(dplyr)
otus_abun_tax_total = merge(re_otus_abun, field_taxonomy, by = "otu", all.x = T)
head(otus_abun_tax_total)
dim(otus_abun_tax_total)
Phylum_tax_total = otus_abun_tax_total %>% group_by(sample,Phylum) %>% summarise(abun = sum(abun))
head(Phylum_tax_total)

library(reshape2)
Phylum_data <- as.data.frame(t(acast(Phylum_tax_total, formula = sample ~ Phylum , value.var = "abun", fill = 0)))

### top leaf
Phylum_data_L = Phylum_data[,subset(field_group, position == "L")$sample]
rownames(Phylum_data_L)
leaf_top = as.data.frame(rowSums(Phylum_data_L)); colnames(leaf_top) = "Abun"; leaf_top$Phylum = rownames(leaf_top) 
top_10_leaf <- leaf_top[order(-leaf_top$Abun), ][1:10, ]

##
Phylum_data_R = Phylum_data[,subset(field_group, position == "R")$sample]
soil_top = as.data.frame(rowSums(Phylum_data_R)); colnames(soil_top) = "Abun"; soil_top$Phylum = rownames(soil_top) 
top_10_soil <- soil_top[order(-soil_top$Abun), ][1:10, ]

library(ggalluvial)
top.col <- c(Acidobacteriota="#2078B3", Actinobacteriota="#AEC6E8", Bacteroidota="#F98307", Chloroflexi="#FFBB79", Cyanobacteria="#2CA129",
             Firmicutes="#97DF88", Gemmatimonadota="#D52728", Myxococcota="#FE9995", Planctomycetota="#9466BC", Proteobacteria="#C2B0D7",
             Verrucomicrobiota = "#8C564B", Unknown="#C49C93")

re_otu_ggalluvial_leaf <- data.frame(Phylum = as.factor(row.names(Phylum_data)), Phylum_data) %>% 
  gather(sample, abun, -Phylum) %>%
  left_join(field_group, by = 'sample') %>%
  filter(position == 'L') %>%
  filter(Phylum %in% top_10_leaf$Phylum) %>%
  #left_join(P_clusters, by='Phylum') %>%
  #left_join(shared_field_taxonomy_genus, by='Phylum') %>%
  group_by(species,Phylum, sampling_week) %>%
  summarise(n_count=sum(abun),
            n_sample=length(unique(sample)),
            norm_ra=n_count/n_sample) 

re_otu_ggalluvial_leaf$Species = factor(re_otu_ggalluvial_leaf$species, 
                                        levels = c("F","P","T"), 
                                        labels = c("F. nilgerrensis", "P. asiatica", "T. repens"))

re_otu_ggalluvial_leaf$Phylum <- gsub("p__", "", re_otu_ggalluvial_leaf$Phylum)
re_otu_ggalluvial_leaf$Phylum <- trimws(re_otu_ggalluvial_leaf$Phylum)
re_otu_ggalluvial_leaf <- re_otu_ggalluvial_leaf %>%
  mutate(Phylum = ifelse(Phylum == "NA", "Others", Phylum))
unique(re_otu_ggalluvial_leaf$Phylum)

prel.ranked_leaf <- re_otu_ggalluvial_leaf %>% 
  ggplot(aes(sampling_week , norm_ra , alluvium = Phylum, fill = Phylum)) +
  geom_alluvium(decreasing = F, width = 0.1,alpha = 0.7, knot.pos = 0.1)+
  labs(title = NULL ,y = 'Ranked Relative Abundance', x = NULL, tag = NULL) +
  scale_fill_manual(values = top.col) +
  theme_bw()+
  theme(legend.position = "none",legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.y = element_text(color = "black", size = 13),
        strip.background = element_rect(fill = "#8C8884"),
        strip.text.x = element_text(color = "white", face = "italic", size = 11),
        legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous("", labels = c("1" = "Sep","2" = "Oct", "3" = "Dec","4" = "Jun","5" = "Aug")) + 
  facet_wrap(~Species); prel.ranked_leaf

### Root 
re_otu_ggalluvial_root <- data.frame(Phylum = as.factor(row.names(Phylum_data)), Phylum_data) %>% 
  gather(sample, abun, -Phylum) %>%
  left_join(field_group, by = 'sample') %>%
  filter(position == 'R') %>%
  filter(Phylum %in% top_10_soil$Phylum) %>%
  #left_join(P_clusters, by='Phylum') %>%
  #left_join(shared_field_taxonomy_genus, by='Phylum') %>%
  group_by(species,Phylum, sampling_week) %>%
  summarise(n_count=sum(abun),
            n_sample=length(unique(sample)),
            norm_ra=n_count/n_sample) 

re_otu_ggalluvial_root$Species = factor(re_otu_ggalluvial_root$species, 
                                        levels = c("F","P","T"), 
                                        labels = c("F. nilgerrensis", "P. asiatica", "T. repens"))

re_otu_ggalluvial_root$Phylum <- gsub("p__", "", re_otu_ggalluvial_root$Phylum)
re_otu_ggalluvial_root$Phylum <- trimws(re_otu_ggalluvial_root$Phylum)
unique(re_otu_ggalluvial_root$Phylum)

prel.ranked_root <- re_otu_ggalluvial_root %>% 
  ggplot(aes(sampling_week , norm_ra , alluvium = Phylum, fill = Phylum)) +
  geom_alluvium(decreasing = F, width = 0.1,alpha = 0.7, knot.pos = 0.1)+
  labs(title = NULL ,y = 'Ranked Relative Abundance', x = NULL, tag = NULL) +
  scale_fill_manual(values = top.col) +
  theme_bw()+
  annotate("segment", x = 1, xend = 3, y = 0.1, yend = 0.1) + 
  annotate("segment", x = 4, xend = 5, y = 0.1, yend = 0.1) + 
  theme(legend.position = "none",legend.title = element_blank(),
        axis.text.x = element_text(color = "black", size = 11),
        axis.text.y = element_text(color = "black", size = 11),
        axis.title.y = element_text(color = "black", size = 13),
        legend.key = element_rect(fill = 'transparent'), 
        legend.background = element_rect(fill = 'transparent'),
        strip.background = element_blank(),strip.text = element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_x_continuous("", labels = c("1" = "Sep","2" = "Oct", "3" = "Dec","4" = "Jun","5" = "Aug")) + 
  facet_wrap(~Species); prel.ranked_root

prel.ranked_leaf/prel.ranked_root

####################FigS4###############
#Take field survey as an example
total <- read.xlsx("field_survey.xlsx", sheet = "ASV1", colNames = TRUE, rowNames = TRUE)
field_group <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)

#Function: Perform EdgeR analysis
perform_edgeR_analysis <- function(species_name, total_data, field_group_data) {
  group1 <- subset(field_group_data, species == species_name)
  total_core_otu <- total_data[, colnames(total_data) %in% rownames(group1),]
  total_core_otu <- total_core_otu[rowSums(total_core_otu) != 0, ]
  
  group <- factor(group1$compartment, levels = c("L", "R"))
  library_sizes <- colSums(total_core_otu)
  total_core_otu <- total_core_otu[, library_sizes > 0]
  group <- group[library_sizes > 0]  
  dgelist <- DGEList(counts = total_core_otu, group = group)
  
  dgelist <- DGEList(counts = total_core_otu, group = group)
  keep <- rowSums(cpm(dgelist) > 1) >= 2
  dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
  dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
  
  design <- model.matrix(~group)
  dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
  fit <- glmFit(dge, design, robust = TRUE)
  lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
  
  edgeR_diff <- lrt$table[order(lrt$table$FDR, lrt$table$logFC, decreasing = c(FALSE, TRUE)), ]
  edgeR_diff[which(edgeR_diff$logFC >= 1 & edgeR_diff$FDR < 0.05), 'sig'] <- paste('Rich in root', sep = '')
  edgeR_diff[which(edgeR_diff$logFC <= -1 & edgeR_diff$FDR < 0.05), 'sig'] <- paste('Rich in leaf', sep = '')
  edgeR_diff[which(abs(edgeR_diff$logFC) <= 1 | edgeR_diff$FDR >= 0.05), 'sig'] <- 'Non-significant'
  
  return(edgeR_diff)
}

# Function: Create a volcano plot
create_volcano_plot <- function(edgeR_diff) {
  edgeR_diff$sig <- as.factor(edgeR_diff$sig)
  
  p <- ggplot(edgeR_diff, aes(logFC, -log(FDR, 10), color = sig)) +
    geom_point(size = 1.5, pch = 21) + 
    scale_color_manual(values = c('Rich in leaf' = "#108b96", 'Non-significant' ='gray40', 'Rich in root' = "#7A8234")) +
    theme_bw() +
    theme(panel.grid = element_blank(), panel.background = element_blank(), plot.title = element_text(size = 12), 
          axis.ticks.y = element_line(), axis.line.y = element_line(), axis.line.x = element_line(color = 'black'),
          axis.text.y = element_text(size = 10, angle = 0, colour = "black"),
          axis.text.x = element_text(size = 10, angle = 0, colour = "black"),
          axis.title.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          legend.position = "none",
          legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
    geom_vline(xintercept = c(-1, 1), color = 'gray', linewidth = 0.25, linetype = 2) + 
    geom_hline(yintercept = -log(0.05, 10), color = 'gray', linewidth = 0.25, linetype = 2) +
    labs(x = expression(Log[2]~(Fold~change)), y = expression(-Log[10]~(P-value)), color = NA)
  
  return(p)
}

# Perform analysis and create graphs
species_list <- c("F", "P", "T")
plots <- list()

for (species in species_list) {
  edgeR_diff <- perform_edgeR_analysis(species, total, field_group)
  plots[[species]] <- create_volcano_plot(edgeR_diff)
  
  filter_up <- subset(edgeR_diff, FDR < 0.05 & logFC >= 1)
  filter_down <- subset(edgeR_diff, FDR < 0.05 & logFC <= -1)
  print(paste('Rich in root for', species, ':', nrow(filter_up)))
  print(paste('Rich in leaf for', species, ':', nrow(filter_down)))
  print(paste('Total for', species, ':', nrow(edgeR_diff)))
}

p1 <- plots[["F"]]
p2 <- plots[["P"]]
p3 <- plots[["T"]]

(p1+p2+p3)


###
total <- read.xlsx("field_translocation.xlsx", sheet = "ASV1", colNames = TRUE, rowNames = TRUE)
field_group <- read.xlsx("field_translocation.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)
field_group <- subset(field_group, Time != "Jun")
total <- total[, colnames(total) %in% rownames(field_group),]

#Function: Perform EdgeR analysis
perform_edgeR_analysis <- function(species_name, total_data, field_group_data) {
  group1 <- subset(field_group_data, species == species_name)
  total_core_otu <- total_data[, colnames(total_data) %in% rownames(group1),]
  total_core_otu <- total_core_otu[rowSums(total_core_otu) != 0, ]
  
  group <- factor(group1$Compartment, levels = c("L", "R"))
  library_sizes <- colSums(total_core_otu)
  total_core_otu <- total_core_otu[, library_sizes > 0]
  group <- group[library_sizes > 0]  
  dgelist <- DGEList(counts = total_core_otu, group = group)
  
  dgelist <- DGEList(counts = total_core_otu, group = group)
  keep <- rowSums(cpm(dgelist) > 1) >= 2
  dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
  dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
  
  design <- model.matrix(~group)
  dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
  fit <- glmFit(dge, design, robust = TRUE)
  lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
  
  edgeR_diff <- lrt$table[order(lrt$table$FDR, lrt$table$logFC, decreasing = c(FALSE, TRUE)), ]
  edgeR_diff[which(edgeR_diff$logFC >= 1 & edgeR_diff$FDR < 0.05), 'sig'] <- paste('Rich in root', sep = '')
  edgeR_diff[which(edgeR_diff$logFC <= -1 & edgeR_diff$FDR < 0.05), 'sig'] <- paste('Rich in leaf', sep = '')
  edgeR_diff[which(abs(edgeR_diff$logFC) <= 1 | edgeR_diff$FDR >= 0.05), 'sig'] <- 'Non-significant'
  
  return(edgeR_diff)
}

# Function: Create a volcano plot
create_volcano_plot <- function(edgeR_diff) {
  edgeR_diff$sig <- as.factor(edgeR_diff$sig)
  
  p <- ggplot(edgeR_diff, aes(logFC, -log(FDR, 10), color = sig)) +
    geom_point(size = 1.5, pch = 21) + 
    scale_color_manual(values = c('Rich in leaf' = "#108b96", 'Non-significant' ='gray40', 'Rich in root' = "#7A8234")) +
    theme_bw() +
    theme(panel.grid = element_blank(), panel.background = element_blank(), plot.title = element_text(size = 12), 
          axis.ticks.y = element_line(), axis.line.y = element_line(), axis.line.x = element_line(color = 'black'),
          axis.text.y = element_text(size = 10, angle = 0, colour = "black"),
          axis.text.x = element_text(size = 10, angle = 0, colour = "black"),
          axis.title.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          legend.position = "none",
          legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent')) +
    geom_vline(xintercept = c(-1, 1), color = 'gray', linewidth = 0.25, linetype = 2) + 
    geom_hline(yintercept = -log(0.05, 10), color = 'gray', linewidth = 0.25, linetype = 2) +
    labs(x = expression(Log[2]~(Fold~change)), y = expression(-Log[10]~(P-value)), color = NA)
  
  return(p)
}

# Perform analysis and create graphs
species_list <- c("F", "P", "T")
plots <- list()

for (species in species_list) {
  edgeR_diff <- perform_edgeR_analysis(species, total, field_group)
  plots[[species]] <- create_volcano_plot(edgeR_diff)
  
  filter_up <- subset(edgeR_diff, FDR < 0.05 & logFC >= 1)
  filter_down <- subset(edgeR_diff, FDR < 0.05 & logFC <= -1)
  print(paste('Rich in root for', species, ':', nrow(filter_up)))
  print(paste('Rich in leaf for', species, ':', nrow(filter_down)))
  print(paste('Total for', species, ':', nrow(edgeR_diff)))
}

p4 <- plots[["F"]]
p5 <- plots[["P"]]
p6 <- plots[["T"]]

(p1+p2+p3)/(p4+p5+p6)
#############FigS5#############
entire_dis <- readRDS("field_unique_bray.rds")
group <- read.xlsx("field_survey.xlsx", sheet = "group", colNames = TRUE, rowNames = TRUE)
entire_dis <- as.matrix(entire_dis)
perform_pcoa_and_plot <- function(data, group_data, compartment1) {
  
  compartment_data <- subset(group_data, compartment == compartment1)
  filtered_data <- data[rownames(data) %in% rownames(compartment_data), ]
  leaf_data <- filtered_data[, colnames(filtered_data) %in% rownames(compartment_data)]
  
  pcoa_result <- cmdscale(leaf_data, k = 2, eig = TRUE)
  points_df <- as.data.frame(pcoa_result$points)
  names(points_df)[1:2] <- c("PCoA1", "PCoA2")
  points_df$SampleID <- rownames(points_df)
  
  
  group_data$SampleID <- rownames(group_data)
  sample_site <- merge(points_df, group_data)
  
  df <- aggregate(cbind(PCoA1, PCoA2) ~ month, sample_site, FUN = mean)
  
  colnames(df)[2:3] <- c("PCoA1_mean", "PCoA2_mean")
  
  sample_site$type <- paste0(sample_site$month)
  df$type <- paste0(df$month)
  sample_site2 <- sample_site %>% left_join(df[, c(2:4)], by = "type")
  
  df_line <- aggregate(cbind(PCoA1, PCoA2) ~ month, sample_site, FUN = mean)
  colnames(df_line)[2:3] <- c("PCoA1_line", "PCoA2_line")
  
  sample_site2$month <- factor(sample_site2$month, levels = c("Sep", "Oct", "Dec", "Jun", "Aug"))
  df$month <- factor(df$month, levels = c("Sep", "Oct", "Dec", "Jun", "Aug"))
  df_line$month <- factor(df_line$month, levels = c("Sep", "Oct", "Dec", "Jun", "Aug"))
  df_line <- df_line %>% arrange(month)

  pcoa_eig <- round(pcoa_result$eig[1:2] / sum(pcoa_result$eig), 3)
  
  ggplot() +
    geom_point(sample_site2, mapping = aes(PCoA1, PCoA2, color = month), shape = 21, size = 1.5, alpha = 0.8) +
    geom_point(df, mapping = aes(PCoA1_mean, PCoA2_mean, fill = month), size = 2.5, alpha = 1, shape = 21) +
    geom_point(df_line, mapping = aes(PCoA1_line, PCoA2_line), color = NA, size = 2.5) +
    scale_color_manual(values = c("#EEE90A", "#32B40E", "#7E549A", "#3B528B", "#21908C")) +
    scale_fill_manual(values = c("#EEE90A", "#32B40E", "#7E549A", "#3B528B", "#21908C")) +
    geom_encircle(data = sample_site2, mapping = aes(PCoA1, PCoA2, group = month, fill = month), expand = 0, spread = 0.5,
                  s_shape = 1, size = 1, linetype = 1, alpha = 0.1) +
    geom_curve(data = df_line, mapping = aes(x = lag(PCoA1_line), y = lag(PCoA2_line), 
                                             xend = PCoA1_line, yend = PCoA2_line),
               curvature = 0.3, arrow = arrow(length = unit(0.3, "cm")), size = 1) +
    labs(x = paste("PCoA1 (", format(100 * pcoa_eig[1], digits = 3), "%)", sep = ""),
         y = paste("PCoA2 (", format(100 * pcoa_eig[2], digits = 3), "%)", sep = "")) +
    theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 12),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid = element_blank(), legend.position = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
    scale_x_continuous(labels = scales::label_comma(accuracy = 0.0001)) +
    scale_y_continuous(labels = scales::label_comma(accuracy = 0.0001))
}
p2a <- perform_pcoa_and_plot(entire_dis, group, "L")
p2b <- perform_pcoa_and_plot(entire_dis, group, "R")
p2c <- perform_pcoa_and_plot(entire_dis, group, "L")
p2d <- perform_pcoa_and_plot(entire_dis, group, "R")
(p2a + p2b)/(p2c + p2d)

