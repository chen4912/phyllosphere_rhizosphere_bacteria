#####################################################################
################Screening of shared microorganisms###################
#####################################################################
Packages <- c("openxlsx", "edgeR","tidyverse")
field_otu <- read.xlsx("field_survey.xlsx", sheet = "ASV",colNames = T, rowNames = T)#加载行为OTU，列为样本的微生物数据
merged_data <- data.frame(rownames(field_otu))
#df <- data.frame(matrix(ncol = 488, nrow = 0))

# Create new column names
new_column_names <- c()
for (i in 1:244) {
  new_column_names <- c(new_column_names, paste0("L", i), paste0("S", i))
}

# Make sure the number of new column names matches the original column names
if (length(new_column_names) == ncol(field_otu)) {
  colnames(field_otu) <- new_column_names
} else {
  print("The number of new column names does not match the number of original column names")
}

for (i in 1:244) {
  col_A <- paste0("L", i)
  col_B <- paste0("S", i)
  merged_data <- cbind(merged_data, 
                       ifelse(field_otu[[col_A]] > 0 & field_otu[[col_B]] > 0, 
                              field_otu[[col_A]], 
                              0),
                       ifelse(field_otu[[col_A]] > 0 & field_otu[[col_B]] > 0, 
                              field_otu[[col_B]], 
                              0))
  colnames(merged_data)[ncol(merged_data) - 1] <- col_A
  colnames(merged_data)[ncol(merged_data)] <- col_B
}
head(merged_data)
rownames(merged_data) <- merged_data[, 1]
microbial_data <- merged_data[, -1]
write.table(microbial_data, file ="microbial_data.csv",sep =",", quote =FALSE) #结果导出

###unique
total_group <- read.xlsx("field_survey.xlsx", sheet = "ASV",colNames = T, rowNames = T)
field_group <- read.xlsx("field_survey.xlsx", sheet = "shared_ASV",colNames = T, rowNames = T)
field_group2 <- total_group - field_group
field_group1 <- field_group[rowSums(field_group != 0) > 0, ]
field_group2 <- field_group2[rowSums(field_group2 != 0) > 0, ]
#write.table(field_group1, file ="SHARED.csv",sep =",", quote =FALSE)
#write.table(field_group2, file ="UNIQUE.csv",sep =",", quote =FALSE)
