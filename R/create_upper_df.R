#
#
#
library(data.table)
library(dplyr)
library(stringr)
# create df_mutation_count matrix


df_mutation = read.csv("C:/BBhack/Lung/Lung_Mutation.txt", sep = "\t") %>% as.data.table


# create unique GENE_CHR combination 
df_mutation$Gene_Chr = str_c(df_mutation$Gene,"_", df_mutation$Chr)

# select unique patients and gene_chr's in the df
patients_mutation = unique(df_mutation$Sample_ID)
gene_chr_mutation = unique(df_mutation$Gene_Chr)

# storage
df_mutation_count = matrix(0, nrow = length(patients_mutation), ncol = length(gene_chr_mutation))
rownames(df_mutation_count) = patients_mutation
colnames(df_mutation_count) = str_c(gene_chr_mutation, "_mut")


# fill by row (for each patient)
for(patient in patients_mutation){
  
  patient_id = match(patient, patients_mutation)
  
  gene_chr_count = 
    df_mutation %>% 
    filter(Sample_ID == patient) %>% 
    group_by(Gene_Chr) %>% 
    summarize(count = n()) %>% 
    as.data.table
    
    
  gene_chr_id = match(gene_chr_count %>% pull(Gene_Chr), gene_chr_mutation)
  
  df_mutation_count[patient_id, gene_chr_id] = (gene_chr_count %>% pull(count))
  
  
}

df_mutation_count = df_mutation_count %>% as.data.table
rownames(df_mutation_count) = patients_mutation


dim(df_mutation)
sum(df_mutation_count)


#  CNV 
df_CNV = read.csv("C:/BBhack/Lung/Lung_CNV.txt", sep = "\t") %>% as.data.table
df_CNV = df_CNV %>% select(-c(Start, Stop, Strand)) %>% as.data.table
colnames(df_CNV) = gsub('\\.', '-', colnames(df_CNV))
df_CNV$Gene_Chr = str_c(df_CNV$Gene, df_CNV$Chr)

# remove duplicates
df_CNV = df_CNV[!duplicated(df_CNV$Gene_Chr, fromLast = TRUE), ] 

# create variable df for patient
df_CNV_count = df_CNV %>% select(-c(Gene, Chr, Gene_Chr)) %>% t() %>% as.data.table()
colnames(df_CNV_count) = str_c(df_CNV$Gene_Chr)
rownames(df_CNV_count) = df_CNV %>% select(-c(Gene, Chr, Gene_Chr)) %>% colnames()

#
df_CNV_count = add_rownames(df_CNV_count) %>% as.data.table
df_mutation_count = add_rownames(df_mutation_count) %>% as.data.table




# create joint dataframe 
df_upper = inner_join(df_CNV_count, df_mutation_count, by = "rowname")

sum(is.na(df_upper))


rm(df_CNV, df_CNV_count, df_mutation, df_mutation_count)


write.csv(df_upper, "")





  


