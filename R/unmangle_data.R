# bring data in normal form
require(purrr)
require(fst)
require(tibble)
require(dplyr)
require(tidyr)
require(stringr)
# txt files to fst

files <- list.files('../../data/Lung', full.names = T, pattern = 'txt$')
files <- set_names(files, str_replace(basename(files), '.txt',''))
files_fst <- map(files,str_c,'.fst')
# walk2(files,files_fst,~fst::write_fst(read_tsv(..1),path=..2))

## Gene expression data

Lung_GeneExpression <- read_fst(files_fst$Lung_GeneExpression) %>% as_data_frame
Lung_GeneExpression_transposed <- Lung_GeneExpression %>% 
  gather(SampleID, value = gene_expression, starts_with("TCGA")) %>% select(Gene,Chr,gene_expression,SampleID) %>%
  tidyr::unite(gene_chr, Gene,Chr) %>% 
  distinct(SampleID, gene_chr,.keep_all=TRUE) %>% Lung_Mutation
  spread(gene_chr,gene_expression) 
write_fst(Lung_GeneExpression_transposed,'../../data/intermediates/Lung_GeneExpression_normal.fst')

## CNV data

Lung_CNV <- read_fst(files_fst$Lung_CNV) %>% as_data_frame
Lung_CNV_transposed <- Lung_CNV %>% 
  gather(SampleID, value = cnv, starts_with("TCGA")) %>% select(Gene,Chr,cnv,SampleID) %>%
  tidyr::unite(gene_chr, Gene,Chr) %>% 
  distinct(SampleID, gene_chr,.keep_all=TRUE) %>% 
  spread(gene_chr,cnv) 
write_fst(Lung_CNV_transposed,'../../data/intermediates/Lung_CNV_transposed.fst')

## mutation

Lung_Mutation <- read_fst(files_fst$Lung_Mutation) %>% as_data_frame
Lung_Mutation_transposed <- Lung_Mutation %>% rename(SampleID = Sample_ID) %>% 
  gather(SampleID, value = Effect, starts_with("TCGA")) %>% select(Gene,Chr,Effect,SampleID) %>%
  tidyr::unite(gene_chr, Gene,Chr) %>% 
  distinct(SampleID, gene_chr,.keep_all=TRUE) %>% 
  spread(gene_chr,Effect) 
write_fst(Lung_Mutation_transposed,'../../data/intermediates/Lung_Mutation_transposed.fst')

## mutation

Lung_Mutation <- read_fst(files_fst$Lung_Methylation) %>% as_data_frame
Lung_Mutation_transposed <- Lung_Mutation %>% rename(SampleID = Sample_ID) %>% 
  gather(SampleID, value = Effect, starts_with("TCGA")) %>% select(Gene,Chr,Effect,SampleID) %>%
  tidyr::unite(gene_chr, Gene,Chr) %>% 
  distinct(SampleID, gene_chr,.keep_all=TRUE) %>% 
  spread(gene_chr,Effect) 
write_fst(Lung_Mutation_transposed,'../../data/intermediates/Lung_Mutation_transposed.fst')



