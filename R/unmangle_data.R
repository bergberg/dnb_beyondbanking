# bring data in normal form
require(purrr)
require(fst)
require(stringr)
# txt files to fst

files <- list.files('Z:/startset/Lung', full.names = T, pattern = 'txt$')
files <- set_names(files, str_replace(basename(files), '.txt',''))
files_fst <- map(files,str_c,'.fst')
# walk2(files,files_fst,~fst::write_fst(read_tsv(..1),path=..2))

## Gene expression data

Lung_GeneExpression <- read_fst(files_fst$Lung_GeneExpression) %>% as_data_frame
Lung_GeneExpression_transposed <- Lung_GeneExpression %>% 
    gather(SampleID, value = gene_expression, starts_with("TCGA")) %>% select(Gene,Chr,gene_expression,SampleID) %>%
    tidyr::unite(gene_chr, Gene,Chr) %>% 
    distinct(SampleID, gene_chr,.keep_all=TRUE) %>% 
    spread(gene_chr,gene_expression) 
write_fst(Lung_GeneExpression_transposed,'Z:/intermediates/Lung_GeneExpression_normal.fst')

