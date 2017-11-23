#' 
#' A script to get all human and mouse genes with theire distance to telomeres
#' 


# load required packages
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(biomaRt)
library(tidyverse)

source("R/get_telomere_dist.R")

#-------------------------------------------------------------------
# A few parameters
#-------------------------------------------------------------------
outPrefix <- "results/telomere_genes_v02"
dir.create("results", showWarnings = FALSE)

#-------------------------------------------------------------------
# get human genes with annotation:
#-------------------------------------------------------------------

# define atributes and download genes from ensemble

ensembl_human <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl")

geneAttributes = c("ensembl_gene_id", "ensembl_transcript_id", 
                   "external_gene_name", "external_gene_source", 
                   "chromosome_name", "transcript_start","transcript_end",  
                   "transcription_start_site", "strand", "gene_biotype")

human_genes = getBM(attributes = geneAttributes, mart = ensembl_human)

# extract seqinfo object
seqInfo <- seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)

# select longest transcript per gene
genesDF_human <- human_genes %>%
  as_tibble() %>% 
  mutate(transcript_size = transcript_end - transcript_start) %>% 
  filter(!is.na(chromosome_name)) %>% 
  group_by(ensembl_gene_id) %>% 
  filter(min_rank(desc(transcript_size)) == 1) %>% 
  # filter for regular chromosomes contained in Bioc object
  filter(paste0("chr", chromosome_name) %in% seqlevels(seqInfo)) %>% 
  ungroup()

# convert into GRanges
tssGR_human <- GRanges(
  paste0("chr", genesDF_human$chromosome_name),
  IRanges(genesDF_human$transcription_start_site, genesDF_human$transcription_start_site),
  strand = ifelse(genesDF_human$strand == 1, "+", "-"),
  seqinfo = seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene)
)
mcols(tssGR_human) <- genesDF_human %>%
  select(-chromosome_name, -transcription_start_site, -strand) %>% 
  as.data.frame()

# add distance to telomere
tssGR_human <- add_telomere_dist(tssGR_human)

humanDF <- mcols(tssGR_human) %>% 
  as.data.frame() %>% 
  as.tibble() %>% 
  select(ensembl_gene_id, external_gene_name, external_gene_source, 
         telomere_dist, telomere, everything())

write_tsv(humanDF, paste0(outPrefix, ".human_genes_with_telomere_distance.tsv"))

#-------------------------------------------------------------------
# get mouse genes with annotation:
#-------------------------------------------------------------------

seqInfo_mouse <- seqinfo(TxDb.Mmusculus.UCSC.mm10.knownGene)

# define atributes and download genes from ensemble
ensembl_mouse <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl")

geneAttributes = c("ensembl_gene_id", "ensembl_transcript_id", 
                   "external_gene_name", "external_gene_source", 
                   "chromosome_name", "transcript_start","transcript_end",  
                   "transcription_start_site", "strand", "gene_biotype")

mouse_genes = getBM(attributes = geneAttributes, mart = ensembl_mouse)

# select longest transcript per gene
genesDF_mouse <- mouse_genes %>%
  as_tibble() %>% 
  mutate(
    transcript_size = transcript_end - transcript_start,
    chr = paste0("chr", chromosome_name)
    ) %>% 
  filter(!is.na(chromosome_name)) %>% 
  group_by(ensembl_gene_id) %>% 
  filter(min_rank(desc(transcript_size)) == 1) %>% 
  # filter for regular chromosomes contained in Bioc object
  filter(chr %in% seqlevels(seqInfo_mouse))

# convert into GRanges
tssGR_mouse <- GRanges(
  genesDF_mouse$chr,
  IRanges(genesDF_mouse$transcription_start_site, genesDF_mouse$transcription_start_site),
  strand = ifelse(genesDF_mouse$strand == 1, "+", "-"),
  seqinfo = seqInfo_mouse)

mcols(tssGR_mouse) <- genesDF_mouse %>%
  select(-chromosome_name, -transcription_start_site, -strand) %>% 
  as.data.frame()

# add distance to telomere
tssGR_mouse <- add_telomere_dist(tssGR_mouse)

mouseDF <- mcols(tssGR_mouse) %>% 
  as.data.frame() %>% 
  as.tibble() %>% 
  select(ensembl_gene_id, external_gene_name, external_gene_source, 
         telomere_dist, telomere, everything())

write_tsv(mouseDF, paste0(outPrefix, ".mouse_genes_with_telomere_distance.tsv"))

#-------------------------------------------------------------------
# get human-mouse ortholog pairs
#-------------------------------------------------------------------

# atributes for orthologs
orthologAttr = c("ensembl_gene_id",
                 paste0("mmusculus", 
                        c("_homolog_ensembl_gene", "_homolog_orthology_type", 
                          "_homolog_subtype", "_homolog_orthology_confidence", 
                          "_homolog_perc_id", "_homolog_perc_id_r1", 
                          "_homolog_dn", "_homolog_ds")))

orthologs = getBM(attributes = orthologAttr, mart = ensembl_human)  


orthologsDF <- orthologs %>% 
  as.tibble() %>% 
  filter(mmusculus_homolog_ensembl_gene != "")

write_tsv(orthologsDF, paste0(outPrefix, ".human_mouse_orthologs.tsv"))

