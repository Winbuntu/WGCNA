#source("http://www.bioconductor.org/biocLite.R")
#biocLite("wiggleplotr")
library("wiggleplotr")
library("dplyr")
library("GenomicRanges")
library("GenomicFeatures")
library("biomaRt")

plotTranscripts(ncoa7_exons, ncoa7_cdss, ncoa7_metadata, rescale_introns = FALSE)

plotTranscripts(ncoa7_exons, ncoa7_cdss, ncoa7_metadata, rescale_introns = TRUE)

###################################
ensembl_mart = useMart("ENSEMBL_MART_ENSEMBL")
ensembl_dataset = useDataset("mmusculus_gene_ensembl",mart=ensembl_mart)
ensembl_dataset
attributes = listAttributes(ensembl_dataset)
head(attributes)
selected_attributes = c("ensembl_transcript_id", "ensembl_gene_id", 
                        "external_gene_name", "strand", 
                        "gene_biotype", "transcript_biotype")
data = getBM(attributes = selected_attributes, mart = ensembl_dataset)
head(data)


data = dplyr::rename(data, 
                     transcript_id = ensembl_transcript_id, 
                     gene_id = ensembl_gene_id, 
                     gene_name = external_gene_name)
head(data)

temporary_file = tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".rds")
saveRDS(data, temporary_file)

transcript_metadata = readRDS(temporary_file)
head(transcript_metadata)

#################################

txdb = makeTxDbFromBiomart(biomart = "ENSEMBL_MART_ENSEMBL", 
                           dataset = "mmusculus_gene_ensembl")

saveDb(txdb, "TranscriptDb_mm10.db")


exons = exonsBy(txdb, by = "tx", use.names = TRUE)
cdss = cdsBy(txdb, by = "tx", use.names = TRUE)


selected_transcripts = transcript_metadata %>%
  dplyr::filter(gene_name == "Bdnf", transcript_biotype == "protein_coding")
tx_ids = selected_transcripts$transcript_id

plotTranscripts(exons[tx_ids], cdss[tx_ids], 
                transcript_metadata, rescale_introns = T)


####################
