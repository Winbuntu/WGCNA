annot = read.table("FPKM/genes.attr_table",header = T)

convert_gene_name <- function(x){
  
  annot$gene_short_name[match(x,annot$tracking_id)]
  
}
convert_gene_name("XLOC_000121")
