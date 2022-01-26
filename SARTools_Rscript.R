setwd('/analysis/projects/RNA-seq/30_SARTools_Ctrl')

annot <- read.delim('../ref/Biomart_ensembl_for_R.txt')
colnames(annot) <- c('gene_id', 'gene_name', 'gene_description')

write_annot <- function(filename) {
    data <- read.delim(filename, check.names=F)
    merge_data <- merge(unique(annot), data, by.x='gene_id', by.y='Id')
    write.table(merge_data, gsub('tables', 'annot_tables', filename), quote=F, row.names=F, sep='\t')
}

write_annot('tables/Ctrl13vsCtrl0.complete.txt')
write_annot('tables/Ctrl13vsCtrl0.down.txt')
write_annot('tables/Ctrl13vsCtrl0.up.txt')

write_annot('tables/Ctrl20vsCtrl0.complete.txt')
write_annot('tables/Ctrl20vsCtrl0.down.txt')
write_annot('tables/Ctrl20vsCtrl0.up.txt')

write_annot('tables/Ctrl20vsCtrl13.complete.txt')
write_annot('tables/Ctrl20vsCtrl13.down.txt')
write_annot('tables/Ctrl20vsCtrl13.up.txt')

write_annot('tables/Ctrl7vsCtrl0.complete.txt')
write_annot('tables/Ctrl7vsCtrl0.down.txt')
write_annot('tables/Ctrl7vsCtrl0.up.txt')

write_annot('tables/Ctrl7vsCtrl13.complete.txt')
write_annot('tables/Ctrl7vsCtrl13.down.txt')
write_annot('tables/Ctrl7vsCtrl13.up.txt')

write_annot('tables/Ctrl7vsCtrl20.complete.txt')
write_annot('tables/Ctrl7vsCtrl20.down.txt')
write_annot('tables/Ctrl7vsCtrl20.up.txt')
