library(GSA)

source('Rcode/cash (Scott\'s stuff)/cash.R')
data=read.table('GSEA/leukemia/expression_data.txt',header=T,sep='\t',comment.char='')
C1=GSA.read.gmt('GSEA/C1/c1.all.v4.0.orig.gmt')$genesets
temp=GSA.read.gmt('GSEA/C1/c1.all.v4.0.entrez.gmt')$genesets