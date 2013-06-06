data=read.table('GSEA/p53_collapsed_symbols.gct',header=T,sep='\t',quote='',skip=2,fill=T)
which(data[,1]=='FLJ90650')