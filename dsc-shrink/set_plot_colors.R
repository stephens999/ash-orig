library("RColorBrewer")
library("ggplot2")
myColors <- brewer.pal(7,"Set1")
names(myColors) <- c("mixfdr.tnull","ash.hu","ash.n","ash.u","qvalue","locfdr","truth")
colScale <- scale_colour_manual(name = "method",values = myColors)

