aggregate(RMSE~method+scenario,res, mean)
library(ggplot2)
library(dplyr)

### Plot RMSE boxplot
ggplot(res,aes(method,RMSE,color=method)) + geom_violin() + facet_grid(.~scenario)

### Plot RMSE boxplot for subset of methods
ggplot(filter(res,scenario %in% c("An","Bn","Cn")),aes(method,RMSE,color=method)) + geom_violin() + facet_grid(.~scenario)

### or boxplot of elapsed time
ggplot(filter(res,!(scenario %in% c("hard-b"))),aes(method,elapsed,color=method)) + geom_boxplot() + facet_grid(.~scenario)


## plot true pi0 vs estimated pi0

p=ggplot(data=filter(res,scenario %in% c("hard","hard-b","easy")),aes(pi0,pi0_est,colour=method)) +geom_point(shape=1) +
  facet_grid(. ~ scenario) +
  geom_abline(colour = "black") +
  xlab("True pi0") +
  ylab("Estimated pi0")
print(p +scale_y_continuous(limits=c(0,1)) +
        scale_x_continuous(limits=c(0,1)) +
        #scale_colour_manual(values=cbbPalette,breaks=breaks,labels=labels) +
        coord_equal(ratio=1))
