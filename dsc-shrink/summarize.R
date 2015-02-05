aggregate(RMSE~method+scenario,res, mean)
library(ggplot2)
ggplot(res,aes(method,RMSE,color=method)) + geom_boxplot() + facet_grid(.~scenario)

p=ggplot(data=res,aes(pi0,pi0_est,colour=method)) +geom_point(shape=1) +
  facet_grid(. ~ scenario) +
  geom_abline(colour = "black") +
  xlab("True pi0") +
  ylab("Estimated pi0")
print(p +scale_y_continuous(limits=c(0,1)) +
        scale_x_continuous(limits=c(0,1)) +
        #scale_colour_manual(values=cbbPalette,breaks=breaks,labels=labels) +
        coord_equal(ratio=1))
