load("res_robust.RData")
quantile((res_robust %>% filter(method=="ash.n"))$diff1,c(0.95,0.99,1))
quantile((res_robust %>% filter(method=="ash.hu"))$diff1,c(0.8,0.9,0.95,0.99,1))
quantile((res_robust %>% filter(method=="ash.u"))$diff1,c(0.8,0.9,0.95,0.99,1))
plot(ecdf((res_robust %>% filter(method %in% c("ash.u","ash.hu")))$diff1))
quantile((res_robust %>% filter(method %in% c("ash.u","ash.hu")))$diff1,c(0.8,0.9,0.95,0.99,1))
mean((res_robust %>% filter(method %in% c("ash.u","ash.hu")))$diff1>1)
plot(ecdf((res_robust %>% filter(method %in% c("ash.n")))$diff1))


