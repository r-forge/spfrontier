
rm(list=ls())
results = data.frame()
for (i in 1:100){
  filename = paste("part_",i,"_20131227_000828.rData",sep="")
  if (file.exists(filename)){
    load(filename)
    obname = paste("ezsim_part_",i,sep="")
    ob = get(obname)
    if (ob$simulation_result[[1]][[1]][1]!=1000) results = rbind(results,ob$simulation_result[[1]][[1]])
  } 
}
colnames(results) <- c("beta0_hat","beta1_hat","beta2_hat","rho","rho2", "sigmaV_hat","sigmaU_hat")
results
save(results,file="results100_sararsfa110.rData")


rm(list=ls())
load("results100_sararsfa110.rData")
plot(density(results$rho))
plot(density(results$rho2))