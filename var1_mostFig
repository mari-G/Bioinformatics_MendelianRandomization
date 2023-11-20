rm(list = ls())
setwd("/MR/")
if (! dir.exists("./01_MR")){
  dir.create("./01_MR")
}
setwd("./01_MR")
# Fasting glucose || id:ieu-b-113
library(TwoSampleMR)
##正向----------
# finn-b-I9_CORATHER  ukb-d-I9_CORATHER
exposure <- extract_instruments(outcomes = "ukb-a-139",p1 = 5e-08,r2 = 0.001,clump = T,kb=10000)
clumped_exposure <- exposure
outcome <-  extract_outcome_data(snps = clumped_exposure$SNP,outcomes = 'ebi-a-GCST90000255',proxies = TRUE)
dat <- harmonise_data(exposure_dat = clumped_exposure,outcome_dat = outcome)
# dat <- dat[!dat$SNP%in%c('rs10830963','rs2524299','rs7644261','rs11708067'),]
# 
write.csv(dat,file = 'dat.csv',row.names = F)
##直接进行mr，并产生OR值
# mr_method_list()
result <- generate_odds_ratios(mr_res = mr(dat,method_list = c("mr_egger_regression","mr_weighted_median","mr_ivw",
                                                               "mr_simple_mode","mr_weighted_mode")))
write.table(result,file = 'result.mr.xls',sep = '\t',row.names = F,quote = F)

##检测异质性
heter.res <- mr_heterogeneity(dat)
heter.res
##当异质性存在与否，如何选取对应的mr(method_list=?)
##计算Cochran Q是为了量化个体因果效应的异质性，P值≤0.05表明存在异质性，因此应使用随机效应IVW MR分析
##Cochran Q分析评估异质性，并认为如果P值高于0.05且没有异质性证据，则固定效应IVW方法是主要方法。如果存在显著的异质性，则采用随机效应IVW方法（P < 0.05）

write.table(heter.res,file = 'heterogeneity.result.xls',sep = '\t',row.names = F,quote = F)

##自己选方法MR分析，并可视化
library(ggplot2)
p <- mr_scatter_plot(dat,mr_results = mr(dat,method_list = c("mr_egger_regression","mr_weighted_median",
                                                             "mr_simple_mode","mr_weighted_mode","mr_ivw")))
p
pdf(file = '01.Scatter.plot.pdf',w=6,h=5)
print(p)
dev.off()
png(file = '01.Scatter.plot.png',w=6,h=5,units = 'in',res = 600)
print(p)
dev.off()

# ##水平多效性
# presso.result <- run_mr_presso(dat,NbDistribution = 1000)
# presso.result <- as.data.frame(presso.result)
# presso.result
# write.table(presso.result,file = 'presso.result.xls',sep = '\t',row.names = F,quote = F)

##水平多效性----
pleiotropy.result <- mr_pleiotropy_test(dat)
pleiotropy.result
write.table(file = 'pleiotropy.result.xls',pleiotropy.result,sep = '\t',row.names = F,quote = F)

##森林图-------
res_single <- mr_singlesnp(dat)
p2 <- mr_forest_plot(res_single)
pdf(file = '02.forest.pdf',w=6,h=6)
print(p2[[1]])
dev.off()
png(file = '02.forest.png',w=6,h=6,units = 'in',res = 600)
print(p2[[1]])
dev.off()


##绘制漏斗图-------
p <- mr_funnel_plot(singlesnp_results = mr_singlesnp(dat,all_method = c('mr_ivw')))
pdf(file = '03.funnel.pdf',w=5,h=5)
print(p)
dev.off()
png(file = '03.funnel.png',w=5,h=5,units = 'in',res = 600)
print(p)
dev.off()


p <- mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
pdf('04.leaveoneout.pdf',w=6,h=6)
print(p)
dev.off()
png(file = '04.leaveoneout.png',w=6,h=6,units = 'in',res = 600)
print(p)
dev.off()


snp <- dat$SNP%>%as.data.frame()
write.csv(snp,file = 'snp.expouse.csv',row.names = F)


