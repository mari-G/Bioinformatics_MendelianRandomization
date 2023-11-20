rm(list = ls()); gc()
setwd("/01_forward_mr/")

outcome_factors <- "ukb-d-30710_raw"

if (!dir.exists(outcome_factors)) {dir.create(outcome_factors)}
setwd(outcome_factors)

# 正向孟德尔 柠檬酸盐-------------------------------------------------------------------
# 00.开始分析：去除LD（连锁不平衡）的SNP，暴露因素为柠檬酸盐----------------------------------------------------------------------
library(TwoSampleMR)
exposure <- TwoSampleMR::extract_instruments(outcomes = outcome_factors,##输入暴露因素ID号，读取暴露，筛选SNP
                                             p1 = 5e-08,##当结果不理想时可调大p值
                                             r2 = 0.001,#
                                             kb = 10000)#结果不理想时可调小，最小为10
clumped_exposure <- exposure
outcome <- TwoSampleMR::extract_outcome_data(snps = clumped_exposure$SNP,outcomes = 'ebi-a-GCST90000255')
# outcome <- readRDS("/data/nas1/zhiyq/project/51_GY0336_MR/ieu-b-108_outcome_dat.rds")
dat <- TwoSampleMR::harmonise_data(exposure_dat = clumped_exposure,outcome_dat = outcome)##统一效应等位与效应量
dat <- dat[!dat$SNP=="rs529565",]
write.csv(dat,file = '00.harmonized_data.csv')
# 01.孟德尔随机化分析 -------------------------------------------------------------
# mr_method_list() %>% as.data.frame %>% View()##查看方法，第一列
result <- TwoSampleMR::generate_odds_ratios(mr_res = mr(dat))##随机化分析，结局为疾病时获取OR值,方法默认为5种，可自行选择方法
write.csv(result,'01.MR_expr_outcome_res.csv',row.names = F)
p1 <- TwoSampleMR::mr_scatter_plot(dat,mr_results = result)##dat为harmonise结果，mr_results为mr结果
height <- 8
width <- 8
png("02.MR_expr_outcome_scatter.png",height = height,width = width,res = 300,units = "in")
p1
dev.off()
pdf("02.MR_expr_outcome_scatter.pdf",height = height,width = width)
p1
dev.off()
res_single <- TwoSampleMR::mr_singlesnp(dat)##对每个SNP分别进行2次MR
p2 <- TwoSampleMR::mr_forest_plot(res_single)
p2
height <- 9
width <- 8
png("03.MR_expr_outcome_forest.png",height = height,width = width,res = 300,units = "in")
p2
dev.off()
pdf("03.MR_expr_outcome_forest.pdf",height = height,width = width)
p2
dev.off()
# 02.异质性和水平多效性检验：判断结果的可靠性 ------------------------------------------------------------
##异质性检验
hete <- TwoSampleMR::mr_heterogeneity(dat)
hete$Q_pval
write.csv(hete,"06.heterogeneity_test.csv",row.names = F)
p3 <- TwoSampleMR::mr_funnel_plot(singlesnp_results = mr_singlesnp(dat,all_method = c('mr_ivw')))
p3
height <- 5
width <- 5
png("04.MR_expr_outcome_funnel.png",height = height,width = width,res = 300,units = "in")
p3
dev.off()
pdf("04.MR_expr_outcome_funnel.pdf",height = height,width = width)
p3
dev.off()
##多效性检测
pleiotropy <- TwoSampleMR::mr_pleiotropy_test(dat)
if (pleiotropy$pval>0.05){print("p>0.05,可继续进行分析")}else{print("p<0.05,不可进行后续分析")}
write.csv(pleiotropy,"07.pleiotropy_test.csv",row.names = F)
p4 <- TwoSampleMR::mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
height <- 9
width <- 8
png("05.MR_expr_outcome_leaveoneout.png",height = height,width = width,res = 300,units = "in")
p4
dev.off()
pdf("05.MR_expr_outcome_leaveoneout.pdf",height = height,width = width)
p4
dev.off()
