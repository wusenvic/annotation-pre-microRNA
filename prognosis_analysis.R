library(survival)
library(survminer)
library(caret)
library(glmnet)
library(pROC)
library(ggplot2)
library(My.stepwise)
library("GDCRNATools")
library("dplyr")
library('rms')
library('GGally')
library('survivalROC')
library('plotROC')
library('tableone')
library(grpreg)
#setwd("D:/WUSEN/usingcount")
#read files
exprS = read.table("selected_MIR_SAM_TMM.txt",header = T,sep = "\t", row.names = 1)
#expr = read.table("../tcga_mir_1043.txt",header = T,sep = "\t",row.names = 1)
expr = t(expr)
rownames(expr) = gsub("-","\\.",rownames(expr))
colnames(expr) = gsub("-","\\.",colnames(expr))
#miR_expr = gdcVoomNormalization(counts = expr,filter = FALSE)
meta = read.table("../tcga_final_clinic.txt",sep = "\t",row.names = 2,header = T)
#cox= c("hsa.mir.3654","hsa.mir.4653","hsa.mir.671","hsa.mir.6850","hsa.mir.29c","hsa.mir.4680","hsa.mir.5691")
#exprS = expr[cox,rownames(meta)]

#test_mir = read.table("../cohort2_miR_norm.txt",header = T,sep = "\t",row.names = 1)
#test_mir = t(test_mir)
#rownames(test_mir) = gsub("-","\\.",rownames(test_mir))
#colnames(test_mir) = gsub("-","\\.",colnames(test_mir))
#miR_expr = gdcVoomNormalization(counts = expr,filter = FALSE)
#test_meta = read.table("../batch2_meta2.txt",sep = "\t",row.names = 2,header = T)
#test_mir = test_mir[,rownames(test_meta)]

rownames(meta) = gsub("-","\\.",rownames(meta))
colnames(meta) = gsub("-","\\.",colnames(meta))
ex_mir = read.table("tcga_mir_exp.txt",header = T,sep = "\t",row.names = 1)
rownames(ex_mir) = gsub("-","\\.",rownames(ex_mir))
colnames(ex_mir) = gsub("-","\\.",colnames(ex_mir))
ex_mir <- t(ex_mir)
ex_meta = read.table("tcga_final_clinic.txt",header = T,sep = "\t",row.names = 1)
rownames(ex_meta) = gsub("-","\\.",rownames(ex_meta))
colnames(ex_meta) = gsub("-","\\.",colnames(ex_meta))
ex_mir = ex_mir[,rownames(ex_meta)]

#tcga_meta =read.table("tcga_final_clinic.txt",header=T,sep = "\t",row.names = 1)
#tcga_mir = read.table("TCGA-BRCA.mirna.tsv",header = T,sep = "\t",row.names = 1)

#rownames(tcga_meta) = gsub("-","\\.",rownames(tcga_meta))
#colnames(tcga_meta) = gsub("-","\\.",colnames(tcga_meta))
#tcga_mir = tcga_mir[,rownames(tcga_meta)]

# first do a cox regression for each microRNA, filter those with p_value less than 0.15
proj = "CRM-BRCA"
coxfile = paste0(proj,"_cox.Rdata")
logrankfile = paste0(proj,"_log_rank_p.Rdata")

#exprSet = t(exprS)
if(!file.exists(logrankfile)){
  log_rank_p <- apply(exprS , 1 , function(gene){
    meta$group=ifelse(gene>median(gene),'high','low')  
    data.survdiff=survdiff(Surv(OS, event)~group,data=meta)
    p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
    return(p.val)
  })
  log_rank_p=sort(log_rank_p)
  save(log_rank_p,file = logrankfile)
}
load(logrankfile)

lr = names(log_rank_p)[log_rank_p<0.15];length(lr)

if(!file.exists(coxfile)){
cox_results <-apply(exprS , 1 , function(gene){
meta$gene = gene
m = coxph(Surv(OS, event) ~ gene, data =  meta)
beta <- coef(m)
se <- sqrt(diag(vcov(m)))
HR <- exp(beta)
HRse <- HR * se
#summary(m)
tmp <- round(cbind(coef = beta,
se = se, z = beta/se,
p = 1 - pchisq((beta/se)^2, 1),
HR = HR, HRse = HRse,
HRz = (HR - 1) / HRse,
HRp = 1 - pchisq(((HR - 1)/HRse)^2, 1),
HRCILL = exp(beta - qnorm(.975, 0, 1) * se),
HRCIUL = exp(beta + qnorm(.975, 0, 1) * se)), 3)
return(tmp['gene',])
#return(tmp['grouphigh',])#二分类变量
})
cox_results=as.data.frame(t(cox_results))
save(cox_results,file = coxfile)
}
# save related file and results
load(coxfile)
table(cox_results$p<0.01)
table(cox_results$p<0.05)
table(cox_results$p<0.2)
#lr = names(log_rank_p)[log_rank_p<0.01];length(lr)
cox = rownames(cox_results)[cox_results$p<0.05];length(cox)
#length(intersect(lr,cox))
save(lr,cox,file = paste0(proj,"_logrank_cox_gene.Rdata"))
#save(cox,file = paste0(proj,"log_cox_gene.Rdata"))
#cox = rownames(cox_results)[cox_results$p<0.1];length(cox)
cox
exprSet = exprS[cox,]
riskscore <- function(survival_cancer_df, candidate_genes_for_cox, cox_report) {
  risk_score_table <- survival_cancer_df[,candidate_genes_for_cox]
  for(each_sig_gene in colnames(risk_score_table)){
    risk_score_table[,each_sig_gene] <- risk_score_table[,each_sig_gene]*(summary(cox_report)$coefficients[each_sig_gene,1])
  }
  risk_score_table$ 'total_risk_score'<- rowSums(risk_score_table)
  risk_score_table <- cbind(survival_cancer_df[,candidate_genes_for_cox], risk_score_table) %>%
    cbind(survival_cancer_df[,c('sample','OS','event')])
  risk_score_table <- risk_score_table[,c('sample','OS','event', candidate_genes_for_cox, 'total_risk_score')]
  risk_score_table
}
#randomly separate to two parts with seed changed
#calculte 1 sample and then to judge the AUC
multi_ROC <- function(time_vector, risk_score_table){
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime = risk_score_table$OS,
                         status = risk_score_table$event,
                         marker = risk_score_table$total_risk_score,
                         predict.time = single_time, method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP, 
             'Cut_values'=for_ROC$cut.values, 'Time_point'=rep(single_time, length(for_ROC$TP)),
             'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
#We evaluate 11 AUCs between 3-5 years.



times <- vector(mode = "numeric",0)
C_index_list <- vector(mode = "numeric",0)
auc_m <- vector(mode = "numeric",0)

set.seed(3849)
sam<- createDataPartition(meta$event, p = .8,list = FALSE)
#sam <- rownames(meta)
#head(sam)

train <- exprSet[,sam]
test <- exprSet[,-sam]
train_meta <- meta[sam,]
test_meta <- meta[-sam,]
#library(glmnet)
x = as.matrix(t(train))
y = train_meta$event
y <- train_meta[,c("OS","event")]
colnames(y) <- c("time","status")
y$time <- as.double(y$time)
y$status <- as.double(y$status)
#y<- data.matrix(survival::Surv(y$time,y$status))
#y<- survival::Surv(y$time,y$status)
cv_fit <- cv.grpreg(x, y$time, group, penalty="grLasso")

y <-as.data.frame(y)
fit <- glmnet(x=x, y=y$status)
#model_lasso_min <- glmnet(x=x, y=y,lambda=cv_fit$lambda.min)
#model_lasso_1se <- glmnet(x=x, y=y,lambda=cv_fit$lambda.1se)
#choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
#choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
lasso.prob <- predict(cv_fit, X=t(test), s=cv_fit$lambda.min)
re=cbind(event = test_meta$event ,lasso.prob)
re=as.data.frame(re)
colnames(re)=c('event','prob_min')
re$event=as.factor(re$event)
m <- roc(test_meta$event, re$prob_min)
#m2 <- roc(test_meta$event, re$prob_1se)
coefficient <- coef(cv_fit,s=cv_fit$lambda.min)
coefficient2 <- as.data.frame(coefficient)
coefficient2 <- subset(coefficient2,coefficient2$coefficient != 0)
#if (auc(m2)-auc(m)>=0){
#coefficient <- coef(cv_fit,s=cv_fit$lambda.1se)}else{
#coefficient <- coef(cv_fit,s=cv_fit$lambda.min)}
#Active.index <- which(as.numeric(coefficient2) != 0)
#active.coefficient <- as.numeric(coefficient)[Active.index]
sig_gene_multi_cox <- rownames(coefficient2)[2:length(rownames(coefficient2))]
if (length(sig_gene_multi_cox)<=3){
times[seed] = seed
C_index_list[seed] = 0
next}else{
train = t(train)
train = cbind(train,train_meta)
train$sample = rownames(train)
test = t(test)
test = cbind(test,test_meta)
test$sample = rownames(test)
ex_test = cbind(t(ex_mir[cox,]),ex_meta)
ex_test$sample <- rownames(ex_test)
f_test = t(exprSet)
f_test = cbind(f_test,meta)
f_test$sample = rownames(f_test)

datasets = list(train,test,ex_test,f_test)

formula_for_multivariate <- as.formula(paste0('Surv(OS,event)~',paste(sig_gene_multi_cox,sep='',collapse='+')))
#multi_variate_cox <- coxph(formula_for_multivariate,data = train)
#ph_hypo_multi <- cox.zph(multi_variate_cox)
#ph_hypo_table <- ph_hypo_multi$table[-nrow(ph_hypo_multi$table),]
#formula_for_multivariate <- as.formula(paste0('Surv(OS,event)~',paste(rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05],sep='',collapse='+')))
multi_variate_cox2 <- coxph(formula_for_multivariate,data = train)

C_index <- multi_variate_cox2$concordance['concordance']
times[seed] = seed
C_index_list[seed] = C_index
auc_m[seed] = auc(m)
if (C_index >= 0.65){

#correction <- cor(train[,rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05]],method= 'pearson')
#p_cor <- ggpairs(train[,rownames(ph_hypo_table)[ph_hypo_table[,3]>0.05]],axisLabels='show')+theme_bw()+theme(panel.background=element_rect(color='black',size=1,fill='white'),panel.grid=element_blank())


#vif <- rms::vif(multi_variate_cox2)
#sqrt(vif < 2)

candidate_genes_for_cox2 <- c(sig_gene_multi_cox)
train_risk_score_table_multi_cox2 <- riskscore(train, candidate_genes_for_cox2, multi_variate_cox2)
train_for_multi_ROC <- multi_ROC(time_vector = c(365*seq(5,7,1)), risk_score_table = train_risk_score_table_multi_cox2)
train_for_multi_ROC$Time_point <- as.factor(train_for_multi_ROC$Time_point)
train_AUC_max <- max(train_for_multi_ROC$AUC)
if (train_AUC_max >= 0.65){

#train_p_km <- ggsurvplot(train_fit_km, conf.int = F,pval = T,legend.title="train_total risk score",
#           legend.labs=train_fit_km$strata, risk.table = T,palette = c('red','blue'), surv.median.line = 'hv')
#ggsave(paste0(new_dir,"/train_p_km.pdf"),plot=train_p_km,width = 10,height = 8)


#train_p_km <- ggsurvplot(train_fit_km, palette = "jco", risk.table =TRUE, pval =TRUE, conf.int =TRUE)
test_risk_score_table_multi_cox2 <- riskscore(test, candidate_genes_for_cox2, multi_variate_cox2)


test_for_multi_ROC <- multi_ROC(time_vector = c(365*seq(5,7,1)), risk_score_table = test_risk_score_table_multi_cox2)
test_for_multi_ROC$Time_point <- as.factor(test_for_multi_ROC$Time_point)

test_AUC_max <- max(test_for_multi_ROC$AUC)


ex_test_risk_score_table_multi_cox2 <- riskscore(ex_test, candidate_genes_for_cox2, multi_variate_cox2)


ex_test_for_multi_ROC <- multi_ROC(time_vector = c(365*seq(5,10,1)), risk_score_table = ex_test_risk_score_table_multi_cox2)
ex_test_for_multi_ROC$Time_point <- as.factor(ex_test_for_multi_ROC$Time_point)

ex_test_AUC_max <- max(ex_test_for_multi_ROC$AUC)
if (ex_test_AUC_max >= 0.59){
train_AUC_max_time <- train_for_multi_ROC$Time_point[which(train_for_multi_ROC$AUC == train_AUC_max)]
train_AUC_max_time <- train_AUC_max_time[!duplicated(train_AUC_max_time)]
train_AUC_max_time <- train_AUC_max_time[length(train_AUC_max_time)]
train_AUC_max_time <- as.numeric(as.character(train_AUC_max_time))

#train_for_multi_ROC$Time_point <- as.factor(train_for_multi_ROC$Time_point)
#find the optimal cutoff value within the ROC curve of the optimal time point.
#train_optimal_time_ROC_df <- train_for_multi_ROC[which(train_for_multi_ROC$Time_point == train_AUC_max_time),]
#train_cut.off <- train_optimal_time_ROC_df$Cut_values[which.max(train_optimal_time_ROC_df$True_positive-train_optimal_time_ROC_df$False_positive)]
train_cut.off <- median(train_risk_score_table_multi_cox2$total_risk_score)
train_high_low <- (train_risk_score_table_multi_cox2$total_risk_score > train_cut.off)
train_high_low[train_high_low == TRUE] <- 'high'
train_high_low[train_high_low == FALSE] <- 'low'
train_risk_score_table_multi_cox2 <- cbind(train_risk_score_table_multi_cox2, train_high_low)

train_risk_score_table_multi_cox2$event[which(train_risk_score_table_multi_cox2$OS > train_AUC_max_time)] <- 0
train_risk_score_table_multi_cox2$OS[which(train_risk_score_table_multi_cox2$OS > train_AUC_max_time)] <- train_AUC_max_time
train_fit_km <- survfit(Surv(OS, event) ~train_high_low, data = train_risk_score_table_multi_cox2)     

test_AUC_max_time <- test_for_multi_ROC$Time_point[which(test_for_multi_ROC$AUC == test_AUC_max)]
test_AUC_max_time <- test_AUC_max_time[!duplicated(test_AUC_max_time)]
test_AUC_max_time <- test_AUC_max_time[length(test_AUC_max_time)]
test_AUC_max_time <- as.numeric(as.character(test_AUC_max_time))

#test_for_multi_ROC$Time_point <- as.factor(test_for_multi_ROC$Time_point)
#find the optimal cutoff value within the ROC curve of the optimal time point.
test_optimal_time_ROC_df <- test_for_multi_ROC[which(test_for_multi_ROC$Time_point == test_AUC_max_time),]
#test_cut.off <- test_optimal_time_ROC_df$Cut_values[which.max(test_optimal_time_ROC_df$True_positive-test_optimal_time_ROC_df$False_positive)]
test_cut.off <- median(test_risk_score_table_multi_cox2$total_risk_score)
test_high_low <- (test_risk_score_table_multi_cox2$total_risk_score > test_cut.off)
test_high_low[test_high_low == TRUE] <- 'high'
test_high_low[test_high_low == FALSE] <- 'low'
test_risk_score_table_multi_cox2 <- cbind(test_risk_score_table_multi_cox2, test_high_low)

test_risk_score_table_multi_cox2$event[which(test_risk_score_table_multi_cox2$OS > test_AUC_max_time)] <- 0
test_risk_score_table_multi_cox2$OS[which(test_risk_score_table_multi_cox2$OS > test_AUC_max_time)] <- test_AUC_max_time
test_fit_km <- survfit(Surv(OS, event) ~test_high_low, data = test_risk_score_table_multi_cox2)  
test_p_val <- surv_pvalue(test_fit_km)
if (test_p_val$pval < 0.05){
new_dir = paste0("seed",seed)
dir.create(new_dir)
pdf(paste0(new_dir,"/cv_fit.pdf"),height=10,width = 10)
plot(cv_fit)
dev.off()
pdf(paste0(new_dir,"/fit_lambda.pdf"),height=10,width = 10)
plot(fit,xvar = "lambda")
dev.off()
p_forest <- ggforest(model = multi_variate_cox2,data = train,main = "Hazard ratios of candidate microRNAs")
train_p_ROC<-ggplot(train_for_multi_ROC, aes(x = False_positive, y = True_positive, label = Cut_values, color = Time_point)) + geom_roc(labels = F, stat = 'identity', n.cuts = 0) + geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2) + theme_bw()+ theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'), panel.grid = element_blank())+ annotate("text",x = 0.75, y = 0.15, label = paste("AUC max = ", round(train_AUC_max, 2), '\n', 'AUC max time = ', train_AUC_max_time, ' days', sep = ''))
test_p_ROC<-ggplot(test_for_multi_ROC, aes(x = False_positive, y = True_positive, label = Cut_values, color = Time_point)) + geom_roc(labels = F, stat = 'identity', n.cuts = 0) + 
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2)+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'), 
        panel.grid = element_blank())+
  annotate("text",x = 0.75, y = 0.15,
           label = paste("AUC max = ", round(test_AUC_max, 2), '\n', 'AUC max time = ', test_AUC_max_time, ' days', sep = ''))
train_p_km <- ggsurvplot(train_fit_km, palette = "jco", risk.table =TRUE, pval =TRUE, conf.int =TRUE)
test_p_km <- ggsurvplot(test_fit_km, palette = "jco", risk.table =TRUE, pval =TRUE, conf.int =TRUE)
#test_p_km <- ggsurvplot(test_fit_km, conf.int = F,pval = T,legend.title="test_total risk score",
#           legend.labs=c(paste0('&gt;',as.character(round(cut.off,2))), paste0('&lt;=',as.character(round(cut.off,2)))), risk.table = T, 
#           palette = c('red','blue'), surv.median.line = 'hv')

#ggsave(paste0(new_dir,"/test_p_km.pdf"),plot=test_p_km,width = 10,height = 8)

ex_test_AUC_max_time <- ex_test_for_multi_ROC$Time_point[which(ex_test_for_multi_ROC$AUC == ex_test_AUC_max)]
ex_test_AUC_max_time <- ex_test_AUC_max_time[!duplicated(ex_test_AUC_max_time)]
ex_test_AUC_max_time <- ex_test_AUC_max_time[length(ex_test_AUC_max_time)]
ex_test_AUC_max_time <- as.numeric(as.character(ex_test_AUC_max_time))
ex_test_p_ROC<-ggplot(ex_test_for_multi_ROC, aes(x = False_positive, y = True_positive, label = Cut_values, color = Time_point)) + geom_roc(labels = F, stat = 'identity', n.cuts = 0) + 
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2)+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'), 
        panel.grid = element_blank())+
  annotate("text",x = 0.75, y = 0.15,
           label = paste("AUC max = ", round(ex_test_AUC_max, 2), '\n', 'AUC max time = ', ex_test_AUC_max_time, ' days', sep = ''))


#find the optimal cutoff value within the ROC curve of the optimal time point.
ex_test_optimal_time_ROC_df <- ex_test_for_multi_ROC[which(ex_test_for_multi_ROC$Time_point == ex_test_AUC_max_time),]
#ex_test_cut.off <- ex_test_optimal_time_ROC_df$Cut_values[which.max(ex_test_optimal_time_ROC_df$True_positive-ex_test_optimal_time_ROC_df$False_positive)]

ex_test_cut.off <- median(ex_test_risk_score_table_multi_cox2$total_risk_score)
ex_test_high_low <- (ex_test_risk_score_table_multi_cox2$total_risk_score > ex_test_cut.off)
ex_test_high_low[ex_test_high_low == TRUE] <- 'high'
ex_test_high_low[ex_test_high_low == FALSE] <- 'low'
ex_test_risk_score_table_multi_cox2 <- cbind(ex_test_risk_score_table_multi_cox2, ex_test_high_low)

ex_test_risk_score_table_multi_cox2$event[which(ex_test_risk_score_table_multi_cox2$OS > ex_test_AUC_max_time)] <- 0
ex_test_risk_score_table_multi_cox2$OS[which(ex_test_risk_score_table_multi_cox2$OS > ex_test_AUC_max_time)] <- ex_test_AUC_max_time
ex_test_fit_km <- survfit(Surv(OS, event) ~ex_test_high_low, data = ex_test_risk_score_table_multi_cox2)     
ex_test_p_km <- ggsurvplot(ex_test_fit_km, palette = "jco", risk.table =TRUE, pval =TRUE, conf.int =TRUE)
#ex_test_p_km <- ggsurvplot(ex_test_fit_km, conf.int = F,pval = T,legend.title="ex_test_total risk score",
#           legend.labs=c(paste0('&gt;',as.character(round(cut.off,2))), paste0('&lt;=',as.character(round(cut.off,2)))), risk.table = T, 
#           palette = c('red','blue'), surv.median.line = 'hv')


#ggsave(paste0(new_dir,"/ex_test_p_km.pdf"),plot=ex_test_p_km,width = 10,height = 8)

f_test_risk_score_table_multi_cox2 <- riskscore(f_test, candidate_genes_for_cox2, multi_variate_cox2)



f_test_for_multi_ROC <- multi_ROC(time_vector = c(365*seq(5,10,1)), risk_score_table = f_test_risk_score_table_multi_cox2)
f_test_for_multi_ROC$Time_point <- as.factor(f_test_for_multi_ROC$Time_point)

f_test_AUC_max <- max(f_test_for_multi_ROC$AUC)
f_test_AUC_max_time <- f_test_for_multi_ROC$Time_point[which(f_test_for_multi_ROC$AUC == f_test_AUC_max)]
f_test_AUC_max_time <- f_test_AUC_max_time[!duplicated(f_test_AUC_max_time)]
f_test_AUC_max_time <- f_test_AUC_max_time[length(f_test_AUC_max_time)]
f_test_AUC_max_time <- as.numeric(as.character(f_test_AUC_max_time))
f_test_p_ROC<-ggplot(f_test_for_multi_ROC, aes(x = False_positive, y = True_positive, label = Cut_values, color = Time_point)) + geom_roc(labels = F, stat = 'identity', n.cuts = 0) + 
  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2)+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'), 
        panel.grid = element_blank())+
  annotate("text",x = 0.75, y = 0.15,
           label = paste("AUC max = ", round(f_test_AUC_max, 2), '\n', 'AUC max time = ', f_test_AUC_max_time, ' days', sep = ''))


#find the optimal cutoff value within the ROC curve of the optimal time point.
f_test_optimal_time_ROC_df <- f_test_for_multi_ROC[which(f_test_for_multi_ROC$Time_point == f_test_AUC_max_time),]
#f_test_cut.off <- f_test_optimal_time_ROC_df$Cut_values[which.max(f_test_optimal_time_ROC_df$True_positive-f_test_optimal_time_ROC_df$False_positive)]
f_test_cut.off <- median(f_test_risk_score_table_multi_cox2$total_risk_score)
f_test_high_low <- (f_test_risk_score_table_multi_cox2$total_risk_score > f_test_cut.off)
f_test_high_low[f_test_high_low == TRUE] <- 'high'
f_test_high_low[f_test_high_low == FALSE] <- 'low'
f_test_risk_score_table_multi_cox2 <- cbind(f_test_risk_score_table_multi_cox2, f_test_high_low)

f_test_risk_score_table_multi_cox2$event[which(f_test_risk_score_table_multi_cox2$OS > f_test_AUC_max_time)] <- 0
f_test_risk_score_table_multi_cox2$OS[which(f_test_risk_score_table_multi_cox2$OS > f_test_AUC_max_time)] <- f_test_AUC_max_time
f_test_fit_km <- survfit(Surv(OS, event) ~f_test_high_low, data = f_test_risk_score_table_multi_cox2)     
f_test_p_km <- ggsurvplot(f_test_fit_km, palette = "jco", risk.table =TRUE, pval =TRUE, conf.int =TRUE)
#f_test_p_km <- ggsurvplot(f_test_fit_km, conf.int = F,pval = T,legend.title="f_test_total risk score",
#           legend.labs=c(paste0('&gt;',as.character(round(cut.off,2))), paste0('&lt;=',as.character(round(cut.off,2)))), risk.table = T, 
#           palette = c('red','blue'), surv.median.line = 'hv')
#          

#ggsave(paste0(new_dir,"/f_test_p_km.pdf"),plot=f_test_p_km,width = 10,height = 8)
#tcga_risk_score_table_multi_cox2 <- riskscore(tcga, candidate_genes_for_cox2, multi_variate_cox2)
#tcga_for_multi_ROC <- multi_ROC(time_vector = c(365*seq(5,10,1)), risk_score_table = tcga_risk_score_table_multi_cox2)
#tcga_for_multi_ROC$Time_point <- as.factor(tcga_for_multi_ROC$Time_point)

#tcga_AUC_max <- max(tcga_for_multi_ROC$AUC)
#tcga_AUC_max_time <- tcga_for_multi_ROC$Time_point[which(tcga_for_multi_ROC$AUC == tcga_AUC_max)]
#tcga_AUC_max_time <- tcga_AUC_max_time[!duplicated(tcga_AUC_max_time)]
#tcga_AUC_max_time <- tcga_AUC_max_time[length(tcga_AUC_max_time)]
#tcga_AUC_max_time <- as.numeric(as.character(tcga_AUC_max_time))
#tcga_p_ROC<-ggplot(tcga_for_multi_ROC, aes(x = False_positive, y = True_positive, label = Cut_values, color = Time_point)) + geom_roc(labels = F, stat = 'identity', n.cuts = 0) + 
#  geom_abline(slope = 1, intercept = 0, color = 'red', linetype = 2)+
#  theme_bw()+
#  theme(panel.background = element_rect(colour = 'black', size=1, fill = 'white'), 
#        panel.grid = element_blank())+
#  annotate("text",x = 0.75, y = 0.15,
#           label = paste("AUC max = ", round(tcga_AUC_max, 2), '\n', 'AUC max time = ', tcga_AUC_max_time, ' days', sep = ''))


#find the optimal cutoff value within the ROC curve of the optimal time point.
#tcga_optimal_time_ROC_df <- tcga_for_multi_ROC[which(tcga_for_multi_ROC$Time_point == tcga_AUC_max_time),]
#tcga_cut.off <- tcga_optimal_time_ROC_df$Cut_values[which.max(tcga_optimal_time_ROC_df$True_positive-tcga_optimal_time_ROC_df$False_positive)]
#tcga_cut.off <- median(tcga_risk_score_table_multi_cox2$total_risk_score)
#tcga_high_low <- (tcga_risk_score_table_multi_cox2$total_risk_score > tcga_cut.off)
#tcga_high_low[tcga_high_low == TRUE] <- 'high'
#tcga_high_low[tcga_high_low == FALSE] <- 'low'
#tcga_risk_score_table_multi_cox2 <- cbind(tcga_risk_score_table_multi_cox2, tcga_high_low)

#tcga_risk_score_table_multi_cox2$event[which(tcga_risk_score_table_multi_cox2$OS > tcga_AUC_max_time)] <- 0
#tcga_risk_score_table_multi_cox2$OS[which(tcga_risk_score_table_multi_cox2$OS > tcga_AUC_max_time)] <- tcga_AUC_max_time
#tcga_fit_km <- survfit(Surv(OS, event) ~tcga_high_low, data = tcga_risk_score_table_multi_cox2)     
#tcga_p_km <- ggsurvplot(tcga_fit_km, palette = "jco", risk.table =TRUE, pval =TRUE, conf.int =TRUE)



#all_ROC <- list(train_p_ROC,test_p_ROC,ex_test_p_ROC,f_test_p_ROC)
#all_p_ROC <- arrange_ggsurvplots(all_ROC,print = TRUE,ncol = 2,nrow = 2)
#ggsave(paste0(new_dir,"/all_p_ROC.pdf"),plot = all_p_ROC,width = 18,height = 8.5,units = 'in')

write.table(multi_variate_cox2$coefficients,paste0(new_dir,"/coefficiant.txt"))
#ggsave(paste0(new_dir,"/p_cor.pdf"),plot = p_cor,width = 10,height = 8)
ggsave(paste0(new_dir,"/p_forest.pdf"),plot = p_forest,width = 10,height = 8)
write.table(train_risk_score_table_multi_cox2,paste0(new_dir,"/train_risk_score.txt"),sep = "\t")
write.table(train_for_multi_ROC,paste0(new_dir,"/train_for_multi_ROC.txt"),sep = "\t")
ggsave(paste0(new_dir,"/train_p_ROC.pdf"),plot=train_p_ROC,width = 10,height = 8)
pdf(paste0(new_dir,"/train_p_km.pdf"))
print(train_p_km, newpage = FALSE)
dev.off()

write.table(test_risk_score_table_multi_cox2,paste0(new_dir,"/test_risk_score_table_multi_cox2.txt"),sep = "\t")
write.table(test_for_multi_ROC,paste0(new_dir,"/test_for_multi_ROC.txt"),sep = "\t")
ggsave(paste0(new_dir,"/test_p_ROC.pdf"),plot=test_p_ROC,width = 10,height = 8)
pdf(paste0(new_dir,"/test_p_km.pdf"))
print(test_p_km, newpage = FALSE)
dev.off()

write.table(ex_test_risk_score_table_multi_cox2,paste0(new_dir,"/ex_test_risk_score_table_multi_cox2.txt"),sep = "\t")
write.table(ex_test_for_multi_ROC,paste0(new_dir,"/ex_test_for_multi_ROC.txt"),sep = "\t")
ggsave(paste0(new_dir,"/ex_test_p_ROC.pdf"),plot=ex_test_p_ROC,width = 10,height = 8)
pdf(paste0(new_dir,"/ex_test_p_km.pdf"))
print(ex_test_p_km, newpage = FALSE)
dev.off()

write.table(f_test_for_multi_ROC,paste0(new_dir,"/f_test_for_multi_ROC.txt"),sep = "\t")
write.table(f_test_risk_score_table_multi_cox2,paste0(new_dir,"/f_test_risk_score_table_multi_cox2.txt"),sep = "\t")
ggsave(paste0(new_dir,"/f_test_p_ROC.pdf"),plot=f_test_p_ROC,width = 10,height = 8)
pdf(paste0(new_dir,"/f_test_p_km.pdf"))
print(f_test_p_km, newpage = FALSE)
dev.off()

#write.table(tcga_for_multi_ROC,paste0(new_dir,"/tcga_for_multi_ROC.txt"),sep = "\t")
#write.table(tcga_risk_score_table_multi_cox2,paste0(new_dir,"/tcga_risk_score_table_multi_cox2.txt"),sep = "\t")
#ggsave(paste0(new_dir,"/tcga_p_ROC.pdf"),plot=tcga_p_ROC,width = 10,height = 8)
#pdf(paste0(new_dir,"/tcga_p_km.pdf"))
#print(tcga_p_km, newpage = FALSE)
#dev.off()

all_km <- list(train_p_km,test_p_km,ex_test_p_km,f_test_p_km)
#all_p_km <- arrange_ggsurvplots(all_km,print = TRUE,ncol = 2,nrow = 2)
all_p_km <- arrange_ggsurvplots(all_km,print = FALSE,ncol = 2,nrow = 2)
ggsave(paste0(new_dir,"/all_p_km.pdf"),plot = all_p_km,width = 18,height = 12,units = 'in')

}else{next}
}
}
}else{next}
}


