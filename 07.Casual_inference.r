#Causal inference for H3K27ac and gene expression
snp <- read.table("potential_casual_snp_from_colocalization.raw",header=T,sep=" ",check.names = F) #0,1 file
gene_expression <- read.table("gene-TMM.bed",header=T)
peak_signal <- read.table("ac-TMM.bed",header=T)
colo <- read.table("colocalization_acQTL_eQTL.txt",header=T,sep="\t")

CausalityTestJM = function(LL,GG,TT){
	no.bootstrap = 50
	sel = (!is.na(LL)) & (!is.na(GG)) & (!is.na(TT))
	dat_f = as.data.frame(cbind(LL,GG,TT),stringsAsFactors=FALSE)
	dat_f = dat_f[sel,]
	names(dat_f) = c("L","G","T")
	Lf = as.factor(dat_f$L)
	dat_f$L = as.integer(Lf) - 1
	llevels = as.integer(levels(as.factor(dat_f$L)))
	dfL = length(llevels) - 1
	pvec = rep(NA,4)
	if(dfL == 2){
	dat_f$L1 = ifelse(dat_f$L == 1,1,0)    
	dat_f$L2 = ifelse(dat_f$L == 2,1,0)    
	fit0 = lm(T ~ 1,data=dat_f)
	fit1 = lm(T ~ L1 + L2,data=dat_f)
	fit2 = lm(G ~ T,data=dat_f)
	fit3 = lm(T ~ G,data=dat_f)
	fit4 = lm(G ~ T + L1 + L2,data=dat_f)
	fit5 = lm(T ~ G + L1 + L2,data=dat_f)
	pvec[1] = anova(fit0,fit1)$"Pr(>F)"[2]
	pvec[2] = anova(fit2,fit4)$"Pr(>F)"[2]
	pvec[3] = anova(fit1,fit5)$"Pr(>F)"[2]
	f_ = anova(fit3,fit5)$F[2]
	fit1G = lm(G ~ L1 + L2,data=dat_f)
	alg = summary(fit1G)$coefficients["(Intercept)",1]
	blg1 = summary(fit1G)$coefficients["L1",1]
	blg2 = summary(fit1G)$coefficients["L2",1]
	alt = summary(fit1)$coefficients["(Intercept)",1]
	blt1 = summary(fit1)$coefficients["L1",1]
	blt2 = summary(fit1)$coefficients["L2",1]
	dat_f$eG = resid(fit1G)
	dat_f$eT = resid(fit1)   
	ss = dim(dat_f)[1]
	fvecr = rep(NA,no.bootstrap)
	fvecr_r = rep(NA,no.bootstrap)
	for(i in 1:no.bootstrap){
		nni <- trunc(1 + ss*runif(ss, 0, 1))
		dat_f$G_ = alg + blg1*dat_f$L1 + blg2*dat_f$L2 + dat_f$eG[nni]
		fit_0 = lm(T ~ G_,data=dat_f)
		fit_1 = lm(T ~ G_ + L1 + L2,data=dat_f)
		fvecr[i] = anova(fit_0,fit_1)$F[2]
		dat_f$T_ = alt + blt1*dat_f$L1 + blt2*dat_f$L2 + dat_f$eT[nni]
		fit_0 = lm(G ~ T_,data=dat_f)
		fit_1 = lm(G ~ T_ + L1 + L2,data=dat_f)
		fvecr_r[i] = anova(fit_0,fit_1)$F[2]}}
	if(dfL == 1){
	dat_f$L1 = ifelse(dat_f$L == 1,1,0)
	fit0 = lm(T ~ 1,data=dat_f)
	fit1 = lm(T ~ L1,data=dat_f)
	fit2 = lm(G ~ T,data=dat_f)
	fit3 = lm(T ~ G,data=dat_f)
	fit4 = lm(G ~ T + L1,data=dat_f)
	fit5 = lm(T ~ G + L1,data=dat_f)
	pvec[1] = anova(fit0,fit1)$"Pr(>F)"[2]
	pvec[2] = anova(fit2,fit4)$"Pr(>F)"[2]
	pvec[3] = anova(fit1,fit5)$"Pr(>F)"[2]
	f_ = anova(fit3,fit5)$F[2]
	fit1G = lm(G ~ L1,data=dat_f)
	alt = summary(fit1)$coefficients["(Intercept)",1]
	blt1 = summary(fit1)$coefficients["L1",1]
	alg = summary(fit1G)$coefficients["(Intercept)",1]
	blg1 = summary(fit1G)$coefficients["L1",1]
	dat_f$eG = resid(fit1G)
	dat_f$eT = resid(fit1)
	ss = dim(dat_f)[1]
	fvecr = rep(NA,no.bootstrap)
	fvecr_r = rep(NA,no.bootstrap)
	for(i in 1:no.bootstrap){
		nni <- trunc(1 + ss*runif(ss, 0, 1))
		dat_f$G_ = alg + blg1*dat_f$L1 + dat_f$eG[nni]
		fit_0 = lm(T ~ G_,data=dat_f)
		fit_1 = lm(T ~ G_ + L1,data=dat_f)
		fvecr[i] = anova(fit_0,fit_1)$F[2]
		dat_f$T_ = alt + blt1*dat_f$L1 + dat_f$eT[nni]
		fit_0 = lm(G ~ T_,data=dat_f)
		fit_1 = lm(G ~ T_ + L1,data=dat_f)
		fvecr_r[i] = anova(fit_0,fit_1)$F[2]}}
	fvecr = fvecr[!is.na(fvecr)]
	df1 = anova(fit3,fit5)$Df[2]
	df2 = anova(fit3,fit5)$Res.Df[2]
	fncp = mean(fvecr,na.rm=TRUE)*(df1/df2)*(df2-df1)-df1
	if(fncp < 0) fncp = 0
	npvals = pf(fvecr,df1,df2,ncp=fncp,lower.tail=TRUE)
	nfvecr = qnorm(npvals)
	npf = pf(f_,df1,df2,ncp=fncp,lower.tail=TRUE)
	zf = qnorm(npf)
	pvec[4] = pnorm(zf,mean=0,sd=sd(nfvecr))
	pvalc = max(pvec)
	fit0G = lm(G ~ 1,data=dat_f)
	pvec1 = rep(NA,4)
	pvec1[1] = anova(fit0G,fit1G)$"Pr(>F)"[2]
	pvec1[2] = anova(fit3,fit5)$"Pr(>F)"[2]
	pvec1[3] = anova(fit1G,fit4)$"Pr(>F)"[2]
	f_ = anova(fit2,fit4)$F[2]
	fvecr_r = fvecr_r[!is.na(fvecr_r)]
	df1 = anova(fit3,fit5)$Df[2]
	df2 = anova(fit3,fit5)$Res.Df[2]
	fncp = mean(fvecr_r,na.rm=TRUE)*(df1/df2)*(df2-df1)-df1
	if(fncp < 0) fncp = 0
	npvals = pf(fvecr_r,df1,df2,ncp=fncp,lower.tail=TRUE)
	nfvecr = qnorm(npvals)
	npf = pf(f_,df1,df2,ncp=fncp,lower.tail=TRUE)
	zf = qnorm(npf)
	pvec1[4] = pnorm(zf,mean=0,sd=sd(nfvecr))
	pvalr = max(pvec1)
	ccall = NA
	ccall = ifelse((pvalc < .05) & (pvalr > .05),1,ccall)
	ccall = ifelse((pvalc > .05) & (pvalr < .05),2,ccall)
	ccall = ifelse((pvalc > .05) & (pvalr > .05),3,ccall)
	ccall = ifelse((pvalc < .05) & (pvalr < .05),0,ccall)
	return(c(pvalc,pvalr,ccall))}

#pvalc = p-value for causal model L -> G -> T
#pvalr = p-value for reactive call L -> T -> G
#ccall = causal call 0:no call, 1:causal, 2:reactive, 3:independent(or other)

for( i in 1:nrow(colo)){ 
L_ID <- colo[i,"lead_snp"]
L <- snp[,L_ID]
ac_ID <- colo[i,"peak_ID"]
AC <- peak_signal[,ac_ID]
rna_ID <- colo[i,"gene_ID"]
G <- gene_expression[,rna_ID]
a <- CausalityTestJM(L,AC,G)
colo[i,c("pvalc")] <- a[1]
colo[i,c("pvalr")] <- a[2]
colo[i,c("ccall")] <- a[3]
}




