## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(echo=T, warning=F, message=F, cache=F, results='hold');
Plot = require(cowplot)&require(ggplot2)&require(reshape2);

## ---- include=T----------------------------------------------------------
library(RNOmni);
# Sample from the chi-1 square distribution
y = rchisq(n=1000,df=1);
# Rank-normalize
z = rankNorm(y);

## ---- echo=F, fig.align="center", fig.width=2*3, eval=Plot---------------
# Data frame
df = data.frame("Original"=y,"INT"=z);
df = suppressMessages(reshape2::melt(df));
colnames(df) = c("Data","x");
q = ggplot(data=df);
q = q+geom_density(aes(x=x,fill=Data),alpha=0.8);
q = q+theme_bw()+scale_fill_manual(labels=c("Original","RankNorm"),values=c("coral","royalblue"));
q = q+xlab("Measurement")+ylab("Density");
q = q+ggtitle(expression(chi[1]^2~Phenotype~Before~and~After~INT));
print(q);

## ---- include=T, echo=T--------------------------------------------------
set.seed(100);
# Sample size
n = 1e3;
## Simulate genotypes
G = replicate(rbinom(n,size=2,prob=0.25),n=1e3);
storage.mode(G) = "numeric";
# Genetic principal components
S = svd(scale(G))$u[,1:4];
S = scale(S);
# Covariates
Z = scale(matrix(rnorm(n*4),nrow=n));
# Overall design
X = cbind(1,Z,S);
# Coefficient
b = c(1,rnorm(n=4,sd=1/sqrt(15)),rnorm(n=4,sd=1/sqrt(60)));
# Linear predictor
h = as.numeric(X%*%b);
# Normal phenotype
yn = h+rnorm(n);
# T-3 phenotype
yt = h+rt(n,df=3)/sqrt(3);

## ---- echo=F, fig.align="center", fig.width=2*3, eval=Plot---------------
# Normal phenotype
q = ggplot(data=data.frame("yn"=yn));
q = q+geom_density(aes(x=yn),alpha=0.8,fill="royalblue");
q = q+theme_bw()+xlab("Phenotype")+ylab("Marginal Density");
q1 = q+ggtitle("Normal Phenotype");
# Log-normal phenotype
q = ggplot(data=data.frame("yt"=yt));
q = q+geom_density(aes(x=yt),alpha=0.8,fill="coral");
q = q+theme_bw()+xlab("Phenotype")+ylab("Marginal Density");
q2 = q+ggtitle("T3 Phenotype");
# Plot
Q = plot_grid(q1,q2,nrow=1);
print(Q);

## ---- include=T, eval=T--------------------------------------------------
cat("Omnibus Test, Normal Phenotype, Average Correaltion Method\n");
pn = RNOmni(y=yn,G=G,X=X,method="AvgCorr");
round(head(pn),digits=3);
cat("\n");
cat("Omnibus Test, Normal Phenotype, Bootstrap Correaltion Method\n");
pn = RNOmni(y=yn,G=G[,1:10],X=X,method="Bootstrap",B=100);
round(head(pn),digits=3);
cat("\n");
cat("Omnibus Test, T3 Phenotype, Average Correaltion Method\n");
pt = RNOmni(y=yt,G=G,X=X,method="AvgCorr");
round(head(pt),digits=3);
cat("\n");
cat("Omnibus Test, T3 Phenotype, Bootstrap Correaltion Method\n");
pt = RNOmni(y=yt,G=G[,1:10],X=X,method="Bootstrap",keep.rho=T,B=100);
round(head(pt),digits=3);
cat("\n");
cat("Replicate the Omnibus Test on the T3 Phenotype, Manually Specifying Correlation\n");
pt = RNOmni(y=yt,G=G,X=X,method="Manual",set.rho=pt[,"Corr"],keep.rho=T);
round(head(pt),digits=3);
cat("\n");

## ---- include=T----------------------------------------------------------
# Basic Association Test, Normal Phenotype
pn = BAT(y=yn,G=G,X=X);
round(head(pn),digits=3);
# Basic Association Test, T3 Phenotype
pt = BAT(y=yt,G=G,X=X);
round(head(pt),digits=3);

## ---- include=T----------------------------------------------------------
# Direct INT Test, Normal Phenotype
pn = DINT(y=yn,G=G,X=X);
round(head(pn),digits=3);
# Direct INT Test, T3 Phenotype
pt = DINT(y=yt,G=G,X=X);
round(head(pt),digits=3);

## ---- include=T----------------------------------------------------------
# Indirect INT Test, Normal Phenotype
pn = IINT(y=yn,G=G,X=X);
round(head(pn),digits=3);
# Indirect INT Test, T3 Phenotype
pt = IINT(y=yt,G=G,X=X);
round(head(pt),digits=3);

