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
maf = runif(n=1e3,min=0.05,max=0.50);
G = sapply(maf,function(x){rbinom(n=n,size=2,prob=x)});
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

## ---- include=T----------------------------------------------------------
# Basic Association Test, Normal Phenotype
Results1 = BAT(y=yn,G=G,X=X);
cat("BAT Applied to Normal Phenotype\n");
round(head(Results1),digits=3);
cat("\n");
# Basic Association Test, T3 Phenotype
cat("BAT Applied to T3 Phenotype\n");
Results2  = BAT(y=yt,G=G,X=X);
round(head(Results2),digits=3);

## ---- include=T----------------------------------------------------------
# Direct INT Test, Normal Phenotype
cat("D-INT Applied to Normal Phenotype\n");
Results1 = DINT(y=yn,G=G,X=X);
round(head(Results1),digits=3);
cat("\n");
# Direct INT Test, T3 Phenotype
cat("D-INT Applied to T3 Phenotype\n");
Results2 = DINT(y=yt,G=G,X=X);
round(head(Results2),digits=3);

## ---- include=T----------------------------------------------------------
# Indirect INT Test, Normal Phenotype
cat("I-INT Applied to Normal Phenotype\n");
Results1 = IINT(y=yn,G=G,X=X);
round(head(Results1),digits=3);
cat("\n");
# Indirect INT Test, T3 Phenotype
cat("I-INT Applied to T3 Phenotype\n");
Results2 = IINT(y=yt,G=G,X=X);
round(head(Results2),digits=3);

## ---- include=T, eval=T--------------------------------------------------
# Omnibus INT Test, Normal Phenotype
cat("O-INT Applied to Normal Phenotype\n");
Results1 = OINT(y=yn,G=G,X=X);
round(head(Results1),digits=3);
cat("\n");
# Omnibus INT Test, T3 Phenotype
cat("O-INT Applied to T3 Phenotype\n");
Results2 = OINT(y=yt,G=G,X=X);
round(head(Results2),digits=3);

