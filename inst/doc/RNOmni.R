## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.width=1.5*3, fig.height=1.5*2, fig.align="center",echo=T, warning=F, message=F, cache=T, results='hold');
library(RNOmni);
C1 = require(ggplot2)&require(reshape2);

## ----A00, include = T----------------------------------------------------
# Chi-1 data
y = rchisq(n=1000,df=1);
# Rank-normalize
z = RNOmni::rankNormal(y);

## ----A01, echo=F, fig.width=2*3, eval=C1---------------------------------
# Data frame
df = data.frame("Original"=y,"INT"=z);
df = suppressMessages(reshape2::melt(df));
colnames(df) = c("Data","x");
q = ggplot(data=df) + geom_density(aes(x=x,fill=Data),alpha=0.8);
q = q + theme_bw() + scale_fill_manual(labels=c("Original","Rank Norm"),values=c("coral","royalblue"));
q = q + xlab("Measurement") + ylab("Density");
q = q + ggtitle(expression(chi[1]^2~Phenotype~Before~and~After~INT));
print(q);

## ----A02, include=T, echo=F----------------------------------------------
## Example data
X = RNOmni::X;
cat("Covariates\n");
round(head(X),digits=2);
cat("\n");
cat("Structure Adjustments\n");
S = RNOmni::S;
round(head(S),digits=2);
cat("\n");
cat("Genotype Matrix\n");
G = RNOmni::G;
G[1:6,1:6];
cat("\n");
cat("Sample Minor Allele Frequency\n");
summary(apply(G,MARGIN=2,FUN=mean)/2);
cat("\n");
cat("Phenotypes\n");
Y = RNOmni::Y;
round(head(Y),digits=2);

## ---- include=T, eval=F--------------------------------------------------
#  cat("Omnibus Test, Normal Phenotype, Average Correaltion Method\n");
#  p1.omni.avg = RNOmni::RNOmni(y=Y[,1],G=G,X=X,S=S,method="AvgCorr");
#  round(head(p1.omni.avg),digits=3);
#  cat("\n");
#  cat("Omnibus Test, Normal Phenotype, Bootstrap Correaltion Method\n");
#  set.seed(100);
#  p1.omni.boot = RNOmni::RNOmni(y=Y[,1],G=G,X=X,S=S,method="Bootstrap",B=100);
#  round(head(p1.omni.boot),digits=3);
#  cat("\n");
#  cat("Omnibus Test, T3 Phenotype, Average Correaltion Method\n");
#  p2.omni.avg = RNOmni::RNOmni(y=Y[,2],G=G,X=X,S=S,method="AvgCorr");
#  round(head(p2.omni.avg),digits=3);
#  cat("\n");
#  cat("Omnibus Test, T3 Phenotype, Bootstrap Correaltion Method\n");
#  p2.omni.boot = RNOmni::RNOmni(y=Y[,2],G=G,X=X,S=S,method="Bootstrap",keep.rho=T,B=100);
#  round(head(p2.omni.boot),digits=3);
#  cat("\n");
#  cat("Replicate the Omnibus Test on the T3 Phenotype, Manually Specifying Correlation\n");
#  p2.omni.manual = RNOmni::RNOmni(y=Y[,2],G=G,X=X,S=S,method="Manual",set.rho=p2.omni.boot[,"Corr"],keep.rho=T);
#  round(head(p2.omni.manual),digits=3);
#  cat("\n");

## ---- include=T----------------------------------------------------------
# Basic Association Test
p.bat = RNOmni::BAT(y=Y[,1],G=G,X=X,S=S);
round(head(p.bat),digits=3);

## ---- include=T----------------------------------------------------------
# Direct INT Test
p.dint = RNOmni::DINT(y=Y[,1],G=G,X=X,S=S);
round(head(p.dint),digits=3);

## ---- include=T----------------------------------------------------------
# Partially Indirect INT Test
p.iint = RNOmni::IINT(y=Y[,1],G=G,X=X,S=S);
round(head(p.iint),digits=3);

