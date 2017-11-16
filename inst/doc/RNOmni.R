## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.width=1.5*3, fig.height=1.5*2, fig.align="center",echo=T, warning=F, message=F, cache=T, results='hold');

## ----A00, include = T----------------------------------------------------
# Chi-1 data
y = rchisq(n=1000,df=1);
# Rank-normalize
z = RNOmni::rankNormal(y);

## ----A01, echo=F, fig.width=2*3------------------------------------------
library(ggplot2);
library(reshape2);
# Data frame
df = data.frame("Original"=y,"INT"=z);
df = suppressMessages(reshape2::melt(df));
colnames(df) = c("Data","x");
q = ggplot(data=df) + geom_density(aes(x=x,fill=Data),alpha=0.8);
q = q + theme_bw() + scale_fill_manual(labels=c("Original","Rank Norm"),values=c("coral","royalblue"));
q = q + xlab("Measurement") + ylab("Density");
q = q + ggtitle(expression(chi[1]^2~Phenotype~Before~and~After~INT));
print(q);

## ----A02, include=T------------------------------------------------------
library(RNOmni);
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

## ----A03, include=T, echo=F, fig.width=2*3, fig.height=2*1.5*2-----------
# Residuals after projection onto X,S
n = length(y);
D = cbind(1,X,S);
M1 = RcppEigen::fastLmPure(X=D,y=Y[,1]);
e1 = sort(resid(M1));
M2 = RcppEigen::fastLmPure(X=D,y=Y[,2]);
e2 = sort(resid(M2));
M3 = RcppEigen::fastLmPure(X=D,y=rankNormal(Y[,2]));
e3 = sort(resid(M3));
M4a = RcppEigen::fastLmPure(X=cbind(1,X),y=Y[,2]);
M4b = RcppEigen::fastLmPure(X=cbind(1,S),y=rankNormal(resid(M4a)));
e4 = sort(resid(M4b));
z = qnorm(p=seq(from=1,to=1000)/(n+1));
df1 = data.frame("Obs"=e1,"Exp"=z);
df2 = data.frame("Obs"=e2,"Exp"=z);
df3 = data.frame("Obs"=e3,"Exp"=z);
df4 = data.frame("Obs"=e4,"Exp"=z);
# QQ plots
q = ggplot(data=df1,aes(x=Exp,y=Obs)) + geom_point(color="royalblue"); 
q = q + theme_bw() + xlab("Expected Quantile") + ylab("Observed Quantile");
q = q + geom_abline(intercept=0,slope=1,color="gray",linetype="dashed");
q1 = q + ggtitle("QQ for Normal Pheno.");

q = ggplot(data=df2,aes(x=Exp,y=Obs)) + geom_point(color="seagreen");
q = q + theme_bw() + xlab("Expected Quantile") + ylab("Observed Quantile");
q = q + geom_abline(intercept=0,slope=1,color="gray",linetype="dashed");
q2 = q + ggtitle("QQ for T3 Pheno.");

q = ggplot(data=df3,aes(x=Exp,y=Obs)) + geom_point(color="darkseagreen");
q = q + theme_bw() + xlab("Expected Quantile") + ylab("Observed Quantile");
q = q + geom_abline(intercept=0,slope=1,color="gray",linetype="dashed");
q3 = q + ggtitle("QQ for DINT of T3 Pheno.");

q = ggplot(data=df4,aes(x=Exp,y=Obs)) + geom_point(color="darkseagreen");
q = q + theme_bw() + xlab("Expected Quantile") + ylab("Observed Quantile");
q = q + geom_abline(intercept=0,slope=1,color="gray",linetype="dashed");
q4 = q + ggtitle("QQ for PIINT of T3 Pheno.");

cowplot::plot_grid(plotlist=list(q1,q2,q3,q4),nrow=2);

## ----B01, include=T------------------------------------------------------
cat("Omnibus Test, Normal Phenotype, Average Correaltion Method\n");
p1.omni.avg = RNOmni::RNOmni(y=Y[,1],G=G,X=X,S=S,method="AvgCorr");
round(head(p1.omni.avg),digits=3);
cat("\n");
cat("Omnibus Test, Normal Phenotype, Bootstrap Correaltion Method\n");
set.seed(100);
p1.omni.boot = RNOmni::RNOmni(y=Y[,1],G=G,X=X,S=S,method="Bootstrap",B=100);
round(head(p1.omni.boot),digits=3);
cat("\n");
cat("Omnibus Test, T3 Phenotype, Average Correaltion Method\n");
p2.omni.avg = RNOmni::RNOmni(y=Y[,2],G=G,X=X,S=S,method="AvgCorr");
round(head(p2.omni.avg),digits=3);
cat("\n");
cat("Omnibus Test, T3 Phenotype, Bootstrap Correaltion Method\n");
p2.omni.boot = RNOmni::RNOmni(y=Y[,2],G=G,X=X,S=S,method="Bootstrap",keep.rho=T,B=100);
round(head(p2.omni.boot),digits=3);
cat("\n");
cat("Replicate the Omnibus Test on the T3 Phenotype, Manually Specifying Correlation\n");
p2.omni.manual = RNOmni::RNOmni(y=Y[,2],G=G,X=X,S=S,method="Manual",set.rho=p2.omni.boot[,"Corr"],keep.rho=T);
round(head(p2.omni.manual),digits=3);
cat("\n");

## ----B02, include=T, echo=F----------------------------------------------
# Subset omnibus p values
P0 = cbind("AvgCorr"=p1.omni.avg[,3],"Bootstrap"=p1.omni.boot[,3]);
P1 = cbind("AvgCorr"=p2.omni.avg[,3],"Bootstrap"=p2.omni.boot[,3]);
# Size
Q0 = apply(X=(P0<=0.05),MARGIN=2,FUN=mean);
Q0 = data.frame("Method"=names(Q0),"Size"=as.numeric(Q0));
Q1 = apply(X=(P1<=0.05),MARGIN=2,FUN=mean);
Q1 = data.frame("Method"=names(Q1),"Size"=as.numeric(Q1));
# CIs
Q0$L = Q0$Size-2*sqrt(Q0$Size*(1-Q0$Size)/1000);
Q0$U = Q0$Size+2*sqrt(Q0$Size*(1-Q0$Size)/1000);
Q1$L = Q1$Size-2*sqrt(Q1$Size*(1-Q1$Size)/1000);
Q1$U = Q1$Size+2*sqrt(Q1$Size*(1-Q1$Size)/1000);
# Results frame
R = rbind(Q0,Q1);
R$Phenotype = rep(c("Normal","T3"),each=2);
cat("Type I Error of Rank Normal Omnibus Test:\n");
show(data.frame(R[,c(5,1)],round(R[,2:4],digits=3)));

## ----C01, include=T------------------------------------------------------
# Basic Association Test
p1.bat = RNOmni::BAT(y=Y[,1],G=G,X=X,S=S);
round(head(p1.bat),digits=3);

## ----C02, include=T------------------------------------------------------
# Direct INT Test
p1.dint = RNOmni::DINT(y=Y[,1],G=G,X=X,S=S);
round(head(p1.dint),digits=3);

## ----C03, include=T------------------------------------------------------
# Partially Indirect INT Test
p1.piint = RNOmni::PIINT(y=Y[,1],G=G,X=X,S=S);
round(head(p1.piint),digits=3);

## ----C04, echo=F---------------------------------------------------------
## Normal phenotype
P1 = cbind("BAT"=p1.bat,p1.omni.avg[,1:2],"Omni.AvgCorr"=p1.omni.avg[,3],"Omni.Boot"=p1.omni.boot[,3]);
# Size
Q1 = apply(X=(P1<=0.05),MARGIN=2,FUN=mean);
Q1 = data.frame("Method"=names(Q1),"y"=as.numeric(Q1));
# CIs
Q1$L = Q1$y-2*sqrt(Q1$y*(1-Q1$y)/1000);
Q1$U = Q1$y+2*sqrt(Q1$y*(1-Q1$y)/1000);
## T-3 phenotype
p2.bat = BAT(y=Y[,2],G,X,S);
P2 = cbind("BAT"=p2.bat,p2.omni.avg[,1:2],"Omni.AvgCorr"=p2.omni.avg[,3],"Omni.Boot"=p2.omni.boot[,3]);
# Size
Q2 = apply(X=(P2<=0.05),MARGIN=2,FUN=mean);
Q2 = data.frame("Method"=names(Q2),"y"=as.numeric(Q2));
# CIs
Q2$L = Q2$y-2*sqrt(Q2$y*(1-Q2$y)/1000);
Q2$U = Q2$y+2*sqrt(Q2$y*(1-Q2$y)/1000);

## ----C05, echo=F, fig.width=2*3, fig.height=2*1.5*2----------------------
## Plotting
Q1$Method = Q2$Method = factor(Q1$Method,levels=c("Omni.Boot","Omni.AvgCorr","PIINT","DINT","BAT"),ordered=T);
q = ggplot(data=Q1) + geom_linerange(aes(x=Method,ymin=L,ymax=U)) + geom_point(aes(x=Method,y=y),color="royalblue",size=1.5);
q = q + geom_hline(yintercept=0.05,linetype="dashed",color="gray");
q = q + theme_bw() + ylab("Type I Error") + coord_flip();
q1 = q + ggtitle(expression(Estimated~Size~against~Normal~Phenotype~at~alpha==0.05));
q = ggplot(data=Q2) + geom_linerange(aes(x=Method,ymin=L,ymax=U)) + geom_point(aes(x=Method,y=y),color="seagreen",size=1.5);
q = q + geom_hline(yintercept=0.05,linetype="dashed",color="gray");
q = q + theme_bw() + ylab("Type I Error") + coord_flip();
q2 = q + ggtitle(expression(Estimated~Size~against~t[3]~Phenotype~at~alpha==0.05));
cowplot::plot_grid(plotlist=list(q1,q2),ncol=1);

## ----D01, include=T------------------------------------------------------
# Subset to 100 loci
H = G[1:100,];
# Time performance
library(microbenchmark);
microbenchmark(BAT(y=Y[,1],G=H,X=X,S=S),DINT(y=Y[,1],G=H,X=X,S=S),PIINT(y=Y[,1],G=H,X=X,S=S),
               RNOmni(y=Y[,1],G=H,X=X,S=S,method="AvgCorr"));
microbenchmark(RNOmni(y=Y[,1],G=H,X=X,S=S,method="Bootstrap",B=100),times=10);

## ----E01, include=T------------------------------------------------------
# Introduce Missingness
y.m = Y[,1];
y.m[sample(length(y.m),size=10,replace=F)] = NA;
G.m = G;
G.m[sample(length(G.m),size=10000,replace=F)] = NA;
X.m = X;
X.m[sample(length(X.m),size=100,replace=F)] = NA;
S.m = S;
S.m[sample(length(S.m),size=10,replace=F)] = NA;
# Association Testing after Missingness
pm.bat = RNOmni::BAT(y=y.m,G=G.m,X=X.m,S=S.m);
pm.dint = RNOmni::DINT(y=y.m,G=G.m,X=X.m,S=S.m);
pm.piint = RNOmni::PIINT(y=y.m,G=G.m,X=X.m,S=S.m);
pm.omni.avg = RNOmni::RNOmni(y=y.m,G=G.m,X=X.m,S=S.m);
pm.omni.boot = RNOmni::RNOmni(y=y.m,G=G.m,X=X.m,S=S.m,method="Bootstrap",B=100);

## ----E02, include=T, echo=F----------------------------------------------
# Display original results
cat("Normal Phenotype in the Absence of Missingness, Combined Results\n");
A = round(Q1[,2:4],digits=3);
colnames(A)[1] = "Size";
show(data.frame("Method"=Q1$Method,A));
cat("\n");
cat("Normal Phenotype in the Presence of Missingness, Combined Results\n");
# Bind missingness results
P.m = cbind("BAT"=pm.bat,"DINT"=pm.dint,"PIINT"=pm.piint,
            "Omni.Avg"=pm.omni.avg[,"RNOmni"],"Omni.Boot"=pm.omni.boot[,"RNOmni"]);
# Size
QM = apply(X=(P.m<=0.05),MARGIN=2,FUN=mean);
QM = data.frame("Method"=names(QM),"Size"=as.numeric(QM));
# CIs
QM$L = QM$Size-2*sqrt(QM$Size*(1-QM$Size)/1000);
QM$U = QM$Size+2*sqrt(QM$Size*(1-QM$Size)/1000);
# Display missingness results
B = round(QM[,2:4],digits=3);
show(data.frame("Method"=QM$Method,B));
cat("\n");

