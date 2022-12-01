## initial population genetic analyses for C. maculatus from the poolseq data
library(data.table)

## read population ids
ids<-read.table("PopIds.txt",header=FALSE)

## read allele depth data
ad1<-as.matrix(fread("cmac_ad1.txt",header=FALSE))
ad2<-as.matrix(fread("cmac_ad2.txt",header=FALSE))

## allele freqs and sample size
p<-ad2/(ad2+ad1)
n<-ad1+ad2

apply(n,2,summary)
               V1        V2        V3        V4        V5         V6         V7
#Min.      0.00000   0.00000   0.00000   0.00000    0.0000    0.00000    0.00000
#1st Qu.  44.00000  59.00000  47.00000  51.00000   52.0000   56.00000   52.00000
#Median   58.00000  77.00000  60.00000  66.00000   65.0000   69.00000   64.00000
#Mean     57.86556  76.47619  59.95771  65.89161   64.5953   68.57422   63.98652
#3rd Qu.  70.00000  93.00000  72.00000  79.00000   76.0000   80.00000   74.00000
#Max.    969.00000 958.00000 970.00000 962.00000 1008.0000 1003.00000 1000.00000
               V8       V9        V10
#Min.      0.00000   0.0000    0.00000
#1st Qu.  58.00000  53.0000   59.00000
#Median   73.00000  66.0000   72.00000
#Mean     72.75885  66.0368   72.49382
#3rd Qu.  86.00000  77.0000   84.00000
#Max.    972.00000 991.0000 1004.00000

boxplot(n)

goodSnsp<-apply(n,1,min) > 25
sum(goodSnsp)
#[1] 18571791

p_g<-p[which(goodSnsp),]

L<-dim(p_g)[1]
rsnps<-sample(1:L,1000000,replace=FALSE)
o<-prcomp(t(p_g[rsnps,]),center=TRUE,scale=FALSE)
summary(o)
#Importance of components:
#                            PC1     PC2      PC3      PC4      PC5      PC6
#Standard deviation     109.2921 70.3789 48.08891 39.80724 33.00031 29.88807
#Proportion of Variance   0.4832  0.2004  0.09356  0.06411  0.04406  0.03614
#Cumulative Proportion    0.4832  0.6836  0.77719  0.84129  0.88535  0.92149

pdf("cmacPCA.pdf",width=5,height=5)
par(mar=c(5,5,1,1))
plot(o$x[,1],o$x[,2],type='n',xlab="PC1 (48.3%)",ylab="PC2 (20.0%)",cex.lab=1.3)
text(o$x[,1],o$x[,2],ids[,1],cex=.8)
dev.off()
## fst
Hs<-2 * p_g * (1-p_g)

comps<-matrix(c(1,2,
	 1,3,
	 1,4,
	 2,3,
	 5,6,
	 5,7,
	 6,7,
	 8,9,
	 8,10,
	 9,10,
	 1,5,
	 1,8,
	 5,8),ncol=2,byrow=TRUE)


win<-rep(1:(2*18571),each=500)
xx<-1:length(win)

Fst<-rep(NA,dim(comps)[1])
FstWin<-matrix(NA,nrow=dim(comps)[1],ncol=2*18571)
for(i in 1:dim(comps)[1]){
	a<-comps[i,1];b<-comps[i,2]
	pbar<-(p_g[,a]+p_g[,b])/2
	Ht<-2*pbar*(1-pbar)
	H<-(Hs[,a] + Hs[,b])/2
	Fst[i]<-mean(Ht-H)/mean(Ht)
	num<-tapply(X=Ht[xx]-H[xx],INDEX=win,mean)
	den<-tapply(X=Ht[xx],INDEX=win,mean)
	FstWin[i,]<-num/den
}

data.frame(p1=ids[comps[,1],1],p2=ids[comps[,2],1],fst=round(Fst,3))
#        p1      p2   fst
#   BF1F1C BF2F20L 0.049
#   BF1F1C BF4F20L 0.058
#   BF1F1C BF5F20L 0.109
#  BF2F20L BF4F20L 0.070
#   BZ1F1C BZ4F20L 0.059
#   BZ1F1C BZ5F20L 0.059
#  BZ4F20L BZ5F20L 0.040
#   CA1F1C CA4F20L 0.035
#   CA1F1C CA5F20L 0.036
# CA4F20L CA5F20L 0.028
#  BF1F1C  BZ1F1C 0.159
#  BF1F1C  CA1F1C 0.154
#  BZ1F1C  CA1F1C 0.083
cors<-round(cor(t(FstWin)),2)

colnames(cors)<-paste(ids[comps[,1],1],ids[comps[,2],1])
rownames(cors)<-paste(ids[comps[,1],1],ids[comps[,2],1])

## scatterplots of various combinations, just one here
plot(FstWin[1,],FstWin[2,],xlab=colnames(cors)[1],ylab=colnames(cors)[2])
## take home is that C to L within source is quite repeatable for high Fst, BZ is crazy
## also repeatable to a lesser extent across sources
## and mostly independent of differences among sources

## SNP Fst and allele freq diff
FstSnp<-matrix(NA,nrow=dim(comps)[1],ncol=L)
DpSnp<-matrix(NA,nrow=dim(comps)[1],ncol=L)
for(i in 1:dim(comps)[1]){
	a<-comps[i,1];b<-comps[i,2]
	pbar<-(p_g[,a]+p_g[,b])/2
	DpSnp[i,]<-abs(p_g[,a]-p_g[,b])
	Ht<-2*pbar*(1-pbar)
	H<-(Hs[,a] + Hs[,b])/2
	FstSnp[i,]<-(Ht-H)/Ht
}


save(list=ls(),file="cmacPopGen.rdat")
