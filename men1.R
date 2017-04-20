## Check MEN1 comutation
## Run after heatmap with "_comutatedGenes"

nameValue=data.frame(name="subset",value="allSamples")
nameValue=data.frame(name="subset",value="jeffVetted")
nameValue=data.frame(name="subset",value="fionaVetted")
write.table(nameValue,"config.tmp", sep="=", col.names=F, row.names=F, quote=F)
k=which(nameValue$name=="subset")
datadir="results/"
if (nameValue$value[k]=="fionaVetted") {
    load(paste(datadir,"tmp_fiona166_T5aT7Assays.RData",sep=""))
    nameValue=read.table("config.tmp",sep="=",h=F,quote="",comment.char="",as.is=T,fill=T,col.names=c("name","value"))
    datadirG="results/comutation/fisherTestGeneLevelPv.05/anyCellSize/samSize0/"
} else {
    load(paste(datadir,"tmp_allAssays.RData",sep=""))
    nameValue=read.table("config.tmp",sep="=",h=F,quote="",comment.char="",as.is=T,fill=T,col.names=c("name","value"))
    gene1=unique(clin1$gene[clin1$assayVersion=="T5a"])
    gene2=unique(clin1$gene[clin1$assayVersion=="T7"])
    subsetFlag="_allSamples_T5aT7Assays"
    samId=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7"))
    samId03=which(clin3$id%in%clin1$id[samId])
    datadirG="results/comutation/fisherTestGeneLevelPv.05/anyCellSize/allSamplesSamSize0/"
    k=which(nameValue$name=="subset")
    if (nameValue$value[k]=="jeffVetted") {
        subsetFlag="_jeffVetted_T5aT7Assays"
        samId03=samId03[which(clin3$disOntTerm3[samId03]=="SCLC" | (clin3$grade[samId03]=="High" & clin3$disOntTerm3[samId03]=="combined"))]
        samId=samId[which(clin1$id[samId]%in%clin3$id[samId03])]
        datadirG="results/comutation/fisherTestGeneLevelPv.05/anyCellSize/jeffVettedSamSize0/"
    }
    samId01=samId
}

genesetFlag="pancreas"

x=datGP[,grep("anyAlt_",colnames(datGP))]
x=x[,match(clin3$id[samId03],sub("anyAlt_","",colnames(x)))]
i1=which(rownames(datGP)=="MEN1")
i2=which(rownames(datGP)=="DAXX")
i1
i2
table(x[i1,],x[i2,])
fisher.test(x[i1,],x[i2,])
dim(clin3)
dim(x)
dim(clin3[samId03,])
table(clin3$id[samId03]==sub("anyAlt_","",colnames(x)))
table(clin3$disOntTerm2[samId03])
j=which(clin3$disOntTerm2[samId03]=="colon")
table(x[i1,j],x[i2,j])
j=which(clin3$disOntTerm2[samId03]=="lungSmall")
table(x[i1,j],x[i2,j])
j=which(clin3$disOntTerm2[samId03]=="pancreas")
table(x[i1,j],x[i2,j])
fisher.test(x[i1,j],x[i2,j])


tbl=read.table(paste(datadirG,ifelse(testFlag=="cosine","cor","comut"),"_",genesetFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
pv=read.table(paste(datadirG,"pv_",genesetFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
y=as.matrix(tbl[,-1])
pv=as.matrix(pv[,-1])
rownames(y)=rownames(pv)=tbl$gene
iM=which(rownames(y)=="MEN1")
y[iM,which(y[iM,]!=0)]

tmp=rep(NA,ncol(pv))
out=data.frame(id=rownames(pv),comut=y[iM,],pv=pv[iM,],qv=tmp,bh=tmp,by=tmp,holm=tmp,bonf=tmp,stringsAsFactors=F)
rownames(out)=NULL
i=which(!is.na(out$pv))
out$qv=NA
res=try(qvalue(out$pv[i])$qvalues)
if ("qvalue"%in%class(res)) out$qv[i]
out$bh=p.adjust(out$pv,method="BH")
out$by=p.adjust(out$pv,method="BY")
out$holm=p.adjust(out$pv,method="holm")
out$bonf=p.adjust(out$pv,method="bonferroni")
cat("P-value < 0.05:\n")
print(out[which(out[,"pv"]<.05),])
cat("Q-value < 0.05:\n")
print(out[which(out[,"qv"]<.05),])
cat("BH < 0.05:\n")
print(out[which(out[,"bh"]<.05),])
cat("BY < 0.05:\n")
print(out[which(out[,"by"]<.05),])
cat("Holm's < 0.05:\n")
print(out[which(out[,"holm"]<.05),])
cat("Bonferroni < 0.05:\n")
print(out[which(out[,"bonf"]<.05),])
print(out[which(out$id%in%c("DAXX","TSC2")),])



"
print(out[which(out[,'pv']<.05),])

# fionaVetted
id comut          pv qv       bh by     holm     bonf
48 DAXX     8 0.007963554 NA 0.589303  1 0.589303 0.589303

# jeffVetted
id comut           pv qv         bh        by       holm       bonf
2   RB1     0 8.694432e-03 NA 0.39559668 1.0000000 0.78249892 0.79119336
34 DAXX    17 6.147538e-05 NA 0.00559426 0.0284947 0.00559426 0.00559426

# allSamples
id comut           pv qv         bh        by       holm       bonf
2   RB1     0 0.0045957545 NA 0.28953253 1.0000000 0.57446931 0.57906506
34 DAXX    19 0.0004636045 NA 0.05841416 0.3164564 0.05841416 0.05841416
40 TSC2     9 0.0431311483 NA 0.80208333 1.0000000 1.00000000 1.00000000


"
