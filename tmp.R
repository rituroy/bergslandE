## version 2

## ---------------------
## Feature-patient-level
## Feature could be multiple for a gene if there are multiple entries for a gene for a patient


out=matrix(NA,nrow=length(geneList),ncol=sum(!duplicated(clin1$id[samId])),dimnames=list(names(gene)[geneList],paste(c("anyAlt"),clin1$id[samId][!duplicated(clin1$id[samId])],sep="_")))
id=sapply(colnames(out),function(x){strsplit(x,"_")[[1]][2]},USE.NAMES=F)
jj=samId[which(!is.na(clin1$fracReads[samId]))]
jj2=which(sapply(colnames(out),function(x){strsplit(x,"_")[[1]][1]},USE.NAMES=F)=="anyAlt")
grp=paste(clin1$id,clin1$gene,clin1$altType)
grpUniq=unique(grp[jj])
#grpUniq=unique(grp[jj][!grp[jj]%in%grp[jj][duplicated(grp[jj])]])
for (gId in 1:length(grpUniq)) {
    j=jj[which(grp[jj]==grpUniq[gId])]
    if (length(j)!=0) {
        if (any(is.na(match(clin1$id[j],id[jj2])))) print(g1)
        j=j[which.max(clin1$fracReads[j])]
        j2=jj2[match(clin1$id[j],id[jj2])]
        k=which(rownames(out)==clin1$gene[j])
        out[k,j2]=clin1$fracReads[j]
    }
}
fracFP=out



source("code/funcs.R")
out=getFamilyLevelInfo()
datGP_m=out$datGP_m
candGene=out$candGene
k=which(candGene$family2=="SWI.SNF_Modifier.Histone_ATM.ATR")
i=which(clin1$gene%in%candGene$gene[k])
x=paste(clin1$id[i],clin1$gene[i])
clin1[i,][which(x%in%x[duplicated(x)][1]),c("id","disOntTerm2","gene","altType","cdsEffect","fracReads")]

k=which(!is.na(clin1$gene) & !duplicated(clin1$gene) & !is.na(clin1$fracReads))
i=which(clin1$gene%in%candGene$gene[k] & !is.na(clin1$fracReads))
x=paste(clin1$id[i],clin1$gene[i])
x1=unique(x[duplicated(x)])
tmp=rep(NA,length(x1))
out=data.frame(gene=x1,num=tmp,low=tmp,high=tmp,sum=tmp,mean=tmp,stringsAsFactors=F)
for (p in 1:length(x1)) {
    ii=which(x==x1[p])
    out$num[p]=length(ii)
    res=range(clin1$fracRead[i][ii])
    out$low[p]=res[1]
    out$high[p]=res[2]
    out$sum[p]=sum(clin1$fracRead[i][ii])
    out$mean[p]=mean(clin1$fracRead[i][ii])
}
## There are some features with sum(fraction read) > 1
out

