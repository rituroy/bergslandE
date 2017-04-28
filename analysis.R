## Can trust all lung entries

####################################################################
####################################################################
## Section 1

computerFlag="cluster"
computerFlag=""

## ----------------------------------------------
if (computerFlag=="cluster") {
	setwd("/home/royr/project/EmilyBergsland")
} else {
	dirSrc="/Users/royr/UCSF/"
	dirSrc2=dirSrc
    dirSrc3="code/"
    dirSrc3="/Users/royr/Downloads/bergslandE-net_2015/"
	setwd(paste(dirSrc2,"EmilyBergsland",sep=""))
}


## ----------------------------------------------
## ----------------------------------------------

capWords=function(s, strict=FALSE) {
    cap=function(s) paste(toupper(substring(s, 1, 1)),
    {s=substring(s, 2); if(strict) tolower(s) else s},
    sep="", collapse=" " )
    sapply(strsplit(s, split=" "), cap, USE.NAMES=!is.null(names(s)))
}

tolowerWords=function(s, strict=FALSE) {
    cap=function(s) paste(tolower(substring(s, 1, 1)),
    {s=substring(s, 2); if(strict) tolower(s) else s},
    sep="", collapse=" " )
    sapply(strsplit(s, split=" "), cap, USE.NAMES=!is.null(names(s)))
}

####################################################################
####################################################################
## Section 2

## ----------------------------------------------
## ----------------------------------------------
## Read files

if (computerFlag=="cluster") {
    datadir="data/"
} else {
    datadir="docs/"
}
clin1=read.table(paste(datadir,"Neuroendocrine_variants_de_identified.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
clin2=read.table(paste(datadir,"Copy of Neuroendocrine_variants_TRFs_GIonly_NoTRFnumbers.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
clinN=read.table(paste(datadir,"Neuroendocrine_likely_passed_variants_grading_de_id_assayVersion(1).txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
clinC=read.table(paste(datadir,"clinical_assays_genes(2).txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
clinK=read.table(paste(datadir,"Neuroendocrine_variants_de_id_assayVersion_known_assay.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
clinV=read.table(paste(datadir,"Neuroendocrine_variants_de_id_assayVersion_vus_assay.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
candGene=read.table(paste(datadir,"targetableGenes.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

candGene=candGene$gene[!duplicated(candGene$gene)]
candGene=gsub("\"","",candGene)
gene=geneAlias=c()
for (i in 1:length(candGene)) {
    x=strsplit(candGene[i],", ")[[1]]
    geneAlias=c(geneAlias,rep(candGene[i],length(x)))
    gene=c(gene,x)
}
candGene=data.frame(gene,geneAlias,stringsAsFactors=F)
##save(candGene,file="candGene.RData")

names(clin1)[match(c("ID","TISSUE","PRIMARYTISSUE","AGE","GENDER","STUDY","DISEASE_ONTOLOGY_TERM","QC_RESULT","GENE","PARTNER_GENE","ALTERATION_TYPE","CDS_EFFECT","PROTEIN_EFFECT","COPY_NUMBER","FRACTION_READS","VAR_STATUS"),names(clin1))]=
c("id","tissue","priTissue","age","gender","study","disOntTerm","qcResult","gene","partnerGene","altType","cdsEffect","proteinEffect","cn","fracReads","varStatus")
names(clin2)[match(c("ID","DISEASE_ONTOLOGY_TERM"),names(clin2))]=c("id","disOntTerm")
names(clinN)[match(c("ID","TISSUE","AGE","GENDER","STUDY","DISEASE_ONTOLOGY_TERM","QC_RESULT","GENE","PARTNER_GENE","ALTERATION_TYPE","CDS_EFFECT","PROTEIN_EFFECT","COPY_NUMBER","FRACTION_READS","VAR_STATUS","assayVersion","GRADE"),names(clinN))]=
c("id","tissue","age","gender","study","disOntTerm","qcResult","gene","partnerGene","altType","cdsEffect","proteinEffect","cn","fracReads","varStatus","assayVersion","grade")
names(clinK)[match(c("ID","TISSUE","PRIMARYTISSUE","AGE","GENDER","STUDY","DISEASE_ONTOLOGY_TERM","QC_RESULT","GENE","PARTNER_GENE","ALTERATION_TYPE","CDS_EFFECT","PROTEIN_EFFECT","COPY_NUMBER","FRACTION_READS","VAR_STATUS","assayVersion"),names(clinK))]=
c("id","tissue","priTissue","age","gender","study","disOntTerm","qcResult","gene","partnerGene","altType","cdsEffect","proteinEffect","cn","fracReads","varStatus","assayVersion")
names(clinV)[match(c("ID","TISSUE","PRIMARYTISSUE","AGE","GENDER","STUDY","DISEASE_ONTOLOGY_TERM","QC_RESULT","GENE","PARTNER_GENE","ALTERATION_TYPE","CDS_EFFECT","PROTEIN_EFFECT","COPY_NUMBER","FRACTION_READS","VAR_STATUS","assayVersion"),names(clinV))]=
c("id","tissue","priTissue","age","gender","study","disOntTerm","qcResult","gene","partnerGene","altType","cdsEffect","proteinEffect","cn","fracReads","varStatus","assayVersion")

for (thisFlag in c("1","2","N","K","V","C")) {
    switch(thisFlag,
    "1"={clin=clin1},
    "2"={clin=clin2},
    "N"={clin=clinN},
    "K"={clin=clinK},
    "V"={clin=clinV},
    "C"={clin=clinC}
    )
    for (k in 1:ncol(clin)) {
        if (class(clin[,k])=="character") {
            clin[,k]=gsub("\"","",clin[,k])
            clin[which(clin[,k]==""),k]=NA
        }
    }
    k=which(names(clin)=="disOntTerm")
    if (length(k)==1) {
        clin[,k]=gsub(",","",clin[,k])
    }
    switch(thisFlag,
    "1"={clin1=clin},
    "2"={clin2=clin},
    "N"={clinN=clin},
    "K"={clinK=clin},
    "V"={clinV=clin},
    "C"={clinC=clin}
    )
}

## Check that clin1 & clinK agree
grp1=paste(clin1$id,clin1$disOntTerm,clin1$gene,clin1$altType,clin1$cdsEffect)
grp2=paste(clinK$id,clinK$disOntTerm,clinK$gene,clinK$altType,clinK$cdsEffect)
table(duplicated(grp1))
table(duplicated(grp2))
table(grp1==grp2)
k=match(names(clin1),names(clinK)); k1=which(!is.na(k)); k2=k[k1]
for (k in 1:length(k1)) {
    if (any(is.na(clin1[,k1[k]])!=is.na(clinK[,k2[k]])) | any(clin1[,k1[k]]!=clinK[,k2[k]],na.rm=T)) print(k)
}
clin1$assayVersion=clinK$assayVersion

## Compare clinK clinN
if (T) {
    grp1=paste(clin1$id,clin1$disOntTerm,clin1$gene,clin1$altType,clin1$cdsEffect)
    grp2=paste(clinN$id,clinN$disOntTerm,clinN$gene,clinN$altType,clinN$cdsEffect)
    table(duplicated(grp1))
    table(duplicated(grp2))
    table(grp2%in%grp1)
    for (k in 1:length(k1)) {
        if (any(is.na(clin1[,k1[k]])!=is.na(clinK[,k2[k]])) | any(clin1[,k1[k]]!=clinK[,k2[k]],na.rm=T)) print(k)
    }
    clin1$grade=NA
    j=match(grp1,grp2); j1=which(!is.na(j)); j2=j[j1]
    clin1$grade[j1]=clinN$grade[j2]
}
if (F) {
    clin1$grade=NA
    j=match(clin1$id,clinN$id); j1=which(!is.na(j)); j2=j[j1]
    table(is.na(j))
    clin1$grade[j1]=clinN$grade[j2]
}

clin1$age[is.infinite(clin1$age)]=NA

## Assay version
tmpC=rep("",nrow(clinC)*ncol(clinC))
out=data.frame(gene=tmpC,assayVersion=tmpC)
out$gene=c(as.matrix(clinC))
out$assayVersion=rep(names(clinC),each=nrow(clinC))
out=out[!is.na(out),]

table(duplicated(clin1$id),duplicated(paste(clin1$id,clin1$disOntTerm)))
table(duplicated(clin2$id))
id=sort(unique(clin1$id))
k=c()
for (j1 in 1:length(id)) {
    j=which(clin1$id==id[j1])
    if (sum(!duplicated(clin1$disOntTerm[j]))>1) k=c(k,j1)
}
## No patient with multiple disease

i=which(clin1$qcResult%in%c("Pass","Qualified"))
clin1=clin1[i,]

table(clin1$disOntTerm[which(clin1$id%in%clin2$id)])
table(clin1$disOntTerm[which(!clin1$id%in%clin2$id)])
#clin1=clin1[which(clin1$id%in%clin2$id),]
clin2=clin2[which(clin2$id%in%clin1$id),]

for (k in 1:ncol(clin1)) {
    cat("\n\n============",k,names(clin1)[k],"\n")
    x=sum(!duplicated(clin1[,k]))
    if (x<26) {
        print(table(clin1[,k],exclude=NULL))
    } else {
        if (x<101) {
            print(sort(unique(clin1[,k])))
        } else if (class(clin1[,k])=="numeric") {
            print(summary(clin1[,k]))
        }
    }
}

clin3=clin1[!duplicated(clin1$id),]
table(clin3$disOntTerm[which(clin3$id%in%clin2$id)])
table(clin3$disOntTerm[which(!clin3$id%in%clin2$id)])

## ----------
grpUniq=sort(unique(clin1$assayVersion))
out=matrix(nrow=nrow(clin3),ncol=length(grpUniq),dimnames=list(clin3$id,grpUniq))
clin3$assayVersion=""
for (k in 1:length(grpUniq)) {
    x=table(clin1$id[which(clin1$assayVersion==grpUniq[k])])
    j=match(names(x),rownames(out))
    out[j,k]=x
    clin3$assayVersion[j]=grpUniq[k]
}
x=apply(out,1,function(x) {sum(x!=0,na.rm=T)})
summary(x)
## No patient on multiple assay
res=table(disOntTerm=clin3$disOntTerm,assayVersion=clin3$assayVersion)
res2=t(apply(res,1,function(x) {round(x/sum(x),2)}))
res
res2

## ----------
res=table(clin3$disOntTerm)
samId13=which(clin1$disOntTerm%in%names(res)[res>=20])

for (subsetFlag in c("","_20pat")) {
    samId=1:nrow(clin1)
    if (subsetFlag=="_20pat") samId=samId13

    samId=samId[which(!is.na(clin1$altType[samId]) & !duplicated(paste(clin1$id,clin1$gene)[samId]))]
    samId=samId[which(clin1$assayVersion[samId]%in%c("T5a","T7"))]
    gene=table(clin1$gene[samId])
    gene[order(gene,decreasing=T)][1:10]
    geneList=order(gene,decreasing=T); geneList=geneList[which(gene[geneList]>9)]
    geneList=order(gene,decreasing=T)

    grp=clin1$disOntTerm[samId]
    grpUniq=sort(unique(clin1$disOntTerm[samId]))

    out=matrix(0,nrow=length(geneList),ncol=length(grpUniq),dimnames=list(names(gene),grpUniq))
    j=samId[which(clin1$assayVersion[samId]=="T5a")]
    x=table(clin1$gene[j],clin1$disOntTerm[j])
    i=match(rownames(x),rownames(out))
    j=match(colnames(x),colnames(out))
    out[i,j]=x
    j=samId[which(clin1$assayVersion[samId]=="T7")]
    x=table(clin1$gene[j],clin1$disOntTerm[j])
    i=match(rownames(x),rownames(out))
    j=match(colnames(x),colnames(out))
    out[i,j]=out[i,j]/x
    out0=out

    out=matrix(0,nrow=length(geneList),ncol=length(grpUniq),dimnames=list(names(gene),grpUniq))
    j=samId[which(clin1$assayVersion[samId]=="T5a")]
    x=table(clin1$gene[j],clin1$disOntTerm[j])
    i=match(rownames(x),rownames(out))
    j=match(colnames(x),colnames(out))
    out[i,j]=x
    out5a=out

    out=matrix(0,nrow=length(geneList),ncol=length(grpUniq),dimnames=list(names(gene),grpUniq))
    j=samId[which(clin1$assayVersion[samId]=="T7")]
    x=table(clin1$gene[j],clin1$disOntTerm[j])
    i=match(rownames(x),rownames(out))
    j=match(colnames(x),colnames(out))
    out[i,j]=x
    out7=out


    out[geneList,][1:50,1:5]

    tbl=apply(round(out0[geneList,],3),c(1,2),as.character)
    tbl[tbl=="NaN"]=""
    tbl=cbind(gene=names(gene)[geneList],tbl)
    write.table(tbl,file=paste("proportion_assayVersion_T5abyT7",subsetFlag,".txt",sep=""),append=F,col.names=T,row.names=F,sep="\t",quote=F)
    
    if (F) {
    if (subsetFlag=="_20pat") {
        out=out0
        summary(c(out0[!is.infinite(c(out0))]))
        lim=c(0,9)
        nm="Lung small cell undifferentiated carcinoma"
        k1=which(colnames(out0)==nm)
        i=geneList[1:100]
        for (k in which(colnames(out0)!=nm)) {
            png(paste("proportionPlot_assayVersion_T5abyT7_",colnames(out0)[k],"vs",nm,sep=""))
            plot(out[i,k1],out[i,k],xlim=lim,ylim=lim,main="Proportion of patients for top 100 genes on assays T5a/T7",xlab=nm,ylab=colnames(out0)[k])
            abline(c(0,1))
            dev.off()
        }
    }
    }

}

## -------------------------------------------------------
## -------------------------------------------------------
## Check data

## Mutations have 2 entries - cdsEffect & proteinEffect
table(cdsEffect=!is.na(clin1$cdsEffect),proteinEffect=!is.na(clin1$proteinEffect))
"
            proteinEffect
cdsEffect   FALSE TRUE
FALSE       2677    0
TRUE            0 5495
"

# An entry is for either mutation or cn, not both
# Some entries do not have either type of alteration
# Upto ~ 315 genes tested for alterations for each patient
table(cdsEffect=!is.na(clin1$cdsEffect),cn=!is.na(clin1$cn))
"
cn
cdsEffect FALSE TRUE
FALSE       433 2244
TRUE       5495    0
"

## Fraction reads for mutations only, not cn
summary(clin1$fracReads[!is.na(clin1$cdsEffect)])
"
Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
0.0100  0.3000  0.4600  0.4882  0.7000  1.0000
"
summary(clin1$fracReads[is.na(clin1$cdsEffect)])
"
Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
NA      NA      NA     NaN      NA      NA    2677
"


table(cdsEffect=!is.na(clin1$cdsEffect),cn=!is.na(clin1$cn),clin1$altType,exclude=NULL)
"
, ,  = amplification

            cn
cdsEffect FALSE TRUE <NA>
FALSE     0 1616    0
TRUE      0    0    0
<NA>      0    0    0

, ,  = loss

            cn
cdsEffect FALSE TRUE <NA>
FALSE     0  628    0
TRUE      0    0    0
<NA>      0    0    0

, ,  = rearrangement

            cn
cdsEffect FALSE TRUE <NA>
FALSE   242    0    0
TRUE      0    0    0
<NA>      0    0    0

, ,  = short variant

            cn
cdsEffect FALSE TRUE <NA>
FALSE     0    0    0
TRUE   5495    0    0
<NA>      0    0    0

, ,  = NA

            cn
cdsEffect FALSE TRUE <NA>
FALSE   191    0    0
TRUE      0    0    0
<NA>      0    0    0
"
## -------------------------------------------------------

varList=c("altType")
k=1
j2=c()
for (j3 in 1:nrow(clin3)) {
    j1=which(clin1$id==clin3$id[j3])
    if (any(duplicated(clin1[j1,varList[k]]))) j2=c(j2,j3)
}
length(j2)

gene=table(clin1$gene)
gene[order(gene,decreasing=T)][1:10]
varList=c("altType")
k=1
for (g1 in order(gene,decreasing=T)[1:10]) {
    cat("\n\n============",g1,names(gene)[g1],gene[g1],"\n")
    j2=j22=c()
    for (j3 in 1:nrow(clin3)) {
        j1=which(clin1$id==clin3$id[j3] & clin1$gene==names(gene)[g1])
        if (any(duplicated(clin1[j1,varList[k]][!is.na(clin1[j1,varList[k]])]))) j2=c(j2,j3)
        if (sum(!duplicated(clin1[j1,varList[k]][!is.na(clin1[j1,varList[k]])]))>1) j22=c(j22,j3)
    }
    print(length(j2))
    print(length(j22))
}

id=unique(clin1$disOntTerm[which(clin1$disOntTerm%in%clin1$disOntTerm[which(is.na(clin1$qcResult))])])
table(clin1$disOntTerm[which(is.na(clin1$qcResult) & clin1$disOntTerm%in%id)])
"
Ovary neuroendocrine carcinoma
4
Soft tissue small cell neuroendocrine carcinoma
4
"
table(clin1$disOntTerm[which(!is.na(clin1$qcResult) & clin1$disOntTerm%in%id)])
"
Ovary neuroendocrine carcinoma
40
Soft tissue small cell neuroendocrine carcinoma
8
"
for (x in c(sort(unique(clin1$qcResult)),NA)) {
    cat("\n------- QC Result: ",x,"\n",sep="")
    j=which(!is.na(clin1$cn))
    j=which(is.na(clin1$cdsEffect) & is.na(clin1$cn))
    j=which(!is.na(clin1$cdsEffect))
    if (is.na(x)) j=j[which(is.na(clin1$qcResult[j]))] else j=j[which(clin1$qcResult[j]==x)]
    print(summary(clin1$fracRead[j]))
}
"
Fraction reads for sequence mutation:

------- QC Result: Fail
Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
0.0200  0.1725  0.2900  0.3217  0.4075  0.9100

------- QC Result: Pass
Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
0.0100  0.3100  0.4600  0.4963  0.7100  1.0000

------- QC Result: Qualified
Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
0.020   0.110   0.340   0.353   0.540   0.990

------- QC Result: NA
Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
0.2700  0.3300  0.3700  0.4875  0.5275  0.9400
"

"
PD-L1 protein: gene CD274
PD-L2 protein: gene PDCD1LG2, CD273
"
gene=c("CD274","PDCD1LG2","CD273")
table(clin1$gene[which(clin1$gene%in%gene)])
"
CD274 PDCD1LG2
9       10
"

id=unique(clin1$id[is.na(clin1$altType)])
length(id)
"191"
table(clin1$id%in%id,alteration=!is.na(clin1$altType))
"
        alteration
        FALSE TRUE
FALSE     0 7981
TRUE    191    0
"
## All diseases have some patients with no mutation
table(clin1$disOntTerm[clin1$id%in%id]%in%clin1$disOntTerm[!clin1$id%in%id])
"
TRUE
191
"
table(clin1$disOntTerm[clin1$id%in%id])

## Genes with varStatus are not in genes with no varStatus
table(clin1$gene[!is.na(clin1$varStatus)]%in%clin1$gene[is.na(clin1$varStatus)])
## Some genes have varStatus known & likely
table(duplicated(clin1$gene[which(clin1$varStatus=="known")][(clin1$gene[which(clin1$varStatus=="known")]%in%clin1$gene[which(clin1$varStatus=="likely")])]))
"
FALSE  TRUE
139  3248
"

## Some genes have multiple types of alterations
x=table(clin1$gene,clin1$altType)
y=apply(x,1,function(x) {sum(x!=0,na.rm=T)})
table(y)

x=table(paste(clin1$disOntTerm,clin1$gene),clin1$altType)
y=apply(x,1,function(x) {sum(x!=0,na.rm=T)})


table(fionaCol1=clin32$fionaCol1,fionaCol2=clin32$fionaCol2,disOntTerm=clin32$disOntTerm)

## ----------------------------------------------
## ----------------------------------------------
## Tabulate

#load("results/tmp_allAssays.RData")

subsetFlag="_fiona102"
subsetFlag="_fiona102_T5aT7Assays"
subsetFlag="_fiona166_T5aT7Assays"
subsetFlag="_allAssays"

tbl=read.table(paste("docs/Sent to Ritu 4-25 Full sample set as reviewed by Fiona - G3 WHO tumors as Y ( 3).txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
names(tbl)[match(c("ID","Gender","Age","Sample.Site","Disease.Ontology","Fiona.Column.1"),names(tbl))]=c("id","gender","age","Sample.Site","disOntTerm","fionaCol1")
kk=ceiling(nrow(clin3)/nrow(tbl))
clin32=NULL
for (k in 1:kk) clin32=rbind(clin32,tbl)
clin32=clin32[1:nrow(clin3),]
for (k in 1:ncol(clin32)) clin32[,k]=NA
clin32$id=clin3$id
j=match(tbl$id,clin32$id); j1=which(!is.na(j)); j2=j[j1]
clin32[j2,]=tbl[j1,]
j=which(!is.na(clin32$disOntTerm))
#clin32$disOntTerm2[j]=clin3$disOntTerm2[j]
#clin32$disOntTerm2[j]=clin3$disOntTerm2[j]
table(clin3$disOntTerm[j],clin32$disOntTerm[j])
if (F) {
    tbl=read.table(paste("docs/Colon and pancreas neuroendocrine passing FGC review.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    names(tbl)[match(c("ID","Gender","Age","Sample.Site","Path.Report.Diagnosis","Disease.Ontology","Jeff.s.Grade","Fiona.Column.1","Fiona.Column.2"),names(tbl))]=c("id","gender","age","Sample.Site","Path.Report.Diagnosis","disOntTerm","gradeJeff","fionaCol1","fionaCol2")
    clin32$cellSize=NA
    j=match(clin3$id,tbl$id); j1=which(!is.na(j)); j2=j[j1]
    clin32$cellSize[j1]=tbl$fionaCol2[j2]
}
tbl=read.table(paste("docs/Sent to Ritu 5-11 JUST Y OR BORDERLINE sample set with LC SC designation as reviewed by Fiona.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
names(tbl)[match(c("ID","Gender","Age","Sample.Site","Disease.Ontology","Fiona.Column.1","CYTO"),names(tbl))]=c("id","gender","age","Sample.Site","disOntTerm","fionaCol1","cyto")
clin32$cellSize=NA
j=match(clin3$id,tbl$id); j1=which(!is.na(j)); j2=j[j1]
clin32$cellSize[j1]=tbl$cyto[j2]
table(disOntTerm=clin3$disOntTerm2,cyto=clin32$cellSize,fionaCol1=clin32$fionaCol1,exclude=NULL)
table(disOntTerm2=clin3$disOntTerm2,cyto=clin32$cellSize,gradeJeff=clin3$grade)
x=table(disOntTerm2=clin3$disOntTerm2,cellSize=clin32$cellSize)
x
x[which(x[,"LC"]!=0),]
tbl=read.table(paste("docs/Sent to Ritu 5-11 JUST Y OR BORDERLINE sample set with TMB MSI.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
names(tbl)[match(c("ID","Gender","Age","Sample.Site","Disease.Ontology","Fiona.Column.1","MSI.Status","TMB..mutations.per.megabase.","TMB.Level"),names(tbl))]=c("id","gender","age","Sample.Site","disOntTerm","fionaCol1","msiStatus","tmbMut","tmbLevel")
j=match(clin3$id,tbl$id); j1=which(!is.na(j)); j2=j[j1]
clin32$msiStatus=clin32$tmbMut=clin32$tmbLevel=NA
clin32$msiStatus[j1]=tbl$msiStatus[j2]
clin32$tmbMut[j1]=as.numeric(tbl$tmbMut[j2])
clin32$tmbLevel[j1]=tbl$tmbLevel[j2]


clin1$disOntTerm2=clin1$disOntTerm
clin1$disOntTerm2[grep("Lung large",clin1$disOntTerm)]="lungLarge"
clin1$disOntTerm2[grep("Lung small",clin1$disOntTerm)]="lungSmall"
clin1$disOntTerm2[grep("Soft tissue primitive",clin1$disOntTerm)]="softTissuePrimitive"
clin1$disOntTerm2[grep("Soft tissue small cell",clin1$disOntTerm)]="softTissueSmall"
clin1$disOntTerm2[grep("Pancreas neuroendocrine",clin1$disOntTerm)]="Pancreas neuroendocrine"
#clin1$disOntTerm2[which(!clin1$disOntTerm2%in%c("Colon neuroendocrine carcinoma","Lung small cell undifferentiated carcinoma","Pancreas neuroendocrine"))]="other"
clin1$disOntTerm2[which(clin1$disOntTerm2%in%c("Duodenum neuroendocrine tumor","Esophagus neuroendocrine carcinoma","Small intestine neuroendocrine carcinoma","Stomach neuroendocrine carcinoma"))]="other"
clin1$disOntTerm2=sapply(clin1$disOntTerm2,function(x) {
    y=strsplit(x," ")[[1]][1]
    y=strsplit(y,"")[[1]]
    y[1]=tolower(y[1])
    y=paste(y,collapse="")
    #y=tolower(y)
    #y=sub("lung","lungSmall",y)
    y
},USE.NAMES=F
)
clin1$disOntTerm3=clin1$disOntTerm2
clin1$disOntTerm3[which(clin1$disOntTerm2%in%c("colon","other","pancreas"))]="combined"

j=match(clin3$id,clin1$id)
clin3$disOntTerm2=clin1$disOntTerm2[j]
clin3$disOntTerm3=clin1$disOntTerm3[j]
#j1=which(!is.na(clin32$disOntTerm))
#j=match(clin32$id[j1],clin1$id); j1=j1[which(!is.na(j))]; j2=j[which(!is.na(j))]
#clin32$disOntTerm2[j1]=clin1$disOntTerm2[j2]
#clin32$disOntTerm3[j1]=clin1$disOntTerm3[j2]

gene1=unique(clin1$gene[clin1$assayVersion=="T5a"])
gene2=unique(clin1$gene[clin1$assayVersion=="T7"])
samId01=which(clin1$gene%in%gene1[gene1%in%gene2])
samId01=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7"))
samId01=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7") & (clin1$grade=="High" | clin1$disOntTerm=="Lung small cell undifferentiated carcinoma"))
switch(subsetFlag,
    "_allAssays"={
        samId01=1:nrow(clin1)
    },
    "_fiona102"={
        j=!is.na(clin32$disOntTerm)
        samId01=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7") & ((clin1$id%in%clin32$id[j]) | (clin1$grade=="High" & !clin1$disOntTerm%in%clin3$disOntTerm[j]) | (clin1$disOntTerm=="Lung small cell undifferentiated carcinoma")))
        samId01=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7") & ((clin1$id%in%clin32$id[j]) | (clin1$grade=="High" & !clin1$disOntTerm2%in%clin3$disOntTerm2[j]) | (clin1$disOntTerm=="Lung small cell undifferentiated carcinoma")))
    },
    "_fiona102_T5aT7Assays"={
        j=!is.na(clin32$disOntTerm)
        samId01=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7") & ((clin1$id%in%clin32$id[j]) | (!clin1$disOntTerm2%in%clin3$disOntTerm2[j])))
    },
    "_fiona166_T5aT7Assays"={
        j=which(clin32$fionaCol1=="Y")
        samId01=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7") & ((clin1$id%in%clin32$id[j]) | (!clin1$disOntTerm2%in%clin3$disOntTerm2[j])))
    }
)
samId03=which(clin3$id%in%clin1$id[samId01])

#samId012=samId011[samId]
#samId013=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7") & (clin1$id%in%clin32$id))


x=table(clin3$disOntTerm[samId03],clin3$assayVersion[samId03])

res=table(clin1$disOntTerm[samId01])
samId11=samId01[which(clin1$disOntTerm[samId01]%in%names(res)[res>=20])]

res=table(clin3$disOntTerm[samId03])
samId13=samId01[which(clin1$disOntTerm[samId01]%in%names(res)[res>=20])]

samId=samId13
samId=samId01

## -----------

if (F) {
    disInfo=sort(unique(clin1$disOntTerm2[samId]))
    disInfo=cbind(disInfo,sapply(disInfo,function(x) {
        y=strsplit(x," ")[[1]][1]
        y=sub("Lung","LungSmall",y)
        y
        },USE.NAMES=F
        ))
    colnames(disInfo)=c("disease","disShort")
    disInfo=as.data.frame(disInfo,stringsAsFactors=F)
}

grpUniq=sort(unique(clin1$disOntTerm2[samId]))
#grpRef="Lung small cell undifferentiated carcinoma"
#grpUniq=matrix(c(rep(grpRef,length(grpUniq)-1),grpUniq[!grpUniq%in%grpRef]),ncol=2,byrow=F)
grpUniq=c("lungSmall","colon","pancreas","other")
n=length(grpUniq)*(length(grpUniq)-1)/2
out=matrix(nrow=n,ncol=2)
k=1
for (k1 in 1:(length(grpUniq)-1)) {
    for (k2 in (k1+1):length(grpUniq)) {
        out[k,1]=grpUniq[k1]
        out[k,2]=grpUniq[k2]
        k=k+1
    }
}
grpUniq=out
grpUniq1=grpUniq
grpUniq2=matrix(c("lungSmall","combined"),ncol=2)

grpUniq=grpUniq1


## ----------------------------------------------

switch(subsetFlag,
    "_allAssays"={
        samIdThis=1:nrow(clin3)
        heading1="All Assays"
    },
    "_fiona166_T5aT7Assays"={
        samIdThis=samId03
        heading1=paste("Fiona certified, ",paste(sort(unique(clin3$assayVersion[samIdThis])),collapse=" & ")," assay samples",sep="")
    }
)
j=samIdThis[which(clin32$tmbLevel[samIdThis]=="High")]

## Boxplots
outFormat="png"
outFormat="pdf"
switch(subsetFlag,
"_allAssays"={
    samIdThis=1:nrow(clin3)
    heading1="All Assays"
},
"_fiona166_T5aT7Assays"={
    samIdThis=samId03
    heading1=paste("Fiona certified, ",paste(sort(unique(clin3$assayVersion[samIdThis])),collapse=" & ")," assay samples",sep="")
}
)
switch(outFormat,
"png"={png(paste("boxplot_tmbMut",subsetFlag,".png",sep=""))},
"pdf"={pdf(paste("boxplot_tmbMut",subsetFlag,".pdf",sep=""))}
)
heading=heading1
x=3-as.integer(as.factor(clin32$tmbLevel[samIdThis]))
res=table(x)
names(res)=clin32$tmbLevel[samIdThis][match(names(res),x)]
boxplot(clin32$tmbMut[samIdThis]~x,names=paste(names(res)," (",res,")",sep=""),main=heading,xlab="TMB level",ylab="TMB mutations per MB")
dev.off()
if (subsetFlag=="_allAssays") {
    samIdThis=1:nrow(clin3)
} else {
    samIdThis=samId03
}
if (any(clin32$tmbMut[samIdThis]>500,na.rm=T)) {
    samIdThis=samIdThis[which(clin32$tmbMut[samIdThis]<500)]
    switch(outFormat,
    "png"={png(paste("boxplot_tmbMut",subsetFlag,"_zoomed.png",sep=""))},
    "pdf"={pdf(paste("boxplot_tmbMut",subsetFlag,"_zoomed.pdf",sep=""))}
    )
    heading=paste(heading1,"\nExcluding TMB outlier",sep="")
    x=3-as.integer(as.factor(clin32$tmbLevel[samIdThis]))
    res=table(x)
    names(res)=clin32$tmbLevel[samIdThis][match(names(res),x)]
    boxplot(clin32$tmbMut[samIdThis]~x,names=paste(names(res)," (",res,")",sep=""),main=heading,xlab="TMB level",ylab="TMB mutations per MB")
    dev.off()
}

if (subsetFlag=="_fiona102_T5aT7Assays") {
    tbl=read.table(paste("docs/Colon and pancreas neuroendocrine passing FGC review.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    names(tbl)[match(c("ID","Gender","Age","Sample.Site","Path.Report.Diagnosis","Disease.Ontology","Jeff.s.Grade","Fiona.Column.1","Fiona.Column.2"),names(tbl))]=c("id","gender","age","Sample.Site","Path.Report.Diagnosis","disOntTerm","gradeJeff","fionaCol1","fionaCol2")
    kk=ceiling(nrow(clin3)/nrow(tbl))
    clin32=NULL
    for (k in 1:kk) clin32=rbind(clin32,tbl)
    clin32=clin32[1:nrow(clin3),]
    for (k in 1:ncol(clin32)) clin32[,k]=NA
    clin32$id=clin3$id
    clin32[match(tbl$id,clin32$id),]=tbl
    j=which(!is.na(clin32$disOntTerm))
    #clin32$disOntTerm2[j]=clin3$disOntTerm2[j]
    #clin32$disOntTerm2[j]=clin3$disOntTerm2[j]
    table(clin3$grade[j],clin32$gradeJeff[j])
    table(clin3$disOntTerm[j],clin32$disOntTerm[j])
}

## ----------------------------------------------
"
* Globally
** proportion of each of the 4 alt type
** Exclude any disease with number >= 20
** enumerate for each disease
** lung small vs. lung large
** lung small vs. other diseases (pair-wise)
"

subsetList=c("_allAssays","")
subsetList="_fiona102"
subsetList="_fiona102_T5aT7Assays"
subsetList="_fiona340_T5aT7Assays"
subsetList="_fiona166_T5aT7Assays"
subsetList="_allAssays"
subsetList=subsetFlag
for (subsetFlag in subsetList) {
    samId=samId01
    heading=paste("High grade & small lung cell, ",paste(sort(unique(clin1$assayVersion[samId01])),collapse=" & ")," assay versions",sep="")
    switch(subsetFlag,
    "_allAssays"={
        heading="All Assays"
        samId=1:nrow(clin1)
    },
    "_t5at7"={
        samId=samId01
    },
    "_fiona102_T5aT7Assays"={
        heading=paste("Fiona certified colon & pancreas + rest, ",paste(sort(unique(clin1$assayVersion[samId01])),collapse=" & ")," assay versions",sep="")
        samId=samId01
    },
    "_fiona340_T5aT7Assays"={
        heading=paste("Fiona checked colon, pancreas & otherGI + rest, ",paste(sort(unique(clin1$assayVersion[samId01])),collapse=" & ")," assay versions",sep="")
        samId=samId01
    },
    "_fiona166_T5aT7Assays"={
        heading=paste("Fiona certified colon, pancreas & otherGI + SCLC, ",paste(sort(unique(clin1$assayVersion[samId01])),collapse=" & ")," assay versions",sep="")
        samId=samId01
    }
    )

    for (disUniq in c("lungSmall+colon+pancreas+other","lungSmall","colon","pancreas","other")) {

        fName=paste("countSummary_",ifelse(disUniq=="","global",disUniq),subsetFlag,".txt",sep="")
        if (disUniq=="") {
            samIdThis=samId
        } else if (disUniq=="lungSmall+colon+pancreas+other") {
            samIdThis=samId[which(clin1$disOntTerm2[samId]%in%c("lungSmall","colon","pancreas","other"))]
        } else {
            samIdThis=samId[which(clin1$disOntTerm2[samId]==disUniq)]
        }

        write.table(heading,file=fName,append=F,col.names=F,row.names=F,sep="\t",quote=F)
        write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
        write.table(toupper(ifelse(disUniq=="","global",disUniq)),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
        write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)

        if (F) {
            res=table(altType=clin1$altType[samIdThis])
            clin4=as.matrix(res)
            res=as.matrix(table(anyAlt=!is.na(clin1$altType[samIdThis])))
            clin4=rbind(noAlt=sum(is.na(clin1$altType[samIdThis])),anyAlt=sum(!duplicated(clin1$id[samIdThis][which(!is.na(clin1$altType[samIdThis]))])),clin4)
            tbl=cbind(rownames(clin4),clin4)
            colnames(tbl)=c("altType","count")
            k=c("noAlt","anyAlt")
            clin4[k,]/sum(clin4[k,],na.rm=T)
            k=which(!colnames(clin4)%in%c("noAlt","anyAlt"))
            clin4[k,]/sum(clin4[k,],na.rm=T)
            nm=c(colnames(tbl),"proportion")
            tbl=cbind(tbl,proportion=clin4)
            colnames(tbl)=nm
            write.table(tbl,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
            write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
        }

        if (F) {
            n=11
            tmp=rep(NA,n)
            tmpC=rep("",n)
            tbl4=data.frame(description=tmpC,count=tmp,proportion=tmp,stringsAsFactors=F)
        }
        j=which(clin3$id%in%clin1$id[samIdThis])
        clin=cbind(clin3,clin32[,which(!names(clin32)%in%names(clin3))])
        names(clin)[match("grade",names(clin))]="gradeJeff"
        clin$ageBi=as.integer(clin$age<=median(clin$age,na.rm=T))
        res=table(clin1$id[samIdThis],clin1$altType[samIdThis])
        res2=apply(res,1,function(x) {sum(x,na.rm=T)})
        summary(res2)
        write.table(paste("Average of ",round(mean(res2),2)," alterations (range ",min(res2,na.rm=T),", ",max(res2,na.rm=T),") per patient",sep=""),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
        write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
        res[res>0]=1
        altTypeUniq1=sort(unique(clin1$altType))
        tmp=matrix(nrow=nrow(clin),ncol=length(altTypeUniq1)+1)
        colnames(tmp)=c("anyAlt",altTypeUniq1)
        clin=cbind(clin,tmp)
        jj=match(rownames(res),clin$id)
        clin$anyAlt[jj]=apply(res,1,function(x) as.integer(any(x==1,na.rm=T)))
        clin[jj,altTypeUniq1]=res
        jj=samIdThis[which(!duplicated(paste(clin1$gene,clin1$varStatus)[samIdThis]))]
        clin4=clin1
        write.table(paste("Median age",round(median(clin$age,na.rm=T),2),sep="\t"),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
        write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
        varList=c("gender","fionaCol1","msiStatus","tmbLevel","gradeJeff","disOntTerm","disOntTerm2","assayVersion","anyAlt",sort(unique(clin1$altType)))
        for (v1 in 1:length(varList)) {
            #res=table(clin[j,varList[v1]])
            res=table(clin[j,varList[v1]],exclude=NULL)
            if (all(is.na(names(res)))) next
            if (res[length(res)]==0) res=res[-length(res)]
            tbl=cbind(names(res),res,round(res/sum(res),2))
            colnames(tbl)=c(varList[v1],"count","proportion")
            if (F) {
            if ("0"%in%tbl[,1]) {
                tbl=tbl[2,]
                tbl[1]=""
                tbl=matrix(tbl,nrow=1)
            }
            }
            if ("0"%in%tbl[,1]) {
                tbl[,1][which(tbl[,1]=="0")]="no mutation"
                tbl[,1][which(tbl[,1]=="1")]="mutation"
            }
            write.table("\tPatient",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
            write.table(tbl,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
            write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
        }
        varList=data.frame(id1=c("disOntTerm2","disOntTerm2","fionaCol1","fionaCol1","msiStatus"),id2=c("fionaCol1","cellSize","msiStatus","tmbLevel","tmbLevel"),stringsAsFactors=F)
        for (v1 in 1:nrow(varList)) {
            res=table(clin[j,varList$id1[v1]],clin[j,varList$id2[v1]],exclude=NULL)
            if (all(is.na(colnames(res))) || all(is.na(rownames(res)))) next
            if (all(res[,ncol(res)]==0)) res=res[,-ncol(res)]
            if (!is.null(dim(res)) && all(res[nrow(res),]==0)) res=res[-nrow(res),]
            if (is.null(dim(res))) next
            tbl=cbind(rownames(res),res)
            colnames(tbl)[1]=varList$id1[v1]
            if ("0"%in%tbl[,1]) {
                tbl[,1][which(tbl[,1]=="0")]="no mutation"
                tbl[,1][which(tbl[,1]=="1")]="mutation"
            }
            if ("0"%in%colnames(tbl[,1])) {
                colnames(tbl)[which(colnames(tbl)=="0")]="no mutation"
                colnames(tbl)[which(colnames(tbl)=="1")]="mutation"
            }
            write.table("\tPatient count",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
            write.table(paste("\t",varList$id2[v1],sep=""),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
            write.table(tbl,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
            write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
        }
    
        if (disUniq=="") {
            if (F) {
                grpRef="lungSmall"
                res=table(disOntTerm2=clin1$disOntTerm2[samId],altType=clin1$altType[samId])
                clin4=as.matrix(res)
                res=table(disOntTerm2=clin1$disOntTerm2[samId],anyAlt=!is.na(clin1$altType[samId]))
                clin4=cbind(noAlt=res[,1],anyAlt=res[,2],clin4)
                clin4=cbind(disOntTerm2=rownames(clin4),as.data.frame(clin4))
                write.table(clin4,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
                write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)

                x=table(clin1[samId,"id"],as.integer(!is.na(clin1$altType[samId])))
                clin5=as.data.frame(cbind(id=as.integer(rownames(x)),noAlt=x[,1],anyAlt=x[,2]))
                clin5=clin5[which(clin5$noAlt!=0 | clin5$anyAlt!=0),]
                clin5$disOntTerm2=clin3$disOntTerm2[match(clin5$id,clin3$id)]
                tmp=rep(NA,nrow(grpUniq))
                out=data.frame(disease=grpUniq[,2],noAlt=tmp,anyAlt=tmp,pvFisher=tmp,coefT=tmp,pvT=tmp,pvWilcox=tmp,stringsAsFactors=F)
                for (k in 1:nrow(out)) {
                    j=which(clin4$disOntTerm2%in%grpUniq[k,2])
                    out[k,c("noAlt","anyAlt")]=clin4[j,c("noAlt","anyAlt")]
                    j=which(clin4$disOntTerm2%in%grpUniq[k,])
                    out$pvFisher[k]=fisher.test(as.matrix(clin4[j,c("noAlt","anyAlt")]))$p.value
                    j=which(clin5$disOntTerm2%in%grpUniq[k,])
                    res=summary(lm(clin5$anyAlt[j]~clin5$disOntTerm2[j]))$coef
                    out[k,c("coefT","pvT")]=res[2,c("Estimate","Pr(>|t|)")]
                    out$pvWilcox[k]=wilcox.test(clin5$anyAlt[j]~clin5$disOntTerm2[j])$p.value
                }
                write.table(paste("Each disease vs. ",grpRef, " (pair-wise) and any vs. no alternation"),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                j=which(clin5$disOntTerm2==grpRef)
                write.table(paste(grpRef,": noAlt = ",sum(clin5$anyAlt[j]==0,na.rm=T), ", anyAlt = ",sum(clin5$anyAlt[j],na.rm=T),sep=""),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                write.table(out,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
                write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
            }
        }
    }
}

## ----------------------------------------------
"
* Gene-level (any alteration)
** proportion of each of the 4 alt type
** Exclude any disease with number >= 20
** enumerate for each disease
** lung small vs. lung large
** lung small vs. other diseases (pair-wise)
"

subsetList=c("_allAssays","")
subsetList="_fiona102"
subsetList="_fiona102_T5aT7Assays"
subsetList="_fiona166_T5aT7Assays"
subsetList="_allAssays"
subsetList=subsetFlag
for (subsetFlag in subsetList) {
    samId=samId13
    samId=samId01
    heading=paste("High grade & small lung cell, ",paste(sort(unique(clin1$assayVersion[samId01])),collapse=" & ")," assay versions",sep="")
    switch(subsetFlag,
    "_allAssays"={
        heading="All assays"
        samId=1:nrow(clin1)
    },
    "_t5at7"={
        samId=samId01
    },
    "_fiona102_T5aT7Assays"={
        heading=paste("Fiona certified colon & pancreas + rest, ",paste(sort(unique(clin1$assayVersion[samId01])),collapse=" & ")," assay versions",sep="")
        samId=samId01
    },
    "_fiona166_T5aT7Assays"={
        heading=paste("Fiona certified colon, pancreas & otherGI + rest, ",paste(sort(unique(clin1$assayVersion[samId01])),collapse=" & ")," assay versions",sep="")
        samId=samId01
    }
    )

    gene=table(clin1$gene[samId])
    gene[order(gene,decreasing=T)][1:10]
    geneList=order(gene,decreasing=T); geneList=geneList[which(gene[geneList]>9)]
    geneList=order(gene,decreasing=T)

    altTypeUniq=c("noAlt","anyAlt",sort(unique(clin1$altType[samId])))
    altTypeUniq1=sort(unique(clin1$altType))
    disOntTerm2Uniq=sort(unique(clin1$disOntTerm2[samId]))

    id=unique(clin1$id[samId])

    ## ---------------------

    fName=paste("countSummary_gene",subsetFlag,".txt",sep="")

    write.table("Gene-level\n",file=fName,append=F,col.names=F,row.names=F,sep="\t",quote=F)

    if (F) {
        out=matrix(0,nrow=length(geneList),ncol=length(altTypeUniq),dimnames=list(names(gene)[geneList],altTypeUniq))
        for (g1 in geneList) {
            geneThis=names(gene)[g1]
            k=which(rownames(out)==geneThis)
            j=samId[which(clin1$gene[samId]==geneThis)]
            res=table(clin1$altType[j])
            out[k,match(names(res),colnames(out))]=res
            out[k,"noAlt"]=sum(id%in%clin1$id[samId] & !id%in%clin1$id[j])
        }
        out[,"anyAlt"]=apply(out[,sort(unique(clin1$altType[samId]))],1,sum,na.rm=T)
        out=data.frame(gene=names(gene)[geneList],out,stringsAsFactors=F)
        write.table(out,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
        write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
    }

    x1=sum(!duplicated(clin1$id))
    out=matrix(0,nrow=length(geneList),ncol=length(altTypeUniq1)+2,dimnames=list(names(gene)[geneList],c("anyAltCount","anyAlt",altTypeUniq1)))
    for (g1 in geneList) {
        geneThis=names(gene)[g1]
        k=which(rownames(out)==geneThis)
        j=samId[which(clin1$gene[samId]==geneThis)]
        res=table(clin1$id[j],clin1$altType[j])
        res[res>0]=1
        res2=apply(res,2,sum,na.rm=T)
        out[k,match(names(res2),colnames(out))]=res2/x1
        out[k,"anyAltCount"]=sum(apply(res,1,function(x) as.integer(any(x==1,na.rm=T))))
        out[k,"anyAlt"]=out[k,"anyAltCount"]/x1
    }
    grp=rep("",length(geneList))
    names(grp)=names(gene)[geneList]
    j=which(names(grp)%in%clin1$gene[which(clin1$varStatus=="known")])
    grp[j]="known"
    j=which(names(grp)%in%clin1$gene[which(clin1$varStatus=="likely")])
    grp[j]=paste(grp[j],"likely",sep="/")
    grp[which(grp=="/likely")]="likely"
    out=data.frame(gene=names(gene)[geneList],varStatus=grp,round(out,3),stringsAsFactors=F)
    write.table(paste("","","Number of patients with","Proportion of patients with",sep="\t"),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
    nm=colnames(out)
    nm[nm=="anyAltCount"]="anyAlt"
    names(out)=nm
    write.table(out,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
    write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)

    ## ---------------------

    fName=paste("countSummary_gene_disease",subsetFlag,".txt",sep="")

    write.table("Gene-disease-level\n",file=fName,append=F,col.names=F,row.names=F,sep="\t",quote=F)

    out=matrix(0,nrow=length(geneList),ncol=length(altTypeUniq)*length(disOntTerm2Uniq),dimnames=list(names(gene)[geneList],paste(altTypeUniq,rep(disOntTerm2Uniq,each=length(altTypeUniq)),sep="_")))
    for (g1 in geneList) {
        geneThis=names(gene)[g1]
        k=which(rownames(out)==geneThis)
        j=samId[which(clin1$gene[samId]==geneThis)]
        res=table(disOntTerm2=clin1$disOntTerm2[j],altType=clin1$altType[j])
        out[k,match(paste(rep(colnames(res),each=nrow(res)),rownames(res),sep="_"),colnames(out))]=c(res)
        out[k,match(paste("anyAlt",rownames(res),sep="_"),colnames(out))]=apply(res,1,sum,na.rm=T)
        id2=id[id%in%clin1$id[samId] & !id%in%clin1$id[j]]
        res=table(clin1$disOntTerm2[samId][clin1$id[samId]%in%id2])
        out[k,match(paste("noAlt",names(res),sep="_"),colnames(out))]=res
    }
    #tbl=data.frame(gene=names(gene)[geneList],out,stringsAsFactors=F)
    #write.table(tbl,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
    #write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
    ann4=out


    out=matrix(0,nrow=length(geneList),ncol=length(disOntTerm2Uniq)*(length(altTypeUniq1)+2),dimnames=list(names(gene)[geneList],paste(c("anyAltCount","anyAlt",altTypeUniq1),rep(disOntTerm2Uniq,each=length(altTypeUniq1)+2),sep="_")))
    for (d1 in 1:length(disOntTerm2Uniq)) {
        jj=samId[which(clin1$disOntTerm2[samId]==disOntTerm2Uniq[d1])]
        x1=sum(!duplicated(clin1$id[jj]))
        for (g1 in geneList) {
            geneThis=names(gene)[g1]
            k=which(rownames(out)==geneThis)
            j=jj[which(clin1$gene[jj]==geneThis)]
            if (length(j)!=0) {
                res=table(clin1$id[j],clin1$altType[j])
                res[res>0]=1
                res2=apply(res,2,sum,na.rm=T)
                out[k,match(paste(names(res2),disOntTerm2Uniq[d1],sep="_"),colnames(out))]=res2/x1
                out[k,paste("anyAltCount",disOntTerm2Uniq[d1],sep="_")]=sum(apply(res,1,function(x) as.integer(any(x==1,na.rm=T))))
                out[k,paste("anyAlt",disOntTerm2Uniq[d1],sep="_")]=out[k,paste("anyAltCount",disOntTerm2Uniq[d1],sep="_")]/x1
            }
        }
    }
    datGD=out
    grp=rep("",length(geneList))
    names(grp)=names(gene)[geneList]
    j=which(names(grp)%in%clin1$gene[which(clin1$varStatus=="known")])
    grp[j]="known"
    j=which(names(grp)%in%clin1$gene[which(clin1$varStatus=="likely")])
    grp[j]=paste(grp[j],"likely",sep="/")
    grp[which(grp=="/likely")]="likely"
    out=data.frame(gene=names(gene)[geneList],varStatus=grp,round(out,3),stringsAsFactors=F)
    write.table(paste(c("","",rep(c("Number of patients with","Proportion of patients with",rep("",length(altTypeUniq1))),length(disOntTerm2Uniq))),collapse="\t"),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
    nm=colnames(out)
    nm=sub("anyAltCount_","anyAlt_",nm)
    names(out)=nm
    write.table(out,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
    write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)

    ## ---------------------
    ## Gene-patient-level

    out=matrix(0,nrow=length(geneList),ncol=sum(!duplicated(clin1$id[samId]))*(length(altTypeUniq1)+1),dimnames=list(names(gene)[geneList],paste(c("anyAlt",altTypeUniq1),rep(clin1$id[samId][!duplicated(clin1$id[samId])],each=length(altTypeUniq1)+1),sep="_")))
    id=sapply(colnames(out),function(x){strsplit(x,"_")[[1]][2]},USE.NAMES=F)
    for (a1 in 1:length(altTypeUniq1)) {
        jj=samId[which(clin1$altType[samId]==altTypeUniq1[a1])]
        jj2=which(sapply(colnames(out),function(x){strsplit(x,"_")[[1]][1]},USE.NAMES=F)==altTypeUniq1[a1])
        for (g1 in geneList) {
            geneThis=names(gene)[g1]
            j=jj[which(clin1$gene[jj]==geneThis)]
            if (length(j)!=0) {
                j2=jj2[match(clin1$id[j],id[jj2])]
                k=which(rownames(out)==geneThis)
                out[k,j2]=1
            }
        }
    }
    j1=which(sapply(colnames(out),function(x){strsplit(x,"_")[[1]][1]},USE.NAMES=F)=="anyAlt")
    for (a1 in 1:length(altTypeUniq1)) {
        j2=which(sapply(colnames(out),function(x){strsplit(x,"_")[[1]][1]},USE.NAMES=F)==altTypeUniq1[a1])
        out[,j1]=out[,j1]+out[,j2]
    }
    out[out>0]=1
    datGP=out
    
    ## ---------------------
    ## Feature-patient-level
    ## Feature could be multiple for a gene if there are multiple entries for a gene for a patient
    ## For multiple entries, we select the one with highest fraction read
    
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

    ## ---------------------

    if (F) {
    fName=paste("pValue_gene_disease",subsetFlag,".txt",sep="")

    write.table("Gene-disease-level\n",file=fName,append=F,col.names=F,row.names=F,sep="\t",quote=F)

    tmp=rep(NA,nrow(grpUniq))
    out=matrix(nrow=nrow(ann4),ncol=2*nrow(grpUniq),dimnames=list(rownames(ann4),paste("pvFisher_",c("anyOrNoAlt","anyOf4Alt"),"_",rep(grpUniq[,2],each=2),sep="")))
    for (g1 in 1:nrow(out)) {
        for (k in 1:nrow(grpUniq)) {
            x=altTypeUniq[altTypeUniq%in%c("noAlt","anyAlt")]
            j1=which(colnames(ann4)%in%paste(x,grpUniq[k,1],sep="_"))
            j2=which(colnames(ann4)%in%paste(x,grpUniq[k,2],sep="_"))
            x=cbind(ann4[g1,j1],ann4[g1,j2])
            if (sum(c(x))>4) {
                out[g1,paste("pvFisher_anyOrNoAlt_",grpUniq[k,2],sep="")]=fisher.test(x)$p.value
            }
            x=altTypeUniq[!altTypeUniq%in%c("noAlt","anyAlt")]
            j1=which(colnames(ann4)%in%paste(x,grpUniq[k,1],sep="_"))
            j2=which(colnames(ann4)%in%paste(x,grpUniq[k,2],sep="_"))
            x=cbind(ann4[g1,j1],ann4[g1,j2])
            if (sum(c(x))>4) {
                out[g1,paste("pvFisher_anyOf4Alt_",grpUniq[k,2],sep="")]=fisher.test(x)$p.value
            }
        }
    }
    tbl=cbind(gene=rownames(out),out)
    write.table(paste("Each disease vs. ",grpRef, " (pair-wise) and 1) any vs. no alternation, 2) ",paste(altTypeUniq[!altTypeUniq%in%c("noAlt","anyAlt")],collapse=" vs. ")),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
    write.table(tbl,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
    write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
    }

    ## ---------------------

    if (F) {
    fName=paste("countSummary_gene_disease_cdsEffect",subsetFlag,".txt",sep="")

    #feat=paste(clin1$disease,clin1$gene,clin1$cdsEffect,sep=" /// ")[samId]
    feat=paste(clin1$gene,clin1$cdsEffect,sep=" /// ")[samId]
    #featUniq=sort(unique(feat))
    featUniq=table(feat)
    featUniq=names(featUniq)[featUniq>1]
    featUniq=featUniq[-grep(" /// NA",featUniq)]

    write.table("Gene-disease-cdsEffect-level\n",file=fName,append=F,col.names=F,row.names=F,sep="\t",quote=F)

    out=matrix(0,nrow=length(featUniq),ncol=length(altTypeUniq)*length(disOntTerm2Uniq),dimnames=list(featUniq,paste(altTypeUniq,rep(disOntTerm2Uniq,each=length(altTypeUniq)),sep="_")))
    xxx=c()
    for (k in 1:length(featUniq)) {
        geneThis=featUniq[k]
        k=which(rownames(out)==geneThis)
        j=samId[which(feat==geneThis)]
        if (length(j)>1) {
            xxx=c(xxx,k)
            res=table(disOntTerm2=clin1$disOntTerm2[j],altType=clin1$altType[j])
            out[k,match(paste(rep(colnames(res),each=nrow(res)),rownames(res),sep="_"),colnames(out))]=c(res)
            out[k,match(paste("anyAlt",rownames(res),sep="_"),colnames(out))]=apply(res,1,sum,na.rm=T)
            id2=id[id%in%clin1$id[samId] & !id%in%clin1$id[j]]
            res=table(clin1$disOntTerm2[samId][clin1$id[samId]%in%id2])
            out[k,match(paste("noAlt",names(res),sep="_"),colnames(out))]=res
        }
    }
    #tbl=out[which(apply(out,1,sum,na.rm=T)>1),]
    tbl=out[,-grep("noAlt_",colnames(out))]
    tbl=data.frame(feature=rownames(tbl),tbl,stringsAsFactors=F)
    write.table(tbl,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
    write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
    }

    ## ---------------------

    fName=paste("countSummary_patient",subsetFlag,".txt",sep="")

    sam=clin1$id[samId]
    samUniq=sort(unique(paste(clin1$disOntTerm2[samId],sam)))
    samUniq=unique(sam[order(paste(clin1$disOntTerm2[samId],formatC(sam,width=nchar(as.character(max(sam))),flag="0")))])

    fName=paste("alterationCount_patient",subsetFlag,".txt",sep="")
    write.table("Patient-level\n",file=fName,append=F,col.names=F,row.names=F,sep="\t",quote=F)

    out=matrix(0,nrow=length(samUniq),ncol=length(altTypeUniq),dimnames=list(samUniq,altTypeUniq))
    for (k in 1:length(samUniq)) {
        samThis=samUniq[k]
        k=which(rownames(out)==samThis)
        j=samId[which(sam==samThis)]
        res=table(clin1$altType[j])
        out[k,match(names(res),colnames(out))]=c(res)
    }
    out[,"anyAlt"]=apply(out[,sort(unique(clin1$altType))],1,sum,na.rm=T)
    tbl=out[,which(colnames(out)!="noAlt")]
    j=match(rownames(tbl),clin3$id)
    tbl=data.frame(patient=rownames(tbl),disOntTerm2=clin3$disOntTerm2[j],tbl,stringsAsFactors=F)
    write.table(paste("Number of alterations"),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
    write.table(tbl,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
    write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)

    fName=paste("alterationStatus_patient",subsetFlag,".txt",sep="")
    write.table("Patient-level\n",file=fName,append=F,col.names=F,row.names=F,sep="\t",quote=F)

    tbl=data.frame(patient=rownames(tbl),disOntTerm2=clin3$disOntTerm2[match(rownames(tbl),clin3$id)],tbl,stringsAsFactors=F)
    write.table(paste("Alteration status (dichotomous)"),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
    write.table(tbl,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
    write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)

    ann=data.frame(id=rep(samUniq,each=length(geneList)),gene=rep(names(gene)[geneList],length(samUniq)))
    out=matrix(0,nrow=length(samUniq)*length(geneList),ncol=length(altTypeUniq),dimnames=list(paste(ann$id,ann$gene,sep="_"),altTypeUniq))
    for (g1 in geneList) {
        geneThis=names(gene)[g1]
        j=samId[which(clin1$gene[samId]==geneThis)]
        res=table(id=clin1$id[j],altType=clin1$altType[j])
        out[match(paste(rownames(res),geneThis,sep="_"),rownames(out)),match(colnames(res),colnames(out))]=res
    }
    out[,"anyAlt"]=apply(out[,sort(unique(clin1$altType))],1,sum,na.rm=T)
    clin6=cbind(ann,clin3[match(ann$id,clin3$id),c("disOntTerm2","disOntTerm3")],as.data.frame(out,stringsAsFactors=F))
    rownames(clin6)=NULL
    tbl=apply(out,c(1,2),function(x) as.integer(any(x!=0,na.rm=T)))
    clin7=cbind(ann,clin3[match(ann$id,clin3$id),c("disOntTerm2","disOntTerm3")],as.data.frame(tbl,stringsAsFactors=F))
    rownames(clin7)=NULL

    ## ---------------------

    if (F) {
    altTypeUniq2=c("anyAlt",altTypeUniq1)

    for (a1 in 1:length(altTypeUniq2)) {
        out=matrix(nrow=length(geneList),ncol=length(sam),dimnames=list(names(gene)[geneList],sam))
        for (g1 in geneList) {
            geneThis=names(gene)[g1]
            k=which(rownames(out)==geneThis)
            j=samId[which(clin1$gene[samId]==geneThis)]
            res=table(disOntTerm2=clin1$disOntTerm2[j],altType=clin1$altType[j])
            out[k,match(paste(rep(colnames(res),each=nrow(res)),rownames(res),sep="_"),colnames(out))]=c(res)
            out[k,match(paste("anyAlt",rownames(res),sep="_"),colnames(out))]=apply(res,1,sum,na.rm=T)
            id2=id[id%in%clin1$id[samId] & !id%in%clin1$id[j]]
            res=table(clin1$disOntTerm2[samId][clin1$id[samId]%in%id2])
            out[k,match(paste("noAlt",names(res),sep="_"),colnames(out))]=res
        }
    }
    }

    ## ---------------------

    library(qvalue)

    altTypeUniq2=c("anyAlt",altTypeUniq1)

    fName=paste("count_gene_disease",subsetFlag,".txt",sep="")
    write.table("Gene-disease-level\n",file=fName,append=F,col.names=F,row.names=F,sep="\t",quote=F)
    varId="disOntTerm2"; grpUniq=grpUniq1
    tmp=rep(NA,nrow(grpUniq))
    out=matrix(nrow=nrow(ann4),ncol=length(altTypeUniq2)*nrow(grpUniq),dimnames=list(rownames(ann4),paste("num_",altTypeUniq2,"_",rep(paste(grpUniq[,2],grpUniq[,1],sep="V"),each=length(altTypeUniq2)),sep="")))
    print(nrow(out))
    timeStamp=Sys.time()
    print(format(timeStamp, "%x %X"))
    for (g1 in 1:nrow(out)) {
        geneThis=rownames(out)[g1]
        if (g1%%20==0) print(g1)
        j1=which(clin7$gene==geneThis)
        for (k in 1:nrow(grpUniq)) {
            j=j1[which(clin7[j1,varId]%in%grpUniq[k,])]
            for (a1 in 1:length(altTypeUniq2)) {
                l=which(colnames(clin7)==altTypeUniq2[a1])
                res=sum(clin7[j,l]==1,na.rm=T)
                out[g1,paste("num_",altTypeUniq2[a1],"_",grpUniq[k,2],"V",grpUniq[k,1],sep="")]=res
            }
        }
    }
    summary(c(out))
    "
    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    0.00    0.00    0.00    3.39    2.00  589.00
    "
    tbl=cbind(gene=rownames(out),out)
    write.table(paste("Number of patients"),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
    write.table(tbl,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
    write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)


    if (F) {
        ## Testing ...
        geneThis="ARID2"
        geneThis="MYCL1"
        j1=which(clin1$disOntTerm2%in%c("colon","lungSmall")); j1=j1[j1%in%samId01]
        j1=samId01[which(clin1$disOntTerm2[samId01]%in%c("colon","lungSmall"))]
        j3=which(clin3$id%in%clin1$id[j1])
        id=clin3$id[j3]
        x=rep(0,length(id))
        for (jj in 1:length(id)) {
            jjj=which(clin1$id[j1]==id[jj]  & clin1$altType[j1]=="amplification" & clin1$gene[j1]==geneThis)
            x[jj]=length(jjj)
        }
        y=as.integer(x!=0)
        fisher.test(clin3$disOntTerm2[j3],y)
        wilcox_test(x~as.factor(clin3$disOntTerm2[j3]), distribution="exact", conf.int=F)

        table(clin7$gene==geneThis & clin7$disOntTerm2%in%c("colon","lungSmall"))
    }

    fName=paste("pValue_fisherTest_gene_disease",subsetFlag,".txt",sep="")
    write.table("Gene-disease-level\n",file=fName,append=F,col.names=F,row.names=F,sep="\t",quote=F)
    varId="disOntTerm2"; grpUniq=grpUniq1
    tmp=rep(NA,nrow(grpUniq))
    out=matrix(nrow=nrow(ann4),ncol=length(altTypeUniq2)*nrow(grpUniq),dimnames=list(rownames(ann4),paste("pv_",altTypeUniq2,"_",rep(paste(grpUniq[,2],grpUniq[,1],sep="V"),each=length(altTypeUniq2)),sep="")))
    print(nrow(out))
    timeStamp=Sys.time()
    print(format(timeStamp, "%x %X"))
    for (g1 in 1:nrow(out)) {
    #for (g1 in 1:10) {
        if (g1%%20==0) print(g1)
        geneThis=rownames(out)[g1]
        j1=which(clin7$gene==geneThis)
        for (k in 1:nrow(grpUniq)) {
            j=j1[which(clin7[j1,varId]%in%grpUniq[k,])]
            for (a1 in 1:length(altTypeUniq2)) {
                l=which(colnames(clin7)==altTypeUniq2[a1])
                x=table(clin7[j,varId],clin7[j,l])
                if (nrow(x)==2 && ncol(x)==2) {
                    pv=fisher.test(x)$p.value
                    out[g1,paste("pv_",altTypeUniq2[a1],"_",grpUniq[k,2],"V",grpUniq[k,1],sep="")]=pv
                }
            }
        }
    }
    out1=out
    varId="disOntTerm3"; grpUniq=grpUniq2
    tmp=rep(NA,nrow(grpUniq))
    out=matrix(nrow=nrow(ann4),ncol=length(altTypeUniq2)*nrow(grpUniq),dimnames=list(rownames(ann4),paste("pv_",altTypeUniq2,"_",rep(paste(grpUniq[,2],grpUniq[,1],sep="V"),each=length(altTypeUniq2)),sep="")))
    print(nrow(out))
    timeStamp=Sys.time()
    print(format(timeStamp, "%x %X"))
    for (g1 in 1:nrow(out)) {
    #for (g1 in 1:10) {
        geneThis=rownames(out)[g1]
        if (g1%%20==0) print(g1)
        j1=which(clin7$gene==geneThis)
        for (k in 1:nrow(grpUniq)) {
            j=j1[which(clin7[j1,varId]%in%grpUniq[k,])]
            for (a1 in 1:length(altTypeUniq2)) {
                l=which(colnames(clin7)==altTypeUniq2[a1])
                x=table(clin7[j,varId],clin7[j,l])
                if (nrow(x)==2 && ncol(x)==2) {
                    pv=fisher.test(x)$p.value
                    out[g1,paste("pv_",altTypeUniq2[a1],"_",grpUniq[k,2],"V",grpUniq[k,1],sep="")]=pv
                }
            }
        }
    }
    grpUniq=grpUniq1
    timeStamp=c(timeStamp,Sys.time())
    print(format(timeStamp[2], "%x %X"))
    print(diff(timeStamp))
    out=cbind(out1,out)
    out[which(out>1)]=1
    tbl=cbind(gene=rownames(out),out)
    write.table(paste("Pair-wise comparison of diseases based on alteration status for a patient"),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
    write.table(paste("P-values from Fisher's exact tests"),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
    write.table(tbl,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
    write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)

    out1=matrix(nrow=nrow(out),ncol=ncol(out),dimnames=list(rownames(out),sub("pv_","qv_",colnames(out))))
    for (k in 1:ncol(out)) {
        i=which(!is.na(out[,k]))
        out1[i,k]=qvalue(out[i,k])$qvalues
        
    }
    fName=sub("pValue","qValue",fName)
    tbl=cbind(gene=rownames(out1),out1)
    write.table(paste("Pair-wise comparison of diseases based on alteration status for a patient"),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
    write.table(paste("Q-values (False discovery rate adjusted p-values from Fisher's exact tests"),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
    write.table(tbl,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
    write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)


    #save.image(paste("tmp",subsetFlag,".RData",sep=""))
}

## --------
load(paste("tmp",subsetFlag,".RData",sep=""))

if (T) {
    grp=paste(clin3$assayVersion," / ",clin3$grade," grade",sep="")
    gene1=unique(clin1$gene[clin1$assayVersion=="T5a"])
    gene2=unique(clin1$gene[clin1$assayVersion=="T7"])
    j1=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7") & (clin1$grade=="High" | clin1$disOntTerm=="Lung small cell undifferentiated carcinoma"))
    j3=which(clin3$id%in%clin1$id[j1])
    cat("\n\nAll samples:\n")
    print(table(disOntTerm=clin3$disOntTerm,asayVersion=clin3$assayVersion))
    j=which(clin1$grade=="High" | clin1$disOntTerm=="Lung small cell undifferentiated carcinoma")
    j=which(clin3$id%in%clin1$id[j])
    cat("\n\nHigh grade & SCLC samples:\n")
    print(table(disOntTerm=clin3$disOntTerm[j],assayVersion=clin3$assayVersion[j]))
    #table(disOntTerm=clin3$disOntTerm,asayVersion.grade=grp)
    #j=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7"))
    j=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7") & (clin1$grade=="High" | clin1$disOntTerm=="Lung small cell undifferentiated carcinoma"))
    j=which(clin3$id%in%clin1$id[j])
    #table(disOntTerm=clin3$disOntTerm[j],asayVersion.grade=grp[j])
    #table(disOntTerm=clin3$disOntTerm[j],assayVersion=clin3$assayVersion[j])
    j=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7") & (clin1$grade=="High" | clin1$disOntTerm=="Lung small cell undifferentiated carcinoma"))
    j=which(clin3$id%in%clin1$id[j])
    cat(paste("\n\nHigh grade & SCLC for genes in common between T5a & T7:\n",sep=""))
    print(table(disOntTerm=clin3$disOntTerm[j],assayVersion=clin3$assayVersion[j]))
    #table(disOntTerm=clin3$disOntTerm2[j],assayVersion=clin3$assayVersion[j])
    #table(disOntTerm=clin3$disOntTerm2[j3],assayVersion=clin3$assayVersion[j3])
    #table(disOntTerm=clin3$disOntTerm2[samId03],assayVersion=clin3$assayVersion[samId03])
}

pThres=0.05

x=table(clin3$disOntTerm2[which(clin3$id%in%clin1$id[samId])])
x=x[order(x,decreasing=T)]
disUniq=names(x)
disUniq=c("lungSmall","pancreas","colon","other")
disUniq=cbind(disUniq,toupper(substr(disUniq,1,1)))
disUniq[,2]=sub("L","S",disUniq[,2])

grpUniq3=grpUniq
k=2
grpUniq3=grpUniq3[which(grpUniq3[,1]==disUniq[k] | grpUniq3[,2]==disUniq[k]),]

fName1="_signifPairwise"
fName1="_signifPairwise_12genes"
fName1=""

propFlag=">="
propFlag=">"


for (thisFlag in c("10%patient","15%patient","signif","targetableGene")) {
    switch(thisFlag,
    "signif"={
        fName=paste("proportionTable_geneByDisease_qv",pThres,fName1,subsetFlag,".txt",sep="")
        heading="Table of significant genes in any of the comparisons pancreas vs. SCLC, other vs. SCLC or colon vs. SCLC"
        k1=grep("qv_anyAlt_",colnames(out1))
        k1=k1[grep("lungSmall",colnames(out1)[k1])]
        k1=k1[-grep("combined",colnames(out1)[k1])]
        nm1=colnames(out1)[k1]
        i=c()
        for (k in k1) {
            i=c(i,which(out1[,k]<pThres))
        }
        i=sort(unique(i))
    },
    "targetableGene"={
        fName=paste("proportionTable_geneByDisease_targetGene",fName1,subsetFlag,".txt",sep="")
        heading="Table of targetable genes"
        k1=grep("qv_anyAlt_",colnames(out1))
        k1=k1[grep("lungSmall",colnames(out1)[k1])]
        k1=k1[-grep("combined",colnames(out1)[k1])]
        nm1=colnames(out1)[k1]
        i=match(toupper(candGene$gene),toupper(rownames(out1)))
        i=i[which(!is.na(i))]
    }
    )
    if (length(grep("%patient",thisFlag))==1) {
        propThres=0.10
        propThres=as.integer(sub("%patient","",thisFlag))/100
        i=c()
        for (k in which(substr(colnames(datGD),1,nchar("anyAlt_"))=="anyAlt_" & sub("anyAlt_","",colnames(datGD))%in%disUniq)) {
            if (propFlag==">") {
                #if (grep("_12genes",fName1)) {
                #i=c(i,which(datGD[,k]>propThres))
                i=c(i,which(round(datGD[,k],2)>propThres))
            } else {
                #i=c(i,which(datGD[,k]>=propThres))
                i=c(i,which(round(datGD[,k],2)>=propThres))
            }
        }
        i=sort(unique(i))
        fName=paste("proportionTable_geneByDisease_",propThres*100,"percPatient",fName1,subsetFlag,".txt",sep="")
        heading=paste("Table of ",length(i)," genes with alterations in ",propFlag," ",propThres*100,"% of patients in any group",sep="")
    }

    #round(datGD[i,grep("anyAlt_",colnames(datGD))],2)
    out2=round(datGD[i,grep("anyAlt_",colnames(datGD))],2)
    i=match(rownames(out1),rownames(out2));
    i1=which(!is.na(i)); i2=i[i1]
    out3=apply(out2,c(1,2),function(x) {
        #y=formatC(x*100,width=2,flag="0")
        #y=paste(" 0.",y," ",sep="")
        y=paste(" ",x*100,"% ",sep="")
        y
    })
    for (d1 in 1:nrow(disUniq)) {
        k1=grep("qv_anyAlt_",colnames(out1))
        k1=k1[-grep("combined",colnames(out1)[k1])]
        k1=k1[grep(disUniq[d1,1],colnames(out1)[k1])]
        nm1=colnames(out1)[k1]
        nm11=nm1
        nm1=sub("qv_anyAlt_","",sub(paste("V",disUniq[d1,1],sep=""),"",sub(paste(disUniq[d1,1],"V",sep=""),"",nm1)))
        nm21=colnames(out2)[grep("anyAlt_",colnames(out2))]
        nm2=sub("anyAlt_","",colnames(out2)[grep("anyAlt_",colnames(out2))])
        x=" "
        for (k in 1:length(k1)) {
            k21=which(nm2==nm1[k])
            i=which(out1[i1,nm11[k]]<pThres)
            if (length(i)!=0) {
                out3[i2[i],k21]=paste(out3[i2[i],k21],disUniq[d1,2],sep="")
            }
        }
    }
    if (F) {
        i=match(rownames(out3),rownames(out1)); k="qv_anyAlt_colonVlungSmall"
        i=match(rownames(out3),rownames(out1)); k="qv_anyAlt_pancreasVlungSmall"
        i=match(rownames(out3),rownames(out1)); k="qv_anyAlt_otherVlungSmall"
        table(out1[i,k]<.05)
        signif(out1[i,k],2)
    }
    x=table(clin3$disOntTerm2[which(clin3$id%in%clin1$id[samId])])
    x=x[order(x,decreasing=T)]
    x=x[match(c("lungSmall","pancreas","colon","other"),names(x))]
    nm=names(x)
    colnames(out3)=sub("anyAlt_","",colnames(out3))
    out3=out3[,nm]
    nm=colnames(out3)
    nm[match("lungSmall",nm)]="SCLC"
    nm=paste(c("gene",paste(nm," (N=",x,")",sep="")),collapse="\t")
    tbl=cbind(gene=rownames(out3),data.frame(out3))
    write.table(heading,file=fName,append=F,col.names=F,row.names=F,sep="\t",quote=F)
    write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
    write.table(nm,file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
    write.table(tbl,file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
}

x1=read.table(paste("proportionTable_geneByDisease_targetGene.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=2)
x2=read.table(paste("results/asco/proportionTable_geneByDisease_targetGene.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,skip=2)

if (F) {
k1=grep("qv_anyAlt_",colnames(out1))
k1=k1[grep("lungSmall",colnames(out1)[k1])]
nm1=colnames(out1)[k1]
nm11=nm1
nm1=sub("qv_anyAlt_","",sub("VlungSmall","",nm1))
for (k in 1:length(nm21)) {
    print(nm11[k])
    print(table(out1[,nm11[k]]<pThres))
    i=which(out1[,nm11[k]]<pThres)
    print(rownames(out1)[i])
    
}
}

#SETD2 - gene


if (T) {
    cat("\n\n========================================================\n",sep="")
    cat("Comparison of mutation status (if a patient has at least 1 mutation vs. none)\n",sep="")
    cat("Computed mutation status for each patient and then compared between two diseases for each gene\n",sep="")
    cat("Performed Fisher's exact tests\n",sep="")
    cat("Q-values were computed to adjust for multiple testing\n",sep="")
    cat("Q-value < 0.05 was considered significant\n",sep="")
    cat("========================================================\n",sep="")
    cat("\n\nTotal no. of comparisons: ",ncol(out1)/length(altTypeUniq2),"\n\n",sep="")
    a1=1
    nm=sapply(colnames(out1)[grep(altTypeUniq2[a1],colnames(out1))],function(x) {
        sub("V"," vs. ",strsplit(x,"_")[[1]][3])
    },USE.NAMES=F)
    cat(nm,sep="\n")
    for (a1 in 1:length(altTypeUniq2)) {
        cat("\n\n========== ",altTypeUniq2[a1],"\n",sep="")
        res=apply(out1[,grep(altTypeUniq2[a1],colnames(out1))],2,function(x) sum(x<0.05,na.rm=T))
        nm=sapply(names(res),function(x) {
            sub("V"," vs. ",strsplit(x,"_")[[1]][3])
        },USE.NAMES=F)
        names(res)=nm
        k=which(res!=0)
        cat("No of comparisons with any significant gene: ",length(k),"\n",sep="")
        if (length(k)!=0) {
            print(res[k])
        }
    }
    
    res1=res
    
    cat("\n\n========================================================\n",sep="")
    cat("\n\nTotal no. of genes: ",nrow(out1),"\n\n",sep="")
    for (a1 in 1:length(altTypeUniq2)) {
        cat("\n\n========== ",altTypeUniq2[a1],"\n",sep="")
        res=apply(out1[,grep(altTypeUniq2[a1],colnames(out1))],1,function(x) sum(x<0.05,na.rm=T))
        k=order(res,decreasing=T); k=k[which(res[k]!=0)]
        cat("No of genes significant for any of the comparisons: ",length(k),"\n",sep="")
        if (length(k)!=0) {
            print(res[k])
        }
    }
}

## ---------------------
load(paste("tmp",subsetFlag,".RData",sep=""))

library(coin)

altTypeUniq2=c("anyAlt",altTypeUniq1)

fName=paste("pValue_wilcoxTest_gene_disease",subsetFlag,".txt",sep="")
write.table("Gene-disease-level\n",file=fName,append=F,col.names=F,row.names=F,sep="\t",quote=F)
varId="disOntTerm3"; grpUniq=grpUniq2
tmp=rep(NA,nrow(grpUniq))
out=matrix(nrow=nrow(ann4),ncol=length(altTypeUniq2)*nrow(grpUniq),dimnames=list(rownames(ann4),paste("pv_",altTypeUniq2,"_",rep(paste(grpUniq[,2],grpUniq[,1],sep="V"),each=length(altTypeUniq2)),sep="")))
print(nrow(out))
timeStamp=Sys.time()
print(format(timeStamp, "%x %X"))
#for (g1 in 1:nrow(out)) {
for (g1 in 1:2) {
    geneThis=rownames(out)[g1]
    if (g1%%20==0) print(g1)
    j1=which(clin7$gene==geneThis)
    for (k in 1:nrow(grpUniq)) {
        j=j1[which(clin7[j1,varId]%in%grpUniq[k,])]
        for (a1 in 1:length(altTypeUniq2)) {
            l=which(colnames(clin7)==altTypeUniq2[a1])
            if (any(clin7[j,l]!=0,na.rm=T)) {
                res=try(wilcox_test(clin6[j,l] ~ as.factor(clin6[j,varId]), distribution="exact", conf.int=F))
                if (class(res)!="try-error") {
                    pv=try(pvalue(res))
                    if (class(res)!="try-error") {
                        out[g1,paste("pv_",altTypeUniq2[a1],"_",grpUniq[k,2],"V",grpUniq[k,1],sep="")]=pv
                    }
                }
            }
        }
    }
}

#save.image(paste("tmp2_1",subsetFlag,".RData",sep=""))
out2=out
varId="disOntTerm2"; grpUniq=grpUniq1
tmp=rep(NA,nrow(grpUniq))
out=matrix(nrow=nrow(ann4),ncol=length(altTypeUniq2)*nrow(grpUniq),dimnames=list(rownames(ann4),paste("pv_",altTypeUniq2,"_",rep(paste(grpUniq[,2],grpUniq[,1],sep="V"),each=length(altTypeUniq2)),sep="")))
print(nrow(out))
timeStamp=Sys.time()
print(format(timeStamp, "%x %X"))
#for (g1 in 1:nrow(out)) {
for (g1 in 1:2) {
    geneThis=rownames(out)[g1]
    if (g1%%20==0) print(g1)
    j1=which(clin7$gene==geneThis)
    for (k in 1:nrow(grpUniq)) {
        j=j1[which(clin7[j1,varId]%in%grpUniq[k,])]
        for (a1 in 1:length(altTypeUniq2)) {
            l=which(colnames(clin7)==altTypeUniq2[a1])
            if (any(clin7[j,l]!=0,na.rm=T)) {
                res=try(wilcox_test(clin6[j,l] ~ as.factor(clin6[j,varId]), distribution="exact", conf.int=F))
                if (class(res)!="try-error") {
                    pv=try(pvalue(res))
                    if (class(res)!="try-error") {
                        out[g1,paste("pv_",altTypeUniq2[a1],"_",grpUniq[k,2],"V",grpUniq[k,1],sep="")]=pv
                    }
                }
            }
        }
    }
}
out1=out
grpUniq=grpUniq1
timeStamp=c(timeStamp,Sys.time())
print(format(timeStamp[2], "%x %X"))
print(diff(timeStamp))
out=cbind(out1,out2)
out[which(out>1)]=1
tbl=cbind(gene=rownames(out),out)
write.table(paste("Pair-wise comparison of diseases based on number of alterations for a patient"),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
write.table(paste("P-values from Wilcoxon's rank sum tests"),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
write.table(tbl,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)

#save.image(paste("tmp2",subsetFlag,".RData",sep=""))

library(qvalue)
out1=matrix(nrow=nrow(out),ncol=ncol(out),dimnames=list(rownames(out),sub("pv_","qv_",colnames(out))))
for (k in 1:ncol(out)) {
    i=which(!is.na(out[,k]))
    out1[i,k]=qvalue(out[i,k])$qvalues
    
}
fName=sub("pValue","qValue",fName)
tbl=cbind(gene=rownames(out1),out1)
write.table(paste("Pair-wise comparison of diseases based on number of alterations for a patient"),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
write.table(paste("Q-values (False discovery rate adjusted p-values from Wilcoxon's rank sum tests"),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
write.table(tbl,file=fName,append=T,col.names=T,row.names=F,sep="\t",quote=F)
write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)

#save.image(paste("tmp3",subsetFlag,".RData",sep=""))


## --------
load(paste("tmp3",subsetFlag,".RData",sep=""))


if (T) {
    cat("\n\n========================================================\n",sep="")
    cat("Comparison of number of mutations\n",sep="")
    cat("Computed number of mutations for each patient and then compared those numbers between two diseases for each gene\n",sep="")
    cat("Performed Wilcoxon's rank sum tests\n",sep="")
    cat("Q-values were computed to adjust for multiple testing\n",sep="")
    cat("Q-value < 0.05 was considered significant\n",sep="")
    cat("========================================================\n",sep="")
    cat("\n\nTotal no. of comparisons: ",ncol(out1)/length(altTypeUniq2),"\n\n",sep="")
    a1=1
    nm=sapply(colnames(out1)[grep(altTypeUniq2[a1],colnames(out1))],function(x) {
        sub("V"," vs. ",strsplit(x,"_")[[1]][3])
    },USE.NAMES=F)
    cat(nm,sep="\n")
    for (a1 in 1:length(altTypeUniq2)) {
        cat("\n\n========== ",altTypeUniq2[a1],"\n",sep="")
        res=apply(out1[,grep(altTypeUniq2[a1],colnames(out1))],2,function(x) sum(x<0.05,na.rm=T))
        nm=sapply(names(res),function(x) {
                sub("V"," vs. ",strsplit(x,"_")[[1]][3])
        },USE.NAMES=F)
        names(res)=nm
        k=which(res!=0)
        cat("No of comparisons with any significant gene: ",length(k),"\n",sep="")
        if (length(k)!=0) {
            print(res[k])
        }
    }

    res1=res

    cat("\n\n========================================================\n",sep="")
    cat("\n\nTotal no. of genes: ",nrow(out1),"\n\n",sep="")
    for (a1 in 1:length(altTypeUniq2)) {
        cat("\n\n========== ",altTypeUniq2[a1],"\n",sep="")
        res=apply(out1[,grep(altTypeUniq2[a1],colnames(out1))],1,function(x) sum(x<0.05,na.rm=T))
        k=order(res,decreasing=T); k=k[which(res[k]!=0)]
        cat("No of genes significant for any of the comparisons: ",length(k),"\n",sep="")
        if (length(k)!=0) {
            print(res[k])
        }
    }
}

## ---------------------
## ---------------------
## NOT USED
## Comparison within a disease

library(qvalue)

datadir=""
datadir="results/"
load(paste(datadir,"tmp_allAssays.RData",sep=""))

subset2Flag="_lungSmallSmallCellGI_T5aT7Assays"
subset2Flag="_tmbLevel_withinVetted_T5aT7Assays"
subset2Flag="_T5aT7Assays"
subset2Flag="_fionaVetted_T5aT7Assays"
subset2Flag="_withinSmallCell_T5aT7Assays"
subset2Flag="_lungSmallLargeCellGI_T5aT7Assays"
subset2Flag="_lungSmallFionaVettedGI_T5aT7Assays"

subset2List=c("_lungSmallSmallCellGI_T5aT7Assays","_tmbLevel_withinVetted_T5aT7Assays","_T5aT7Assays","_fionaVetted_T5aT7Assays","_withinSmallCell_T5aT7Assays","_lungSmallLargeCellGI_T5aT7Assays")
subset2List="_lungSmallFionaVettedGI_T5aT7Assays"
subset2List="_tmbLevel_T5aT7Assays"

for (subset2Flag in subset2List) {
    
    pThres=0.05
    
    outFormat="png"
    outFormat="pdf"
    
    colList=c("red","blue","purple","cyan","skyblue","orange","yellow","darkgreen")
    
    propFlag=">="
    propFlag=">"
    
    altTypeUniq1=sort(unique(clin1$altType))
    altTypeUniq2=cbind(c(altTypeUniq1,paste(altTypeUniq1[1],"+",altTypeUniq1[3]),paste(altTypeUniq1[1:3],"+",altTypeUniq1[4])),as.character(c(1:length(altTypeUniq1),10*(1)+3,10*(1:3)+4)))
    
    gene1=unique(clin1$gene[clin1$assayVersion=="T5a"])
    gene2=unique(clin1$gene[clin1$assayVersion=="T7"])
    j=which(!is.na(clin32$fionaCol1))
    j=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7") & ((clin1$id%in%clin32$id[j]) | (!clin1$disOntTerm2%in%clin3$disOntTerm2[j])))
    clinThis=cbind(clin3,clin32[,which(!names(clin32)%in%names(clin3))])[which(clin3$id%in%clin1$id[j]),]
    clinThis$fionaCol1[which(!clinThis$fionaCol1%in%c("N","Y","Y/N"))]=NA
    clinThis$fionaCol1=paste("fiona",sub("/","", clinThis$fionaCol1),sep="")
    clinThis$fionaCol1[clinThis$fionaCol1=="fionaNA"]=NA
    clinThis$cellSize[which(clinThis$cellSize%in%c("MC","SC/LC","SC/MC"))]="SMLC"
    clinThis$cellSize2=clinThis$cellSize
    clinThis$cellSize2[which(clinThis$disOntTerm2=="lungSmall")]="SC"
    clinThis$disease=clinThis$disOntTerm2
    clinThis$disease[which(clinThis$disease=="lungSmall")]="SCLC"
    clinThis$disOntTerm3[which(clinThis$disOntTerm3=="lungSmall")]="SCLC"
    clinThis$disOntTerm3[which(clinThis$disOntTerm3=="combined")]="pancreas+colon+otherGI"
    
    if (length(grep("_withinSmallCell",subset2Flag)==1)) {
        clinThis=clinThis[which(clinThis$cellSize2=="SC"),]
        heading1=paste("Small cell Fiona checked + SCLC, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
        #disUniq=cbind(c("pancreas","colon","other"),c("Pancreas","Colon","Other GI"))
        disUniq=cbind(c(""),c("SCLC, pancreas, colon, other GI"))
        varList=c("disease")
        varName=c("disease")
    } else if (length(grep("_lungSmallFionaVettedGI",subset2Flag)==1)) {
        clinThis=clinThis[which(clinThis$disOntTerm3=="SCLC" | clinThis$fionaCol1=="fionaY"),]
        heading1=paste("Fiona vetted GI + SCLC, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
        disUniq=cbind(c(""),c("SCLC, pancreas, colon, other GI"))
        varList=c("disease")
        varName=c("disease")
    } else if (length(grep("_lungSmallSmallCellGI",subset2Flag)==1)) {
        clinThis=clinThis[which(clinThis$disOntTerm3=="SCLC" | clinThis$cellSize2=="SC"),]
        heading1=paste("Small cell Fiona checked + SCLC, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
        disUniq=cbind(c(""),c("SCLC, small cell pancreas, colon, other GI pooled"))
        varList=c("disOntTerm3")
        varName=c("disease")
    } else if (length(grep("_lungSmallLargeCellGI",subset2Flag)==1)) {
        clinThis=clinThis[which(clinThis$disOntTerm3=="SCLC" | clinThis$cellSize2=="LC"),]
        heading1=paste("Large cell Fiona checked + SCLC, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
        disUniq=cbind(c(""),c("SCLC, large cell pancreas, colon, other GI pooled"))
        varList=c("disOntTerm3")
        varName=c("disease")
    } else if (length(grep("_tmbLevel",subset2Flag)==1)) {
        heading1=paste("Fiona checked, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
        disUniq=cbind(c(""),c("pancreas, colon, other GI pooled"))
        varList=c("tmbLevel")
        varName=c("tmbLevel")
    } else {
        clinThis=clinThis[which(!is.na(clinThis$fionaCol1)),]
        heading1=paste("Fiona checked, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
        disUniq=cbind(c("","pancreas","colon","other"),c("Pancreas + colon + other GI","Pancreas","Colon","Other GI"))
        varList=c("fionaCol1","cellSize")
        varName=c("vetted","cellSize")
    }
    if (length(grep("_withinVetted",subset2Flag)==1)) {
        clinThis=clinThis[which(clinThis$fionaCol1=="fionaY"),]
        heading1=sub("Fiona checked","Fiona vetted",heading1)
    }
    
    for (dId in 1:nrow(disUniq)) {
        cat("\n\n================",disUniq[dId,1],"===============\n\n")
        if (length(grep("_withinSmallCell",subset2Flag)==1)) {
            if (disUniq[dId,1]=="") {
                samIdThis=1:nrow(clinThis)
            } else {
                samIdThis=which(clinThis$disOntTerm2%in%c("lungSmall",disUniq[dId,1]))
            }
        } else {
            if (disUniq[dId,1]=="") {
                samIdThis=1:nrow(clinThis)
            } else {
                samIdThis=which(clinThis$disOntTerm2==disUniq[dId,1])
            }
        }
        for (varId in 1:length(varList)) {
            cat("\n\n================",varList[varId],"\n\n")
            if (varList[varId]=="fionaCol1") {
                if (disUniq[dId,1]=="") next
                grpUniq3=paste("fiona",sort(c("N","Y","YN")),sep="")
                grpUniq3=cbind(grpUniq3,sub("fiona","vetted ",sub("YN","Borderline",grpUniq3)),sub("fiona","",sub("YN","B",grpUniq3)))
                grpUniq3[,2]=paste("Vetted ",c("Yes","No","Borderline"),sep="")
                #grpUniq3=cbind(grpUniq3,paste("fiona",sub("/","", grpUniq3),sep=""))
            } else if (varList[varId]=="cellSize") {
                grpUniq3=c("SC","LC","SMLC")
                grpUniq3=cbind(grpUniq3,c("Small Cell","Large Cell","Small/Medium/Large Cell"),c("S","L","N"))
            } else if (length(grep("_lungSmallSmallCellGI",subset2Flag)==1) | length(grep("_lungSmallLargeCellGI",subset2Flag)==1)) {
                grpUniq3=c("SCLC","pancreas+colon+otherGI")
                grpUniq3=cbind(grpUniq3,capWords(grpUniq3),toupper(substr(grpUniq3,1,1)))
            } else if (varList[varId]=="tmbLevel") {
                grpUniq3=c("Low","Intermediate")
                grpUniq3=cbind(grpUniq3,paste("TMB level ",capWords(grpUniq3),sep=""),toupper(substr(grpUniq3,1,1)))
            } else {
                grpUniq3=c("SCLC","pancreas","colon","other")
                grpUniq3=cbind(grpUniq3,capWords(sub("other","other GI",grpUniq3)),toupper(substr(grpUniq3,1,1)))
            }
            n=nrow(grpUniq3)*(nrow(grpUniq3)-1)/2
            grpUniq32=matrix(nrow=n,ncol=2)
            k=1
            for (k1 in 1:(nrow(grpUniq3)-1)) {
                for (k2 in (k1+1):nrow(grpUniq3)) {
                    grpUniq32[k,1]=grpUniq3[k1,1]
                    grpUniq32[k,2]=grpUniq3[k2,1]
                    k=k+1
                }
            }
            outAnyAlt=datGP[,grep("anyAlt_",colnames(datGP))]
            #colnames(outAnyAlt)=sapply(colnames(outAnyAlt),function(x) strsplit(x,"_")[[1]][2],USE.NAMES=F)
            j=match(paste("anyAlt_",clinThis$id[samIdThis],sep=""),colnames(outAnyAlt))
            if (any(is.na(j))) {
                cat("Mismatch in patients !!!\n")
                break
            }
            outAnyAlt=outAnyAlt[,j]
            
            outPropAll=matrix(0,nrow=nrow(outAnyAlt),ncol=nrow(grpUniq3),dimnames=list(rownames(outAnyAlt),grpUniq3[,1]))
            outPropCombAll=rep(0,nrow(outAnyAlt))
            j=which(clinThis[samIdThis,varList[varId]]%in%grpUniq3[,1])
            for (g1 in 1:nrow(outPropAll)) {
                res=table(grp=clinThis[samIdThis[j],varList[varId]],anyAlt=outAnyAlt[g1,j])
                x1=apply(res,1,sum)
                if ("1"%in%colnames(res)) {
                    outPropAll[g1,match(rownames(res),colnames(outPropAll))]=res[,"1"]/x1
                    outPropCombAll[g1]=sum(res[,"1"])/sum(x1)
                }
            }
            
            #out4Alt=datGP[,grep("anyAlt_",colnames(datGP))]
            j=grep("anyAlt_",colnames(datGP))
            out4Alt=matrix(0,nrow=nrow(datGP),ncol=ncol(datGP[,j]),dimnames=list(rownames(datGP),colnames(datGP)[j]))
            id1=sapply(colnames(datGP),function(x) {strsplit(x,"_")[[1]][2]},USE.NAMES=F)
            id2=sapply(colnames(out4Alt),function(x) {strsplit(x,"_")[[1]][2]},USE.NAMES=F)
            jj=match(sapply(colnames(datGP),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F),altTypeUniq1)
            for (g1 in 1:nrow(out4Alt)) {
                for (a1 in 1:length(altTypeUniq1)) {
                    j=which(jj==a1 & datGP[g1,]==1)
                    j2=match(id1[j],id2)
                    j21=j2[which(out4Alt[g1,j2]==0)]
                    j22=j2[which(out4Alt[g1,j2]!=0)]
                    out4Alt[g1,j21]=a1
                    out4Alt[g1,j22]=10*out4Alt[g1,j22]+a1
                }
            }
            j=match(paste("anyAlt_",clinThis$id[samIdThis],sep=""),colnames(out4Alt))
            out4Alt=out4Alt[,j]
            
            out4AltPropCombAll=matrix(0,nrow=nrow(outAnyAlt),ncol=nrow(altTypeUniq2)+1,dimnames=list(rownames(outAnyAlt),c("noAlt",altTypeUniq2[,1])))
            colId=c("0",altTypeUniq2[,2])
            x=as.character(sort(unique(c(out4Alt))))
            k=which(!x%in%colId)
            if (length(k)!=0) cat("out4AltPropCombAll: ",x[k]," not in list !!!",sep="")
            for (g1 in 1:nrow(outPropAll)) {
                res=table(out4Alt[g1,])
                out4AltPropCombAll[g1,match(names(res),colId)]=res/sum(res)
            }
            out4AltPropCombAll=out4AltPropCombAll[,altTypeUniq2[,1]]
            
            if (F) {
                ## -------------------
                ## Long tail plot
                
                colList=c("red","blue","purple","cyan","skyblue","orange","yellow")
                fName=paste("longTailPlot_",capWords(varList[varId]),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),sep="")
                switch(outFormat,
                "png"={png(paste(fName,".png",sep=""))},
                "pdf"={pdf(paste(fName,".pdf",sep=""),width=2*7,height=7)}
                )
                x=apply(out4AltPropCombAll,1,sum)
                k=which(x!=0)
                k=which(x>=0.01)
                k=k[order(x[k],decreasing=T)]
                barplot(100*t(out4AltPropCombAll[k,]),names.arg=rownames(out4AltPropCombAll)[k],ylab="% of alteration",col=colList,las=3,cex.axis=.7)
                dev.off()
                ## -------------------
            }
            
            outPvAll=matrix(nrow=nrow(outPropAll),ncol=nrow(grpUniq32),dimnames=list(rownames(outPropAll),paste("pv_anyAlt_",rep(paste(grpUniq32[,2],grpUniq32[,1],sep="V"),each=1),sep="")))
            print(nrow(outPvAll))
            timeStamp=Sys.time()
            print(format(timeStamp, "%x %X"))
            for (k in 1:nrow(grpUniq32)) {
                j=which(clinThis[samIdThis,varList[varId]]%in%grpUniq32[k,])
                if (length(j)!=0) {
                    for (a1 in which(altTypeUniq2[,1]=="anyAlt")) {
                        for (g1 in 1:nrow(outPvAll)) {
                            x=table(clinThis[samIdThis[j],varList[varId]],outAnyAlt[g1,j])
                            if (nrow(x)==2 && ncol(x)==2) {
                                pv=fisher.test(x)$p.value
                                outPvAll[g1,paste("pv_",altTypeUniq2[a1,1],"_",grpUniq32[k,2],"V",grpUniq32[k,1],sep="")]=pv
                            }
                        }
                    }
                }
            }
            outPvAll[which(outPvAll>1)]=1
            cat("No. of comparisons with pv < 0.05",sum(c(outPvAll)<0.05,na.rm=T),"\n")
            outQvAll=matrix(nrow=nrow(outPvAll),ncol=ncol(outPvAll),dimnames=list(rownames(outPvAll),sub("pv_","qv_",colnames(outPvAll))))
            for (k in 1:ncol(outPvAll)) {
                i=which(!is.na(outPvAll[,k]))
                if (length(i)!=0) {
                    outQvAll[i,k]=qvalue(outPvAll[i,k])$qvalues
                }
            }
            cat("No. of comparisons with qv < 0.05",sum(c(outQvAll)<0.05,na.rm=T),"\n")
            
            #for (thisFlag in c("","10%patient","15%patient")) {
            for (thisFlag in c("10%patient","15%patient")) {
                heading=paste(heading1,"\n")
                if (length(grep("%patient",thisFlag))==1) {
                    propThres=0.10
                    propThres=as.integer(sub("%patient","",thisFlag))/100
                    i=c()
                    for (k in 1:ncol(outPropAll)) {
                        if (propFlag==">") {
                            i=c(i,which(round(outPropAll[,k],2)>propThres))
                        } else {
                            i=c(i,which(round(outPropAll[,k],2)>=propThres))
                        }
                    }
                    i=sort(unique(i))
                    fName=paste("proportionTable_geneBy",capWords(varList[varId]),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),"_",propThres*100,"percPatient",fName1,subset2Flag,".txt",sep="")
                    heading=paste(heading,"Table of ",length(i)," genes with alterations in ",propFlag," ",propThres*100,"% of patients in any group within ",disUniq[dId,2],sep="")
                } else {
                    i=1:nrow(outPropAll)
                    fName=paste("proportionTable_geneBy",capWords(varList[varId]),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),".txt",sep="")
                    heading=paste(heading,"Table of ",length(i)," genes with any alteration in a group within ",disUniq[dId,2],sep="")
                }
                if (length(i)!=0) {
                    nmR=rownames(outPropAll)[i]
                    outPropThis=round(outPropAll[i,colnames(outPropAll)],2)
                    if (length(i)>1) {
                        i=match(rownames(outQvAll),rownames(outPropThis)); i1=which(!is.na(i)); i2=i[i1]
                        outPercThis=apply(outPropThis,c(1,2),function(x) {
                            y=paste(" ",x*100,"% ",sep="")
                            y
                        })
                    } else if (length(i)==1) {
                        outPercThis=matrix(paste(" ",outPropThis*100,"% ",sep=""),nrow=1,dimnames=list(nmR,colnames(outPropAll)))
                    }
                    colnames(outPercThis)=paste("anyAlt_",colnames(outPercThis),sep="")
                    for (grpId in 1:nrow(grpUniq3)) {
                        k1=grep("qv_anyAlt_",colnames(outQvAll))
                        k1=k1[grep(grpUniq3[grpId,1],colnames(outQvAll)[k1],fixed=T)]
                        nm1=colnames(outQvAll)[k1]
                        nm11=nm1
                        nm1=sub("qv_anyAlt_","",sub(paste("V",grpUniq3[grpId,1],sep=""),"",sub(paste(grpUniq3[grpId,1],"V",sep=""),"",nm1)))
                        nm21=colnames(outPercThis)[grep("anyAlt_",colnames(outPercThis))]
                        nm2=sub("anyAlt_","",colnames(outPercThis)[grep("anyAlt_",colnames(outPercThis))])
                        x=" "
                        for (k in 1:length(k1)) {
                            k21=which(nm2==nm1[k])
                            i=which(outQvAll[i1,nm11[k]]<pThres)
                            if (length(i)!=0) {
                                outPercThis[i2[i],k21]=paste(outPercThis[i2[i],k21],grpUniq3[grpId,3],sep="")
                            }
                        }
                    }
                    x=table(clinThis[samIdThis,varList[varId]][which(clinThis[samIdThis,varList[varId]]%in%grpUniq3[,1])])
                    x=x[order(x,decreasing=T)]
                    nm=names(x)
                    colnames(outPercThis)=sub("anyAlt_","",colnames(outPercThis))
                    outPercThis=outPercThis[,nm]
                    if (!"matrix"%in%class(outPercThis)) {
                        outPercThis=matrix(outPercThis,nrow=1,dimnames=list(nmR,colnames(outPropAll)))
                    }
                    nm=grpUniq3[match(colnames(outPercThis),grpUniq3[,1]),2]
                    nm=paste(c("gene",paste(nm," (N=",x,")",sep="")),collapse="\t")
                    tbl=cbind(gene=rownames(outPercThis),data.frame(outPercThis))
                    write.table(heading,file=fName,append=F,col.names=F,row.names=F,sep="\t",quote=F)
                    write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                    write.table(nm,file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                    write.table(tbl,file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                    
                    ## -------------------
                    ## Long tail plot
                    
                    fName=sub("proportionTable_","longTailPlot_",sub(".txt","",fName,fixed=T))
                    heading=sub("Table of ","Long tail plot of ",heading)
                    switch(outFormat,
                    "png"={png(paste(fName,".png",sep=""))},
                    "pdf"={pdf(paste(fName,".pdf",sep=""),width=2*7,height=7)}
                    )
                    x=apply(out4AltPropCombAll,1,sum)
                    k=match(rownames(outPropThis),rownames(out4AltPropCombAll))
                    k=k[order(x[k],decreasing=T)]
                    barplot(100*t(out4AltPropCombAll[k,]),names.arg=rownames(out4AltPropCombAll)[k],main=heading,ylab="% of alteration",col=colList,las=3,cex.axis=.7)
                    dev.off()
                    ## -------------------
                    
                }
            }
        }
    }
    
}

fName="_forLongTailPlot"
width = 480; height = 480
if (outFormat=="png") {
    png(paste("altTypeColorLegend",fName,".png",sep=""),width=width,height=height)
} else {
    pdf(paste("altTypeColorLegend",fName,".pdf",sep=""))
}
grpUniq=altTypeUniq2[,1]
cexThis=NULL
if (outFormat=="pdf") cexThis=1
sampleColorLegend(tls=grpUniq,col=colList,legendTitle="alteration type",cex=cexThis)
dev.off()


if (F) {
    ## Testing
    outPropAll[(outPropAll[,"fionaY"]>.1),]
    
    geneThis="MEN1"
    geneThis="DAXX"
    x2=clinThis[samIdThis,]
    x1=outAnyAlt
    colnames(x1)=sub("anyAlt_","",colnames(x1))
    i=which(rownames(x1)==geneThis)
    j=which(x2$fionaCol1%in%c("fionaN","fionaY","fionaYN"))
    j=which(x2$fionaCol1%in%c("fionaY","fionaYN"))
    j=which(x2$fionaCol1%in%c("fionaN","fionaY"))
    table(x1[i,j],x2$fionaCol1[j])
    fisher.test(x1[i,j],x2$fionaCol1[j])
}

rm(outAnyAlt,outPropAll,outPropThis,outPercThis,outPvAll,outQvAll,out4Alt,out4AltPropCombAll)


if (F) {
    ## Check if any patient has > 1 type of alteration for a gene
    ## There are some
    k1=grep("anyAlt_",colnames(datGP))
    print(max(k1))
    y=c()
    for (k in k1) {
        if (k%%1000==0) print(k)
        x=apply(datGP[,k+(1:4)],1,sum)
        if (any(x>1)) y=c(y,k)
    }
}

## ---------------------
## ---------------------
## Distance of each sample to disease centroid

getDistance=function(x,y,method="euclidean") {
    switch(method,
    "euclidean"={
        z=sqrt(sum((x-y)^2,na.rm=T))
    }
    )
    z
}

arrayData=datGP[,which(substr(colnames(datGP),1,nchar("anyAlt_"))=="anyAlt_")]
arrayData2=arrayData

grpUniq=c("colon","lungSmall","pancreas","skin")
tmpC=rep("",ncol(arrayData))
dat=data.frame(id=clin3$id[samId03],disOntTerm2=tmpC,stringsAsFactors=F)
j1=samId03[which(clin3$disOntTerm2[samId03]%in%grpUniq)]
j=match(clin3$id[j1],dat$id); j1=j1[which(!is.na(j))]; j2=j[which(!is.na(j))]
dat$disOntTerm2[j2]=clin3$disOntTerm2[j1]
disDistMat=matrix(nrow=nrow(dat),ncol=length(grpUniq),dimnames=list(dat$id,grpUniq))
jj=which(dat$disOntTerm2=="")
jj=1:nrow(dat)
for (gId in 1:length(grpUniq)) {
    j=which(dat$disOntTerm2==grpUniq[gId])
    centr=apply(arrayData2[,j],1,mean,na.rm=T)
    x=apply(arrayData2[,jj],2,getDistance,y=centr)
    disDistMat[jj,gId]=x
}

j=which(clin3$disOntTerm3[samId03]=="unknown")
j=1:length(samId03)
tmp=disDistMat[j,]
colnames(tmp)=paste("dist2",colnames(disDistMat),"_euclidean",sep="")
tbl=cbind(clin3[samId03[j],],tmp)
tmp=t(apply(disDistMat[j,],1,function(x) {
    y=order(order(x))
    y[is.na(x)]=NA
    y
}))
colnames(tmp)=paste("dist2",colnames(disDistMat),"_rank",sep="")
tbl=cbind(tbl,tmp)
write.table(tbl,file=paste("dist2disease",subsetFlag,".txt",sep=""),append=F,col.names=T,row.names=F,sep="\t",quote=F)

## ---------------------
## ---------------------
* Co-mutations
** Consider genes with high mutation
** Heatmap?

if (F) {
load(paste("tmp3",subsetFlag,".RData",sep=""))

x=datGD[,grep("anyAltCount_",colnames(datGD))]
out=matrix()
cutoff=10
y=x
y[y<cutoff]=NA
y=y[apply(y,1,function(x) {sum(!is.na(x))>1}),]
dim(y)
}

## ---
load(paste("tmp",subsetFlag,".RData",sep=""))

x=table(clin3$disOntTerm2[which(clin3$id%in%clin1$id[samId])])
x=x[order(x,decreasing=T)]
disUniq=names(x)
disUniq=cbind(disUniq,toupper(substr(disUniq,1,1)))
disUniq[,2]=sub("L","S",disUniq[,2])
disUniq[,2][match(c("lungSmall","unknown","lungLarge","pancreas","skin","colon","other","prostate","head","softTissuePrimitive","breast","cervix","brain","ovary","bladder","uterus","softTissueSmall","vagina"),disUniq[,1])]=c("Ls","Un","Ll","Pa","S","Co","Ot","Pr","H","Sp","Brt","Ce","Brn","Ov","Bl","Ut","Sm","V")


## --------------
heading1="All assays, any alteration"
pThres=0.05
genesetFlag="_mmr"
genesetFlag="_notch"
genesetFlag="_setd2"
genesetFlag="_brca"

## --------------
switch(genesetFlag,
    "_signifPairwise_12genes"={
        fName1=genesetFlag
        fName=paste("proportionTable_geneByDisease_10percPatient",fName1,subsetFlag,".txt",sep="")
        heading="Table of 15 genes with alterations in at least 10% of patients in any group"
        i=c()
        propThres=0.10
        for (k in grep("anyAlt_",colnames(datGD))) {
            if (grep("_12genes",fName1)) {
                i=c(i,which(round(datGD[,k],2)>propThres))
            } else {
                i=c(i,which(round(datGD[,k],2)>=propThres))
            }
        }
        i=sort(unique(i))

        n=length(i)*(length(i)-1)/2
        tmpC=rep("",n)
        annThis=data.frame(gene1=tmpC,gene2=tmpC,stringsAsFactors=F)
        k=1
        for (i1 in 1:(length(i)-1)) {
            for (i2 in (i1+1):length(i)) {
                annThis$gene1[k]=rownames(datGD)[i[i1]]
                annThis$gene2[k]=rownames(datGD)[i[i2]]
                k=k+1
            }
        }
    },
    "_mmr"={
        fName1=genesetFlag
        heading=paste(heading1,", MMR genes",sep="")
        tbl=read.table(paste("docs/candGene.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        tbl=tbl[which(tbl$family=="Mismatch excision repair (MMR)"),]
        x=unique(c(tbl$gene,unlist(sapply(tbl$alias,function(x) {strsplit(x,", ")[[1]]},USE.NAMES=F))))
        i=which(toupper(rownames(datGD))%in%x)

        n=length(i)*(length(i)-1)/2
        tmpC=rep("",n)
        annThis=data.frame(gene1=tmpC,gene2=tmpC,stringsAsFactors=F)
        k=1
        for (i1 in 1:(length(i)-1)) {
            for (i2 in (i1+1):length(i)) {
                annThis$gene1[k]=rownames(datGD)[i[i1]]
                annThis$gene2[k]=rownames(datGD)[i[i2]]
                k=k+1
            }
        }
    },
    "_notch"={
        fName1=genesetFlag
        heading=paste(heading1,", NOTCH genes",sep="")
        i=grep("NOTCH",toupper(rownames(datGP)))
        n=length(i)*(length(i)-1)/2
        tmpC=rep("",n)
        annThis=data.frame(gene1=tmpC,gene2=tmpC,stringsAsFactors=F)
        k=1
        for (i1 in 1:(length(i)-1)) {
            for (i2 in (i1+1):length(i)) {
                annThis$gene1[k]=rownames(datGD)[i[i1]]
                annThis$gene2[k]=rownames(datGD)[i[i2]]
                k=k+1
            }
        }
    },
    "_brca"={
        fName1=genesetFlag
        heading=paste(heading1,", BRCA genes",sep="")
        x=c("BRCA1","BRCA2","ATM","ATR","PTEN","PARP1","MRE11A","PALB2","CDK12","BAP1","PRKDC","BARD1","BACH1","FANCD2","FANCE","CHEK1","CHEK2","PMS2","FANCA","FANCC","FANCF","XPA","NBN","BLM","BRIP1","RAD51B","RAD51C","RAD51D","RAD50")
        i=which(toupper(rownames(datGD))%in%x)
        n=length(i)*(length(i)-1)/2
        tmpC=rep("",n)
        annThis=data.frame(gene1=tmpC,gene2=tmpC,stringsAsFactors=F)
        k=1
        for (i1 in 1:(length(i)-1)) {
            for (i2 in (i1+1):length(i)) {
                annThis$gene1[k]=rownames(datGD)[i[i1]]
                annThis$gene2[k]=rownames(datGD)[i[i2]]
                k=k+1
            }
        }
    },
    "_setd2"={
        fName1=genesetFlag
        heading=paste(heading1,", SETD2 vs. other genes (pair-wise comparison)",sep="")
        i=which(toupper(rownames(datGD))=="SETD2")
        annThis=data.frame(gene1=rep("SETD2",nrow(datGD)-1),gene2=rownames(datGD)[-i],stringsAsFactors=F)
    }
)

## --------------
id=sapply(colnames(datGP),function(x){strsplit(x,"_")[[1]][2]},USE.NAMES=F)
disOntTerm2=clin3$disOntTerm2[match(id,clin3$id)]
n=length(i)*(length(i)-1)/2
for (a1 in which(altTypeUniq2=="anyAlt")) {
    cat("\n\n============ ",altTypeUniq2[a1]," =============\n",sep="")
    varId=grep(altTypeUniq2[a1],colnames(datGP))
    out1=matrix(nrow=nrow(annThis),ncol=nrow(disUniq))
    colnames(out1)=disUniq[,1]
    #out2=matrix(nrow=n,ncol=nrow(disUniq))
    #colnames(out2)=disUniq[,1]
    for (k in 1:nrow(annThis)) {
        for (d1 in 1:nrow(disUniq)) {
            j=varId[which(disOntTerm2[varId]==disUniq[d1,1])]
            if (F) {
                out=nm=c()
                out22=c()
                for (i1 in 1:(length(i)-1)) {
                    for (i2 in (i1+1):length(i)) {
                        out=c(out,mean(c(datGP[i1,j]==1 & datGP[i2,j]==1),na.rm=T))
                        nm=c(nm,paste(rownames(datGP)[i1],"_",rownames(datGP)[i2],sep=""))
                        out22=c(out22,sum(c(datGP[i1,j]==datGP[i2,j]),na.rm=T))
                    }
                }
                out1[,d1]=out
                out2[,d1]=out22
            }
            i1=which(rownames(datGP)==annThis$gene1[k])
            i2=which(rownames(datGP)==annThis$gene2[k])
            out1[k,d1]=mean(c(datGP[i1,j]==1 & datGP[i2,j]==1),na.rm=T)
            
            #out2[,d1]=sum(c(datGP[i1,j]==datGP[i2,j]),na.rm=T)
        }
    }
    nm=paste(annThis$gene1,annThis$gene2,sep="_")
    rownames(out1)=nm
    #out1=round(out1,2)
    out1=signif(out1,2)
    cutoff=.25
    cutoff=.10
    k=apply(out1,1,function(x) any(x>=cutoff))
    
    cat("No. of gene-pairs with >= ",cutoff," proportion of patients with co-mutation in any disease group: ",sum(k),sep="")
    out1[k,]
}

id0=unique(clin1$id[samId][which(clin1$disOntTerm2[samId]=="pancreas")])
id1=unique(clin1$id[samId][which(!is.na(clin1$altType[samId]) & clin1$disOntTerm2[samId]=="pancreas" & clin1$gene[samId]=="TP53")])
id2=unique(clin1$id[samId][which(!is.na(clin1$altType[samId]) & clin1$disOntTerm2[samId]=="pancreas" & clin1$gene[samId]=="APC")])
mean(id0%in%id1 & id0%in%id2)
out1[2,]

library(qvalue)
p=nrow(disUniq)*(nrow(disUniq)-1)/2
out1=apply(out1,c(1,2),function(x) paste(x," ",sep=""))
if (nrow(out1)<10) multTestFlag="bonferroni" else multTestFlag="qvalue"
#multTestFlag="bonferroni"
switch(multTestFlag,
"bonferroni"={
    heading=paste(heading,", significance determined as bonferroni corrected p-value < ",pThres,sep="")
},
"qvalue"={
    heading=paste(heading,", significance determined as q-value < ",pThres,sep="")
}
)
for (a1 in which(altTypeUniq2=="anyAlt")) {
    cat("\n\n============ ",altTypeUniq2[a1]," =============\n",sep="")
    varId=grep(altTypeUniq2[a1],colnames(datGP))
    out3=matrix(nrow=nrow(annThis),ncol=p)
    p=1
    nm2=c()
    for (d1 in 1:(nrow(disUniq)-1)) {
        j1=varId[which(disOntTerm2[varId]==disUniq[d1,1])]
        for (d2 in (d1+1):nrow(disUniq)) {
            j2=varId[which(disOntTerm2[varId]==disUniq[d2,1])]
            if (F) {
                out=nm=c()
                for (i1 in 1:(length(i)-1)) {
                    for (i2 in (i1+1):length(i)) {
                        j=c(j1,j2)
                        x=table(rep(disOntTerm2[j],each=2),c(datGP[i1,j],datGP[i2,j]))
                        if (nrow(x)==2 && ncol(x)==2) {
                            pv=fisher.test(x)$p.value
                        } else {
                            pv=NA
                        }
                        out=c(out,pv)
                        nm=c(nm,paste(rownames(datGP)[i1],"_",rownames(datGP)[i2],sep=""))
                    }
                }
            }
            
            out=rep(NA,nrow(annThis))
            nm=rep("",nrow(annThis))
            j=c(j1,j2)
            for (k in 1:nrow(annThis)) {
                i1=which(rownames(datGP)==annThis$gene1[k])
                i2=which(rownames(datGP)==annThis$gene2[k])
                #x=table(rep(disOntTerm2[j],2),c(datGP[i1,j],datGP[i2,j]))
                x=table(disOntTerm2[j],as.integer(datGP[i1,j]==1 & datGP[i2,j]==1))
                if (nrow(x)==2 && ncol(x)==2) {
                    pv=fisher.test(x)$p.value
                } else {
                    pv=NA
                }
                out[k]=pv
                nm[k]=paste(rownames(datGP)[i1],"_",rownames(datGP)[i2],sep="")
            }

            nm2=c(nm2,paste(disUniq[d1,1],"_",disUniq[d2,1],sep=""))
            out[which(out>1)]=1
            out3[,p]=out
            ii=which(!is.na(out))
            pv=rep(NA,length(ii))
            if (multTestFlag=="bonferroni") {
                pv=p.adjust(out[ii],method=multTestFlag)
            } else {
                res=try(qvalue(out[ii]))
                if (class(res)=="qvalue") {
                    pv=res$qvalues
                }
            }
            out3[ii,p]=pv
            ii=which(out3[,p]<0.05)
            if (length(ii)!=0) {
                k1=which(colnames(out1)%in%disUniq[d1,1])
                out1[ii,k1]=paste(out1[ii,k1],disUniq[d2,2],sep="")
                k2=which(colnames(out1)%in%disUniq[d2,1])
                out1[ii,k2]=paste(out1[ii,k2],disUniq[d1,2],sep="")
            }
            p=p+1
        }
    }
    rownames(out3)=nm
    colnames(out3)=nm2
    out
    cutoff=.25
    k=apply(out1,1,function(x) any(x>=cutoff))
    cutoff=0.05
    k=apply(out3,1,function(x) any(x<pThres,na.rm=T))
    cat("No. of gene-pairs with p-value/q-value < ",pThres," for at least 1 pair of diseases: ",sum(k),sep="")
    if (any(k)) {
        out1[k,]
        k3=apply(out3,2,function(x) {any(!is.na(x))})
        k3=apply(out3,2,function(x) {any(x<cutoff,na.rm=T)})
        out3[k,k3]
        k1=unique(c(sapply(colnames(out3)[k3],function(x) {strsplit(x,"_")[[1]]},USE.NAMES=F)))
        out1[k,k1]
        nm=t(sapply(rownames(out3)[k],function(x) {strsplit(x,"_")[[1]]},USE.NAMES=F))
        colnames(nm)=c("gene1","gene2")
        tbl=cbind(nm,out1[k,k1])
        fName=paste("frequencyTable_anyAlt_disease",subsetFlag,fName1,sep="")
        x=table(disOntTerm2[varId])
        x=x[match(k1,names(x))]
        x=paste(paste(rep("\t",1),collapse=""),"No. of patients\t",paste(x,collapse="\t"),sep="")
        "Proportion of patients with mutation in both genes"
        write.table(heading,file=paste(fName,".txt",sep=""),append=F,col.names=F,row.names=F,sep="\t",quote=F)
        write.table(x,file=paste(fName,".txt",sep=""),append=T,col.names=F,row.names=F,sep="\t",quote=F)
        write.table(tbl,file=paste(fName,".txt",sep=""),append=T,col.names=T,row.names=F,sep="\t",quote=F)
        for (kk in 1:length(nm)) {
            i1=which(rownames(datGP)==nm[kk,1])
            i2=which(rownames(datGP)==nm[kk,2])
            #j1=varId[which(disOntTerm2[verId]==colnames())
        }
    }
}
"
Co-mutation of genes for the 12 genes where > 10% of patients have a mutation for any of the 4 disease groups:

============ anyAlt =============
No. of gene-pairs with >= 0.25 proportion of patients with a mutation in any disease group: 30

============ amplification =============
No. of gene-pairs with >= 0.25 proportion of patients with a mutation in any disease group: 0

============ loss =============
No. of gene-pairs with >= 0.25 proportion of patients with a mutation in any disease group: 0

============ rearrangement =============
No. of gene-pairs with >= 0.25 proportion of patients with a mutation in any disease group: 0

============ short variant =============
No. of gene-pairs with >= 0.25 proportion of patients with a mutation in any disease group: 27


Proportion of patients with any mutation for a pair of genes
Tabulated only those with >= 0.25 proportion of patients with a gene-pair mutation in any disease group
--------------------------------------------
            colon lungSmall other pancreas
--------------------------------------------
lungSmall pancreas colon other
TP53_RB1           0.65     0.04  0.27  0.25
TP53_APC           0.02     0.01  0.33  0.05
TP53_MLL2          0.12     0.00  0.04  0.02
TP53_LRP1B         0.11     0.02  0.03  0.03
TP53_KRAS          0.04     0.04  0.25  0.02
TP53_PTEN          0.08     0.02  0.03  0.07
TP53_CDKN2A        0.03     0.04  0.02  0.10
TP53_MYCL1         0.06     0.00  0.00  0.00
TP53_MYC           0.06     0.02  0.08  0.08
TP53_MEN1          0.00     0.06  0.01  0.00
TP53_PIK3CA        0.05     0.00  0.03  0.00
RB1_APC            0.01     0.00  0.21  0.03
RB1_MLL2           0.08     0.00  0.04  0.00
RB1_LRP1B          0.08     0.00  0.01  0.02
RB1_KRAS           0.02     0.02  0.17  0.00
RB1_PTEN           0.06     0.02  0.03  0.03
RB1_CDKN2A         0.01     0.00  0.02  0.03
RB1_MYCL1          0.05     0.00  0.00  0.00
RB1_MYC            0.04     0.00  0.01  0.05
RB1_MEN1           0.00     0.00  0.00  0.00
RB1_PIK3CA         0.04     0.00  0.01  0.00
APC_MLL2           0.01     0.00  0.03  0.00
APC_LRP1B          0.00     0.01  0.03  0.00
APC_KRAS           0.00     0.01  0.25  0.00
APC_PTEN           0.00     0.00  0.02  0.02
APC_CDKN2A         0.01     0.01  0.04  0.03
APC_MYCL1          0.00     0.00  0.00  0.00
APC_MYC            0.00     0.00  0.04  0.02
APC_MEN1           0.00     0.00  0.00  0.00
APC_PIK3CA         0.00     0.00  0.04  0.00
MLL2_LRP1B         0.01     0.00  0.00  0.00
MLL2_KRAS          0.01     0.00  0.04  0.00
MLL2_PTEN          0.01     0.00  0.00  0.00
MLL2_CDKN2A        0.01     0.00  0.00  0.00
MLL2_MYCL1         0.01     0.00  0.00  0.00
MLL2_MYC           0.01     0.00  0.00  0.00
MLL2_MEN1          0.00     0.01  0.00  0.00
MLL2_PIK3CA        0.00     0.00  0.00  0.00
LRP1B_KRAS         0.01     0.01  0.02  0.00
LRP1B_PTEN         0.01     0.00  0.00  0.02
LRP1B_CDKN2A       0.01     0.01  0.00  0.00
LRP1B_MYCL1        0.01     0.00  0.00  0.00
LRP1B_MYC          0.01     0.01  0.01  0.02
LRP1B_MEN1         0.00     0.02  0.00  0.00
LRP1B_PIK3CA       0.01     0.00  0.00  0.00
KRAS_PTEN          0.01     0.00  0.02  0.00
KRAS_CDKN2A        0.01     0.03  0.03  0.00
KRAS_MYCL1         0.01     0.00  0.00  0.00
KRAS_MYC           0.00     0.02  0.05  0.03
KRAS_MEN1          0.00     0.01  0.00  0.00
KRAS_PIK3CA        0.00     0.00  0.03  0.00
PTEN_CDKN2A        0.00     0.02  0.01  0.07
PTEN_MYCL1         0.01     0.00  0.00  0.00
PTEN_MYC           0.01     0.00  0.00  0.00
PTEN_MEN1          0.00     0.02  0.00  0.00
PTEN_PIK3CA        0.01     0.00  0.01  0.00
CDKN2A_MYCL1       0.00     0.00  0.00  0.00
CDKN2A_MYC         0.01     0.02  0.00  0.00
CDKN2A_MEN1        0.00     0.07  0.00  0.00
CDKN2A_PIK3CA      0.01     0.00  0.00  0.02
MYCL1_MYC          0.00     0.00  0.00  0.00
MYCL1_MEN1         0.00     0.01  0.00  0.00
MYCL1_PIK3CA       0.00     0.00  0.00  0.00
MYC_MEN1           0.00     0.01  0.00  0.00
MYC_PIK3CA         0.01     0.00  0.00  0.00
MEN1_PIK3CA        0.00     0.00  0.00  0.00
--------------------------------------------
"

## ---


## ---------------------
"
* PD-L protein related genes
PD-L1 protein: gene CD274
PD-L2 protein: gene PDCD1LG2, CD273
There are 19 entries for PD-L proteins
CD274 PDCD1LG2
9       10

"

g1=which(rownames(ann4)%in%c("CD274","PDCD1LG2","CD273"))
summary(c(ann4[g1,]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#0.0     0.0     0.0   106.5     0.0  2625.0


## ---------------------
# SETD2 - gene
for (a1 in which(altTypeUniq2=="anyAlt")) {
}

## ---------------------
## ---------------------
* Other projects: Would be great to assess frequency of SETD2 mutations by disease group (focused on samples tested with platform that included SETD2) and characterize type of mutation— this relates to another project with ashworth lab
* Other projects: Would be great to look at DNA repair alterations by disease group: BRCA1/2, ATM, PALB2, BRIP1, BAP1, CDK12, FANCA, CHEK2, RAD51B, and RAD51C

i=which(toupper(rownames(datGP))%in%c("SETD2"))
j1=sub("anyAlt_","",colnames(datGP))
j=match(,clin1$id)
table(datGP[i])

## ----------------------------------------------
## Tabulate mutations
## NOT USED. Code should be fine. Results not used

fName=""
fName=subsetFlag

grpUniq=sort(unique(clin1$disOntTerm))

gene=table(clin1$gene)
gene[order(gene,decreasing=T)][1:10]
varList=c("altType")

k=1
geneList=order(gene,decreasing=T); geneList=geneList[which(gene[geneList]>9)]
geneList=order(gene,decreasing=T)

nm=c("numTot","numNoMut","numMut","propNoMut","propMut")
out=matrix(0,nrow=length(geneList),ncol=length(nm)*length(grpUniq),dimnames=list(names(gene)[geneList],paste(rep(nm,length(grpUniq)),"_",rep(grpUniq,each=length(nm)),sep="")))
j0=1:nrow(clin1)
j1=j0[which(is.na(clin1$gene[j0]))]
x2=table(disOntTerm=as.factor(clin1$disOntTerm)[j1])
for (g1 in geneList) {
    j1=j0[which(clin1$gene[j0]==names(gene)[g1] & !is.na(clin1$altType[j0]))]
    j1=j1[!duplicated(clin1$id[j1])]
    k1=which(rownames(out)==names(gene)[g1])
    x1=table(as.factor(clin1$disOntTerm)[j1])
    out[k1,match(paste("numNoMut_",names(x2),sep=""),colnames(out))]=x2
    out[k1,match(paste("numMut_",names(x1),sep=""),colnames(out))]=x1
    out[k1,match(paste("numTot_",names(x1),sep=""),colnames(out))]=out[k1,match(paste("numNoMut_",names(x1),sep=""),colnames(out))]+out[k1,match(paste("numMut_",names(x1),sep=""),colnames(out))]
    out[k1,match(paste("propNoMut_",rownames(x1),sep=""),colnames(out))]=out[k1,match(paste("numNoMut_",rownames(x1),sep=""),colnames(out))]/out[k1,match(paste("numTot_",rownames(x1),sep=""),colnames(out))]
    out[k1,match(paste("propMut_",rownames(x1),sep=""),colnames(out))]=out[k1,match(paste("numMut_",rownames(x1),sep=""),colnames(out))]/out[k1,match(paste("numTot_",rownames(x1),sep=""),colnames(out))]
}
out[is.na(out)]=NA
write.table(cbind(gene=rownames(out),out),file=paste("frequencyTable",fName,".txt",sep=""),append=F,col.names=T,row.names=F,sep="\t",quote=F)
propMat=out

####################################################################
####################################################################
# NOT USED
library(marray)
#source(paste(dirSrc,"functions/heatmap.5.R",sep=""))
#source(paste(dirSrc,"functions/heatmapAcgh.7.R",sep=""))
source(paste(dirSrc,"functions/heatmap.5.3.R",sep=""))
source(paste(dirSrc,"functions/heatmapAcgh.7.1.R",sep=""))

asc <- function(x) { strtoi(charToRaw(x),16L) }
chr <- function(n) { rawToChar(as.raw(n)) }

colId2="etoeprt"
colId2=c("gene","partnerGene","altType","cdsEffect","proteinEffect","cn","fracReads")

datadir="results/"
load(paste(datadir,"tmp_allAssays.RData",sep=""))
annColThis=clin3[,which(!names(clin3)%in%colId2)]

datadir="results/"
#load(paste(datadir,"tmp_fiona102_T5aT7Assays.RData",sep=""))
load(paste(datadir,"tmp_fiona166_T5aT7Assays.RData",sep=""))
#load(paste(datadir,"tmp_allAssays.RData",sep=""))
#load(paste(datadir,"tmp_T5aT7Assays.RData",sep="")); subsetFlag="_T5aT7Assays"

annCol=annColThis
rm(annColThis)
if (exists("clin32")) {
    annCol=cbind(annCol,clin32[,which(!names(clin32)%in%names(annCol))])
}
if (exists("disDistMat")) {
    j=match(annCol$id,rownames(disDistMat)); j1=which(!is.na(j)); j2=j[j1]
    tmp=matrix(nrow=nrow(annCol),ncol=ncol(disDistMat),dimnames=list(annCol$id,paste("dist2",colnames(disDistMat),"_euclidean",sep="")))
    tmp[j1,]=disDistMat[j2,]
    annCol=cbind(annCol,tmp)
    tmp=matrix(nrow=nrow(annCol),ncol=ncol(disDistMat),dimnames=list(annCol$id,paste("dist2",colnames(disDistMat),"_rank",sep="")))
    tmp[j1,]=disDistMat[j2,]
    tmp[j1,]=t(apply(disDistMat[j2,],1,function(x) {
        y=order(order(x))
        y[is.na(x)]=NA
        y
    }))
    annCol=cbind(annCol,tmp)
    
}
annColAll=annCol

altTypeFlag="_anyAlt"
altTypeFlag="_4Alt"
altTypeList=c("_anyAlt","_anyAlt1minus1","_4Alt")
altTypeList=c("_4Alt")

#subsetFlag=""
#subsetFlag="_allAssays"

subset2Flag=""
subset2Flag="_giUnknown"
subset2Flag="_tmbHigh"

datTypeFlag="_prop"
datTypeFlag="_binary"
datTypeFlag="_patient"

subset2Flag=""; genesetFlag="_10percPatient"
subset2Flag="_giUnknown"; genesetFlag="_10percPatient"
subset2Flag=""; genesetFlag="_setd2"
subset2Flag=""; genesetFlag="_notch"
subset2Flag=""; genesetFlag="_mmr"
subset2Flag=""; genesetFlag="_brca"
subset2Flag=""; genesetFlag="_allGenes"

subsetName=subsetFlag
if (subsetFlag=="") {
    heading1=paste("T5a/T7 assays",sep="")
    subsetName="_T5aT7Assays"
    subset2Flag="_giUnknown"
    subset2Flag="_noLung"
    genesetList=c("_10percPatient","_setd2","_notch","_mmr","_brca")
    genesetList=c("_setd2","_notch","_mmr","_brca")
} else {
    genesetList=c("_setd2","_notch","_mmr","_brca")
    genesetList=c("_notch")
    genesetList=c("_allGenes")
    genesetList=c("_setd2","_notch","_mmr","_brca","_allGenes")
    switch(subsetFlag,
    "_allAssays"={
        heading1=paste("All assays",sep="")
        subset2Flag=""
        subset2Flag="_noLung"
    },
    "_T5aT7Assays"={
        heading1=paste("T5a & T7 assays",sep="")
        subset2Flag="_fiona"
        #genesetList=c("_allGenes")
    },
    "_fiona102"={
        heading1=paste("Fiona certified colon & pancreas + SCLC & otherGI",sep="")
        subset2Flag=""
    },
    "_fiona102_T5aT7Assays"={
        heading1=paste("Fiona certified colon & pancreas + rest",sep="")
        subset2Flag=""
        subset2Flag="_unknown"
        #genesetList=c("_allGenes")
    },
    "_fiona166_T5aT7Assays"={
        heading1=paste("Fiona certified colon, pancreas & otherGI + rest",sep="")
        subset2Flag="_unknown"
        subset2Flag=""
        #genesetList=c("_allGenes")
        
        heading1=paste("Fiona certified, TMB high",sep="")
        subset2Flag="_tmbHigh"
        
        heading1=paste("Fiona certified colon, pancreas & otherGI",sep="")
        subset2Flag="_gi"
        
        heading1=paste("Fiona certified colon, pancreas & otherGI + SCLC",sep="")
        subset2Flag="_giLungSmall"
    }
    )
}

subset3Flag=rep("",2)
subset3Flag[1]="_no0AltGene"
subset3Flag[1]=""

getCorFlag=T
getCorFlag=F

outFormat="png"
outFormat="pdf"

sampleBar=""
sampleBar="cluster"
geneBar=""
geneBar="clusterPr"

distMethod="kendall"
distMethod="pearson"
distMethod="spearman"
distMethod="euclidean"
distMethod="cosine"
distMethodList=c("euclidean","cosine")
distMethodList=c("euclidean")

linkMethod="ward.D2"

colList=c("skyblue","blue","yellow","purple","red")
colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","magenta","purple","darkgreen")
colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","magenta","purple","darkgreen","limegreen","salmon","gray","gold","antiquewhite","steelblue","aquamarine","lightcyan","turquoise","hotpink","black")
png("tmp.png")
barplot(1:length(colList),col=colList)
dev.off()
colList2=c("skyblue","blue")
colList2=c("skyblue","blue","yellow","purple","red","green","cyan")
colHMCon=list("red","blue","grey")
colHMCat=list(c("red","blue","purple","cyan"),NULL,NULL)

for (altTypeFlag in altTypeList) {
    for (distMethod in distMethodList) {
        
        for (genesetFlag in genesetList) {
            
            nClust=c(NA,NA)
            nClust=c(3,3)
            nClust=c(1,5)
            nClust=c(5,5)
            nClust=c(NA,5)
            
            varList=varName=NULL
            varFList=varFName=NULL
            
            fName1=""
            fNameOut=paste(datTypeFlag,subsetName,subset2Flag,subset3Flag[1],subset3Flag[2],"_",distMethod,altTypeFlag,sep="")
            fNameOut2=paste(datTypeFlag,sep="")
            limit1=c(0,1)
            
            cat("\n\n======================",fNameOut,genesetFlag,"==========================\n\n",sep="")
            
            samId3=which(clin3$id%in%clin1$id[samId])
            heading=heading1
            switch(subset2Flag,
            "_giUnknown"={
                samId3=samId3[which(clin3$disOntTerm3[samId3]%in%c("combined","unknown"))]
                heading=paste(heading,", gi+unkown",sep="")
            },
            "_gi"={
                samId3=samId3[which(clin3$disOntTerm3[samId3]%in%c("combined"))]
                heading=paste(heading,sep="")
                #heading=paste(heading,", gi",sep="")
            },
            "_giLungSmall"={
                samId3=samId3[which(clin3$disOntTerm3[samId3]%in%c("lungSmall","combined"))]
                heading=paste(heading,sep="")
            },
            "_noLung"={
                samId3=samId3[-grep("lung",clin3$disOntTerm3[samId3])]
                heading=paste(heading,", no lung",sep="")
            },
            "_fiona"={
                samId3=samId3[which(!is.na(clin32$fionaCol1[samId3]))]
                heading=paste(heading,", checked by Fiona",sep="")
            },
            "_unknown"={
                samId3=samId3[which(clin3$disOntTerm3[samId3]=="unknown")]
                heading=paste(heading,", unknown primary",sep="")
                heading=paste("Unknown primary",sep="")
            },
            "_tmbHigh"={
                samId3=samId3[which(clin32$tmbLevel[samId3]=="High")]
                heading=paste(heading,", high TMB",sep="")
                heading=paste("High TMB",sep="")
            }
            )
            
            switch(genesetFlag,
            "_10percPatient"={
                fNameOut=paste(fNameOut,genesetFlag,sep="")
                k1=grep("anyAlt_",colnames(datGD))
                nm=sub("anyAlt_","",colnames(datGD)[k1])
                k1=k1[which(nm%in%clin3$disOntTerm2[samId3])]
                i=c()
                propThres=0.10
                for (k in grep("anyAlt_",colnames(datGD))) {
                    #if (grep("_12genes",fName1)) {
                    i=c(i,which(round(datGD[,k],2)>propThres))
                    #} else {
                    #i=c(i,which(round(datGD[,k],2)>=propThres))
                    #}
                }
                i=sort(unique(i))
                heading=paste(heading,", genes with alterations in > 10% of patients in any group",sep="")
            },
            "_signif"={
                fName=paste("proportionTable_geneByDisease_qv",pThres,fName1,subsetName,".txt",sep="")
                heading="Table of significant genes in any of the comparisons pancreas vs. SCLC, other vs. SCLC or colon vs. SCLC"
                k1=grep("qv_anyAlt_",colnames(out1))
                k1=k1[grep("lungSmall",colnames(out1)[k1])]
                k1=k1[-grep("combined",colnames(out1)[k1])]
                nm1=colnames(out1)[k1]
                i=c()
                for (k in k1) {
                    i=c(i,which(out1[,k]<pThres))
                }
                i=sort(unique(i))
            },
            "_targetableGene"={
                fName=paste("proportionTable_geneByDisease_targetGene",fName1,subsetName,".txt",sep="")
                heading="Table of targetable genes"
                k1=grep("qv_anyAlt_",colnames(out1))
                k1=k1[grep("lungSmall",colnames(out1)[k1])]
                k1=k1[-grep("combined",colnames(out1)[k1])]
                nm1=colnames(out1)[k1]
                i=match(toupper(candGene$gene),toupper(rownames(out1)))
                i=i[which(!is.na(i))]
            },
            "_mmr"={
                fNameOut=paste(fNameOut,genesetFlag,sep="")
                tbl=read.table(paste("docs/candGene.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                tbl=tbl[which(tbl$family=="Mismatch excision repair (MMR)"),]
                x=unique(c(tbl$gene,unlist(sapply(tbl$alias,function(x) {strsplit(x,", ")[[1]]},USE.NAMES=F))))
                i=which(toupper(rownames(datGD))%in%x)
                heading=paste(heading,", MMR genes",sep="")
            },
            "_notch"={
                fNameOut=paste(fNameOut,genesetFlag,sep="")
                i=grep("NOTCH",toupper(rownames(datGP)))
                heading=paste(heading,", NOTCH genes",sep="")
            },
            "_brca"={
                fNameOut=paste(fNameOut,genesetFlag,sep="")
                x=c("BRCA1","BRCA2","ATM","ATR","PTEN","PARP1","MRE11A","PALB2","CDK12","BAP1","PRKDC","BARD1","BACH1","FANCD2","FANCE","CHEK1","CHEK2","PMS2","FANCA","FANCC","FANCF","XPA","NBN","BLM","BRIP1","RAD51B","RAD51C","RAD51D","RAD50")
                i=which(toupper(rownames(datGD))%in%x)
                heading=paste(heading,", BRCA genes",sep="")
            },
            "_setd2"={
                fNameOut=paste(fNameOut,genesetFlag,sep="")
                i=1:nrow(datGD)
                heading=paste(heading,", genes with alt, SETD2 marked",sep="")
            },
            "_allGenes"={
                fNameOut=paste(fNameOut,genesetFlag,sep="")
                i=1:nrow(datGD)
                heading=paste(heading,", all genes",sep="")
            }
            )
            
            colHM=colHMCon
            switch(datTypeFlag,
            "_prop"={
                arrayData=propMat[,grep("propMut_",colnames(propMat))]
                
                x1=apply(arrayData,1,function(x) mean(!is.na(x)))
                x2=apply(arrayData,2,function(x) mean(!is.na(x)))
                j=which(x2>=.5)
                arrayData=arrayData[,j]
                annRow=annCol=NULL
                
                x=apply(arrayData,1,function(x) sum(!duplicated(x[!is.na(x)])))
                i=which(x>1)
                arrayData=arrayData[i,]
                colnames(arrayData)=sub("carcinoma","cancer",sub("undifferentiated","ud",sub("primitive neuroectoderm tumor (PNET)","PNET",sub("neuroendocrine tumor","NET",sub("neuroendocrine carcinoma","NEC",sub("propMut_","",colnames(arrayData)))),fixed=T)))
            },
            "_binary"={
                distMethod="euclidean"
                arrayData=matrix(0,nrow=nrow(propMat),ncol=nrow(clin3),dimnames=list(rownames(propMat),clin3$id))
                tmpC=rep("",ncol(arrayData))
                annCol=data.frame(id=colnames(arrayData),disOntTerm=tmpC,stringsAsFactors=F)
                for (i in 1:nrow(arrayData)) {
                    k=which(clin1$gene==rownames(arrayData)[i] & !is.na(clin1$altType))
                    j=match(colnames(arrayData),clin1$id[k]); j2=which(!is.na(j)); j1=j[j2]
                    arrayData[i,j2]=1
                    annCol$disOntTerm[j2]=clin1$disOntTerm[k][j1]
                }
                annRow=NULL
                
                varList="disOntTerm"
                grpUniq=sort(unique(annCol$disOntTerm))
                nm=names(annCol)
                nm2=c()
                k1=round(seq(65,90,length=5)); k2=c(k1[2:(length(k1)-1)]-1,k1[length(k1)]); k1=k1[1:(length(k1)-1)]
                for (k in 1:length(k1)) {
                    x=annCol$disOntTerm
                    x[which(!substr(annCol$disOntTerm,1,1)%in%strsplit(chr(k1[k]:k2[k]),"")[[1]])]=NA
                    if (any(!is.na(x))) {
                        annCol=cbind(annCol,x)
                        nm2=c(nm2,paste("disOntTerm_",paste(chr(k1[k]),chr(k2[k]),sep="_"),sep=""))
                    }
                }
                names(annCol)=c(nm,nm2)
                varList=c(varList,nm2)
                varName=paste(varList," ",sep="")
            },
            "_patient"={
                #distMethod="euclidean"
                #distMethod="cosine"
                #arrayData=datGP[i,grep("anyAlt_",colnames(datGP))]
                #colnames(arrayData)=sub("anyAlt_","",colnames(arrayData))
                altTypeUniq1=sort(unique(clin1$altType))
                arrayData=datGP[i,grep(altTypeUniq1[1],colnames(datGP))]
                colnames(arrayData)=sub(paste(altTypeUniq1[1],"_",sep=""),"",colnames(arrayData))
                if (length(altTypeUniq1)==1) {
                    
                } else {
                    for (a1 in 2:length(altTypeUniq1)) {
                        arrayData[datGP[i,grep(altTypeUniq1[a1],colnames(datGP))]==1]=a1
                    }
                }
                #j=match(colnames(arrayData),clin3$id); j1=which(!is.na(j)); j2=j[j1]
                j=match(colnames(arrayData),clin3$id[samId3]); j1=which(!is.na(j)); j2=samId3[j[j1]]
                #annCol=data.frame(id=colnames(arrayData),disOntTerm2=clin3$disOntTerm2[j],stringsAsFactors=F)
                arrayData=arrayData[,j1]
                annCol=clin3[j2,which(!names(clin3)%in%colId2)]
                if (exists("clin32")) {
                    annCol=cbind(annCol,clin32[j2,which(!names(clin32)%in%names(annCol))])
                }
                annRow=NULL
                annRow=data.frame(id=rownames(arrayData),gene=rownames(arrayData),stringsAsFactors=F)
                
                #annCol$disOntTerm2[which(annCol$disOntTerm2=="other" & !annCol$disOntTerm%in%c("Duodenum neuroendocrine tumor","Esophagus neuroendocrine carcinoma","Small intestine neuroendocrine carcinoma","Stomach neuroendocrine carcinoma"))]=NA
                #annCol$disOntTerm2[which(annCol$disOntTerm=="Unknown primary undifferentiated neuroendocrine carcinoma")]="Unknown primary undifferentiated neuroendocrine carcinoma"
                #annCol$disOntTerm2[which(!annCol$disOntTerm3%in%c("combined","unknown"))]=NA
                
                if (exists("disDistMat")) {
                    j=match(annCol$id,rownames(disDistMat)); j1=which(!is.na(j)); j2=j[j1]
                    tmp=matrix(nrow=nrow(annCol),ncol=ncol(disDistMat),dimnames=list(annCol$id,paste("dist2",colnames(disDistMat),"_euclidean",sep="")))
                    tmp[j1,]=disDistMat[j2,]
                    annCol=cbind(annCol,tmp)
                    tmp=matrix(nrow=nrow(annCol),ncol=ncol(disDistMat),dimnames=list(annCol$id,paste("dist2",colnames(disDistMat),"_rank",sep="")))
                    tmp[j1,]=disDistMat[j2,]
                    tmp[j1,]=t(apply(disDistMat[j2,],1,function(x) {
                        y=order(order(x))
                        y[is.na(x)]=NA
                        y
                    }))
                    annCol=cbind(annCol,tmp)
                }
                
                if (subsetFlag=="_allAssays") {
                    varList=c("disOntTerm2","grade","gender","assayVersion")
                    varName=paste(c("disOntTerm","grade","gender","assay")," ",sep="")
                } else {
                    varList="disOntTerm2"
                    varName=paste("disOntTerm ",sep="")
                    varList=c("disOntTerm2","grade","gender","assayVersion")
                    varName=paste(c("disOntTerm","grade","gender","assay")," ",sep="")
                }
                if ("fionaCol1"%in%names(annCol)) {
                    varList=c(varList,"fionaCol1")
                    varName=c(varName,"vetted ")
                }
                x=names(annCol)[grep("_rank",names(annCol))]
                if (length(x!=0)) {
                    varList=c(varList,x)
                    varName=c(varName,paste(x," ",sep=""))
                }
                colHM=colHMCat
                limit1[2]=max(c(arrayData),na.rm=T)
            }
            )
            
            i=1:nrow(annRow)
            j=1:nrow(annCol)
            switch(subset2Flag,
            "_giUnknown"={
                j=which(annCol$disOntTerm3%in%c("combined","unknown"))
            },
            "_noLung"={
                j=(1:nrow(annCol))[-grep("Lung",annCol$disOntTerm)]
            },
            "_unknown"={
                j=which(annCol$disOntTerm3=="unknown")
            }
            )
            switch(genesetFlag,
            "_setd2"={
                #j=which(apply(arrayData,2,function(x) {mean(x!=0,na.rm=T)})!=0)
                i=which(annRow$gene=="SETD2")
                j=which(arrayData[i,]!=0)
                i=which(apply(arrayData[,j],1,function(x) {mean(x!=0,na.rm=T)})!=0)
            }
            )
            if (length(i)==0 | length(j)==0) {
                cat("No genes or samples !!!\n")
                next
            }
            arrayData=arrayData[i,j]
            annRow=annRow[i,]
            annCol=annCol[j,]
            
            if (genesetFlag%in%c("_setd2","_notch","_mmr","_brca","_allGenes")) {
                if (subset3Flag[1]=="_no0AltGene") {
                    heading=paste(heading,", genes with > 0 alt",sep="")
                    heading=sub("genes, genes","genes",heading)
                    i=apply(arrayData,1,function(x) {mean(x!=0,na.rm=T)!=0})
                    arrayData=arrayData[i,]
                    if (!is.null(annRow)) annRow=annRow[i,]
                }
                heading=paste(heading,", patients with > 0 alt",sep="")
                j=apply(arrayData,2,function(x) {mean(x!=0,na.rm=T)!=0})
                if (sum(j)<5) j=rep(T,ncol(arrayData))
                if (sum(j)==1) {
                    arrayData=matrix(arrayData[,j],ncol=1)
                    annCol=matrix(annCol[j,],nrow=1)
                } else {
                    arrayData=arrayData[,j]
                    annCol=annCol[j,]
                }
            }
            if (F) {
                switch(altTypeFlag,
                "_4Alt"={
                    heading=paste(heading,", 4 alt types",sep="")
                },
                "_anyAlt"={
                    heading=paste(heading,", alt types 0, 1",sep="")
                },
                "_anyAlt1minus1"={
                    heading=paste(heading,", alt types -1 (for no alt), 1",sep="")
                }
                )
            }
            
            arrayData2=arrayData
            if (altTypeFlag=="_anyAlt") {
                arrayData2[arrayData2!=0]=1
            }
            if (altTypeFlag=="_anyAlt1minus1") {
                arrayData2[arrayData2!=0]=1; arrayData2[arrayData2==0]=-1
            }
            arrayData[arrayData==0]=NA
            
            if (!is.na(nClust[1]) && nrow(arrayData)<nClust[1]) {
                nClust[1]=2
            }
            if (!is.na(nClust[2]) && ncol(arrayData)<20) {
                nClust[2]=NA
            }
            
            
            annRowAll=annRow
            #annColAll=annCol
            varListAll=varList
            varNameAll=varName
            varFListAll=varFList
            varFNameAll=varFName

            ## -------------------
            #if (sampleBar=="cluster") {
            if (sampleBar=="cluster" & nrow(arrayData2)>1 & ncol(arrayData2)>1) {
                switch(distMethod,
                "pearson"={
                    distMat=as.dist(1 - cor(t(arrayData2),method=distMethod,use="complete.obs"))
                },
                "spearman"={
                    distMat=as.dist(1 - cor(t(arrayData2),method=distMethod,use="complete.obs"))
                },
                "kendall"={
                    if (getCorFlag) {
                        corMatSam2=cor(t(arrayData2),method=distMethod,use="complete.obs")
                        save(corMatSam2,file=paste("corMatSam",fNameOut,".RData",sep=""))
                    } else {
                        load(file=paste(datadir,"corMatSam",fNameOut,".RData",sep=""))
                    }
                    distMat=as.dist(1 - corMatSam2)
                },
                "euclidean"={
                    distMat=dist(t(arrayData2), method=distMethod)
                },
                "cosine"={
                    distMat=getCosineDist(t(arrayData2))
                }
                )
                clustC=hclust(distMat, method=linkMethod)
            } else {
                clustC=NA
                nClust[2]=NA
            }
            
            #if (geneBar=="clusterPr") {
            if (geneBar=="clusterPr" & nrow(arrayData2)>1 & ncol(arrayData2)>1) {
                switch(distMethod,
                "pearson"={
                    distMat=as.dist(1 - cor(arrayData2,method=distMethod,use="complete.obs"))
                },
                "spearman"={
                    distMat=as.dist(1 - cor(arrayData2,method=distMethod,use="complete.obs"))
                },
                "kendall"={
                    if (getCorFlag) {
                        corMat2=cor(arrayData2,method=distMethod,use="complete.obs")
                        save(corMat2,file=paste("corMat",fNameOut,".RData",sep=""))
                    } else {
                        load(file=paste(datadir,"corMat",fNameOut,".RData",sep=""))
                    }
                    distMat=as.dist(1 - corMat2)
                },
                "euclidean"={
                    distMat=dist(arrayData2, method=distMethod)
                },
                "cosine"={
                    distMat=getCosineDist(arrayData2)
                }
                )
                clustR=try(hclust(distMat, method=linkMethod))
                if (class(clustR)=="try-error") {
                    clustR=NA
                    nClust[1]=NA
                }
            } else {
                clustR=NA
                nClust[1]=NA
            }
            
            ## -------------------
            if (nrow(arrayData)<=20) nameRow=rownames(arrayData) else nameRow=rep("",nrow(arrayData))
            if (nrow(arrayData)<=100) nameRow=rownames(arrayData) else nameRow=rep("",nrow(arrayData))
            switch(genesetFlag,
            "_setd2"={
                i=which(annRow$gene=="SETD2")
                nameRow[i]=paste(annRow$gene[i]," ",paste(rep("*",10),collapse=""),sep="")
            }
            )
            colRow=NULL
            if (is.null(varFList)) {
                colRow=NULL
            } else {
                colRow=matrix(nrow=length(varFList),ncol=nrow(annRow))
                for (varId in 1:length(varFList)) {
                    if (varFList[varId]%in%c("sd")) {
                        j=match(annRow$feature,annRowAll$feature)
                        x=round(100*annRowAll[,varFList[varId]])+1
                        lim=range(x,na.rm=T)
                        x[x<lim[1]]=lim[1]; x[x>lim[2]]=lim[2]
                        grpUniq=lim[1]:lim[2]
                        colRowUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                        colRow[varId,]=colRowUniq[x[j]]
                    } else {
                        x=as.character(annRowAll[,varFList[varId]])
                        x[x==""]=NA; x=as.integer(as.factor(x))
                        grpUniq=sort(unique(x))
                        x=x[match(annRow$feature,annRowAll$feature)]
                        if (length(grpUniq)<=length(colList2)) {
                            colRow[varId,]=colList2[x]
                        } else if (length(grpUniq)<=length(colList)) {
                            colRow[varId,]=colList[x]
                        } else {
                            colRow[varId,]=rainbow(length(grpUniq))[x]
                        }
                    }
                }
                rownames(colRow)=varFName
            }
            
            ## -------------------
            if (datTypeFlag=="_prop") nameCol=colnames(arrayData) else nameCol=rep("",ncol(arrayData))
            if (ncol(arrayData)<6) nameCol=colnames(arrayData)
            if (is.null(varList)) {
                colCol=NULL
            } else {
                colCol=matrix(nrow=length(varList),ncol=nrow(annCol))
                for (varId in 1:length(varList)) {
                    if (varList[varId]%in%c("sd")) {
                        j=match(annCol$id,annColAll$id)
                        x=round(100*annColAll[,varList[varId]])+1
                        lim=range(x,na.rm=T)
                        if (varList[varId]==c("sd")) lim=limSdSam
                        x[x<lim[1]]=lim[1]; x[x>lim[2]]=lim[2]
                        grpUniq=lim[1]:lim[2]
                        colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                        colCol[varId,]=colColUniq[x[j]]
                    } else if (length(grep("dist2",varList[varId]))==1) {
                        j=match(annCol$id,annColAll$id)
                        x=round(annColAll[,varList[varId]])
                        lim=range(x,na.rm=T)
                        grpUniq=lim[1]:lim[2]
                        colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                        colCol[varId,]=colColUniq[x[j]]
                    } else {
                        x=as.character(annColAll[,varList[varId]])
                        x[x==""]=NA; x=as.integer(as.factor(x))
                        grpUniq=sort(unique(x))
                        x=x[match(annCol$id,annColAll$id)]
                        if (length(grpUniq)<=length(colList2)) {
                            colCol[varId,]=colList2[x]
                        } else if (length(grpUniq)<=length(colList)) {
                            colCol[varId,]=colList[x]
                        } else {
                            colCol[varId,]=rainbow(length(grpUniq))[x]
                        }
                    }
                }
                rownames(colCol)=varName
            }
            
            ## -------------------
            if (F) {
                print("summary(c(arrayData))")
                print(summary(c(arrayData)))
                print("quantile(abs(c(arrayData)),probs=seq(0,1,by=.1),na.rm=T)")
                print(quantile(abs(c(arrayData)),probs=seq(0,1,by=.1),na.rm=T))
            }
            print("table(c(arrayData2))")
            print(table(c(arrayData2)))
            main=heading
            
            if (F) {
                subDir <- paste(fNameOut,sep="")
                if (!file.exists(subDir)){
                    dir.create(file.path(subDir))
                }
                subDir=paste(subDir,"/",sep="")
            }
            subDir=""
            if (outFormat=="png") {
                margins=c(5,20)
                if (datTypeFlag=="_prop") margins=c(20,1) else margins=c(1,1)
                margins[2]=5
                png(paste(subDir,"heatmap",fNameOut,".png",sep=""),width=480*2,height=480*2)
            } else {
                margins=c(12,5)
                pdf(paste(subDir,"heatmap",fNameOut,".pdf",sep=""))
                par(cex.main=0.7)
            }
            hcc=heatmap3(x=arrayData, Rowv=clustR, Colv=clustC, distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=colCol, RowSideColors=colRow, labCol=nameCol, labRow=nameRow, ncr=nClust[1], ncc=nClust[2], scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit1, cexCol=2, , high=colHM[[1]], low=colHM[[2]], mid=colHM[[3]])
            dev.off()
            
            
            #xxx1=list(x=arrayData, Rowv=clustR, Colv=clustC,ColSideColors=colCol, RowSideColors=colRow, labCol=nameCol, labRow=nameRow)
            
            
            ## -------------------
            #if (!is.na(nClust[1])) {
            fNameThis=paste(subDir,"clusterInfoFeature",fNameOut,".txt",sep="")
            
            x=paste("Gene cluster information for heatmap of samples from\n",heading,"\nGenes are ordered as in heatmap from bottom to top. See columns 'clustId' and 'order'",sep="")
            write.table(x, fNameThis, sep="\t", col.names=F, row.names=F, quote=F)
            
            if (F) {
                pdf(paste(subDir,"clusterFeatures",fNameOut,".pdf",sep=""))
                plot(clustR,main=paste("Feature clusters with ",nClust[1]," main clusters marked in red",sep=""),xlab="",sub="",ylab=NULL,axes=F, cex=.2); rect.hclust(clustR,k=nClust[1])
                dev.off()
            }
            
            nClustThis=ifelse(is.na(nClust[1]),1,nClust[1])
            clustId=cutree(clustR,k=nClustThis)[clustR$order]
            k1=which(!duplicated(clustId))
            for (k in 1:length(k1)) {
                clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
            }
            
            tbl=cbind(annRow[clustR$order,],clustId,order=1:nrow(annRow))
            
            write.table(tbl, fNameThis, sep="\t", col.names=T, row.names=F, quote=F,append=T)
            #}
            
            #if (!is.na(nClust[2])) {
            fNameThis=paste(subDir,"clusterInfoSample",fNameOut,".txt",sep="")
            
            x=paste("Sample cluster information for heatmap of samples from\n",heading,"\nSamples are ordered as in heatmap from left to right. See columns 'clustId' and 'order'",sep="")
            write.table(x, fNameThis, sep="\t", col.names=F, row.names=F, quote=F)
            
            tbl=annCol
            #write.table("Global",file=fNameThis,append=F,col.names=F,row.names=F,sep="\t",quote=F)
            #write.table(heading,file=fNameThis,append=T,col.names=F,row.names=F,sep="\t",quote=F)
            write.table("",file=fNameThis,append=T,col.names=F,row.names=F,sep="\t",quote=F)
            write.table(paste("No. of patients",nrow(tbl),sep="\t"),file=fNameThis,append=T,col.names=F,row.names=F,sep="\t",quote=F)
            write.table("",file=fNameThis,append=T,col.names=F,row.names=F,sep="\t",quote=F)
            tbl$ageBi=as.integer(tbl$age<=median(tbl$age,na.rm=T))
            res=matrix(0,nrow=nrow(tbl),ncol=length(altTypeUniq1))
            colnames(res)=altTypeUniq1
            for (k in 1:length(altTypeUniq1)) {
                res[apply(arrayData2,2,function(x){any(x==k,na.rm=T)}),k]=1
            }
            tbl=cbind(tbl,anyAlt=apply(arrayData2,2,function(x) {as.integer(any(x!=0,na.rm=T))}),res)
            res2=apply(arrayData2,2,function(x) {sum(x!=0,na.rm=T)})
            summary(res2)
            write.table(paste("Average of ",round(mean(res2),2)," alterations (range ",min(res2,na.rm=T),", ",max(res2,na.rm=T),") per patient",sep=""),file=fNameThis,append=T,col.names=F,row.names=F,sep="\t",quote=F)
            write.table("",file=fNameThis,append=T,col.names=F,row.names=F,sep="\t",quote=F)
            #res[res>0]=1
            #altTypeUniq1=sort(unique(tbl$altType))
            #tmp=matrix(nrow=nrow(tbl),ncol=length(altTypeUniq1)+1)
            #colnames(tmp)=c("anyAlt",altTypeUniq1)
            #tbl=cbind(tbl,tmp)
            #jj=match(rownames(res),tbl$id)
            #tbl$anyAlt[jj]=apply(res,1,function(x) as.integer(any(x==1,na.rm=T)))
            #k=match(colnames(res),altTypeUniq1)
            #tbl[jj,altTypeUniq1[k]]=res
            #jj=which(!duplicated(paste(tbl$gene,tbl$varStatus)))
            write.table(paste("Median age",round(median(tbl$age,na.rm=T),2),sep="\t"),file=fNameThis,append=T,col.names=F,row.names=F,sep="\t",quote=F)
            write.table("",file=fNameThis,append=T,col.names=F,row.names=F,sep="\t",quote=F)
            varList2=c("gender","grade","disOntTerm","disOntTerm2","assayVersion","anyAlt",sort(unique(tbl$altType)))
            varList2=c("gender","disOntTerm2","assayVersion","anyAlt",sort(unique(tbl$altType)))
            varList2=c("gender","disOntTerm2","assayVersion","anyAlt",altTypeUniq1)
            for (v1 in 1:length(varList2)) {
                res=table(tbl[,varList2[v1]])
                tbl2=cbind(names(res),res,round(res/sum(res),2))
                colnames(tbl2)=c(varList2[v1],"count","proportion")
                if (any(c("0","1")%in%tbl2[,1],na.rm=T)) {
                    tbl2[,1][which(tbl2[,1]=="0")]="no mutation"
                    tbl2[,1][which(tbl2[,1]=="1")]="mutation"
                }
                write.table("\tPatient",file=fNameThis,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                write.table(tbl2,file=fNameThis,append=T,col.names=T,row.names=F,sep="\t",quote=F)
                write.table("",file=fNameThis,append=T,col.names=F,row.names=F,sep="\t",quote=F)
            }
            
            if (F) {
                pdf(paste(subDir,"clusterSamples",fNameOut,".pdf",sep=""))
                plot(clustC,main=paste("Sample clusters with ",nClust[2]," main clusters marked in red",sep=""),xlab="",sub="",ylab=NULL,axes=F, cex=.2); rect.hclust(clustC,k=nClust[2])
                dev.off()
            }
            
            nClustThis=ifelse(is.na(nClust[2]),1,nClust[2])
            clustId=cutree(clustC,k=nClustThis)[clustC$order]
            k1=which(!duplicated(clustId))
            for (k in 1:length(k1)) {
                clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
            }
            
            tbl=cbind(annCol[clustC$order,],clustId,order=1:nrow(annCol))
            write.table(tbl, fNameThis, sep="\t", col.names=T, row.names=F, quote=F,append=T)
            #}
            
            if (!is.null(colRow)) {
                for (varId in 1:length(varFListAll)) {
                    if (varFListAll[varId]%in%c("sd")) {
                        width = 480; height = 140
                    } else {
                        width = 480; height = 480
                    }
                    if (outFormat=="png") {
                        png(paste("heatmapFeatureColorBarLegend_",varFListAll[varId],".png",sep=""),width=width,height=height)
                    } else {
                        pdf(paste("heatmapFeatureColorBarLegend_",varFListAll[varId],".pdf",sep=""))
                    }
                    if (varFListAll[varId]%in%c("sd")) {
                        x=round(100*annRowAll[,varFListAll[varId]])+1
                        lim=range(x,na.rm=T)
                        grpUniq=lim[1]:lim[2]
                        colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                        lim=(lim-1)/100
                        heatmapColorBar(limit=lim,cols=c(colColUniq[c(length(colColUniq),1,median(1:length(colColUniq)))]))
                    } else {
                        x=as.character(annRowAll[,varFListAll[varId]]); x[x==""]=NA
                        grpUniq=table(x)
                        grpUniq=names(grpUniq)
                        k=1:length(grpUniq)
                        if (length(grpUniq)<=length(colList2)) {
                            sampleColorLegend(tls=grpUniq[k],col=colList2,legendTitle=varFNameAll[varId])
                        } else if (length(grpUniq)<=length(colList)) {
                            cexThis=NULL
                            if (length(grpUniq)>10) cexThis=1
                            sampleColorLegend(tls=grpUniq[k],col=colList,legendTitle=varFNameAll[varId],cex=cexThis)
                        } else {
                            sampleColorLegend(tls=grpUniq[k],col=rainbow(length(grpUniq)),legendTitle=varFNameAll[varId])
                        }
                    }
                    dev.off()
                }
            }
            
            
            ## -------------------
            if (!is.null(colCol)) {
                for (varId in 1:length(varListAll)) {
                    if (length(grep("dist2",varList[varId]))==1) {
                        if (outFormat=="png") {
                            png(paste("heatmapSampleColorBarLegend_",varListAll[varId],fNameOut2,".png",sep=""),width=480,height=140)
                        } else {
                            pdf(paste("heatmapSampleColorBarLegend_",varListAll[varId],fNameOut2,".pdf",sep=""))
                        }
                        x=round(annColAll[,varListAll[varId]])+1
                        lim=range(x,na.rm=T)
                        #if (length(grep("dist2class",varList[varId]))==1) lim=limDist2classSam
                        grpUniq=lim[1]:lim[2]
                        colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                        lim=lim-1
                        heatmapColorBar(limit=lim,cols=c(colColUniq[c(length(colColUniq),1,median(1:length(colColUniq)))]))
                    } else {
                        if (varListAll[varId]=="disOntTerm2") {
                            #x=as.character(annColAll[samId03,varListAll[varId]]); x[x==""]=NA
                            x=as.character(annColAll[,varListAll[varId]]); x[x==""]=NA
                        } else {
                            x=as.character(annColAll[,varListAll[varId]]); x[x==""]=NA
                        }
                        grp=table(x)
                        grpUniq=names(grp)
                        k=1:length(grpUniq)
                        ttl=paste(grpUniq[k]," (",grp[k],")",sep="")
                        ttl=grpUniq[k]
                        if (length(grpUniq)<6) {
                            width = 480; height = 480
                        } else {
                            width = 560; height = 960
                        }
                        if (outFormat=="png") {
                            png(paste("heatmapSampleColorBarLegend_",varListAll[varId],fNameOut2,".png",sep=""),width=width,height=height)
                        } else {
                            pdf(paste("heatmapSampleColorBarLegend_",varListAll[varId],fNameOut2,".pdf",sep=""))
                        }
                        cexThis=NULL
                        if (outFormat=="pdf" & (length(grpUniq)>15 | max(nchar(grpUniq))>20)) cexThis=1
                        if (outFormat=="pdf") cexThis=1
                        if (length(grpUniq)<=length(colList2)) {
                            sampleColorLegend(tls=ttl,col=colList2,legendTitle=varNameAll[varId],cex=cexThis)
                        } else if (length(grpUniq)<=length(colList)) {
                            sampleColorLegend(tls=ttl,col=colList,legendTitle=varNameAll[varId],cex=cexThis)
                        } else {
                            sampleColorLegend(tls=ttl,col=rainbow(length(grpUniq)),legendTitle=varNameAll[varId],cex=cexThis)
                        }
                    }
                    dev.off()
                }
            }
            
            if (length(colHM[[1]])==1) {
                if (outFormat=="png") {
                    png(paste("heatmapColorRange",fNameOut2,".png",sep=""),width=480,height=140)
                } else {
                    pdf(paste("heatmapColorRange",fNameOut2,".pdf",sep=""))
                }
                heatmapColorBar=function(limit,cols=c("green","red","black")) {
                    try <- maPalette(high=cols[1], low=cols[2], mid=cols[3],k=15)
                    maColorBar(try, scale=limit,k=3)
                }
                heatmapColorBar(limit=limit1,cols=unlist(colHM))
                dev.off()
            } else {
                width = 480; height = 480
                if (outFormat=="png") {
                    png(paste("heatmapColorLegend",fNameOut2,".png",sep=""),width=width,height=height)
                } else {
                    pdf(paste("heatmapColorLegend",fNameOut2,".pdf",sep=""))
                }
                grpUniq=sort(unique(clin1$altType))
                cexThis=NULL
                if (outFormat=="pdf") cexThis=1
                sampleColorLegend(tls=grpUniq,col=colHM[[1]],legendTitle="alteration type",cex=cexThis)
                dev.off()
                
            }
            
        }
    }
}

"

Skin sample in setd2 heatmap is closest to SCLC, skin is a close second
> cbind(clin3[samId03,],disDistMat)[which(rownames(disDistMat)%in%xxx),]
id tissue priTissue age gender       study                 disOntTerm qcResult  gene
6897 1690  Block      Skin  80   Male CLINICAL-T7 Skin Merkel cell carcinoma     Pass SPTA1
partnerGene       altType cdsEffect proteinEffect cn fracReads varStatus assayVersion grade
6897        <NA> short variant   3067C>T        Q1023* NA      0.59    likely           T7  <NA>
disOntTerm2 disOntTerm3    colon lungSmall pancreas     skin
6897        skin        skin 2.350162  2.185624 2.542175 2.229089


These are SCLC samples in setd2 heatmap which are not closest to SCLC centroid but it is a close second
> cbind(clin3[samId03,],disDistMat)[which(rownames(disDistMat)%in%annCol$id[which(annCol$disOntTerm2=="lungSmall" & annCol$dist2lung_rank!=1)]),]
    id tissue   priTissue age gender       study                                 disOntTerm
2832  774  Block  Lymph Node  67 Female CLINICAL-T7 Lung small cell undifferentiated carcinoma
3591  939  Block Soft Tissue  61 Female CLINICAL-T7 Lung small cell undifferentiated carcinoma
4725 1219  Slide       Liver  49 Female CLINICAL-T7 Lung small cell undifferentiated carcinoma
    qcResult  gene partnerGene       altType cdsEffect proteinEffect cn fracReads varStatus
2832     Pass RBM10        <NA> short variant   1753C>T         Q585* NA      0.82    likely
3591     Pass   WT1       OPCML rearrangement      <NA>          <NA> NA        NA    likely
4725     Pass SETD2        <NA> short variant   3918G>A        W1306* NA      0.17    likely
    assayVersion grade disOntTerm2 disOntTerm3    colon lungSmall pancreas     skin
2832           T7  <NA>   lungSmall   lungSmall 2.813435  2.846475 2.938174 2.878718
3591           T7  <NA>   lungSmall   lungSmall 2.472144  2.485958 2.230088 2.335313
4725           T7  <NA>   lungSmall   lungSmall 1.348102  1.237234 1.374112 1.218190

"

## ----------------------------------------------
## Comparisons of diseases

library(coin)

clin4=clin1

altTypeUniq=c("",sort(clin4$altType[!duplicated(clin4$altType)])); altTypeUniq=altTypeUniq[!is.na(altTypeUniq)]

grpUniq=c("Lung large cell neuroendocrine carcinoma","Lung small cell undifferentiated carcinoma")
grpName=c("Large cell","Small cell")
grpName2=c("largeLung","smallLung")

compFlag=paste("_",grpName2[2],"V",grpName2[1],sep="")

table(clin3$disOntTerm)

clin4$grp=clin4$disOntTerm
for (k in 1:length(grpUniq)) clin4$grp[which(clin4$grp==grpUniq[k])]=grpName[k]
fName=compFlag
heading="Lung carcinoma: Small cell vs. large cell"

varList2=c("fracReads","cn")
varName2=c("fracReads","cn")

varList3=c("cdsEffect","proteinEffect")
varName3=c("cdsEffect","proteinEffect")



geneList=order(gene,decreasing=T); geneList=geneList[which(gene[geneList]>9)]
nm=c("pvFisher","pvBinom")
out=matrix(NA,nrow=length(geneList),ncol=(length(grpUniq)+length(nm))*length(altTypeUniq),dimnames=list(names(gene)[geneList],paste(c(paste("num_",grpName2,sep=""),paste(nm,compFlag,sep="")),"_",rep(altTypeUniq,each=length(grpUniq)+length(nm)),sep="")))
j0=1:nrow(clin4)
j1=j0[!duplicated(clin4$id[j0])]
x0=table(disOntTerm=as.factor(clin4$disOntTerm)[j1])
x0[names(x0)==grpUniq[2]]/(x0[names(x0)==grpUniq[1]]+x0[names(x0)==grpUniq[2]])
prob=x0[names(x0)==grpUniq[2]]/(x0[names(x0)==grpUniq[1]]+x0[names(x0)==grpUniq[2]])
# 0.7551515
j1=j0[which(is.na(clin4$gene[j0]))]
x2=table(disOntTerm=as.factor(clin4$disOntTerm)[j1])
varThis=as.integer(!is.na(clin4$altType))
for (g1 in geneList) {
    k1=match(names(gene)[g1],rownames(out))
    for (a1 in 1:length(altTypeUniq)) {
        if (altTypeUniq[a1]=="") {
            j11=j0[which(clin4$disOntTerm[j0]==grpUniq[1] & clin4$gene[j0]==names(gene)[g1])]
            j21=j0[which(clin4$disOntTerm[j0]==grpUniq[2] & clin4$gene[j0]==names(gene)[g1])]
            j1=j0[which(clin4$disOntTerm[j0]==grpUniq[1] & (varThis[j0]==0 | clin4$gene[j0]==names(gene)[g1]))]
            j2=j0[which(clin4$disOntTerm[j0]==grpUniq[2] & (varThis[j0]==0 | clin4$gene[j0]==names(gene)[g1]))]
        } else {
            j11=j0[which(clin4$disOntTerm[j0]==grpUniq[1] & clin4$gene[j0]==names(gene)[g1] & clin4$altType[j0]==altTypeUniq[a1])]
            j21=j0[which(clin4$disOntTerm[j0]==grpUniq[2] & clin4$gene[j0]==names(gene)[g1] & clin4$altType[j0]==altTypeUniq[a1])]
            j1=j0[which(clin4$disOntTerm[j0]==grpUniq[1] & (varThis[j0]==0 | (clin4$gene[j0]==names(gene)[g1] & clin4$altType[j0]==altTypeUniq[a1])))]
            j2=j0[which(clin4$disOntTerm[j0]==grpUniq[2] & (varThis[j0]==0 | (clin4$gene[j0]==names(gene)[g1] & clin4$altType[j0]==altTypeUniq[a1])))]
        }
        j11=j11[!duplicated(clin4$id[j11])]
        j21=j21[!duplicated(clin4$id[j21])]
        j1=j1[!duplicated(clin4$id[j1])]
        j2=j2[!duplicated(clin4$id[j2])]
        j=c(j1,j2)
        x1=table(varThis[j],clin4$grp[j])
        out[k1,match(paste("num_",grpName2,"_",altTypeUniq[a1],sep=""),colnames(out))]=c(length(j11),length(j21))
        if (length(j11)>4 & length(j21)>4) {
            pv=try(fisher.test(varThis[j],clin4$grp[j])$p.value)
            if (class(pv)=="try-error") {
                pv=NA
            }
            out[k1,match(paste("pvFisher",compFlag,"_",altTypeUniq[a1],sep=""),colnames(out))]=pv
            pv=try(binom.test(x=length(j21),n=length(j11)+length(j21),p=prob)$p.value)
            if (class(pv)=="try-error") {
                pv=NA
            } else if (class(pv)=="logical") {
                pv=2.2^-16
            }
            out[k1,match(paste("pvBinom",compFlag,"_",altTypeUniq[a1],sep=""),colnames(out))]=pv
        }
    }
}
colnames(out)=sub("_$","",colnames(out))
write.table(cbind(gene=rownames(out),out),file=paste("pValue",fName,".txt",sep=""),append=F,col.names=T,row.names=F,sep="\t",quote=F)
png("plot_pv_BinomVsFisher.png")
par(mfrow=c(3,3))
k1=grep("pvFisher",colnames(out))
k2=grep("pvBinom",colnames(out))
lim=c(0,1)
for (k in 1:length(k1)) {
    plot(out[,k1[k]],out[,k2[k]],xlim=lim,ylim=lim,main=sub("pvFisher_","",colnames(out)[k1[k]]),xlab="PV (Fisher's test)",ylab="PV (Binomial test)"); abline(c(0,1))
}
dev.off()

gene=table(clin4$gene)
gene[order(gene,decreasing=T)][1:10]
varList=c("altType")
k=1
out=NULL
geneList=order(gene,decreasing=T); geneList=geneList[which(gene[geneList]>9)]
for (g1 in geneList) {
    #j0=which(clin4$id==clin3$id[j3] & clin4$gene==names(gene)[g1])
    j0=1:nrow(clin4)
    for (a1 in 1:length(altTypeUniq)) {
        if (altTypeUniq[a1]=="") {
            j1=j0[which(clin4$gene[j0]==names(gene)[g1] & clin4$disOntTerm[j0]==grpUniq[1])]
            j2=j0[which(clin4$gene[j0]==names(gene)[g1] & clin4$disOntTerm[j0]==grpUniq[2])]
        } else {
            j1=j0[which(clin4$gene[j0]==names(gene)[g1] & clin4$altType[j0]==altTypeUniq[a1] & clin4$disOntTerm[j0]==grpUniq[1])]
            j2=j0[which(clin4$gene[j0]==names(gene)[g1] & clin4$altType[j0]==altTypeUniq[a1] & clin4$disOntTerm[j0]==grpUniq[2])]
        }
        j1=j1[!duplicated(clin4$id[j1])]
        j2=j2[!duplicated(clin4$id[j2])]
        if (length(j1)>4 & length(j2)>4) {
            j=c(j1,j2)
            for (vId2 in 1:length(varList2)) {
                if (sum(!is.na(clin4[j,varList2[vId2]]))>9) {
                    cat(g1,a1,"\n")
                    pv=wilcox.test(clin4[j,varList2[vId2]]~clin4$grp[j])$p.value
                    out=rbind(out,c(fName,varName2[vId2],names(gene)[g1],altTypeUniq[a1],"Wilcoxon test",pv))
                    ttl=paste(grpName," (",c(length(j1),length(j2)),")",sep="")
                    png(paste("boxplot_",names(gene)[g1],"_",altTypeUniq[a1],"_",varName2[vId2],fName,".png",sep=""))
                    boxplot(clin4[j,varList2[vId2]]~clin4$grp[j],names=ttl,main=paste(heading,": ",names(gene)[g1],", ",altTypeUniq[a1],"\npv ",signif(pv,2),sep=""),ylab=varName2[vId2])
                    dev.off()
                }
            }
            for (vId3 in 1:length(varList3)) {
                if (sum(!is.na(clin4[j,varList3[vId3]]))>9 & sum(!duplicated(clin4[j,varList3[vId3]]))<11) {
                    cat(g1,a1,"\n")
                    pv=try(fisher.test(clin4[j,varList3[vId3]],clin4$grp[j])$p.value)
                    nm="Fisher test"
                    if (class(pv)=="try-error") {
                        pv=try(chisq.test(clin4[j,varList3[vId3]],clin4$grp[j])$p.value)
                        nm="Chisq test"
                    }
                    out=rbind(out,c(fName,varName3[vId3],names(gene)[g1],altTypeUniq[a1],nm,pv))
                }
            }
        }
    }
}
colnames(out)=c("part","variable","gene","altType","test","pv")
write.table(out,file=paste("summary",fName,".txt",sep=""),append=F,col.names=T,row.names=F,sep="\t",quote=F)


RPTOR_short variant
j=c(4210,4483,5396,5515,6255,7783,1310,3346,4052,4329,4476,7541)
table(clin4[j,"proteinEffect"],clin4$grp[j])
            Large cell Small cell
A727fs*24           1          0
C1298Y              0          1
E1973fs*6           0          1
F519fs*14           0          1
H1904R              0          1
L532*               0          1
N1574fs*22          0          1
P1962L              0          2
Q1445*              0          1
S1431fs*47          1          0
T1066fs*16          0          1

proteinEffect:
"P1962L"
"Residue change:  From Proline (P) to Leucine (L) at position 1962 (P1962L, p.Pro1962Leu)"


####################################################################
####################################################################
