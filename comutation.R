## Co-mutation

# Actionable genes
# grpUniq=c("BRAF","SETD2","MMR","ATRX","BRCA1","BRCA2")
#
# Run genesetFlag="" & genesetFlag="_actionableGene"
# Run subset2Flag="_lungSmallJeffVettedHighGradeGI_T5aT7Assays" first then subset2Flag="_lungSmallFionaVettedGI_T5aT7Assays"

## ---------------------
## ---------------------
## Comparison within a disease

library(qvalue)
source(paste(dirSrc,"functions/heatmapAcgh.7.1.R",sep=""))

nameValueList=c("allSamples")
nameValueList=c("allSamples","jeffVetted","fionaVetted")
for (nameValueFlag in nameValueList) {
switch(nameValueFlag,
    "allSamples"={nameValue=data.frame(name="subset",value="allSamples")
    },
    "jeffVetted"={nameValue=data.frame(name="subset",value="jeffVetted")
    },
    "fionaVetted"={nameValue=data.frame(name="subset",value="fionaVetted")
    }
)
write.table(nameValue,"config.tmp", sep="=", col.names=F, row.names=F, quote=F)
k=which(nameValue$name=="subset")
datadir="results/"
if (nameValue$value[k]=="fionaVetted") {
    #load(paste(datadir,"tmp_fiona166_T5aT7Assays.RData",sep=""))
    load(paste(datadir,"tmp_allAssays.RData",sep=""))
    nameValue=read.table("config.tmp",sep="=",h=F,quote="",comment.char="",as.is=T,fill=T,col.names=c("name","value"))
    subset2List=c("_lungSmallFionaVettedGI_T5aT7Assays")
} else {
    load(paste(datadir,"tmp_allAssays.RData",sep=""))
    nameValue=read.table("config.tmp",sep="=",h=F,quote="",comment.char="",as.is=T,fill=T,col.names=c("name","value"))
    subset2List=c("_allSamples_T5aT7Assays")
    k=which(nameValue$name=="subset")
    if (nameValue$value[k]=="jeffVetted") {
        subset2List=c("_lungSmallJeffVettedHighGradeGI_T5aT7Assays")
    }
}

## -------------------
geneCorTh=0.9
geneCorTh=0
geneCorTh=0.5
pThres=0.001
pThres=0.01
pThres=0.05
pThres=99
samSize2=5
samSize2=10

testFlag="cosine"
testFlag="fisher"

propFlag=">="
propFlag=">"

outFormat="png"
outFormat="pdf"

source(paste(dirSrc3,"funcs.R",sep=""))
out=getFamilyLevelInfo()
datGP_m=out$datGP_m
candGene=out$candGene

geneFamilyFlag=F
geneFamilyFlag=T

#candGeneThis=candGene

## -------------------

confIntFlag=F
confIntFlag=T

samSize=c(0,0)
samSize=c(50,10)
samSize=c(50,30)

datadir2=""
geneFlag=F
geneFlag=T

genesetFlag=""; geneset2List=c("10%patient","15%patient","100topGene")
genesetFlag=""; geneset2List=c("10%patientGrp1Grp2")
genesetFlag=""; geneset2List=c("10%patientDisease")
genesetFlag=""; geneset2List=c("10%patient")
genesetFlag="_actionableGene"; geneset2List=c("actionableGene")
genesetFlag=""; geneset2List=c("")
genesetFlag="_swiSnfCompEtc"; geneset2List=c("swiSnfCompEtc")
genesetFlag="_swiSnfEtc"; geneset2List=c("swiSnfEtc")
genesetFlag="_swiSnfHisModAtmAtr"; geneset2List=c("swiSnfHistoneModifierAtmAtr")
    
#subset2Flag="_T5aT7Assays"

#subset2List=c("_lungSmallJeffVettedHighGradeGI_T5aT7Assays","_lungSmallFionaVettedGI_T5aT7Assays","_T5aT7Assays","_fionaVetted_T5aT7Assays","_withinSmallCell_T5aT7Assays","_lungSmallFionaVettedSmallCellGI_T5aT7Assays","_lungSmallSmallCellGI_T5aT7Assays","_lungSmallLargeCellGI_T5aT7Assays","_tmbLevel_T5aT7Assays","_tmbLevel_withinVetted_T5aT7Assays")
#subset2List=c("_lungSmallJeffVettedHighGradeGI_T5aT7Assays","_lungSmallFionaVettedGI_T5aT7Assays","_fionaVetted_T5aT7Assays","_lungSmallFionaVettedSmallCellGI_T5aT7Assays","_tmbLevel_withinVetted_T5aT7Assays")
#subset2List=sub("_T5aT7Assays","_noOtherGI_T5aT7Assays",subset2List)

if (F) {
    genesetFlag=""; geneset2List=c("100topGene"); subset2List=c("_fionaVetted_noOtherGI_T5aT7Assays_fionaVetted_noOtherGI_T5aT7Assays")
    geneFlag=F
    samSize=c(50,30)
}

if (F) {
    subset2List=c("_lungSmallFionaVettedGI_noOtherGI_T5aT7Assays","_fionaVetted_noOtherGI_T5aT7Assays","_lungSmallFionaVettedSmallCellGI_noOtherGI_T5aT7Assays")
    subset2List=c("_lungSmallFionaVettedSmallCellGI_noOtherGI_T5aT7Assays")
    subset2List=c("_lungSmallFionaVettedGI_noOtherGI_T5aT7Assays")
    geneFlag=F
    samSize=c(50,30)

    subset2List=c("_lungSmallFionaVettedGI_T5aT7Assays")
    geneFlag=F
    samSize=c(0,0)

    subset2List=c("_lungSmallJeffVettedHighGradeGI_T5aT7Assays")
    geneFlag=F
    samSize=c(0,0)

    subset2List=c("_lungSmall_T5aT7Assays")
    geneFlag=F
    samSize=c(0,0)
}

geneFlag=F
samSize=c(50,30)
samSize=c(0,0)

if (geneFamilyFlag) {
    datGPThis=datGP_m
} else {
    datGPThis=datGP
}

candGeneThis=getCandidateGenes(genesetFlag)

if (genesetFlag!="") {
    if (geneFamilyFlag) {
        i=match(toupper(unique(candGeneThis$family2)),toupper(rownames(datGPThis)))
    } else {
        i=match(toupper(candGeneThis$gene),toupper(rownames(datGPThis)))
    }
    datGPThis=datGPThis[i,]
}

datGPName=rep("",nrow(datGPThis)*(nrow(datGPThis)-1)/2)
k=1
for (i1 in 1:(nrow(datGPThis)-1)) {
    for (i2 in (i1+1):nrow(datGPThis)) {
        datGPName[k]=paste(rownames(datGPThis)[i1],rownames(datGPThis)[i2],sep="_")
        k=k+1
    }
}


for (subset2Flag in subset2List) {
    fName1=paste(genesetFlag,subset2Flag,sep="")
    altTypeUniq1=sort(unique(clin1$altType))
    #altTypeUniq2=cbind(c(altTypeUniq1,paste(altTypeUniq1[1],"+",altTypeUniq1[3]),paste(altTypeUniq1[1:3],"+",altTypeUniq1[4])),as.character(c(1:length(altTypeUniq1),10*(1)+3,10*(1:3)+4)))
    altTypeUniq2=cbind(c(altTypeUniq1,paste(altTypeUniq1[1],"+",altTypeUniq1[2]),paste(altTypeUniq1[1],"+",altTypeUniq1[3]),paste(altTypeUniq1[1:3],"+",altTypeUniq1[4]),paste(altTypeUniq1[1],"+",altTypeUniq1[2],"+",altTypeUniq1[4])),as.character(c(1:length(altTypeUniq1),10*(1)+2:3,10*(1:3)+4,124)))
    altTypeUniq2=cbind(c(altTypeUniq1,"multiple"),as.character(c(1:length(altTypeUniq1),10)))
    altTypeUniq3="anyAlt"
    
    if ("multiple"%in%altTypeUniq2[,1]) {
        colList=c("red","blue","orange","cyan","yellow","indianred","yellow2","skyblue","bisque3","indianred4")
    } else {
        colList=c("red","blue","purple","cyan","skyblue","orange","yellow","darkgreen")
        colList=c("red","blue","indianred","cyan","orange","yellow","skyblue","bisque3")
        colList=c("red","blue","indianred","cyan","yellow2","orange","yellow","skyblue","bisque3","indianred4")
    }
    
    gene1=unique(clin1$gene[clin1$assayVersion=="T5a"])
    gene2=unique(clin1$gene[clin1$assayVersion=="T7"])
    if (subset2Flag%in%c("_allSamples_T5aT7Assays","_lungSmallJeffVettedHighGradeGI_T5aT7Assays","_lungSmall_T5aT7Assays")) {
        j=1:nrow(clin3)
    } else {
        j=which(!is.na(clin32$fionaCol1))
    }
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
    
    if (subset2Flag=="_lungSmallFionaVettedGI_noOtherGI_T5aT7Assays") {
        clinThis=clinThis[which(clinThis$disOntTerm3=="SCLC" | clinThis$fionaCol1=="fionaY"),]
        heading1=paste("Fiona vetted, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
        disUniq=cbind(c("lungSmall","pancreas","colon"),c("SCLC","Pancreas","Colon"))
        varList=c("disease")
        varName=c("disease")
        geneCorList=list(lungSmall=NULL,pancreas=NULL,colon=NULL)
    } else if (subset2Flag%in%c("_lungSmallFionaVettedGI_T5aT7Assays")) {
        clinThis=clinThis[which(clinThis$disOntTerm3=="SCLC" | clinThis$fionaCol1=="fionaY"),]
        heading1=paste("Fiona vetted, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
        disUniq=cbind(c("lungSmall","pancreas","colon","other"),c("SCLC","Pancreas","Colon","Other GI"))
        disUniq=cbind(c("lungSmall","pancreas","colon","other"),c("SCLC","Pancreas","Colon","Other GI"))
        varList=c("disease")
        varName=c("disease")
        geneCorList=list(lungSmall=NULL,pancreas=NULL,colon=NULL)
    } else if (length(grep("_allSamples_T5aT7Assays",subset2Flag))==1) {
        heading1=paste("All samples: GI + SCLC, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
        disUniq=cbind(c(""),c("SCLC, pancreas, colon, other GI"))
        disUniq=cbind(c("lungSmall","pancreas","colon","other"),c("SCLC","Pancreas","Colon","Other GI"))
        varList=c("disease")
        varName=c("disease")
        geneCorList=list(lungSmall=NULL,pancreas=NULL,colon=NULL,other=NULL)
    } else if (subset2Flag%in%c("_lungSmallJeffVettedHighGradeGI_T5aT7Assays")) {
        clinThis=clinThis[which(clinThis$disOntTerm3=="SCLC" | (clinThis$grade=="High" & clinThis$disOntTerm3=="pancreas+colon+otherGI")),]
        heading1=paste("Jeff vetted, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
        disUniq=cbind(c(""),c("Pancreas + colon + other GI"))
        disUniq=cbind(c("lungSmall","pancreas","colon","other"),c("SCLC","Pancreas","Colon","Other GI"))
        varList=c("disease")
        varName=c("disease")
        geneCorList=list(lungSmall=NULL,pancreas=NULL,colon=NULL,other=NULL)
    } else if (subset2Flag%in%c("_lungSmall_T5aT7Assays")) {
        clinThis=clinThis[which(clinThis$disOntTerm3=="SCLC" | clinThis$disOntTerm3=="pancreas+colon+otherGI"),]
        heading1=paste("Jeff vetted, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
        disUniq=cbind(c(""),c("Pancreas + colon + other GI"))
        disUniq=cbind(c("lungSmall","pancreas","colon","other"),c("SCLC","Pancreas","Colon","Other GI"))
        varList=c("disease")
        varName=c("disease")
        geneCorList=list(lungSmall=NULL,pancreas=NULL,colon=NULL,other=NULL)
    } else if (subset2Flag=="_lungSmallFionaVettedSmallCellGI_noOtherGI_T5aT7Assays") {
        clinThis=clinThis[which(clinThis$disOntTerm3=="SCLC" | (clinThis$fionaCol1=="fionaY" & clinThis$cellSize=="SC")),]
        heading1=paste("Fiona vetted, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
        disUniq=cbind(c("lungSmall","pancreas","colon"),c("SCLC","Pancreas","Colon"))
        varList=c("disease")
        varName=c("disease")
        geneCorList=list(lungSmall=NULL,pancreas=NULL,colon=NULL)
    } else {
        clinThis=clinThis[which(!is.na(clinThis$fionaCol1)),]
        heading1=paste("Fiona checked, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
        disUniq=cbind(c("","pancreas","colon","other"),c("Pancreas + colon + other GI","Pancreas","Colon","Other GI"))
        varList=c("fionaCol1","cellSize")
        varName=c("vetted","cellSize")
        geneCorList=list(lungSmall=NULL,pancreas=NULL,colon=NULL)
    }
    geneMutList=list(lungSmall=NULL,pancreas=NULL,colon=NULL)
    
    cat("\n\n================ ",subset2Flag," ===============\n\n",sep="")
    print(table(clinThis$tmbLevel,exclude=NULL))
    
    for (dId in 1:nrow(disUniq)) {
        cat("\n\n================ Disease: ",disUniq[dId,1]," ===============\n\n",sep="")
        if (length(grep("_withinSmallCell",subset2Flag))==1) {
            if (disUniq[dId,1]=="") {
                samIdAll=1:nrow(clinThis)
            } else {
                samIdAll=which(clinThis$disOntTerm2%in%c("lungSmall",disUniq[dId,1]))
            }
        } else {
            if (disUniq[dId,1]=="") {
                samIdAll=1:nrow(clinThis)
            } else {
                samIdAll=which(clinThis$disOntTerm2==disUniq[dId,1])
            }
        }
        for (varId in 1:length(varList)) {
            cat("\n\n================ Variable: ",varList[varId],"\n\n",sep="")

            samIdThis=samIdAll
            outAnyAlt=datGPThis[,grep("anyAlt_",colnames(datGPThis))]
            #colnames(outAnyAlt)=sapply(colnames(outAnyAlt),function(x) strsplit(x,"_")[[1]][2],USE.NAMES=F)
            j=match(paste("anyAlt_",clinThis$id[samIdThis],sep=""),colnames(outAnyAlt))
            if (any(is.na(j))) {
                cat("Mismatch in patients !!!\n")
                break
            }
            outAnyAlt=outAnyAlt[,j]
            
            arrayData2=outAnyAlt
            distMat=getCosineDist(arrayData2)
            geneCorMat=1-as.matrix(distMat)
            genePvMat=geneOrMat=geneComutMat=matrix(nrow=nrow(arrayData2),ncol=nrow(arrayData2))
            for (i1 in 1:(nrow(outAnyAlt)-1)) {
                for (i2 in (i1+1):nrow(outAnyAlt)) {
                    x=table(arrayData2[i1,],arrayData2[i2,])
                    geneComutMat[i1,i2]=geneComutMat[i2,i1]=sum(arrayData2[i1,]==1 & arrayData2[i2,]==1,na.rm=T)
                    if (min(dim(x))>1) {
                        res=fisher.test(x)
                        genePvMat[i1,i2]=genePvMat[i2,i1]=res$p.value
                        geneOrMat[i1,i2]=geneOrMat[i2,i1]=res$estimate
                    }
                    if (testFlag=="cosine") {
                        if (sum(outAnyAlt[i1,],na.rm=T)<samSize2 | sum(outAnyAlt[i2,],na.rm=T)<samSize2) geneCorMat[i1,i2]=geneCorMat[i2,i1]=NA
                    }
                }
            }
            if (testFlag=="cosine") {
                geneCorMat[abs(geneCorMat)<geneCorTh]=NA
                geneOrMat[abs(geneCorMat)<geneCorTh]=NA
                geneComutMat[abs(geneCorMat)<geneCorTh]=NA
            } else {
                x=genePvMat[lower.tri(genePvMat)]
                tmp=rep(NA,length(x))
                out=data.frame(id=datGPName,pv=x,qv=tmp,bh=tmp,by=tmp,holm=tmp,bonf=tmp,stringsAsFactors=F)
                i=which(!is.na(out$pv))
                out$qv=NA
                res=try(qvalue(out$pv[i]))
                if ("qvalue"%in%class(res)) out$qv[i]=res$qvalues
                out$by=p.adjust(out$pv,method="BY")
                out$holm=p.adjust(out$pv,method="holm")
                out$bonf=p.adjust(out$pv,method="bonferroni")
                cat("Q-value < 0.05:\n")
                print(out[which(out[,"qv"]<.05),])
                cat("BY < 0.05:\n")
                print(out[which(out[,"by"]<.05),])
                cat("Holm's < 0.05:\n")
                print(out[which(out[,"holm"]<.05),])
                cat("Bonferroni < 0.05:\n")
                print(out[which(out[,"bonf"]<.05),])
                
                geneCorMat[genePvMat>=pThres]=NA
                geneOrMat[genePvMat>=pThres]=NA
                geneComutMat[genePvMat>=pThres]=NA
            }
            if (F) {
                fNameOut=paste("cor_",disUniq[dId],fName1,sep="")
                tbl=geneCorMat; colnames(tbl)=rownames(arrayData2)
                tbl=data.frame(gene=rownames(arrayData2),tbl)
                write.table(tbl, paste(fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
                fNameOut=paste("lor_",disUniq[dId],fName1,sep="")
                tbl=geneOrMat; colnames(tbl)=rownames(arrayData2)
                tbl=data.frame(gene=rownames(arrayData2),tbl)
                write.table(tbl, paste(fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
            }
            fNameOut=paste("pv_",disUniq[dId],fName1,sep="")
            tbl=genePvMat; colnames(tbl)=rownames(arrayData2)
            tbl=data.frame(gene=rownames(arrayData2),tbl)
            write.table(tbl, paste(fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
            fNameOut=paste("comut_",disUniq[dId],fName1,sep="")
            tbl=geneComutMat; colnames(tbl)=rownames(arrayData2)
            tbl=data.frame(gene=rownames(arrayData2),tbl)
            write.table(tbl, paste(fNameOut,".txt",sep=""), sep="\t", col.names=T, row.names=F, quote=F)
            if (testFlag=="cosine") {
                y=geneCorMat[lower.tri(geneCorMat)]
            } else {
                y=geneOrMat[lower.tri(geneOrMat)]
            }
            names(y)=datGPName
            geneCorList[[disUniq[dId,1]]]=y[!is.na(y)]
            nm=unique(c(sapply(names(geneCorList[[disUniq[dId,1]]]),function(x) {strsplit(x,"_")[[1]]},USE.NAMES=F)))
            geneMutList[[disUniq[dId,1]]]=apply(arrayData2[which(rownames(arrayData2)%in%nm),],1,function(x) {mean(x!=0,na.rm=T)})
        }
    }
}

## -------------------

if ("other"%in%names(geneCorList)) {
    geneList=list(lungSmall=NULL,pancreas=NULL,colon=NULL,other=NULL)
} else {
    geneList=list(lungSmall=NULL,pancreas=NULL,colon=NULL)
}
for (k in 1:length(geneList)) {
    cat("\n\n=============== ",names(geneList)[[k]],"\n\n",sep="")
    geneList[[k]]=unique(c(sapply(names(geneCorList[[k]]),function(x) {strsplit(x,"_")[[1]]},USE.NAMES=F)))
    geneMutList[[k]]=geneMutList[[k]][match(geneList[[k]],names(geneMutList[[k]]))]
    cat("No. of co-mutated genes: ",length(geneList[[k]]),"\n",sep="")
}
cat("Total no. of co-mutated genes: ",sum(!duplicated(unlist(geneList))),"\n",sep="")
save(geneMutList,geneList,file=paste("geneList",fName1,".RData",sep=""))


table(p.adjust(genePvMat[lower.tri(genePvMat)],method="BH")<.05)

if (F) {
    for (k in 1:length(geneMutList)) print(sum("MEN1"%in% geneList[[k]]))
    geneId=c("MEN1","DAXX","TSC2")
    sampleId=NULL
    save(geneId,sampleId,file=paste("geneSampleId_pancreas.RData",sep=""))
}

}
