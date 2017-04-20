## Run analysis.R section 1 first
#
# Actionable genes
# grpUniq=c("BRAF","SETD2","MMR","ATRX","BRCA1","BRCA2")
#
# Run genesetFlag="" & genesetFlag="_actionableGene"
# Run subset2Flag="_lungSmallJeffVettedHighGradeGI_T5aT7Assays" first then subset2Flag="_lungSmallFionaVettedGI_T5aT7Assays"
#
# Set as significant if q-value>pThres and samSize>=samSize[1]
# Show % altered if samSize>=samSize[2]

## ---------------------
## ---------------------
## Comparison within a disease

library(qvalue)
library(coin)

nameValueList=c("allSamples","jeffVetted","fionaVetted")
nameValueList=c("allSamples")
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

    source("code/funcs.R")
    out=getFamilyLevelInfo()
    datGP_m=out$datGP_m
    candGene=out$candGene

    geneFamilyFlag=F
    geneFamilyFlag=T
    
    testFlag="_3tmbLevels"
    testFlag="_2tmbLevels"
    testFlag="_tmbMut"

    ## -------------------

    confIntFlag=T
    confIntFlag=F

    samSize=c(0,0)
    samSize=c(50,10)
    samSize=c(50,30)

    datadir2=""
    geneFlag=F
    geneFlag=T


    geneFlag=F
    samSize=c(50,0)
    samSize=c(0,0)

    ## All genesets
    #genesetList=c("_swiSnfComp","_histoneMod","_swiSnfHisMod","_swiSnfPlusHisMod","_atmAtr","_rbEtc","_swiSnfEtc")
    ## For component level
    genesetList=c("_swiSnfComp","_histoneMod","_swiSnfHisMod","_swiSnfPlusHisMod","_atmAtr","_rbEtc")
    ## For gene level
    genesetList=c("_swiSnfComp","_histoneMod","_swiSnfPlusHisMod","_atmAtr","_rbEtc")
    genesetList=c("_swiSnfEtc")
    genesetList=c("_swiSnf")
    genesetList=c("_swiSnfComp")
    genesetList=c("_swiSnf","_swiSnfComp")
    genesetList=c("_swiSnfHisModAtmAtr")
    
    for (genesetFlag in genesetList) {
        geneset2List=c("30topGene")
        switch(genesetFlag,
        "_swiSnf"={geneset2List=c("swiSnf")
        },
        "_histoneMod"={geneset2List=c("histoneModifier")
        },
        "_atmAtr"={geneset2List=c("atmAtr")
        },
        "_swiSnfHisMod"={geneset2List=c("swiSnfHistoneModifier")
        },
        "_swiSnfPlusHisMod"={geneset2List=c("swiSnfPlusHistoneModifier")
        },
        "_swiSnfHisModAtmAtr"={geneset2List=c("swiSnfHistoneModifierAtmAtr")
        },
        "_swiSnfPlusHisModPlusAtmAtr"={geneset2List=c("swiSnfPlusHistoneModifierPlusAtmAtr")
        },
        "_swiSnfComp"={geneset2List=c("swiSnfComponent")
        },
        "_swiSnfCompEtc"={geneset2List=c("swiSnfCompEtc")
        },
        "_swiSnfEtc"={geneset2List=c("swiSnfEtc")
        },
        "_rbEtc"={geneset2List=c("rbEtc")
        }
        )
        
        if (geneFamilyFlag) {
            datGPThis=datGP_m
        } else {
            datGPThis=datGP
        }
        
        for (subset2Flag in subset2List) {
            
            cat("\n\n================",genesetFlag,subset2Flag,"===============\n\n")
            pThres=0.05
            
            altTypeUniq1=sort(unique(clin1$altType))
            
            #altTypeUniq2=cbind(c(altTypeUniq1,paste(altTypeUniq1[1],"+",altTypeUniq1[3]),paste(altTypeUniq1[1:3],"+",altTypeUniq1[4])),as.character(c(1:length(altTypeUniq1),10*(1)+3,10*(1:3)+4)))
            altTypeUniq2=cbind(c(altTypeUniq1,paste(altTypeUniq1[1],"+",altTypeUniq1[2]),paste(altTypeUniq1[1],"+",altTypeUniq1[3]),paste(altTypeUniq1[1:3],"+",altTypeUniq1[4]),paste(altTypeUniq1[1],"+",altTypeUniq1[2],"+",altTypeUniq1[4])),as.character(c(1:length(altTypeUniq1),10*(1)+2:3,10*(1:3)+4,124)))
            altTypeUniq2=cbind(c(altTypeUniq1,"multiple"),as.character(c(1:length(altTypeUniq1),10)))
            altTypeUniq3="anyAlt"
            
            outFormat="png"
            outFormat="pdf"
            
            if ("multiple"%in%altTypeUniq2[,1]) {
                colList=c("red","blue","orange","cyan","yellow","indianred","yellow2","skyblue","bisque3","indianred4")
            } else {
                colList=c("red","blue","purple","cyan","skyblue","orange","yellow","darkgreen")
                colList=c("red","blue","indianred","cyan","orange","yellow","skyblue","bisque3")
                colList=c("red","blue","indianred","cyan","yellow2","orange","yellow","skyblue","bisque3","indianred4")
            }
            
            propFlag=">="
            propFlag=">"
            
            candGeneThis=getCandidateGenes(genesetFlag)
            
            gene1=unique(clin1$gene[clin1$assayVersion=="T5a"])
            gene2=unique(clin1$gene[clin1$assayVersion=="T7"])
            if (subset2Flag%in%c("_allSamples_T5aT7Assays","_lungSmallJeffVettedHighGradeGI_T5aT7Assays")) {
                j=1:nrow(clin3)
            } else {
                j=which(!is.na(clin32$fionaCol1))
            }
            j=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7") & ((clin1$id%in%clin32$id[j]) | (!clin1$disOntTerm2%in%clin3$disOntTerm2[j])))
            clinThis=cbind(clin3,clin32[,which(!names(clin32)%in%names(clin3))])[which(clin3$id%in%clin1$id[j]),]
            clinThis$fionaCol1[which(!clinThis$fionaCol1%in%c("N","Y","Y/N"))]=NA
            clinThis$fionaCol1=paste("fiona",sub("/","", clinThis$fionaCol1),sep="")
            clinThis$fionaCol1[clinThis$fionaCol1=="fionaNA"]=NA
            if (F) {
                clinThis$cellSize[which(clinThis$cellSize%in%c("MC","SC/LC","SC/MC"))]="SMLC"
                clinThis$cellSize2=clinThis$cellSize
                clinThis$cellSize2[which(clinThis$disOntTerm2=="lungSmall")]="SC"
            }
            clinThis$cellSize2=clinThis$cellSize
            clinThis$cellSize2[which(clinThis$disOntTerm2=="lungSmall")]="SC"
            clinThis$cellSize2[which(clinThis$cellSize%in%c("MC","SC/LC","SC/MC"))]="other"
            clinThis$cellSize2[which(!clinThis$cellSize%in%c("SC","LC","other"))]=NA
            clinThis$disease=clinThis$disOntTerm2
            clinThis$disease[which(clinThis$disease=="lungSmall")]="SCLC"
            clinThis$disOntTerm3[which(clinThis$disOntTerm3=="lungSmall")]="SCLC"
            clinThis$disOntTerm3[which(clinThis$disOntTerm3=="combined")]="pancreas+colon+otherGI"
            
            if (length(grep("_withinSmallCell",subset2Flag))==1) {
                clinThis=clinThis[which(clinThis$cellSize2=="SC"),]
                heading1=paste("Small cell Fiona checked + SCLC, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
                #disUniq=cbind(c("pancreas","colon","other"),c("Pancreas","Colon","Other GI"))
                disUniq=cbind(c(""),c("SCLC, pancreas, colon, other GI"))
                varList=c("disease")
                varName=c("disease")
            } else if (length(grep("_allSamples_T5aT7Assays",subset2Flag))==1) {
                heading1=paste("All samples: GI + SCLC, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
                disUniq=cbind(c(""),c("SCLC, pancreas, colon, other GI"))
                varList=c("disease")
                varName=c("disease")
            } else if (length(grep("_lungSmallJeffVettedHighGradeGI_T5aT7Assays",subset2Flag))==1) {
                clinThis=clinThis[which(clinThis$disOntTerm3=="SCLC" | (clinThis$grade=="High" & clinThis$disOntTerm3=="pancreas+colon+otherGI")),]
                heading1=paste("Jeff vetted high grade GI + SCLC, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
                disUniq=cbind(c(""),c("SCLC, pancreas, colon, other GI"))
                varList=c("disease")
                varName=c("disease")
            } else if (length(grep("_lungSmallFionaVettedGI",subset2Flag))==1) {
                clinThis=clinThis[which(clinThis$disOntTerm3=="SCLC" | clinThis$fionaCol1=="fionaY"),]
                heading1=paste("Fiona vetted GI + SCLC, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
                #disUniq=cbind(c("lungSmall","pancreas","colon"),c("SCLC","pancreas","colon"))
                disUniq=cbind(rep("",3),rep("SCLC, pancreas, colon, other GI",3),c("SCLC","pancreas","colon"))
                disUniq=cbind(c(""),c("SCLC, pancreas, colon, other GI"))
                varList=c("cellSize2"); varName=c("cellSize")
                varList=varName=c("disease")
            } else if (length(grep("_lungSmallFionaVettedSmallCellGI",subset2Flag))==1) {
                clinThis=clinThis[which(clinThis$disOntTerm3=="SCLC" | (clinThis$fionaCol1=="fionaY" & clinThis$cellSize2=="SC")),]
                heading1=paste("Small cell Fiona vetted + SCLC, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
                disUniq=cbind(c(""),c("SCLC, pancreas, colon, other GI"))
                varList=c("disOntTerm3","disease")
                varName=c("disOntTerm3","disease")
            } else if (length(grep("_lungSmallSmallCellGI",subset2Flag))==1) {
                clinThis=clinThis[which(clinThis$disOntTerm3=="SCLC" | clinThis$cellSize2=="SC"),]
                heading1=paste("Small cell Fiona checked + SCLC, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
                disUniq=cbind(c(""),c("SCLC, small cell pancreas, colon, other GI pooled"))
                varList=c("disOntTerm3")
                varName=c("disease")
            } else if (length(grep("_lungSmallLargeCellGI",subset2Flag))==1) {
                clinThis=clinThis[which(clinThis$disOntTerm3=="SCLC" | clinThis$cellSize2=="LC"),]
                heading1=paste("Large cell Fiona checked + SCLC, ",paste(sort(unique(clinThis$assayVersion)),collapse=" & ")," assay samples",sep="")
                disUniq=cbind(c(""),c("SCLC, large cell pancreas, colon, other GI pooled"))
                varList=c("disOntTerm3")
                varName=c("disease")
            } else if (length(grep("_tmbLevel",subset2Flag))==1) {
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
            if (length(grep("_withinVetted|_fionaVetted",subset2Flag))==1) {
                clinThis=clinThis[which(clinThis$fionaCol1=="fionaY"),]
                heading1=sub("Fiona checked","Fiona vetted",heading1)
                k=which(!varList%in%"fionaCol1")
                varList=varList[k]
                varName=varName[k]
            }
            if (length(grep("_noOtherGI",subset2Flag))==1) {
                clinThis=clinThis[which(clinThis$disease!="other"),]
            }
            if (testFlag!="") {
                clinThis=clinThis[which(!is.na(clinThis$tmbMut)),]
                heading1=paste(heading1,"\nConsidered only samples with non-missing tumor burgen information",sep="")
            }
            if (ncol(disUniq)==2) disUniq=cbind(disUniq,rep("",nrow(disUniq)))
            
            if (genesetFlag=="_swiSnfHisModAtmAtr") {
                j=match(paste("anyAlt_",clinThis$id,sep=""),colnames(datGPThis))
                if (any(is.na(j))) {
                    cat("Mismatch in patients !!!\n")
                    break
                }
                tbl=t(datGPThis[which(rownames(datGPThis)%in%candGeneThis$family2),j])
                clinThis=cbind(clinThis,tbl)
            }
            
            
            for (dId in 1:nrow(disUniq)) {
                #cat("\n\n================",disUniq[dId,1],"===============\n\n")
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
                    #cat("\n\n================",varList[varId],"\n\n")
                    if (varList[varId]=="fionaCol1") {
                        if (disUniq[dId,1]=="") next
                        grpUniq3=paste("fiona",sort(c("N","Y","YN")),sep="")
                        grpUniq3=cbind(grpUniq3,sub("fiona","vetted ",sub("YN","Borderline",grpUniq3)),sub("fiona","",sub("YN","B",grpUniq3)))
                        grpUniq3[,2]=paste("Vetted ",c("Yes","No","Borderline"),sep="")
                        #grpUniq3=cbind(grpUniq3,paste("fiona",sub("/","", grpUniq3),sep=""))
                    } else if (varList[varId]=="cellSize") {
                        grpUniq3=c("SC","LC","SMLC")
                        grpUniq3=cbind(grpUniq3,c("Small Cell","Large Cell","Small/Medium/Large Cell"),c("S","L","N"))
                        #} else if (length(grep("_lungSmallSmallCellGI",subset2Flag))==1 | length(grep("_lungSmallLargeCellGI",subset2Flag))==1) {
                    } else if (varList[varId]=="cellSize2") {
                        grpUniq3=c("SC","LC","other")
                        grpUniq3=cbind(grpUniq3,c("Small Cell","Large Cell","Small/Medium/Large Cell"),c("S","L","N"))
                        #} else if (length(grep("_lungSmallSmallCellGI",subset2Flag))==1 | length(grep("_lungSmallLargeCellGI",subset2Flag))==1) {
                    } else if (varList[varId]=="disOntTerm3") {
                        grpUniq3=c("SCLC","pancreas+colon+otherGI")
                        grpUniq3=cbind(grpUniq3,capWords(grpUniq3),toupper(substr(grpUniq3,1,1)))
                    } else if (varList[varId]=="tmbLevel") {
                        grpUniq3=c("Low","Intermediate")
                        grpUniq3=cbind(grpUniq3,paste("TMB level ",capWords(grpUniq3),sep=""),toupper(substr(grpUniq3,1,1)))
                    } else {
                        grpUniq3=c("SCLC","pancreas","colon","other")
                        grpUniq3=cbind(grpUniq3,capWords(sub("other","other GI",grpUniq3)),toupper(substr(grpUniq3,1,1)))
                    }
                    x=table(clinThis[samIdAll,varList[varId]])
                    k1=match(grpUniq3[,1],names(x))
                    #k=match(grpUniq3[,1],names(x)); k1=k[!is.na(k)]
                    k=which(x[k1]>=samSize[2])
                    if (length(k)<2) next
                    grpUniq3=grpUniq3[k,]
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
                    samIdThis=samIdAll
                    outAnyAlt=datGPThis[,grep("anyAlt_",colnames(datGPThis))]
                    rownames(outAnyAlt)=rownames(datGPThis)
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
                    
                    #out4Alt=datGPThis[,grep("anyAlt_",colnames(datGPThis))]
                    j=grep("anyAlt_",colnames(datGPThis))
                    #out4Alt=matrix(0,nrow=nrow(datGPThis),ncol=ncol(datGPThis[,j]),dimnames=list(rownames(datGPThis),colnames(datGPThis)[j]))
                    out4Alt=matrix(0,nrow=nrow(datGPThis),ncol=ncol(datGPThis[,j]),dimnames=list(rownames(datGPThis),colnames(datGPThis)[j]))
                    id1=sapply(colnames(datGPThis),function(x) {strsplit(x,"_")[[1]][2]},USE.NAMES=F)
                    id2=sapply(colnames(out4Alt),function(x) {strsplit(x,"_")[[1]][2]},USE.NAMES=F)
                    jj=match(sapply(colnames(datGPThis),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F),altTypeUniq1)
                    for (g1 in 1:nrow(out4Alt)) {
                        for (a1 in 1:length(altTypeUniq1)) {
                            j=which(jj==a1 & datGPThis[g1,]==1)
                            j2=match(id1[j],id2)
                            j21=j2[which(out4Alt[g1,j2]==0)]
                            j22=j2[which(out4Alt[g1,j2]!=0)]
                            out4Alt[g1,j21]=a1
                            if ("multiple"%in%altTypeUniq2[,1]) {
                                out4Alt[g1,j22]=10
                            } else {
                                out4Alt[g1,j22]=10*out4Alt[g1,j22]+a1
                            }
                        }
                    }
                    j=match(paste("anyAlt_",clinThis$id[samIdThis],sep=""),colnames(out4Alt))
                    out4Alt=out4Alt[,j]
                    
                    out4AltPropCombAll=matrix(0,nrow=nrow(outAnyAlt),ncol=nrow(altTypeUniq2)+1,dimnames=list(rownames(outAnyAlt),c("noAlt",altTypeUniq2[,1])))
                    out4AltPropAll=matrix(0,nrow=nrow(outAnyAlt),ncol=nrow(grpUniq3)*(nrow(altTypeUniq2)+1),dimnames=list(rownames(outAnyAlt),paste(rep(c("noAlt",altTypeUniq2[,1]),each=nrow(grpUniq3)),grpUniq3[,1],sep="_")))
                    colId=c("0",altTypeUniq2[,2])
                    x=as.character(sort(unique(c(out4Alt))))
                    k=which(!x%in%colId)
                    if (length(k)!=0) cat("out4AltPropCombAll: ",x[k]," not in list !!!",sep="")
                    for (g1 in 1:nrow(outPropAll)) {
                        res=table(out4Alt[g1,])
                        out4AltPropCombAll[g1,match(names(res),colId)]=res/sum(res)
                        for (k in 1:nrow(grpUniq3)) {
                            j=which(clinThis[samIdThis,varList[varId]]==grpUniq3[k,1])
                            kk=grep(paste("_",grpUniq3[k,1],sep=""),colnames(out4AltPropAll),fixed=T)
                            res=table(out4Alt[g1,j])
                            out4AltPropAll[g1,kk[match(names(res),colId)]]=res/sum(res)
                        }
                    }
                    out4AltPropCombAll=out4AltPropCombAll[,altTypeUniq2[,1]]
                    out4AltPropAll=out4AltPropAll[,paste(rep(altTypeUniq2[,1],each=nrow(grpUniq3)),grpUniq3[,1],sep="_")]
                    
                    outPvAll=matrix(nrow=nrow(outPropAll),ncol=nrow(grpUniq3),dimnames=list(rownames(outPropAll),paste("pv_anyAlt_",grpUniq3[,1],sep="")))
                    outEstAll=matrix(nrow=nrow(outPropAll),ncol=nrow(grpUniq3),dimnames=list(rownames(outPropAll),paste("qv_anyAlt_",grpUniq3[,1],sep="")))
                    #print(nrow(outPvAll))
                    timeStamp=Sys.time()
                    #print(format(timeStamp, "%x %X"))
                    grpThis=clinThis[samIdThis,"tmbLevel"]
                    grpThis[which(grpThis%in%c("High","Intermediate"))]="Intermediate/High"
                    for (k in 1:nrow(grpUniq3)) {
                        j=which(clinThis[samIdThis,varList[varId]]%in%grpUniq3[k,1])
                        print(nrow(clinThis[j,]))
                        print(table(swiSnf=clinThis$SWI.SNF[j],histMod=clinThis$Modifier.Histone[j],atmAtr=clinThis$ATM.ATR[j]))

                        library(MASS)
                        cat("======== loglm\n")
                        res=loglm(~SWI.SNF+Modifier.Histone+ATM.ATR, data=clinThis[j,])
                        x=table(SWI.SNF=clinThis$SWI.SNF[j],Modifier.Histone=clinThis$Modifier.Histone[j],ATM.ATR=clinThis$ATM.ATR[j])
                        res=loglm(~SWI.SNF+Modifier.Histone+ATM.ATR, data=x)
                        print(res)
                        res=loglm(~SWI.SNF*Modifier.Histone*ATM.ATR, data=x)
                        print(res)

                        j=which(clinThis[samIdThis,varList[varId]]%in%grpUniq3[k,1] & !is.na(clinThis[samIdThis,"tmbLevel"]))
                        if (length(j)!=0) {
                            colId=paste("qv_anyAlt_",grpUniq3[k,1],sep="")
                            outEstAll[,colId]=""
                            for (a1 in which(altTypeUniq3=="anyAlt")) {
                                for (g1 in 1:nrow(outPvAll)) {
                                    if (testFlag=="_tmbMut") {
                                        x=table(outAnyAlt[g1,j]); x=matrix(c(x,x),nrow=length(x))
                                        if (nrow(x)>1 && ncol(x)>1) {
                                            res=wilcox_test(clinThis[samIdThis[j],"tmbMut"]~as.factor(outAnyAlt[g1,j]),distribution="exact")
                                            colId=paste("qv_",altTypeUniq3[a1],"_",grpUniq3[k,1],sep="")
                                            outEstAll[g1,colId]=ifelse(statistic(res)>0,"+","-")
                                            colId=paste("pv_",altTypeUniq3[a1],"_",grpUniq3[k,1],sep="")
                                            outPvAll[g1,colId]=pvalue(res)
                                        }
                                    } else {
                                        if (testFlag=="_2tmbLevels") {
                                            x=table(grpThis[j],outAnyAlt[g1,j])
                                        } else {
                                            x=table(clinThis[samIdThis[j],"tmbLevel"],outAnyAlt[g1,j])
                                        }
                                        if (nrow(x)>1 && ncol(x)>1) {
                                            res=fisher.test(x)
                                            #colId=paste("qv_",altTypeUniq3[a1],"_",grpUniq3[k,1],sep="")
                                            #outEstAll[g1,colId]=paste(round(c(res$estimate,res$conf.int),1),collapse=",")
                                            colId=paste("pv_",altTypeUniq3[a1],"_",grpUniq3[k,1],sep="")
                                            outPvAll[g1,colId]=res$p.value
                                        }
                                    }
                                }
                            }
                        }
                    }
                    outPvAll[which(outPvAll>1)]=1
                    #cat("No. of comparisons with pv < 0.05",sum(c(outPvAll)<0.05,na.rm=T),"\n")
                    outQvAll=matrix(nrow=nrow(outPvAll),ncol=ncol(outPvAll),dimnames=list(rownames(outPvAll),sub("pv_","qv_",colnames(outPvAll))))
                    for (k in 1:ncol(outPvAll)) {
                        i=which(!is.na(outPvAll[,k]))
                        if (length(i)!=0) {
                            outQvAll[i,k]=qvalue(outPvAll[i,k])$qvalues
                        }
                    }
                    #cat("No. of comparisons with qv < 0.05",sum(c(outQvAll)<0.05,na.rm=T),"\n")
                    
                    #for (geneset2Flag in c("","10%patient","15%patient")) {
                    for (geneset2Flag in geneset2List) {
                        heading=paste(heading1,"\n")
                        if (length(grep("%patient",geneset2Flag))==1) {
                            propThres=0.10
                            if (length(grep("%patientGrp1Grp2",geneset2Flag))==1) {
                                propThres=as.integer(sub("%patientGrp1Grp2","",geneset2Flag))/100
                                i=match(c("TP53","RB1","APC","CDKN2A","KRAS","MEN1","CDKN2B","CCNE1","DAXX","FBXW7"),rownames(outPropAll))
                                fName=paste("proportionTable_geneBy",capWords(varList[varId]),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),"_",propThres*100,"percPatientGrp1Grp2",fName1,subset2Flag,".txt",sep="")
                                heading=paste(heading,"Table of ",length(i)," genes with alterations in ",propFlag," ",propThres*100,"% of patients for group 1 or 2",sep="")
                            } else if (length(grep("%patientDisease",geneset2Flag))==1) {
                                propThres=as.integer(sub("%patientDisease","",geneset2Flag))/100
                                i=match(c("TP53","RB1","MLL2","APC","CDKN2A","LRP1B","KRAS","MEN1","FBXW7","DAXX"),rownames(outPropAll))
                                fName=paste("proportionTable_geneBy",capWords(varList[varId]),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),"_",propThres*100,"percPatientInAnyDisease",fName1,subset2Flag,".txt",sep="")
                                heading=paste(heading,"Table of ",length(i)," genes with alterations in ",propFlag," ",propThres*100,"% of patients for SCLC, pancreas or colon",sep="")
                            }  else if (subset2Flag=="_lungSmallFionaVettedGI_T5aT7Assays" & geneFlag) {
                                propThres=as.integer(sub("%patient","",geneset2Flag))/100
                                fName2=paste("geneSampleId_geneBy",capWords(varList[varId]),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),"_",propThres*100,"percPatient",fName1,"_lungSmallJeffVettedHighGradeGI_T5aT7Assays.RData",sep="")
                                load(paste(datadir2,fName2,sep=""))
                                i=match(geneId,rownames(outPropAll))
                                fName=paste("proportionTable_geneBy",capWords(varList[varId]),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),"_",propThres*100,"percPatient",fName1,subset2Flag,".txt",sep="")
                                heading=paste(heading,"Table of ",length(i)," genes with alterations in ",propFlag," ",propThres*100,"% of patients in any group within ",disUniq[dId,2],sep="")
                            } else {
                                propThres=as.integer(sub("%patient","",geneset2Flag))/100
                                i=c()
                                if (disUniq[dId,3]=="") {
                                    k1=1:ncol(outPropAll)
                                } else {
                                    k1=which(colnames(outPropAll)==disUniq[dId,3])
                                }
                                for (k in k1) {
                                    if (propFlag==">") {
                                        i=c(i,which(round(outPropAll[,k],2)>propThres))
                                    } else {
                                        i=c(i,which(round(outPropAll[,k],2)>=propThres))
                                    }
                                }
                                i=sort(unique(i))
                                fName=paste("proportionTable_geneBy",capWords(varList[varId]),ifelse(disUniq[dId,3]=="","",paste("_within",capWords(disUniq[dId,3]),sep="")),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),"_",propThres*100,"percPatient",fName1,subset2Flag,".txt",sep="")
                                heading=paste(heading,"Table of ",length(i)," genes with alterations in ",propFlag," ",propThres*100,"% of patients in ",ifelse(disUniq[dId,3]=="",paste("any group within ",disUniq[dId,2],sep=""),disUniq[dId,3]),sep="")
                            }
                        } else if (length(grep("patient",geneset2Flag))==1) {
                            numThres=as.integer(strsplit(geneset2Flag,"patient")[[1]][1])
                            j=1:nrow(clinThis)
                            x2=strsplit(geneset2Flag,"patientIn")[[1]]
                            if (length(x2)==2) {
                                j=which(tolower(clinThis$disOntTerm2)==tolower(x2[2]))
                                nm=clinThis$disease[j][1]
                            }
                            x=apply(out4Alt[,j],1,function(x) {sum(x!=0,na.rm=T)})
                            i=order(x,decreasing=T)
                            i=i[which(x[i]>=numThres)]
                            if (length(i)==0) {cat("No genes !!!\n\n"); next}
                            fName=paste("proportionTable_geneBy",capWords(varList[varId]),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),"_",geneset2Flag,fName1,subset2Flag,".txt",sep="")
                            heading=paste(heading,"Table of ",length(i)," genes with >= ",numThres," patients with alterations",ifelse(length(x2)==2,paste(" in ",nm,sep=""),""),sep="")
                        } else if (length(grep("topGene",geneset2Flag))==1) {
                            #numThres=as.integer(sub("topGene","",geneset2Flag))
                            numThres=as.integer(strsplit(geneset2Flag,"topGene")[[1]][1])
                            j=1:nrow(clinThis)
                            x2=strsplit(geneset2Flag,"topGeneIn")[[1]]
                            if (length(x2)==2) {
                                j=which(tolower(clinThis$disOntTerm2)==tolower(x2[2]))
                                nm=clinThis$disease[j][1]
                            }
                            x=apply(out4Alt[,j],1,function(x) {sum(x!=0,na.rm=T)})
                            i=order(x,decreasing=T)[1:numThres]
                            fName=paste("proportionTable_geneBy",capWords(varList[varId]),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),"_",geneset2Flag,fName1,subset2Flag,".txt",sep="")
                            heading=paste(heading,"Table of ",length(i)," genes with most alterations",ifelse(length(x2)==2,paste(" in ",nm,sep=""),""),sep="")
                        } else if (geneset2Flag%in%c("actionableGene","swiSnf","swiSnfComponent","histoneModifier","atmAtr","swiSnfHistoneModifier","swiSnfPlusHistoneModifier","swiSnfHistoneModifierAtmAtr","swiSnfPlusHistoneModifierPlusAtmAtr","rbEtc","swiSnfCompEtc","swiSnfEtc")) {
                            print("------------------ 20 ---------------")
                            if (geneFamilyFlag) {
                                i=match(toupper(unique(candGeneThis$family2)),toupper(rownames(outPropAll)))
                            } else {
                                i=match(toupper(candGeneThis$gene),toupper(rownames(outPropAll)))
                            }
                            fName=paste("proportionTable_geneBy",capWords(varList[varId]),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),genesetFlag,fName1,subset2Flag,".txt",sep="")
                            heading=paste(heading,"Table of ",length(i)," genes with any alteration in a group within ",disUniq[dId,2],sep="")
                        } else {
                            i=1:nrow(outPropAll)
                            fName=paste("proportionTable_geneBy",capWords(varList[varId]),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),".txt",sep="")
                            heading=paste(heading,"Table of ",length(i)," genes with any alteration in a group within ",disUniq[dId,2],sep="")
                        }
                        fName=sub("proportionTable_","tmbAssociation_",fName)
                        if (length(i)!=0) {
                            nmR=rownames(outPropAll)[i]
                            outPropThis=round(outPropAll[i,colnames(outPropAll)],2)
                            if (length(i)>1) {
                                outPercThis=apply(outPropThis,c(1,2),function(x) {
                                    y=paste(" ",x*100,"% ",sep="")
                                    y
                                })
                            } else if (length(i)==1) {
                                outPropThis=matrix(outPropThis,nrow=1,dimnames=list(nmR,colnames(outPropAll)))
                                outPercThis=matrix(paste(" ",c(outPropThis)*100,"% ",sep=""),nrow=1,dimnames=list(nmR,colnames(outPropAll)))
                            }
                            i=match(rownames(outQvAll),rownames(outPropThis)); i1=which(!is.na(i)); i2=i[i1]
                            
                            cat("No. of comparisons with pv < 0.05",sum(c(outPvAll[i1,])<0.05,na.rm=T),"\n")
                            if (length(i1)==1) {
                                outQvAll=outPvAll
                                colnames(outQvAll)=sub("pv_","qv_",colnames(outQvAll))
                            } else {
                                outQvAll=matrix(nrow=nrow(outPvAll),ncol=ncol(outPvAll),dimnames=list(rownames(outPvAll),sub("pv_","qv_",colnames(outPvAll))))
                                x=outQvAll
                                for (k in 1:ncol(outPvAll)) {
                                    i=i1
                                    if (length(i)!=0) {
                                        x[i,k]=p.adjust(outPvAll[i,k],method="BH")
                                        if (F) {
                                        res=try(qvalue(outPvAll[i,k]))
                                        if (class(res)=="try-error") {
                                            cat("Could not compute q-values!!!\n")
                                        } else {
                                            outQvAll[i,k]=qvalue(outPvAll[i,k])$qvalues
                                        }
                                        }
                                        outQvAll=x
                                    }
                                }
                                print(table(bh=c(x[i1,]<0.05),qv=c(outQvAll[i1,])<0.05,exclude=NULL))
                            }
                            cat("No. of comparisons with qv < 0.05",sum(c(outQvAll[i1,])<0.05,na.rm=T),"\n")

                            colnames(outPercThis)=paste("anyAlt_",colnames(outPercThis),sep="")
                            res=table(clinThis[samIdThis,varList[varId]][which(clinThis[samIdThis,varList[varId]]%in%grpUniq3[,1])])
                            for (grpId in 1:nrow(grpUniq3)) {
                                if (res[which(names(res)==grpUniq3[grpId,1])]<samSize[1]) next
                                k1=grep("qv_anyAlt_",colnames(outQvAll),fixed=T)
                                k1=k1[grep(grpUniq3[grpId,1],colnames(outQvAll)[k1],fixed=T)]
                                nm1=colnames(outQvAll)[k1]
                                nm11=nm1
                                nm1=sub("qv_anyAlt_","",nm1,fixed=T)
                                nm21=colnames(outPercThis)[grep("anyAlt_",colnames(outPercThis),fixed=T)]
                                nm2=sub("anyAlt_","",colnames(outPercThis)[grep("anyAlt_",colnames(outPercThis),fixed=T)],fixed=T)
                                #x=" "
                                for (k in 1:length(k1)) {
                                    if (res[which(names(res)==nm1[k])]<samSize[1]) next
                                    k21=which(nm2==nm1[k])
                                    i=which(outQvAll[i1,nm11[k]]<pThres)
                                    if (length(i)!=0) {
                                        x=rep("",length(i))
                                        if (confIntFlag) {
                                            x=paste(x,grpUniq3[grpId,3],"(",outEstAll[i1[i],nm11[k]],"),",sep="")
                                        } else {
                                            #x=paste(x,grpUniq3[grpId,3],sep="")
                                            #x=paste(x,"**",sep="")
                                            x=paste(x,outEstAll[i1[i],nm11[k]],sep="")
                                        }
                                        outPercThis[i2[i],k21]=paste(outPercThis[i2[i],k21],x,sep="")
                                    }
                                    i=which(is.na(outEstAll[i1,nm11[k]]))
                                    if (length(i)!=0) {
                                        outPercThis[i2[i],k21]=""
                                    }
                                }
                            }
                            outPercThis=sub(",$","",outPercThis)
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
                            if (length(grep("_noOtherGI",subset2Flag))==1) {
                                heading=sub(", other GI","",heading)
                                nm=sub("+otherGI","",nm,fixed=T)
                            }
                            write.table(heading,file=fName,append=F,col.names=F,row.names=F,sep="\t",quote=F)
                            write.table(paste("N = ",sum(x),sep=""),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                            write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                            write.table(nm,file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                            write.table(tbl,file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                            
                            ## -------------------
                            ## Box plot
                            
                            if (F) {
                                x=apply(out4AltPropCombAll,1,sum)
                                gId=match(rownames(outPropThis),rownames(out4AltPropCombAll))
                                fName2=sub("tmbAssociation_","boxplot_tmb_",sub(".txt","",fName,fixed=T))
                                heading2=sub("Table of ","Boxplots of ",heading)
                                heading2=strsplit(heading2," \n")[[1]][1]
                                #par(mar=c(5, 4, 4, 2) + 0.1)
                                switch(outFormat,
                                "png"={
                                    png(paste(fName2,".png",sep=""),width=2.75*240, height=1.25*240)
                                    cexAxis=0.7; cexAxis2=0.8; cexLab=2
                                },
                                "pdf"={
                                    cexAxis=1.4; cexAxis2=1.2; cexLab=2
                                    if (nrow(grpUniq3)>3) {
                                        cexAxis=1
                                    }
                                    if (length(gId)>13) {
                                        pdf(paste(fName2,".pdf",sep=""),width=3*7,height=7)
                                        cexAxis2=.7
                                    } else {
                                        pdf(paste(fName2,".pdf",sep=""),width=2*7,height=7)
                                        cexAxis2=1.2
                                    }
                                }
                                )
                                
                                cexAxis2=1
                                par(mar=c(8, 5, 0, 0) + 0.1)
                                par(mar=c(8, 5, 1, 0) + 0.1)
                                out=NULL
                                nm=c()
                                for (kk in gId) {
                                    for (grpId in 1:nrow(grpUniq3)) {
                                        ll=grep(paste("_",grpUniq3[grpId,1],sep=""),colnames(outAnyAlt),fixed=T)
                                        out=rbind(out,outAnyAlt[kk,ll])
                                        nm=c(nm,paste(grpUniq3[grpId,2],"    ",sep=""))
                                    }
                                    out=rbind(out,rep(0,ncol(out)))
                                    nm=c(nm,"")
                                }
                                if (length(grep("_noOtherGI",subset2Flag))==1) {
                                    nm=sub("+otherGI","",nm,fixed=T)
                                }
                                nm=sub("+","    \n+ ",nm,fixed=T)
                                rownames(out)=nm
                                grpId=1
                                boxplot(100*t(out),names.arg=rownames(out),main=ifelse(grpId==1,heading2,""),ylab="% alteration",col=colList,las=3,cex.names=cexAxis,cex.axis=cexAxis,cex.lab=cexLab,space=0)
                                x=rownames(outAnyAlt)[gId]
                                axis(side=1,at=seq(1,nrow(out),by=nrow(grpUniq3)+1),labels=x,tick=F,cex.axis=cexAxis2,cex.lab=cexLab)
                                dev.off()
                                
                                i=match(nmR,rownames(outAnyAlt))
                                j=which(clinThis[,varList[varId]]%in%grpUniq3[,1])
                                if (varList[varId]%in%c("disease")) {
                                    j=j[which(!clinThis[j,varList[varId]]%in%"SCLC")]
                                }
                                x=outAnyAlt[i,j]
                                #x=paste(c(x))
                                y=clinThis$tmbMut[j]
                                #y[which(x>200)]=200
                                y=log(y+1)
                                ylim=range(x,na.rm=T)
                                i=1
                                boxplot(y~x[i,],xlim=1:length(nmR),ylim=ylim)
                                for (i in 2:length(nmR)) {
                                    
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
