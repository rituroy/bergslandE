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


# Actionable genes
# grpUniq=c("BRAF","SETD2","MMR","ATRX","BRCA1","BRCA2")
#
# Run genesetFlag="" & genesetFlag="_actionableGene"
# Run subset2Flag="_lungSmallJeffVettedHighGradeGI_T5aT7Assays" first then subset2Flag="_lungSmallFionaVettedGI_T5aT7Assays"
#
# Set as significant if q-value>pThres and samSize>=samSize[1]
# Show % altered if samSize>=samSize[2]
#
# Set geneFamilyFlag, genesetList, subset2List, datadirG

## ---------------------
## ---------------------
## Comparison within a disease

library(qvalue)

nameValueList=c("allSamples","jeffVetted","fionaVetted")
nameValueList=c("fionaVettedUcsf500")
nameValueList=c("allSamples")
nameValueList=c("fionaVetted")
nameValueList=c("ucsf500")
nameValueList=c("ucsf500Fmi")
for (nameValueFlag in nameValueList) {
    switch(nameValueFlag,
    "allSamples"={nameValue=data.frame(name="subset",value="allSamples")
    },
    "jeffVetted"={nameValue=data.frame(name="subset",value="jeffVetted")
    },
    "fionaVetted"={nameValue=data.frame(name="subset",value="fionaVetted")
    },
    "ucsf500"={nameValue=data.frame(name="subset",value="ucsf500")
    },
    "ucsf500Fmi"={nameValue=data.frame(name="subset",value="ucsf500Fmi")
    },
    "fionaVettedUcsf500"={nameValue=data.frame(name="subset",value="fionaVettedUcsf500")
    }
    )
    write.table(nameValue,"config.tmp", sep="=", col.names=F, row.names=F, quote=F)
    datadir="results/"
    if (nameValue$value[which(nameValue$name=="subset")]=="fionaVetted") {
        #load(paste(datadir,"tmp_fiona166_T5aT7Assays.RData",sep=""))
        load(paste(datadir,"tmp_allAssays.RData",sep=""))
        nameValue=read.table("config.tmp",sep="=",h=F,quote="",comment.char="",as.is=T,fill=T,col.names=c("name","value"))
        subset2List=c("_lungSmallFionaVettedGI_T5aT7Assays")
    } else if (nameValue$value[which(nameValue$name=="subset")]%in%c("ucsf500","ucsf500Fmi")) {
        load(paste(datadir,"tmp_allAssays.RData",sep=""))
        datadir="results/"
        load(paste(datadir,"tmp_ucsf500.RData",sep=""))
        nameValue=read.table("config.tmp",sep="=",h=F,quote="",comment.char="",as.is=T,fill=T,col.names=c("name","value"))
        subset2List=paste("_",c("primarySiteEB"),nameValue$value[which(nameValue$name=="subset")],sep="")
        # For genesetList = _30topGeneInUcsf500, _30topGeneInUcsf500Fmi
        subset2List=paste("_",c("poorDiff_","wellDiff_"),nameValue$value[which(nameValue$name=="subset")],sep="")
        # For genesetList = _30topGene
        subset2List=paste("_",c(""),nameValue$value[which(nameValue$name=="subset")],sep="")
    } else if (nameValue$value[which(nameValue$name=="subset")]=="fionaVettedUcsf500") {
        load(paste(datadir,"tmp_allAssays.RData",sep=""))
        datadir="results/"
        load(paste(datadir,"tmp_ucsf500.RData",sep=""))
        names(clinU)[match(c("primarySite","cellSizeEB"),names(clinU))]=c("disOntTerm2","cellSize")
        clinU$cellSizeEB[which(clinU$cellSizeEB=="small")]="SC"
        clinU$cellSize[which(clinU$cellSizeEB=="large")]="LC"
        nameValue=read.table("config.tmp",sep="=",h=F,quote="",comment.char="",as.is=T,fill=T,col.names=c("name","value"))
        subset2List=c("_fionaVettedUcsf500")
    } else {
        load(paste(datadir,"tmp_allAssays.RData",sep=""))
        nameValue=read.table("config.tmp",sep="=",h=F,quote="",comment.char="",as.is=T,fill=T,col.names=c("name","value"))
        subset2List=c("_allSamples_T5aT7Assays")
        if (nameValue$value[which(nameValue$name=="subset")]=="jeffVetted") {
            subset2List=c("_lungSmallJeffVettedHighGradeGI_T5aT7Assays")
        }
    }
    
    source(paste(dirSrc3,"funcs.R",sep=""))
    if (nameValue$value[which(nameValue$name=="subset")]%in%c("ucsf500","ucsf500Fmi")) {
        samSizeAll=50
        samSizeAll=60
        if (F) {
            x=datGPU[,which(clinU$testUCSF500orFMI=="UCSF500")]
            gene1=rownames(x)[apply(x,1,function(x) {mean(x,na.rm=T)!=0})]
            x=datGPU[,which(clinU$testUCSF500orFMI=="FMI")]
            gene2=rownames(x)[apply(x,1,function(x) {mean(x,na.rm=T)!=0})]
            i=which(rownames(datGPU)%in%gene1[gene1%in%gene2])
        }
        gene1=unique(clin1$gene[clin1$assayVersion=="T5a"])
        gene2=unique(clin1$gene[clin1$assayVersion=="T7"])
        j=1:nrow(clin3)
        j=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7") & ((clin1$id%in%clin32$id[j]) | (!clin1$disOntTerm2%in%clin3$disOntTerm2[j])))
        i=which(rownames(datGPU)%in%clin1$gene[j])
        if (nameValue$value[which(nameValue$name=="subset")]=="ucsf500") {
            #i=which(rownames(datGPU)%in%gene1)
            i=1:nrow(datGPU)
            j=which(clinU$testUCSF500orFMI%in%c("UCSF500"))
        } else {
            j=which(clinU$testUCSF500orFMI%in%c("UCSF500","FMI"))
        }
        clinU=clinU[j,]
        datGP=datGPU[i,j]
        samSizeAll=ceiling(nrow(clinU)/10)*10
        
        tbl1=read.table("docs/ucsf500/gene/diff-v1-dropped.txt",sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
        tbl2=read.table("docs/ucsf500/gene/diff-v2-added.txt",sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
        tbl3=read.table("docs/ucsf500/gene/bad-genes.txt",sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
        x=unique(c(tbl1[,1],tbl2[,1],tbl3[,1]))
        datGP=datGP[which(!rownames(datGP)%in%x),]
    } else {
        gene1=unique(clin1$gene[clin1$assayVersion=="T5a"])
        gene2=unique(clin1$gene[clin1$assayVersion=="T7"])
        #if (subset2Flag%in%c("_allSamples_T5aT7Assays","_lungSmallJeffVettedHighGradeGI_T5aT7Assays")) {
        if (nameValue$value[which(nameValue$name=="subset")]%in%c("allSamples","jeffVetted")) {
            j=1:nrow(clin3)
        } else {
            j=which(!is.na(clin32$fionaCol1))
        }
        j=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7") & ((clin1$id%in%clin32$id[j]) | (!clin1$disOntTerm2%in%clin3$disOntTerm2[j])))
        #j=1:nrow(clin1)
        datGP=datGP[which(rownames(datGP)%in%clin1$gene[j]),]
        out=getFamilyLevelInfo()
        datGP_m=out$datGP_m
        candGene=out$candGene
    }
    
    geneFamilyFlag=T
    geneFamilyFlag=F
    
    #candGeneThis=candGene
    
    ## -------------------
    
    confIntFlag=T
    confIntFlag=F
    
    samSize=c(0,0)
    samSize=c(50,10)
    samSize=c(50,30)
    
    datadirG=""
    geneFlag=F
    geneFlag=T
    
    
    if (F) {
        genesetFlag=""; geneset2List=c("10%patient","15%patient","100topGene")
        genesetFlag=""; geneset2List=c("10%patientGrp1Grp2")
        genesetFlag="_actionableGene"; geneset2List=c("actionableGene")
        genesetFlag=""; geneset2List=c("10%patientGrp1Grp2")
        genesetFlag=""; geneset2List=c("10%patientDisease")
        genesetFlag=""; geneset2List=c("15%patient")
        genesetFlag=""; geneset2List=c("10%patient")
        genesetFlag=""; geneset2List=c("10topGene")
        genesetFlag=""; geneset2List=c("3patient")
        genesetFlag=""; geneset2List=c("20topGene")
        
        subset2Flag="_T5aT7Assays"
        
        subset2List=c("_lungSmallJeffVettedHighGradeGI_T5aT7Assays","_lungSmallFionaVettedGI_T5aT7Assays","_T5aT7Assays","_fionaVetted_T5aT7Assays","_withinSmallCell_T5aT7Assays","_lungSmallFionaVettedSmallCellGI_T5aT7Assays","_lungSmallSmallCellGI_T5aT7Assays","_lungSmallLargeCellGI_T5aT7Assays","_tmbLevel_T5aT7Assays","_tmbLevel_withinVetted_T5aT7Assays")
        subset2List=c("_lungSmallJeffVettedHighGradeGI_T5aT7Assays","_lungSmallFionaVettedGI_T5aT7Assays","_fionaVetted_T5aT7Assays","_lungSmallFionaVettedSmallCellGI_T5aT7Assays","_tmbLevel_withinVetted_T5aT7Assays")
        subset2List=sub("_T5aT7Assays","_noOtherGI_T5aT7Assays",subset2List)
        
        genesetFlag=""; geneset2List=c("100topGene"); subset2List=c("_fionaVetted_noOtherGI_T5aT7Assays_fionaVetted_noOtherGI_T5aT7Assays")
        geneFlag=F
        samSize=c(50,30)
    }
    
    if (F) {
        genesetFlag=""; geneset2List=paste("30topGeneIn",capWords(c("lungSmall","pancreas","colon")),sep="")
        subset2List=c("_lungSmallFionaVettedGI_T5aT7Assays")
        geneFlag=T
        geneFlag=F
        samSize=c(0,0)
        samSize=c(50,0)
        
        genesetFlag=""; geneset2List=paste("20topGeneIn",capWords(c("pancreas")),sep="")
        
        genesetFlag=""; geneset2List=paste("30topGeneIn",capWords(c("lungSmall","pancreas","colon","other")),sep="")
        genesetFlag="_actionableGene"; geneset2List=c("actionableGene")
        genesetFlag=""; geneset2List=c("10%patientGrp1Grp2")
        genesetFlag=""; geneset2List=c("10%patientDisease")
        genesetFlag=""; geneset2List=c("15%patient")
        genesetFlag=""; geneset2List=c("10%patient")
        genesetFlag=""; geneset2List=paste("20topGeneIn",capWords(c("lungSmall","pancreas","colon","other")),sep="")
        genesetFlag=""; geneset2List=paste("10topGeneIn",capWords(c("lungSmall","pancreas","colon","other")),sep="")
        genesetFlag=""; geneset2List=paste("3patientIn",capWords(c("lungSmall","pancreas","colon","other")),sep="")
        
        subset2List=c("_lungSmallJeffVettedHighGradeGI_T5aT7Assays")
        geneFlag=F
        samSize=c(0,0)
        subset2List=c("_lungSmallFionaVettedGI_noOtherGI_T5aT7Assays","_fionaVetted_noOtherGI_T5aT7Assays","_lungSmallFionaVettedSmallCellGI_noOtherGI_T5aT7Assays")
        subset2List=c("_lungSmallFionaVettedGI_noOtherGI_T5aT7Assays")
        subset2List=c("_lungSmallFionaVettedGI_T5aT7Assays")
        geneFlag=F
        samSize=c(50,30)
        samSize=c(50,0)
    }
    
    geneFlag=F
    samSize=c(50,0)
    if (nameValue$value[which(nameValue$name=="subset")]%in%c("ucsf500","ucsf500Fmi")) {
        samSize=c(0,0)
    }
    
    subset1Flag=""
    subset1Flag="_withLung"
    
    ## All genesets
    #genesetList=c("_swiSnfComp","_histoneMod","_swiSnfHisMod","_swiSnfPlusHisMod","_atmAtr","_rbEtc","_swiSnfEtc")
    ## For component level
    ## For gene level
    #genesetList=c("_swiSnf")
    genesetList=c("_swiSnfComp","_histoneMod","_swiSnfPlusHisMod","_atmAtr","_rbEtc")
    genesetList=c("_swiSnfComp","_histoneMod","_atmAtr","_swiSnfHisMod","_swiSnfPlusHisMod","_swiSnfHisModAtmAtr","_swiSnfPlusHisModPlusAtmAtr","_rbEtc")
    genesetList=c("_swiSnfHisMod","_swiSnfPlusHisMod","_swiSnfHisModAtmAtr","_swiSnfPlusHisModPlusAtmAtr")
    genesetList=c("_swiSnfHisModAtmAtr","_swiSnfPlusHisModPlusAtmAtr")
    genesetList=c("_20topGene")
    genesetList=c("_swiSnfHisModAtmAtr")
    genesetList=c("_swiSnfPlusHisModPlusAtmAtr")
    genesetList=c("_swiSnfEtc")
    genesetList=c("")
    genesetList=c("_30topGeneInUcsf500")
    genesetList=c("_30topGeneInUcsf500Fmi")
    genesetList=c("_30topGene")
    genesetList=c("_p53Etc")
    
    yLim=NULL
    yLim=c(0,60) # Used for genesetList=_p53Etc, UCSF 500, FMI with lung: Poorly differentiated, well differentiated
    
    
    for (genesetFlag in genesetList) {
        geneset2List=sub("_","",genesetFlag)
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
        },
        "_p53Etc"={geneset2List=c("p53Etc")
        }
        )
        
        for (subset2Flag in subset2List) {
            
            pThres=0.05
            
            if (length(grep("_ucsf500",subset2Flag))==1) {
                altTypeUniq1="anyAlt"
                altTypeUniq2=cbind(altTypeUniq1,"1")
                colList=c("grey")
                if (genesetFlag%in%c("_p53Etc")) {
                    candGeneThis=getCandidateGenes(genesetFlag)
                    write.table(candGeneThis,file=paste("geneList",genesetFlag,".txt",sep=""),append=F,col.names=T,row.names=F,sep="\t",quote=F)
                } else {
                    candGeneThis=NULL
                }
            } else {
                altTypeUniq1=sort(unique(clin1$altType))
                #altTypeUniq2=cbind(c(altTypeUniq1,paste(altTypeUniq1[1],"+",altTypeUniq1[3]),paste(altTypeUniq1[1:3],"+",altTypeUniq1[4])),as.character(c(1:length(altTypeUniq1),10*(1)+3,10*(1:3)+4)))
                altTypeUniq2=cbind(c(altTypeUniq1,paste(altTypeUniq1[1],"+",altTypeUniq1[2]),paste(altTypeUniq1[1],"+",altTypeUniq1[3]),paste(altTypeUniq1[1:3],"+",altTypeUniq1[4]),paste(altTypeUniq1[1],"+",altTypeUniq1[2],"+",altTypeUniq1[4])),as.character(c(1:length(altTypeUniq1),10*(1)+2:3,10*(1:3)+4,124)))
                altTypeUniq2=cbind(c(altTypeUniq1,"multiple"),as.character(c(1:length(altTypeUniq1),10)))
                if ("multiple"%in%altTypeUniq2[,1]) {
                    colList=c("red","blue","orange","cyan","yellow","indianred","yellow2","skyblue","bisque3","indianred4")
                } else {
                    colList=c("red","blue","purple","cyan","skyblue","orange","yellow","darkgreen")
                    colList=c("red","blue","indianred","cyan","orange","yellow","skyblue","bisque3")
                    colList=c("red","blue","indianred","cyan","yellow2","orange","yellow","skyblue","bisque3","indianred4")
                }
                candGeneThis=getCandidateGenes(genesetFlag)
                write.table(candGeneThis,file=paste("geneList",genesetFlag,".txt",sep=""),append=F,col.names=T,row.names=F,sep="\t",quote=F)
            }
            colListThis=colList
            
            altTypeUniq3="anyAlt"
            
            outFormat="png"
            outFormat="pdf"
            
            propFlag=">="
            propFlag=">"
            
            
            if (geneFamilyFlag) {
                datGPThis=datGP_m
            } else {
                datGPThis=datGP
            }
            
            if (length(grep("_ucsf500",subset2Flag))==1) {
                if (length(grep("_ucsf500Fmi",subset2Flag))==1) {
                    heading1=paste("UCSF 500, FMI",sep="")
                } else {
                    heading1=paste("UCSF 500",sep="")
                }
                if (length(grep("Diff",subset2Flag))==1) {
                    j=which(clinU$diffEB==sub("Diff","",strsplit(subset2Flag,"_")[[1]][2]))
                    clinThis=clinU[j,]
                    heading1=paste(heading1,": ",ifelse(clinThis$diffEB[1]=="poor","Poorly","Well")," differentiated",sep="")
                    clinThis$disease=clinThis$diffEB[1]
                    disUniq=cbind(c(""),paste(clinThis$diffEB[1],ifelse(clinThis$diffEB[1]=="poor","ly",""),c(" differentiated"),sep=""))
                    clinThis$disease=sub("_","",strsplit(subset2Flag,"Diff")[[1]][1])
                    disUniq=cbind(c(""),c("All samples"))
                    varList=c("disease")
                    varName=varList
                    if (!genesetFlag%in%c("_p53Etc")) {
                        if (any(!is.na(clinThis$cellSizeEB))) {
                            varList=c(varList,"cellSizeEB")
                        }
                        varName=varList
                        if (any(!is.na(clinThis$grade))) {
                            varList=c(varList,"gradeThis")
                            varName=c(varName,"grade")
                        }
                    }
                } else {
                    heading1=paste(heading1,": All samples",sep="")
                    clinThis=clinU
                    clinThis$disease="allSamples"
                    disUniq=cbind(c(""),c("All samples"))
                    varList=c("disease","diffEB","cellSizeEB")
                    varName=c("disease","diffEB","cellSizeEB")
                    #varList=c("disease","primarySiteEB2","diffEB","cellSizeEB")
                    #varName=c("disease","primarySiteEB","diffEB","cellSizeEB")
                    varList=c("disease","primarySiteEB2","diffEB","cellSizeEB","gradeThis")
                    varName=c("disease","primarySiteEB","diffEB","cellSizeEB","grade")
                }
                clinThis$primarySiteEB2=tolower(clinThis$primarySiteEB)
                clinThis$primarySiteEB2[which(clinThis$primarySiteEB2=="other(lung)")]="lung"
                clinThis$primarySiteEB2[grep("othergi(", clinThis$primarySiteEB2,fixed=T)]="otherGI"
                clinThis$primarySiteEB2[grep("other(", clinThis$primarySiteEB2,fixed=T)]="other"
                clinThis$primarySiteEB2[grep("unknown", clinThis$primarySiteEB2,fixed=T)]="unknown"
                clinThis$gradeThis=clinThis$grade
                clinThis$gradeThis[which(clinThis$grade=="grade2/3")]="gradeT"
            } else {
                gene1=unique(clin1$gene[clin1$assayVersion=="T5a"])
                gene2=unique(clin1$gene[clin1$assayVersion=="T7"])
                if (subset2Flag%in%c("_allSamples_T5aT7Assays","_lungSmallJeffVettedHighGradeGI_T5aT7Assays")) {
                    j=1:nrow(clin3)
                } else {
                    j=which(!is.na(clin32$fionaCol1))
                }
                j=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7") & ((clin1$id%in%clin32$id[j]) | (!clin1$disOntTerm2%in%clin3$disOntTerm2[j])))
                #j=1:nrow(clin1)
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
                } else if (subset2Flag=="_fionaVettedUcsf500") {
                    clinThis=clinThis[which(!is.na(clinThis$fionaCol1)),]
                    clinThis$id=paste("F",clinThis$id,sep="")
                    clinThis2=clinU; clinThis2$id=paste("U",clinThis2$id,sep="")
                    k=match(names(clinThis),names(clinThis2)); k1=which(!is.na(k)); k2=k[k1]
                    clinThis=rbind(clinThis[,k1],clinThis2[,k2])
                    rm(clinThis2)
                    heading1=paste("Fiona vettes & UCSF 500: Pancreas",sep="")
                    disUniq=cbind(c("pascreas"),c("Pancreas"))
                    varList=c("cellSize")
                    varName=c("cellSize")
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
            }
            if (ncol(disUniq)==2) disUniq=cbind(disUniq,rep("",nrow(disUniq)))
            
            for (dId in 1:nrow(disUniq)) {
                cat("\n\n================ disease: ",disUniq[dId,1]," ===============\n\n",sep="")
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
                if (subset1Flag=="") {
                    samIdAll=samIdAll[!clinThis$caseId[samIdAll]%in%samExclId]
                }
                for (varId in 1:length(varList)) {
                    cat("\n\n================ variable: ",varList[varId],"\n\n",sep="")
                    if (length(grep("_ucsf500",subset2Flag))==1) {
                        switch(varList[varId],
                        "disease"={
                            x=strsplit(subset2Flag,"Diff")[[1]]
                            if (length(x)==1) {
                                grpUniq3=c("allSamples")
                                grpUniq3=cbind(grpUniq3,"All samples","A")
                            } else {
                                grpUniq3=sub("_","",x[1])
                                grpUniq3=cbind(grpUniq3,paste(ifelse(grpUniq3=="poor","Poorly","Well")," diff",sep=""),toupper(substr(grpUniq3,1,1)))
                            }
                        },
                        "primarySiteEB2"={
                            grpUniq3=sort(unique(clinThis[,varList[varId]]))
                            grpUniq3=c("colorectal","other gi","other non-gi","pancreas","unknown")
                            grpUniq3=cbind(grpUniq3,grpUniq3,toupper(substr(grpUniq3,1,1)))
                        },
                        "diffEB"={grpUniq3=c("poor","well","nr")
                            grpUniq3=cbind(grpUniq3,c("Poor","Well","NR"),toupper(substr(grpUniq3,1,1)))
                        },
                        "cellSizeEB"={grpUniq3=c("small","large","nr")
                            grpUniq3=cbind(grpUniq3,c("Small","Large","NR"),toupper(substr(grpUniq3,1,1)))
                        },
                        "gradeThis"={
                            #grpUniq3=paste("grade",c("2","3"),sep="")
                            grpUniq3=paste("grade",c("1","2","T","3"),sep="")
                            grpUniq3=cbind(grpUniq3,paste("grade",c("1","2","2/3","3"),sep=""),c("1","2","T","3"))
                        }
                        )
                    } else {
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
                    }
                    x=table(clinThis[samIdAll,varList[varId]])
                    k1=match(grpUniq3[,1],names(x))
                    #k=match(grpUniq3[,1],names(x)); k1=k[!is.na(k)]
                    k=which(x[k1]>=samSize[2])
                    if (length(grep("_ucsf500",subset2Flag))==1) {
                    } else {
                        if (length(k)<2) next
                    }
                    if (length(k)==1) {
                        grpUniq3=matrix(grpUniq3[k,],nrow=1)
                        n=1
                        grpUniq32=matrix(nrow=n,ncol=2)
                        k1=k2=1
                        grpUniq32[k,1]=grpUniq3[k1,1]
                        grpUniq32[k,2]=grpUniq3[k2,1]
                    } else {
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
                    if (nrow(altTypeUniq2)==1) {
                        out4AltPropCombAll=matrix(out4AltPropCombAll,ncol=1,dimnames=list(names(out4AltPropCombAll),altTypeUniq2[,1]))
                    }
                    if (nrow(grpUniq3)==1) {
                        out4AltPropAll=matrix(out4AltPropAll,ncol=1,dimnames=list(names(out4AltPropAll),paste(rep(altTypeUniq2[,1],each=nrow(grpUniq3)),grpUniq3[,1],sep="_")))
                    }
                    
                    outPvAll=matrix(nrow=nrow(outPropAll),ncol=nrow(grpUniq32),dimnames=list(rownames(outPropAll),paste("pv_anyAlt_",rep(paste(grpUniq32[,2],grpUniq32[,1],sep="V"),each=1),sep="")))
                    outEstAll=matrix(nrow=nrow(outPropAll),ncol=nrow(grpUniq32),dimnames=list(rownames(outPropAll),paste("qv_anyAlt_",rep(paste(grpUniq32[,2],grpUniq32[,1],sep="V"),each=1),sep="")))
                    print(nrow(outPvAll))
                    timeStamp=Sys.time()
                    print(format(timeStamp, "%x %X"))
                    for (k in 1:nrow(grpUniq32)) {
                        j=which(clinThis[samIdThis,varList[varId]]%in%grpUniq32[k,])
                        if (length(j)!=0) {
                            for (a1 in which(altTypeUniq3=="anyAlt")) {
                                for (g1 in 1:nrow(outPvAll)) {
                                    x=table(clinThis[samIdThis[j],varList[varId]],outAnyAlt[g1,j])
                                    if (nrow(x)==2 && ncol(x)==2) {
                                        res=fisher.test(x)
                                        colId=paste("qv_",altTypeUniq3[a1],"_",grpUniq32[k,2],"V",grpUniq32[k,1],sep="")
                                        outEstAll[g1,colId]=paste(round(c(res$estimate,res$conf.int),1),collapse=",")
                                        #outEstAll[g1,colId]=paste(c(res$estimate,res$conf.int),collapse=",")
                                        colId=paste("pv_",altTypeUniq3[a1],"_",grpUniq32[k,2],"V",grpUniq32[k,1],sep="")
                                        outPvAll[g1,colId]=res$p.value
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
                                load(paste(datadirG,fName2,sep=""))
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
                            #j=1:nrow(clinThis)
                            j=1:nrow(clinThis[samIdThis,])
                            x2=strsplit(geneset2Flag,"patientIn")[[1]]
                            if (length(x2)==2) {
                                #j=which(tolower(clinThis$disOntTerm2)==tolower(x2[2]))
                                #nm=clinThis$disease[j][1]
                                j=which(tolower(clinThis$disOntTerm2[samIdThis])==tolower(x2[2]))
                                nm=clinThis$disease[samIdThis][j][1]
                            }
                            x=apply(out4Alt[,j],1,function(x) {sum(x!=0,na.rm=T)})
                            i=order(x,decreasing=T)
                            i=i[which(x[i]>=numThres)]
                            if (length(i)==0) {cat("No genes !!!\n\n"); next}
                            fName=paste("proportionTable_geneBy",capWords(varList[varId]),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),"_",geneset2Flag,fName1,subset2Flag,".txt",sep="")
                            heading=paste(heading,"Table of ",length(i)," genes with >= ",numThres," patients with alterations",ifelse(length(x2)==2,paste(" in ",nm,sep=""),""),sep="")
                        } else if (length(grep("topGeneIn",geneset2Flag))==1) {
                            fName2=strsplit(genesetFlag,"In")[[1]]; fName2=paste(fName2[1],"_",tolowerWords(fName2[2]),sep="")
                            fName2=paste("geneSampleId_geneBy",capWords(varList[varId]),fName2,subset1Flag,".RData",sep="")
                            load(paste(datadirG,fName2,sep=""))
                            x2=strsplit(geneset2Flag,"topGeneIn")[[1]]
                            i=match(geneId,rownames(outPropAll))
                            fName=paste("proportionTable_geneBy",capWords(varList[varId]),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),"_",geneset2Flag,fName1,subset2Flag,".txt",sep="")
                            heading=paste(heading,"Table of ",length(i)," genes with most alterations",ifelse(length(x2)==2,paste(" in ",toupper(x2[2]),sep=""),""),sep="")
                        } else if (length(grep("topGene",geneset2Flag))==1) {
                            #numThres=as.integer(sub("topGene","",geneset2Flag))
                            numThres=as.integer(strsplit(geneset2Flag,"topGene")[[1]][1])
                            #j=1:nrow(clinThis)
                            j=1:nrow(clinThis[samIdThis,])
                            x2=strsplit(geneset2Flag,"topGeneIn")[[1]]
                            if (length(x2)==2) {
                                #j=which(tolower(clinThis$disOntTerm2)==tolower(x2[2]))
                                #nm=clinThis$disease[j][1]
                                j=which(tolower(clinThis$disOntTerm2[samIdThis])==tolower(x2[2]))
                                nm=clinThis$disease[samIdThis][j][1]
                            }
                            x=apply(out4Alt[,j],1,function(x) {sum(x!=0,na.rm=T)})
                            i=order(x,decreasing=T)[1:numThres]
                            fName=paste("proportionTable_geneBy",capWords(varList[varId]),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),"_",geneset2Flag,fName1,subset2Flag,".txt",sep="")
                            heading=paste(heading,"Table of ",length(i)," genes with most alterations",ifelse(length(x2)==2,paste(" in ",nm,sep=""),""),sep="")
                        } else if (geneset2Flag%in%c("actionableGene","swiSnf","swiSnfComponent","histoneModifier","atmAtr","swiSnfHistoneModifier","swiSnfPlusHistoneModifier","swiSnfHistoneModifierAtmAtr","swiSnfPlusHistoneModifierPlusAtmAtr","p53Etc","rbEtc","swiSnfCompEtc","swiSnfEtc")) {
                            if (geneFamilyFlag) {
                                i=match(toupper(unique(candGeneThis$family2)),toupper(rownames(outPropAll)))
                            } else {
                                i=match(toupper(candGeneThis$gene),toupper(rownames(outPropAll)))
                            }
                            i=i[!is.na(i)]
                            fName=paste("proportionTable_geneBy",capWords(varList[varId]),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),genesetFlag,fName1,subset2Flag,".txt",sep="")
                            #heading=paste(heading,"Table of ",length(i)," genes with any alteration in a group within ",disUniq[dId,2],sep="")
                            heading=paste(heading,"Table of ",length(i)," gene",ifelse(geneFamilyFlag," family",""),ifelse(length(i)==1,"","s"),sep="")
                        } else {
                            i=1:nrow(outPropAll)
                            fName=paste("proportionTable_geneBy",capWords(varList[varId]),ifelse(disUniq[dId,1]=="","",paste("_within",capWords(disUniq[dId,1]),sep="")),ifelse(geneset2Flag=="","","_"),geneset2Flag,fName1,subset2Flag,".txt",sep="")
                            heading=paste(heading,"Table of ",length(i)," genes with any alteration",ifelse(disUniq[dId,1]=="","",paste(" in a group within ",disUniq[dId,2],sep="")),sep="")
                            
                        }
                        if (subset1Flag=="_withLung") {
                            fName=sub(".txt",paste(subset1Flag,".txt",sep=""),fName,fixed=T)
                            heading=sub(":"," with lung:",heading)
                        }
                        if (length(i)!=0) {
                            colListThis=colList
                            if (length(grep("_ucsf500",subset2Flag))==1) {
                                x=table(clinThis[samIdThis,varList[varId]][which(clinThis[samIdThis,varList[varId]]%in%grpUniq3[,1])])
                                x=x[match(grpUniq3[,1],names(x))]
                                colListThis=rev(gray((1:samSizeAll)/samSizeAll))[c(x,1)]
                            }
                            nmR=rownames(outPropAll)[i]
                            outPropThis=round(outPropAll[i,colnames(outPropAll)],2)
                            if (ncol(outPropAll)==1) outPropThis=matrix(outPropThis,ncol=1,dimnames=list(rownames(outPropAll)[i],colnames(outPropAll)))
                            if (length(i)>1) {
                                outPercThis=apply(outPropThis,c(1,2),function(x) {
                                    y=paste(" ",x*100,"% ",sep="")
                                    y
                                })
                                if (ncol(outPropAll)==1) outPercThis=matrix(outPercThis,ncol=1,dimnames=list(rownames(outPropThis),colnames(outPropThis)))
                            } else if (length(i)==1) {
                                outPropThis=matrix(outPropThis,nrow=1,dimnames=list(nmR,colnames(outPropAll)))
                                outPercThis=matrix(paste(" ",c(outPropThis)*100,"% ",sep=""),nrow=1,dimnames=list(nmR,colnames(outPropAll)))
                            }
                            i=match(rownames(outQvAll),rownames(outPropThis)); i1=which(!is.na(i)); i2=i[i1]
                            colnames(outPercThis)=paste("anyAlt_",colnames(outPercThis),sep="")
                            res=table(clinThis[samIdThis,varList[varId]][which(clinThis[samIdThis,varList[varId]]%in%grpUniq3[,1])])
                            for (grpId in 1:nrow(grpUniq3)) {
                                if (res[which(names(res)==grpUniq3[grpId,1])]<samSize[1]) next
                                k1=grep("qv_anyAlt_",colnames(outQvAll),fixed=T)
                                k1=k1[grep(grpUniq3[grpId,1],colnames(outQvAll)[k1],fixed=T)]
                                nm1=colnames(outQvAll)[k1]
                                nm11=nm1
                                nm1=sub("qv_anyAlt_","",sub(paste("V",grpUniq3[grpId,1],sep=""),"",sub(paste(grpUniq3[grpId,1],"V",sep=""),"",nm1,fixed=T),fixed=T))
                                nm21=colnames(outPercThis)[grep("anyAlt_",colnames(outPercThis),fixed=T)]
                                nm2=sub("anyAlt_","",colnames(outPercThis)[grep("anyAlt_",colnames(outPercThis),fixed=T)],fixed=T)
                                #x=" "
                                for (k in 1:length(k1)) {
                                    if (res[which(names(res)==nm1[k])]<samSize[1]) next
                                    k21=which(nm2==nm1[k])
                                    ii2=which(outQvAll[i1,nm11[k]]<pThres)
                                    ii=which(outPvAll[i1,sub("qv_","pv_",nm11[k])]<pThres)
                                    if (length(ii)!=0) {
                                        x=rep("",length(ii))
                                        if (confIntFlag) {
                                            x=paste(x,grpUniq3[grpId,3],"(",ifelse(length(ii2)!=0,outEstAll[i1[ii],nm11[k]],tolower(outEstAll[i1[ii],nm11[k]])),"),",sep="")
                                        } else {
                                            x=paste(x,ifelse(length(ii2)!=0,grpUniq3[grpId,3],tolower(grpUniq3[grpId,3])),sep="")
                                        }
                                        outPercThis[i2[ii],k21]=paste(outPercThis[i2[ii],k21],x,sep="")
                                    }
                                }
                            }
                            outPercThis=sub(",$","",outPercThis)
                            x=table(clinThis[samIdThis,varList[varId]][which(clinThis[samIdThis,varList[varId]]%in%grpUniq3[,1])])
                            x=x[order(x,decreasing=T)]
                            nm=names(x)
                            colnames(outPercThis)=sub("anyAlt_","",colnames(outPercThis))
                            outPercThis=outPercThis[,nm]
                            if (!is.matrix(outPercThis)) {
                                if (ncol(outPropAll)==1) {
                                    outPercThis=matrix(outPercThis,ncol=1,dimnames=list(nmR,colnames(outPropAll)))
                                } else {
                                    outPercThis=matrix(outPercThis,nrow=1,dimnames=list(nmR,colnames(outPropAll)))
                                }
                            }
                            nm=grpUniq3[match(colnames(outPercThis),grpUniq3[,1]),2]
                            nm=paste(c("gene",paste(nm," (N=",x,")",sep="")),collapse="\t")
                            tbl=cbind(gene=rownames(outPercThis),data.frame(outPercThis))
                            if (length(grep("_noOtherGI",subset2Flag))==1) {
                                heading=sub(", other GI","",heading)
                                nm=sub("+otherGI","",nm,fixed=T)
                            }
                            write.table(heading,file=fName,append=F,col.names=F,row.names=F,sep="\t",quote=F)
                            #write.table(paste("N = ",sum(x),sep=""),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                            write.table(paste(ifelse(varList[varId]=="gradeThis" & "gradeT"%in%grpUniq3[,1],"T = grade 2/3, ",""),"N = ",sum(x),sep=""),file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                            write.table("",file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                            write.table(nm,file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                            write.table(tbl,file=fName,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                            
                            ## -------------------
                            ## Long tail plot
                            
                            x=apply(out4AltPropCombAll,1,sum)
                            gId=match(rownames(outPropThis),rownames(out4AltPropCombAll))
                            if (F) {
                                ## For combined groups
                                fName2=sub("proportionTable_","longTailPlot_",sub(".txt","",fName,fixed=T))
                                heading2=sub("Table of ","Long tail plot of ",heading)
                                switch(outFormat,
                                "png"={png(paste(fName2,".png",sep=""))},
                                "pdf"={pdf(paste(fName2,".pdf",sep=""),width=2*7,height=7)}
                                )
                                #gId=gId[order(x[gId],decreasing=T)]
                                barplot(100*t(out4AltPropCombAll[gId,]),names.arg=rownames(out4AltPropCombAll)[gId],main=heading2,ylab="% alteration",col=colListThis,las=3,cex.axis=1.5,cex.names=2)
                                dev.off()
                            }
                            
                            fName2=sub("proportionTable_","longTailPlot2_",sub(".txt","",fName,fixed=T))
                            heading2=sub("Table of ","Long tail plot of ",heading)
                            heading2=sub(" \n"," - ",heading2)
                            #heading2=strsplit(heading2," \n")[[1]][1]
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
                                } else if (nrow(grpUniq3)>2 & length(gId)>10) {
                                    cexAxis=1
                                }
                                if (length(gId)>10) {
                                    pdf(paste(fName2,".pdf",sep=""),width=3*7,height=7)
                                    cexAxis2=.7
                                } else {
                                    pdf(paste(fName2,".pdf",sep=""),width=2*7,height=7)
                                    cexAxis2=1.2
                                }
                            }
                            )
                            
                            #cexAxis2=1
                            par(mar=c(8, 5, 0, 0) + 0.1)
                            par(mar=c(8.5, 5, 1, 0) + 0.1)
                            out=NULL
                            nm=c()
                            for (kk in gId) {
                                for (grpId in 1:nrow(grpUniq3)) {
                                    ll=grep(paste("_",grpUniq3[grpId,1],sep=""),colnames(out4AltPropAll),fixed=T)
                                    out=rbind(out,out4AltPropAll[kk,ll])
                                    #nm=c(nm,paste(ifelse(grpId==1,paste(rownames(out4AltPropAll)[kk]," ",sep=""),""),grpUniq3[grpId,2],sep=""))
                                    nm=c(nm,paste(grpUniq3[grpId,2],"    ",sep=""))
                                    #nm=c(nm,paste(c(grpUniq[grpId,2],rep("",length(gId)-1)),rownames(out4AltPropAll)[gId],sep=""))
                                    #barplot(100*t(out4AltPropAll[gId,ll]),names.arg=rownames(out4AltPropAll)[gId],main=ifelse(grpId==1,heading2,""),sub=grpUniq3[grpId,2],ylab="% alteration",col=colListThis,las=3,cex.axis=.7)
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
                            x=table(clinThis[samIdThis,varList[varId]][which(clinThis[samIdThis,varList[varId]]%in%grpUniq3[,1])])
                            if (ncol(out)==1) out2=100*c(out) else out2=100*t(out)
                            barplot(out2,names.arg=rownames(out),ylim=yLim,main=ifelse(grpId==1,heading2,""),ylab="% alteration",col=colListThis,las=3,cex.names=cexAxis,cex.axis=cexAxis,cex.lab=cexLab,space=0)
                            #axis(side=1,at=1:length(gId),labels=rownames(out4AltPropAll)[gId])
                            #axis(side=1,at=seq(1,nrow(out),by=nrow(grpUniq3)+1),labels=rownames(out4AltPropAll)[gId],tick=F,cex.axis=cexAxis2,cex.lab=cexLab)
                            #x=rownames(out4AltPropAll)[gId]; k=seq(1,length(x),by=2); x[k]=paste(x[k],"\n")
                            x=rownames(out4AltPropAll)[gId]
                            axis(side=1,at=seq(1,nrow(out),by=nrow(grpUniq3)+1),labels=x,tick=F,cex.axis=cexAxis2,cex.lab=cexLab)
                            dev.off()
                            
                            fName2=sub("proportionTable_","geneSampleId_",sub(".txt","",fName,fixed=T))
                            geneId=rownames(outPropThis)
                            sampleId=as.integer(sapply(colnames(outAnyAlt),function(x) strsplit(x,"_")[[1]][2],USE.NAMES=F))
                            #sampleId=clinThis$id[which(clinThis$id%in%sampleId & clinThis[,varList[varId]]%in%colnames(outPercThis))]
                            sampleId=clinThis$id[samIdThis][which(clinThis$id[samIdThis]%in%sampleId & clinThis[samIdThis,varList[varId]]%in%colnames(outPercThis))]
                            save(geneId,sampleId,file=paste(fName2,".RData",sep=""))
                        }
                    }
                }
            }
        }
    }
}

library(marray)
fName="_forLongTailPlot"
width = 480; height = 480
if (outFormat=="png") {
    png(paste("samSizeColorLegend",fName,".png",sep=""),width=width,height=height)
} else {
    pdf(paste("samSizeColorLegend",fName,".pdf",sep=""))
}
grpUniq=altTypeUniq2[,1]
cexThis=NULL
if (outFormat=="pdf") cexThis=1
colColUniq=rev(gray((1:samSizeAll)/samSizeAll))
lim=c(1,samSizeAll)
heatmapColorBar(limit=lim,cols=c(colColUniq[c(length(colColUniq),1,median(1:length(colColUniq)))]))
dev.off()

if (F) {
    if (F) {
        y="colon"
        x=outAnyAlt; nm=as.integer(sub("anyAlt_","",colnames(x)))
        x=x[match(geneId,rownames(x)),which(clinThis$disOntTerm2==y)]
        z1=as.integer(apply(x,1,function(x) {any(!is.na(x)) & any(x!=0,na.rm=T)}))
        z2=as.integer(apply(x,2,function(x) {any(!is.na(x)) & any(x!=0,na.rm=T)}))
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

    #rm(outAnyAlt,outPropAll,outPropThis,outPercThis,outPvAll,outQvAll,out4Alt,out4AltPropCombAll,geneId,sampleId)


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

    if (F) {
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
        j=which(clinThis$fionaCol1=="fionaY")
        table(clinThis$disOntTerm2[j],clinThis$cellSize[j])
        "
        BLANK LC SC SMLC
        colon        0 17 32    2
        other        0  1 10    6
        pancreas     1 13 70    7
        "
    }

}

if (F) {
    j=which(clin1$assayVersion%in%c("T5a","T7"))
    dim(datGPThis)
    table(rownames(datGPThis)%in%clin1$gene[j])
    dim(clinThis)
    table(clinThis$id%in%clin1$id[j])
    table(clinThis$disOntTerm2)

    i=which(rownames(datGP)%in%candGeneThis$gene)
    j=match(paste("anyAlt_",clinThis$id,sep=""),colnames(datGP))

    grpUniq=unique(candGeneThis$family2)
    for (k in 1:length(grpUniq)) {
        print(grpUniq[k])
        i=which(rownames(datGP)%in%candGeneThis$gene[which(candGeneThis$family2==grpUniq[k])])
        x=apply(datGP[i,],2,function(x) {as.integer(any(x!=0,na.rm=T))})
        i=which(rownames(datGP_m)%in%grpUniq[k])
        print(table(datGP_m[i,]==x,exclude=NULL))
    }
}
