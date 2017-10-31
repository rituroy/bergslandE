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


## Run analysis.R section 1 first
## ordFlag, clusterFlag, subset2List - set these parameters
## clusterGene - list of genes to cluster samples by. NULL - cluster based on all genes. Set annCol$grp
## ordFlag - sample groups to order by
## fracRead - not working, CHECK

genesetFlag="geneByFionaCol1_withinColon_15percPatient_fionaVetted_T5aT7Assays"
genesetFlag="geneByTmbLevel_15percPatient_tmbLevel_withinVetted_T5aT7Assays"
genesetFlag="geneByTmbLevel_10percPatient_tmbLevel_T5aT7Assays"
genesetFlag="geneByFionaCol1_withinPancreas_15percPatient_fionaVetted_T5aT7Assays"
genesetFlag="geneByDisOntTerm3_10percPatient_lungSmallSmallCellGI_T5aT7Assays"
genesetFlag="geneByCellSize_actionableGene_fionaVetted_noOtherGI_T5aT7Assays"

####################################################################
####################################################################

nameValueList=c("allSamples","jeffVetted","fionaVetted")
nameValueList=c("jeffVetted")
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
    }
    )
    write.table(nameValue,"config.tmp", sep="=", col.names=F, row.names=F, quote=F)
    
    colId2="etoeprt"
    colId2=c("gene","partnerGene","altType","cdsEffect","proteinEffect","cn","fracReads")
    
    datadir="results/"
    if (nameValue$value[which(nameValue$name=="subset")]%in%c("ucsf500","ucsf500Fmi")) {
        load(paste(datadir,"tmp_allAssays.RData",sep=""))
        clinF1=clin1
        clinF31=clin3
        clinF32=clin32
        rm(clin1,clin3,clin32)
        datadir="results/"
        load(paste(datadir,"tmp_ucsf500.RData",sep=""))
        nameValue=read.table("config.tmp",sep="=",h=F,quote="",comment.char="",as.is=T,fill=T,col.names=c("name","value"))
        #colId2=unique(c(colId2,c("currPatStatus","ageAtDiag","survSinceInitDiag","pathReview","typeOfPath","stageAtDiag","survSinceStageIV","diff","cellSize","grade","mitoticRate","ki67","note","testUCSF500orFMI","dateOfUCSFsampleUsed","chromosomalCopyChange","noMutDelRearrange","dead")))
        colId2=unique(c(colId2,c("currPatStatus","ageAtDiag","survSinceInitDiag","pathReview","typeOfPath","stageAtDiag","survSinceStageIV","diff","cellSize","mitoticRate","ki67","note","testUCSF500orFMI","dateOfUCSFsampleUsed","chromosomalCopyChange","noMutDelRearrange","dead")))
        clin3=clinU
        #names(clin3)[match(c("primarySite"),names(clin3))]="disOntTerm2"
        clin3$primarySiteEB2=tolower(clin3$primarySiteEB)
        clin3$primarySiteEB2[which(clin3$primarySiteEB2=="other(lung)")]="lung"
        clin3$primarySiteEB2[grep("othergi(", clin3$primarySiteEB2,fixed=T)]="otherGI"
        clin3$primarySiteEB2[which(clin3$primarySiteEB2=="other gi")]="otherGI"
        clin3$primarySiteEB2[grep("other(", clin3$primarySiteEB2,fixed=T)]="other"
        clin3$primarySiteEB2[which(clin3$primarySiteEB2=="other non-gi")]="other"
        clin3$primarySiteEB2[grep("unknown", clin3$primarySiteEB2,fixed=T)]="unknown"
        #clin3$disOntTerm2=clin3$primarySiteEB2
        if (T) {
            clin3$cellSizeEB[which(clin3$cellSizeEB=="small")]="SC"
            clin3$cellSizeEB[which(clin3$cellSizeEB=="large")]="LC"
            clin3$cellSizeEB[which(clin3$cellSizeEB=="nr")]="NR"
        }
        if (F) {
            clin3$cellSize=clin3$cellSizeEB
            clin3$cellSize[which(clin3$cellSize=="small")]="SC"
            clin3$cellSize[which(clin3$cellSize=="large")]="LC"
            clin3$cellSize[which(clin3$cellSize=="nr")]="NR"
        }
    } else {
        load(paste(datadir,"tmp_allAssays.RData",sep=""))
        clinF1=clin1
        clinF31=clin3
        clinF32=clin32
    }
    annColThis=clin3[,which(!names(clin3)%in%colId2)]
    
    candGene1=c("TP53","RB1","MLL2","APC","CDKN2A","RAS-MAPK","mTOR","MMR","Other")
    
    clusterGene=NULL
    
    clusterFlag=""; lineFlag=F
    
    clusterFlag="_comutated"; lineFlag=T
    clusterFlag="_comutatedGenes"; lineFlag=F
    
    clusterFlag="_supervised"; lineFlag=T
    #clusterFlag=c(clusterFlag,"_topGenes"); lineFlag=T
    
    #k=which(nameValue$name=="subset")
    datadir="results/"
    if (nameValue$value[which(nameValue$name=="subset")]%in%c("ucsf500","ucsf500Fmi")) {
        load(paste(datadir,"tmp_ucsf500.RData",sep=""))
        nameValue=read.table("config.tmp",sep="=",h=F,quote="",comment.char="",as.is=T,fill=T,col.names=c("name","value"))
        subsetFlag=paste("_",nameValue$value[which(nameValue$name=="subset")],sep="")
        clin3=clinU
        #names(clin3)[match(c("primarySite"),names(clin3))]="disOntTerm2"
        clin3$primarySiteEB2=tolower(clin3$primarySiteEB)
        clin3$primarySiteEB2[which(clin3$primarySiteEB2=="other(lung)")]="lung"
        clin3$primarySiteEB2[grep("othergi(", clin3$primarySiteEB2,fixed=T)]="otherGI"
        clin3$primarySiteEB2[which(clin3$primarySiteEB2=="other gi")]="otherGI"
        clin3$primarySiteEB2[grep("other(", clin3$primarySiteEB2,fixed=T)]="other"
        clin3$primarySiteEB2[which(clin3$primarySiteEB2=="other non-gi")]="other"
        clin3$primarySiteEB2[grep("unknown", clin3$primarySiteEB2,fixed=T)]="unknown"
        #clin3$disOntTerm2=clin3$primarySiteEB2
        if (F) {
            clin3$diffEB[which(clin3$diffEB=="nr")]="1:NR"
            clin3$diffEB[which(clin3$diffEB=="poor")]="2:poor"
            clin3$diffEB[which(clin3$diffEB=="well")]="3:well"
        }
        clin3$diffEB[which(clin3$diffEB=="nr")]="NR"
        
        if (T) {
            clin3$cellSizeEB[which(clin3$cellSizeEB=="nr")]="NR"
            clin3$cellSizeEB[which(clin3$cellSizeEB=="small")]="SC"
            clin3$cellSizeEB[which(clin3$cellSizeEB=="large")]="LC"
        }
        if (F) {
            clin3$cellSize=clin3$cellSizeEB
            clin3$cellSize[which(clin3$cellSize=="nr")]="NR"
            clin3$cellSize[which(clin3$cellSize=="small")]="SC"
            clin3$cellSize[which(clin3$cellSize=="large")]="LC"
        }
    } else {
        if (nameValue$value[which(nameValue$name=="subset")]=="fionaVetted") {
            #load(paste(datadir,"tmp_fiona166_T5aT7Assays.RData",sep=""))
            load(paste(datadir,"tmp_allAssays.RData",sep=""))
            nameValue=read.table("config.tmp",sep="=",h=F,quote="",comment.char="",as.is=T,fill=T,col.names=c("name","value"))
            subsetFlag=c("_fiona166_T5aT7Assays","_lungSmallFionaVettedGI_T5aT7Assays")
        } else {
            load(paste(datadir,"tmp_allAssays.RData",sep=""))
            nameValue=read.table("config.tmp",sep="=",h=F,quote="",comment.char="",as.is=T,fill=T,col.names=c("name","value"))
            gene1=unique(clin1$gene[clin1$assayVersion=="T5a"])
            gene2=unique(clin1$gene[clin1$assayVersion=="T7"])
            subsetFlag="_allSamples_T5aT7Assays"
            samId=which(clin1$gene%in%gene1[gene1%in%gene2] & clin1$assayVersion%in%c("T5a","T7"))
            samId03=which(clin3$id%in%clin1$id[samId])
            #k=which(nameValue$name=="subset")
            if (nameValue$value[which(nameValue$name=="subset")]=="jeffVetted") {
                subsetFlag=c("_jeffVetted_T5aT7Assays","_lungSmallJeffVettedHighGradeGI_T5aT7Assays")
                samId03=samId03[which(clin3$disOntTerm3[samId03]=="lungSmall" | (clin3$grade[samId03]=="High" & clin3$disOntTerm3[samId03]=="combined"))]
                samId=samId[which(clin1$id[samId]%in%clin3$id[samId03])]
            }
            samId01=samId
        }
    }
    if (length(subsetFlag)==1) subsetFlag=rep(subsetFlag,2)
    if (F) {
        #if (clusterFlag[1]=="") {
        #if (clusterFlag[1]%in%c("","_comutated","_comutatedGenes")) {
        if (clusterFlag[1]%in%c("","_supervised","_comutated","_comutatedGenes")) {
            datadir="results/"
            ##load(paste(datadir,"tmp_fiona102_T5aT7Assays.RData",sep=""))
            ##load(paste(datadir,"tmp_T5aT7Assays.RData",sep="")); subsetFlag[1]="_T5aT7Assays"
            #load(paste(datadir,"tmp_allAssays.RData",sep=""))
            load(paste(datadir,"tmp_fiona166_T5aT7Assays.RData",sep=""))
        }
    }
    
    datTypeFlag=c("_fracRead","")
    datTypeFlag=c("","")
    
    subset1Flag=""
    subset1Flag="_withLung"
    
    ordFlag=""
    ordFlag="_swiSnfEtc"
    ordFlag="_tmbSwiSnfEtc"
    ordFlag="_diffEBki67EB"
    ordFlag="_diffEBpriSiteEBki67EB"
    ordFlag="_priSiteEBki67EB"
    
    outFormat="png"
    outFormat="pdf"
    
    source(paste(dirSrc3,"funcs.R",sep=""))
    if (nameValue$value[which(nameValue$name=="subset")]%in%c("ucsf500","ucsf500Fmi")) {
        if (F) {
            x=datGPU[,which(clinU$testUCSF500orFMI=="UCSF500")]
            gene1=rownames(x)[apply(x,1,function(x) {mean(x,na.rm=T)!=0})]
            x=datGPU[,which(clinU$testUCSF500orFMI=="FMI")]
            gene2=rownames(x)[apply(x,1,function(x) {mean(x,na.rm=T)!=0})]
            i=which(rownames(datGPU)%in%gene1[gene1%in%gene2])
        }
        gene1=unique(clinF1$gene[clinF1$assayVersion=="T5a"])
        gene2=unique(clinF1$gene[clinF1$assayVersion=="T7"])
        j=1:nrow(clinF31)
        j=which(clinF1$gene%in%gene1[gene1%in%gene2] & clinF1$assayVersion%in%c("T5a","T7") & ((clinF1$id%in%clinF32$id[j]) | (!clinF1$disOntTerm2%in%clinF31$disOntTerm2[j])))
        i=which(rownames(datGPU)%in%clinF1$gene[j])
        if (nameValue$value[which(nameValue$name=="subset")]=="ucsf500") {
            #i=which(rownames(datGPU)%in%gene1)
            i=1:nrow(datGPU)
            j=which(clinU$testUCSF500orFMI%in%c("UCSF500"))
        } else {
            j=which(clinU$testUCSF500orFMI%in%c("UCSF500","FMI"))
        }
        clinU=clinU[j,]
        clin3=clin3[j,]
        datGP=datGPU[i,j]
        #rm(datGP_m,candGene)
        if (subset1Flag=="") {
            j=which(!clinU$caseId%in%samExclId)
            clinU=clinU[j,]
            clin3=clin3[j,]
            datGP=datGP[,j]
        }
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
        if (subset1Flag=="") {
            j=which(!clinU$caseId%in%samExclId)
            clinU=clinU[j,]
            clin3=clin3[j,]
            datGP=datGP[,j]
            datGP_m=datGP_m[,j]
        }
    }
    
    library(marray)
    #source(paste(dirSrc,"functions/heatmap.5.R",sep=""))
    #source(paste(dirSrc,"functions/heatmapAcgh.7.R",sep=""))
    #source(paste(dirSrc,"functions/heatmap.5.4.R",sep=""))
    #source(paste(dirSrc,"functions/heatmapAcgh.7.1.R",sep=""))
    #source(paste(dirSrc,"functions/heatmapAcgh.7.2.R",sep=""))
    source(paste(dirSrc,"functions/heatmapAcgh.7.3.R",sep=""))
    #source(paste("heatmap.5.5.R",sep=""))
    #source(paste(dirSrc,"functions/heatmap.5.5.R",sep=""))
    source(paste(dirSrc,"functions/heatmap.5.6.R",sep=""))
    
    datadir2="results/"
    
    asc <- function(x) { strtoi(charToRaw(x),16L) }
    chr <- function(n) { rawToChar(as.raw(n)) }
    limPerc=c(0,100)
    
    colId2="etoeprt"
    colId2=c("gene","partnerGene","altType","cdsEffect","proteinEffect","cn","fracReads")
    if (subsetFlag[1]%in%c("_ucsf500","_ucsf500Fmi")) {
        #colId2=unique(c(colId2,c("currPatStatus","ageAtDiag","survSinceInitDiag","pathReview","typeOfPath","stageAtDiag","survSinceStageIV","diff","cellSize","grade","mitoticRate","ki67","note","testUCSF500orFMI","dateOfUCSFsampleUsed","chromosomalCopyChange","noMutDelRearrange","dead")))
        colId2=unique(c(colId2,c("currPatStatus","ageAtDiag","survSinceInitDiag","pathReview","typeOfPath","stageAtDiag","survSinceStageIV","diff","cellSize","mitoticRate","ki67","note","testUCSF500orFMI","dateOfUCSFsampleUsed","chromosomalCopyChange","noMutDelRearrange","dead")))
    }
    
    if (datTypeFlag[1]=="_fracRead" | nameValue$value[which(nameValue$name=="subset")]%in%c("ucsf500","ucsf500Fmi","fionaVettedUcsf500")) {
        altTypeUniq1="anyAlt"
        altTypeUniq2=cbind(altTypeUniq1,"1")
        #colList=c("grey")
    } else {
        altTypeUniq1=sort(unique(clin1$altType))
        altTypeUniq2=cbind(c(altTypeUniq1,paste(altTypeUniq1[1],"+",altTypeUniq1[3]),paste(altTypeUniq1[1:3],"+",altTypeUniq1[4])),as.character(c(1:length(altTypeUniq1),10*(1)+3,10*(1:3)+4)))
        altTypeUniq2=cbind(c(altTypeUniq1,"multiple"),as.character(c(1:length(altTypeUniq1),10)))
    }
    colList=c("skyblue","blue","yellow","purple","red")
    colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","magenta","purple","darkgreen")
    colList=c("brown","red","orange","yellow","green","cyan","skyblue","blue","pink","peachpuff","purple","darkgreen","limegreen","salmon","gray","gold","antiquewhite","steelblue","aquamarine","lightcyan","turquoise","hotpink","black")
    corListGene=gray((0:4)/4)
    colList2=c("skyblue","blue")
    colList2=c("skyblue","blue","yellow","purple","red","green","cyan")
    colList2=c("gray","black","green","yellow","purple","red","cyan")
    colList2=c("gray","black","green","yellow","purple","red","cyan")
    colHMCon=list("red","blue","grey")
    
    distMethod="kendall"
    distMethod="pearson"
    distMethod="spearman"
    distMethod="euclidean"
    distMethod="cosine"
    distMethodList=c("euclidean","cosine")
    distMethodList=c("euclidean") # Used this for FMI data
    distMethodList=c("cosine") # Used this for UCSF500 data because CDKN2A & CDKN2B were not clustering together with euclidean distance
    
    linkMethod="ward.D2"
    
    png("tmp.png")
    barplot(1:length(colList),col=colList)
    dev.off()
    if ("multiple"%in%altTypeUniq2[,1]) {
        colHMCat=list(c("red","blue","orange","cyan","yellow","indianred","yellow2","skyblue","bisque3","indianred4"),NULL,NULL)
    } else if ("anyAlt"%in%altTypeUniq2[,1]) {
        colHMCat=list(c("grey"),NULL,NULL)
    } else {
        #colHMCat=list(c("red","blue","purple","cyan"),NULL,NULL)
        #colHMCat=list(c("red","blue","purple","cyan","skyblue","orange","yellow","darkgreen"),NULL,NULL)
        colHMCat=list(c("red","blue","indianred","cyan","yellow2","orange","yellow","skyblue","bisque3","indianred4"),NULL,NULL)
    }
    
    annCol=annColThis
    rm(annColThis)
    if (!subsetFlag[1]%in%c("_ucsf500","_ucsf500Fmi")) {
        if (exists("clin32")) {
            annCol=cbind(annCol,clin32[,which(!names(clin32)%in%names(annCol))])
            annCol$cellSize[which(annCol$disOntTerm2=="lungSmall")]="SC"
            annCol$cellSize[which(annCol$cellSize%in%c("MC","SC/LC","SC/MC"))]="other"
            annCol$cellSize[which(!annCol$cellSize%in%c("SC","LC","other"))]=NA
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
    }
    annColAll=annCol
    
    colListD=colList
    k=which(names(annColAll)%in%c("disOntTerm2","primarySiteEB2"))[1]
    x=sort(unique(annColAll[,k]))
    colListD[which(x%in%c("lungSmall","lung"))]="blue"
    colListD[which(x=="pancreas")]="magenta"
    colListD[which(x%in%c("colon","colorectal"))]="purple"
    colListD[which(x%in%c("otherGI","other gi"))]="pink"
    colListD[which(x%in%c("other","other non-gi"))]="lightcyan"
    colListD[which(x=="unknown")]="gold"
    cat("x: ",x,"\n")
    cat("colListD: ",colListD,"\n")
    
    colListC=colList
    if (subsetFlag[1]%in%c("_ucsf500","_ucsf500Fmi")) {
        x=annColAll$cellSizeEB
        #x=annColAll$cellSize
    } else {
        x=annColAll$cellSize
    }
    x=sort(unique(x))
    colListC[which(x=="SC")]="palegreen"
    colListC[which(x=="LC")]="darkgreen"
    colListC[which(x=="other")]="yellowgreen"
    #colListC[which(x=="NR")]="yellow"
    colListC[which(x=="NR")]="grey"
    
    colListT=colList
    if (F) {
        if (("tmbLevel")%in%names(annColAll)) {
            x=sort(unique(annColAll$tmbLevel))
            colListT[which(x=="Low")]="darkorange1"
            colListT[which(x=="Intermediate")]="orange1"
            colListT[which(x=="High")]="darkorange4"
        }
    }
    
    colListDi=c("grey","skyblue","blue")
    x=annColAll$diffEB
    x=sort(unique(x))
    colListDi[which(x=="poor")]="skyblue"
    colListDi[which(x=="well")]="blue"
    colListDi[which(x=="NR")]="grey"
    
    colListGr=c("yellow","orange","darkorange1","darkorange3")
    
    colListGene=grey((0:4)/4)
    colListGene=grey((0:3)/3)
    
    if (datTypeFlag[1]%in%c("_fracRead")) {
        altTypeList="_anyAlt"
    } else {
        altTypeFlag="_4Alt"
        altTypeList=c("_anyAlt","_anyAlt1minus1","_4Alt")
        altTypeList=c("_4Alt")
    }
    
    #subsetFlag=""
    #subsetFlag="_allAssays"
    #if (length(subsetFlag)==1) subsetFlag=rep(subsetFlag,2)
    
    subset2Flag=""
    subset2Flag="_giUnknown"
    subset2Flag="_tmbHigh"
    
    datTypeFlag[2]="_prop"
    datTypeFlag[2]="_binary"
    datTypeFlag[2]="_patient"
    
    subset2Flag=""; genesetFlag="_10percPatient"
    subset2Flag="_giUnknown"; genesetFlag="_10percPatient"
    subset2Flag=""; genesetFlag="_setd2"
    subset2Flag=""; genesetFlag="_notch"
    subset2Flag=""; genesetFlag="_mmr"
    subset2Flag=""; genesetFlag="_brca"
    subset2Flag=""; genesetFlag="_allGenes"
    
    subset3Flag=rep("",2)
    subset3Flag[1]="_no0AltGene"
    subset3Flag[1]=""
    
    subset4Flag=""
    
    subsetFFlag=""
    
    mutFlag=""
    mutFlag="_numMut"
    
    subsetName=subsetFlag[1]
    if (subsetFlag[1]=="") {
        heading1=paste("T5a/T7 assays",sep="")
        subsetName="_T5aT7Assays"
        subset2Flag="_giUnknown"
        subset2Flag="_noLung"
        genesetList=c("_10percPatient","_setd2","_notch","_mmr","_brca")
        genesetList=c("_setd2","_notch","_mmr","_brca")
    } else {
        genesetList=c("_setd2","_notch","_mmr","_brca")
        genesetList=c("_notch")
        genesetList=c("_setd2","_notch","_mmr","_brca","_allGenes")
        genesetList=c("_allGenes")
        switch(subsetFlag[1],
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
            
            heading1=paste("Fiona certified pancreas, colon",sep="")
            subset2Flag="_gi"
            
            heading1=paste("Fiona vetted colon, pancreas & otherGI + SCLC",sep="")
            heading1=paste("Fiona vetted",sep="")
            subset2Flag="_giLungSmall"
        },
        "_jeffVetted_T5aT7Assays"={
            heading1=paste("Jeff vetted colon, pancreas & otherGI + SCLC",sep="")
            heading1=paste("Jeff vetted",sep="")
            subset2Flag="_giLungSmall"
        },
        "_allSamples_T5aT7Assays"={
            heading1=paste("All colon, pancreas & otherGI + SCLC",sep="")
            heading1=paste("All samples",sep="")
            subset2Flag="_giLungSmall"
        },
        "_ucsf500"={
            heading1=paste("UCSF 500: All samples",sep="")
            subset2Flag=subsetFlag[1]
        },
        "_ucsf500Fmi"={
            heading1=paste("UCSF 500, FMI: All samples",sep="")
            subset2Flag=subsetFlag[1]
        }
        )
    }
    #fName1=""
    fName1=datTypeFlag[1]
    switch(datTypeFlag[1],
    "_fracRead"={
        heading1=paste("FracRead: ",heading1,sep="")
    }
    )
    geneFamilyList=NULL
    if (clusterFlag[1]=="_supervised") {
        if (F) {
            datadirG="geneSampleId/"
            #genesetList=genesetList[grep(subsetFlag[1],genesetList)]
            genesetList=gsub("geneSampleId_","",sub(".RData","",genesetList,fixed=T))
            k=grep("geneByTmbLevel|_withinPancreas|_withinColon|_withinOther",genesetList)
            if (length(k)!=0) {
                genesetList=genesetList[-k]
            }
            #heading1=paste("Fiona certified colon, pancreas & otherGI + SCLC",sep="")
            subset2Flag=""
            subset3Flag[1]="_no0AltGene"
            genesetList="geneByCellSize_actionableGene_fionaVetted_noOtherGI_T5aT7Assays"
            geneFamilyFlag=F
            geneFamilyFlag=T
            
            subset2Flag=c("_other")
            subset2Flag=c("_colon")
            subset2Flag=c("_pancreas")
            subset2Flag=c("_lungSmall")
            subset3Flag[1]=""
            datadirG="results/3pat/geneSampleId/"
            genesetList=dir(datadirG); genesetList=sub(".RData","",sub("geneSampleId_","",genesetList),fixed=T)
            genesetList=genesetList[grep(paste("patientIn",capWords(sub("_","",subset2Flag)),sep=""),genesetList)]
            geneFamilyFlag=T
            geneFamilyFlag=F
            clusterGene=NULL
            #heading1=""
            
            fName1=paste("_men1",subsetFlag[1],sep="")
            #heading1=paste(heading1,", genes co-mutated with men1")
            subset2Flag=c("_pancreas")
            subset3Flag[1]=""
            datadirG="results/comutation/fisherTestGeneLevelPv.05/anyCellSize/allSamplesSamSize0/men1Comut/"
            genesetList="pancreas"
            geneFamilyFlag=T
            geneFamilyFlag=F
            clusterGene=NULL
            
            fName1=""
            heading1=paste(heading1,", SWI/SNF...",sep="")
            datadirG="results/15percPatInAnyDisease/geneSampleId/"
            datadirG="results/15percPatInADisease/geneSampleId/"
            datadirG="results/20topGene/geneSampleId/"
            datadirG="results/30topGene/30topGene/geneSampleId/"
            datadirG="results/30topGene/genesCommonForDataset/geneSampleId/"
            datadirG="results/30topGene/geneSampleId/"
            #k=which(nameValue$name=="subset")
            genesetList="geneByDisease_20topGeneInPancreas_lungSmallJeffVettedHighGradeGI_T5aT7Assays"
            genesetList=dir(datadirG); genesetList=sub(".RData","",sub("geneSampleId_","",genesetList),fixed=T)
            genesetList=dir(datadirG); genesetList=genesetList[grep(tolower(nameValue$value[which(nameValue$name=="subset")]),tolower(genesetList))]; genesetList=sub(".RData","",sub("geneSampleId_","",genesetList),fixed=T)
            subset2Flag=c("_lungSmall")
            subset2Flag=c("_other")
            subset2Flag=c("_colon")
            subset2Flag=c("_pancreas")
            subset3Flag[1]=""
            geneFamilyFlag=T
            geneFamilyFlag=F
            clusterGene=NULL
            
            fName1=fName1
            heading1=paste(heading1,", SWI/SNF, Histone modifier, ATM/ATR",sep="")
            if (datTypeFlag[1]%in%c("_fracRead")) heading1=sub("Histone modifier","Hist mod",heading1)
            geneFamilyFlag=F
            datadirG=paste("results/swiSnf_modifier/geneSampleId/",ifelse(geneFamilyFlag,"componentLevel","geneLevel"),"/",sep="")
            #k=which(nameValue$name=="subset")
            genesetList=dir(datadirG,pattern="geneSampleId_geneByDisease"); genesetList=genesetList[grep(tolower(nameValue$value[which(nameValue$name=="subset")]),tolower(genesetList))]; genesetList=sub(".RData","",sub("geneSampleId_","",genesetList),fixed=T)
            genesetList=genesetList[grep("_swiSnfHisModAtmAtr",genesetList)]
            geneFamilyList=c("SWI.SNF","Modifier.Histone","ATM.ATR")
            subset2Flag=c("_lungSmall")
            subset2Flag=c("_colon")
            subset2Flag=c("_other")
            subset2Flag=c("_pancreas")
            subset3Flag[1]=""
            geneFamilyFlag=T
            geneFamilyFlag=F
            clusterGene=NULL
            
            geneFamilyFlag=F
            fName1=""
            heading1=paste(heading1,", SWI/SNF...",sep="")
            datadirG=paste("results/swiSnf_modifier/extra/",ifelse(geneFamilyFlag,"componentLevel","geneLevel"),"/",sep="")
            #k=which(nameValue$name=="subset")
            genesetList=dir(datadirG,pattern="geneSampleId_geneByDisease"); genesetList=genesetList[grep(tolower(nameValue$value[which(nameValue$name=="subset")]),tolower(genesetList))]; genesetList=sub(".RData","",sub("geneSampleId_","",genesetList),fixed=T)
            subset2Flag=c("_lungSmall")
            subset2Flag=c("_colon")
            subset2Flag=c("_other")
            subset2Flag=c("_pancreas")
            subset3Flag[1]=""
            geneFamilyFlag=T
            geneFamilyFlag=F
            clusterGene=NULL
        }
        
        fName1=""
        datadirG=paste("results/ucsf500/",ifelse(subset1Flag=="_withLung","withLung","noLung"),"/",sub("_","",subsetFlag[1]),"/geneSampleId/",sep="")
        genesetList=dir(datadirG,pattern="geneSampleId_geneByDisease"); genesetList=genesetList[grep(tolower(nameValue$value[which(nameValue$name=="subset")]),tolower(genesetList))]; genesetList=sub(".RData","",sub("geneSampleId_","",genesetList),fixed=T)
        datadirG=paste("results/ucsf500/p53Etc/",ifelse(subset1Flag=="_withLung","withLung","noLung"),"/",sub("_","",subsetFlag[1]),"/geneSampleId/",sep="")
        genesetList=dir(datadirG,pattern="geneSampleId_geneByDisease"); genesetList=genesetList[grep(tolower(nameValue$value[which(nameValue$name=="subset")]),tolower(genesetList))]; genesetList=sub(".RData","",sub("geneSampleId_","",genesetList),fixed=T)
        subset2Flag=c("")
        subset3Flag[1]=""
        geneFamilyFlag=F
        clusterGene=NULL
        
    } else if (clusterFlag[1]%in%c("_comutated")) {
        subset2Flag=c("_lungSmall","_pancColon") ## Heatmap appears grey for _lungSmall since there are too many cases with no mutation
        subset2Flag=c("_pancColon")
        subset2Flag=c("_lungSmall")
        subset2Flag=c("_colon")
        subset2Flag=c("_pancreas")
        switch(subsetFlag[1],
        "_fiona166_T5aT7Assays"={
            subset4Flag="_smallCell"
            subset4Flag=""
            subset4Flag="_fionaVetted" ## same as subset4Flag=""
        },
        "_jeffVetted_T5aT7Assays"={
            subset4Flag="_jeffVetted"
        },
        "_allSamples_T5aT7Assays"={
            subset4Flag="_allSamples"
        }
        )
        subset3Flag[1]="_no0AltGene"
        if (subset4Flag[1]=="_jeffVetted") {
            datadirG="results/comutation/fisherTestGeneLevelPv.05/anyCellSize/jeffVettedSamSize0/"
        } else if (subset4Flag[1]=="_allSamples") {
            datadirG="results/comutation/fisherTestGeneLevelPv.05/anyCellSize/allSamplesSamSize0/"
        } else if (subset4Flag[1]=="_smallCell") {
            datadirG="results/comutation/fisherTestGeneLevelPv.05/smallCell/"
            datadirG="results/comutation/fisherTestGeneLevelPv.01/smallCell/"
        } else {
            #datadirG="results/comutation/wrong/"
            datadirG="results/comutation/cosine0.5/"
            datadirG="results/comutation/fisherTestGeneLevelPv.05/anyCellSize/"
            datadirG="results/comutation/fisherTestGeneLevelPv.01/anyCellSize/"
            datadirG="results/comutation/fisherTestGeneLevelPv.05/anyCellSize/samSize0/"
        }
        #subset4Flag=c(subset4Flag,"_smallLargeCell")
        subset4Flag=c(subset4Flag,"")
        subsetFFlag="_men1"
        subsetFFlag=""
        load(file=paste(datadirG,"geneList",subset1Flag,".RData",sep=""))
        
        geneFamilyFlag=T
        geneFamilyFlag=F
        
        datadirG=paste("results/swiSnf_modifier/comutation/",ifelse(geneFamilyFlag,"componentLevel","geneLevel"),"/",sep="")
        geneset2Flag="_swiSnfEtc"
        load(file=paste(datadirG,"geneList",geneset2Flag,subsetFlag[2],subset1Flag,".RData",sep=""))
        
        genesetList=names(geneList)
        genesetList=genesetList[which(genesetList==paste("_",sep=""))]
        
        testFlag="cosine"
        testFlag="fisher"
    } else if (clusterFlag[1]%in%c("_comutatedGenes")) {
        altTypeList="_anyAlt"
        datTypeFlag[2]="_gene"
        subset2Flag=""
        switch(subsetFlag[1],
        "_fiona166_T5aT7Assays"={
            subset4Flag="_smallCell"
            subset4Flag=""
            subset4Flag="_fionaVetted" ## same as subset4Flag=""
        },
        "_jeffVetted_T5aT7Assays"={
            subset4Flag="_jeffVetted"
        },
        "_allSamples_T5aT7Assays"={
            subset4Flag="_allSamples"
        }
        )
        if (subset4Flag[1]=="_jeffVetted") {
            datadirG="results/comutation/fisherTestGeneLevelPv.05/anyCellSize/jeffVettedSamSize0/"
        } else if (subset4Flag[1]=="_allSamples") {
            datadirG="results/comutation/fisherTestGeneLevelPv.05/anyCellSize/allSamplesSamSize0/"
        } else if (subset4Flag[1]=="_smallCell") {
            datadirG="results/comutation/fisherTestGeneLevelPv.05/smallCell/"
            datadirG="results/comutation/fisherTestGeneLevelPv.01/smallCell/"
        } else {
            #datadirG="results/comutation/wrong/"
            datadirG="results/comutation/fisherTestGeneLevelPv.05/anyCellSize/"
            datadirG="results/comutation/cosine0.5/"
            datadirG="results/comutation/fisherTestGeneLevelPv.01/anyCellSize/"
            datadirG="results/comutation/fisherTestGeneLevelPv.05/anyCellSize/samSize0/"
        }
        subsetFFlag="_men1"
        subsetFFlag=""
        load(file=paste(datadirG,"geneList",subset1Flag,".RData",sep=""))
        
        geneFamilyFlag=F
        geneFamilyFlag=T
        
        datadirG=paste("results/swiSnf_modifier/genesCommonForDataset/comutation/",ifelse(geneFamilyFlag,"componentLevel","geneLevel"),"/",sep="")
        geneset2Flag="_swiSnfEtc"
        datadirG=paste("results/swiSnf_modifier/comutation/",ifelse(geneFamilyFlag,"componentLevel","geneLevel"),"/",sep="")
        geneset2Flag="_swiSnfCompEtc"
        geneset2Flag="_swiSnfEtc"
        geneset2Flag="_swiSnfHisModAtmAtr"
        load(file=paste(datadirG,"geneList",geneset2Flag,subsetFlag[2],subset1Flag,".RData",sep=""))
        
        genesetList=names(geneList)
        geneset2List=geneset2Flag
        
        testFlag="cosine"
        testFlag="fisher"
        colHMCon=list("darkred","indianred","red")
    }
    if (subset1Flag=="_withLung") {
        fName1=paste(fName1,subset1Flag,sep="")
        heading1=sub(":"," with lung:",heading1)
    }
    subset2List=c("_lungSmall","_pancreas","_colon","_other")
    subset2List=subset2Flag
    subset2List=c("_poorDiff","_wellDiff")
    
    getCorFlag=T
    getCorFlag=F
    
    if (length(subset4Flag)==1) subset4Flag=c(subset4Flag,subset4Flag)
    
    if (clusterFlag[1]=="_supervised") {
        #if (length(grep("percPatient",genesetList))!=0) geneBar="clusterPr"
        geneBar=""
        geneBar="clusterPr"
        sampleBar=""
        sampleBar="cluster"
    } else if (clusterFlag[1]=="_comutated") {
        sampleBar=""
        sampleBar="cluster"
        geneBar="clusterPr"
    } else if (clusterFlag[1]=="_comutatedGenes") {
        sampleBar=""
        geneBar=""
    } else {
        sampleBar=""
        sampleBar="cluster"
        geneBar=""
        geneBar="clusterPr"
    }
    
    
    #genesetList="colon"
    #genesetList="pancreas"
    #genesetList="lungSmall"
    
    
    
    for (altTypeFlag in altTypeList) {
        for (distMethod in distMethodList) {
            for (subset2Flag in subset2List) {
                for (genesetFlag in genesetList) {
                    candGeneThis=getCandidateGenes()
                    if (!geneFamilyFlag & length(grep("_actionableGene",genesetFlag))==1) {
                        candGeneThis=getCandidateGenes("_actionableGene")
                    }
                    if (clusterFlag[1]=="_comutated") {
                        if (subset2Flag=="_lungSmall") {
                            geneset2List=paste("_",candGene1,sep="")
                            geneset2List=""
                        } else {
                            geneset2List=paste("_",candGene1,sep="")
                            geneset2List=""
                        }
                    } else if (clusterFlag[1]=="_comutatedGenes") {
                        #geneset2List=c("_swiSnfEtc")
                    } else {
                        if (clusterFlag[1]=="_supervised") {
                            if (length(grep("InPancreas",genesetFlag))==1) {
                                clusterGene=list("TP53","RB1","CDKN2A","MEN1","ATRX","DAXX",c("MEN1","ATRX","DAXX"))
                            } else if (length(grep("InColon",genesetFlag))==1) {
                                if (geneFamilyFlag) {
                                    clusterGene=list("TP53","RB1","APC","RAS-MAPK",c("APC","RAS-MAPK"))
                                } else {
                                    clusterGene=list("TP53","RB1","APC","KRAS",c("APC","KRAS"))
                                }
                            } else if (length(grep("InLungSmall",genesetFlag))==1) {
                                clusterGene=list("TP53","RB1","MLL2","LRP1B")
                            } else if (length(grep("_swiSnfEtc",genesetFlag))==1) {
                                geneset2List=c("_swiSnfEtc")
                                candGeneThis=getCandidateGenes(geneset2List)
                                if (geneFamilyFlag) {
                                    clusterGene=list("TP53","RB1","APC","RAS-MAPK",c("APC","RAS-MAPK"))
                                    clusterGene=NULL
                                } else {
                                    clusterGene=list(candGeneThis$gene[which(candGeneThis$family2=="SWI.SNF")],candGeneThis$gene[which(candGeneThis$family2=="Modifier.Histone")],c("ATRX","DAXX"),"RB1")
                                }
                            } else {
                                clusterGene=list("RB1","TP53")
                                clusterGene=NULL
                            }
                            geneset2List=""
                            if (!is.null(clusterGene)) sampleBar="cluster"
                        }
                    }
                    for (geneset2Flag in geneset2List) {
                        
                        
                        annCol=clin3[,which(!names(clin3)%in%colId2)]
                        if (!subsetFlag[1]%in%c("_ucsf500","_ucsf500Fmi")) {
                            if (exists("clin32")) {
                                annCol=cbind(annCol,clin32[,which(!names(clin32)%in%names(annCol))])
                                annCol$cellSize[which(annCol$disOntTerm2=="lungSmall")]="SC"
                                annCol$cellSize[which(annCol$cellSize%in%c("MC","SC/LC","SC/MC"))]="other"
                                annCol$cellSize[which(!annCol$cellSize%in%c("SC","LC","other"))]=NA
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
                        }
                        #annColAll=annCol
                        
                        if (datTypeFlag[1]=="_fracRead") {
                            if (geneFamilyFlag) {
                                datGPThis=datGDThis=NULL
                            } else {
                                datGPThis=fracFP
                                datGDThis=NULL
                            }
                        } else {
                            if (geneFamilyFlag) {
                                datGPThis=datGP_m
                                datGDThis=NULL
                            } else {
                                datGPThis=datGP
                                if (exists("datGD")) datGDThis=datGD
                            }
                        }
                        
                        nClust=c(1,5)
                        nClust=c(5,5)
                        nClust=c(NA,5)
                        nClust=c(NA,NA)
                        nClust=c(NA,1)
                        nClust=c(3,3)
                        
                        varList=varName=NULL
                        varFList=varFName=NULL
                        
                        fNameOut=paste(datTypeFlag[2],subsetName,subset2Flag,subset4Flag[2],subset3Flag[1],subset3Flag[2],mutFlag,"_",distMethod,altTypeFlag,sep="")
                        fNameOut2=paste(datTypeFlag[1],datTypeFlag[2],sep="")
                        limit1=c(0,1)
                        getClustInfoFlag=T
                        
                        cat("\n\n======================",fNameOut,genesetFlag,geneset2Flag,"==========================\n\n",sep="")
                        
                        if (clusterFlag[1]=="_supervised") {
                            x=strsplit(genesetFlag,"_")[[1]]
                            #fNameOut=paste("_",paste(x[2:length(x)],collapse="_"),sep="")
                            fNameOut=paste("_",genesetFlag,fName1,sep="")
                            if (subset1Flag!="" & length(grep(subset1Flag,fNameOut))==2) fNameOut=sub(subset1Flag,"",fNameOut)
                            if (subset2Flag!="" ) {
                                #if (length(grep(sub("_","",subset2Flag),fNameOut))==0) {
                                if (length(grep(paste(subset2Flag,"_",sep=""),fNameOut))==0) {
                                    fNameOut=sub("_geneByDisease",paste("_geneByDisease",subset2Flag,sep=""),fNameOut)
                                }
                            }
                            varList=tolowerWords(sub("geneBy","",x[1]))
                            varList=varList[varList%in%names(annCol)]
                            if (subsetFlag[1]%in%c("_ucsf500","_ucsf500Fmi")) {
                                varList=c("primarySiteEB2","diffEB","ki67EB2","cellSizeEB","grade")
                                varName=unique(names(clin3))
                            } else {
                                varList[which(varList=="disOntTerm3")]="disOntTerm2"
                                varList[which(varList=="disease")]="disOntTerm2"
                                if (!any(c("disease","disOntTerm","disOntTerm2","disOntTerm3")%in%varList)) varList=cbind(varList,"disOntTerm2")
                                #if (subset2Flag%in%c("_pancreas","_colon","_pancColon","_other")) varList=c(varList,"cellSize")
                                if (subset2Flag%in%c("_lungSmall","_pancreas","_colon","_pancColon","_other")) varList=c(varList,"cellSize","tmbLevel")
                                #if (subset2Flag%in%c("_lungSmall","_pancreas","_colon","_pancColon","_other")) varList=c(varList,"cellSize")
                                if (exists("clin32")) {
                                    varName=unique(c(names(clin3),names(clin32)))
                                } else {
                                    varName=unique(names(clin3))
                                }
                            }
                            varName=varName[match(varList,varName)]
                            k=which(varName%in%c("disease","disOntTerm","disOntTerm2","disOntTerm3"))
                            if (length(k)==1) {
                                varName[k]="disease"
                            }
                            k=match(c("fionaCol1","tmbLevel","primarySiteEB2","ki67EB2"),varName); k2=which(!is.na(k)); k1=k[k2]
                            varName[k1]=c("vetted","TMB","primarySiteEB","ki67EB")[k2]
                            if (subsetFlag[1]%in%c("_ucsf500","_ucsf500Fmi")) {
                                k=match(c("disease"),varName); k2=which(!is.na(k)); k1=k[k2]
                                varName[k1]=c("primarySiteEB")[k2]
                            }
                            varName=paste(varName," ",sep="")
                            varList[which(varList=="disease")]="disOntTerm2"
                            
                            if (!geneFamilyFlag & length(grep("_actionableGene",genesetFlag))==1) {
                                varFList="geneFamily"
                                varFName="family "
                            }
                            if (!is.null(geneFamilyList)) {
                                varFList="geneFamily"; varFName="family "
                            }
                            
                            load(paste(datadirG,"geneSampleId_",genesetFlag,".RData",sep=""))
                            if(length(grep("_fionaVetted|_withinVetted",genesetFlag))==1) {
                                heading=paste("Fiona vetted ",ifelse(length(grep("within",x[2]))==1,sub("within","",tolower(x[2])),"colon, pancreas & otherGI"),sep="")
                            } else {
                                heading=paste("Fiona checked ",ifelse(length(grep("within",x[2]))==1,sub("within","",tolower(x[2])),"colon, pancreas & otherGI + SCLC"),sep="")
                            }
                            if (F) {
                                heading=ifelse(length(grep("FionaVetted",genesetFlag))==1,"Group 2","Group 1")
                                heading=paste(heading,ifelse(length(grep("SmallCellGI",genesetFlag))==1,", small cell",""),sep="")
                                heading=paste(heading," pancreas + colon",ifelse(length(grep("_noOtherGI",genesetFlag))==1,""," + otherGI"),sep="")
                            }
                            heading=""
                            heading=heading1
                            x=strsplit(genesetFlag,"_")[[1]]
                            k=grep("percPatient",x)
                            if (length(k)!=0) {
                                heading=x[k]
                            }
                            k=grep("topGene",x)
                            if (length(k)!=0) {
                                fNameOut=paste("_patient_fiona166_T5aT7Assays_lungSmall_no0AltGene_euclidean_4Alt_genesComutatedInColon",genesetFlag,sep="")
                                fNameOut=paste("_patient",subset2Flag,"_",paste(x[2:5],collapse="_"),sep="")
                                fNameOut=paste("_patient",subset2Flag,"_",paste(x[2:length(x)],collapse="_"),sep="")
                                heading=x[k]
                                if (subsetFlag[1]%in%c("_ucsf500","_ucsf500Fmi")) heading=paste(heading1,", ",heading,sep="")
                            } else {
                                k=grep("patient",x)
                                if (length(k)!=0) {
                                    #fNameOut=paste("_patient",subset2Flag,"_",paste(x[2:5],collapse="_"),sep="")
                                    heading=x[k]
                                }
                            }
                            g0=match(geneId,rownames(datGPThis))
                            g0=g0[!is.na(g0)]
                            if (length(clusterFlag)>1 && clusterFlag[2]=="_topGenes") {
                                g0=match(unique(unlist(clusterGene)),rownames(datGPThis))
                            }
                            if (subsetFlag[1]%in%c("_ucsf500","_ucsf500Fmi")) {
                                samId3=which(clin3$id%in%sampleId)
                            } else {
                                if (subset2Flag=="") {
                                    #samId3=which(clin3$id%in%sampleId)
                                    samId3=which(clin3$id%in%sampleId & clin3$disOntTerm3=="combined")
                                } else {
                                    samId3=which(clin3$id%in%sampleId)
                                    #samId3=which(clin3$id%in%clin1$id[samId])
                                }
                            }
                            switch(subset2Flag,
                            "_pancColon"={
                                samId3=samId3[which(clin3$disOntTerm2[samId3]%in%c("pancreas","colon") & clin32$fionaCol1[samId3]=="Y")]
                                heading=paste(heading,": Pancreas+colon samples",sep="")
                            },
                            "_pancreas"={
                                samId3=samId3[which(clin3$disOntTerm2[samId3]%in%c("pancreas"))]
                                heading=paste(heading,": Pancreas samples",sep="")
                            },
                            "_colon"={
                                samId3=samId3[which(clin3$disOntTerm2[samId3]%in%c("colon"))]
                                heading=paste(heading,": Colon samples",sep="")
                            },
                            "_other"={
                                samId3=samId3[which(clin3$disOntTerm2[samId3]%in%c("other"))]
                                heading=paste(heading,": Other samples",sep="")
                            },
                            "_lungSmall"={
                                samId3=samId3[which(clin3$disOntTerm2[samId3]%in%c("lungSmall"))]
                                heading=paste(heading,": SCLC samples",sep="")
                            }
                            )
                            if (length(grep("Diff",subset2Flag))==1) {
                                samId3=samId3[which(clin3$diffEB[samId3]==sub("Diff","",strsplit(subset2Flag,"_")[[1]][2]))]
                                heading=paste(sub(": All samples","",heading),": ",ifelse(clin3$diffEB[1]=="poor","Poorly","Well")," differentiated",sep="")
                            }
                            if (subset4Flag[2]=="_smallCell") {
                                samId3=samId3[which(clin3$disOntTerm2[samId3]=="lungSmall" | clin32$cellSize[samId3]=="SC")]
                            } else if (subset4Flag[2]=="_smallLargeCell") {
                                samId3=samId3[which(clin3$disOntTerm2[samId3]=="lungSmall" | clin32$cellSize[samId3]%in%c("SC","LC"))]
                            }
                            if (length(g0)==0 | length(samId3)==0) {
                                cat("(length(g0)==0 | length(samId3)==0)\n")
                                next
                            }
                            #varList="disOntTerm2"
                            #varName="disease "
                            
                            if (is.na(varList) || varList=="NA") stop()
                        } else if (clusterFlag[1]=="_comutated") {
                            fNameOut=paste(fNameOut,"_genesComutatedIn",capWords(genesetFlag),capWords(sub("_","",subset4Flag[1])),geneset2Flag,sep="")
                            heading=paste(sub("_","",subset4Flag[1]),": Genes co-mutated",ifelse(geneset2Flag=="","",paste(" with ",sub("_","",geneset2Flag),sep="")),ifelse(subsetFFlag=="","",paste(" with ",toupper(sub("_","",subsetFFlag)),sep=""))," in ",ifelse(genesetFlag=="lungSmall","SCLC",genesetFlag),sep="")
                            g0=match(geneList[[genesetFlag]],rownames(datGPThis))
                            if (subsetFFlag!="") {
                                nm=toupper(sub("_","",subsetFFlag))
                                #tbl=read.table(paste(datadirG,ifelse(testFlag=="cosine","cor","comut"),"_",genesetFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                                tbl=read.table(paste(datadirG,ifelse(testFlag=="cosine","cor","comut"),"_",genesetFlag,geneset2Flag,subsetFlag[2],".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                                x=as.matrix(tbl[,-1])
                                rownames(x)=tbl$gene
                                x=x[which(rownames(x)==nm),which(colnames(x)!=nm)]
                                if (length(x)==0) {
                                    cat("No ",nm," gene!!!\n")
                                    next
                                } else {
                                    x=names(x)[which(x!=0)]
                                    if (length(x)==0) {
                                        cat("No genes co-mutated with ",nm,"!!!\n")
                                        next
                                    }
                                    g0=match(c(nm,x),rownames(datGPThis))
                                }
                            }
                            samId3=which(clin3$id%in%clin1$id[samId])
                            switch(subset2Flag,
                            "_pancColon"={
                                samId3=samId3[which(clin3$disOntTerm2[samId3]%in%c("pancreas","colon") & clin32$fionaCol1[samId3]=="Y")]
                                heading=paste(heading,": Pancreas+colon samples",sep="")
                            },
                            "_pancreas"={
                                samId3=samId3[which(clin3$disOntTerm2[samId3]%in%c("pancreas"))]
                                heading=paste(heading,": Pancreas samples",sep="")
                            },
                            "_colon"={
                                samId3=samId3[which(clin3$disOntTerm2[samId3]%in%c("colon"))]
                                heading=paste(heading,": Colon samples",sep="")
                            },
                            "_lungSmall"={
                                samId3=samId3[which(clin3$disOntTerm2[samId3]%in%c("lungSmall"))]
                                heading=paste(heading,": SCLC samples",sep="")
                            }
                            )
                            if (subset4Flag[2]=="_smallCell") {
                                samId3=samId3[which(clin3$disOntTerm2[samId3]=="lungSmall" | clin32$cellSize[samId3]=="SC")]
                            } else if (subset4Flag[2]=="_smallLargeCell") {
                                samId3=samId3[which(clin3$disOntTerm2[samId3]=="lungSmall" | clin32$cellSize[samId3]%in%c("SC","LC"))]
                            }
                            if (length(g0)==0 | length(samId3)==0) {
                                cat("(length(g0)==0 | length(samId3)==0)\n")
                                next
                            }
                            varList="disOntTerm2"
                            varName="disease "
                            varList=c("disOntTerm2")
                            #if (subset2Flag%in%c("_lungSmall","_pancreas","_colon","_pancColon","_other")) varList=c(varList,"cellSize","tmbLevel")
                            if (subset2Flag%in%c("_lungSmall","_pancreas","_colon","_pancColon","_other")) varList=c(varList,"cellSize")
                            varName=paste(varList," ",sep="")
                            varName[which(varList=="disOntTerm2")]="disease "
                        } else if (clusterFlag[1]=="_comutatedGenes") {
                            fNameOut=paste(fNameOut,"_genesComutatedIn",capWords(genesetFlag),capWords(sub("_","",subset4Flag[1])),sep="")
                            #heading=paste(sub("_","",subset4Flag[1]),": Correlation coefficient for genes co-mutated",ifelse(subsetFFlag=="","",paste(" with ",toupper(sub("_","",subsetFFlag)),sep=""))," in ",genesetFlag," ",sep="")
                            heading=paste(sub("_","",subset4Flag[1]),": Genes co-mutated",ifelse(subsetFFlag=="","",paste(" with ",toupper(sub("_","",subsetFFlag)),sep=""))," in ",genesetFlag," ",sep="")
                            g0=match(geneList[[genesetFlag]],rownames(datGPThis))
                            if (subsetFFlag!="") {
                                nm=toupper(sub("_","",subsetFFlag))
                                #tbl=read.table(paste(datadirG,ifelse(testFlag=="cosine","cor","comut"),"_",genesetFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                                tbl=read.table(paste(datadirG,ifelse(testFlag=="cosine","cor","comut"),"_",genesetFlag,geneset2Flag,subsetFlag[2],".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                                x=as.matrix(tbl[,-1])
                                rownames(x)=tbl$gene
                                x=x[which(rownames(x)==nm),which(colnames(x)!=nm)]
                                if (length(x)==0) {
                                    cat("No ",nm," gene!!!\n")
                                    next
                                } else {
                                    x=names(x)[which(x!=0)]
                                    if (length(x)==0) {
                                        cat("No genes co-mutated with ",nm,"!!!\n")
                                        next
                                    }
                                    g0=match(c(nm,x),rownames(datGPThis))
                                }
                            }
                            samId3=which(clin3$id%in%clin1$id[samId])
                            switch(paste("_",genesetFlag,sep=""),
                            "_pancColon"={
                                samId3=samId3[which(clin3$disOntTerm2[samId3]%in%c("pancreas","colon") & clin32$fionaCol1[samId3]=="Y")]
                            },
                            "_pancreas"={
                                samId3=samId3[which(clin3$disOntTerm2[samId3]%in%c("pancreas"))]
                            },
                            "_colon"={
                                samId3=samId3[which(clin3$disOntTerm2[samId3]%in%c("colon"))]
                            },
                            "_lungSmall"={
                                samId3=samId3[which(clin3$disOntTerm2[samId3]%in%c("lungSmall"))]
                            }
                            )
                            if (subset4Flag[2]=="_smallCell") {
                                samId3=samId3[which(clin3$disOntTerm2[samId3]=="lungSmall" | clin32$cellSize[samId3]=="SC")]
                            } else if (subset4Flag[2]=="_smallLargeCell") {
                                samId3=samId3[which(clin3$disOntTerm2[samId3]=="lungSmall" | clin32$cellSize[samId3]%in%c("SC","LC"))]
                            }
                            if (length(g0)==0 | length(samId3)==0) {
                                cat("(length(g0)==0 | length(samId3)==0)\n")
                                next
                            }
                        } else {
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
                            if (subset4Flag[2]=="_smallCell") {
                                samId3=samId3[which(clin3$disOntTerm2[samId3]=="lungSmall" | clin32$cellSize[samId3]=="SC")]
                            } else if (subset4Flag[2]=="_smallLargeCell") {
                                samId3=samId3[which(clin3$disOntTerm2[samId3]=="lungSmall" | clin32$cellSize[samId3]%in%c("SC","LC"))]
                            }
                            if (length(g0)==0 | length(samId3)==0) {
                                cat("(length(g0)==0 | length(samId3)==0)\n")
                                next
                            }
                            
                            x=strsplit(genesetFlag,"_")[[1]]
                            k=grep("percPatient",x)
                            if (length(x)==1) {
                                fNameOut=paste(fNameOut,genesetFlag,sep="")
                                k1=grep("anyAlt_",colnames(datGDThis))
                                nm=sub("anyAlt_","",colnames(datGDThis)[k1])
                                k1=k1[which(nm%in%clin3$disOntTerm2[samId3])]
                                g0=c()
                                propThres=as.integer(gsub("_|percPatient","",x[k]))/100
                                for (k in grep("anyAlt_",colnames(datGDThis))) {
                                    g0=c(g0,which(round(datGDThis[,k],2)>propThres))
                                }
                                g0=sort(unique(g0))
                                heading=paste(heading,", genes with alterations in > ",propThres*100,"% of patients in any group",sep="")
                                
                            } else {
                                switch(genesetFlag,
                                "_10percPatient"={
                                    fNameOut=paste(fNameOut,genesetFlag,sep="")
                                    k1=grep("anyAlt_",colnames(datGDThis))
                                    nm=sub("anyAlt_","",colnames(datGDThis)[k1])
                                    k1=k1[which(nm%in%clin3$disOntTerm2[samId3])]
                                    g0=c()
                                    propThres=0.10
                                    for (k in grep("anyAlt_",colnames(datGDThis))) {
                                        g0=c(g0,which(round(datGDThis[,k],2)>propThres))
                                    }
                                    g0=sort(unique(g0))
                                    heading=paste(heading,", genes with alterations in > 10% of patients in any group",sep="")
                                },
                                "_signif"={
                                    fName=paste("proportionTable_geneByDisease_qv",pThres,subsetName,".txt",sep="")
                                    heading="Table of significant genes in any of the comparisons pancreas vs. SCLC, other vs. SCLC or colon vs. SCLC"
                                    k1=grep("qv_anyAlt_",colnames(out1))
                                    k1=k1[grep("lungSmall",colnames(out1)[k1])]
                                    k1=k1[-grep("combined",colnames(out1)[k1])]
                                    nm1=colnames(out1)[k1]
                                    g0=c()
                                    for (k in k1) {
                                        g0=c(g0,which(out1[,k]<pThres))
                                    }
                                    g0=sort(unique(g0))
                                },
                                "_targetableGene"={
                                    fName=paste("proportionTable_geneByDisease_targetGene",subsetName,".txt",sep="")
                                    heading="Table of targetable genes"
                                    k1=grep("qv_anyAlt_",colnames(out1))
                                    k1=k1[grep("lungSmall",colnames(out1)[k1])]
                                    k1=k1[-grep("combined",colnames(out1)[k1])]
                                    nm1=colnames(out1)[k1]
                                    g0=match(toupper(candGene$gene),toupper(rownames(out1)))
                                    g0=g0[which(!is.na(g0))]
                                },
                                "_mmr"={
                                    fNameOut=paste(fNameOut,genesetFlag,sep="")
                                    tbl=read.table(paste("docs/candGene.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                                    tbl=tbl[which(tbl$family=="Mismatch excision repair (MMR)"),]
                                    x=unique(c(tbl$gene,unlist(sapply(tbl$alias,function(x) {strsplit(x,", ")[[1]]},USE.NAMES=F))))
                                    g0=which(toupper(rownames(datGPThis))%in%x)
                                    heading=paste(heading,", MMR genes",sep="")
                                },
                                "_notch"={
                                    fNameOut=paste(fNameOut,genesetFlag,sep="")
                                    g0=grep("NOTCH",toupper(rownames(datGPThis)))
                                    heading=paste(heading,", NOTCH genes",sep="")
                                },
                                "_brca"={
                                    fNameOut=paste(fNameOut,genesetFlag,sep="")
                                    x=c("BRCA1","BRCA2","ATM","ATR","PTEN","PARP1","MRE11A","PALB2","CDK12","BAP1","PRKDC","BARD1","BACH1","FANCD2","FANCE","CHEK1","CHEK2","PMS2","FANCA","FANCC","FANCF","XPA","NBN","BLM","BRIP1","RAD51B","RAD51C","RAD51D","RAD50")
                                    g0=which(toupper(rownames(datGPThis))%in%x)
                                    heading=paste(heading,", BRCA genes",sep="")
                                },
                                "_setd2"={
                                    fNameOut=paste(fNameOut,genesetFlag,sep="")
                                    g0=1:nrow(datGP)
                                    heading=paste(heading,", genes with alt, SETD2 marked",sep="")
                                },
                                "_allGenes"={
                                    fNameOut=paste(fNameOut,genesetFlag,sep="")
                                    g0=1:nrow(datGP)
                                    heading=paste(heading,", all genes",sep="")
                                }
                                )
                            }
                        }
                        if (length(g0)==0) {
                            cat("(length(g0)==0)\n")
                            next
                        }
                        colHM=colHMCon
                        genePvMat=geneQvMat=NULL
                        switch(datTypeFlag[2],
                        "_prop"={
                            arrayData=propMat[,grep("propMut_",colnames(propMat))]
                            
                            x1=apply(arrayData,1,function(x) mean(!is.na(x)))
                            x2=apply(arrayData,2,function(x) mean(!is.na(x)))
                            j=which(x2>=.5)
                            arrayData=arrayData[,j]
                            annRow=annCol=NULL
                            
                            x=apply(arrayData,1,function(x) sum(!duplicated(x[!is.na(x)])))
                            g0=which(x>1)
                            arrayData=arrayData[g0,]
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
                            #arrayData=datGPThis[g0,grep("anyAlt_",colnames(datGPThis))]
                            #colnames(arrayData)=sub("anyAlt_","",colnames(arrayData))
                            if (T) {
                                jj=grep(altTypeUniq1[1],colnames(datGPThis))
                                arrayData=matrix(0,nrow=length(g0),ncol=length(jj))
                                colnames(arrayData)=sub(paste(altTypeUniq1[1],"_",sep=""),"",colnames(datGPThis)[jj])
                                rownames(arrayData)=rownames(datGPThis)[g0]
                                #if (length(altTypeUniq1)==1) {
                                #} else {
                                if (F) {
                                    for (a1 in 2:length(altTypeUniq1)) {
                                        arrayData[datGPThis[g0,jj]==1]=a1
                                    }
                                }
                                if (F) {
                                    arrayData1=arrayData
                                    for (a1 in 2:length(altTypeUniq1)) {
                                        arrayData1[datGPThis[g0,jj]==1]=a1
                                    }
                                }
                                if (datTypeFlag[1]%in%c("_fracRead")) {
                                    arrayData=datGPThis[g0,]
                                    colnames(arrayData)=sub(paste(altTypeUniq1[1],"_",sep=""),"",colnames(datGPThis)[jj])
                                } else {
                                    id1=sapply(colnames(datGPThis),function(x) {strsplit(x,"_")[[1]][2]},USE.NAMES=F)
                                    id2=colnames(arrayData)
                                    jj=match(sapply(colnames(datGPThis),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F),altTypeUniq1)
                                    for (g1 in 1:nrow(arrayData)) {
                                        for (a1 in 1:length(altTypeUniq1)) {
                                            i1=g0[g1]
                                            j=which(jj==a1 & datGPThis[i1,]==1)
                                            j2=match(id1[j],id2)
                                            j21=j2[which(arrayData[g1,j2]==0)]
                                            j22=j2[which(arrayData[g1,j2]!=0)]
                                            arrayData[g1,j21]=a1
                                            if ("multiple"%in%altTypeUniq2[,1]) {
                                                arrayData[g1,j22]=10
                                            } else {
                                                arrayData[g1,j22]=10*arrayData[g1,j22]+a1
                                            }
                                        }
                                    }
                                }
                                #}
                            } else {
                                altTypeUniq1=sort(unique(clin1$altType))
                                arrayData=datGPThis[g0,grep(altTypeUniq1[1],colnames(datGPThis))]
                                colnames(arrayData)=sub(paste(altTypeUniq1[1],"_",sep=""),"",colnames(arrayData))
                                if (length(altTypeUniq1)==1) {
                                    
                                } else {
                                    for (a1 in 2:length(altTypeUniq1)) {
                                        arrayData[datGPThis[g0,grep(altTypeUniq1[a1],colnames(datGPThis))]==1]=a1
                                    }
                                }
                            }
                            #j=match(colnames(arrayData),clin3$id); j1=which(!is.na(j)); j2=j[j1]
                            j=match(colnames(arrayData),clin3$id[samId3]); j1=which(!is.na(j)); j2=samId3[j[j1]]
                            #annCol=data.frame(id=colnames(arrayData),disOntTerm2=clin3$disOntTerm2[j],stringsAsFactors=F)
                            arrayData=arrayData[,j1]
                            annCol=clin3[j2,which(!names(clin3)%in%colId2)]
                            if (!subsetFlag[1]%in%c("_ucsf500","_ucsf500Fmi")) {
                                if (exists("clin32")) {
                                    annCol=cbind(annCol,clin32[j2,which(!names(clin32)%in%names(annCol))])
                                    annCol$cellSize[which(annCol$disOntTerm2=="lungSmall")]="SC"
                                    annCol$cellSize[which(annCol$cellSize%in%c("MC","SC/LC","SC/MC"))]="other"
                                    annCol$cellSize[which(!annCol$cellSize%in%c("SC","LC","other"))]=NA
                                }
                            }
                            annRow=NULL
                            annRow=data.frame(id=rownames(arrayData),gene=rownames(arrayData),stringsAsFactors=F)
                            if (!is.null(geneFamilyList)) {
                                annRow$geneFamily=NA
                                i1=which(candGene$family2%in%geneFamilyList)
                                i=match(candGene$gene[i1],annRow$gene); i1=i1[!is.na(i)]; i2=i[!is.na(i)]
                                if (length(i1)!=0) {
                                    annRow$geneFamily[i2]=candGene$family2[i1]
                                }
                            }
                            
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
                            
                            if (clusterFlag[1]=="") {
                                if (subsetFlag[1]=="_allAssays") {
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
                                if (length(x)!=0) {
                                    varList=c(varList,x)
                                    varName=c(varName,paste(x," ",sep=""))
                                }
                            }
                            if (datTypeFlag[1]%in%c("_fracRead")) colHM=colHMCon else colHM=colHMCat
                        },
                        "_gene"={
                            tbl=read.table(paste(datadirG,"pv_",genesetFlag,geneset2Flag,subsetFlag[2],".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                            datGPName=rep("",nrow(tbl)*(nrow(tbl)-1)/2)
                            k=1
                            for (i1 in 1:(nrow(tbl)-1)) {
                                for (i2 in (i1+1):nrow(tbl)) {
                                    datGPName[k]=paste(tbl$gene[i1],tbl$gene[i2],sep="_")
                                    k=k+1
                                }
                            }
                            genePvMat=as.matrix(tbl[,-1])
                            x=genePvMat[lower.tri(genePvMat)]
                            tmp=rep(NA,length(x))
                            out=data.frame(id=datGPName,pv=x,qv=tmp,bh=tmp,by=tmp,holm=tmp,bonf=tmp,stringsAsFactors=F)
                            i=which(!is.na(out$pv))
                            out$qv=NA
                            res=try(qvalue(out$pv[i]))
                            if (class(res)=="qvalue") {
                                out$qv[i]=res$qvalues
                            }
                            geneQvMat=matrix(nrow=nrow(genePvMat),ncol=nrow(genePvMat))
                            k=1
                            for (i1 in 1:(nrow(genePvMat)-1)) {
                                for (i2 in (i1+1):nrow(genePvMat)) {
                                    geneQvMat[i1,i2]=geneQvMat[i2,i1]=out$qv[k]
                                    k=k+1
                                }
                            }
                            #tbl=read.table(paste(datadirG,ifelse(testFlag=="cosine","cor","cor"),"_",genesetFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                            #tbl=read.table(paste(datadirG,ifelse(testFlag=="cosine","cor","comut"),"_",genesetFlag,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                            tbl=read.table(paste(datadirG,ifelse(testFlag=="cosine","cor","comut"),"_",genesetFlag,geneset2Flag,subsetFlag[2],".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
                            arrayData=as.matrix(tbl[,-1])
                            rownames(arrayData)=tbl$gene
                            if (subsetFFlag!="") {
                                nm=toupper(sub("_","",subsetFFlag))
                                x=arrayData
                                x=x[which(rownames(x)==nm),which(colnames(x)!=nm)]
                                if (length(x)==0) {
                                    cat("No ",nm," gene!!!\n")
                                    next
                                } else {
                                    x=names(x)[which(x!=0)]
                                    if (length(x)==0) {
                                        cat("No genes co-mutated with ",nm,"!!!\n")
                                        next
                                    }
                                    j=match(c(nm,x),rownames(arrayData))
                                    arrayData=arrayData[j,j]
                                    genePvMat=genePvMat[j,j]
                                    geneQvMat=geneQvMat[j,j]
                                }
                            }
                            diag(arrayData)=NA
                            diag(genePvMat)=NA
                            diag(geneQvMat)=NA
                            #rownames(arrayData)=colnames(arrayData)
                            x2=apply(arrayData,2,function(x) mean(!is.na(x)))
                            j=which(x2>0)
                            arrayData=arrayData[j,j]
                            genePvMat=genePvMat[j,j]
                            geneQvMat=geneQvMat[j,j]
                            x=geneMutList[[genesetFlag]]; names(x)=sub("-",".",names(x))
                            if (subsetFFlag!="") {
                                x=geneMutList[[genesetFlag]]; names(x)=sub("-",".",names(x))
                                x=x[match(rownames(arrayData),names(x))]
                            }
                            j=match(names(x),rownames(arrayData))
                            arrayData=arrayData[j,j]
                            genePvMat=genePvMat[j,j]
                            geneQvMat=geneQvMat[j,j]
                            annRow=annCol=data.frame(id=colnames(arrayData),gene=colnames(arrayData),geneFamily=rep("",ncol(arrayData)),mutPerc=x[match(colnames(arrayData),names(x))],stringsAsFactors=F)
                            arrayData[lower.tri(arrayData,diag=T)]=NA
                            genePvMat[lower.tri(arrayData,diag=T)]=NA
                            geneQvMat[lower.tri(arrayData,diag=T)]=NA
                            
                            jj=grep(altTypeUniq1[1],colnames(datGPThis))
                            x=matrix(0,nrow=length(g0),ncol=length(jj))
                            colnames(x)=sub(paste(altTypeUniq1[1],"_",sep=""),"",colnames(datGPThis)[jj])
                            rownames(x)=rownames(datGPThis)[g0]
                            j=match(colnames(x),clin3$id[samId3]); j1=which(!is.na(j)); j2=samId3[j[j1]]
                            x=x[,j1]
                            annRow$numSample=annCol$numSample=apply(x,1,function(x) {sum(!is.na(x))})
                            if (mutFlag=="_numMut") {
                                annRow$mutPerc=annCol$mutPerc=round(annRow$mutPerc*annRow$numSample[1])
                            } else {
                                annRow$mutPerc=annCol$mutPerc==100*round(annRow$mutPerc,2)
                                if (testFlag=="fisher") {
                                    arrayData=arrayData/annRow$numSample[1]
                                }
                            }
                            
                            if (nrow(arrayData)>2) {
                                i=(nrow(arrayData)-1):1
                                arrayData=arrayData[i,]
                                genePvMat=genePvMat[i,]
                                geneQvMat=geneQvMat[i,]
                                annRow=annCol[i,]
                                
                                #if (F) {
                                #j=ncol(arrayData):1
                                j=2:ncol(arrayData)
                                arrayData=arrayData[,j]
                                genePvMat=genePvMat[,j]
                                geneQvMat=geneQvMat[,j]
                                annCol=annCol[j,]
                                #}
                            } else {
                                getClustInfoFlag=F
                                fNameThis=paste(subDir,"clusterInfoFeature",fNameOut,".txt",sep="")
                                fNameThis=paste(subDir,"clusterInfoFeature",fNameOut,".txt",sep="")
                                x=paste("Gene cluster information for heatmap of samples from\n",heading,sep="")
                                #x=paste(x,"\nGene co-mutation: ",100*round(arrayData[1,2],2),"%",sep="")
                                write.table(x, fNameThis, sep="\t", col.names=F, row.names=F, quote=F)
                                tbl=cbind(annRow)
                                write.table(tbl, fNameThis, sep="\t", col.names=T, row.names=F, quote=F,append=T)
                                i=nrow(arrayData)
                                rownames(arrayData)[i]=annRow$id[i]=""; annRow[i,c("gene","mutPerc")]=NA
                                i=1
                                colnames(arrayData)[i]=annCol$id[i]=""; annCol[i,c("gene","mutPerc")]=NA
                                
                                #next
                            }
                            
                            #varList="mutPerc"; varName=paste(varList," ",sep="")
                            varList="mutPerc"; varName=""
                            varFList="mutPerc"; varFName=paste(varFList," ",sep="")
                            #annColAll=annCol; annRowAll=annRow
                            annRowAll=annRow
                            
                            colHM=colHMCon
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
                            ##j=which(apply(arrayData,2,function(x) {mean(x!=0,na.rm=T)})!=0)
                            #j=apply(arrayData,2,function(x) {y=mean(x!=0,na.rm=T); !is.na(y) && y!=0})
                            i=which(annRow$gene=="SETD2")
                            j=which(arrayData[i,]!=0)
                            #i=which(apply(arrayData[,j],1,function(x) {mean(x!=0,na.rm=T)})!=0)
                            i=apply(arrayData[,j],1,function(x) {y=mean(x!=0,na.rm=T); !is.na(y) && y!=0})
                        }
                        )
                        if (length(i)==0 | length(j)==0) {
                            cat("No genes or samples !!!\n")
                            next
                        }
                        arrayData=arrayData[i,j]
                        annRow=annRow[i,]
                        annCol=annCol[j,]
                        if (!is.null(genePvMat)) {
                            genePvMat=genePvMat[i,j]
                        }
                        if (!is.null(geneQvMat)) {
                            geneQvMat=geneQvMat[i,j]
                        }
                        
                        #if (nrow(arrayData)<3) {
                        if (clusterFlag[1]%in%c("_supervised","_comutated") | (clusterFlag[1]=="" & genesetFlag%in%c("_setd2","_notch","_mmr","_brca","_allGenes"))) {
                            if (subset3Flag[1]=="_no0AltGene") {
                                heading=paste(heading,", genes with >0 alt",sep="")
                                heading=sub("genes, genes","genes",heading)
                                #i=apply(arrayData,1,function(x) {mean(x!=0,na.rm=T)!=0})
                                i=apply(arrayData,1,function(x) {y=mean(x!=0,na.rm=T); !is.na(y) && y!=0})
                                arrayData=arrayData[i,]
                                if (!is.null(annRow)) annRow=annRow[i,]
                                if (!is.null(genePvMat)) {
                                    genePvMat=genePvMat[i,]
                                }
                                if (!is.null(geneQvMat)) {
                                    geneQvMat=geneQvMat[i,]
                                }
                            }
                            #if (clusterFlag[1]%in%c("_comutated") & subset2Flag=="_lungSmall") {
                            if (geneset2Flag!="") {
                                if (F) {
                                    heading=paste(heading,"\npatients with > 0 alt in >= 50% genes",sep="")
                                    i=which(annRow$id%in%candGene1)
                                    j=apply(arrayData[i,],2,function(x) {y=mean(x!=0,na.rm=T); !is.na(y) && y!=0})
                                    i=which(!annRow$id%in%candGene1)
                                    j=j | apply(arrayData[i,],2,function(x) {y=mean(x!=0,na.rm=T); !is.na(y) && y>=.5})
                                }
                                heading=paste(heading,", pats with >0 alt",sep="")
                                ii=which(annRow$id%in%candGene1)
                                ii=which(annRow$id%in%sub("_","",geneset2Flag))
                                if (length(ii)==0) {
                                    cat("(length(ii)==0)\n")
                                    next
                                }
                                j=rep(F,nrow(annCol))
                                for (i in ii) {
                                    x1=rep(T,nrow(annCol))
                                    x2=rep(T,nrow(annCol))
                                    x1=!is.na(arrayData[i,]); x1[x1][arrayData[i,x1]==0]=F
                                    if (nrow(arrayData)<=2) {
                                        x2=sapply(arrayData[-i,],function(x) {!is.na(x) && x!=0},USE.NAMES=F)
                                    } else {
                                        x2=apply(arrayData[-i,],2,function(x) {y=mean(x!=0,na.rm=T); !is.na(y) && y!=0})
                                    }
                                    j[x1 & x2]=T
                                }
                            } else {
                                j=rep(T,nrow(annCol))
                                #if (clusterFlag[1]=="") {
                                if (clusterFlag[1]%in%c("","_comutated") | (length(clusterFlag)==1 & clusterFlag[1]%in%c("_supervised") & length(grep("topGene",genesetFlag))==1)) {
                                    #heading=paste(heading,", pats with >0 alt",sep="")
                                }
                                #if (!clusterFlag[1]%in%c("_supervised") | (length(clusterFlag)==1 & clusterFlag[1]%in%c("_supervised"))) {
                                #if (!clusterFlag[1]%in%c("_supervised")) {
                                if (!clusterFlag[1]%in%c("_supervised") & !subsetFlag[1]%in%c("_ucsf500","_ucsf500Fmi")) {
                                    #j=apply(arrayData,2,function(x) {mean(x!=0,na.rm=T)!=0})
                                    j=apply(arrayData,2,function(x) {y=mean(x!=0,na.rm=T); !is.na(y) && y!=0})
                                }
                            }
                            if (sum(j)==0) {
                                cat("(sum(j)==0)\n")
                                next
                            }
                            if (!clusterFlag[1]%in%c("_comutated")) {
                                if (sum(j)<5) j=rep(T,ncol(arrayData))
                            }
                            if (sum(j)==1) {
                                arrayData=matrix(arrayData[,j],ncol=1)
                                #annCol=matrix(annCol[j,],nrow=1)
                                if (!is.null(genePvMat)) {
                                    genePvMat=matrix(genePvMat[,j],ncol=1)
                                }
                                if (!is.null(geneQvMat)) {
                                    geneQvMat=matrix(geneQvMat[,j],ncol=1)
                                }
                            } else {
                                arrayData=arrayData[,j]
                                #annCol=annCol[j,]
                                if (!is.null(genePvMat)) {
                                    genePvMat=genePvMat[,j]
                                }
                                if (!is.null(geneQvMat)) {
                                    geneQvMat=geneQvMat[,j]
                                }
                            }
                            annCol=annCol[j,]
                            if (subset3Flag[1]=="_no0AltGene") {
                                i=apply(arrayData,1,function(x) {y=mean(x!=0,na.rm=T); !is.na(y) && y!=0})
                                if (ncol(arrayData)==1) {
                                    arrayData=matrix(arrayData[i,],ncol=1)
                                    if (!is.null(genePvMat)) {
                                        genePvMat=matrix(genePvMat[i,],ncol=1)
                                    }
                                    if (!is.null(geneQvMat)) {
                                        geneQvMat=matrix(geneQvMat[i,],ncol=1)
                                    }
                                } else {
                                    if (sum(i)==1) {
                                        arrayData=matrix(arrayData[i,],nrow=1)
                                        if (!is.null(genePvMat)) {
                                            genePvMat=matrix(genePvMat[i,],nrow=1)
                                        }
                                        if (!is.null(geneQvMat)) {
                                            geneQvMat=matrix(geneQvMat[i,],nrow=1)
                                        }
                                    } else {
                                        arrayData=arrayData[i,]
                                        if (!is.null(genePvMat)) {
                                            genePvMat=genePvMat[i,]
                                        }
                                        if (!is.null(geneQvMat)) {
                                            geneQvMat=geneQvMat[i,]
                                        }
                                    }
                                }
                                
                                if (!is.null(annRow)) annRow=annRow[i,]
                            }
                        }
                        if (clusterFlag[1]%in%c("_comutated") & any(dim(arrayData)==1)) {
                            cat("(clusterFlag[1]%in%c(_comutated) & any(dim(arrayData)==1))\n")
                            next
                        }
                        
                        if (!geneFamilyFlag & length(grep("_actionableGene",genesetFlag))==1) {
                            annRow$geneFamily=""
                            i=match(toupper(candGeneThis$gene),toupper(annRow$gene)); i2=which(!is.na(i)); i1=i[i2]
                            annRow$geneFamily[i1]=candGeneThis$family2[i2]
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
                        if (datTypeFlag[1]%in%c("_fracRead")) {
                            arrayData2[is.na(arrayData)]=-99
                        } else {
                            if (altTypeFlag=="_anyAlt") {
                                arrayData2[arrayData2!=0]=1
                            }
                            if (altTypeFlag=="_anyAlt1minus1") {
                                arrayData2[arrayData2!=0]=1; arrayData2[arrayData2==0]=-1
                            }
                            arrayData[arrayData==0]=NA
                        }
                        
                        if (!is.na(nClust[1]) && nrow(arrayData)>=100) {
                            nClust[1]=10
                        }
                        if (!is.na(nClust[2]) && ncol(arrayData)>=100) {
                            #nClust[2]=10
                        }
                        if (!is.na(nClust[1]) && nrow(arrayData)<(nClust[1]+1)) {
                            nClust[1]=2
                        }
                        if (!is.na(nClust[2]) && ncol(arrayData)<(nClust[2]+1)) {
                            nClust[2]=2
                        }
                        if (geneBar=="" || (!is.na(nClust[1]) && nrow(arrayData)<10)) {
                            nClust[1]=NA
                        }
                        if (sampleBar=="" || (!is.na(nClust[2]) && ncol(arrayData)<10)) {
                            nClust[2]=NA
                        }
                        
                        annRowAll=annRow
                        #annColAll=annCol
                        varListAll=varList
                        varNameAll=varName
                        varFListAll=varFList
                        varFNameAll=varFName
                        colHMAll=colHM
                        
                        ## -------------------
                        #if (clusterFlag[1]=="_supervised") {
                        if ((length(clusterFlag)==1 & clusterFlag[1]%in%c("_supervised")) | (clusterFlag[1]%in%c("_comutated") & sampleBar=="")) {
                            var2List=varList
                            varId=NULL
                            varId=which(varList=="cellSize")
                            
                            annCol$grp=""
                            switch(ordFlag,
                            "_tmbLevel"={
                                #annCol$grp=annCol$tmbLevel
                                #annCol$grp[is.na(annCol$grp)]=""
                                annCol$grp=as.integer(as.factor(annCol$tmbLevel))
                                annCol$grp[is.na(annCol$grp)]=99
                                annCol$grp=as.character(annCol$grp)
                            },
                            "_rb1P53"={
                                annCol$grp=arrayData[which(annRow$id=="RB1"),]
                                annCol$grp=paste(annCol$grp,arrayData[which(annRow$id=="TP53"),])
                            },
                            "_swiSnfEtc"={
                                grpUniq=unlist(clusterGene)
                                annCol$grp=arrayData[which(annRow$id==grpUniq[1]),]
                                for (grpId in 2:length(grpUniq)) {
                                    annCol$grp=paste(annCol$grp,arrayData[which(annRow$id==grpUniq[grpId]),])
                                }
                            },
                            "_tmbSwiSnfEtc"={
                                grpUniq=unlist(clusterGene)
                                annCol$grp=annCol$tmbLevel
                                for (grpId in 1:length(grpUniq)) {
                                    annCol$grp=paste(annCol$grp,arrayData[which(annRow$id==grpUniq[grpId]),])
                                }
                            },
                            "_diffEBki67EB"={
                                annCol$grp=10*annCol$ki67EB2; annCol$grp[which(annCol$ki67EB=="NR")]=2*max(annCol$grp,na.rm=T)
                                annCol$grp=paste(annCol$diffEB,formatC(annCol$grp,width=nchar(as.character(max(annCol$grp,na.rm=T))),flag="0"))
                            },
                            "_priSiteEBki67EB"={
                                annCol$grp=10*annCol$ki67EB2; annCol$grp[which(annCol$ki67EB=="NR")]=2*max(annCol$grp,na.rm=T)
                                annCol$grp=paste(annCol$primarySiteEB2,formatC(annCol$grp,width=nchar(as.character(max(annCol$grp,na.rm=T))),flag="0"))
                            },
                            "_diffEBpriSiteEBki67EB"={
                                annCol$grp=10*annCol$ki67EB2; annCol$grp[which(annCol$ki67EB=="NR")]=2*max(annCol$grp,na.rm=T)
                                annCol$grp=paste(annCol$diffEB,annCol$primarySiteEB2,formatC(annCol$grp,width=nchar(as.character(max(annCol$grp,na.rm=T))),flag="0"))
                            }
                            )
                            
                            annColAll$grp=NA
                            annColAll$grp[match(annCol$id,annColAll$id)]=annCol$grp
                            var2List="grp"
                            varId=which(var2List=="grp")
                            if (length(varId)==1) {
                                grp=annCol[,var2List[varId]]
                            } else {
                                varId=1
                                if (any(c("disease","disOntTerm","disOntTerm2","disOntTerm3","primarySiteEB2")%in%varList[varId])) {
                                    grp=annCol[,varList[varId]]
                                } else {
                                    grp=paste(annCol[,varList[varId]],annCol$disOntTerm2)
                                }
                            }
                            grpUniq=sort(unique(grp))
                            j=c()
                            for (grpId in 1:length(grpUniq)) {
                                jj=which(grp==grpUniq[grpId])
                                if (length(jj)==1) {
                                    dat=matrix(arrayData2[,jj],ncol=1)
                                } else {
                                    dat=arrayData2[,jj]
                                }
                                if (!is.null(clusterGene)) {
                                    grpC=nm=c()
                                    for (k in 1:length(clusterGene)) {
                                        p=1
                                        kk=which(rownames(dat)==clusterGene[[k]][p])
                                        if (length(kk)!=0) {
                                            out=as.integer(dat[kk,]!=0)
                                            nm1=clusterGene[[k]][p]
                                            if (length(clusterGene[[k]])>1) {
                                                for (p in 2:length(clusterGene[[k]])) {
                                                    kk=which(rownames(dat)==clusterGene[[k]][p])
                                                    if (length(kk)!=0) {
                                                        out=as.integer(out!=0 & dat[kk,]!=0)
                                                        nm1=paste(nm1,clusterGene[[k]][p],sep="_")
                                                    }
                                                }
                                            }
                                            grpC=cbind(grpC,out)
                                            nm=c(nm,nm1)
                                        }
                                    }
                                    if (length(grpC)!=0) {
                                        nm=gsub("-",".",nm)
                                        colnames(grpC)=nm
                                        dat=dat[which(rownames(dat)%in%unlist(clusterGene)),]
                                        
                                        k=match(colnames(grpC),names(annCol)); k=colnames(grpC)[!is.na(k)]
                                        if (length(k)==0) {
                                            tmp=matrix(nrow=nrow(annCol),ncol=ncol(grpC),dimnames=list(rownames(annCol),nm))
                                            annCol=cbind(annCol,tmp)
                                            if (F) {
                                                varList=c(varList,nm)
                                                varName=c(varName,paste(tolower(nm)," ",sep=""))
                                                if (!"gene"%in%varListAll) {
                                                    varListAll=c(varListAll,"gene")
                                                    varNameAll=c(varNameAll,paste("gene"," ",sep=""))
                                                }
                                            }
                                        }
                                        for (k in colnames(grpC)) {
                                            annCol[jj,k]=grpC[,k]
                                        }
                                        #}
                                        k=match(colnames(grpC),names(annColAll)); k=colnames(grpC)[!is.na(k)]
                                        if (length(k)==0) {
                                            tmp=matrix(nrow=nrow(annColAll),ncol=ncol(grpC),dimnames=list(rownames(annColAll),nm))
                                            annColAll=cbind(annColAll,tmp)
                                        }
                                        k=match(colnames(grpC),names(annColAll)); k=colnames(grpC)[!is.na(k)]
                                        tmp=annColAll[,k]
                                        tmp[match(annCol$id[jj],annColAll$id),]=grpC[,k]
                                        for (k in colnames(tmp)) {
                                            annColAll[,k]=tmp[,k]
                                        }
                                    }
                                }
                                if (sampleBar=="cluster" & nrow(dat)>1 & ncol(dat)>1) {
                                    switch(distMethod,
                                    "pearson"={
                                        distMat=as.dist(1 - cor(t(dat),method=distMethod,use="complete.obs"))
                                    },
                                    "spearman"={
                                        distMat=as.dist(1 - cor(t(dat),method=distMethod,use="complete.obs"))
                                    },
                                    "kendall"={
                                        if (getCorFlag) {
                                            corMatSam2=cor(t(dat),method=distMethod,use="complete.obs")
                                            save(corMatSam2,file=paste("corMatSam",fNameOut,".RData",sep=""))
                                        } else {
                                            load(file=paste(datadir2,"corMatSam",fNameOut,".RData",sep=""))
                                        }
                                        distMat=as.dist(1 - corMatSam2)
                                    },
                                    "euclidean"={
                                        distMat=dist(t(dat), method=distMethod)
                                    },
                                    "cosine"={
                                        if (any(apply(dat,2,sum,na.rm=T)==0)) dat=dat+1
                                        distMat=getCosineDist(t(dat))
                                    }
                                    )
                                    clustC=hclust(distMat, method=linkMethod)
                                    j=c(j,jj[clustC$order])
                                } else {
                                    j=c(j,jj)
                                }
                            }
                            if (is.na(nClust[2])) {
                                x=F
                            } else {
                                x=length(clustC$labels)==nrow(annCol)
                                if (x) x=all(clustC$labels==annCol$id)
                            }
                            if (x) {
                                j=1:ncol(arrayData)
                            } else {
                                clustC=NA
                                nClust[2]=NA
                            }
                            i=nrow(annRow):1
                            arrayData=arrayData[i,j]
                            arrayData2=arrayData2[i,j]
                            annRow=annRow[i,]
                            annCol=annCol[j,]
                            if (!is.null(genePvMat)) {
                                genePvMat=genePvMat[i,j]
                            }
                            if (!is.null(geneQvMat)) {
                                geneQvMat=geneQvMat[i,j]
                            }
                        } else {
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
                                        load(file=paste(datadir2,"corMatSam",fNameOut,".RData",sep=""))
                                    }
                                    distMat=as.dist(1 - corMatSam2)
                                },
                                "euclidean"={
                                    distMat=dist(t(arrayData2), method=distMethod)
                                },
                                "cosine"={
                                    dat=arrayData2
                                    if (any(apply(dat,2,sum,na.rm=T)==0)) dat=dat+1
                                    distMat=getCosineDist(t(dat))
                                }
                                )
                                clustC=hclust(distMat, method=linkMethod)
                            } else {
                                clustC=NA
                                nClust[2]=NA
                            }
                        }
                        
                        
                        
                        nClust[1]=NA
                        #nClust[2]=NA
                        
                        
                        
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
                                    load(file=paste(datadir2,"corMat",fNameOut,".RData",sep=""))
                                }
                                distMat=as.dist(1 - corMat2)
                            },
                            "euclidean"={
                                distMat=dist(arrayData2, method=distMethod)
                            },
                            "cosine"={
                                dat=arrayData2
                                if (any(apply(dat,1,sum,na.rm=T)==0)) dat=dat+1
                                distMat=getCosineDist(dat)
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
                        #thres=20
                        thres=100
                        thres=150
                        if (nrow(arrayData)<=thres) nameRow=annRow$id else nameRow=rep("",nrow(arrayData))
                        switch(genesetFlag,
                        "_setd2"={
                            i=which(annRow$gene=="SETD2")
                            nameRow[i]=paste(annRow$gene[i]," ",paste(rep("*",10),collapse=""),sep="")
                        }
                        )
                        #if (F) {
                        if (clusterFlag[1]%in%c("_comutatedGenes")) {
                            i=which(nameRow!="")
                            #nameRow[i]=paste(annRow$id," ",annRow$numSample,sep="")[i]
                            nameRow[i]=paste(annRow$id," ",annRow$mutPerc,sep="")[i]
                        }
                        #}
                        colRow=NULL
                        if (clusterFlag[1]=="_supervised") {
                            colRow=matrix(nrow=3,ncol=nrow(annRow))
                            rownames(colRow)=rep(" ",nrow(colRow))
                        }
                        if (is.null(varFList)) {
                            colRow=NULL
                        } else {
                            colRow=matrix(nrow=length(varFList),ncol=nrow(annRow))
                            for (varId in 1:length(varFList)) {
                                if (varFList[varId]%in%c("sd")) {
                                    j=match(annRow$id,annRowAll$id)
                                    x=round(100*annRowAll[,varFList[varId]])+1
                                    lim=range(x,na.rm=T)
                                    x[x<lim[1]]=lim[1]; x[x>lim[2]]=lim[2]
                                    grpUniq=lim[1]:lim[2]
                                    colRowUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                    colRow[varId,]=colRowUniq[x[j]]
                                } else if (varFList[varId]%in%c("mutPerc")) {
                                    j=match(annRow$id,annRowAll$id)
                                    x=round(annRowAll[,varList[varId]])+1
                                    lim=range(x,na.rm=T)
                                    if (varList[varId]==c("mutPerc")) lim=limPerc
                                    grpUniq=lim[1]:lim[2]
                                    colRowUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                    colRow[varId,]=colRowUniq[x[j]]
                                } else {
                                    x=as.character(annRowAll[,varFList[varId]])
                                    x[x==""]=NA; x=as.integer(as.factor(x))
                                    grpUniq=sort(unique(x))
                                    x=x[match(annRow$id,annRowAll$id)]
                                    if (varFList[varId]%in%c("geneFamily")) {
                                        colRow[varId,]=colListGene[x]
                                    } else {
                                        if (length(grpUniq)<=length(colList2)) {
                                            colRow[varId,]=colList2[x]
                                        } else if (length(grpUniq)<=length(colList)) {
                                            colRow[varId,]=colList[x]
                                        } else {
                                            colRow[varId,]=rainbow(length(grpUniq))[x]
                                        }
                                    }
                                }
                            }
                            rownames(colRow)=varFName
                        }
                        
                        ## -------------------
                        #thres=6
                        thres=31
                        thres=51
                        if (datTypeFlag[2]%in%c("_prop","_gene") | ncol(arrayData)<thres) nameCol=annCol$id else nameCol=rep("",ncol(arrayData))
                        #if (F) {
                        if (clusterFlag[1]%in%c("_comutatedGenes")) {
                            i=which(nameCol!="")
                            #nameCol[i]=paste(annCol$id," ",annCol$numSample,sep="")[i]
                            nameCol[i]=paste(annCol$id," ",annCol$mutPerc,sep="")[i]
                        }
                        #}
                        if (is.null(varList)) {
                            colCol=NULL
                        } else {
                            colCol=matrix(nrow=length(varList),ncol=nrow(annCol))
                            for (varId in 1:length(varList)) {
                                if (varList[varId]%in%c("sd","ki67EB2")) {
                                    j=match(annCol$id,annColAll$id)
                                    if (varList[varId]==c("ki67EB2")) {
                                        x=round(annColAll[,varList[varId]])+1
                                        lim=c(0,200)+1
                                        lim=c(0,100)+1
                                    } else {
                                        x=round(100*annColAll[,varList[varId]])+1
                                        lim=range(x,na.rm=T)
                                        if (varList[varId]==c("sd")) lim=limSdSam
                                    }
                                    x[x<lim[1]]=lim[1]; x[x>lim[2]]=lim[2]
                                    grpUniq=lim[1]:lim[2]
                                    if (varList[varId]==c("ki67EB2")) {
                                        colColUniq <- maPalette(low="indianred1",high="darkred",k=100)
                                        #maColorBar(1:length(pal), col=colColUniq, horizontal=FALSE, k=length(pal))
                                    } else {
                                        #colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                        #colColUniq=rev(gray(1:(length(grpUniq))/length(grpUniq)))
                                        colColUniq=rev(rgb((0:length(grpUniq))/length(grpUniq), green = 0, blue = 0, names = paste("red", 0:length(grpUniq), sep = ".")))
                                        #x=colors(); x=x[grep("red",x)]
                                    }
                                    colCol[varId,]=colColUniq[x[j]]
                                } else if (varList[varId]%in%c("mutPerc") | length(grep("dist2",varList[varId]))==1) {
                                    j=match(annCol$id,annColAll$id)
                                    x=round(annColAll[,varList[varId]])+1
                                    lim=range(x,na.rm=T)
                                    if (varList[varId]==c("mutPerc")) lim=limPerc
                                    grpUniq=lim[1]:lim[2]
                                    colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                    colCol[varId,]=colColUniq[x[j]]
                                } else {
                                    x=as.character(annColAll[,varList[varId]])
                                    x[x==""]=NA; x=as.integer(as.factor(x))
                                    grpUniq=sort(unique(x))
                                    x=x[match(annCol$id,annColAll$id)]
                                    if (varList[varId]%in%c("disOntTerm2","primarySiteEB2")) {
                                        colCol[varId,]=colListD[x]
                                    } else if (varList[varId]%in%c("cellSize","cellSizeEB")) {
                                        colCol[varId,]=colListC[x]
                                    } else if (varList[varId]%in%c("diffEB")) {
                                        colCol[varId,]=colListDi[x]
                                    } else if (varList[varId]%in%c("grade")) {
                                        colCol[varId,]=colListGr[x]
                                    } else if (varList[varId]%in%c("tmbLevel")) {
                                        colCol[varId,]=colListT[x]
                                    } else {
                                        if (length(grpUniq)<=length(colList2)) {
                                            colCol[varId,]=colList2[x]
                                        } else if (length(grpUniq)<=length(colList)) {
                                            colCol[varId,]=colList[x]
                                        } else {
                                            colCol[varId,]=rainbow(length(grpUniq))[x]
                                        }
                                    }
                                }
                            }
                            rownames(colCol)=varName
                        }
                        
                        ## -------------------
                        if (datTypeFlag[1]%in%c("_fracRead") | datTypeFlag[2]%in%c("_gene")) {
                            print("summary(c(arrayData))")
                            print(summary(c(arrayData)))
                            print("quantile(abs(c(arrayData)),probs=seq(0,1,by=.1),na.rm=T)")
                            print(quantile(abs(c(arrayData)),probs=seq(0,1,by=.1),na.rm=T))
                        } else {
                            print("table(c(arrayData2))")
                            print(table(c(arrayData2)))
                        }
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
                            if (datTypeFlag[2]=="_prop") margins=c(20,1) else margins=c(1,1)
                            margins[2]=5
                            png(paste(subDir,"heatmap",fNameOut,".png",sep=""),width=480*2,height=480*2)
                        } else {
                            margins=c(12,5)
                            margins=c(12,7)
                            margins=c(2,7)
                            #if (clusterFlag[1]%in%c("_comutated")) {
                            if (clusterFlag[1]%in%c("_supervised","_comutated")) {
                                margins=c(4,7)
                            } else if (clusterFlag[1]%in%c("_comutatedGenes")) {
                                margins=c(7,7)
                            } else {
                                if (nrow(arrayData)>50) {
                                    margins=c(2,7)
                                } else {
                                    margins=c(round(2+(50-2*nrow(arrayData))),7)
                                }
                            }
                            pdf(paste(subDir,"heatmap",fNameOut,".pdf",sep=""),width=7,height=7)
                            par(cex.main=0.7)
                        }
                        x=arrayData
                        if (!datTypeFlag[1]%in%c("_fracRead")) {
                            for (j in 1:ncol(x)) {
                                for (k in 1:nrow(altTypeUniq2)) {
                                    x[which(arrayData[,j]==as.integer(altTypeUniq2[k,2])),j]=k
                                }
                                #x[which(arrayData2[,j]==4),j]=7
                                #x[which(arrayData2[,j]>4),j]=NA
                            }
                        }
                        lim=60
                        lim=150
                        cexThis=c(2,2,0.5)
                        if (clusterFlag[1]%in%c("_comutated","_comutatedGenes")) lim=150
                        if (nrow(arrayData)>lim) {
                            nameRow=rep("",length(nameRow))
                            cexThis[1]=0.75
                        } else if (nrow(arrayData)>30) {
                            i=seq(1,length(nameRow),by=2)
                            nameRow[i]=paste(nameRow[i],"                  ",sep="")
                            cexThis[1]=0.75
                        } else {
                            cexThis[1]=1.5
                            cexThis[3]=1.5
                            cexThis[1]=1
                            cexThis[3]=1
                        }
                        if (ncol(arrayData)>lim) {
                            nameCol=rep("",length(nameCol))
                            cexThis[2]=0.75
                        } else if (ncol(arrayData)>30) {
                            i=seq(1,length(nameCol),by=2)
                            nameCol[i]=paste("                  ",nameCol[i],sep="")
                            cexThis[2]=0.75
                        } else {
                            cexThis[2]=1.5
                            cexThis[3]=1.5
                            cexThis[2]=1
                            cexThis[3]=1
                        }
                        nameCol=rep("",length(nameCol))
                        if ("gene"%in%varListAll) cexThis[2]=1
                        #if (length(colHM[[1]])>1) {
                        if (is.null(colHM[[2]])) {
                            y=unique(c(x)); y=sort(y[!is.na(y)])
                            #colHM[[1]]=colHM[[1]][which((1:nrow(altTypeUniq2))%in%as.character(y))]
                            colHM[[1]]=colHM[[1]][1:max(y)]
                            #colHM[[1]]=colHM[[1]][y]
                            if (length(colHM[[1]])==1) colHM[[1]]=rep(colHM[[1]],2)
                            limit1[2]=max(c(x),na.rm=T)
                        }
                        if (lineFlag) {
                            lim=100
                            lim=150
                            if (max(dim(arrayData))>lim) {
                                lineList=list(row=NULL,col=NULL,color=NULL)
                            } else {
                                lineList=list(row=seq(0,(nrow(x)))+0.5,col=seq(0,(ncol(x)))+0.5,color="grey")
                            }
                        } else {
                            lineList=list(row=NULL,col=NULL,color=NULL)
                        }
                        margins[which(margins<0)]=0
                        #hcc=heatmap3(x=x, Rowv=clustR, Colv=clustC, distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=colCol, RowSideColors=colRow, labCol=nameCol, labRow=nameRow, ncr=nClust[1], ncc=nClust[2], scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit1, cexCol=cexThis[2], cexRow=cexThis[1], high=colHM[[1]], low=colHM[[2]], mid=colHM[[3]],lineRow=lineList$row, lineCol=lineList$col, lineColor=lineList$color, addText=ifelse(clusterFlag[1]=="_comutatedGenes",apply(100*round(x,2),c(1,2),as.character),NULL))
                        if (clusterFlag[1]=="_comutatedGenes") {
                            if (mutFlag=="_numMut") {
                                addTextThis=apply(x,c(1,2),as.character)
                                if (testFlag=="fisher") {
                                    x=x/annRow$numSample[1]
                                }
                            } else {
                                addTextThis=apply(100*round(x,2),c(1,2),as.character)
                            }
                            for (i1 in 1:nrow(x)) {
                                for (i2 in 1:nrow(x)) {
                                    if (!is.na(geneQvMat[i1,i2]) && geneQvMat[i1,i2]<.05) {
                                        addTextThis[i1,i2]=ifelse(is.na(addTextThis[i1,i2]),"**",paste(addTextThis[i1,i2],"**"))
                                    } else if (!is.na(genePvMat[i1,i2]) && genePvMat[i1,i2]<.05) {
                                        addTextThis[i1,i2]=ifelse(is.na(addTextThis[i1,i2]),"*",paste(addTextThis[i1,i2],"*"))
                                    }
                                }
                            }
                        } else {
                            addTextThis=NULL
                        }
                        hcc=heatmap3(x=x, Rowv=clustR, Colv=clustC, distfun=distMethod, hclustfun=hclust, symm=F, ColSideColors=colCol, RowSideColors=colRow, labCol=nameCol, labRow=nameRow, ncr=nClust[1], ncc=nClust[2], scale="none", na.rm=F, margins=margins, main=main, xlab=NULL, ylab=NULL, zlm=limit1, cexCol=cexThis[2], cexRow=cexThis[1], high=colHM[[1]], low=colHM[[2]], mid=colHM[[3]],lineRow=lineList$row, lineCol=lineList$col, lineColor=lineList$color, addText=addTextThis, cexText=cexThis[3])
                        dev.off()
                        
                        
                        ## -------------------
                        if (datTypeFlag[2]%in%c("_gene") & getClustInfoFlag) {
                            #if (!is.na(nClust[1])) {
                            fNameThis=paste(subDir,"clusterInfoFeature",fNameOut,".txt",sep="")
                            x=paste("Gene cluster information for heatmap of samples from\n",heading,"\nGenes are ordered as in heatmap from top to bottom. See column 'order'",sep="")
                            write.table(x, fNameThis, sep="\t", col.names=F, row.names=F, quote=F)
                            nClustThis=1
                            tbl=rbind(annCol,annRow)
                            tbl=tbl[nrow(annCol):nrow(tbl),]
                            tbl=tbl[nrow(tbl):1,]
                            tbl=cbind(tbl,order=1:nrow(tbl))
                            write.table(tbl, fNameThis, sep="\t", col.names=T, row.names=F, quote=F,append=T)
                        }
                        
                        if (datTypeFlag[2]%in%c("_patient")) {
                            #if (!is.na(nClust[1])) {
                            fNameThis=paste(subDir,"clusterInfoFeature",fNameOut,".txt",sep="")
                            x=paste("Gene cluster information for heatmap of samples from\n",heading,"\nGenes are ordered as in heatmap from bottom to top. See columns 'clustId' and 'order'",sep="")
                            write.table(x, fNameThis, sep="\t", col.names=F, row.names=F, quote=F)
                            if (F) {
                                pdf(paste(subDir,"clusterFeatures",fNameOut,".pdf",sep=""))
                                plot(clustR,main=paste("Feature clusters with ",nClust[1]," main clusters marked in red",sep=""),xlab="",sub="",ylab=NULL,axes=F, cex=.2); rect.hclust(clustR,k=nClust[1])
                                dev.off()
                            }
                            if (is.na(nClust[1])) {
                                nClustThis=1
                                tbl=cbind(annRow,order=1:nrow(annRow))
                            } else {
                                nClustThis=nClust[1]
                                clustId=cutree(clustR,k=nClustThis)[clustR$order]
                                k1=which(!duplicated(clustId))
                                for (k in 1:length(k1)) {
                                    clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                                }
                                tbl=cbind(annRow[clustR$order,],clustId,order=1:nrow(annRow))
                            }
                            write.table(tbl, fNameThis, sep="\t", col.names=T, row.names=F, quote=F,append=T)
                        }
                        
                        if (datTypeFlag[2]%in%c("_patient")) {
                            #if (!is.na(nClust[2])) {
                            colId=names(annCol)[which(!names(annCol)%in%c("noteEB","ki67EB2","grp"))]
                            colId=cbind(colId,colId)
                            if (F) {
                                if (subsetFlag[1]%in%c("_ucsf500","_ucsf500Fmi")) {
                                    colId[match(c("disOntTerm2"),colId[,1]),2]=c("primarySite")
                                }
                            }
                            fNameThis=paste(subDir,"clusterInfoSample",fNameOut,".txt",sep="")
                            x=paste("Sample cluster information for heatmap of samples from\n",heading,"\nSamples are ordered as in heatmap from left to right. See columns 'clustId' and 'order'",sep="")
                            write.table(x, fNameThis, sep="\t", col.names=F, row.names=F, quote=F)
                            #tbl=annCol
                            tbl=annCol[,colId[,1]]; names(tbl)=colId[,2]
                            #write.table("Global",file=fNameThis,append=F,col.names=F,row.names=F,sep="\t",quote=F)
                            #write.table(heading,file=fNameThis,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                            write.table("",file=fNameThis,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                            write.table(paste("No. of patients",nrow(tbl),sep="\t"),file=fNameThis,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                            write.table("",file=fNameThis,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                            if ("age"%in%names(tbl)) {
                                tbl$ageBi=as.integer(tbl$age<=median(tbl$age,na.rm=T))
                            }
                            if (!"anyAlt"%in%altTypeUniq1) {
                                tbl=cbind(tbl,anyAlt=apply(arrayData2,2,function(x) {as.integer(any(x!=0,na.rm=T))}))
                            }
                            res=matrix(0,nrow=nrow(tbl),ncol=length(altTypeUniq1))
                            colnames(res)=altTypeUniq1
                            for (k in 1:length(altTypeUniq1)) {
                                res[apply(arrayData2,2,function(x){any(x==k,na.rm=T)}),k]=1
                            }
                            tbl=cbind(tbl,res)
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
                            if ("age"%in%names(tbl)) {
                                write.table(paste("Median age",round(median(tbl$age,na.rm=T),2),sep="\t"),file=fNameThis,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                            }
                            write.table("",file=fNameThis,append=T,col.names=F,row.names=F,sep="\t",quote=F)
                            if (subsetFlag[1]%in%c("_ucsf500","_ucsf500Fmi")) {
                                #names(tbl)[match(c("disOntTerm2"),names(tbl))]=c("primarySite")
                                #varList2=c("sex","primarySite","diffEB","cellSize","anyAlt")
                                varList2=c("primarySiteEB2","diffEB","cellSizeEB","anyAlt")
                            } else {
                                #varList2=c("gender","grade","disOntTerm","disOntTerm2","assayVersion","anyAlt",sort(unique(tbl$altType)))
                                #varList2=c("gender","disOntTerm2","assayVersion","anyAlt",sort(unique(tbl$altType)))
                                varList2=c("gender","disOntTerm2","assayVersion","anyAlt",altTypeUniq1)
                            }
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
                            if (is.na(nClust[2])) {
                                nClustThis=1
                                tbl=annCol[,colId[,1]]; names(tbl)=colId[,2]
                                tbl=cbind(tbl,order=1:nrow(annCol))
                            } else {
                                nClustThis=nClust[2]
                                clustId=cutree(clustC,k=nClustThis)[clustC$order]
                                k1=which(!duplicated(clustId))
                                for (k in 1:length(k1)) {
                                    clustId[which(clustId==clustId[k1[k]])]=paste("cluster",k,sep="")
                                }
                                tbl=annCol[clustC$order,colId[,1]]; names(tbl)=colId[,2]
                                tbl=cbind(tbl,clustId,order=1:nrow(annCol))
                            }
                            write.table(tbl, fNameThis, sep="\t", col.names=T, row.names=F, quote=F,append=T)
                        }
                        
                        if (!is.null(colRow) & !is.null(varFListAll)) {
                            for (varId in 1:length(varFListAll)) {
                                if (varFListAll[varId]%in%c("sd","mutPerc")) {
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
                                    if (varFList[varId]==c("sd")) lim=limSdSam
                                    grpUniq=lim[1]:lim[2]
                                    colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                    lim=(lim-1)/100
                                    heatmapColorBar(limit=lim,cols=c(colColUniq[c(length(colColUniq),1,median(1:length(colColUniq)))]))
                                } else if (varFListAll[varId]%in%c("mutPerc")) {
                                    x=round(annRowAll[,varFListAll[varId]])+1
                                    lim=range(x,na.rm=T)
                                    if (varFList[varId]==c("mutPerc")) lim=limPerc+1
                                    #if (length(grep("dist2class",varFList[varId]))==1) lim=limDist2classSam
                                    grpUniq=lim[1]:lim[2]
                                    colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                    lim=lim-1
                                    heatmapColorBar(limit=lim,cols=c(colColUniq[c(length(colColUniq),1,median(1:length(colColUniq)))]))
                                } else {
                                    x=as.character(annRowAll[,varFListAll[varId]]); x[x==""]=NA
                                    grpUniq=table(x)
                                    grpUniq=names(grpUniq)
                                    k=1:length(grpUniq)
                                    if (varFListAll[varId]%in%c("geneFamily")) {
                                        sampleColorLegend(tls=grpUniq[k],col=colListGene,legendTitle=varFNameAll[varId])
                                    } else {
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
                                }
                                dev.off()
                            }
                        }
                        
                        
                        ## -------------------
                        if (!is.null(colCol) & !is.null(varListAll)) {
                            for (varId in 1:length(varListAll)) {
                                if (varListAll[varId]=="mutPerc") {
                                    cat("(varListAll[varId]==mutPerc)\n")
                                    next
                                }
                                if (length(grep("dist2",varList[varId]))==1 | varList[varId]%in%c("ki67EB2")) {
                                    if (outFormat=="png") {
                                        png(paste("heatmapSampleColorBarLegend_",varListAll[varId],fNameOut2,".png",sep=""),width=480,height=140)
                                    } else {
                                        pdf(paste("heatmapSampleColorBarLegend_",varListAll[varId],fNameOut2,".pdf",sep=""))
                                    }
                                    x=round(annColAll[,varListAll[varId]])+1
                                    lim=range(x,na.rm=T)
                                    if (varList[varId]==c("mutPerc")) lim=limPerc+1
                                    if (varList[varId]==c("ki67EB2")) {
                                        lim=c(0,200)+1
                                        lim=c(0,100)+1
                                    }
                                    #if (length(grep("dist2class",varList[varId]))==1) lim=limDist2classSam
                                    grpUniq=lim[1]:lim[2]
                                    if (varList[varId]==c("ki67EB2")) {
                                        colColUniq <- maPalette(low="indianred1",high="darkred",k=100)
                                    } else {
                                        #colColUniq=gray(0:(length(grpUniq)-1)/length(grpUniq))
                                        #colColUniq=rev(gray(1:(length(grpUniq))/length(grpUniq)))
                                        colColUniq=rev(rgb((0:length(grpUniq))/length(grpUniq), green = 0, blue = 0, names = paste("red", 0:length(grpUniq), sep = ".")))
                                    }
                                    lim=lim-1
                                    heatmapColorBar(limit=lim,cols=c(colColUniq[c(length(colColUniq),1,median(1:length(colColUniq)))]))
                                } else {
                                    if (varListAll[varId]%in%c("disOntTerm2","primarySiteEB2")) {
                                        #x=as.character(annColAll[samId03,varListAll[varId]]); x[x==""]=NA
                                        x=as.character(annColAll[,varListAll[varId]]); x[x==""]=NA
                                        x2=as.character(annCol[,varList[varId]]); x2[x2==""]=NA
                                    } else if (varListAll[varId]=="gene") {
                                        x=x2=c("0: no alt","1: any alt")
                                    } else {
                                        x=as.character(annColAll[,varListAll[varId]]); x[x==""]=NA
                                        x2=as.character(annCol[,varList[varId]]); x2[x2==""]=NA
                                    }
                                    if (all(is.na(x2))) {
                                        cat("(all(is.na(x2)))\n")
                                        next
                                    }
                                    x2=table(x2)
                                    x2=names(x2)
                                    grp=table(x)
                                    grpUniq=names(grp)
                                    k=1:length(grpUniq)
                                    ttl=paste(grpUniq[k]," (",grp[k],")",sep="")
                                    ttl=grpUniq[k]
                                    if (clusterFlag[1]=="_supervised") {
                                        k=match(x2,grpUniq)
                                    } else {
                                        k=1:length(grpUniq)
                                    }
                                    if (!is.null(varFListAll) && varFListAll[varId]%in%c("sd","mutPerc")) {
                                        width = 480; height = 140
                                    } else {
                                        if (length(grpUniq)<6) {
                                            width = 480; height = 480
                                        } else {
                                            width = 560; height = 960
                                        }
                                    }
                                    if (outFormat=="png") {
                                        png(paste("heatmapSampleColorBarLegend_",varListAll[varId],fNameOut2,".png",sep=""),width=width,height=height)
                                    } else {
                                        pdf(paste("heatmapSampleColorBarLegend_",varListAll[varId],fNameOut2,".pdf",sep=""))
                                    }
                                    cexThis=NULL
                                    if (outFormat=="pdf" & (length(grpUniq)>15 | max(nchar(grpUniq))>20)) cexThis=1
                                    if (outFormat=="pdf") cexThis=1
                                    if (varList[varId]%in%c("disOntTerm2","primarySiteEB2")) {
                                        sampleColorLegend(tls=ttl[k],col=colListD[k],legendTitle=varNameAll[varId],cex=cexThis)
                                    } else if (varList[varId]%in%c("cellSize","cellSizeEB")) {
                                        sampleColorLegend(tls=ttl[k],col=colListC[k],legendTitle=varNameAll[varId],cex=cexThis)
                                    } else if (varList[varId]%in%c("diffEB")) {
                                        sampleColorLegend(tls=ttl[k],col=colListDi[k],legendTitle=varNameAll[varId],cex=cexThis)
                                    } else if (varList[varId]%in%c("grade")) {
                                        sampleColorLegend(tls=ttl[k],col=colListGr[k],legendTitle=varNameAll[varId],cex=cexThis)
                                    } else if (varList[varId]%in%c("tmbLevel")) {
                                        sampleColorLegend(tls=ttl[k],col=colListT[k],legendTitle=varNameAll[varId],cex=cexThis)
                                    } else {
                                        if (length(grpUniq)<=length(colList2)) {
                                            sampleColorLegend(tls=ttl[k],col=colList2[k],legendTitle=varNameAll[varId],cex=cexThis)
                                        } else if (length(grpUniq)<=length(colList)) {
                                            sampleColorLegend(tls=ttl[k],col=colList[k],legendTitle=varNameAll[varId],cex=cexThis)
                                        } else {
                                            sampleColorLegend(tls=ttl[k],col=rainbow(length(grpUniq))[k],legendTitle=varNameAll[varId],cex=cexThis)
                                        }
                                    }
                                }
                                dev.off()
                            }
                        }
                        
                        if (length(colHMAll[[1]])==1) {
                            if (outFormat=="png") {
                                png(paste("heatmapColorRange",fNameOut2,".png",sep=""),width=480,height=140)
                            } else {
                                pdf(paste("heatmapColorRange",fNameOut2,".pdf",sep=""))
                            }
                            heatmapColorBar=function(limit,cols=c("green","red","black")) {
                                try <- maPalette(high=cols[1], low=cols[2], mid=cols[3],k=15)
                                maColorBar(try, scale=limit,k=3)
                            }
                            heatmapColorBar(limit=limit1,cols=unlist(colHMAll))
                            dev.off()
                        } else {
                            width = 480; height = 480
                            if (outFormat=="png") {
                                png(paste("heatmapColorLegend",fNameOut2,".png",sep=""),width=width,height=height)
                            } else {
                                pdf(paste("heatmapColorLegend",fNameOut2,".pdf",sep=""))
                            }
                            ##grpUniq=c(sort(unique(clin1$altType)),"multiple")
                            #grpUniq=c(sort(unique(clin1$altType[!is.na(clin1$altType)])),"multiple")
                            grpUniq=sort(altTypeUniq2[,1])
                            cexThis=NULL
                            if (outFormat=="pdf") cexThis=1
                            sampleColorLegend(tls=grpUniq,col=colHMAll[[1]],legendTitle="alteration type",cex=cexThis)
                            dev.off()
                            
                        }
                        
                    }
                }
            }
        }
    }
}

if (F) {
    tbl=candGene0[,which(!names(candGene0)%in%c("family2"))]
    tbl$altInAnyPatient=0
    tbl$altInAnyPatient[which(tbl$gene%in%rownames(arrayData))]=1
    write.table(tbl,file=paste("actionableGene_for",fNameOut,".txt",sep=""),append=F,col.names=T,row.names=F,sep="\t",quote=F)
}
