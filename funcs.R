getFamilyLevelInfo=function(verbose=FALSE) {

    ## -------------------
    if (F) {
        #geneFamilyFlag=T
        candGene=c("BRAF","SETD2","MMR","ATRX","BRCA1","BRCA2")

        genesetFlag="_mmr"
        datGP_m=datGP
        switch(genesetFlag,
        "_mmr"={
            tbl=read.table(paste("docs/candGene.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            tbl=tbl[which(tbl$family=="Mismatch excision repair (MMR)"),]
            x=unique(c(tbl$gene,unlist(sapply(tbl$alias,function(x) {strsplit(x,", ")[[1]]},USE.NAMES=F))))
            i=which(toupper(rownames(datGP))%in%x)
            y=apply(datGP[i,],2,function(x) {as.integer(any(x==1,na.rm=T))})
            for (ii in i) {
                datGP_m[ii,]=y
            }
            rownames(datGP_m)[i]=toupper(sub("_","",genesetFlag))
        }
        )
        datGP_m=datGP_m[!duplicated(rownames(datGP_m)),]
        if (F) {
            k=which(rownames(datGP_m)=="MMR")
            x=which(datGP_m[k,]==1)[1:5]
            x=which(datGP_m[k,]==0)[1:5]
            datGP[i,x]
            datGP_m[k,x]
        }
    }

    ## -------------------
    #geneFamilyFlag=F
    #geneFamilyFlag=T

    candGene=read.table(paste("docs/candGene2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    candGene$family2=sapply(candGene$family,function(x) {strsplit(x," ")[[1]][1]},USE.NAMES=F)
    k=which(candGene$family2=="mTOR"); k=k[1]
    x2=candGene[k,]; x2$gene="MTOR"
    candGene=rbind(candGene[1:(k-1),],x2,candGene[k:nrow(candGene),])
    tbl=read.table(paste("docs/candGene_swiSnfMod.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    #tbl=read.table(paste("docs/candGene_swiSnfMod_v2.txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    #tbl$family2=tbl$family
    tbl=tbl[,match(names(candGene),names(tbl))]
    candGene=rbind(candGene,tbl)
    x=c("ATM","ATR")
    x=x[which(x%in%rownames(datGP))]
    tbl=candGene[1:length(x),]; tbl$gene=x
    tbl$family=tbl$family2="ATM.ATR"
    candGene=rbind(candGene,tbl)
    candGene=candGene[-which(candGene$family2=="Other" & candGene$gene=="ATM"),]
    x=c("TP53","RB1","CDKN2A","CDKN2B","MEN1","ATRX","DAXX","APC","RAS-MAPK")
    x=c("TP53","RB1","CDKN2A","CDKN2B","MEN1","ATRX","DAXX","APC")
    x=x[which(x%in%rownames(datGP))]
    tbl=candGene[1:length(x),]; tbl$gene=x
    #tbl=rbind(tbl,candGene[which(candGene$family=="RAS-MAPK Signaling Pathway"),])
    tbl$family=tbl$family2="RB, etc."
    tbl$family=tbl$family2=tbl$gene
    candGene=rbind(candGene,tbl)
    #tbl=candGene[grep("SWI.SNF",candGene$family2),]
    #tbl$family=tbl$family2="SWI/SNF"
    #candGene=rbind(candGene,tbl)
    candGene$family2=gsub("-|/|: |_",".",candGene$family2)
    x1=toupper(rownames(datGP))
    x2=toupper(candGene$gene)
    #x2=unique(toupper(candGene$gene))
    out=NULL
    k1=1
    for (k in which(!x2%in%x1)) {
        if (verbose) cat("\n\n=========== ",x2[k],":\n",sep="")
        i=grep(x2[k],x1)
        ii=which(x1[i]%in%x2)
        if (length(ii)!=0) {
            i=i[-ii]
        }
        if (length(i)==0) {
            if (verbose) cat("In our data: None")
        } else {
            if (verbose) cat("In our data:",paste(x1[i],collapse=", "))
            out2=candGene[i,]
            for (ll in 1:ncol(out2)) {
                out2[,ll]=NA
            }
            out2$family=candGene$family[k]
            out2$family2=candGene$family2[k]
            out2$gene=x1[i]
            out=rbind(out,candGene[k1:(k-1),],out2)
            k1=k+1
        }
    }
    if (k1<=nrow(candGene)) {
        out=rbind(out,candGene[k1:nrow(candGene),])
    }
    out$altType[which(is.na(out$altType))]=""
    candGene=out[which(toupper(out$gene)%in%toupper(rownames(datGP))),]
    candGene0=out
    candGene0$inData=0
    candGene0$inData[which(toupper(out$gene)%in%toupper(rownames(datGP)))]=1
    
    ## -------------------
    grpUniq3=unique(candGene$family2)
    datGP_m=datGP
    k=c()
    for (gId in 1:length(grpUniq3)) {
        i=which(toupper(rownames(datGP))%in%toupper(candGene$gene[which(candGene$family2==grpUniq3[gId])]))
        if (length(i)==0) k=c(k,gId)
        if (length(i)==1) {
            y=as.integer(datGP[i,]==1)
        } else {
            y=apply(datGP[i,],2,function(x) {as.integer(any(x==1,na.rm=T))})
        }
        for (ii in i) {
            datGP_m[ii,]=y
        }
        rownames(datGP_m)[i]=grpUniq3[gId]
    }
    datGP_m=datGP_m[!duplicated(rownames(datGP_m)),]
    
    tbl=candGene[grep("SWI.SNF",candGene$family2),]
    tbl$family="SWI/SNF"; tbl$family2="SWI.SNF"
    tbl2=candGene[c(grep("SWI.SNF",candGene$family2),which(candGene$family2=="Modifier.Histone")),]
    tbl2$family="SWI/SNF, Modifier: Histone"; tbl2$family2="SWI.SNF_Modifier.Histone"
    tbl=rbind(tbl,tbl2)
    tbl2=candGene[c(grep("SWI.SNF",candGene$family2),which(candGene$family2=="Modifier.Histone"),which(candGene$family2=="ATM.ATR")),]
    tbl2$family="SWI/SNF, Modifier: Histone, ATM/ATR"; tbl2$family2="SWI.SNF_Modifier.Histone_ATM.ATR"
    tbl=rbind(tbl,tbl2)
    grpUniq3=c("SWI.SNF","SWI.SNF_Modifier.Histone","SWI.SNF_Modifier.Histone_ATM.ATR")
    #out=datGP
    #k=c()
    out=NULL
    for (gId in 1:length(grpUniq3)) {
        switch(grpUniq3[gId],
        "SWI.SNF"={i=which(toupper(rownames(datGP))%in%toupper(candGene$gene[which(substr(candGene$family2,1,nchar("SWI.SNF."))=="SWI.SNF.")]))
        },
        "SWI.SNF_Modifier.Histone"={i=which(toupper(rownames(datGP))%in%toupper(candGene$gene[c(which(substr(candGene$family2,1,nchar("SWI.SNF."))=="SWI.SNF."),which(candGene$family2=="Modifier.Histone"))]))
        },
        "SWI.SNF_Modifier.Histone_ATM.ATR"={i=which(toupper(rownames(datGP))%in%toupper(candGene$gene[c(which(substr(candGene$family2,1,nchar("SWI.SNF."))=="SWI.SNF."),which(candGene$family2=="Modifier.Histone"),which(candGene$family2=="ATM.ATR"))]))
        }
        )
        i=unique(i)
        #if (length(i)==0) k=c(k,gId)
        if (length(i)==1) {
            y=as.integer(datGP[i,]==1)
        } else {
            y=apply(datGP[i,],2,function(x) {as.integer(any(x==1,na.rm=T))})
        }
        if (F) {
            for (ii in i) {
                out[ii,]=y
            }
            rownames(out)[i]=grpUniq3[gId]
        }
        out=rbind(out,y)
    }
    #k=match(grpUniq3,rownames(out))
    #nm=c(rownames(datGP_m),rownames(out)[k])
    #datGP_m=rbind(datGP_m,out[k,])
    nm=c(rownames(datGP_m),grpUniq3)
    datGP_m=rbind(datGP_m,out)
    rownames(datGP_m)=nm
    candGene=rbind(candGene,tbl)

    if (F) {
        i=which(toupper(rownames(datGP))%in%toupper(candGene$gene[which(candGene$family2=="MMR")]))
        k=which(rownames(datGP_m)=="MTOR")
        x=which(datGP_m[k,]==1)[1:5]
        x=which(datGP_m[k,]==0)[1:5]
        datGP[i,x]
        datGP_m[k,x]
    }
    
    list(candGene=candGene,datGP_m=datGP_m)
}


getCandidateGenes=function(genesetFlag="") {
    candGeneThis=NULL
    switch(genesetFlag,
    "_actionableGene"={
        candGeneThis=candGene[which(substr(candGene$family2,1,nchar("SWI.SNF"))!="SWI.SNF" & candGene$family2!="Modifier.Histone"),]
    },
    "_swiSnf"={
        candGeneThis=candGene[which(candGene$family2=="SWI.SNF"),]
    },
    "_swiSnfComp"={
        candGeneThis=candGene[which(substr(candGene$family2,1,nchar("SWI.SNF."))=="SWI.SNF."),]
    },
    "_histoneMod"={
        candGeneThis=candGene[which(candGene$family=="Modifier: Histone methylase/demethylase"),]
    },
    "_atmAtr"={
        x=c("ATM","ATR")
        x=x[which(x%in%rownames(datGP))]
        candGeneThis=candGene[1:length(x),]; candGeneThis$gene=x
        candGeneThis$family=candGeneThis$family2="ATM.ATR"
        candGeneThis=candGene[which(candGene$family2%in%c("ATM.ATR")),]
    },
    "_swiSnfHisMod"={
        candGeneThis=candGene[which(candGene$family2=="SWI.SNF" | candGene$family2%in%c("Modifier.Histone")),]
    },
    "_swiSnfPlusHisMod"={
        candGeneThis=candGene[which(candGene$family2=="SWI.SNF_Modifier.Histone"),]
    },
    "_swiSnfHisModAtmAtr"={
        candGeneThis=candGene[which(candGene$family2=="SWI.SNF" | candGene$family2%in%c("Modifier.Histone") | candGene$family2%in%c("ATM.ATR")),]
    },
    "_swiSnfPlusHisModPlusAtmAtr"={
        candGeneThis=candGene[which(candGene$family2=="SWI.SNF_Modifier.Histone_ATM.ATR"),]
    },
    "_rbEtc"={
        x=c("TP53","RB1","CDKN2A","CDKN2B","MEN1","ATRX","DAXX","APC","RAS-MAPK")
        x=x[which(x%in%rownames(datGP))]
        candGeneThis=candGene[1:length(x),]; candGeneThis$gene=x
        candGeneThis=rbind(candGeneThis,candGene[which(candGene$family=="RAS-MAPK Signaling Pathway"),])
        candGeneThis$family=candGeneThis$family2="RB, etc."
        k=match(c("TP53","RB1","CDKN2A","CDKN2B","MEN1","ATRX","DAXX","APC"),candGene$family2); k=k[!is.na(k)]
        k=c(k,which(candGene$family2=="RAS-MAPK"))
        candGeneThis=candGene[k,]
    },
    "_swiSnfCompEtc"={
        k=match(c("TP53","RB1","CDKN2A","CDKN2B","MEN1","ATRX","DAXX","APC"),candGene$family2); k=k[!is.na(k)]
        k=c(k,which(candGene$family2%in%c("RAS-MAPK")),which(substr(candGene$family2,1,nchar("SWI.SNF."))=="SWI.SNF."),which(candGene$family2%in%c("Modifier.Histone")),which(candGene$family2%in%c("ATM.ATR")))
        candGeneThis=candGene[k,]
    },
    "_swiSnfEtc"={
        k=match(c("TP53","RB1","CDKN2A","CDKN2B","MEN1","ATRX","DAXX","APC"),candGene$family2); k=k[!is.na(k)]
        k=c(k,which(candGene$family2%in%c("RAS-MAPK")),which(candGene$family2=="SWI.SNF"),which(candGene$family2%in%c("Modifier.Histone")),which(candGene$family2%in%c("ATM.ATR")))
        candGeneThis=candGene[k,]
    }
    )
    candGeneThis
}
