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
    setwd(paste(dirSrc2,"EmilyBergsland",sep=""))
}

####################################################################
####################################################################
fName1="_ucsf500"
fName1=""

if (F) {
    #clin2=read.table("docs/ucsf500/ritu.1.3.31.17_noTherapyTypes.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    clin2=read.table("docs/ucsf500/ritu.1.3.31.17_noTherapy_noMutComment.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    clin2[1:5,1]

    clinU=read.table("docs/ucsf500/HighGradeOutcomes.3.6.17.v3.ekb.deidentified. reviewed.17_RRedit_noTherapyTypes.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

    sort(names(clinU)[!tolower(names(clinU))%in%tolower(names(clin2))])
    sort(names(clin2)[!tolower(names(clin2))%in%tolower(names(clinU))])

    k=match(tolower(names(clinU)),tolower(names(clin2))); k1=which(!is.na(k)); k2=k[k1]
    n=nrow(clinU)
    for ( k in 1:length(k1)) {
        if (any(clinU[1:n,k1[k]]!=clin2[1:n,k2[k]],na.rm=T) | any(is.na(clinU[1:n,k1[k]])!=is.na(clin2[1:n,k2[k]]))) {
            cat("\n=============== ",k,": ",names(clinU)[k1[k]],", ",names(clin2)[k2[k]],"\n",sep="")
            if (sum(!duplicated(clinU[1:n,k1[k]]))<10) {
                print(table(clinU[1:n,k1[k]],clin2[1:n,k2[k]],exclude=NULL))
            } else {
                print(summary(clinU[1:n,k1[k]]))
                print(summary(clin2[1:n,k2[k]]))
            }
        }
    }

    clinU=read.table("docs/ucsf500/HighGradeOutcomes.3.6.17.v3.ekb.deidentified. reviewed.17..txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    clinU[1:5,1]
}

## Remove therapy types
clinU=read.table("docs/ucsf500/HighGradeOutcomes.3.6.17.v3.ekb.deidentified. reviewed.17_RRedit_noTherapyTypes.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
## Remove therapy types, no mutation comment
clinU=read.table("docs/ucsf500/ritu.1.3.31.17_noTherapy_noMutComment.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

## For HighGradeOutcomes.3.6.17.v3.ekb.deidentified. reviewed.17_RRedit_noTherapyTypes.txt
varInfo=data.frame(variable1=c("Sex","Current.Patient.Status","Age.at.Diagnosis","Survival.since.Initial.Diagnosis..mo.","Primary.Site","Path.Review..UCSF.or.outside.","Type.of.Path","Stage.at.Diagnosis","Survival.since.Stage.IV..mo.","EB.GRADE","EB.differentiation","EB.ki.67","EB.MR","EB.cell.type","EB.NOTEs","Differentiation","Large.Small.Cell","Grade","Mitotic.Rate","Ki67","Comments","Test..ucsf500.or.FMI.","Date.of.UCSF.Sample.Used","TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","chromosomal.copy.change.Y.N","X.NO.mutations..deletions..rearrangement.","CDKN2A","CDKN2B","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFbR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM6A","FGFR1","VHL","SMARCA4","X.CDK4..MDM2..FRS2.","NFE2L2","MYC","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB2","PTEN","BCORL1","CHD4","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3"),
    variable2=c("sex","currPatStatus","ageAtDiag","survSinceInitDiag","primarySite","pathReview","typeOfPath","stageAtDiag","survSinceStageIV","gradeEB","diffEB","ki67EB","mitoticRateEB","cellSizeEB","noteEB","diff","cellSize","grade","mitoticRate","ki67","note","testUCSF500orFMI","dateOfUCSFsampleUsed","TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","chromosomalCopyChange","noMutDelRearrange","CDKN2A","CDKN2B","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFbR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM6A","FGFR1","VHL","SMARCA4","CDK4_MDM2_FRS2","NFE2L2","MYC","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB2","PTEN","BCORL1","CHD4","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3"),stringsAsFactors=F)
## For ritu.1.3.31.17_noTherapy_noMutComment.txt
varInfo=data.frame(variable1=c("Path.Review..UCSF.or.outside.","EB.PRIMARY.SITE","EB.DIFFERENTIATION","EB.KI67","EB.MITOTIC.RATE","EB.CELL.TYPE","EB.COMMENTS","Stage.at.Diagnosis","Survival.since.Stage.IV..mo.","Grade","Comments","Test..ucsf500.or.FMI.","TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","chromosomal.copy.change.Y.N","X.NO.mutations..deletions..rearrangement.","CDKN1B","CDKN2A","CDKN2B","CDKN2C","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFBR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM5A","KDM6A","FGFR1","FGFR3","VHL","SMARCA4","X.CDK4..MDM2..FRS2.","NFE2L2","MYC","ERBB2","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB22","PTEN","BCORL1","CHD4","IDH1","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3","MYCN","BRCA2","FH","NBN","ARID2","AKT2","AXL","FGF19","FGF23","FGF3","FGF4","FGF6","AR","DNMT3A","AURKA","GNAS","ARFRP1","FLT3","LZTR1","MAP2K4","ZNF217","NF1","FANCA","RET","MCL1","BRCA1","EGFR","JAK3","TERT","FANCC","PTPRD","STK11","RNF43","BRIP1","BRIP2","PTPN11"),
    variable2=c("pathReview","primarySiteEB","diffEB","ki67EB","mitoticRateEB","cellSizeEB","noteEB","stageAtDiag","survSinceStageIV","grade","note","testUCSF500orFMI","TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","chrCopyChng","noMutDelRearr","CDKN1B","CDKN2A","CDKN2B","CDKN2C","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFBR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM5A","KDM6A","FGFR1","FGFR3","VHL","SMARCA4","CDK4_MDM2_FRS2","NFE2L2","MYC","ERBB2","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB22","PTEN","BCORL1","CHD4","IDH1","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3","MYCN","BRCA2","FH","NBN","ARID2","AKT2","AXL","FGF19","FGF23","FGF3","FGF4","FGF6","AR","DNMT3A","AURKA","GNAS","ARFRP1","FLT3","LZTR1","MAP2K4","ZNF217","NF1","FANCA","RET","MCL1","BRCA1","EGFR","JAK3","TERT","FANCC","PTPRD","STK11","RNF43","BRIP1","BRIP2","PTPN11"),stringsAsFactors=F)
k=match(varInfo[,1],names(clinU)); k1=which(!is.na(k)); k2=k[k1]
names(clinU)[k1]=varInfo[,2]

clinU$id=1:nrow(clinU)
for (k in 1:ncol(clinU)) {
    if (is.character(clinU[,k])) {
        clinU[,k]=gsub("\"","",clinU[,k])
        clinU[,k]=gsub(" ","",clinU[,k])
        clinU[which(clinU[,k]%in%c("N.R.","NR")),k]="NR"
    }
}

for (k in 1:ncol(clinU)) {
    if (is.character(clinU[,k])) {
        if (sum(grep("REF",clinU[,k]))!=0) {
            cat("\n\n================ ",k,", ",names(clinU)[k],": ",class(clinU[,k]),", ",sum(!duplicated(clinU[,k])),"\n",sep="")
        }
    } else if (is.numeric(clinU[,k])) {
        cat("\n\n================ ",k,", ",names(clinU)[k],": ",class(clinU[,k]),", ",sum(!duplicated(clinU[,k])),"\n",sep="")
    } else {
        #cat("\n\n================ ",k,", ",names(clinU)[k],": ",class(clinU[,k]),", ",sum(!duplicated(clinU[,k])),"\n",sep="")
    }
}


colId=which(names(clinU)%in%c("TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","CDKN2A","CDKN2B","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFbR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM6A","FGFR1","VHL","SMARCA4","CDK4_MDM2_FRS2","NFE2L2","MYC","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB2","PTEN","BCORL1","CHD4","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3"))
colId=which(names(clinU)%in%c("TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","CDKN1B","CDKN2A","CDKN2B","CDKN2C","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFBR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM5A","KDM6A","FGFR1","FGFR3","VHL","SMARCA4","CDK4_MDM2_FRS2","NFE2L2","MYC","ERBB2","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB22","PTEN","BCORL1","CHD4","IDH1","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3","MYCN","BRCA2","FH","NBN","ARID2","AKT2","AXL","FGF19","FGF23","FGF3","FGF4","FGF6","AR","DNMT3A","AURKA","GNAS","ARFRP1","FLT3","LZTR1","MAP2K4","ZNF217","NF1","FANCA","RET","MCL1","BRCA1","EGFR","JAK3","TERT","FANCC","PTPRD","STK11","RNF43","BRIP1","BRIP2","PTPN11"))
kk=c()
for (k in colId) {
    x=mean(is.na(clinU[,k]))
    if (x>=1) kk=c(kk,k)
}
table(colId%in%kk)

datGPU=t(apply(as.matrix(clinU[,colId]),c(1,2),function(x) {ifelse(is.na(x),0,as.integer(x!=""))}))
colnames(datGPU)=paste("anyAlt_",clinU$id,sep="")
sort(unique(c(as.matrix(clinU[,colId]))),na.last=T)
table(c(datGPU),c(as.matrix(clinU[,colId]))!="",exclude=NULL)
clinU=clinU[,-colId]
if ("dead"%in%names(clinU)) {
    clinU$dead=NA
    clinU$dead[which(substr(clinU$currPatStatus,1,nchar("Alive"))=="Alive")]=0
    clinU$dead[which(substr(clinU$currPatStatus,1,nchar("Deceased"))=="Deceased")]=1
}
if ("primarySite"%in%names(clinU)) {
    clinU$primarySite=tolower(clinU$primarySite)
}

if ("ki67EB"%in%names(clinU)) {
    sort(unique(clinU$ki67EB))
    clinU$ki67EB2=as.numeric(gsub("%|>| ","",clinU$ki67EB))
    j=grep(">",clinU$ki67EB)
    clinU$ki67EB2[j]=clinU$ki67EB2[j]+1
}

colId=c("diffEB","ki67EB2","cellSizeEB")
for (k in colId) {
    cat("\n\n================ ",k,": ",class(clinU[,k]),", ",sum(!duplicated(clinU[,k])),"\n",sep="")
    if (sum(!duplicated(clinU[,k]))<=10) {
        if (is.character(clinU[,k])) clinU[,k]=tolower(clinU[,k])
        print(table(clinU[,k],exclude=NULL))
        if (is.character(clinU[,k])) clinU[,k]=tolower(clinU[,k])
    } else {
        print(summary(as.numeric(clinU[,k])))
    }
}

res=apply(datGPU,1,function(x) mean(x==1,na.rm=T))
res=apply(datGPU,1,function(x) sum(x==1,na.rm=T))

save(clinU,datGPU,fName1,file="tmp_ucsf500.RData")

## -------------------------------------
library(coin)

samSize=5

varList1=c("diffEB","cellSizeEB"); varName1=varList1
varList2=c("ki67EB2"); varName2=c("ki67EB")
pdf("boxplot_ucsf500.pdf")
for (subsetFlag in c("_ucsf500","_ucsf500Fmi")) {
    j=1:nrow(clinU)
    heading=paste("All samples",sep="")
    switch(subsetFlag,
    "_ucsf500"={
        j=which(clinU$testUCSF500orFMI%in%c("UCSF500"))
        heading=paste("UCSF 500: All samples",sep="")
    },
    "_ucsf500Fmi"={
        j=which(clinU$testUCSF500orFMI%in%c("UCSF500","FMI"))
        heading=paste("UCSF 500, FMI: All samples",sep="")
    }
    )
    #par(mfrow=c(2,1))
    lim=c(0,100)
    for (vId2 in 1:length(varList2)) {
        for (vId1 in 1:length(varList1)) {
            x=clinU[j,varList1[vId1]]
            ttl=sort(unique(x))
            if (varList1[vId1]=="diffEB") {
                x[which(x=="poor")]="1:Poor"
                x[which(x=="well")]="2:Well"
                x[which(x=="nr")]="99:NR"
                ttl=sapply(sort(unique(x)),function(x) strsplit(x,":")[[1]][2],USE.NAMES=F)
            } else if (varList1[vId1]=="cellSizeEB") {
                x[which(x=="small")]="1:SC"
                x[which(x=="large")]="2:LC"
                x[which(x=="nr")]="99:NR"
                ttl=sapply(sort(unique(x)),function(x) strsplit(x,":")[[1]][2],USE.NAMES=F)
            }
            ttl=paste(ttl," (",table(x),")",sep="")
            grpUniq=sort(unique(x))
            pv=rep(NA,length(grpUniq))
            heading2="pv: "
            k=1
            for (grpId1 in 1:(length(grpUniq)-1)) {
                for (grpId2 in (grpId1+1):length(grpUniq)) {
                    jj=which(x%in%grpUniq[c(grpId1,grpId2)])
                    if (all(table(x[jj])>=samSize)) {
                        pv[k]=pvalue(wilcox_test(clinU[j[jj],varList2[vId2]]~as.factor(x[jj])))
                    }
                    heading2=paste(heading2,signif(pv[k],2)," (",substr(ttl[grpId1],1,1)," vs ",substr(ttl[grpId2],1,1),"), ",sep="")
                    k=k+1
                }
            }
            heading2=substr(heading2,1,nchar(heading2)-1)
            boxplot(clinU[j,varList2[vId2]]~x,names=ttl,main=paste(heading,"\n",heading2,sep=""),ylim=lim,xlab=varName1[vId1],ylab=varName2[vId2])
        }
    }
}
dev.off()


## -------------------------------------

library(survival)

time=clinU$survSinceStageIV; event=clinU$dead

varList=c("diffEB","cellSizeEB")
pdf("kmPlot.pdf")
for (varId in 1:length(varList)) {
    x=clinU[,varList[varId]]
    fit=survfit(Surv(time, event) ~ x)
    plot(fit, xlab="Survival since stage IV (mo)", lty = 2:3)
    legend(30, .9, sort(unique(x)), title=varList[varId], lty = 2:3)
    title("Kaplan-Meier Curves\nfor UCSF500 data")
}
dev.off()
## -------------------------------------
