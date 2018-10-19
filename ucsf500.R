####################################################################
####################################################################
## Section 1

computerFlag="cluster"
computerFlag=""

## ----------------------------------------------
if (computerFlag=="cluster") {
    setwd("/home/royr/project/bergslandE")
} else {
    dirSrc="/Users/royr/UCSF/"
    dirSrc2=dirSrc
    setwd(paste(dirSrc2,"bergslandE",sep=""))
}

####################################################################
####################################################################
datVerFlag="_20170306"
datVerFlag="_20170331"
datVerFlag="_20170623"
datVerFlag="_20171019"
datVerFlag="_20180925"
datVerFlag="_20180929"

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

if (datVerFlag=="_20170306") {
    ## Remove columns therapy types, giving error when reading table
    clinU=read.table("docs/ucsf500/HighGradeOutcomes.3.6.17.v3.ekb.deidentified. reviewed.17_RRedit_noTherapyTypes.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    varInfo=data.frame(variable1=c("Sex","Current.Patient.Status","Age.at.Diagnosis","Survival.since.Initial.Diagnosis..mo.","Primary.Site","Path.Review..UCSF.or.outside.","Type.of.Path","Stage.at.Diagnosis","Survival.since.Stage.IV..mo.","EB.GRADE","EB.differentiation","EB.ki.67","EB.MR","EB.cell.type","EB.NOTEs","Differentiation","Large.Small.Cell","Grade","Mitotic.Rate","Ki67","Comments","Test..ucsf500.or.FMI.","Date.of.UCSF.Sample.Used","TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","chromosomal.copy.change.Y.N","X.NO.mutations..deletions..rearrangement.","CDKN2A","CDKN2B","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFbR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM6A","FGFR1","VHL","SMARCA4","X.CDK4..MDM2..FRS2.","NFE2L2","MYC","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB2","PTEN","BCORL1","CHD4","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3"),
        variable2=c("sex","currPatStatus","ageAtDiag","survSinceInitDiag","primarySite","pathReview","typeOfPath","stageAtDiag","survSinceStageIV","gradeEB","diffEB","ki67EB","mitoticRateEB","cellSizeEB","noteEB","diff","cellSize","grade","mitoticRate","ki67","note","testUCSF500orFMI","dateOfUCSFsampleUsed","TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","chromosomalCopyChange","noMutDelRearrange","CDKN2A","CDKN2B","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFbR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM6A","FGFR1","VHL","SMARCA4","CDK4_MDM2_FRS2","NFE2L2","MYC","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB2","PTEN","BCORL1","CHD4","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3"),stringsAsFactors=F)
} else if (datVerFlag=="_20170331") {
    ## Remove columns therapy types, no mutation comment, giving error when reading table
    clinU=read.table("docs/ucsf500/ritu.1.3.31.17_noTherapy_noMutComment.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    varInfo=data.frame(variable1=c("Path.Review..UCSF.or.outside.","EB.PRIMARY.SITE","EB.DIFFERENTIATION","EB.KI67","EB.MITOTIC.RATE","EB.CELL.TYPE","EB.COMMENTS","Stage.at.Diagnosis","Survival.since.Stage.IV..mo.","Grade","Comments","Test..ucsf500.or.FMI.","TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","chromosomal.copy.change.Y.N","X.NO.mutations..deletions..rearrangement.","CDKN1B","CDKN2A","CDKN2B","CDKN2C","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFBR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM5A","KDM6A","FGFR1","FGFR3","VHL","SMARCA4","X.CDK4..MDM2..FRS2.","NFE2L2","MYC","ERBB2","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB22","PTEN","BCORL1","CHD4","IDH1","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3","MYCN","BRCA2","FH","NBN","ARID2","AKT2","AXL","FGF19","FGF23","FGF3","FGF4","FGF6","AR","DNMT3A","AURKA","GNAS","ARFRP1","FLT3","LZTR1","MAP2K4","ZNF217","NF1","FANCA","RET","MCL1","BRCA1","EGFR","JAK3","TERT","FANCC","PTPRD","STK11","RNF43","BRIP1","BRIP2","PTPN11"),
        variable2=c("pathReview","primarySiteEB","diffEB","ki67EB","mitoticRateEB","cellSizeEB","noteEB","stageAtDiag","survSinceStageIV","grade","note","testUCSF500orFMI","TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","chrCopyChng","noMutDelRearr","CDKN1B","CDKN2A","CDKN2B","CDKN2C","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFBR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM5A","KDM6A","FGFR1","FGFR3","VHL","SMARCA4","CDK4_MDM2_FRS2","NFE2L2","MYC","ERBB2","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB22","PTEN","BCORL1","CHD4","IDH1","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3","MYCN","BRCA2","FH","NBN","ARID2","AKT2","AXL","FGF19","FGF23","FGF3","FGF4","FGF6","AR","DNMT3A","AURKA","GNAS","ARFRP1","FLT3","LZTR1","MAP2K4","ZNF217","NF1","FANCA","RET","MCL1","BRCA1","EGFR","JAK3","TERT","FANCC","PTPRD","STK11","RNF43","BRIP1","BRIP2","PTPN11"),stringsAsFactors=F)
} else if (datVerFlag=="_20170623") {
    #clinU=read.table("docs/ucsf500/Ritu_HGN1_DATA_2017-06-23_1442.csv",sep=",",h=T,quote="",comment.char="",as.is=T,fill=T)
    #clinU=read.table("docs/ucsf500/Ritu (2) HGN1_DATA_2017-06-23_1637.csv",sep=",",h=T,quote="",comment.char="",as.is=T,fill=T)
    clinU=read.table("docs/ucsf500/Ritu_Final_HGN1_ImportTemplate_2017-06-24.csv",sep=",",h=T,quote="",comment.char="",as.is=T,fill=T)
    varInfo=data.frame(variable1=c("caseid","sex","status","dateofdiagnosis","survivalsincedx","agedx","ebprimarysite","metsites___0","metsites___1","metsites___2","metsites___3","pathreviewcenter","pathdate","ebdifferentiation","ebki67","eb_celltype","diagnosis_stage","grade","test","treatment_firstline___0","treatment_firstline___1","treatment_firstline___2","treatment_firstline___3","treatment_firstline___4","treatment_firstline___5","treatment_firstline___6","treatment_firstline___7","treatment_firstline___8","treatment_firstline___9","treatment_firstline___10","treatment_firstline___11","treatment_secondline___0","treatment_secondline___1","treatment_secondline___2","treatment_secondline___3","treatment_secondline___4","treatment_secondline___5","treatment_secondline___6","treatment_secondline___7","treatment_secondline___8","treatment_secondline___9","treatment_secondline___10","treatment_secondline___11","treatment_secondline___12","treatment_secondline___13","treatment_secondline___14","treatment_secondline___15","treatment_secondline___16","treatment_secondline___17","treatment_secondline___18","treatment_thirdline___0","treatment_thirdline___1","treatment_thirdline___2","treatment_thirdline___3","treatment_thirdline___4","treatment_thirdline___5","treatment_thirdline___6","treatment_thirdline___7","treatment_thirdline___8","treatment_thirdline___9","treatment_thirdline___10","tp53","rb1","crc_side","comments","variables_complete"),
        variable2=c("caseId","sex","status","dateDx","survSinceDx","ageDx","primarySiteEB","metsites0","metsites1","metsites2","metsites3","pathReviewCenter","pathDate","diffEB","ki67EB","cellSizeEB","stageAtDiag","grade","testUCSF500orFMI","treatFirstline0","treatFirstline1","treatFirstline2","treatFirstline3","treatFirstline4","treatFirstline5","treatFirstline6","treatFirstline7","treatFirstline8","treatFirstline9","treatFirstline10","treatFirstline11","treatSecondline0","treatSecondline1","treatSecondline2","treatSecondline3","treatSecondline4","treatSecondline5","treatSecondline6","treatSecondline7","treatSecondline8","treatSecondline9","treatSecondline10","treatSecondline11","treatSecondline12","treatSecondline13","treatSecondline14","treatSecondline15","treatSecondline16","treatSecondline17","treatSecondline18","treatThirdline0","treatThirdline1","treatThirdline2","treatThirdline3","treatThirdline4","treatThirdline5","treatThirdline6","treatThirdline7","treatThirdline8","treatThirdline9","treatThirdline10","TP53","RB1","crcSide","comments","variablesComplete"),stringsAsFactors=F)
} else if (datVerFlag=="_20171019") {
    ## Remove columns therapy types, no mutation comment, giving error when reading table
    #clinU1=read.table("docs/ucsf500/RITU_HighGradeOutcomes_JW_10Oct2017_noTherapy_noMutNote_noMetSite.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    clinU=read.table("docs/ucsf500/Ritu_HighGradeOutcomes_JW_19Oct2017_noTherapy_noMutNote_noMetSite_noNote.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    varInfo=data.frame(variable1=c("Case.ID","Sex","DOB","Diagnosis.Date","Current.Patient.Status","Last.Contact.Date","Age.at.Diagnosis","Survival.since.Initial.Diagnosis..mo.","Primary.Site","Path.Review..UCSF.or.outside.","Date.of.Path.Used","Type.of.Path","EB.PRIMARY.SITE","EB.DIFFERENTIATION","EB.KI67","EB.MITOTIC.RATE","EB.CELL.TYPE","EB.COMMENTS","Stage.at.Diagnosis","Survival.since.Stage.IV..mo.","Differentiation","Large.Small.Cell","Grade","Mitotic.Rate","Ki67","Metastatic.Site.s.","Comments","Therapy.Types","Test..ucsf500.or.FMI.","Date.of.UCSF.Sample.Used","TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","chromosomal.copy.change.Y.N","X.NO.mutations..deletions..rearrangement.","CDKN1B","CDKN2A","CDKN2B","CDKN2C","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFBR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM5A","KDM6A","FGFR1","FGFR3","VHL","SMARCA4","X.CDK4..MDM2..FRS2.","NFE2L2","MYC","ERBB2","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB22","PTEN","BCORL1","CHD4","IDH1","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3","MYCN","BRCA2","FH","NBN","ARID2","AKT2","AXL","FGF19","FGF23","FGF3","FGF4","FGF6","AR","DNMT3A","AURKA","GNAS","ARFRP1","FLT3","LZTR1","MAP2K4","ZNF217","NF1","FANCA","RET","MCL1","BRCA1","EGFR","JAK3","TERT","FANCC","PTPRB","PTPRD","STK11","RNF43","PTPN11","BRIP1","CREBBP","AXIN2","MGA","PBRM1","SDHD","STAG2","HRAS","ASXL1","BRD4","FAM123B","Genetic.Mutations.Comments"),
    variable2=c("caseId","sex","dob","dateDx","status","dateLastContact","ageDx","survSinceDx","primarySite","pathReviewCenter","pathDate","typeOfPath","primarySiteEB","diffEB","ki67EB","mitoticRateEB","cellSizeEB","noteEB","stageAtDiag","survSinceStageIV","diff","cellSize","grade","mitoticRate","ki67","Metastatic.Site.s.","note","therapyType","testUCSF500orFMI","dateOfUCSFsampleUsed","TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","chromosomalCopyChange","noMutDelRearrange","CDKN1B","CDKN2A","CDKN2B","CDKN2C","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFBR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM5A","KDM6A","FGFR1","FGFR3","VHL","SMARCA4","CDK4_MDM2_FRS2","NFE2L2","MYC","ERBB2","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB22","PTEN","BCORL1","CHD4","IDH1","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3","MYCN","BRCA2","FH","NBN","ARID2","AKT2","AXL","FGF19","FGF23","FGF3","FGF4","FGF6","AR","DNMT3A","AURKA","GNAS","ARFRP1","FLT3","LZTR1","MAP2K4","ZNF217","NF1","FANCA","RET","MCL1","BRCA1","EGFR","JAK3","TERT","FANCC","PTPRB","PTPRD","STK11","RNF43","PTPN11","BRIP1","CREBBP","AXIN2","MGA","PBRM1","SDHD","STAG2","HRAS","ASXL1","BRD4","FAM123B","noteMut"),stringsAsFactors=F)
} else if (datVerFlag=="_20180925") {
    nm=read.table("docs/ucsf500/Ritu - HG Data with Diff Updates (25SEP18).txt",sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=1)
    nm=unlist(nm[1,]); names(nm)=NULL
    nm1=nm
    #nm=gsub(" +","_",nm)
    nm=gsub(" +","",nm)
    clinU=read.table("docs/ucsf500/Ritu - HG Data with Diff Updates (25SEP18).txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    names(clinU)=nm
    varInfo=data.frame(variable1=c("CaseID","Sex","DOB","Diagnosis.Date","CurrentPatientStatus(Updated:)","Last.Contact.Date","AgeatDiagnosis","Survival.since.Initial.Diagnosis..mo.","Primary.Site","Path.Review..UCSF.or.outside.","Date.of.Path.Used","Type.of.Path","EBPRIMARYSITE","EB.DIFFERENTIATION","EBKI67","EBMITOTICRATE","EBCELLTYPE","EB.COMMENTS","NJDIFFERENTIATION","Stage.at.Diagnosis","Survival.since.Stage.IV..mo.","Differentiation","Large.Small.Cell","Grade","Mitotic.Rate","Ki67","Metastatic.Site.s.","Comments","Therapy.Types","Test(ucsf500orFMI)","Date.of.UCSF.Sample.Used","TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","chromosomal.copy.change.Y.N","X.NO.mutations..deletions..rearrangement.","CDKN1B","CDKN2A","CDKN2B","CDKN2C","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFBR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM5A","KDM6A","FGFR1","FGFR3","VHL","SMARCA4","X.CDK4..MDM2..FRS2.","NFE2L2","MYC","ERBB2","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB22","PTEN","BCORL1","CHD4","IDH1","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3","MYCN","BRCA2","FH","NBN","ARID2","AKT2","AXL","FGF19","FGF23","FGF3","FGF4","FGF6","AR","DNMT3A","AURKA","GNAS","ARFRP1","FLT3","LZTR1","MAP2K4","ZNF217","NF1","FANCA","RET","MCL1","BRCA1","EGFR","JAK3","TERT","FANCC","PTPRB","PTPRD","STK11","RNF43","PTPN11","BRIP1","CREBBP","AXIN2","MGA","PBRM1","SDHD","STAG2","HRAS","ASXL1","BRD4","FAM123B","MYB","KDM5C","CDK4","MDM2","FRS2","FGFR2","CDK12","EMSY","SOX2","ELF3","BAP1","CBFB","DCC","FUPB1","EP300","NTRK1","SMAD2","TSHR","TET2","FGF10","HNF1A","JAK2","SMARCB1","MSH6","MSH7","GATA2","FUBP1","TCF7L2","Genetic.Mutations.Comments"),
    variable2=c("caseId","sex","dob","dateDx","status","dateLastContact","ageDx","survSinceDx","primarySite","pathReviewCenter","pathDate","typeOfPath","primarySiteEB","diffEB","ki67EB","mitoticRateEB","cellSizeEB","noteEB","diffNJ","stageAtDiag","survSinceStageIV","diff","cellSize","grade","mitoticRate","ki67","Metastatic.Site.s.","note","therapyType","testUCSF500orFMI","dateOfUCSFsampleUsed","TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","chromosomalCopyChange","noMutDelRearrange","CDKN1B","CDKN2A","CDKN2B","CDKN2C","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFBR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM5A","KDM6A","FGFR1","FGFR3","VHL","SMARCA4","CDK4_MDM2_FRS2","NFE2L2","MYC","ERBB2","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB22","PTEN","BCORL1","CHD4","IDH1","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3","MYCN","BRCA2","FH","NBN","ARID2","AKT2","AXL","FGF19","FGF23","FGF3","FGF4","FGF6","AR","DNMT3A","AURKA","GNAS","ARFRP1","FLT3","LZTR1","MAP2K4","ZNF217","NF1","FANCA","RET","MCL1","BRCA1","EGFR","JAK3","TERT","FANCC","PTPRB","PTPRD","STK11","RNF43","PTPN11","BRIP1","CREBBP","AXIN2","MGA","PBRM1","SDHD","STAG2","HRAS","ASXL1","BRD4","FAM123B","MYB","KDM5C","CDK4","MDM2","FRS2","FGFR2","CDK12","EMSY","SOX2","ELF3","BAP1","CBFB","DCC","FUPB1","EP300","NTRK1","SMAD2","TSHR","TET2","FGF10","HNF1A","JAK2","SMARCB1","MSH6","MSH7","GATA2","FUBP1","TCF7L2","noteMut"),stringsAsFactors=F)
} else if (datVerFlag=="_20180929") {
    nm=read.table("docs/ucsf500/Ritu - HG Data with Diff Updates (29SEP18).txt",sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=1)
    nm=unlist(nm[1,]); names(nm)=NULL
    nm1=nm
    #nm=gsub(" +","_",nm)
    nm=gsub(" +","",nm)
    clinU=read.table("docs/ucsf500/Ritu - HG Data with Diff Updates (29SEP18).txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    names(clinU)=nm
    varInfo=data.frame(variable1=c("CaseID","Sex","DOB","Diagnosis.Date","CurrentPatientStatus(Updated:)","Last.Contact.Date","AgeatDiagnosis","Survival.since.Initial.Diagnosis..mo.","Primary.Site","Path.Review..UCSF.or.outside.","Date.of.Path.Used","Type.of.Path","EBPRIMARYSITE","EB.DIFFERENTIATION","EBKI67","EBMITOTICRATE","EBCELLTYPE","EB.COMMENTS","NJDIFFERENTIATION","Stage.at.Diagnosis","Survival.since.Stage.IV..mo.","Differentiation","Large.Small.Cell","Grade","Mitotic.Rate","Ki67","Metastatic.Site.s.","Comments","Therapy.Types","Test(ucsf500orFMI)","Date.of.UCSF.Sample.Used","TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","chromosomal.copy.change.Y.N","X.NO.mutations..deletions..rearrangement.","CDKN1B","CDKN2A","CDKN2B","CDKN2C","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFBR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM5A","KDM6A","FGFR1","FGFR3","VHL","SMARCA4","X.CDK4..MDM2..FRS2.","NFE2L2","MYC","ERBB2","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB22","PTEN","BCORL1","CHD4","IDH1","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3","MYCN","BRCA2","FH","NBN","ARID2","AKT2","AXL","FGF19","FGF23","FGF3","FGF4","FGF6","AR","DNMT3A","AURKA","GNAS","ARFRP1","FLT3","LZTR1","MAP2K4","ZNF217","NF1","FANCA","RET","MCL1","BRCA1","EGFR","JAK3","TERT","FANCC","PTPRB","PTPRD","STK11","RNF43","PTPN11","BRIP1","CREBBP","AXIN2","MGA","PBRM1","SDHD","STAG2","HRAS","ASXL1","BRD4","FAM123B","MYB","KDM5C","CDK4","MDM2","FRS2","FGFR2","CDK12","EMSY","SOX2","ELF3","BAP1","CBFB","DCC","FUPB1","EP300","NTRK1","SMAD2","TSHR","TET2","FGF10","HNF1A","JAK2","SMARCB1","MSH6","MSH7","GATA2","FUBP1","TCF7L2","Genetic.Mutations.Comments"),
    variable2=c("caseId","sex","dob","dateDx","status","dateLastContact","ageDx","survSinceDx","primarySite","pathReviewCenter","pathDate","typeOfPath","primarySiteEB","diffEB","ki67EB","mitoticRateEB","cellSizeEB","noteEB","diffNJ","stageAtDiag","survSinceStageIV","diff","cellSize","grade","mitoticRate","ki67","Metastatic.Site.s.","note","therapyType","testUCSF500orFMI","dateOfUCSFsampleUsed","TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","chromosomalCopyChange","noMutDelRearrange","CDKN1B","CDKN2A","CDKN2B","CDKN2C","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFBR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM5A","KDM6A","FGFR1","FGFR3","VHL","SMARCA4","CDK4_MDM2_FRS2","NFE2L2","MYC","ERBB2","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB22","PTEN","BCORL1","CHD4","IDH1","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3","MYCN","BRCA2","FH","NBN","ARID2","AKT2","AXL","FGF19","FGF23","FGF3","FGF4","FGF6","AR","DNMT3A","AURKA","GNAS","ARFRP1","FLT3","LZTR1","MAP2K4","ZNF217","NF1","FANCA","RET","MCL1","BRCA1","EGFR","JAK3","TERT","FANCC","PTPRB","PTPRD","STK11","RNF43","PTPN11","BRIP1","CREBBP","AXIN2","MGA","PBRM1","SDHD","STAG2","HRAS","ASXL1","BRD4","FAM123B","MYB","KDM5C","CDK4","MDM2","FRS2","FGFR2","CDK12","EMSY","SOX2","ELF3","BAP1","CBFB","DCC","FUPB1","EP300","NTRK1","SMAD2","TSHR","TET2","FGF10","HNF1A","JAK2","SMARCB1","MSH6","MSH7","GATA2","FUBP1","TCF7L2","noteMut"),stringsAsFactors=F)
    if (F) {
        altInfo=read.table("docs/ucsf500/Ritu - Mutation Categories for High Grades.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        names(altInfo)=c("caseId","testUCSF500orFMI","alt","altType")
        altInfo$alt=gsub("\"","",altInfo$alt)
        altInfo$alt=gsub("−","-",altInfo$alt)
        x=t(sapply(altInfo$alt,function(x) {
            y=gsub("\"","",x)
            y=strsplit(y," ",fixed=T)[[1]]
            y=c(y[1],paste(y[2:length(y)],collapse=" "))
            y
        },USE.NAMES=F))
        altInfo$geneSym=x[,1]
        altInfo$alt1=altInfo$alt
        altInfo$alt=x[,2]
        }
    
    altInfo=read.table("docs/ucsf500/Ritu - Mutations & Categories with Levels of Evidence (29SEPT18).txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
    names(altInfo)=c("caseId","testUCSF500orFMI","alt","altType","evidLevel","indDrug")
    if (F) {
        x=t(sapply(altInfo$alt,function(x) {
            y=gsub("\"","",x)
            y=strsplit(y," ",fixed=T)[[1]]
            y=c(y[1],paste(y[2:length(y)],collapse=" "))
            y
        },USE.NAMES=F))
    }
    altInfo$altType=tolower(altInfo$altType)
    altInfo$altType[which(altInfo$altType=="deep deletion")]="deletion"
    altInfo$altType[which(altInfo$altType=="loss??")]="loss"
    altInfo$altType[which(altInfo$altType=="loss")]="deletion"
    tmpC=rep("",nrow(altInfo))
    out=data.frame(gene=tmpC,alt=tmpC,stringsAsFactors=F)
    for (i in 1:nrow(altInfo)) {
        x=altInfo$alt[i]
        y=gsub("\"","",x)
        y=try(strsplit(y," ",fixed=T)[[1]])
        y=c(y[1],paste(y[2:length(y)],collapse=" "))
        if (y[1]%in%c("CDKN2A,B","CDKN2A/B")) y[1]="CDKN2AB"
        if (y[1]=="none") y[1]=NA
        y[1]=sub(",","",y[1])
        y[1]=sub("/|-","_",y[1])
        out[i,]=y
    }
    altInfo=cbind(altInfo,out)
    altInfo[which(is.na(altInfo$gene)),]
    id=paste(clinU$CaseID,sep="")
    datGPU=matrix("",nrow=sum(!duplicated(altInfo$gene)),ncol=nrow(clinU),dimnames=list(unique(altInfo$gene),paste("anyAlt_",id,sep="")))
    j1=which(!is.na(altInfo$caseId)); j2=c(j1[2:length(j1)]-1,nrow(altInfo))
    id=paste("anyAlt_",altInfo$caseId,sep="")
    for (j in 1:length(j1)) {
        i=match(altInfo$gene[j1[j]:j2[j]],rownames(datGPU)); i1=which(!is.na(i)); i2=i[i1]
        jj=which(colnames(datGPU)==id[j1[j]])
        datGPU[i2,jj]=altInfo$altType[j1[j]:j2[j]][i1]
    }
    jj=c()
    for (j in 1:ncol(datGPU)) if (any(datGPU[,j]=="",na.rm=T)) {
        jj=c(jj,j)
        i=which(datGPU[,20]=="")
        #print(rownames(datGPU)[i])
    }
    datGPU[datGPU==""]=NA
    
    altTypeUniq1=sort(unique(altInfo$altType))
    altTypeUniq1=altTypeUniq1[which(!altTypeUniq1%in%c(""))]
    id=paste(clinU$CaseID,sep="")
    out=matrix(0,nrow=sum(!duplicated(altInfo$gene)),ncol=nrow(clinU)*(length(altTypeUniq1)+1),dimnames=list(unique(altInfo$gene),paste(c("anyAlt",altTypeUniq1),rep(id,each=length(altTypeUniq1)+1),sep="_")))
    j1=which(!is.na(altInfo$caseId)); j2=c(j1[2:length(j1)]-1,nrow(altInfo))
    for (a1 in 1:length(altTypeUniq1)) {
        id=paste(altTypeUniq1[a1],"_",altInfo$caseId,sep="")
        for (j in 1:length(j1)) {
            jj=j1[j]:j2[j]
            jj=jj[which(altInfo$altType[jj]==altTypeUniq1[a1])]
            i=which(rownames(out)%in%altInfo$gene[jj])
            jj=which(colnames(out)==id[j1[j]])
            out[i,jj]=1
        }
    }
    j1=which(sapply(colnames(out),function(x){strsplit(x,"_")[[1]][1]},USE.NAMES=F)=="anyAlt")
    for (a1 in 1:length(altTypeUniq1)) {
        j2=which(sapply(colnames(out),function(x){strsplit(x,"_")[[1]][1]},USE.NAMES=F)==altTypeUniq1[a1])
        out[,j1]=out[,j1]+out[,j2]
    }
    out[out>0]=1
    datGPU=out
}

k=match(varInfo[,1],names(clinU)); k1=which(!is.na(k)); k2=k[k1]
names(clinU)[k2]=varInfo[k1,2]

if (datVerFlag=="_20180929") {
    clinU$id=clinU$caseId
    clinU$primarySiteEB[which(clinU$primarySiteEB=="OTHER GI (GE)")]="OTHER GI"
    #clinU$cellSizeEB=sub(" & ","/",clinU$cellSizeEB)
} else {
    clinU$id=1:nrow(clinU)
}
if ("caseId"%in%names(clinU)) {
    if (is.integer(clinU$caseId)) clinU$caseId=paste("X",clinU$id,sep="")
} else {
    clinU$caseId=paste("X",clinU$id,sep="")
}
for (k in 1:ncol(clinU)) {
    if (is.character(clinU[,k])) {
        clinU[,k]=gsub("\"","",clinU[,k])
        clinU[,k]=gsub("−","-",clinU[,k])
        #clinU[,k]=gsub(" ","",clinU[,k])
        clinU[,k]=sub(" +$","",clinU[,k])
        clinU[which(clinU[,k]%in%c("N.R.","NR")),k]="NR"
    }
}
k1=match(c("diff","cellSize","ki67","primarySite","mitoticRate"),names(clinU))
if (any(!is.na(k1))) {
    if (datVerFlag%in%c("_20171019")) {k2=match(c("diffEB","cellSizeEB","ki67EB","primarySiteEB","mitoticRateEB"),names(clinU))
    } else {k2=match(c("diffNJ","cellSizeEB","ki67EB","primarySiteEB","mitoticRateEB"),names(clinU))}
    for (k in 1:length(k1)) {
        cat("\n\n============= ",names(clinU)[k1[k]],"\n")
        #print(table(is.na(clinU[,k2[k]])))
        #print(table(clinU[,k1[k]],exclude=NULL))
        #print(table(clinU[,k2[k]],exclude=NULL))
        #print(table(clinU[,k1[k]],clinU[,k2[k]],exclude=NULL))
        #j=which(is.na(clinU[,k1[k]])!=is.na(clinU[,k2[k]]) | (clinU[,k1[k]]=="" & clinU[,k2[k]]!="") | (clinU[,k1[k]]!="" & clinU[,k2[k]]==""))
        jj=63:nrow(clinU)
        j=which(clinU[jj,k2[k]]=="")
        #clinU[jj[j],k2[k]]=clinU[jj[j],k1[k]]
        jj=1:nrow(clinU)
        jj=63:nrow(clinU)
        j=jj[which(clinU[jj,k2[k]]=="")]
        #j=jj
        if (length(j)!=0) print(table(notEB=clinU[j,k1[k]],EB=clinU[j,k2[k]],exclude=NULL))
    }
}
colId=sort(names(clinU)[grep("diff|cellSize|ki67|grade|primarySite|testUCSF500orFMI",names(clinU))])
colId=c("diffEB","cellSizeEB","ki67EB","grade","primarySiteEB","testUCSF500orFMI")
colId=c("diffNJ","cellSizeEB","ki67EB","grade","primarySiteEB","testUCSF500orFMI")
for (k in colId) {
    cat("\n\n============= ",k,"\n")
    print(names(table(tolower(clinU[,k]),exclude=NULL)))
}

if (datVerFlag%in%c("_20171019")) {
    clinU$diffEB=tolower(as.character(clinU$diffEB))
    clinU$diffEB[which(clinU$diffEB%in%c("0","poorly-differentiated"))]="poor"
    clinU$diffEB[which(clinU$diffEB%in%c("1","well-differentiated"))]="well"
    clinU$diffEB[which(clinU$diffEB%in%c("2"))]="nr"
    clinU$diffEB[which(clinU$diffEB%in%c(""))]=NA
} else {
    clinU$diffNJ=tolower(as.character(clinU$diffNJ))
    clinU$diffNJ[which(clinU$diffNJ%in%c(""))]=NA
}
clinU$cellSizeEB=tolower(as.character(clinU$cellSizeEB))
clinU$cellSizeEB[which(clinU$cellSizeEB%in%c("0","smallcell"))]="small"
clinU$cellSizeEB[which(clinU$cellSizeEB%in%c("1","largecell"))]="large"
clinU$cellSizeEB[which(clinU$cellSizeEB%in%c("large/manec"))]="large-manec"
clinU$cellSizeEB[which(clinU$cellSizeEB%in%c("small & large"))]="small-large"
clinU$cellSizeEB[which(clinU$cellSizeEB=="3")]="nr"
clinU$cellSizeEB[which(clinU$cellSizeEB%in%c("2",""))]=NA
clinU$ki67EB=tolower(as.character(clinU$ki67EB))
clinU$ki67EB=paste(sub("%","",clinU$ki67EB),"%",sep="")
clinU$ki67EB[which(clinU$ki67EB%in%c("NA%","%"))]=NA
clinU$ki67EB[which(clinU$ki67EB%in%c("nr%"))]="NR"
clinU$grade[which(clinU$grade=="3-Feb")]="2/3"
clinU$grade=paste("grade",clinU$grade,sep="")
clinU$grade[which(clinU$grade=="grade0")]="NR"
clinU$grade[which(clinU$grade%in%c("gradeNA","grade"))]=NA
clinU$primarySiteEB=as.character(clinU$primarySiteEB)
clinU$primarySiteEB[which(clinU$primarySiteEB%in%c("OTHER(lung)"))]="LUNG"
clinU$primarySiteEB=sapply(as.character(clinU$primarySiteEB),function(x) {
    y=strsplit(strsplit(toupper(x),"(",fixed=T)[[1]][1],"-")[[1]][1]
    y
},USE.NAMES=F)
clinU$primarySiteEB[which(clinU$primarySiteEB%in%c("0","COLON","ANUS"))]="COLORECTAL"
clinU$primarySiteEB[which(clinU$primarySiteEB=="1")]="PANCREAS"
clinU$primarySiteEB[which(clinU$primarySiteEB%in%c("2","OTHERGI","STOMACH"))]="OTHER GI"
clinU$primarySiteEB[which(clinU$primarySiteEB%in%c("3","OTHER","BLADDER"))]="OTHER non-GI"
clinU$primarySiteEB[which(clinU$primarySiteEB=="4")]="UNKNOWN"
clinU$primarySiteEB=tolower(clinU$primarySiteEB)
clinU$testUCSF500orFMI=as.character(clinU$testUCSF500orFMI)
clinU$testUCSF500orFMI[which(clinU$testUCSF500orFMI=="0")]="UCSF500"
clinU$testUCSF500orFMI[which(clinU$testUCSF500orFMI=="1")]="FMI"
clinU$testUCSF500orFMI[which(clinU$testUCSF500orFMI=="2")]="Other"

## Exclude lung cancer primaries
samExclId=c("HGO-031","HGO-028","HGO-107","HGO-023")
if (F) {
    samId=c("HGO-031","HGO-028","HGO-107","HGO-023")
    clinU[clinU$caseId%in%samId,c("caseId","primarySite","primarySiteEB")]
    clinU=clinU[which(!clinU$caseId%in%samId),]
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


if (datVerFlag%in%c("_20171019")) {
    colId=which(names(clinU)%in%c("TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","CDKN2A","CDKN2B","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFbR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM6A","FGFR1","VHL","SMARCA4","CDK4_MDM2_FRS2","NFE2L2","MYC","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB2","PTEN","BCORL1","CHD4","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3"))
    colId=which(names(clinU)%in%c("TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","CDKN1B","CDKN2A","CDKN2B","CDKN2C","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFBR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM5A","KDM6A","FGFR1","FGFR3","VHL","SMARCA4","CDK4_MDM2_FRS2","NFE2L2","MYC","ERBB2","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB22","PTEN","BCORL1","CHD4","IDH1","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3","MYCN","BRCA2","FH","NBN","ARID2","AKT2","AXL","FGF19","FGF23","FGF3","FGF4","FGF6","AR","DNMT3A","AURKA","GNAS","ARFRP1","FLT3","LZTR1","MAP2K4","ZNF217","NF1","FANCA","RET","MCL1","BRCA1","EGFR","JAK3","TERT","FANCC","PTPRD","STK11","RNF43","BRIP1","BRIP2","PTPN11"))
    colId=which(names(clinU)%in%c("TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","CDKN1B","CDKN2A","CDKN2B","CDKN2C","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFBR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM5A","KDM6A","FGFR1","FGFR3","VHL","SMARCA4","CDK4_MDM2_FRS2","NFE2L2","MYC","ERBB2","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB22","PTEN","BCORL1","CHD4","IDH1","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3","MYCN","BRCA2","FH","NBN","ARID2","AKT2","AXL","FGF19","FGF23","FGF3","FGF4","FGF6","AR","DNMT3A","AURKA","GNAS","ARFRP1","FLT3","LZTR1","MAP2K4","ZNF217","NF1","FANCA","RET","MCL1","BRCA1","EGFR","JAK3","TERT","FANCC","PTPRB","PTPRD","STK11","RNF43","PTPN11","BRIP1","CREBBP","AXIN2","MGA","PBRM1","SDHD","STAG2","HRAS","ASXL1","BRD4","FAM123B"))
} else {
    colId=which(names(clinU)%in%c("TP53","RB1","MEN1","ATRX","DAXX","ARID1A","MLL2","SETD2","APC","CDKN1B","CDKN2A","CDKN2B","CDKN2C","POT1","KRAS","MYCL1","SMAD4","CCNE1","RICTOR","FAM58A","MLH1","CTNNB1","ACVR2A","TGFBR2","PIK3CA","NRAS","POLE","FBXW7","PDPK1","TSC1","TSC2","MUTYH","USP9X","BRAF","SPTA1","ERCC2","BLM","PTCH1","SYNE1","NOTCH1","NOTCH2","KMT2D","ARID1A2","CIC","LRP1B","GNAQ","KIT","PRKDC","CCND1","CCND2","RBM10","CYLD","NFKBIE","SRSF2","MET","PPP2R1A","KDM5A","KDM6A","FGFR1","FGFR3","VHL","SMARCA4","CDK4_MDM2_FRS2","NFE2L2","MYC","ERBB2","ERBB3","BCL2","FAT1","SOX9","BCOR","ATM","ERBB22","PTEN","BCORL1","CHD4","IDH1","CTCF","CTNNA1","MSH2","NOTCH3","SMAD3","TBX3","MYCN","BRCA2","FH","NBN","ARID2","AKT2","AXL","FGF19","FGF23","FGF3","FGF4","FGF6","AR","DNMT3A","AURKA","GNAS","ARFRP1","FLT3","LZTR1","MAP2K4","ZNF217","NF1","FANCA","RET","MCL1","BRCA1","EGFR","JAK3","TERT","FANCC","PTPRB","PTPRD","STK11","RNF43","PTPN11","BRIP1","CREBBP","AXIN2","MGA","PBRM1","SDHD","STAG2","HRAS","ASXL1","BRD4","FAM123B","MYB","KDM5C","CDK4","MDM2","FRS2","FGFR2","CDK12","EMSY","SOX2","ELF3","BAP1","CBFB","DCC","FUPB1","EP300","NTRK1","SMAD2","TSHR","TET2","FGF10","HNF1A","JAK2","SMARCB1","MSH6","MSH7","GATA2","FUBP1","TCF7L2"))
}
kk=c()
for (k in colId) {
    x=mean(is.na(clinU[,k]))
    if (x>=1) kk=c(kk,k)
}
table(colId%in%kk)

if (F) {
    x=names(table(as.matrix(clinU[,colId]))); x=x[x!=""]; x=c("alterationType",x)
    write.table(x,file=paste("alterationType.txt",sep=""),append=F,col.names=F,row.names=F,sep="\t",quote=F)
}

if (datVerFlag=="_20180929") {
} else {
    datGPU=t(apply(as.matrix(clinU[,colId]),c(1,2),function(x) {y=as.character(x); ifelse(is.na(y),0,as.integer(!y%in%c("","0")))}))
    colnames(datGPU)=paste("anyAlt_",clinU$id,sep="")
    print("sort(unique(c(as.matrix(clinU[,colId]))),na.last=T)")
    print(sort(unique(c(as.matrix(clinU[,colId]))),na.last=T))
    print("table(c(datGPU),c(as.matrix(clinU[,colId]))!='',exclude=NULL)")
    print(table(c(datGPU),c(as.matrix(clinU[,colId]))!="",exclude=NULL))
}
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

if (datVerFlag%in%c("_20171019")) {colId=c("diffEB","ki67EB2","cellSizeEB")
} else {colId=c("diffNJ","ki67EB2","cellSizeEB")}
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

if (datVerFlag=="_20180929") {
    save(clinU,datGPU,samExclId,fName1,datVerFlag,altInfo,file="tmp_ucsf500.RData")
} else {
    save(clinU,datGPU,samExclId,fName1,datVerFlag,file="tmp_ucsf500.RData")
}

table(clinU$cellSizeEB,clinU$grade)
table(altInfo$altType)

## -------------------------------------
library(coin)

pThres=0.05
samSize=5

dat=as.data.frame(t(datGPU),stringsAsFactors=F)
for (k in 1:ncol(dat)) {
    dat[,k]=as.character(dat[,k])
    dat[which(dat[,k]=="0"),k]="no mutation"
    dat[which(dat[,k]=="1"),k]="mutation"
}
dat=cbind(clinU,dat)

if (datVerFlag%in%c("_20171019")) {
    varList1=c("diffEB","cellSizeEB"); varName1=varList1
    varList1=c("primarySiteEB","diffEB","cellSizeEB","TP53","RB1"); varName1=varList1
} else {
    varList1=c("primarySiteEB","diffNJ","cellSizeEB","TP53","RB1"); varName1=varList1
}
varList2=c("ki67EB2"); varName2=c("ki67EB")
tbl=NULL
pdf("ki67_boxplot_ucsf500.pdf")
for (subsetFlag in paste(rep(c("","_grade3"),each=2),c("_ucsf500","_ucsf500Fmi"),sep="")) {
    #par(mfrow=c(2,1))
    lim=c(0,100)
    for (vId2 in 1:length(varList2)) {
        #j=1:nrow(dat)
        j=which(!is.na(dat[,varList2[vId2]]))
        heading=paste("All samples",sep="")
        subset2Flag=strsplit(subsetFlag,"_")[[1]]
        subset2Flag=paste("_",subset2Flag[(length(subset2Flag)-1):length(subset2Flag)],sep="")
        subset2Flag[subset2Flag=="_"]=""
        names(subset2Flag)=c("grade","test")
        switch(subset2Flag["grade"],
        "_grade3"={
            j=j[which(dat$grade[j]%in%c("grade3"))]
            heading=paste("Grade 3",sep="")
        }
        )
        switch(subset2Flag["test"],
        "_ucsf500"={
            j=j[which(dat$testUCSF500orFMI[j]%in%c("UCSF500"))]
            heading=paste("UCSF 500: ",heading,sep="")
        },
        "_ucsf500Fmi"={
            j=j[which(dat$testUCSF500orFMI[j]%in%c("UCSF500","FMI"))]
            heading=paste("UCSF 500, FMI: ",heading,sep="")
        }
        )
        for (vId1 in 1:length(varList1)) {
            x=dat[j,varList1[vId1]]
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
            } else if (varList1[vId1]%in%c("TP53","RB1")) {
                x[which(x=="no mutation")]="0:No mutation"
                x[which(x=="mutation")]="1:Mutation"
                ttl=sapply(sort(unique(x)),function(x) strsplit(x,":")[[1]][2],USE.NAMES=F)
            }
            ttl=paste(ttl," (",table(x),")",sep="")
            grpUniq=sort(unique(x))
            ttl2=substr(ttl,1,1)
            if (varList1[vId1]=="primarySiteEB") {
                ttl2[which(grpUniq=="OTHER GI")]="OG"
                ttl2[which(grpUniq=="OTHER non-GI")]="ON"
            }
            pv=rep(NA,length(grpUniq))
            heading2="pv: "
            k=1
            k2=0
            for (grpId1 in 1:(length(grpUniq)-1)) {
                for (grpId2 in (grpId1+1):length(grpUniq)) {
                    jj=which(x%in%grpUniq[c(grpId1,grpId2)])
                    if (all(table(x[jj])>=samSize)) {
                        pv[k]=pvalue(wilcox_test(dat[j[jj],varList2[vId2]]~as.factor(x[jj]),distribution="exact"))
                    }
                    if (!is.na(pv[k])) {
                        if (k==6) heading2=paste(heading2,"\n",sep="")
                        heading2=paste(heading2,signif(pv[k],2)," (",ttl2[grpId1]," vs ",ttl2[grpId2],"), ",sep="")
                        k2=k2+1
                    }
                    #tbl=rbind(tbl,c(subsetFlag,varList2[vId2],paste(varList1[vId1]," (",grpUniq[grpId1]," vs ",grpUniq[grpId2],") ",sep=""),signif(pv[k],2)))
                    tbl=rbind(tbl,c(subsetFlag,varList2[vId2],paste(varList1[vId1]," (",ttl[grpId1]," vs ",ttl[grpId2],") ",sep=""),signif(pv[k],2)))
                    k=k+1
                }
            }
            heading2=substr(heading2,1,nchar(heading2)-2)
            if (heading2=="pv") heading2=""
            boxplot(dat[j,varList2[vId2]]~x,names=ttl,main=paste(heading,"\n",heading2,sep=""),ylim=lim,xlab=varName1[vId1],ylab=varName2[vId2],las=0,cex.main=.7,cex.axis=.7)
        }
    }
}
dev.off()
tbl=as.data.frame(tbl,stringsAsFactors=F)
names(tbl)=c("subset","variable1","variable2","pv")
tbl$subset=sub("_","",tbl$subset)
tbl$variable1[which(tbl$variable1=="ki67EB2")]="ki67EB"
for (k in c("pv")) tbl[,k]=as.numeric(tbl[,k])
k=which(!is.na(tbl$pv))
tbl$pv[k]=paste(tbl$pv[k],ifelse(tbl$pv[k]<pThres,"*",""))
write.table(tbl,file=paste("ki67_association_ucsf500.txt",sep=""),append=F,col.names=T,row.names=F,sep="\t",quote=F)


## -------------------------------------

library(survival)

time=clinU$survSinceStageIV; event=clinU$dead

if (datVerFlag%in%c("_20171019")) {varList=c("diffEB","cellSizeEB")
} else {varList=c("diffNJ","cellSizeEB")}
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
# Are CDKN2A & RB1 mutually exclusive?
# Seen this in FMI data. Test in UCSF 500 data

candGene=c("RB1","CDKN2A")
i=which(rownames(datGPU)%in%candGene)
j=which(clinU$testUCSF500orFMI%in%c("UCSF500"))
dat=datGPU[i,j]

for (i1 in 1:(nrow(dat)-1)) {
    for (i2 in (i1+1):nrow(dat)) {
        x=table(dat[i1,],dat[i2,],dnn=list(rownames(dat)[i1],rownames(dat)[i2]))
        print(x)
        print(fisher.test(x))
    }
}
"
    CDKN2A
RB1  0  1
0   29  5
1    5  0

Fisher's Exact Test for Count Data
data:  x
p-value = 1
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
0.000000 8.559462
sample estimates:
odds ratio
0
"


## -------------------------------------
## NOT USED

k=match(names(clinU1),names(clinU2)); k1=which(!is.na(k)); k2=k[k1]
j=match(clinU1[,1],clinU2[,2]); j1=which(!is.na(j)); j2=j[j1]
for (k in 1:length(k1)) {
    if (any(clinU1[j1,k1[k]]!=clinU2[j2,k2[k]],na.rm=T) | any(is.na(clinU1[j1,k1[k]])!=is.na(clinU2[j2,k2[k]]))) print(k)
}

k=match(names(clinU1),names(clinU2)); k1=which(!is.na(k)); k2=k[k1]
j=match()
for (k in 1:length(k1)) {
    cat("\n\n================ ",names(clinU1)[k1[k]],"\n",sep="")
    cat("1: \n")
    print(table(clinU1[,k1[k]],exclude=NULL))
    cat("2: \n")
    print(table(clinU2[,k2[k]],exclude=NULL))
}

k=match(rownames(datGPU1),rownames(datGPU2)); k1=which(!is.na(k)); k2=k[k1]
for (k in 1:length(k1)) {
    cat("\n\n================ ",rownames(datGPU1)[k1[k]],"\n",sep="")
    cat("1: \n")
    print(table(datGPU1[k1[k],],exclude=NULL))
    cat("2: \n")
    print(table(datGPU2[k2[k],],exclude=NULL))
}

if (datVerFlag%in%c("_20171019")) {colId=c("diffEB","ki67EB2","cellSizeEB","grade","primarySiteEB")
} else {colId=c("diffNJ","ki67EB2","cellSizeEB","grade","primarySiteEB")}
for (k in colId) {
    cat("\n\n================\n",sep="")
    for (gId in 1:nrow(datGPU)) {
        j=which(!is.na(datGPU[gId,]) & !is.na(clinU[,k]))
        if (sum(!duplicated(clinU[j,k]))<11) {
            print(table(datGPU[gId,j],clinU[j,k],dnn=c(rownames(datGPU)[gId],k)))
        } else {
            cat("Cases having ",k," (continuous) info:\n",sep="")
            print(table(datGPU[gId,j],dnn=c(rownames(datGPU)[gId])))
        }
        cat("\n\n-------------\n",sep="")
    }
}



"
================
diffEB
TP53 nr poor well
0  3   18   16
1  1   21    1


-------------
diffEB
RB1 nr poor well
0  2   29   17
1  2   10    0


-------------


================
Cases having ki67EB2 (continuous) info:
TP53
0  1
33 20


-------------
Cases having ki67EB2 (continuous) info:
RB1
0  1
45  8


-------------


================
cellSizeEB
TP53 large nr small
0    12  4     5
1    10  5     6


-------------
cellSizeEB
RB1 large nr small
0    18  6     6
1     4  3     5


-------------


================
grade
TP53  1  2  3
0  2 14 21
1  0  2 18


-------------
grade
RB1  1  2  3
0  2 16 28
1  0  0 11


-------------


================
primarySiteEB
TP53  4 COLORECTAL OTHER PANCREAS UNKNOWN
0 11          7     2       15       2
1  3          8     4        5       3


-------------
primarySiteEB
RB1  4 COLORECTAL OTHER PANCREAS UNKNOWN
0 11         12     4       17       4
1  3          3     2        3       1


-------------
"
