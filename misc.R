######################################################################
## Download from ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/

fName1="candGene_swiSnfMod"
tbl1=read.table(paste("/Users/royr/Downloads/Homo_sapiens.gene_info",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
tbl2=read.table(paste("docs/",fName1,".txt",sep=""),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)

table(toupper(tbl2$gene)%in%toupper(tbl1$Symbol))
x2=toupper(unique(tbl2$gene[!toupper(tbl2$gene)%in%toupper(tbl1$Symbol)]))

if (F) {
    ## NOT USED
    x1=sapply(toupper(tbl1$Synonyms),function(x,target) {
        #y=strsplit(x,"|",fixed=T)[[1]]
        if (x=="-") i=c() else i=grep(x,target)
        ifelse(length(i)==0,NA,paste(i,collapse=","))
    },target=toupper(x2),USE.NAMES=F)
    i=which(!is.na(x1) & !duplicated(x1))
    x21=data.frame(id1=i,id2=paste(",",x1[i],",",sep=""),stringsAsFactors=F)
    i=match(x2,toupper(tbl2$gene))
    x22=data.frame(i1=rep(NA,length(i)),id2=i,stringsAsFactors=F)
    for (i2 in 1:length(x2)) {
        i=grep(paste(",",i2,",",sep=""),x21$id2)
        if (length(i)!=0) {
            
            i1=strsplit(x21$id2,",")[[1]]
            i1=i1[2:length(i1)]
            x=c()
            for (i in 1:length(i1)) {
                x=c(x,strsplit(tbl1$Synonyms,"|",fixed=T)[[1]])
            }
            x=unique(toupper(x))
        }
    }
}

x21=paste("|",toupper(tbl1$Symbol),"|",toupper(tbl1$Synonyms),"|",sep="")
x22=paste("|",toupper(x2),"|",sep="")
tbl21=tbl2
for (i22 in 1:length(x2)) {
    i=grep(x22[i22],x21,fixed=T)
    if (length(i)!=0) {
        x=strsplit(x21[i],"|",fixed=T)[[1]]; x=x[2:length(x)]
        i23=which(toupper(tbl2$gene)==x2[i22])
        n=length(x)
        tbl22=data.frame(family=rep(tbl2$family[i23],n),gene=x,altType=rep(tbl2$altType[i23],n),alias=rep(tbl2$alias[i23],n),family2=rep(tbl2$family2[i23],n),stringsAsFactors=F)
        tbl21=rbind(tbl21,tbl22)
    }
}
tbl21=tbl21[!duplicated(tbl21$gene),]

write.table(tbl21,file=paste(fName1,"_v2.txt",sep=""),append=F,col.names=T,row.names=F,sep="\t",quote=F)

######################################################################
## Create genelist

candGene=read.table("geneList_2018104.txt",sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
load("results/ucsf500/30topGene/geneSampleId/geneSampleId_geneByDisease_30topGene_ucsf500Fmi_grade3_noAmp.RData")
i=match(candGene$gene,geneId)
table(is.na(i))
geneId=geneId[i]
save(geneId,sampleId,file="geneSampleId_geneByDisease_27topGene_ucsf500Fmi_grade3_noAmp.RData")

######################################################################
