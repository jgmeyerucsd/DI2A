#setwd("P:/JGM_DI2A/20190402_FAIMS_boudica/DI2A_conditionGrid/FAIMS/")

#f<-list.files(pattern=".txt")

### gives error if 1% FDR not reached
peplvlfdr=function(msplitresults="20190402_MCF7_FAIMS_17_1.txt", fdrlevel=0.01) {
  res.tab<-read.delim(msplitresults,stringsAsFactors = F, header=T)
  sorted.results<-res.tab[order(-res.tab[,"cosine"]),]
  #peptide.factors<-as.factor(sorted.results[,"Peptide"])
  t.first <- sorted.results[match(unique(sorted.results$Peptide), sorted.results$Peptide),]
  maxlines<-nrow(t.first)
  decoylines<-grep(t.first[,"Name"], pattern="DECOY")
  n.decoylines<-length(decoylines)
  fdr=0
  #i = 1
  lastdecoy<-c()
  ### loop through the decoy lines
  if(n.decoylines==0){lastdecoy<-NULL
  } else { 
    for(i in 1:n.decoylines){
      fdr<-i/decoylines[i]
      if(fdr>fdrlevel){
        lastdecoy<-c(lastdecoy,i)
      }
    }
  }
  
  if(is.null(lastdecoy)==TRUE){
    print(paste("not enough decoys, accept all",maxlines-n.decoylines,  "peptides @ FDR=", n.decoylines/maxlines, sep=" "))
          #print(n.decoylines/maxlines)
          #print(paste("peptides =",maxlines-n.decoylines))
          pep.output<-t.first[-decoylines,]
          return(pep.output)
  }
  
  if(is.null(lastdecoy)==FALSE & (min(lastdecoy)-1) !=0){
    i=min(lastdecoy)-1
    fdr<-i/decoylines[i]
    print(paste("fdr", round(length(decoylines[1:i])/decoylines[i], digits = 4)))
    print("score cutoff")
    print(t.first[decoylines[i],"cosine"])
    print(paste("peptide hits=", decoylines[i]-i))
    pep.output<-t.first[1:max(decoylines),]
    return(pep.output)
  }
  if(is.null(lastdecoy)==FALSE & (min(lastdecoy)-1) ==0){
    print(paste("FDR over", fdrlevel))
  }
  ### make output
  #pep.output
}
