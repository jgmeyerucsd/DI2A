names(unlist(vmatchPattern(subject=fas, pattern=pepvec.cleaned[230], fixed=TRUE)))
pepvec.cleaned[230]


f<-list.files(pattern="txt")
f<-list.files(pattern="10p")

getwd()
setwd("P:/JGM_DI2A/20190405/FAIMS")
setwd("P:/JGM_DI2A/20190405/noFAIMS")


f<-list.files(pattern="txt")

f

for(x in f){
  print(substr(x,start=21,stop=22))
  
  peplvlfdr(msplitresults = x,fdrlevel = 0.01)
  
}

f

for(x in f){
  print(unlist(strsplit(x,split="_"))[4])
  peplvlfdr(msplitresults = x,fdrlevel = 0.01)
  
}
x
