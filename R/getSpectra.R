#source("https://bioconductor.org/biocLite.R")
#biocLite("mzR")
library(mzR)
#library(msdata)

#### gets a single spectra
###     to implement whole spectra matching process in R
####    change to get all spectra as R object

mgf.lib<-readLines(con="P:/JGM_DI2A/MSPLIT-DIAv1.0/human.faims.fixed.decoy.mgf")

get.mgf.spec=function(mgf=mgf.lib, 
                      specID="TITLE=MSfragger1.45180.45180")
{
  tmp.index<-grep(pattern=specID, mgf)[1]+5
  mgf[tmp.index]#tmp.index
  first.line<-tmp.index+1
  while(mgf[tmp.index]!="END IONS")  {tmp.index=tmp.index+1}
  last.line<-tmp.index-1
  spectra <- as.numeric(unlist(strsplit(c(mgf[first.line:last.line]),split=" ")))
  mz<-spectra[seq(from=1, to=length(spectra), by=2)]
  int<-spectra[seq(from=2, to=length(spectra), by=2)]
  #plot(mz,int,type="h",lwd=1, xlim=c(200, 1200))  
  return(data.frame(mz,int))
}


### feed in output from spectra
filter.dia=function(spec=rawspec,
                    lib_spec=secondpep,
                    tol=30,
                    tol.type="ppm"){
  ### define the tolerances
  if(tol.type=="ppm"){
    low.mz=lib_spec[,"mz"]-lib_spec[,"mz"]*(tol/1000000)
    high.mz=lib_spec[,"mz"]+lib_spec[,"mz"]*(tol/1000000)
  }
  if(tol.type=="da"){
    low.mz=lib_spec[,"mz"]-tol
    high.mz=lib_spec[,"mz"]+tol
  }
  
  ### start building filtered matched spectra
  btw<-matrix(spec[spec[,1]>=low.mz[1] & spec[,1]<=high.mz[1]],byrow=F, ncol=2)
  
  ### build output of ppm error
  ppm<-c()
  #rbind(btw,btw)
  for(i in 2:length(low.mz)){
    #print(i)
    tmpmatch<-spec[spec[,1]>=low.mz[i] & spec[,1]<=high.mz[i]]
    #tmpmatch
    btw<-rbind(btw,tmpmatch)
    if(length(tmpmatch)==0){
      btw<-rbind(btw,c(lib_spec[i,"mz"],0))
    }
    ### if matched, compute ppm error
    if(length(tmpmatch)>0 && tol.type=="ppm"){
      ppm<-c(ppm,((lib_spec[i,"mz"]-tmpmatch[1])/lib_spec[i,"mz"])*1000000)
    }
    #print(nrow(btw))
  }
  
  return(list(btw,ppm))
}

#filtered.rawspec<-filter.dia()


#peaks(ms1da)
#s18560<-spectra(ms1da, scans=18560)

df<-read.delim("P:/JGM_DI2A/SILAC/2da10ppm_1to16_n1b.txt", sep="\t", stringsAsFactors = F)

head(df)

df <- df[order(-df[,"cosine"]),]



### plot all 3 ---- pt 4
### spec lib
#mgf.lib<-readLines(con="C:/Users/jmeyer/Documents/msfragger_decoys_fixed.mgf")
### IT pt4
#ms1da<-openMSfile(filename="D:/20180816_FIA_DIA/201808aug22_JGM_FDIA2_IT_pt4_ol.mzXML")
### IT pt 8 
#ms1da<-openMSfile(filename="D:/20180816_FIA_DIA/201808aug15_JGM_FDIA2_IT_pt8_ol.mzXML")
### OT 1 Th
### 50k OT
ms1da<-openMSfile(filename="P:/JGM_DI2A/SILAC/20190411_DI2A_1to16_n1b.mzXML")

line <- 7077
head(df)
tail(df)
nrow(df)
df[2,]
df[,"Scan."][line]
rawspec<-spectra(ms1da, scans=df[,"Scan."][line])
secondpep<-get.mgf.spec(mgf=mgf.lib, specID=unlist(strsplit(df[,"Name"][line], split = " "))[1])

### IT, 0.3 Da
#f.secondpep<-filter.dia(spec=rawspec, lspec=secondpep, tol=0.3, tol.type="da")
### OT, 10ppm
f.secondpep<-filter.dia(spec=rawspec, lib_spec=secondpep, tol=10, tol.type="ppm")
#f.secondpep
#f.secondpep[[2]]
hist(f.secondpep[[2]])

ps<-f.secondpep[[1]]

#dev.off()
#par(mfcol=c(3,1),cex.lab=1.2, cex.axis=1.2, mai=c(0.5,0.5,0.5,0))
#plot(rawspec[,1],rawspec[,2],type="h",lwd=1, xlim=c(200, 1200), xlab="mz", ylab="int", main="raw.dia")
#plot(secondpep[,1],secondpep[,2],type="h",lwd=1, xlim=c(200, 1200),xlab="mz", ylab="int", main="library spec")  
#plot(f.secondpep[,1],f.secondpep[,2],type="h",lwd=1, xlim=c(200, 1200),xlab="mz", ylab="int", main="filtered raw")  
range(rawspec[,1])
# make mirrored library spectra and projected spectra
dev.off()
?tiff
tiff("P:/JGM_DI2A/R/outputs/spectra/lastpepRaw.tiff", width= 4, height=2, units="in", res=600)
par(mai=c(0.5,0.5,0.5,0.5))
#layout(matrix(1,1,1,2,2,2,2,2,2), 3, 3)
plot(rawspec[,1],rawspec[,2],type="h",lwd=1, 
     xlim=c(range(rawspec[,1])[1]-25, range(rawspec[,1])[2]+50), 
     xlab="m/z",
     ylab="int",
     main="raw.dia")
dev.off()



### normalize heights
normheights1<-secondpep[,2]/max(secondpep[,2])
normheights2<-ps[,2]/max(ps[,2])

tiff("P:/JGM_DI2A/R/outputs/spectra/lastpep_projected.tiff", width= 4, height=3, units="in", res=600)
par(mai=c(0.5,0.5,0.5,0.5))
plot(secondpep[,1],-normheights1,type="h",col="red",lwd=1, 
     xlim=c(range(rawspec[,1])[1]-25, range(rawspec[,1])[2]+50),
     ylim=c(-1,1),xlab="mz", ylab="int", main="projected(top) vs. library(bot)")  
lines(ps[,1],normheights2,type="h",col="black",lwd=1, xlim=c(200, 1300),xlab="m/z", ylab="int")  
abline(h=0)
dev.off()






library(OrgMassSpecR)
dev.off()

SpectrumSimilarity
SpectrumSimilarity(spec.top=data.frame(secondpep,normheights1), 
                   spec.bottom=data.frame(f.secondpep, normheights2), 
                   t = 0.11, b = 0.0001, 
                   top.label = NULL, bottom.label = "Library", 
                   xlim = c(200, 1200))


#### compute the cosine similarity


cos.sim <- function(ix)
{
  A = ix[,1]
  B = ix[,2]
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
} 

length(secondpep[,1])
length(f.secondpep[,1])

ix<-matrix(c(secondpep[,1],normheights1),
           c(f.secondpep[,1],normheights2), ncol=2, nrow=length(secondpep[,1]))

ix
cos.sim(ix)


### compute the angle between the 2 spectra

theta <- acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )