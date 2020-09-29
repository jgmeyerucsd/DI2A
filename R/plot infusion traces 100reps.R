
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("mzR")
BiocManager::install("MSnbase")

### plot infusion traces and XICs

library(mzR)
#library(msdata)
library(MSnbase)
library(viridis)

setwd("C:/Users/jesse/Documents/MAGIC/")
files <- list.files(pattern="mzXML")

files
## Read the data as an MSnExp
msd <- readMSData(files, msLevel = 2)
msd
## Extract the total ion chromatogram for each file:
tic <- chromatogram(msd, msLevel=2)



fileNames(msd)
## Extract the TIC for the second file:
tic[1, 2]
#> Object of class: Chromatogram
#> Intensity values aggregated using: sum 
#> length of object: 198
#> from file: 2
#> mz range: [95.51765, 1005.043]
#> rt range: [0.486, 66.7818]
#> MS level: 1
dev.off()
svg("P:/JGM_DI2A/Manuscript/revision/figures/infusion_MS1_traces_wider.svg",height = 6, width=12 )
par(cex=1.5)

## Plot the TIC for the first file
plot(rtime(tic[1, 1]), intensity(tic[1, 1]), type = "l",
     xlab = "rtime", ylab = "intensity", main = "TIC", ylim=c(0,2e7), lwd=2, col="darkgrey")


COLORS<-viridis(8)
COLORS<-c(COLORS[2:7],"black")
n=1

for(i in file_indexes){
  print(fileNames(msd)[i])
  print(COLORS[n])
  lines(rtime(tic[1, i]), intensity(tic[1, i]), type = "l",
        xlab = "rtime", ylab = "intensity", main = "TIC", col=COLORS[n], lwd=3)
  n=n+1
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

cvtext <- substr(substrRight(fileNames(msd)[file_indexes], 12), 1, 3)
COLORS
cvtext<- paste("CV=", cvtext)
cvtext[7]<- "no FAIMS"
cvtext[8] <-"blank"
#?svg
legend(5, 2.5e9, legend=cvtext, fill = c(COLORS, "darkgrey"))

dev.off()



### now generate random XICs from the noFAIMS data
## masses: 710.3770, 660.8396
z
xic_masses <-c(710.3770, 660.8396, 872.4123)
xic_matrix <- matrix(c(xic_masses-xic_masses*1e-5, xic_masses+xic_masses*1e-5), byrow = F, nrow = 3)
xic_matrix
xics <- chromatogram(msd, mz = xic_matrix)
fileNames(msd)
xics

### plot four panel

dev.off()
svg("P:/JGM_DI2A/Manuscript/revision/figures/infusion_MS1_xics.svg",height = 6, width=8 )
#par(mfcol=c(4,1))
maxintensity_TIC <- max(intensity(tic[1,9]))
plot(rtime(tic[1, 9]), (intensity(tic[1, 9])/maxintensity_TIC)+3, type = "l",
     xlab = "retention time (seconds)", ylab = "intensity", main = "TIC", ylim=c(0,4), lwd=3, col="black")
maxint_xic1 <- max(na.omit(intensity(xics[1,9])))
lines(rtime(xics[1,9]), (intensity(xics[1,9])/maxint_xic1)+2, lwd=3,type = "l",
     xlab = "retention time (seconds)", ylab = "intensity", main = paste("XIC, m/z", xic_masses[1]), ylim=c(0,1), col="darkgrey")

## plot xic 2
maxint_xic2 <- max(na.omit(intensity(xics[2,9])))

lines(rtime(xics[2,9]), (intensity(xics[2,9])/maxint_xic2)+1, lwd=3,type = "l",
     xlab = "retention time (seconds)", ylab = "intensity", main = paste("XIC, m/z", xic_masses[2]), ylim=c(0,1), col="darkgrey")

## plot xic 3
maxint_xic3 <- max(na.omit(intensity(xics[3,9])))
lines(rtime(xics[3,9]), (intensity(xics[3,9])/maxint_xic3)+0, lwd=3,type = "l",
     xlab = "retention time (seconds)", ylab = "intensity", main = paste("XIC, m/z", xic_masses[3]), ylim=c(0,1), col="darkgrey")
dev.off()



### angiotensin trace

dev.off()
svg("P:/JGM_DI2A/Manuscript/revision/figures/infusion_angio_trace.svg",height = 4, width=8 )
par(cex=1.5)

## Plot the TIC for the first file
plot(rtime(tica[1, 1]), intensity(tica[1, 1]), type = "l",
     xlab = "retention time (seconds)", ylab = "intensity", main = "TIC", lwd=3, col="black")
dev.off()



plot()
