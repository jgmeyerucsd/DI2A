### This script takes the unfiltered list of peptide IDs from the python script:
# "DI2A_make_quant_targets_file_for_R.ipynb"

# determines the protein-level FDR
# and outputs a list of quantitative assays to use for computing the quantitative SILAC ratios
# this output is used with the python script:
# "DI2A_targeted_protein_quantification.ipynb"
# these lists are used to generate the targeted data collection method on the Orbitrap Fusion Lumos


require(seqinr)
require(Biostrings)
df <-read.delim("P:/JGM_DI2A/Python/outputs/peptide_target_table_unfiltered.txt", sep="\t", stringsAsFactors = F)

head(df) 

## add column for stripped sequences to use for searching for protein hits
ppepcln<-gsub(x=df$Peptide, pattern="+[0-9]*.[0-9]", replacement = "")
df <- cbind(df, ppepcln)

sequences <- as.character(df$ppepcln)
length(sequences)

fasta= "P:/JGM_DI2A/MSPLIT-DIAv1.0/2019-03-14-td-UP000005640.fasta"
findpep=function(X, fas){names(unlist(vmatchPattern(subject=fas, pattern=X, fixed=TRUE)))}
fas_ss<-readAAStringSet(filepath=fasta, format="fasta",
                        nrec=-1L, skip=0L, seek.first.rec=FALSE,
                        use.names=TRUE, with.qualities=FALSE)

op <- lapply(FUN=findpep, X=sequences, fas=fas_ss)

nprots<-c()
for( x in op){
  print(length(x))
  nprots <-c(nprots, length(x))
}

decoy_positions = which(df$Name=="DECOY_null")

### combine protein names, add to the table as new column
combineProteins = function(proteinstostring = op[[10]]){
  return(paste(proteinstostring, collapse=", "))
}

combinedproteinslist <- lapply(FUN=combineProteins, X = op)

#unlist(combinedproteinslist)
#combinedproteinslist
df_wprot <- cbind(df, unlist(combinedproteinslist), nprots)
head(df_wprot)
nrow(df_wprot)
#df_wprot[1810,]
# remove the lines with 2 or more protein hits
df_wprot_onlysingles <- df_wprot[which(nprots==0 | nprots==1),]

nrow(df_wprot_onlysingles)



df_wprot_onlysingles$`unlist(combinedproteinslist)`<- as.character(df_wprot_onlysingles$`unlist(combinedproteinslist)`)
for( x in which(df_wprot_onlysingles$Name=="DECOY_null") ){
  df_wprot_onlysingles$`unlist(combinedproteinslist)`[x]<- paste("Decoy_null", x, sep="",collapse="")
}



which(df_wprot_onlysingles$Name=="DECOY_null")

df_wprot_onlysingles[1810,]

df_wprot_onlysingles$cosine


## get first occurance of the protein ##
onlyfirstprot  <- df_wprot_onlysingles[match(unique(df_wprot_onlysingles$`unlist(combinedproteinslist)`), df_wprot_onlysingles$`unlist(combinedproteinslist)`),]
head(onlyfirstprot)
onlyfirstprot

# positions of decoys in the new df
protdecoyindex<-which(onlyfirstprot$Name=="DECOY_null")

## loop to determine FDR

FDR = 0.01
FDRs<-c()

#while(FDR<0.01)
lastdecoy<-c()
for( i in 1:length(protdecoyindex)){
  FDRs <- c(FDRs, i/protdecoyindex[i])
  if(FDR>fdr_level){
    lastdecoy<-c(protdecoyindex[i])
  }
}


FDRs
lastdecoy
protdecoyindex


#write table of proteins for downstream analysis
nrow(onlyfirstprot)
unique(onlyfirstprot$`unlist(combinedproteinslist)`)
length(protdecoyindex)
13/565

### all proteins would be FDR < 2.5%

########    alternatively, get the most intense peptides instead of best scoring peptides  ############
?order
intsort <- df_wprot_onlysingles[order(df_wprot_onlysingles$IonCount, decreasing=T ),]

only1stint  <- intsort[match(unique(intsort $`unlist(combinedproteinslist)`), intsort $`unlist(combinedproteinslist)`),]
head(only1stint)
only1stint$IonCount
hist(only1stint$IonCount, breaks=c(seq(0,1e7, 100000)))

#remove the lines that are decoy
finaltargets <- onlyfirstprot[-protdecoyindex,]

# same for only1stintense
mostintdecoyindex<-which(only1stint$Name=="DECOY_null")

finalmostint <- only1stint[-mostintdecoyindex,]

nrow(finaltargets)
nrow(finalmostint)
### distribution of total ion count for final targets

par(mfcol=c(2,1))
hist(finaltargets$IonCount, breaks=c(seq(0,8e6, 10000)))
hist(only1stint$IonCount, breaks=c(seq(0,8e6, 10000)))

### before writing table, separate into several separate dataframes where each contains one set of CVs

cv_factors<-as.factor(finaltargets$CV)
hist(finaltargets$CV)

cv_levels<- levels(as.factor(finaltargets$CV))


#### making loop to separate into list of dfs
cvslist = list()

for(x in cv_levels){
  cvslist[[x]] = finaltargets[cv_factors==x,]
}

cvslist[[1]]


#######   write table of all the assays based on best scoring:  ##########

cv_levels<-c("-30", "-40", "-50", "-60", "-70", "-80")
scan_counter = 1
cv_counter = 0

for(x in cv_levels){
  print(x)
  df<-data.frame(Peptide=cvslist[[x]]$Peptide[1], 
                    prec_light_mz= cvslist[[x]]$Mz.1[1],
                    prec_heavy_mz = cvslist[[x]]$Heavymz[1],
                    z = cvslist[[x]]$z.1[1], 
                    scan=scan_counter,
                    ylight = cvslist[[x]]$ylight_quant_frags[1],
                    yheavy = cvslist[[x]]$yheavy_quant_frags[1],
                    yordinals = cvslist[[x]]$fragment_ordinals[1], 
                    CV=x)
  scan_counter=scan_counter+1
  for( i in 2:nrow(cvslist[[x]])){
    print(i)
    print(cvslist[[x]]$Peptide[i])
    df<-rbind(df, data.frame(Peptide=cvslist[[x]]$Peptide[i], 
                             prec_light_mz= cvslist[[x]]$Mz.1[i],
                             prec_heavy_mz = cvslist[[x]]$Heavymz[i],
                             z = cvslist[[x]]$z.1[i], 
                             scan=scan_counter,
                             ylight = cvslist[[x]]$ylight_quant_frags[i],
                             yheavy = cvslist[[x]]$yheavy_quant_frags[i],
                             yordinals = cvslist[[x]]$fragment_ordinals[i], 
                             CV=x))
      scan_counter=scan_counter+1
  }
  cv_counter = cv_counter + 1
  if(cv_counter==1)totaldf<-df
  if(cv_counter>=2)totaldf<-rbind(totaldf, df)
}
head(totaldf)
nrow(totaldf)

write.table(totaldf, file="P:/JGM_DI2A/R/outputs/bestscore_quantmzlist.txt", row.names = F, quote=F, col.names = TRUE, sep="\t")



########################################################################################
########################################################################################
########################################################################################
########################################################################################
############    write new tables with most intense pep per protein


### before writing table, separate into several separate dataframes where each contains one set of CVs

cv_factors<-as.factor(finalmostint$CV)
hist(finalmostint$CV)

cv_levels<- levels(as.factor(finalmostint$CV))

#### making loop to separate into list of dfs
cvslist = list()
totalrows = 0
for(x in cv_levels){
  cvslist[[x]] = finalmostint[cv_factors==x,]
  print(nrow(cvslist[[x]]))
  totalrows = totalrows+nrow(cvslist[[x]])
}
totalrows

cvslist[[1]]

44+179+130+84+63+52
#######   write table of all the assays:  ##########

cv_levels<-c("-30", "-40", "-50", "-60", "-70", "-80")
scan_counter = 1
cv_counter = 0

for(x in cv_levels){
  print(x)
  df<-data.frame(Peptide=cvslist[[x]]$Peptide[1], 
                 prec_light_mz= cvslist[[x]]$Mz.1[1],
                 prec_heavy_mz = cvslist[[x]]$Heavymz[1],
                 z = cvslist[[x]]$z.1[1], 
                 scan=scan_counter,
                 ylight = cvslist[[x]]$ylight_quant_frags[1],
                 yheavy = cvslist[[x]]$yheavy_quant_frags[1],
                 yordinals = cvslist[[x]]$fragment_ordinals[1], 
                 CV=x)
  scan_counter=scan_counter+1
  for( i in 2:nrow(cvslist[[x]])){
    #print(i)
    #print(cvslist[[x]]$Peptide[i])
    df<-rbind(df, data.frame(Peptide=cvslist[[x]]$Peptide[i], 
                             prec_light_mz= cvslist[[x]]$Mz.1[i],
                             prec_heavy_mz = cvslist[[x]]$Heavymz[i],
                             z = cvslist[[x]]$z.1[i], 
                             scan=scan_counter,
                             ylight = cvslist[[x]]$ylight_quant_frags[i],
                             yheavy = cvslist[[x]]$yheavy_quant_frags[i],
                             yordinals = cvslist[[x]]$fragment_ordinals[i], 
                             CV=x))
    scan_counter=scan_counter+1
  }
  cv_counter = cv_counter + 1
  if(cv_counter==1)totaldf<-df
  if(cv_counter>=2)totaldf<-rbind(totaldf, df)
}

scan_counter


head(totaldf)
nrow(totaldf)
#totaldf[555,]
write.table(totaldf, file="P:/JGM_DI2A/R/outputs/mostintense_quantmzlist.txt",
            row.names = F, quote=F, col.names = TRUE, sep="\t")



