#Check TSR data for Paralog coverage

tsrset = read.table("/N/dc2/scratch/tlicknac/test_STRIPE_Dec/tsr/TSRsetMerged-1.txt", header=T)

tsrset = tsrset[-which(is.na(tsrset$featureID)),]
tsrsetP = gsub(".G0", ".P0", tsrset$featureID)
tsrset$featureID = tsrsetP


#Recent WGD1 Paralogs
recentWGD = read.table("~/Paramecium_POFF/all_aurelias-cut-mod.poff", header=T)
recentWGD = recentWGD$pdec
wgdPairs = strsplit(as.character(recentWGD), ",")
wgdPairs = lapply(wgdPairs, as.character)
wgdPairs = wgdPairs[-which(sapply(wgdPairs, FUN=function(X) "." %in% X ))]
wgdPairs = wgdPairs[-which(sapply(wgdPairs, FUN=function(X) "NA" %in% X ))]
wgdPairs = do.call("rbind", wgdPairs)

finalWGD1df = data.frame(matrix(ncol=6))
colnames(finalWGD1df) = c("ParalogA", "ntagsA", "msiA", "ParalogB", "ntagsB", "msiB")

tmpdf =  data.frame(matrix(ncol=6))
colnames(tmpdf) = c("ParalogA", "ntagsA", "msiA", "ParalogB", "ntagsB", "msiB")

multHits = c()
first = T

for(q in 1:nrow(tsrset)){        
  hitPair =  wgdPairs[max( c( which(wgdPairs[,2] == tsrset$featureID[q]), which(wgdPairs[,1] == tsrset$featureID[q]) ) ),]
   
  if(all(is.na(hitPair)) == F){ 
    
    if(length(hitPair) > 0){ #----- removing misses -----
      firstRow = tsrset[which(tsrset$featureID == hitPair[1]),]
      secondRow = tsrset[which(tsrset$featureID == hitPair[2]),]
    
      if(length(firstRow$seq) > 0 & length(secondRow$seq) > 0){         #-----
        paralogs = c(hitPair[1], hitPair[2])
        TAGs = c(firstRow$nTAGs, secondRow$nTAGs)
        MSI = c(firstRow$tsrMSI, secondRow$tsrMSI)
      }
      if(length(firstRow$seq) == 0 & length(secondRow$seq) > 0){
        paralogs = c(NA, hitPair[2])
        TAGs = c(NA, secondRow$nTAGs)
        MSI = c(NA, secondRow$tsrMSI)
      }
      if(length(firstRow$seq) > 0 & length(secondRow$seq) == 0){
      paralogs = c(hitPair[1], NA)
        TAGs = c(firstRow$nTAGs, NA)
        MSI = c(firstRow$tsrMSI, NA)
      }                                                                 #-----
      
      if(first == T){ #-------------------- for the first row --------------------
        finalWGD1df$ParalogA= paralogs[1]
        finalWGD1df$ntagsA = TAGs[1]
        finalWGD1df$msiA = MSI[1]
        finalWGD1df$ParalogB= paralogs[2]
        finalWGD1df$ntagsB = TAGs[2]
        finalWGD1df$msiB= MSI[2]
      } #-------------------- for the first row --------------------
      if(first == F){ #------------------------- for the first row -------------------------
        tmpdf$ParalogA= paralogs[1]
        tmpdf$ntagsA = TAGs[1]
        tmpdf$msiA = MSI[1]
        tmpdf$ParalogB= paralogs[2]
        tmpdf$ntagsB = TAGs[2]
        tmpdf$msiB= MSI[2]
        finalWGD1df = rbind(finalWGD1df, tmpdf)
      } #------------------------- for the first row -------------------------
      first = F
    } #------ nothing after this
  }
}
row.names(finalWGD1df) = NULL
finalWGD1df = finalWGD1df[!duplicated(finalWGD1df),]  #cleanup duplciated rows

nrow(finalWGD1df)   #1803 total paralogs with default settings... 1785 with seed length = 60 ... 1767 with bwa mem... 

missingWGD = finalWGD1df[rowSums(is.na(finalWGD1df)) > 0,]
nrow(missingWGD)    #1138 are missing expression of one copy... 1130 with l=60 ... 1071 with bwa mem
bothWGD = finalWGD1df[rowSums(is.na(finalWGD1df)) == 0,]
nrow(bothWGD)       #665 are expressing both ... 655 with l=60 ... 696 with bwa mem

write.table(finalWGD1df, file = "pdec_tsr-nTAGs-MSI_wgd1-ALL-paralogs-mem.tab", sep = "\t", quote = F, row.names = F)
write.table(missingWGD, file = "pdec_tsr-nTAGs-MSI_wgd1-ONE-paralogs-mem.tab", sep = "\t", quote = F, row.names = F)
write.table(bothWGD, file = "pdec_tsr-nTAGs-MSI_wgd1-BOTH-paralogs-mem.tab", sep = "\t", quote = F, row.names = F)

summary(tsrset$nTAGs)
#Default
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#5       8      23    1027      71 5420500
#bwa mem
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#5      10      28    1005      99 4873032 

summary(tsrset$tsrMSI)
#Default
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.2500  1.0000  0.6726  1.0000  1.0000
#bwa mem
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.0800  0.2300  0.4418  1.0000  1.0000

summary(c(finalWGD1df$ntagsA, finalWGD1df$ntagsB))  
#Default
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#5.0     8.0    21.0   122.5    65.0 51676.0    1115
#bwa mem
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#5.0    10.0    27.0   168.8    96.0 46629.0    1049 

summary(c(finalWGD1df$msiA, finalWGD1df$msiB))  
#Default
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.00    0.23    1.00    0.67    1.00    1.00    1115 
#bwa mem
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#0.0000  0.0800  0.2300  0.4403  1.0000  1.0000    1049 

summary(c(bothWGD$ntagsA, bothWGD$ntagsB))  
#Default
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#5.0    13.0    36.5   153.2   104.0  9636.0
#bwa mem
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#5.0    15.0    49.0   225.8   148.2 20835.0 
summary(c(bothWGD$msiA, bothWGD$msiB)) 
#Default
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.1800  0.5600  0.5787  1.0000  1.0000
#bwa mem
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.0800  0.1900  0.3618  0.5375  1.0000

library(ggplot2)
ggplot(data=bothWGD, aes(ntagsB)) + geom_histogram()                                                                 #look at distribution of nTAGs (expression)
ggplot(data=bothWGD, aes(ntagsB)) + geom_histogram() + xlim(0,2000)
ggplot(data=bothWGD, aes(ntagsA, ntagsB)) + geom_point() + geom_smooth(method = 'lm') + stat_smooth_(method = "lm")                              #Look at relationship
ggplot(data=bothWGD, aes(ntagsA, ntagsB)) + geom_point() + geom_smooth(method = 'lm') + xlim(0,2500) + ylim(0,2500)  #throw away outliers
cor.test(bothWGD$ntagsA, bothWGD$ntagsB, method = "spearman", continuity = T)  #rho=0.28 with default... rho=0.288

ggplot(data=bothWGD, aes(msiA)) + geom_histogram()
ggplot(data=bothWGD, aes(msiA, msiB)) + geom_point() + geom_smooth(method = 'lm')
cor.test(bothWGD$msiA, bothWGD$msiB, method = "spearman", continuity = T)  #rho=0.15 with default... rho=0.1

bothWGD[bothWGD$ntagsA > 10000,]







#---------------------------------------------------------------------------------------------
#Old Paralog Comparisons

wgd = read.table("~/Paramecium_WGD3_Trees/pdecaurelia_223_annotation_v1.0.WGD.tree", header=T)
wgd$NB = NULL

lhits = list()

for(i in 1:nrow(tsrset)){
  tsrGene = tsrset$featureID[i]
  
  for(j in 1:nrow(wgd)){
    wgdrow = as.character(unlist(as.list(wgd[j,])))
    tmp = which(tsrGene == wgdrow)
    
    if(length(tmp > 0)){
      cat("Line: ", j, "\n")
      geneFamily = wgdrow[which(wgdrow != ".")]
      
      hits = tsrset[which(tsrset$featureID == geneFamily),]
      lhits = append(lhits, hits)
    }
  }
}


#------------------------
#Prepare gff of small number of paralogs from bothWGD above
#bash --> head -n12 Pdec-Paralog_TSRinfo.tab > Pdec-Paralog_TSRinfo-HEAD.tab
gff = read.table("/N/u/tlicknac/Carbonate/Paramecium_GFF/pdec-full.gff", header=F)

sample = read.table("/N/dc2/scratch/tlicknac/Data-TSR/pdec_tsr-nTAGs-MSI_wgd1-paralogs-HEAD.tab", header=T)

name = "Row"

for(z in 1:nrow(sample)){
  p1 = gsub(".P", ".G", sample[z,1] )
  p2 = gsub(".P", ".G", sample[z,4] )
  
  gffRow1 = gff[grep(p1, gff$V9),][1,]
  gffRow2 = gff[grep(p2, gff$V9),][1,]
  out = rbind(gffRow1, gffRow2)
  outname = paste(name, z, sep = "")
  write.table(out, file=outname, row.names = F, sep="\t", col.names = F, quote = F)
}



