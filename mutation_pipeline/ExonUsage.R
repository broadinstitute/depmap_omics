## reads rna junction read count file and generates exon usage ratios
# Mahmoud Ghandi 7/30/16
## inputs: 
# juncReadFN: example: 'fh_A253_SALIVARY_GLAND-Tumor.J.txt'
# exonsFN: exons: exons definition. Needs to have the following four fields: chr, ExonUsage3.pos, ExonUsage5.pos, enam (which is just an identifier for the exon). example: '/xchip/cle/Resources/ExonUsage/exons.RData'
# outFN: txt output FN  #'output.txt'
# outRobjFN: R object output FN  #'output.Rdata'
# minCov: minimum coverage needed (junctions with less than minCov will be replaced by NA --- no robust ratio estimated). example: 10
## 


readOneFile = function(infile, exons, minCov){
  
  #infile= 'fh_A253_SALIVARY_GLAND-Tumor.J.txt'
  x = read.table(infile, header=TRUE)
  #dim(x)
  #x[1:3,]
  
  schr=c(1:22,'X','Y')
  chr =match(x$chr,schr)
  fpos=x$start
  tpos=x$end
  
  #juncName=paste(paste(chr,fpos, sep='_'), tpos, paste='_')
  juncName=exons$enam
  
  incF = rep(NA, length(juncName)); 
  names(incF)= juncName; 
  excF = incF;
  incT = incF;
  excT = incF;
  
  
  ichr=22
  for(ichr in 24:1){
    ii =which(chr==ichr);
    
    ii_exons = which(exons$chr==ichr)
    
    
    datii=x[ii, ]
    fpi=datii$start
    tpi=datii$end
    npi=datii[,4] 
    
    print(paste('chr',ichr))
    for(i in ii_exons){
      fpii = exons$ExonUsage3.pos[i]
      tpii = exons$ExonUsage5.pos[i]
      
      intp = which(tpi==tpii)     # every edges from upstream that goes to tpii (aprox. expr of the exon that starts at tpii minus its expr. when it is the first exon or intron retained)
      outfp = which(fpi==fpii)    # every edges downstream that goes from fpii (aprox. expr of the exon that ends at fpii minus its expr. when it is the last exon or intron retained)
      skiptp = which((fpi<tpii)&(tpi>tpii)) # every edges that skip over tpii (aprox. isoforms in which exon starting at tpii is skipped)
      skipfp = which((fpi<fpii)&(tpi>fpii)) # every edges that skip over fpii (aprox. isoforms in which exon ending at fpii is skipped)
      
      if(length(intp)>0) {incT[i]=sum(npi[intp], na.rm=TRUE)}else{incT[i]=0}
      if(length(outfp)>0){incF[i]=sum(npi[outfp], na.rm=TRUE)}else{incF[i]=0}
      if(length(skiptp)>0){excT[i]=sum(npi[skiptp], na.rm=TRUE)}else{excT[i]=0}
      if(length(skipfp)>0){excF[i]=sum(npi[skipfp], na.rm=TRUE)}else{excF[i]=0}
      
    }
  }
  
  rat5 = incT/(incT+excT)  #ratio of usage for each exon 
  rat3 = incF/(incF+excF)  #ratio of usage for each exon 
  cnt5 = incT+excT
  cnt3 = incF+excF
  
  rat5[which(cnt5<minCov)]=NA; 
  rat3[which(cnt3<minCov)]=NA; 
  
  res = list(rat5=rat5, cnt5=cnt5, rat3=rat3, cnt3=cnt3)
  return(res)  
}

args<-commandArgs(TRUE)

juncReadFN = as.character(args[[1]])  #'fh_A253_SALIVARY_GLAND-Tumor.J.txt'
exonsFN = as.character(args[[2]])  # '/xchip/cle/Resources/ExonUsage/exons.RData'
outFN =  as.character(args[[3]])  #'output.txt'
outRobjFN = as.character(args[[4]])  #'output.Rdata'
minCov = as.numeric(args[[5]])  #10; 


load(exonsFN)
res = readOneFile(juncReadFN, exons=exons, minCov=minCov);
save(res,file = outRobjFN);

dat = cbind(res$rat5, res$cnt5, res$rat3, res$cnt3  )
colnames(dat)=c('ExonUsage5','coverage5', 'ExonUsage3','coverage3')
write.table(cbind(exons, dat), quote = FALSE, row.names = FALSE, sep='\t', file = outFN)

