## reads selected fields from maf files and aggregates them 
args<-commandArgs(TRUE)
print(args)

outfn = args[1];
inFNs = args[2];

SelectFields = as.character(args[[3]])  ## use SelectFields = 'Chromosome,Start_position,t_alt_count,t_ref_count,Variant_Type' for PoN filter
selFields=unlist(strsplit(SelectFields,','))


#inFNs = '/xchip/cle/analysis/mghandi/CCLE/R/ExonUsage/fns.txt'

fns = as.character(read.csv(inFNs, header=FALSE)[,])
fe = sapply(fns, file.exists)
samples = fns[which(fe)]

nsamples = length(samples)

readsel = function(fn, selFields){
	print(fn)
  f = read.delim(fn, sep='\t', comment.char = '#', header = TRUE, quote='')
  ii = match(selFields,colnames(f))
  if(length(which(is.na(ii)))>0){
    cat('\n Warning: some fields not found: ', fn,' ', selFields[which(is.na(ii))])
  }
  res = f[,ii]
  print(length(res))
}

for(i in 1:nsamples){
  ri = readsel(samples[i],selFields)
  write.table(ri,  file=outfn,  sep='\t', row.names = FALSE, quote = FALSE, append = (i>1), col.names = !(i>1))
}




