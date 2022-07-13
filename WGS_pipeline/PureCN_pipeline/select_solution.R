library(VariantAnnotation)
####################### 

This is a copy of a piece of code from PureCN to update results and plots from a solution to the duplication value.

#############
callAlterations <- function(res, id = 1, cutoffs = c(0.5, 6, 7),
                            log.ratio.cutoffs = c(-0.9, 0.9),
                            failed = NULL, all.genes = FALSE) {
  
  if (!is(res$results[[id]]$gene.calls, "data.frame")) {
    .stopUserError("This function requires gene-level calls.\n",
                   "Please add a column 'Gene' containing gene symbols to the ",
                   "interval.file.")
  }
  
  amp.ids <- (res$results[[id]]$gene.calls$focal &
                res$results[[id]]$gene.calls$C >= cutoffs[2]) |
    res$results[[id]]$gene.calls$C >= cutoffs[3]
  
  del.ids <- res$results[[id]]$gene.calls$C < cutoffs[1]
  
  if (is.null(failed)) failed <- res$results[[id]]$failed
  
  if (failed) {
    amp.ids <- res$results[[id]]$gene.calls$gene.mean >=
      log.ratio.cutoffs[2]
    del.ids <- res$results[[id]]$gene.calls$gene.mean <
      log.ratio.cutoffs[1]
  }
  
  calls <- res$results[[1]]$gene.calls
  calls$type <- NA
  calls$type[amp.ids] <- "AMPLIFICATION"
  calls$type[del.ids] <- "DELETION"
  
  bm <- res$results[[id]]$SNV.posterior
  if (!is.null(bm)) {
    segids <- bm$posteriors$seg.id
    calls$num.snps <- sapply(calls$seg.id, function(i)
      sum(segids == i, na.rm = TRUE))
    calls$M <- bm$posteriors$ML.M.SEGMENT[match(calls$seg.id, segids)]
    calls$M.flagged <- bm$posteriors$M.SEGMENT.FLAGGED[match(calls$seg.id, segids)]
    calls$loh <- bm$posteriors$ML.M.SEGMENT[match(calls$seg.id, segids)] == 0
  }
  
  calls <- calls[, !grepl("^\\.", colnames(calls))]
  
  if (!all.genes) {
    return(calls[!is.na(calls$type), ])
  }
  calls
}

predictSomatic <- function(res, id = 1, return.vcf = FALSE) {
  pp   <- .addSymbols(res$results[[id]])
  if (return.vcf) {
    vcf <- res$input$vcf[
      res$results[[id]]$SNV.posterior$vcf.ids]
    return(.annotatePosteriorsVcf(pp, vcf))
  }
  pp
}

.calcMultisamplePosteriors <- function(ret.list) {
  pool <- Reduce(union, lapply(ret.list, function(r)
    as.character(rowRanges(r$input$vcf))))
  pos <- lapply(ret.list, function(r) match(as.character(GRanges(r$results[[1]]$SNV.posterior$posteriors)), pool))
  pos.align <- lapply(seq_along(pool), function(i)
    lapply(pos, function(j) which(j == i)))
  likelihoods <- lapply(pos.align, function(i) do.call(rbind, lapply(seq_along(i), function(j)
    ret.list[[j]]$results[[1]]$SNV.posterior$likelihoods[i[[j]],])))
  idx <- grep("SOMATIC.M", colnames(likelihoods[[1]]))
  pp <- sapply(likelihoods, function(l) 
    sum(apply(l,1,function(x) sum(x[idx]))) /
      sum(apply(l,1,function(x) sum(x))))
  names(pp) <- pool
  lapply(ret.list, function(x) {
    ps <- predictSomatic(x)
    keys <- as.character(GRanges(ps))
    data.frame(Sampleid = x$input$sampleid,
               ps,
               MULTI.POSTERIOR.SOMATIC = pp[keys])
  })
}       

.addSymbols <- function(result) {
  if (is(result$gene.calls, "data.frame")) {
    g.gr <- GRanges(result$gene.calls)
    p.gr <- GRanges(result$SNV.posterior$posteriors)
    ov <- findOverlaps(p.gr, g.gr)
    
    result$SNV.posterior$posteriors$gene.symbol <- NA
    result$SNV.posterior$posteriors$gene.symbol[queryHits(ov)] <- 
      rownames(result$gene.calls)[subjectHits(ov)]
  }    
  result$SNV.posterior$posteriors
}

.annotatePosteriorsVcf <- function(pp, vcf) {
  if (nrow(pp) != nrow(vcf)) {
    .stopRuntimeError("Posteriors and filtered VCF do not align.")
  }
  if (is.null(pp$gene.symbol)) pp$gene.symbol <- "."
  
  idxColsPp <- grep("^SOMATIC|^GERMLINE", colnames(pp))
  colnames(pp)[idxColsPp] <- gsub("OMATIC\\.|ERMLINE\\.", "",
                                  colnames(pp)[idxColsPp])
  colnames(pp)[colnames(pp) == "CN.SUBCLONAL"] <- "CS"
  colnames(pp)[colnames(pp) == "CELLFRACTION"] <- "CF"
  colnames(pp)[colnames(pp) == "POSTERIOR.SOMATIC"] <- "PS"
  colnames(pp)[colnames(pp) == "log.ratio"] <- "LR"
  colnames(pp)[colnames(pp) == "gene.symbol"] <- "GS"
  
  descriptionPp <- paste0(ifelse(grepl("^S", colnames(pp)[idxColsPp]),
                                 "Somatic", "Germline"), gsub("^.*M(\\d)", 
                                                              " posterior probability, multiplicity \\1", 
                                                              colnames(pp)[idxColsPp]))
  
  descriptionPp <- gsub("GCONTHIGH", 
                        " homozygous, reference allele contamination", descriptionPp)
  descriptionPp <- gsub("GCONTLOW", 
                        " alt allele contamination", descriptionPp)
  descriptionPp <- gsub("GHOMOZYGOUS", 
                        " homozygous", descriptionPp)
  idxCols <- grep("^ML", colnames(pp))
  prefix <- .getPureCNPrefixVcf(vcf)
  colnames(pp) <- paste0(prefix, colnames(pp))
  newInfoPosterior <- DataFrame(
    Number = 1, 
    Type = "Float", 
    Description = descriptionPp, 
    row.names = colnames(pp)[idxColsPp]
  )
  infoType <- as.character(sapply(sapply(pp[idxCols], class), function(x)
    ifelse(x=="logical", "Flag", ifelse(x=="integer", "Integer", "Float"))))
  newInfoMl <- DataFrame(
    Number = ifelse(infoType == "Flag", 0, 1),
    Type = infoType, 
    Description = c(
      "Maximum likelihood state is a somatic state",
      "Maximum likelihood multiplicity",
      "Maximum likelihood integer copy number",
      "Maximum likelihood minor segment copy number",
      "Expected allelic fraction of the maximum likelihood state",
      "TRUE if segment is most likely in LOH"
    ),
    row.names = colnames(pp)[idxCols]
  )
  
  idxColsMisc <- paste0(prefix, c("CS", "CF", "PS", "LR", "GS", "FLAGGED"))
  newInfoMisc <- DataFrame(
    Number = c(0, 1, 1, 1, 1, 0),
    Type = c("Flag", "Float", "Float", "Float", "String", "Flag"),
    Description = c("Sub-clonal copy number gain", "Cellular fraction",
                    "Posterior probability variant is somatic mutation.",
                    "Copy number log-ratio", "Gene Symbol", "QC Flagged"),
    row.names = idxColsMisc
  )
  info(header(vcf)) <- rbind(info(header(vcf)), newInfoPosterior, newInfoMl,
                             newInfoMisc)
  idxColsMisc <- match(idxColsMisc, colnames(pp))
  pp[, idxColsPp] <- round(pp[, idxColsPp], digits = 4)
  pp[, idxColsMisc[2:4]] <- round(pp[, idxColsMisc[2:4]], digits = 4)
  
  info(vcf) <- cbind(info(vcf), pp[, c(idxColsPp, idxCols, idxColsMisc)])
  vcf
}

.getSexFromRds <- function(rds) {
  # if run without VCF, then we don't have sex information from VCF
  if (is.null(rds$input$sex.vcf)) return(rds$input$sex)
  
  # conflict of coverage and snp based sex genotyper?
  if (!is.na(rds$input$sex) && !is.na(rds$input$sex.vcf)) {
    if (rds$input$sex == rds$input$sex.vcf) return(rds$input$sex)
    return(paste("Coverage:", rds$input$sex, "VCF:", rds$input$sex.vcf))
  }
  # believe coverage based more than VCF in case we have only limited
  # number of SNPs on chrX
  if (!is.na(rds$input$sex)) {
    return(rds$input$sex)
  }
  return(rds$input$sex.vcf)
}

args = commandArgs(trailingOnly=TRUE)

#file.rds <- "~/Data/Copy_Number/example_3/CDS-00rz9N.rds"
#file.rds <- "~/Data/Copy_Number/example_4/PureCN_CDS-051M8l.rds"
#i <- 3
file.rds <- args[1]
print(file.rds)
out <- args[2]
print(out)
i <- as.numeric(args[3])
print(i)

ret <- readRDS(file.rds)
print(names(ret))

### Create output files -------------------------------------------------------

res <- ret$results[[i]]
sampleid <- res$seg$ID[1]
print(sampleid)
contamination <- res$SNV.posterior$posterior.contamination
contamination <- if (is.null(contamination)) 0 else contamination
d.f.curation <- data.frame(
  Sampleid = sampleid,
  Purity = res$purity,
  Ploidy = res$ploidy,
  Sex = .getSexFromRds(ret),
  Contamination = contamination,
  Flagged = res$flag,
  Failed = FALSE,
  Curated = TRUE,
  Comment = res$flag_comment
)

write.csv(d.f.curation, file = paste0(out, ".csv"), row.names = FALSE)

file.seg <- paste0(out, "_dnacopy.seg")
seg <- ret$results[[i]]$seg
seg <- seg[, c(1:6, match("C", colnames(seg)))]
write.table(seg, file = file.seg, sep = "\t", quote = FALSE,
            row.names = FALSE)

if (is(ret$results[[i]]$gene.calls, "data.frame")) {
  file.genes <- paste0(out, "_genes.csv")
  allAlterations <- callAlterations(ret, id = i, all.genes = TRUE)
  
  write.csv(cbind(Sampleid = sampleid, gene.symbol = rownames(allAlterations),
                  allAlterations), row.names = FALSE, file = file.genes, quote = FALSE)
} else {
  flog.warn("--intervals does not contain gene symbols. Not generating gene-level calls.")
}

if (!is.null(ret$input$vcf) &&
    !is.null(ret$results[[i]]$SNV.posterior)) {

  file.csv <- paste0(out, "_variants.csv")
  write.csv(cbind(Sampleid = sampleid, predictSomatic(ret, id = i)), file = file.csv,
            row.names = FALSE, quote = FALSE)
  
  file.loh <- paste0(out, "_loh.csv")
  write.csv(cbind(Sampleid = sampleid, PureCN::callLOH(ret,id = i)), file = file.loh,
            row.names = FALSE, quote = FALSE)
}

