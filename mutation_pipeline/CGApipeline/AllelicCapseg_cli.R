## Dump Frames if error
options(error = quote({dump.frames(to.file = TRUE); q()}))

library(optparse)
option.list <- list(
  make_option("--SID", dest="SID"),
  make_option("--capseg.probe.fn", dest="capseg.probe.fn"),
  make_option("--capseg.seg.fn", dest="capseg.seg.fn"),
  make_option("--germline.het.fn", dest="germline.het.fn"),
  make_option("--drop.x", dest="drop.x"),
  make_option("--drop.y", dest="drop.y"),
  make_option("--seg.merge.thresh", dest="seg.merge.thresh"),
  make_option("--min.seg.size", dest="min.seg.size"),
  make_option("--verbose", dest="verbose"),
  make_option("--base.output.dir", dest="base.output.dir"),
  make_option("--initial.merge", dest="initial.merge"),
  make_option("--split.merge", dest="split.merge"),
  make_option("--allelic.outlier.thresh", dest="allelic.outlier.thresh"),
  make_option("--working.dir", dest="working.dir"))

opt <- parse_args(OptionParser(option_list=option.list))
print(opt)
save(opt, file="debug.RData")

drop.x = as.logical(opt[["drop.x"]])
drop.y = as.logical(opt[["drop.y"]])
seg.merge.thresh = as.numeric(opt[["seg.merge.thresh"]])
min.seg.size = as.numeric(opt[["min.seg.size"]])
verbose = as.logical(opt[["verbose"]])
merge = as.logical(opt[["initial.merge"]])
split.merge = as.logical(opt[["split.merge"]])
allelic.outlier.thresh = as.numeric(opt[["allelic.outlier.thresh"]])

CODE_DIR = opt[["working.dir"]]
files_list=list("AllelicCapseg.R", "capture_model_fit.R", "seq_segfit_plots_allchr.R", "SegLabelCNLOH.R", "densities.R", "seq_segfit_plots_range.R", "Split.R", "file_parsers.R", "seq_segfit_plots_range_frac.R", "SplitFunctions.R", "seg_smooth_capture.R", "utils.R", "capseg_error_model.R", "seq_segfit_plots.R")
lapply(paste(CODE_DIR, "/", files_list, sep=""), source)

suppressMessages(library(doMC))
registerDoMC(1)

library(numDeriv)  ## needed for 'hessian()'

library("DNAcopy") ## needed for segmentation

AllelicCapseg( opt[["capseg.probe.fn"]], opt[["capseg.seg.fn"]], opt[["germline.het.fn"]], opt[["SID"]], opt[["base.output.dir"]], min.seg.size, drop.x, drop.y, seg.merge.thresh, verbose, merge, split.merge, allelic.outlier.thresh )
