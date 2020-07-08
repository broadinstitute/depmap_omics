library(optparse)
option.list <- list(
  make_option("--seg_dat_fn", dest="seg_dat_fn"),
  make_option("--maf_fn", dest="maf_fn"),
  make_option("--indelmaf_fn", dest="indelmaf_fn"),
  make_option("--sample_name", dest="sample_name"),
  make_option("--results_dir", dest="results_dir",
              default=getwd()),
  make_option("--ssnv_skew", dest="ssnv_skew",
              default=0.9883274),
  make_option("--abs_lib_dir", dest="abs_lib_dir"),
  make_option("--platform", dest="platform", default="Illumina_WES"),
  make_option("--force_alpha", dest="f_a",type="double",default=NA),
  make_option("--force_tau", dest="f_t",type="double",default=NA))
#--seg_dat_fn --maf_fn --indelmaf_fn --sample_name --results_dir --ssnv_skew

#"""Rscript /run/ABSOLUTE_cli_start.R \
# --seg_dat_fn {seg_file} --maf_fn {snv_maf} --indelmaf_fn {indel_maf} --sample_name {sample_name} --results_dir results \
#--ssnv_skew {skew} """

# --abs_lib_dir /xchip/tcga/Tools/absolute/releases/v1.5/

opt <- parse_args(OptionParser(option_list=option.list))
print(opt)

min.ploidy= 1.1
max.ploidy= 6
sigma.h= 0.01
max.non.clonal= 0.99
min_probes= 1
max.as.seg.count= 5000
primary.disease= NA
platform= opt[["platform"]]
copy_num_type= "allelic"
genome_build= "hg19"
#CGA_DIR= "~/CGA/R/"
max.neg.genome= 0.05
maf.fn= opt[["maf_fn"]]
indel.maf.fn= opt[["indelmaf_fn"]]
min.mut.af= 0
output.fn.base= opt[["sample_name"]]
max_sd= 100
SSNV_skew= opt[["ssnv_skew"]]
filter_segs= TRUE
verbose= TRUE
results.dir= opt[["results_dir"]]
sample.name= opt[["sample_name"]]
seg.dat.fn= opt[["seg_dat_fn"]]
gender= NA
force.alpha= opt[["f_a"]]
force.tau= opt[["f_t"]]

#library(ABSOLUTE)


require(GenomicRanges)
require(gplots)
require(RColorBrewer)

#suppressPackageStartupMessages(require(gplots))
#suppressPackageStartupMessages(require(RColorBrewer))
CGA_DIR_ABS=opt[["abs_lib_dir"]] ##"/soft/local/absolute" #/xchip/tcga/Tools/absolute/releases/v1.5/

   print( paste("sourcing files in ", CGA_DIR_ABS, sep=""))
   #rr = dir( file.path( CGA_DIR_ABS, "ABSOLUTE/sandbox"), full.names=TRUE )
   rr = dir(CGA_DIR_ABS,full.names=TRUE )
   for( i in 1:length(rr) ) { 
   source(rr[i])}


# moved into the absolute root dir
#   source( file.path( CGA_DIR_ABS, "Phylogic", "log_DP_post_fixed.R" ) )
#   source( file.path( CGA_DIR_ABS, "Phylogic", "log_DP_post.R" ) )
#   source( file.path( CGA_DIR_ABS, "Phylogic", "run_DP.R" ) )
#   source( file.path( CGA_DIR_ABS, "Phylogic", "DP_utils.R" ) )
#   source( file.path( CGA_DIR_ABS, "Phylogic", "summarize_DP.R" ) )
#   source( file.path( CGA_DIR_ABS, "Phylogic", "CCF_DP_plots.R" ) )



#debug (RunAbsolute)
load(file.path( CGA_DIR_ABS,"data","ChrArmsDat.RData"))
load(file.path( CGA_DIR_ABS,"data","ChrArmPriorDb.RData"))
load(file.path( CGA_DIR_ABS,"data","diseaseMap.RData"))
RunAbsolute( seg.dat.fn, primary.disease, platform, sample.name, results.dir, copy_num_type,
                        genome_build, gender, min.ploidy, max.ploidy,
                        max.as.seg.count, max.non.clonal, max.neg.genome,
                        maf.fn, indel.maf.fn, min.mut.af,
                        output.fn.base, min_probes, max_sd, sigma.h, SSNV_skew,
			filter_segs, force.alpha, force.tau, verbose )


file.base = paste( output.fn.base, ".ABSOLUTE", sep = "")
absolute.files= file.path( results.dir, paste(file.base, "RData", sep = "."))
print (absolute.files)
CreateReviewObject( sample.name, absolute.files, ".", copy_num_type, plot.modes=TRUE, num_solutions_plotted=num_solutions_plotted, verbose=TRUE)
