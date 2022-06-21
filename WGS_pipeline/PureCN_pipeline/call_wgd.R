#suppressMessages(library(tidyverse))

args = commandArgs(trailingOnly=TRUE)

in_loh <- args[1]
in_selected <- args[2]
#in_loh <- "~/Data/Copy_Number/example/CDS-00rz9N_loh.csv"
#in_selected <- "~/Data/Copy_Number/example/CDS-00rz9N.csv"

loh <- read.csv(in_loh)
loh$len  <- loh$end - loh$start

genome_len <-  sum(loh$len)
loh_len <- sum(loh[loh$type %in% c("LOH","COPY-NEUTRAL LOH","WHOLE ARM COPY-NEUTRAL LOH"),]$len)
loh_frac = loh_len/genome_len

selected <- read.csv(in_selected)
ploidy <- selected[1,"Ploidy"]

if (-2*loh_frac + 3 < ploidy) {
  writeLines(paste("TRUE\n",as.character(loh_frac),sep = ""),"out.txt")
} else {
  writeLines(paste("FALSE\n",as.character(loh_frac),sep = ""),"out.txt")
}
