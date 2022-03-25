#Function from PureCN
callCIN <- function(loh,allele.specific = TRUE, reference.state =
                      c("dominant", "normal")) {
  loh$size <- loh$end - loh$start + 1
  # should not happen
  loh <- loh[!is.na(loh$size), ]
  if (allele.specific) loh <- loh[!is.na(loh$M), ]
  reference.state <- match.arg(reference.state)
  loh$state <- if (allele.specific) paste0(loh$C, "/", loh$M) else loh$C
  dominant.state <-  sort(sapply(split(loh$size, loh$state), sum),
                          decreasing = TRUE)[1]
  reference.state.cn <- names(dominant.state)
  if (reference.state == "normal") {
    reference.state.cn <- if (allele.specific) "2/1" else "2"
  }
  loh$is.reference <- loh$state == reference.state.cn
  sum(loh$size[!loh$is.reference]) / sum(loh$size)
}

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

cin <- callCIN(loh,allele.specific = FALSE, reference.state = "normal")
cin.allele.specific = callCIN(loh, reference.state = "normal")
cin.ploidy.robust = callCIN(loh, allele.specific = FALSE)
cin.allele.specific.ploidy.robust = callCIN(loh)

wgd <-  -2*loh_frac + 3 < ploidy
writeLines(paste(as.character(wgd),
                 as.character(loh_frac),
                 as.character(cin),
                 as.character(cin.allele.specific),
                 as.character(cin.ploidy.robust),
                 as.character(cin.allele.specific.ploidy.robust),
                  sep = "\n"),"out.txt")

