## ---------------------------------------------------------------------------
## What the buffergen() fixes gained: three versions compared on the same grid.
##
##   v3.4.1   the released function
##   patch 1  sign fixes for fluoride and sulfate, explicit HSO4 / HF only on the
##            pH scales where they are alkalinity terms
##   patch 2  patch 1 plus NA for derivatives that are not defined on a pH scale
##
## Versions are read from git so nothing has to be copied around. Edit REFS if your
## commits are named differently. If patch 2 is committed, use "HEAD~1" for patch 1
## and "HEAD" for patch 2 instead of "working".
##
## Run from this directory:  Rscript buffergen_before_after.R
## It reuses GRID, bg() and indep() from ../tests/testthat/helper-buffergen.R, which
## evaluate the fluoride and sulfate derivatives independently of buffergen().
## ---------------------------------------------------------------------------
source("setup.R")
source("../tests/testthat/helper-buffergen.R")

REFS <- c("v3.4.1" = "v3.4.1", "patch 1" = "HEAD", "patch 2" = "working")

get_fun <- function(ref) {
  txt <- if (ref == "working") readLines("../R/buffergen.R")
         else system2("git", c("show", paste0(ref, ":R/buffergen.R")), stdout = TRUE)
  e <- new.env(); eval(parse(text = txt), e); e$buffergen
}

FLU <- c("F", "HF", "FLUO"); SUL <- c("SO4", "HSO4", "ST")
Q8  <- c("dALK.dH", "dtotX.dH", "dALK.dX", "dtotX.dX", "dH.dALK", "dH.dtotX", "dX.dALK", "dX.dtotX")

cat(sprintf("\n%-9s %14s %14s %15s %16s\n", "version", "RF vs buffsun", "RF spread over", "defined entries", "undefined entries"))
cat(sprintf("%-9s %14s %14s %15s %16s\n", "", "(max, T scale)", "F, T, SWS (max)", "wrong (>1e-6)", "not NA"))
for (v in names(REFS)) {
  f <- get_fun(REFS[[v]])
  ## RF against buffsun() on the total scale (diprotic silicate forced with a tiny NH4t)
  dev <- 0; spread <- 0
  for (i in seq_len(nrow(GRID))) {
    g <- GRID[i, ]
    rf <- as.numeric(bg(i, "T", fun = f)$RF)
    ref <- buffsun(flag = 15, var1 = g$ALK, var2 = g$DIC, S = g$S, T = g$T, P = g$P, Pt = g$Pt,
      Sit = g$Sit, NH4t = g$NH4 + 1e-20, HSt = g$HSt, k1k2 = OPT$k1k2, kf = OPT$kf, ks = OPT$ks,
      pHscale = "T", b = OPT$b)$Rf
    dev <- max(dev, abs(rf / ref - 1))
    r <- sapply(c("F", "T", "SWS"), function(sc) as.numeric(bg(i, sc, kf = "pf", fun = f)$RF))
    spread <- max(spread, diff(range(r)) / r[["T"]])
  }
  ## derivatives of the fluoride and sulfate systems, against the alkalinity definition
  nwrong <- 0; ndef <- 0; nundef <- 0; nnum <- 0
  for (i in seq_len(nrow(GRID))) for (sc in c("F", "T", "SWS")) {
    ind <- indep(i, sc, fun = f)
    defined <- if (sc == "F") c(FLU, SUL) else if (sc == "T") FLU else character(0)
    for (n in c(FLU, SUL)) for (q in Q8) {
      b <- as.numeric(ind$g[[q]][n])
      if (n %in% defined) {
        ndef <- ndef + 1
        if (!is.finite(b) || relerr(b, ind$truth[[n]][[q]]) > 1e-6) nwrong <- nwrong + 1
      } else {
        nundef <- nundef + 1
        if (!is.na(b)) nnum <- nnum + 1
      }
    }
  }
  cat(sprintf("%-9s %14.1e %14.1e %9d of %3d %11d of %3d\n", v, dev, spread, nwrong, ndef, nnum, nundef))
}
cat("\nRF vs buffsun and RF spread are relative. Entries are the eight derivative outputs for the\n",
    "fluoride and sulfate species, over the nine states of GRID on the F, T and SWS scales.\n", sep = "")
