## ---------------------------------------------------------------------------
## buffsun(), buffzwg(), buffesm() and buffer() on the free and seawater pH scales.
##
## These four routines build the noncarbonate alkalinity slope w = -dA_nc/dh from
## the acid-base systems that are SEPARATE members of the total alkalinity on the
## chosen pH scale (same rule as buffergen(); see test-buffergen.R):
##   free (F): sulfate and fluoride;  total (T): fluoride only;  SWS: neither.
## Until seacarb 3.4.1 they always used the total-scale rule, so on "F" and "SWS"
## the Revelle factor was wrong by up to about 13 ppm. The Revelle factor is a
## dimensionless ratio, so it cannot depend on the pH scale used to evaluate it.
##
## What is tested
##   1. each routine gives the same Rf on F, T and SWS;
##   2. each routine agrees with buffergen() (an independent implementation of the
##      derivatives) on every scale;
##   3. the total-scale result, which is the default, is unchanged.
## kf = "pf" is used because with kf = "dg" carbfull() itself returns a slightly
## different state on the SWS scale (see test-buffergen.R).
## A tiny NH4t (1e-20 mol/kg) is added so that every routine uses carbfull()
## (diprotic silicate), as buffergen() does. That changes nothing else.
## Tolerances are not to be widened to make a failing test pass.
## ---------------------------------------------------------------------------

rf_routines <- function(i, sc) {
  g <- GRID[i, ]
  a <- list(flag = 15, var1 = g$ALK, var2 = g$DIC, S = g$S, T = g$T, P = g$P,
            Pt = g$Pt, Sit = g$Sit, NH4t = g$NH4 + 1e-20, HSt = g$HSt,
            k1k2 = "l", kf = "pf", ks = "d", pHscale = sc, b = "u74")
  c(buffsun = as.numeric(do.call(buffsun, c(a, warn = "n"))$Rf),
    buffzwg = as.numeric(do.call(buffzwg, c(a, warn = "n"))$RF),
    buffesm = as.numeric(do.call(buffesm, c(a, warn = "n"))$R),
    buffer  = as.numeric(do.call(buffer,  c(a, warn = "n"))$BetaD),
    buffergen = as.numeric(suppressWarnings(do.call(buffergen, a))$RF))
}

test_that("buffsun, buffzwg, buffesm and buffer give the same Rf on every pH scale", {
  for (i in seq_len(nrow(GRID))) {
    ref <- rf_routines(i, "T")[["buffergen"]]
    for (sc in c("F", "T", "SWS")) {
      r <- rf_routines(i, sc)
      expect_equal(unname(r[c("buffsun", "buffzwg", "buffesm", "buffer")]),
                   rep(ref, 4), tolerance = 1e-8)
    }
  }
})

test_that("the four routines agree with buffergen() on each pH scale", {
  for (i in seq_len(nrow(GRID))) {
    for (sc in c("F", "T", "SWS")) {
      r <- rf_routines(i, sc)
      expect_equal(unname(r[c("buffsun", "buffzwg", "buffesm", "buffer")]),
                   rep(r[["buffergen"]], 4), tolerance = 1e-9)
    }
  }
})

test_that("total-scale Rf at the reference state is unchanged", {
  ## S=35, T=18, ALK=2300, DIC=2000 umol/kg, no nutrients (kf = "pf").
  a <- list(flag = 15, var1 = 2300e-6, var2 = 2000e-6, S = 35, T = 18,
            k1k2 = "l", kf = "pf", ks = "d", pHscale = "T", b = "u74", warn = "n")
  r <- c(do.call(buffsun, a)$Rf, do.call(buffzwg, a)$RF, do.call(buffesm, a)$R,
         do.call(buffer, a)$BetaD)
  expect_equal(as.numeric(r), rep(9.6913805691, 4), tolerance = 1e-9)
})
