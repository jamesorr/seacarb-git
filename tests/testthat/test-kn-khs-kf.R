## ---------------------------------------------------------------------------
## Ammonium and sulfide constants (Kn, Khs) must use the same Kf as every other
## constant, for every kf option and every pH scale.
##
## Until seacarb 3.4.1, calculate_carb() and the four buffer routines called Kn()
## and Khs() without the SWS-to-chosen-scale factor. Both then computed their own
## factor with kf = "x" (Perez and Fraga in range, Dickson and Riley outside), which
## differs from an explicit kf = "dg". With NH4t or HSt above zero, the carbonate
## system returned by carbfull() then depended on the pH scale (pCO2 by up to about
## 4e-6 relative) and so did Rf (about 1.6 ppm), on the seawater scale only.
## kf = "pf" hid the problem because the fallback agrees with it in range.
##
## GRID comes from helper-buffergen.R (states 4, 5, 8 and 9 have NH4t and/or HSt).
## ---------------------------------------------------------------------------

carb_kf <- function(i, sc, kf) {
  g <- GRID[i, ]
  suppressWarnings(carbfull(flag = 15, var1 = g$ALK, var2 = g$DIC, S = g$S, T = g$T,
    P = g$P, Pt = g$Pt, Sit = g$Sit, NH4t = g$NH4, HSt = g$HSt, k1k2 = "l",
    kf = kf, ks = "d", pHscale = sc, b = "u74"))
}

rf_kf <- function(i, sc, kf) {
  g <- GRID[i, ]
  a <- list(flag = 15, var1 = g$ALK, var2 = g$DIC, S = g$S, T = g$T, P = g$P,
            Pt = g$Pt, Sit = g$Sit, NH4t = g$NH4 + 1e-20, HSt = g$HSt,
            k1k2 = "l", kf = kf, ks = "d", pHscale = sc, b = "u74")
  c(buffsun = as.numeric(do.call(buffsun, c(a, warn = "n"))$Rf),
    buffzwg = as.numeric(do.call(buffzwg, c(a, warn = "n"))$RF),
    buffesm = as.numeric(do.call(buffesm, c(a, warn = "n"))$R),
    buffer  = as.numeric(do.call(buffer,  c(a, warn = "n"))$BetaD))
}

test_that("carbfull() returns the same pCO2 on every pH scale, for kf = pf and dg", {
  for (kf in c("pf", "dg")) {
    for (i in seq_len(nrow(GRID))) {
      ref <- as.numeric(carb_kf(i, "T", kf)$pCO2)
      for (sc in c("F", "SWS"))
        expect_equal(as.numeric(carb_kf(i, sc, kf)$pCO2), ref, tolerance = 1e-9)
    }
  }
})

test_that("buffsun, buffzwg, buffesm and buffer give the same Rf on every pH scale, for kf = pf and dg", {
  for (kf in c("pf", "dg")) {
    for (i in seq_len(nrow(GRID))) {
      ref <- rf_kf(i, "T", kf)
      for (sc in c("F", "SWS"))
        expect_equal(unname(rf_kf(i, sc, kf)), unname(ref), tolerance = 1e-8)
    }
  }
})

test_that("buffergen() gives the same Rf on every pH scale, for kf = pf and dg", {
  for (kf in c("pf", "dg")) {
    for (i in seq_len(nrow(GRID))) {
      ref <- as.numeric(bg(i, "T", kf = kf)$RF)
      for (sc in c("F", "SWS"))
        expect_equal(as.numeric(bg(i, sc, kf = kf)$RF), ref, tolerance = 1e-8)
    }
  }
})
