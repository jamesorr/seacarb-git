## ---------------------------------------------------------------------------
## Helpers for test-buffergen.R and notebooks/buffergen_before_after.R.
## GRID of states, a wrapper around buffergen(), and indep(): an INDEPENDENT
## numerical evaluation of the fluoride and sulfate derivatives from the
## alkalinity definition of each pH scale (see test-buffergen.R for the rationale).
## Deliberately not simplified or refactored to call buffergen() internals.
## ---------------------------------------------------------------------------
OPT <- list(k1k2 = "l", kf = "dg", ks = "d", b = "u74")

## S, T, P, ALK, DIC, Pt, Sit, NH4t, HSt: coastal and polar states, with and
## without nutrients, and with pressure.
GRID <- data.frame(
  S   = c(35, 35, 35, 38, 20, 35, 35, 33, 36),
  T   = c(18, 18, -1.5, 30, 10, 25, 2, 15, 22),
  P   = c(0, 300, 0, 0, 0, 3000, 0, 0, 500),
  ALK = c(2300, 2300, 2300, 2380, 2200, 2300, 2250, 2400, 2500) * 1e-6,
  DIC = c(2000, 2000, 2150, 1950, 2100, 2000, 2100, 2150, 2100) * 1e-6,
  Pt  = c(0, 0, 1.5, 0.1, 2, 0, 1, 3, 5) * 1e-6,
  Sit = c(0, 0, 60, 2, 50, 0, 20, 80, 100) * 1e-6,
  NH4 = c(0, 0, 0, 0.5, 3, 0, 0, 5, 10) * 1e-6,
  HSt = c(0, 0, 0, 0, 2, 0, 0, 5, 10) * 1e-6)

bg <- function(i, sc, kf = OPT$kf, fun = buffergen) {
  g <- GRID[i, ]
  suppressWarnings(fun(flag = 15, var1 = g$ALK, var2 = g$DIC, S = g$S, T = g$T,
    P = g$P, Pt = g$Pt, Sit = g$Sit, NH4t = g$NH4, HSt = g$HSt,
    k1k2 = OPT$k1k2, kf = kf, ks = OPT$ks, pHscale = sc, b = OPT$b))
}

## Relative error that is safe for exact zeros
relerr <- function(a, b) abs(a - b) / pmax(abs(b), 1e-300)

## Central difference with one Richardson step. 'rel' is the relative step.
d1 <- function(f, x, rel) {
  h1 <- x * rel; h2 <- h1 / 2
  D1 <- (f(x + h1) - f(x - h1)) / (2 * h1)
  D2 <- (f(x + h2) - f(x - h2)) / (2 * h2)
  (4 * D2 - D1) / 3
}

## ---------------------------------------------------------------------------
## Independent derivative check
## ---------------------------------------------------------------------------
indep <- function(i, sc, fun = buffergen) {
  g <- bg(i, sc, fun = fun)
  cf <- as.list(as.data.frame(g$Carbfull))
  H0 <- cf$H
  eH <- (sc == "F")                 # HSO4 is a separate alkalinity term
  eF <- sc %in% c("T", "F")         # HF is a separate alkalinity term
  K1 <- H0 * cf$HCO3 / cf$CO2; K2 <- H0 * cf$CO3 / cf$HCO3
  Kb <- H0 * cf$BOH4 / cf$BOH3; Kw <- H0 * cf$OH
  Kf <- H0 * cf$F / cf$HF; Ks <- H0 * cf$SO4 / cf$HSO4
  DIC <- cf$DIC; BT <- cf$BOR; FT0 <- cf$FLUO; ST0 <- cf$ST
  Pt <- cf$Pt; Sit <- cf$Sit; Nt <- cf$NH4t; Ht <- cf$HSt
  if (Pt > 0)  { P1 <- H0 * cf$H2PO4 / cf$H3PO4; P2 <- H0 * cf$HPO4 / cf$H2PO4; P3 <- H0 * cf$PO4 / cf$HPO4 }
  if (Sit > 0) { Q1 <- H0 * cf$SiOOH3 / cf$SiOH4; Q2 <- H0 * cf$SiO2OH2 / cf$SiOOH3 }
  if (Nt > 0)  Kn <- H0 * cf$NH3 / cf$NH4
  if (Ht > 0)  Khs <- H0 * cf$HS / cf$H2S

  ## everything except the fluoride and sulfate systems, at fixed totals
  rest <- function(H) {
    a <- DIC * (K1 * H + 2 * K1 * K2) / (H^2 + K1 * H + K1 * K2) + BT * Kb / (Kb + H) + Kw / H - H
    if (Pt > 0)  a <- a + Pt * (P1 * P2 * H + 2 * P1 * P2 * P3 - H^3) /
                          (H^3 + P1 * H^2 + P1 * P2 * H + P1 * P2 * P3)
    if (Sit > 0) a <- a + Sit * (Q1 * H + 2 * Q1 * Q2) / (H^2 + Q1 * H + Q1 * Q2)
    if (Nt > 0)  a <- a + Nt * Kn / (Kn + H)
    if (Ht > 0)  a <- a + Ht * Khs / (Khs + H)
    a
  }
  ALKf <- function(H, FT = FT0, ST = ST0) {
    a <- rest(H)
    if (eF) a <- a - FT * H / (Kf + H)
    if (eH) a <- a - ST * H / (Ks + H)
    a
  }
  spec <- function(H, FT, ST) list(F = FT * Kf / (Kf + H), HF = FT * H / (Kf + H),
    FLUO = FT * Kf / (Kf + H), SO4 = ST * Ks / (Ks + H), HSO4 = ST * H / (Ks + H),
    ST = ST * Ks / (Ks + H))

  ## Step sizes: the carbonate and borate terms vary on the scale of [H+] itself, so
  ## the alkalinity itself is differentiated with a small relative step (1e-3). The
  ## fluoride and sulfate species vary on the scale of Kf and Ks, which are orders
  ## of magnitude larger than [H+], so a larger relative step (0.1) is exact for them
  ## and avoids round-off noise.
  dAdH <- d1(function(h) ALKf(h), H0, 1e-3)

  out <- list()
  for (n in c("F", "HF", "FLUO", "SO4", "HSO4", "ST")) {
    fl <- n %in% c("F", "HF", "FLUO")
    Kg <- if (fl) Kf else Ks
    expl <- if (fl) eF else eH
    base <- c(F = "F", FLUO = "F", HF = "HF", SO4 = "SO4", ST = "SO4", HSO4 = "HSO4")[[n]]
    acid <- base %in% c("HF", "HSO4")
    X <- spec(H0, FT0, ST0)[[base]]
    acidc <- function(H, x) if (acid) x else x * H / Kg
    totx  <- function(H, x) x + if (acid) x * Kg / H else x * H / Kg
    other <- function(H) {
      o <- rest(H)
      if (fl && eH)  o <- o - ST0 * H / (Ks + H)
      if (!fl && eF) o <- o - FT0 * H / (Kf + H)
      o
    }
    r <- list()
    r$dALK.dH  <- d1(function(h) other(h) - if (expl) acidc(h, X) else 0, H0, 1e-3)
    r$dtotX.dH <- d1(function(h) totx(h, X), H0, if (acid) 1e-3 else 0.1)   # acid form varies as 1/[H+]
    r$dALK.dX  <- d1(function(x) other(H0) - if (expl) acidc(H0, x) else 0, X, 0.1)
    r$dtotX.dX <- d1(function(x) totx(H0, x), X, 0.1)
    r$dH.dALK  <- 1 / dAdH
    sp1 <- function(H) spec(H, FT0, ST0)[[n]]
    dXdH <- d1(sp1, H0, 0.1)
    if (fl) {
      dAdT <- d1(function(t) ALKf(H0, FT = t), FT0, 0.1)
      dXdT <- d1(function(t) spec(H0, t, ST0)[[n]], FT0, 0.1)
    } else {
      dAdT <- d1(function(t) ALKf(H0, ST = t), ST0, 0.1)
      dXdT <- d1(function(t) spec(H0, FT0, t)[[n]], ST0, 0.1)
    }
    r$dH.dtotX <- -dAdT / dAdH
    r$dX.dALK  <- dXdH * r$dH.dALK
    r$dX.dtotX <- dXdT + dXdH * r$dH.dtotX
    out[[n]] <- r
  }
  list(truth = out, g = g)
}

OUTS8 <- c("dALK.dH", "dtotX.dH", "dALK.dX", "dtotX.dX",
           "dH.dALK", "dH.dtotX", "dX.dALK", "dX.dtotX")
OUTS12 <- c(OUTS8, "dALK.dpH", "dtotX.dpH", "dpH.dALK", "dpH.dtotX")
FLUORIDE <- c("F", "HF", "FLUO")
SULFATE  <- c("SO4", "HSO4", "ST")
