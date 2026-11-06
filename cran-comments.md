# v3.4.2
## Submission summary

seacarb 3.4.2 is a patch release. It corrects the Revelle factor
routines buffergen(), buffsun(), buffzwg(), buffesm(), and buffer() on
the free and seawater pH scales, and the ammonium and sulfide constants
(Kn and Khs) that were converted to the chosen pH scale with a different
Kf than all other constants when kf was set explicitly to "dg". Results
on the total pH scale with the default settings are unchanged. In
buffergen(), outputs that are not defined on a given pH scale now
return NA. The Kf help pages (21 files) were corrected to state that
both Kf formulations are on the free scale. Regression tests were added
(test-buffergen.R, test-buff-scales.R, test-kn-khs-kf.R).

## Test environments

* local Ubuntu 24.04.5, R 4.3.3, with NOT_CRAN=true and --as-cran

## R CMD check results

0 errors | 0 warnings | 3 notes

* CRAN incoming feasibility: maintainer note. The URL
  https://systemco2.net/t/27 (used in the Kf and related help pages)
  redirects to the full thread address and returns status 200. The short
  form is kept because the full address contains the thread title.
* checking for future file timestamps: unable to verify current time
  (the local machine could not check the time).
* checking HTML version of manual: math rendering skipped because
  package V8 is unavailable locally.

## Reverse dependencies

seacarb has three reverse dependencies on CRAN (ecolMod, marelac,
respirometry). We checked all three against this version (R 4.3.3,
Ubuntu 24.04) and found no problems. None of them uses the functions
changed in this release (the buffer routines, carbfull(), or explicit
kf settings): marelac uses rho(), respirometry uses carb() with default
arguments, Kn(), and rho(), and ecolMod only suggests seacarb.


# v3.4.1
## Submission summary

seacarb 3.4.1 is a patch release. It corrects the pH scale of the
Perez and Fraga (1987) hydrogen fluoride dissociation constant (Kf,
argument kf = "pf"), which was previously treated as total-scale and
is in fact free-scale (confirmed by F. Perez, 2026). The unwarranted
total-to-free conversion has been removed, the unused Ks_p0 argument
has been dropped from Kf(), and the internal callers of Kf() and the
Kf documentation have been updated accordingly. The Dickson and Riley
(1979) option (kf = "dg") is unaffected.

## Test environments

* local: Ubuntu 24.04, R 4.4.1
* GitHub Actions: ubuntu-latest (release, devel), macos-latest
  (release), windows-latest (release) — all passing

## R CMD check results

Status: 2 NOTEs, 0 errors, 0 warnings.

* checking for future file timestamps ... NOTE
  unable to verify current time

  This is a network note: the check could not reach a time server to
  verify timestamps. It is unrelated to the package.

* checking HTML version of manual ... NOTE
  Skipping checking math rendering: package 'V8' unavailable

  V8 is not installed in the local check environment. This note
  concerns the check environment, not the package, and does not appear
  on systems where V8 is available.

## Reverse dependencies

seacarb has three reverse dependencies on CRAN (ecolMod, marelac,
respirometry). We checked all three against this version and found no
problems. The kf = "pf" default now returns slightly different Kf
values (altering downstream carbonate-system results only in the
fourth significant figure), and none of the three relies on the
previous values.

