# Copyright (C) 2007 Karline Soetaert (K.Soetaert@nioo.knaw.nl)
#
# This file is part of seacarb.
#
# Seacarb is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or any later version.
#
# Seacarb is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with seacarb; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

"Khs"              <- function (S=35,T=25,P=0,pHscale="T",warn="y")

#--------------------------------------------------------------
# Dissociation constant of H2S
#--------------------------------------------------------------

{

nK <- max(length(S), length(T), length(P), length(pHscale))

##-------- Creation de vecteur pour toutes les entrees (si vectorielles)

if(length(S)!=nK){S <- rep(S[1], nK)}
if(length(T)!=nK){T <- rep(T[1], nK)}
if(length(P)!=nK){P <- rep(P[1], nK)}
if(length(pHscale)!=nK){pHscale <- rep(pHscale[1], nK)}


tk = 273.15           # [K] (for conversion [deg C] <-> [K])
TC = T + tk           # TC [C]; T[K]

	#--------------------------------------------------------------
  # Dissociation constant of H2S on total scale - Millero 1995
 	#--------------------------------------------------------------
lnK <- 225.838 + 0.3449*sqrt(S) - 0.0274*S -13275.3/TC  -34.6435*log(TC)

Khs      <- exp(lnK)

# ---- Conversion from Total scale to seawater scale before pressure corrections
factor <- kconv(S=S, T=T, P=rep(0,nK))$ktotal2SWS
Khs <- Khs * factor

# ----------------- Pressure Correction ------------------	
Khs <- Pcorrect(Kvalue=Khs, Ktype="Khs", T=T, S=S, P=P, pHscale="SWS", warn=warn)


###----------------pH scale corrections
factor <- rep(NA,nK)
pHsc <- rep(NA,nK)
for(i in (1:nK)){   
 if(pHscale[i]=="T"){factor[i] <- kconv(S=S[i], T=T[i], P=P[i])$kSWS2total ; pHsc[i] <- "total scale"}
 if(pHscale[i]=="F"){factor[i] <- kconv(S=S[i], T=T[i], P=P[i])$kSWS2free ; pHsc[i] <- "free scale"}
 if(pHscale[i]=="SWS"){factor[i] <- 1 ; pHsc[i] <- "seawater scale"}
Khs[i] <- Khs[i]*factor[i]
}

##------------Warnings

    is_w <- warn == "y"

    if (any (is_w & (T>45 | T<0 | S<0 | S>45)) ) {warning("S and/or T is outside the range of validity of the formulation chosen for Khs.")}

  attr(Khs,"unit")     = "mol/kg-soln"
  attr(Khs,"pH scale") = pHsc
  return(Khs)

}  ## END Khs
