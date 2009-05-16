# Copyright (C) 2009 
# Sébastien Déjean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio González, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Lê Cao, French National Institute for Agricultural Research and 
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


# ---------------------------print for summary with (s)PLS object or rcc ------------------------------

'print.summary' <-
function(x, ...){

print.gap = 4
what = x$what
digits = x$digits


if(x$method == "pls" | x$method == "spls"){

if (any(what == "all") || any(what == "summarised") || any(what == "communalities") || any(what == "redundancy") || any(what == "VIP")) {
#cat("\n\n", "Data:   X dimension:", c(n, p), "\n")
#cat("         Y dimension:", c(n, q), "\n\n")
#cat(" Number of components considered:", ncomp, "\n\n")

if (x$method == "pls")
cat(" PLS mode:", x$mode, "\n")
else
cat(" sPLS mode:", x$mode, "\n")

}


#-------------------------------------------affichage communauté --#
#--------------------------#
if (any(what == "all") || any(what == "communalities")) { 
cat("\n\n Communalities Analysis:\n",
"----------------------")

cat("\n X-Variables vs their own Variates: see object$CM.X$own   \n")

##print(round(x$Cm.X$own, digits = digits), print.gap = print.gap)

cat("\n X-Variables vs the opposite Variates: see object$CM.X$opp   \n")
##print(round(x$Cm.X$opp, digits = digits), print.gap = print.gap)

cat("\n Y-Variables vs their own Variates: see object$CM.Y$opp  \n")
##print(round(x$Cm.Y$own, digits = digits), print.gap = print.gap)

cat("\n Y-Variables vs the opposite Variates:  see object$CM.Y$opp \n")
##print(round(x$Cm.Y$opp, digits = digits), print.gap = print.gap)
}

#----------------------------------------------affichage redondance --#
#--------------------------#
if (any(what == "all") || any(what == "redundancy")) {
cat("\n\n Redundancy Analysis:\n",
"-------------------\n")

cat(" X-Variables vs:                                             \n",
"            Their own Variates         The opposite Variates \n",
"         ------------------------    ------------------------\n")

affich = cbind(x$Rd.X$own, x$Rd.X$opp) 
print(round(affich, digits = digits), print.gap = print.gap)

cat(                                                           "\n\n",
"Y-Variables vs:                                              \n",
"            Their own Variates         The opposite Variates \n",
"         ------------------------    ------------------------\n")

affich = cbind(x$Rd.Y$own, x$Rd.Y$opp)
print(round(affich, digits = digits), print.gap = print.gap)
}


#----------------------------------------tableau VIP --#
#-----------------#
if (any(what == "all") || any(what == "VIP")) {

cat("\n\n", "Variable Importance in the Projection (VIP): see object$VIP \n",
"-------------------------------------------\n\n")
##print(x$VIP, digits = digits, print.gap = print.gap)

}

}  #end if pls


# --------------------------- output rcc ---------------------------

if(x$method == "rcc" ){

print.gap = 4
if (any(what == "all") || any(what == "summarised")) {

#cat("\n\n Data:   X dimension:", c(n, p), "\n")
#cat("         Y dimension:", c(n, q), "\n\n")
cat(" Number of canonical variates considered:", x$ncomp, "\n")
cat("\n Canonical correlations:",
"\n ----------------------\n")
print(round(x$can.cor, digits = digits), print.gap = print.gap)
}

#-- affichage communauté --#
#--------------------------#
if (any(what == "all") || any(what == "communalities")) { 
cat("\n\n Canonical Communalities Analysis:\n",
"--------------------------------")

cat("\n X-Variables vs their own Canonical Variates: see object$Cm.X$own   \n")
#print(round(x$Cm.X$own, digits = digits), print.gap = print.gap)

cat("\n X-Variables vs the opposite Canonical Variates: see object$Cm.X$opp    \n")
#print(round(x$Cm.X$opp, digits = digits), print.gap = print.gap)

cat("\n Y-Variables vs their own Canonical Variates: see object$Cm.Y$own  \n")
#print(round(x$Cm.Y$own, digits = digits), print.gap = print.gap)

cat("\n Y-Variables vs the opposite Canonical Variates: see object$Cm.Y$opp  \n")
#print(round(x$Cm.Y$opp, digits = digits), print.gap = print.gap)
}

#-- affichage redondance --#
#--------------------------#
if (any(what == "all") || any(what == "redundancy")) {
cat("\n\n Canonical Redundancy Analysis:\n",
"-----------------------------\n")

cat(" X-Variables vs:                                             \n",
"                 Their own                 The opposite      \n",
"            Canonical Variates          Canonical Variates   \n",
"         ------------------------    ------------------------\n")

affich = cbind(x$Rd.X$own, x$Rd.X$opp)
print(round(affich, digits = digits), print.gap = print.gap)

cat(                                                           "\n\n",
"Y-Variables vs:                                              \n",
"                 Their own                 The opposite      \n",
"            Canonical Variates          Canonical Variates   \n",
"         ------------------------    ------------------------\n")

affich = cbind(x$Rd.Y$own, x$Rd.Y$opp)
print(round(affich, digits = digits), print.gap = print.gap)
}

}  #end rcc



}

