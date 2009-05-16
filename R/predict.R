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



## includes predict.pls, and predict.spls

#`predict` <-
#function(object, ...) UseMethod("predict")


# ----------------------------PLS or sPLS ---------------------------
`predict.pls` <-
function(object, newdata, ...)
{
if (missing(newdata))
stop("No new data available.")

X = object$X
Y = object$Y
q = ncol(Y)
p = ncol(X)
mode = object$mode

# ---warnings --------------
if (length(dim(newdata)) == 2) {
if (ncol(newdata) != p)
stop("'newdata' must be a numeric matrix with ncol = ", p, 
" or a vector of length = ", p, ".")
}

if (length(dim(newdata)) == 0) {
if (length(newdata) != p)
stop("'newdata' must be a numeric matrix with ncol = ", p, 
" or a vector of length = ", p, ".")
}

if (mode == 'canonical') stop("Only regression, classic or invariant mode are allowed !")




ncomp = object$ncomp
a = object$loadings$X
b = object$loadings$Y
c = object$mat.c

means.X = attr(X, "scaled:center")
sigma.X = attr(X, "scaled:scale")


#coeff de regression 
B.hat <- array(0, dim = c(p, q, ncomp))
#prediction
Y.hat <- array(0, dim = c(nrow(newdata), q, ncomp))
#variates
t.pred = array(0, dim = c(nrow(newdata), ncomp))


#-- center and scale the new data --#
#newdata = scale(newdata, center = TRUE, scale = TRUE)
newdata = scale(newdata, center = means.X, scale = sigma.X)

for(h in 1:ncomp){
	W = a[, 1:h] %*% solve(t(c[, 1:h]) %*% a[, 1:h]) 
	B.hat[,,h] = W %*% t(b[, 1:h])   
	Y.hat[,,h] = newdata %*% B.hat[,,h] + colMeans(as.matrix(Y))
	t.pred[,h] = newdata %*% W[,h]
#	t.pred[,h] = scale(newdata, center = means.X, scale = sigma.X) %*% W[,h]
#	t.pred[,h] = t.pred[,h] / drop(crossprod(t.pred[,h]))
}  #end h


rownames(t.pred) = rownames(newdata)
colnames(t.pred) = paste('dim', c(1:ncomp), sep=' ')
colnames(Y.hat) = colnames(Y)

return(invisible(list(Y.hat = Y.hat, variates = t.pred)))
}



