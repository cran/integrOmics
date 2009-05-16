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




`cim.rcc` <-
function(object, dim, X.names = NULL, Y.names = NULL, ...) 
{

p = ncol(object$X)
q = ncol(object$Y)

if (!is.numeric(dim) || dim < 1)
    stop("invalid number for the dimensionality 'dim'.")
    else if(dim > min(p, q)) {
warning("Reset the dimensionality 'dim' to min(ncol(X), ncol(Y)) = ", min(p, q))
dim = min(p, q)
}

dim = round(dim)

if (is.null(X.names)) X.names = object$names$X
if (is.null(Y.names)) Y.names = object$names$Y

bisect = object$variates$X[, 1:dim] + object$variates$Y[, 1:dim]
cord.X = cor(object$X, bisect, use = "pairwise")
cord.Y = cor(object$Y, bisect, use = "pairwise")
sim.XY = as.matrix(cord.X %*% t(cord.Y))

result = cim(sim.XY, ...)

return(invisible(result))
}

