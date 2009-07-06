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




network.spls <-
function(object, X.names = NULL, Y.names = NULL,
        color.node = c("white", "white"),
        shape.node = c("circle", "rectangle"), color.edge = c("blue", "red"),
        lty.edge = c("solid", "solid"), lwd.edge = c(1, 1), 
        show.edge.labels = FALSE, ...) 		
{

	# validation des arguments #
	#--------------------------#

	if (length(color.node) != 2) 
		stop("'color.node' must be a vector of length.")

	if (length(shape.node) != 2) 
		stop("'shape.node' must be a vector of length.")

	if (length(color.edge) != 2) 
		stop("'color.edge' must be a vector of length.")
		
	if (length(lty.edge) != 2) 
		stop("'lty.edge' must be a vector of length.")

	if (length(lwd.edge) != 2) 
		stop("'lwd.edge' must be a vector of length.")

	p = ncol(object$X)
	q = ncol(object$Y)

	if (is.null(X.names)) X.names = object$names$X
	if (is.null(Y.names)) Y.names = object$names$Y

	# Calcul de la matrice des associations entre les variables X et Y #
	#------------------------------------------------------------------#
	simMat = predict(object, object$X[1, ])$B.hat[, , object$ncomp]
	simMat = as.vector(t(simMat))

	# Définition des sommets #
	#------------------------#	
	nodes = data.frame(name = c(X.names, Y.names), group = c(rep("x", p), rep("y", q)))

	node.X = rep(X.names, each = q)
	node.Y = rep(Y.names, p)

	# Définition des arêtes #
	#-----------------------#
	relations = data.frame(from = node.X, to = node.Y, weight = simMat)

	# Décide quels sont les arêtes à incluir dans le réseau #
	#-------------------------------------------------------#
	idx = (simMat != 0)
	relations = relations[idx, ]

	# Génère un graphe avec toutes les arêtes signifiantes #
	#------------------------------------------------------#
	gR = graph.data.frame(relations, directed = FALSE, vertices = nodes)
	
	# Attributs des sommets #
	#-----------------------#
	V(gR)$label = V(gR)$name
	
	V(gR)$label.color = "black"
	
	V(gR)$color = color.node[1]
	V(gR)$color[V(gR)$group == "y"] = color.node[2]

	V(gR)$shape = shape.node[1]
	V(gR)$shape[V(gR)$group == "y"] = shape.node[2]
	
	# Attributs des arêtes #
	#----------------------#
	if (show.edge.labels) E(gR)$label = round(E(gR)$weight, 2)
	
	E(gR)$label.color = "black"
	
	E(gR)$color = color.edge[1]
	E(gR)$color[E(gR)$weight < 0] = color.edge[2]
	
	E(gR)$lty = lty.edge[1]
	E(gR)$lty[E(gR)$weight < 0] = lty.edge[2]
	
	E(gR)$width = lwd.edge[1]
	E(gR)$width[E(gR)$weight < 0] = lwd.edge[2]
	
	gR = delete.vertices(gR, which(degree(gR) == 0) - 1)

	#----------------------------------#
	# Construction du graphe de départ #
	#----------------------------------#
	nn = vcount(gR)
	V(gR)$label.cex = min(2/log(nn), 1)
	E(gR)$label.cex = min(2.25/log(nn), 1)
	cex0 = 2*V(gR)$label.cex
	
	op = par(no.readonly = TRUE)
	
	par(pty = "s", mar = c(0, 0, 0, 0))
	plot(1:100, 1:100, type = "n", xaxt = "n")
	cha = V(gR)$label
	cha = paste(" ", cha, " ")
	xh = strwidth(cha, cex = cex0)
	yh = strheight(cha, cex = cex0) * 5/2.75
	dev.off()
	
	V(gR)$size = xh
	V(gR)$size2 = yh
	
	if (nn < 40 ) {
		l = layout.fruchterman.reingold(gR)
	}
	else {
		weights = apply(cbind(abs(E(gR)$weight)/max(abs(E(gR)$weight)), 
						rep(0.6, length(E(gR)$weight))), 1, max)	
		l = layout.fruchterman.reingold(gR, weights = weights)
	}
		
	par(pty = "s", mar = c(0, 0, 0, 0))
	plot(gR, layout = l)

	par(op)
	return(invisible(gR))
}	
		
