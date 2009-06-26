\name{cim}
\encoding{latin1}
\alias{cim}
\alias{cim.default}
\alias{cim.rcc}

\title{Clustered Image Maps (CIMs) ("heat maps")}

\description{
This function generates color-coded Clustered Image Maps (CIMs) ("heat maps") to 
represent "high-dimensional" data sets.
}

\usage{
\method{cim}{default}(mat, breaks, col = jet.colors, 
    labRow = NULL, labCol = NULL, symkey = TRUE, zoom = FALSE, 
    main = NULL, xlab = NULL, ylab = NULL, keysize = 1, 
    cexRow = min(1, 0.2 + 1/log10(nr)), 
    cexCol = min(1, 0.2 + 1/log10(nc)), 
    margins = c(5, 5), lhei = NULL, lwid = NULL, ...)
			
\method{cim}{rcc}(object, dim, X.names = NULL, Y.names = NULL, ...)
}

\arguments{
  \item{mat}{numeric matrix of values to be plotted.}
  \item{object}{object returned by \code{rcc}.}
  \item{dim}{the dimensionality to adequately account for the data association.}
  \item{X.names, Y.names}{character vector containing the names of \eqn{X}- and  
    \eqn{Y}-variables.}
  \item{breaks}{(optional) either a numeric vector indicating the splitting points for 
    binning \code{mat} into colors, or a integer number of break points to be used, in 
	which case the break points will be spaced equally between \code{min(mat)} 
	and \code{max(mat)}.}
  \item{col}{a character string specifying the colors function to use: 
    \code{\link{terrain.colors}}, \code{\link{topo.colors}}, 
	\code{\link{rainbow}} or similar functions. Defaults to \code{\link{jet.colors}}.}
  \item{labRow}{character vectors with row labels to use. 
	Defaults to \code{rownames(mat)}.}
  \item{labCol}{character vectors with column labels to use. 
	Defaults to \code{colnames(mat)}.}	
  \item{symkey}{boolean indicating whether the color key should be made symmetric about 0. 
	Defaults to \code{TRUE}.}
  \item{zoom}{logical. Whether to use zoom for interactively zooming-out. See Details.}
  \item{main, xlab, ylab}{main, \eqn{x}- and \eqn{y}-axis titles; defaults to none.}
  \item{cexRow, cexCol}{positive numbers, used as \code{cex.axis} in for the row or column 
	axis labeling. The defaults currently only use number of rows or columns, respectively.}
  \item{keysize}{numeric value indicating the size of the color key.}
  \item{margins}{numeric vector of length two containing the margins (see \code{\link{par}(mar)}) 
	for co\-lumn and row names respectively.}
  \item{lhei, lwid}{arguments passed to \code{layout} to divide the device up into two 
    rows and two columns, with the row-heights \code{lhei} and the column-widths \code{lwid}.}
  \item{\ldots}{arguments passed to \code{cim.default}.}
}

\details{ 
One matrix Clustered Image Map (default method) is a 2-dimensional visualization of a 
real-valued matrix (basically \code{\link{image}(t(mat))}) with a dendrogram added to the 
left side and to the top. The rows and columns are reordered according to some hierarchical
clustering method to identify interesting patterns. By default the used clustering method for
rows and columns is the \emph{complete linkage} method. 
 
In \code{cim.rcc} method, the matrix \code{mat} is created where element \eqn{(i,j)} is the
scalar product value between every pairs of vectors in dimension \code{dim} representing 
the variables \eqn{X_i} and \eqn{Y_j} on the axis defined by \eqn{(Z_1,\ldots ,Z_\code{dim})}, where 
\eqn{Z_i}, \eqn{i=1, \ldots ,}\code{dim}, is the bisector vector between the \eqn{i}-th \eqn{X} and 
\eqn{Y} canonical variate.

For visualization of "high-dimensional" data sets, a nice zooming tool was created. 
\code{zoom=TRUE} open a new device, one for CIM, one for zoom-out region and
define an interactive `zoom' process: click two points at imagen map region by pressing 
the first mouse button. It then draws a rectangle around the selected region and zoom-out 
this at new device. The process can be repeated to zoom-out other regions of interest. 

The zoom process is terminated by clicking the second button and selecting 'Stop' 
from the menu, or from the 'Stop' menu on the graphics window.
}

\value{
A object of class \code{"rcc"}, a list containing the following components:
  \item{rowInd}{row index permutation vectors as returned 
    by \code{\link{order.dendrogram}}.}
  \item{colInd}{column index permutation vectors as returned 
    by \code{\link{order.dendrogram}}.}	
  \item{ddr, ddc}{object of class \code{"dendrogram"} which describes the row and 
	column trees produced by \code{cim}.}
  \item{labRow, labCol}{character vectors with row and column labels used.}
}

\references{
Eisen, M. B., Spellman, P. T., Brown, P. O. and Botstein, D. (1998). 
Cluster analysis and display of genome-wide expression patterns. 
\emph{Proceeding of the National Academy of Sciences of the USA} \bold{95}, 14863-14868.

Weinstein, J. N., Myers, T. G., O'Connor, P. M., Friend, S. H., Fornace Jr., A. J., 
Kohn, K. W., Fojo, T., Bates, S. E., Rubinstein, L. V., Anderson, N. L., 
Buolamwini, J. K., van Osdol, W. W., Monks, A. P., Scudiero, D. A., Sausville, E. A., 
Zaharevitz, D. W., Bunow, B., Viswanadhan, V. N., Johnson, G. S., Wittes, R. E. 
and Paull, K. D. (1997). An information-intensive approach to the molecular 
pharmacology of cancer. \emph{Science} \bold{275}, 343-349.
}

\author{Sébastien Déjean and Ignacio González.}

\seealso{\code{\link{image}}, \code{\link{heatmap}}, \code{\link{hclust}}.}

\examples{
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
nutri.res <- rcc(X, Y, lambda1 = 0.064, lambda2 = 0.008)

dends <- cim(nutri.res, dim = 3, xlab = "genes", 
             ylab = "lipids", margins = c(5, 6))		
plot(dends$ddr, axes = FALSE)

## interactive 'zoom' 
\dontrun{
cim(nutri.res, dim = 3, zoom = TRUE)
## select the region and "see" the zoom-out region
}
}

\keyword{multivariate}
\keyword{iplot}
\keyword{hplot}
\keyword{graphs}
\keyword{cluster}

