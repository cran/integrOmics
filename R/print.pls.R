'print.pls' <-
function(x, ...){

# -------------print pls ----------------------------------------
#if(any(class(x)=="pls")){
  cat("\nCall:\n", deparse(x$call), "\n")
  cat(" PLS with a", x$mode, "mode with", x$ncomp, "PLS components"  ,"\n",sep=" ")
  cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n",sep=" ")
  cat(" You entered data Y of dimensions:", nrow(x$Y), ncol(x$Y), "\n \n",sep=" ")

  cat(" No variable selection", "\n \n",sep=" ")
  cat(" Available components: latent variables, loading vectors and names for each data set", "\n",sep=" ")

#}
}


# ----------------print spls ----------------------------------

'print.spls' <-
function(x, ...){

#if(any(class(liver.spls) ==  c('spls'))){
  cat("\nCall:\n", deparse(x$call), "\n")
  cat(" sPLS with a", x$mode, "mode with", x$ncomp, "sPLS components"  ,"\n",sep=" ")
  cat(" You entered data X of dimensions :", nrow(x$X), ncol(x$X), "\n",sep=" ")
  cat(" You entered data Y of dimensions :", nrow(x$Y), ncol(x$Y), "\n \n",sep=" ")
  cat(" Selection of", x$keepX, "variables on each of the sPLS components on the X data set", "\n",sep=" ")
  cat(" Selection of", x$keepY, "variables on each of the sPLS components on the Y data set", "\n \n",sep=" ")

  cat(" Available components: latent variables, loading vectors and names for each data set", "\n",sep=" ")
#  cat(" :","\n",sep=" ")
#  cat(  as.character(attributes(x)$names[10:12]), "\n",sep=" ")
#  cat(c("variates", "loadings", "names"), "\n",sep=" ")

#}


}




