\name{newNorm-package}
\alias{newNorm-package}
\alias{newNorm}
\docType{package}
\title{
What the package does (short line)
NewNorm
}
\description{
More about what it does (maybe more than one line)
 A concise (1-5 lines) description of the package
}
\details{
\tabular{ll}{
Package: \tab newNorm\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-04-29\cr
License: \tab What license is it under?\cr
}
An overview of how to use the package, including the most important functions
}
\author{
Who wrote it

Maintainer: Who to complain to <yourfault@somewhere.net>
The author and/or maintainer of the package 
}
\references{
Literature or other references for background information
}

\keyword{ package }

\examples{
    data("data")
    ncmp <- 4
    newnormout <- newnorm(ifquant=FALSE, sigA=sigAsample, sigB=sigBsample, Annot=Annotsample, 
                      controlred=matred, controlgrn=matgrn, 
                      cp.types=cp.types, cell_type = cell_type,
                      ncmp = ncmp, save.loess=TRUE, applyloess=TRUE, logit.quant=TRUE)
    origBeta <- newnormout[[1]]
    newBeta <- newnormout[[2]]
    rownames(newBeta) <- rownames(origBeta)
    colnames(newBeta) <- colnames(origBeta)
}
