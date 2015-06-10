\name{newnorm}
\alias{newnorm}

\title{
Function to normalize Illumina Infinium Human Methylation 450 BeadChip (Illumina 450K)
}
\description{
Function takes Illumina Infinium Human Methylation 450 BeadChip data ans performs normalization 
correcting for quantile 
}
\usage{
newnorm(sigA, sigB, Annot = NULL, quantiledat = NULL, controlred, controlgrn, cp.types = NULL, cell_type, ncmp = 4, save.quant = TRUE, save.loess = TRUE, apply.loess = TRUE, logit.quant = TRUE, validate = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sigA, sigB}{Matrices with signals A and B.

}
%%  \item{sigB}{
%%     ~~Describe \code{sigB} here~~
%%}
  \item{Annot}{
%%     ~~Describe \code{Annot} here~~
}
  \item{quantiledat}{
%%     ~~Describe \code{quantiledat} here~~
}
  \item{controlred}{
%%     ~~Describe \code{controlred} here~~
}
  \item{controlgrn}{
%%     ~~Describe \code{controlgrn} here~~
}
  \item{cp.types}{
%%     ~~Describe \code{cp.types} here~~
}
  \item{cell_type}{
%%     ~~Describe \code{cell_type} here~~
}
  \item{ncmp}{
%%     ~~Describe \code{ncmp} here~~
}
  \item{save.quant}{
%%     ~~Describe \code{save.quant} here~~
}
  \item{save.loess}{
%%     ~~Describe \code{save.loess} here~~
}
  \item{apply.loess}{
%%     ~~Describe \code{apply.loess} here~~
}
  \item{logit.quant}{
%%     ~~Describe \code{logit.quant} here~~
}
  \item{validate}{
%%     ~~Describe \code{validate} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# to normalize methylation data:
    data("data")
    ncmp <- 4
   newnormout <- newnorm(sigA=sigAsample, sigB=sigBsample, Annot=Annotsample, 
                      controlred=matred, controlgrn=matgrn, 
                      cp.types=cp.types, cell_type = cell_type,
                       ncmp=ncmp, save.quant=TRUE, save.loess=TRUE, apply.loess=TRUE, logit.quant=TRUE, validate=FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
