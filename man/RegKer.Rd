\name{RegKer}
\alias{RegKer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Function to do kernel regression
}
\description{Performs kernel regression using the Naraya-Watson estimator with the specified kernel and bandwidth}

\usage{
RegKer(x,X,Y,ker,h)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the point (or points) where we want to estimate Y value}
  \item{X}{sample of X}
  \item{Y}{sample of Y}
  \item{Ker}{kernel function to use}
  \item{h}{bandwith}
}
\details{
}
\value{
vector containing the predictions for x
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Pablo Ramsés Alonso Martín
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as

{x1<-rnorm(10,0,1); y1<-exp(x1) + rnorm(100,0,0.2)
RegKer(0,x1,y1,function(x){x^2},0.1)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
