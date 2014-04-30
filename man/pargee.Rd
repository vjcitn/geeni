\name{pargee}
\alias{pargee}
\alias{pargee,gdManager-method}
\alias{pargee,ANY-method}

\title{Parallel GEE solver}

\description{
  A parallelized, out-of-memory based GEE solver
}

\usage{
  \S4method{pargee}{gdManager}(gd, family, binit, maxit=20, 
            tol=1e-06, sandwich=TRUE, ...)
}

\arguments{
  \item{gd}{
    A \code{\link{gdManager-class}}.
  }
  \item{family}{
    A GLM stats \code{family} instance, e.g., \code{\link{gaussian}}. See
    ?\code{family} for details.
  }
  \item{binit}{
    Initial values for the regression parameter vector.
  }
  \item{maxit}{
    Maximum number of iterations allowed.
  }
  \item{tol}{
    Convergence criterion. The absolute value of the incremental change in
    parameter estimates is measured in each iteration step. Iterations will 
    stop when the largest of these changes is smaller than \code{tol}.
  }
  \item{sandwich}{
    Logical indicating if the robust variance should be computed.
  }
  \item{\dots}{
  Arguments passed to other methods.
  }
}

\details{
  Preliminary demonstration, mildly validated on orthodontic data for 
  gaussian case and respiratory data for binomal case.
}

\value{
  A \code{list} with two named elements: `coefficients` and `robust.variance`.
}

\author{
  Vince Carey <stvjc@channing.harvard.edu>
}

\seealso{
}

\examples{

## Create a 'groupedData' object from 'Orthondont' data set:
library(nlme)
  orthGD <-
       groupedData(distance ~ age | Subject,
                   data = as.data.frame(Orthodont),
                   FUN = mean,
                   outer = ~ Sex,
                   labels = list(x = "Age",
                     y = "Distance from pituitary to pterygomaxillary fissure"),
                   units = list( x = "(yr)", y = "(mm)"))
dim(orthGD)

## Create a 'gdManager' that points to the ff files in tempdir():
orthMgr <- groupedData2ff(gd=orthGD, gcol=3, prefix=tempfile())

## In the ff file 'Sex' is coded as `1` for male and `2` for female. 
## We adjust the coding to `0` and `1`.
numdat(orthMgr)[,"Sex"] <- numdat(orthMgr)[,"Sex"] - 1
 
## Call the GEE solver:
res <- pargee(orthMgr, gaussian, c(0, 0, 0))

## Coefficients and standard errors:
round(res$coefficients, 4)
round(sqrt(diag(res$robust.variance)), 4)

## Compare results with those from geeglm() from the geepack package:
library(geepack)
geeglm <- geeglm(distance ~ age + Sex, data=orthGD, id=Subject)
round(summary(geeglm)$coefficients, 4)
}

\keyword{models}
\keyword{methods}
