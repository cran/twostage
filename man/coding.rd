
\name{coding}
\alias{coding}
\title{combines two or more surrogate/auxiliary variables into a vector} 


\description{

   recodes a matrix of categorical variables into a vector which takes 
   a unique value for each combination \cr


\bold{BACKGROUND}

From the matrix Z of first-stage covariates, this function creates 
a vector which takes a unique value for each combination as follows:

	\tabular{rrrr}{
	z1 \tab z2 \tab z3 \tab new.z \cr
	0 \tab 0 \tab 0 \tab 1 \cr
	1 \tab 0 \tab 0 \tab 2 \cr
	0 \tab 1 \tab 0 \tab 3 \cr
	1 \tab 1 \tab 0 \tab 4 \cr
	0 \tab 0 \tab 1 \tab 5 \cr
	1 \tab 0 \tab 1 \tab 6 \cr
	0 \tab 1 \tab 1 \tab 7 \cr
	1 \tab 1 \tab 1 \tab 8 \cr
	}

If some of the combinations do not exist, the function will adjust
accordingly: for example if the combination (0,1,1) is absent above,
then (1,1,1) will be coded as 7. \cr

The values of this new.z are reported as \code{new.z} in the printed output 
(see \code{value} below) \cr

This function should be run on second stage data prior to using
the ms.nprev function, as it illustrates the order in which the call
to ms.nprev expects the first-stage sample sizes to be provided.
}


\usage{
coding(y=y,x=x,z=z,return=F)
}

\arguments{
REQUIRED ARGUMENTS

\item{y}{response variable (should be binary 0-1)}
\item{x}{matrix of predictor variables for regression model} 
\item{z}{matrix of any surrogate or auxiliary variables which must be categorical \cr


OPTIONAL ARGUMENTS}

\item{return}{logical value; if it's TRUE(T) the original surrogate
	  or auxiliary variables and the re-coded auxilliary 
	  variables will be returned.   
	  The default is False (F). 
}
}
\value{
This function does not return any values \bold{except} if \code{return}=T. \cr

If used with only second stage (i.e. complete) data, it will print the 
following:
\item{ylevel}{the distinct values (or levels) of response variable}
\item{\eqn{\bold{z}1 \dots \bold{z}i}}{the distinct values of first stage variables 
\eqn{\bold{z}1 \dots \bold{z}i}}
\item{new.z}{recoded first stage variables. Each value represents a unique combination of 
first stage variable values.}
\item{n2}{second stage sample sizes in each (\code{ylevel},\code{new.z}) stratum. \cr

If used with combined first and second stage data (i.e. with NA for 
missing values), in addition to the above items, the function will also print the following:}

\item{n1}{first-stage sample sizes in each (\code{ylevel},\code{new.z}) stratum.}

}

\examples{

\dontrun{The CASS2 data set in Reilly (1996) has 2 categorical first-stage 
variables in columns 2 (sex) and 3 (categorical weight). The predictor 
variables are  column 2 (sex) and columns 4-9 and the response variable 
is in column 1 (mort). See help(cass2) for further details. 

The commands}
data(cass2)
coding(x=cass2[,c(2,4:9)],y=cass2[,1], z=cass2[,2:3])

\dontrun{give the following coding scheme and first-stage and second-stage 
sample sizes (n1 and n2 respectively)

[1] "For calls requiring n1 or prev as input, use the following order"
      ylevel sex wtcat new.z n2
 [1,]      0   0     1     1 10
 [2,]      0   1     1     2 10
 [3,]      0   0     2     3 10
 [4,]      0   1     2     4 10
 [5,]      0   0     3     5 10
 [6,]      0   1     3     6 10
 [7,]      1   0     1     1  8
 [8,]      1   1     1     2 10
 [9,]      1   0     2     3 10
[10,]      1   1     2     4 10
[11,]      1   0     3     5 10
[12,]      1   1     3     6 10
}
}
\references{
	Reilly,M. 1996. Optimal sampling strategies for 
		        two-stage studies. \emph{Amer. J. Epidemiol.} 
		        \bold{143:}92-100
}
\seealso{
\code{\link[twostage]{ms.nprev}},\code{\link[twostage]{fixed.n}},\code{\link[twostage]{budget}}
\code{\link[twostage]{precision}}, \code{\link[twostage]{cass1}},\code{\link[twostage]{cass2}}
}

\keyword{utilities}