\name{ms.nprev}
\alias{ms.nprev}
\title{Logistic regression of two-stage data using second stage sample 
       and first stage sample sizes or proportions (prevalences) as input}
	   

\description{Weighted logistic regression using the Mean Score method \cr

\bold{BACKGROUND}

This algorithm will analyse the second stage data from a two-stage
design, incorporating as appropriate weights the first stage sample
sizes in each of the strata defined by the first-stage variables.
If the first-stage sample sizes are unknown, you can still get
estimates (but not standard errors) using estimated relative 
frequencies (prevalences)of the strata. To ensure that the sample
sizes or prevalences are provided in the correct order, it is 
advisable to first run the \code{\link[twostage]{coding}} function.
}

\usage{

ms.nprev(y=y,x=x,z=z,n1="option",prev="option",factor=NULL,print.all=F)
}

\arguments{

REQUIRED ARGUMENTS
\item{y}{response variable (should be binary 0-1)}
\item{x}{matrix of predictor variables for regression model} 
\item{z}{matrix of any surrogate or auxiliary variables which must be categorical, \cr 

and one of the following:}
\item{n1}{vector of the first stage sample sizes 
 for each (y,z) stratum: must be provided
 in the correct order (see \code{\link[twostage]{coding}} function) \cr
OR}

\item{prev}{vector of the first-stage or population
 	  proportions (prevalences) for each (y,z) stratum:
          must be provided in the correct order 
          (see \code{\link[twostage]{coding}} function) \cr 
	  

OPTIONAL ARGUMENTS}

\item{print.all}{logical value determining all output to be printed. 
		 The default is False (F).} 
\item{factor}{factor variables; if the columns of the matrix of
	  predictor variables have names, supply these names, 
	  otherwise supply the column numbers. MS.NPREV will fit 
	  separate coefficients for each level of the factor variables.}

}

\value{

If called with \code{prev} will return only:

	  \bold{A list called \code{table} containing the following:}

\item{ylevel}{the distinct values (or levels) of y}
\item{zlevel}{the distinct values (or levels) of z}
\item{prev}{the prevalences for each (y,z) stratum}
\item{n2}{the sample sizes at the second stage in each stratum 
	  defined by (y,z) \cr

	  \bold{and a list called \code{parameters} containing:}}

\item{est}{the Mean score estimates of the coefficients in the
	  logistic regression model \cr \cr
	
If called with \code{n1} it will return:

	  \bold{a list called \code{table} containing:}}

\item{ylevel}{the distinct values (or levels) of y}
\item{zlevel}{the distinct values (or levels) of z}
\item{n1}{the sample size at the first stage in each (y,z) stratum}
\item{n2}{the sample sizes at the second stage in each stratum 
	  defined by (y,z) \cr

	  \bold{and a list called \code{parameters} containing:}}

\item{est}{the Mean score estimates of the coefficients in the
	  logistic regression model}	
\item{se}{the standard errors of the Mean Score estimates}
\item{z}{Wald statistic for each coefficient}
\item{pvalue}{2-sided p-value (H0: coeff=0) \cr \cr

If print.all=T, the following lists will also be returned:}

\item{Wzy}{the weight matrix used by the mean score algorithm,
	   for each Y,Z stratum: this will be in the same order 
	   as n1 and prev} 	
\item{varsi}{the variance of the score in each Y,Z stratum}
\item{Ihat}{the Fisher information matrix} 
}		   
               
\details{

	The response, predictor and surrogate variables 
	have to be numeric. If you have multiple columns of 
	z, say (z1,z2,..zn), these will be recoded into
      a single vector \code{new.z}

\tabular{rrrr}{
	z1 \tab z2 \tab z3 \tab new.z \cr
	0 \tab 0 \tab 0 \tab 1 \cr
	1 \tab	0 \tab 0 \tab 2 \cr
	0 \tab 1 \tab 0 \tab 3 \cr
	1 \tab 1 \tab 0 \tab 4 \cr
	0 \tab	0 \tab 1 \tab 5 \cr
	1 \tab 0 \tab 1 \tab 6 \cr
	0 \tab	1 \tab 1 \tab 7 \cr
	1 \tab 1 \tab 1 \tab 8 \cr
	}

	If some of the value combinations do not exist 
	in your data, the function will adjust accordingly. 
	For example if the combination (0,1,1) is absent,
	then (1,1,1) will be coded as 7.
}

\examples{
\dontrun{As an illustrative example, we use the CASS pilot data,"cass1"
from Reilly (1996)

Use }
data(cass1)        #to load the data
\dontrun{and}
help(cass1)        #for details


\dontrun{The first-stage sample sizes are:

	Y	Z	n
	0	0	6666
	0	1	1228
	1	0	144
	1	1	58	

An analysis of the pilot data using Mean Score}

# supply the first stage sample sizes in the correct order
n1_c(6666, 1228, 144, 58)
ms.nprev(y=cass1[,1], x=cass1[,2:3],z=cass1[,3],n1=n1)

\dontrun{gives the results:

[1] "please run coding function to see the order in which you"
[1] "must supply the first-stage sample sizes or prevalences"
[1] " Type ?coding for details!"
[1] "For calls requiring n1 or prev as input, use the following order"
     ylevel z new.z n2
[1,]      0 0     0 25
[2,]      0 1     1 25
[3,]      1 0     0 25
[4,]      1 1     1 25
[1] "Check sample sizes/prevalences"
$table
     ylevel zlevel   n1 n2
[1,]      0      0 6666 25
[2,]      0      1 1228 25
[3,]      1      0  144 25
[4,]      1      1   58 25

$parameters
                    est         se         z       pvalue
(Intercept) -5.06286163 1.46495235 -3.455991 0.0005482743
age          0.02166536 0.02584049  0.838427 0.4017909402
sex          0.67381300 0.21807878  3.089769 0.0020031236
}

\dontrun{Note that the Mean Score algorithm produces smaller 
standard errors of estimates than the complete-case
analysis, due to the additional information in the
incomplete cases.}
}

\seealso{
\code{\link[twostage]{fixed.n}},\code{\link[twostage]{budget}},\code{\link[twostage]{precision}}
\code{\link[twostage]{coding}},\code{\link[twostage]{cass1}},\code{\link[twostage]{cass2}}
}

\references{
	Reilly,M and M.S. Pepe. 1995. A mean score method for 
	missing and auxiliary covariate data in 
	regression models. \emph{Biometrika} \bold{82:}299-314 \cr

	Reilly,M. 1996. Optimal sampling strategies for 
		two-stage studies. \emph{Amer. J. Epidemiol.} 
		\bold{143:}92-100

}
	
\keyword{regression, missing data}

	
