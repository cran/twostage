\name{fixed.n}
\alias{fixed.n}
\title{Optimal second stage sampling fractions, subject to
	fixed sample sizes at the first and second stage}



\description{
Optimal second stage sampling fractions (and sample sizes) using
mean score method in \bold{logistic regression setting}, 
based on first-stage sample sizes and pilot second-stage data as input.

Optimality is with respect to the standard error of a coefficient 
of interest, specified in the call to the function. \cr

\bold{BACKGROUND} \cr

This function gives the optimal second stage sampling fractions
(and sample sizes) for applications where a first-stage sample
(of size n) has already been gathered and the size of sample to
be gathered at the second stage is also fixed. Such a situation
might arise where outcome data (Y) and some covariates (Z) are
available on a database, and it is decided to pursue additional
variables on a subsample of subjects, where the size of the 
subsample is determined by time/cost considerations (an example 
would be the testing of stored bloods for a new marker which has 
been discovered since an initial case-control study was done). \cr
Since the first-stage data is available, the count (or proportion)
of first-stage observations in each (Z,Y) stratum can be computed,
and one of these vectors must be provided in the call to the "fixed.n" 
function. \cr
The optimal second-stage sampling fractions can also be found
for the situation where the first-stage data is NOT available
provided we specify the ratio of second stage sample size 
to first stage sample size (i.e the overall sampling fraction
at the second stage), and estimates of prevalences of the
(Z,Y) strata in the population. However, this situation is 
likely to be rare compared to the first scenario above.\cr

Before running the \code{fixed.n} function you should run the \code{coding}
function, to see in which order you must supply the vector of
prevalences. (see help (\code{\link[twostage]{coding}}) for details) 
}

\usage{

fixed.n (x=x,y=y,z=z,factor=NULL,n2=n2,var="var",n1="option",prev="option",frac="option")

}

\arguments{

REQUIRED ARGUMENTS

\item{x}{matrix of predictor variables} 
\item{y}{response variable (binary 0-1)}
\item{z}{matrix of the first stage variables which must be categorical
	  (can be more than one column)}
\item{n2}{size of second stage sample}
\item{var}{The name of the predictor variable whose coefficient is to be optimised. 
	    See \bold{DETAILS} if this is a factor variable \cr	 
\bold{and one of the following:}}
\item{n1}{vector of the first stage sample sizes for each (y,z) stratum \cr

OR}

\item{prev}{vector of estimated prevalences for each 
    	  (y,z) stratum, AND}
\item{frac}{the second stage sampling fraction i.e., the ratio of second stage sample 
		size to first stage sample size 
(NOTE: if \code{prev} is given, \code{frac} will also be required) \cr


OPTIONAL ARGUMENTS}
\item{factor}{the names of any factor variables in the predictor matrix}
}

\value{

\bold{A list called \code{design} consisting of the following items:}

\item{ylevel}{the different levels of response variable}
\item{zlevel}{the different levels of first stage variables z.}
\item{n1}{the first stage sample size for each (\code{ylevel},\code{zlevel}) stratum}
\item{n2}{the sample size of pilot observations for each (\code{ylevel},\code{zlevel}) stratum}
\item{prop}{optimal 2nd stage sampling proportion for each (\code{ylevel},\code{zlevel}) stratum}
\item{samp.2nd}{optimal 2nd stage sample size for each (\code{ylevel},\code{zlevel}) stratum \cr

\bold{and a list called \code{se} containing:}}
	
\item{se}{the standard errors of estimates achieved by the optimal design.}
}

\details{

	The response, predictor and first stage variables 
	have to be numeric. If you have multiple columns of 
	z, say (z1,z2,..zn), these will be recoded into
      a single vector \code{new.z}. These \code{new.z} values are
	reported as \code{zlevel} in the output (see \code{value}).


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
	then (1,1,1) will be coded as 7. \cr

	If you wish to optimise the coefficient of a factor variable, 
      you need to specify which level of the variable to optimise. 
	For example, if "weight" is a factor variable with 3 categories
	1,2 and 3 then var="weight2" will optimise the estimation of the
	coefficient which measures the difference between weight=2 and
	the baseline (weight=1). By default the baseline is always 
	the category with the smallest value. \cr
}

\examples{
\dontrun{This example of computing second stage sampling fractions subject
to a fixed total second-stage sample size uses the CASS data 
(Reilly, 1996). Once the TWOSTAGE library has been attached,
this data can be made available by:
}
data(cass1)

\dontrun{and a detailed description of the data can be obtained by} 

help (cass1) 

\dontrun{In this example, we suppose that the CASS registry only has available
the mortality(Y) and sex(Z) for the 8096 "first-stage" subjects. The pilot
data consists of 25 observations from each (Y,Z) stratum, where the sizes of
the strata are (see Reilly 1996):
	Y	Z	N
	0	0	6666
	0	1	1228
	1	0	144
	1	1	58
We wish to use this pilot information to compute the optimal design to 
minimise the variance of the sex coefficient in a logistic model 
with Sex and Age as predictors . Assume that we wish to sample a total
of 1000 subjects at the second stage.

The following commands give the output below:
}

data(cass1)
y_cass1[,1]     #--- the response variable is mortality
z_cass1[,3]     #--- the auxiliary variable is sex
x_cass1[,2:3] #--- the variables in the model are sex and age

# run CODING function to see in which order we should enter n1
coding(x=x,y=y,z=z)	
#supplying the first stage sample sizes
n1_c(6666, 1228, 144, 58)
 
# variable to be optimised (in our case sex)
fixed.n(x=x,y=y,z=z,n1=n1,var="sex",n2=1000)

\dontrun{
will give us the following output 
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
$design
     ylevel zlevel   n1 n2   prop samp.2nd
[1,]      0      0 6666 25 0.1128      752
[2,]      0      1 1228 25 0.0375       46
[3,]      1      0  144 25 1.0000      144
[4,]      1      1   58 25 1.0000       58

$se
                  [,1]
(Intercept) 0.55496070
age         0.00956422
sex         0.16472156
}
}

\seealso{
\code{\link[twostage]{ms.nprev}},\code{\link[twostage]{budget}},
\code{\link[twostage]{precision}},\code{\link[twostage]{cass1}},
\code{\link[twostage]{cass2}},\code{\link[twostage]{coding}}
}


\references{
	Reilly,M and M.S. Pepe. 1995. A mean score method for 
	missing and auxiliary covariate data in 
	regression models. \emph{Biometrika} \bold{82:}299-314 \cr

	Reilly,M. 1996. Optimal sampling strategies for 
		two-stage studies. \emph{Amer. J. Epidemiol.} 
		\bold{143:}92-100

}

\keyword{design}
