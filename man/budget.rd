\name{budget}
\alias{budget}
\title{Optimal sampling design for 2-stage studies with fixed budget}

\description{
Optimal design for two-stage-study with budget constraint
using the Mean Score method \cr

\bold{BACKGROUND} \cr

This function calculates the total number of study observations
and the second-stage sampling fractions that will maximise precision 
subject to an available budget. The user must also supply the unit
cost of observations at the first and second stage, and the vector
of prevalences in each of the strata defined by different levels 
of dependent variable and first stage covariates . 

Before running the \code{budget} function you should run the \code{coding}
function, to see in which order you must supply the vector of
prevalences. (see help (\code{\link[twostage]{coding}}) for details) 

}


\usage{

budget (x=x,y=y,z=z,factor=NULL,prev=prev,var="var",b=b,c1=c1,c2=c2)
}

\arguments{

REQUIRED ARGUMENTS
\item{x}{matrix of predictor variables} 
\item{y}{response variable (binary 0-1)}
\item{z}{matrix of the first stage variables which must be categorical
	 (can be more than one column)}
\item{prev}{vector of estimated prevalences for each (y,z) stratum}
\item{var}{The name of the predictor variable whose coefficient is to be optimised. 
	    See \bold{DETAILS} if this is a factor variable}	 
\item{b}{the total budget available}
\item{c1}{the cost per first stage observation}
\item{c2}{the cost per second stage observation \cr

OPTIONAL ARGUMENTS}
\item{factor}{the names of any factor variables in the predictor matrix}
}

\value{
The following lists will be returned:

\item{\bold{n}}{the optimal number of observations (first stage sample size)}

\item{\bold{se}}{the standard error of estimates achieved by the optimal design \cr

\bold{and a list called \code{design} consisting of the following items:}}

\item{ylevel}{the different levels of the response variable}
\item{zlevel}{the different levels of first stage covariates z.}
\item{prev}{the prevalence of each (\code{ylevel},\code{zlevel}) stratum}
\item{n2}{the sample size of pilot observations for each (\code{ylevel},\code{zlevel}) stratum}
\item{prop}{optimal 2nd stage sampling proportion for each (\code{ylevel},\code{zlevel}) stratum}
\item{samp.2nd}{optimal 2nd stage sample size for each (\code{ylevel},\code{zlevel}) stratum}

}

\details{
	The response, predictor and first stage variables 
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

\dontrun{We give an example using the pilot subsample from the CASS data 
discussed in Reilly(1996). The data are in the cass2 matrix, which can be loaded using}

data(cass2)

\dontrun{and a description of the dataset can be seen using}

help(cass2)

\dontrun{In our examples below, we use sex and weight as auxiliary variables.

Given an available budget of £10,000, a first-stage cost of £ 1/unit
and second-stage cost £ 0.5/unit, the codes below will calculate
the sampling strategy that optimises the precision of the 
coefficient for SEX : see output below.}

data(cass2)
y_cass2[,1]            #response variable
z_cass2[,10]           #auxiliary variable
x_cass2[,c(2,4:9)]     #predictor variables

# run CODING function to see in which order we should enter prevalences
coding(x=x,y=y,z=z)	
# supplying the prevalence (from Table 5, Reilly 1996)

prev_c(0.0197823937,0.1339020772,0.6698813056,0.0544015826,
+ 0.0503214639,0.0467359050,0.0009891197,0.0040801187,0.0127349159,
+ 0.0022255193,0.0032146390,0.0017309594)

# optimise sex coefficient

budget(x=x,y=y,z=z,var="sex",prev=prev,b=10000,c1=1,c2=0.5)


\dontrun{OUTPUT

[1] "please run coding function to see the order in which you"
[1] "must supply the first-stage sample sizes or prevalences"
[1] " Type ?coding for details!"
[1] "For calls requiring n1 or prev as input, use the following order"
      ylevel z new.z n2
 [1,]      0 1     1 10
 [2,]      0 2     2 10
 [3,]      0 3     3 10
 [4,]      0 4     4 10
 [5,]      0 5     5 10
 [6,]      0 6     6 10
 [7,]      1 1     1  8
 [8,]      1 2     2 10
 [9,]      1 3     3 10
[10,]      1 4     4 10
[11,]      1 5     5 10
[12,]      1 6     6 10
[1] "Check sample sizes/prevalences"
$n
[1] 9166

$design
      ylevel zlevel         prev n2   prop samp.2nd
 [1,]      0      1 0.0197823937 10 0.5230       95
 [2,]      0      2 0.1339020772 10 0.2841      349
 [3,]      0      3 0.6698813056 10 0.0726      446
 [4,]      0      4 0.0544015826 10 0.4488      224
 [5,]      0      5 0.0503214639 10 0.2480      114
 [6,]      0      6 0.0467359050 10 0.4922      211
 [7,]      1      1 0.0009891197  8 1.0000        9
 [8,]      1      2 0.0040801187 10 1.0000       37
 [9,]      1      3 0.0127349159 10 1.0000      117
[10,]      1      4 0.0022255193 10 1.0000       20
[11,]      1      5 0.0032146390 10 1.0000       29
[12,]      1      6 0.0017309594 10 1.0000       16

$se
                   [,1]
(Intercept) 1.193504705
sex         0.217235702
weight      0.006718422
age         0.014588813
angina      0.245831383
chf         0.077039239
lve         0.010071151
surg        0.179887419}
}
\seealso{
\code{\link[twostage]{ms.nprev}},\code{\link[twostage]{fixed.n}},
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
