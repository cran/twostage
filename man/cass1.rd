\name{cass1}
\alias{cass1}
\title{The CASS pilot dataset with sex as auxiliary covariate}

\description{ 

This is a pilot dataset from the Coronary Artery Surgery Study (CASS) 
and is discussed in Reilly (1996). It consists of a random sample
of 25 subjects from each stratum defined by different levels of
mortality and sex. The "first-stage" data from which this pilot
sample was chosen had 8096 observations. There are three variables
in the pilot dataset (see Data Description).
}

\format{
The dataset has 100 observations with 3 variables arranged in the
following columns: \cr

Column 1 (Mortality) \cr
The operative mortality of the patients \cr
(0 = Alive, 1 = Dead) \cr

Column 2 (age) \cr

Column 3 (sex) \cr
(0 = Male, 1 = Female) \cr

}

\source{
	Vliestra, et.al. (1980)
}

\seealso{
\code{\link[twostage]{ms.nprev}},\code{\link[twostage]{fixed.n}},
\code{\link[twostage]{budget}},\code{\link[twostage]{precision}},
\code{\link[twostage]{cass2}},\code{\link[twostage]{coding}}
}

\references{

Reilly,M. 1996. Optimal sampling strategies for two-stage studies. 
\emph{Amer. J. Epidemiol.} \bold{143:}92-100


Vliestra,R.E.,Frye R.L., Kromnal R.A.,et.al.(1980). Risk factors and angiographic
coronary artery disease: a report from the Coronary Artery Surgery Study (CASS).
\emph{Circulation}\bold{ 62:254-61}

}

\keyword{datasets}






