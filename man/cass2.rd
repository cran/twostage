\name{cass2}
\alias{cass2}
\title{The CASS pilot dataset with sex and categorical weight 
	as auxiliary variables}

\description{

This is a pilot dataset from the Coronary Artery Surgery Study (CASS) 
register (Reilly, 1996). It consists of a random sample of 10
observations from each of the strata defined by different levels
of mortality, sex and weight category. The "first-stage" data from
which this pilot sample was chosen had 8096 observations. There are 
ten variables in the pilot dataset (see FORMAT).
}

\format{
The dataset has 118 observations with 10 variables arranged in the
following columns: \cr

Column 1 (Mortality) \cr
The operative mortality of the patients \cr
(0 = Alive, 1 = Dead) \cr

Column 2 (sex) \cr
(0 = Male, 1 = Female) \cr

Column 3 (Categorical Weight) \cr
1 = weight less than 60 kg \cr
2 = weight between 60-70 kg \cr
3 = weight 70 kg or more \cr


Column 4 (weight) \cr
The actual weight (in Kg) of the patients at the time of bypass surgery \cr

Column 5 (age) \cr
Age of patients at the time of bypass surgery \cr

Column 6 (Unstable angina) \cr
The angina status of the patients at the time of bypass surgery \cr
(0 = Stable, 1 = Unstable) \cr

Column 7 (Congestive Heart Failure Score (CHF) score) \cr

Column 8 (Left ventricular end diastolic blood pressure (LVEDBP)) \cr

Column 9 (Urgency of Surgery) \cr
(0 = Not urgent, 1 = Urgent) \cr

Column 10 (the auxiliary variables, Z) \cr
Different values describe the different levels of sex and weight 
category \cr

1 = Male; weight less than 60 kg \cr
2 = Male; weight between 60-70 kg \cr
3 = Male; weight 70 kg or more\cr
4 = Female; weight less than 60 kg \cr
5 = Female; weight between 60-70 kg \cr
6 = Female; weight 70 kg or more\cr

}

\source{
	Vliestra, et.al. (1980)
}

\seealso{
\code{\link[twostage]{msnprev}},\code{\link[twostage]{fixed.n}},
\code{\link[twostage]{budget}},\code{\link[twostage]{precision}},
\code{\link[twostage]{cass1}},\code{\link[twostage]{cass2}},\code{\link[twostage]{coding}}
}

\references{

Reilly,M. 1996. Optimal sampling strategies for two-stage studies. 
\emph{Amer. J. Epidemiol.} \bold{143:}92-100

Vliestra,R.E.,Frye R.L., Kromnal R.A.,et.al.(1980). Risk factors and angiographic
coronary artery disease: a report from the Coronary Artery Surgery Study (CASS).
\emph{Circulation}\bold{ 62:254-61}

}

\keyword{datasets}



