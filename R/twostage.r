
fixed.n<-function(x=x,y=y,z=z,factor=NULL,n2=n2,var=NULL,n1="option",prev="option",frac="option")
{
 #--function to calculate the optimal second stage 
 #--sample once you already get your first stage sample and want
 #--to sample certain amount of them for the second stage. 

 #--The objective is to minimize the variance of kth variable
 #--on the logistic regression model, therefore your response
 #--variable have to be binary.

 Wzykk<-1 

 #---if you only know the prevalence of the first stage sample-----
 #---you also need to suply the relative fraction of second stage
 #---sample to the first stage sample

 chkvar<-is.character(var)
 if(!chkvar)
	stop("Please enter the variable name to be optimised!")

 #--if you know the prevalences -----------
 if (prev[1]!="option"){

 	#--- check the prevalences
 	if (abs(sum(prev)-1)>0.01) 
		stop("Please check the prevalences it does not some up to 1")

 	out<-ms.nprev(x=x,y=y,z=z,factor=factor,prev=prev,print.all=TRUE)
 	table<-out$table
 	label<-rownames(out$parameters)

 	Wzy<-out$Wzy
 	invI<-solve(out$Ihat)

 	nvar<-ncol(invI)
 	k<-match(var,label)

 	if(is.na(k))
		stop(paste(var,"IS NOT FOUND, PLS CHECK THE VAR NAME"),call.=FALSE)

 	for (i in 1:length(prev))
		Wzykk[i]<-Wzy[k,k,i]


 	numer<-frac*sqrt(Wzykk)
 	denom<-sum(prev*sqrt(Wzykk))
 	prop<-numer/denom

 	#-- if some of the strata have zero sampling fraction the program will terminate

 	if (sum(round(prop,4)==0)>0) {
		table.zero<-cbind(table,prop)
		id<-1:nrow(table)
		zrpost<-id[round(prop,4)==0]
		table.zero<-round(table.zero[zrpost,],4)
		print(table.zero)
		cat("These strata have 0 (zero) 2nd samp.fraction","\n") 
		stop("please redefined your strata OR gather more pilot obs",call.=FALSE)
  	}	

 	# correcting the proportion of second stage sample size <=100%
 	id<-1:length(prop)
 	index1<-NULL
	while (sum(prop>1)>0) {
		index1<-c(index1,id[rank(prop)==length(prop)])
		frac.1<-frac-sum(prev[index1])
		numer1<-(frac.1)*sqrt(Wzykk)
		denom1<-sum(sqrt(Wzykk[-index1])*prev[-index1])
		prop<-numer1/denom1
		prop[index1]<-1
	}
 	n2<-round(prop*prev*n2/frac)
 	prop<-round(prop,4)
	final.est<-cbind(prop,n2)
 	colnames(final.est)<-c("prop","samp.2nd")
	list(prop=cbind(table,final.est))
 }

#--if you know the sample size of the first stage study-----------
else {
	out<-ms.nprev(x=x,y=y,z=z,factor=factor,n1=n1,print.all=TRUE)
	table<-out$table

	Wzy<-out$Wzy 
	invI<-solve(out$Ihat)
	label<-rownames(out$parameters)
	nvar<-ncol(invI)
	k<-match(var,label)

	if(is.na(k))
		stop(paste(var,"IS NOT FOUND, PLS CHECK THE VAR NAME"),call.=FALSE)

	for (i in 1:length(n1))
		Wzykk[i]<-Wzy[k,k,i]

	prev<-n1/sum(n1)
	numer<-(n2/sum(n1))*sqrt(Wzykk)
	denom<-sum(prev*sqrt(Wzykk))
	prop<-numer/denom

#-- if some of the strata have zero sampling fraction the program will terminate

	if (sum(round(prop,4)==0)>0) {
		table.zero<-cbind(table,prop)
		id<-1:nrow(table)
		zrpost<-id[round(prop,4)==0]
		table.zero<-round(table.zero[zrpost,],4)
		print(table.zero)
		cat("These strata have 0 (zero) 2nd samp.fraction","\n") 
		stop("please redefined your strata OR gather more pilot obs",call.=FALSE)
	}	

# correcting the proportion of second stage sample size <=100%
	id<-1:length(prop)
	index1<-NULL
	while (sum(prop>1)>0){
		index1<-c(index1,id[rank(prop)==length(prop)])
		n2.1<-n2-sum(n1[index1])
		numer1<-(n2.1/sum(n1))*sqrt(Wzykk)
		denom1<-sum(sqrt(Wzykk[-index1])*prev[-index1])
		prop<-numer1/denom1
		prop[index1]<-1
	}

#---calculate the minimum variance obtained by the design
	wgt<-prev/prop*(1-prop)
	result <- 0
	varsi<-ar<-array(0,dim=c(nrow(invI),nrow(invI),length(wgt)))
	for(i in 1:length(wgt)) {
		varsi <-Wzy[,,i]
		result <- result + varsi * wgt[i]
		}
	
	Vhat <- result	
	
	V <- (invI + Vhat)/sum(n1)


#--the function will return the number of sample you should
#--take from each stratum for the second stage in order to
#--minimize the variance

	n2<-round(n1*prop)
	prop <- round(prop,4)
	final.est<-cbind(prop,n2)
	colnames(final.est)<-c("prop","samp.2nd")
	se<-as.matrix(sqrt(diag(V)))
	rownames(se)<-label
	list(design=cbind(table,final.est),se=se)
 }
# close function
}

budget<-
function(x="covariate",y="response",z="auxiliary",factor=NULL,prev=prev,var=NULL,b="budget", c1="first cost",c2="second cost")

{
 #---function to compute optimal design for two-stage-studies
 #---subject to budget restriction. The function will calculate    
 #   optimal' number of sample  
 #---and the second stage proportion in order to minimize the variance  
 #   of the kth variable (the variable we want to optimise) on the    
 #   logistic regression model.

 chkvar<-is.character(var)
 if(!chkvar)
	stop("Please enter the variable name to be optimised!")

 Wzykk<-cost<-1
 out<-ms.nprev(x=x,z=z,y=y,factor=factor,prev=prev,print.all=TRUE)
 table<-out$table
 label<-rownames(out$parameters)

 Wzy<-out$Wzy
 n2<-out$table[,4]
 invI<-solve(out$Ihat)
 nvar<-ncol(invI)
 k<-match(var,label)
 if(is.na(k))
	stop(paste(var,"IS NOT FOUND, PLS CHECK THE VAR NAME"),call.=FALSE)

 #--- check the prevalences
 if (abs(sum(prev)-1)>0.01) 
	stop("Please check the prevalences it does not some up to 1")

 for (i in 1:length(prev))
	Wzykk[i]<-Wzy[k,k,i]

 id<-1:length(prev)
 index1<-NULL
 subs<-Wzykk
 denom<-invI[k,k]-sum(prev*Wzykk)

 #---check if the denominator is less than zero (we will sample 100%
 #---from stratum causing the denominator be negative. The optimization
 #---will take place on the rest of the strata 
 while (denom<=0){
	index1<-c(index1,id[rank(subs[id])==(length(subs)-length(index1))])
	denom<-invI[k,k]-sum(prev[-index1]*Wzykk[-index1])
	}

 if (!is.null(index1)){
	above<-sqrt((c1+sum(prev[index1]*c2))*c2)*sum(prev[-index1]*sqrt(Wzykk)[-index1])
	below<-sqrt(invI[k,k]-sum(prev[-index1]*Wzykk[-index1]))
	n<-b*((c1+sum(prev[index1]*c2))+above/below)^(-1)
	prop<-(b-n*(c1+sum(prev[index1]*c2)))/(n*c2)*(sqrt(Wzykk)/sum(prev[-index1]*sqrt(Wzykk[-index1])))
	prop[index1]<-1
	}

 else {
	above<-sqrt(c1*c2)*sum(prev*sqrt(Wzykk))
	below<-sqrt(invI[k,k]-sum(prev*Wzykk))
	n<-b*(c1+above/below)^(-1)
	prop<-(b-n*c1)/(n*c2)*(sqrt(Wzykk)/sum(prev*sqrt(Wzykk)))
	}

 #-- if some of the strata have zero sampling fraction the program will terminate

 if (sum(round(prop,4)==0)>0) {
	table.zero<-cbind(table,prop)
	id<-1:nrow(table)
	zrpost<-id[round(prop,4)==0]
	table.zero<-round(table.zero[zrpost,],4)
	print(table.zero)
	cat("These strata have 0 (zero) 2nd samp.fraction","\n") 
	stop("please redefined your strata OR gather more pilot obs",call.=FALSE)
	}	

 #---correcting the proportion of second stage sample size <=100%
 while (sum(prop>1)>0){
	index1<-c(index1,id[rank(prop)==length(prop)])
	above<-sqrt((c1+sum(prev[index1]*c2))*c2)*sum(prev[-index1]*sqrt(Wzykk)[-index1])
	below<-sqrt(invI[k,k]-sum(prev[-index1]*Wzykk[-index1]))
	n<-b*((c1+sum(prev[index1]*c2))+above/below)^(-1)
	prop<-(b-n*(c1+sum(prev[index1]*c2)))/(n*c2)*(sqrt(Wzykk)/sum(prev[-index1]*sqrt(Wzykk[-index1])))
	prop[index1]<-1
	cost<-(c1*round(n)+sum(round((prop*prev*round(n)*c2))))
 	}

 #---calculate the minimum variance obtained by the design
 wgt<-prev/prop*(1-prop)
 result <- 0
 varsi<-ar<-array(0,dim=c(nrow(invI),nrow(invI),length(wgt)))
 for(i in 1:length(wgt)) {
	varsi <-Wzy[,,i]
	result <- result + varsi * wgt[i]
	}
 Vhat <- result
 V <- (invI + Vhat)/n

 #---the function will return the total number of sample and the second stage proportion
 #---for each stratum to achieve minimum variance.
 n<-round(n)
 n2<-round(prev*prop*n)
 prop<-round(prop,4)
 final.est<-cbind(prop,n2)
 colnames(final.est)<-c("prop","samp.2nd")
 se<-as.matrix(sqrt(diag(V)))
 rownames(se)<-label

 return(n=n,design=cbind(table,final.est),se=se)
 }


precision<-function(x=x,y=y,z=z,factor=NULL,var=NULL,prev=prev,prc="precision",c1="first cost",c2="second cost")

{
 #---functionto compute optimal sampling design for two-stage-studies
 #---subject to the fixed variance of the kth variable. The function will calculate 'optimal' 
 #---number of sample and the second stage proportion in order to minimize the cost
 #---to achieve the wanted variance.

 Wzykk<-cost<-1

 chkvar<-is.character(var)
 if(!chkvar)
	stop("Please enter the variable name to be optimised!")

 out<-ms.nprev(x=x,z=z,y=y,factor=factor,prev=prev,print.all=TRUE)
 table<-out$table
 label<-rownames(out$parameters)

 Wzy<-out$Wzy
 n2<-out$table[,4]
 invI<-solve(out$Ihat)

 nvar<-ncol(invI)
 k<-match(var,label)
 if(is.na(k))
	stop(paste(var,"IS NOT FOUND, PLS CHECK THE VAR NAME"),call.=FALSE)

 #--- check the prevalences

 if (abs(sum(prev)-1)>0.01) 
	stop("Please check the prevalences it does not some up to 1")

 for (i in 1:length(prev))
	Wzykk[i]<-Wzy[k,k,i]

 id<-1:length(prev)
 index1<-NULL
 subs<-prev*Wzykk
 denom<-invI[k,k]-sum(subs)

 #--check if the denominator is less than zero (we decide to take 100% of
 #--the first stage sample from the stratum causing negative denominator

 while (denom<=0){
	index1<-c(index1,id[rank(subs[id])==(length(subs)-length(index1))])
	denom<-invI[k,k]-sum(subs[-index1])
	}

 # if one or more strata has a big value of prev*Wzykk, we will
 # sample all strata members for the second stage

 if (!is.null(index1)){
	n<-denom/prc+sqrt(c2/(c1+sum(prev[index1]*c2)))/prc*sum(prev[-index1]*
	sqrt(Wzykk[-index1])*sqrt(denom))	
	prop<-sqrt((c1+sum(prev[index1]*c2))/c2)*sqrt(Wzykk/denom)
	prop[index1]<-1
	}

 else {
	n<-denom/prc+sqrt(c2/c1)/prc*sum(prev*
	sqrt(Wzykk)*sqrt(denom))	
	prop<-sqrt(c1/c2)*sqrt(Wzykk/denom)
	prop[index1]<-1
	}

 #-- if some of the strata have zero sampling fraction the program will terminate
 if (sum(round(prop,4)==0)>0) {
	table.zero<-cbind(table,prop)
	id<-1:nrow(table)
	zrpost<-id[round(prop,4)==0]
	table.zero<-round(table.zero[zrpost,],4)
	print(table.zero)
	cat("These strata have 0 (zero) 2nd samp.fraction","\n") 
	stop("please redefined your strata OR gather more pilot obs",call.=FALSE)
 	}	

# correcting the proportion of second stage sample size <=100%
 while (sum(prop>1)>0){
	index1<-c(index1,id[rank(prop)==length(prop)])
	denom<-invI[k,k]-sum(subs[-index1])
	n<-denom/prc+sqrt(c2/(c1+sum(prev[index1]*c2)))/prc*sum(prev[-index1]*
	sqrt(Wzykk[-index1])*sqrt(denom))	
	prop<-sqrt((c1+sum(prev[index1]*c2))/c2)*sqrt(Wzykk/denom)
	prop[index1]<-1
	}

 n<-round(n)
 n2<-round(prev*prop*n)
 prop<-round(prop,4)

 wgt<-prev/prop*(1-prop)
 result <- 0
 varsi<-ar<-array(0,dim=c(nrow(invI),nrow(invI),length(wgt)))
 for(i in 1:length(wgt)) {
	varsi <-Wzy[,,i]
	result <- result + varsi * wgt[i]
	}
 Vhat <- result	
 V <- (invI + Vhat)/n
 cost<-(c1*round(n)+sum(round((prop*prev*round(n)*c2))))

 #---The function will return the total number of sample, the proportion of
 #---the second stage sample and the minimum cost for achieving the
 #---wanted variance.
 final.est<-cbind(prop,n2)
 colnames(final.est)<-c("prop","samp.2nd")
 var<-as.matrix(diag(V))
 rownames(var)<-label

 return(n=n,design=cbind(table,final.est),cost=cost,var=var)
 }


#   @@@@@@@@@@@@@@@@@@@@@@@@@   MS.NPREV  FUNCTION  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
"ms.nprev"<-function(x="complete data",y=y,z=z,n1="option",prev="option",factor=NULL,print.all=FALSE)
{
# this function uses the second-stage (i.e. complete ) data and the first-stage
# sample sizes (or prevalences) to compute Mean Score estimates of the coefficients
# in a logistic regression model. The function requires the following input:
#
# x= the matrix of predictor variables in the regression model
#    defined as a data.frame (before calling the function)
#
# y= the outcome variable vector
#
# z= the surrogate variable vectors, defined as data.frame
#
# n1 OR prev where
#              n1=  the vector of first-stage sample sizes
#                   in the Y,Z strata in the same order as
#                   given by the table(y,z) command i.e.sorted by Y and Z(within Y)
#             prev= prevalence of the Y,Z strata in the same
#                   order as specified for n1
#
#
# The function called with "prev" returns only:
#          ylevel=  the distinct values (or levels) of y
#          zlevel=  the distinct values (or levels) of z
#              n2=  the sample sizes at the second stage at each stratum 
#                   defined by (ylevel,zlevel)
#             est=  the mean score estimates
#  and if called with n1 also returns:
#              se=  the standard errors of the MS estimates
#             Wzy=  the Wzy matrix  for each Y,Z stratum
#		    in the same order as n1 and prev 	
#           varsi=  the variance of score in each Y,Z stratum
#            Ihat=  the estimated information matrix
#              n2=  the second-stage sample sizes in each (Z,Y) stratum
#
#
 print("please run coding function to see the order in which you")
 print("must supply the first-stage sample sizes or prevalences")
 print (" Type ?coding for details!")

 stop1<-c("ARE NOT FOUND PLEASE CHECK COL NAMES OR ENTER COL NUMBER IN THE PREDICTOR MATRIX")         
 z1<-data.frame(z)
 z.old<-as.matrix(z) 

 z<-coding(x=x,y=y,z=z,return=TRUE)$z
  
 ylev<-as.numeric(levels(factor(y)))
 zlev<-as.numeric(levels(factor(z)))
 ylevel<-rep(ylev,rep(length(zlev),length(ylev)))
 zlevel<-rep(zlev,length(ylev))
 n2<-c(t(table(y,z)))

 if(min(n2)<2) {
	stop("WARNING: One or more strata with less than 2 obs!")  
	}

 w.MS <<- rep(1, length(y))
		
 if (prev[1]!="option") {
	print("Check sample sizes/prevalences")
	wt<-prev/n2
	for(i in 1:length(wt)) {
		w.MS<<- ifelse(y == ylevel[i] & z == zlevel[i], wt[i], w.MS)
	}
	rdata<-data.frame(y,x)

	# recode the factor variables
	if (length(factor) > 0) {
		if(is.character(factor)) {
		   for (i in 1:length(factor)) {
			ind<-ifelse(colnames(rdata)==factor[i],1,0)
			if (sum(ind)==0) {
			stop(paste(factor[i],stop1),call.=FALSE)
			}	
		varpost<-order(ind)[ncol(rdata)]
		ff<-factor(rdata[,varpost])
		factlev<-levels(ff)
		dummy<-as.data.frame(model.matrix(~ ff - 1)[,-1])
		colnames(dummy)<-paste(colnames(rdata)[varpost], factlev[-1], sep = "")
		rdata<-rdata[,-varpost]
		rdata<-cbind(rdata,dummy)
			
		}
 	}

	else if (is.numeric(factor)){
		varpost<-factor+1
		ind<-ifelse(varpost>ncol(rdata),1,0)
		if (sum(ind)>0) {
			stop("COLUMN NUMBER OF FACTOR VARIABLES IS OUT OF BOUND,PLEASE CHECK!")			  	
			}
		}
		
		for (i in 1:length(factor)) {
			ff<-factor(rdata[,varpost[i]])
			factlev<-levels(ff)
			dummy<-as.data.frame(model.matrix(~ ff - 1)[,-1])
			colnames(dummy)<-paste(colnames(rdata)[varpost[i]], factlev[-1], sep = "")
			rdata<-rdata[,-varpost[i]]
			rdata<-cbind(rdata,dummy)
			}
	}

	glm.MS <- glm(y ~ ., family = "binomial",weight=w.MS, data = rdata)			
	
	X <- cbind(1, rdata[,-1])
	pi.hat <- glm.MS$fitted.values
	Ihat <- (t(as.matrix(X)) %*% (as.matrix(X) * w.MS * pi.hat * (1 - pi.hat)))
	S <- X * (y - pi.hat)
	
	invI <- solve(Ihat)
	varsi<-ar<-array(0,dim=c(nrow(Ihat),nrow(Ihat),length(wt)))
	for(i in 1:length(wt)) {
		si <- S[y == ylevel[i] & z == zlevel[i],  ]
		varsi[,,i]<-var(si)
		ar[,,i]<-invI%*%varsi[,,i]%*%invI
		}

	if (print.all)
	list(table=cbind(ylevel=ylevel, zlevel=zlevel, prev=prev, n2=n2),
             parameters=cbind(est=glm.MS$coef),Wzy=ar,Ihat=Ihat,varsi=varsi)
	else list(table=cbind(ylevel=ylevel, zlevel=zlevel,prev=prev,n2=n2),
		parameters=cbind(est=glm.MS$coef))
    }
		

else 
   { 
	print("Check sample sizes/prevalences")
	N1<-sum(n1)
	prev<-n1/N1
	rdata<-data.frame(y,x)

	wt<-n1/n2
	for(i in 1:length(wt)) {
		w.MS<<- ifelse(y == ylevel[i] & z == zlevel[i], wt[i], w.MS)
				}


	# recode the factor variables
	if (length(factor) > 0) {
		
		if(is.character(factor)) {
		   for (i in 1:length(factor)) {
			ind<-ifelse(colnames(rdata)==factor[i],1,0)
			if (sum(ind)==0) {
			 stop(paste(factor[i],stop1),call.=FALSE)
			 	}	
			varpost<-order(ind)[ncol(rdata)]
			ff<-factor(rdata[,varpost])
			factlev<-levels(ff)
			dummy<-as.data.frame(model.matrix(~ ff - 1)[,-1])
			colnames(dummy)<-paste(colnames(rdata)[varpost], factlev[-1], sep = "")
			rdata<-rdata[,-varpost]
			rdata<-cbind(rdata,dummy)
			
			}
		    }

		else if (is.numeric(factor)){
		  varpost<-factor+1
		  ind<-ifelse(varpost>ncol(rdata),1,0)
		  if (sum(ind)>0) {
			stop("COLUMN NUMBER OF FACTOR VARIABLES IS OUT OF BOUND,PLEASE CHECK!")			  	}
		
		  for (i in 1:length(factor)) {
			ff<-factor(rdata[,varpost[i]])
			factlev<-levels(ff)
			dummy<-as.data.frame(model.matrix(~ ff - 1)[,-1])
			colnames(dummy)<-paste(colnames(rdata)[varpost[i]], factlev[-1], sep = "")
			rdata<-rdata[,-varpost[i]]
			rdata<-cbind(rdata,dummy)
			
			}
		    }

		}
 	
	glm.MS <- glm(y ~ ., family = "binomial", weights = w.MS, data = rdata)	

	X <- cbind(1, rdata[,-1])
	pi.hat <- glm.MS$fitted.values
	Ihat <- (t(as.matrix(X)) %*% (as.matrix(X) * w.MS * pi.hat * (1 - pi.hat)))/N1
	wgt<-n1/n2*(n1-n2)
	S <- X * (y - pi.hat)
	
	invI <- solve(Ihat)
	result <- 0
	varsi<-ar<-array(0,dim=c(nrow(Ihat),nrow(Ihat),length(wgt)))
	for(i in 1:length(wgt)) {
		si <- S[y == ylevel[i] & z == zlevel[i],  ]
		varsi[,,i]<-var(si)
		ar[,,i]<-invI%*%varsi[,,i]%*%invI
		result <- result + var(si) * wgt[i]
		}
	
	Vhat <- result/N1	
	
	V <- (invI + invI %*% Vhat %*% invI)/N1
	
	z.value <- glm.MS$coef/sqrt(diag(V))
	p.value <- 2*(1-pnorm(abs(z.value)))

	if (print.all)
	list(table=cbind(ylevel=ylevel, zlevel=zlevel, n1=n1, n2=n2),
        parameters=cbind(est=glm.MS$coef,se=sqrt(diag(V)),z=z.value,pvalue=p.value),
			Wzy=ar,Ihat=Ihat,varsi=varsi)
	else list(table=cbind(ylevel=ylevel, zlevel=zlevel,n1=n1,n2=n2),
	  parameters=cbind(est=glm.MS$coef,se=sqrt(diag(V)),z=z.value,pvalue=p.value))
  }

# close function call
}




 # @@@@@@@@@@@@@@@@@@@@@@@     CODING FUNCTION @@@@@@@@@@@@@@@@@@@@@@@@@@@@

 coding<-function(x=x,y=y,z=z,return=FALSE) {
 ### This function is used to combine multiple columns of z into one column
 ### If used with combined first and second stage data (i.e. with NA for missing
 ### values), it will return sample sizes for the first and second stage 
 ### for each (Y,Z) stratum. If used with only second stage (i.e. complete) data
 ### it will return the second stage sample sizes in each (Y,Z) stratum.
 ### This function should be run on second stage data prior to calling 
 ### optimal sampling functions that require sample sizes or prevalences
 ### as input, as it illustrates the order in which sample sizes should
 ### be provided.
 

 z1<-data.frame(z)
 z.old<-as.matrix(z)

 if (ncol(z1)>1){
 	ncz<-ncol(z1)
  	nrz<-nrow(z1)
  	zlst<-leve<-NULL
  	for (i in 1:ncz){
      	zlst<-c(zlst,list(z1[,i]))
      	leve<-c(leve,length(levels(as.factor(z1[,i]))))
      	zlst[[i]]<-as.factor(zlst[[i]])
      	levels(zlst[[i]])<-c(1:leve[i])
      }
      z1<-matrix(unlist(zlst),nrz,ncz)
        
	m<-max(leve)
	m1<-m^c(1:ncz)
	nz<-z1%*%m1	
	nz<-as.factor(nz)
	nlev<-length(as.numeric(levels(nz)))
	levels(nz)<-c(1:nlev)
	list(nz=as.numeric(nz),z=z1)
	z<-as.numeric(nz)
 }

 levels(z)<-1:length(levels(as.factor(z))) 
        
 ########  now prepare levels of new Z for printing        

 id<-1:length(z)
 index<-NULL
 nlev<-length(levels(as.factor(z))) 
       
 for (i in 1:nlev){ 
 	if (ncol(z1)>1){
      	id1<-id[z==levels(z)[i]]
            id1<-sample(id1,1)
            index<-c(index,id1)
	}
            
      else {
            id1<-id[z==levels(as.factor(z))[i]]
            id1<-sample(id1,1)
            index<-c(index,id1)}
      }

 data<-data.frame(y,z,x)
 n1<-c(t(table(y,z)))
 Cdata<-na.omit(data)
 n2<-c(t(table(Cdata[,1],Cdata[,2])))

 ylev<-as.numeric(levels(factor(y)))
 zlev<-as.numeric(levels(factor(z)))
 ylevel<-rep(ylev,rep(length(zlev),length(ylev)))
 zlevel<-rep(zlev,length(ylev))
 index<-rep(index,length(ylev))         

 ### now label the columns of z for printing

 if (is.null(colnames(z.old)))
 	colnames(z.old)<-paste("z",1:ncol(z.old),sep="") 

 if (ncol(z1)>1){
 	if (sum(n1==n2)<length(n1)){ 
      	print("For calls requiring n1 or prev as input, use the following order!!")
        	print(cbind(ylevel=ylevel,z.old[index,],new.z=zlevel,n1=n1,n2=n2))
	}

      else {
		print("For calls requiring n1 or prev as input, use the following order")
           	print(cbind(ylevel=ylevel,z.old[index,],new.z=zlevel,n2=n2))
	}
 }

 else{ 
     if (sum(n1==n2)<length(n1)){
 		print("For calls requiring n1 or prev as input, use the following order")
           	print(cbind(ylevel=ylevel,z=z.old[index],new.z=zlevel,n1=n1,n2=n2))}
     else {
   		print("For calls requiring n1 or prev as input, use the following order")
         	print(cbind(ylevel=ylevel,z=z.old[index],new.z=zlevel,n2=n2))}
     }

 ######## end of code for printing original and recoded Z

 #### if return=TRUE, return the original z variable(s)
 if (return)
 return(z=z,z.old=z.old)
 }


