#From Wu M, Diez-Roux A, Raghunathan TE, Sánchez BN. 
#FPCA-based method to select optimal sampling schedules that capture between-subject variability in longitudinal studies. 
#Biometrics. 2018;74(1):229-238. doi:10.1111/biom.12714
library(orthogonalsplinebasis)
library(MASS)


get_Btheta=function(i, B, theta)
{
	B[[i]] %*% theta
}

get_alphaVar=function(i, s2, Dinv, Btheta)
{
	solve((Dinv + t(Btheta[[i]]) %*% Btheta[[i]] / s2))
}

get_alpha=function(i,  alphaVar, s2, Btheta, y_mu)
{
	alphaVar[,,i] %*% (t(Btheta[[i]]) %*% y_mu[[i]]) / s2
}

get_aa=function(i, alpha, alphaVar)
{
	outer(alpha[,i], alpha[,i])+alphaVar[,,i]
}

get_aaBB=function(n_basis, m, aa, BB, penalty_matrix)
{
	aaBB=array(NA, dim=c(n_basis, n_basis, m))
	for (j in 1:m)
	{
		a2=aa[j,j,]
		aaBB[,,j]=apply(rep(a2, each=n_basis^2) * BB, c(1,2), sum) + penalty_matrix
	}
	aaBB
}

get_partial_random_effect=function(i, j, y_mu, B, theta, alpha, aa)
{
	aa_leave_one_out=c(aa[j,-j,i])
	t(B[[i]]) %*% (alpha[j,i]*y_mu[[i]]-B[[i]] %*% theta[,-j] %*% aa_leave_one_out)
}

get_partial_random_effect_for_1pc=function(i, y_mu, B, alpha)
{
	t(B[[i]]) %*% y_mu[[i]]  * alpha[1,i]
}
get_theta=function(j, y_mu, aaBB, B, theta, alpha, aa, null_space)
{
	n=length(y_mu)
	null_space %*%
		solve(t(null_space) %*% aaBB[,,j] %*% null_space, 
			t(null_space) %*% 
				apply(sapply(1:n, get_partial_random_effect, j, y_mu, B, theta, alpha, aa),
					1, sum)
		)
}

get_s2=function(i, y_mu, alpha, alphaVar, Btheta)
{
	c(sum((y_mu[[i]]-Btheta[[i]] %*% alpha[,i])^2)+
		tr(Btheta[[i]] %*% alphaVar[,,i] %*% t(Btheta[[i]]))
		)
}

get_y_random=function(i, y, B, theta, alpha)
{
	y[[i]]-B[[i]] %*% theta %*% alpha[,i]
}

normalize=function(x)
{
	x/sqrt(sum(x^2))
}

tr=function(x)
{
	sum(diag(x))
}


####################################################################################################
fpca_em=function(formula, data, grouping, m, knots, lambda, 
		mean_parm_start, theta_start, D_start, s2_start, tol, status=F, showlik=F,
		response_name=NA, time_name=NA, maxiter=5000)
# m is the number of components
{
# figure the names
	if(length(formula)==3) 
	# if formula==NA, then length(formula)==1
	# Yes. This is an ugly hack. Hope to find a better solution.
	{
		response_name=formula[[2]]
		mean_parm_name=names(mean_parm_start)
		time_name=setdiff(all.vars(formula[[3]]), c(response_name, mean_parm_name))
	}
	
# Translate formula
	Y=data[[as.character(response_name)]]
	TIME=data[[as.character(time_name)]]
	SID=data[[as.character(grouping)]]


# Restructure the data
	y=tapply(Y, SID, list)
	time=tapply(TIME, SID, list)
	n=length(unique(SID))
	n_obs=length(Y)


# Pre calculate the bases
	obase=OrthogonalSplineBasis(expand.knots(knots))
	B=list()
	n_basis=dim(obase)[2]
	BB=array(NA, dim=c(n_basis, n_basis, n))
	for (i in 1:n)
	{
		B[[i]]=evaluate(obase, time[[i]])
		BB[,,i]=t(B[[i]]) %*% B[[i]]
	}


# Initialization
	id_idx=1:n
	mean_parm=mean_parm_start
	y_random=y # to begin with, random effects are zero.
	s2=s2_start
	if (any(is.na(theta_start)))
		theta=matrix(rnorm(m*n_basis, mean=0, sd=0.1), nrow=n_basis, ncol=m)
	else
		theta=theta_start
	if (any(is.na(D_start)))
		D=diag(1, m)
	else
		D=D_start
	Dinv=solve(D)

	penalty_matrix_original=OuterProdSecondDerivative(obase)
	
	parm_diff=2*tol
	iter_counter=0
	while (parm_diff > tol && iter_counter<=maxiter)
	{
		iter_counter=iter_counter+1

		last_mean_parm=mean_parm
		last_theta=theta
		last_D=D
		last_s2=s2
		
# M-step for mean parameters
		if(length(formula)==3) # again, the same ugly hack.
		{
			working_data=data
			working_data[[response_name]]=unlist(y_random)
			working_fit=nls(formula, data=working_data, start=mean_parm)
			mean_parm=coefficients(working_fit)
			y_mu=tapply(Y-predict(working_fit), SID, list)
		}
		else
			y_mu=tapply(Y, SID, list)

		penalty_matrix=penalty_matrix_original*lambda*2*s2

# E-step
		Btheta=lapply(id_idx, get_Btheta, B, theta)
		alphaVar=sapply(id_idx, get_alphaVar, s2, Dinv, Btheta)
		alphaVar=array(c(alphaVar), dim=c(m, m, n))
		alpha=sapply(id_idx, get_alpha, alphaVar, s2, Btheta, y_mu)
		alpha=array(c(alpha), dim=c(m, n))
		aa=sapply(id_idx, get_aa, alpha, alphaVar)
		aa=array(c(aa), dim=c(m, m, n))

# M-step for principal components: theta
		aaBB=get_aaBB(n_basis, m, aa, BB, penalty_matrix)
		if (m>=2)
		{
			old_theta=theta+2*tol
			#while (max(abs(theta-old_theta))>tol)
			{
				old_theta=theta
				for (j in 1:m)
				{
					null_space=svd(t(theta[,-j]), nv=n_basis)$v[,m:n_basis]

					theta[,j]=get_theta(j, y_mu, aaBB, B, theta, alpha, aa, null_space)
					theta[,j]=theta[,j]/sqrt(sum(theta[,j]^2))

				}
			}
		}
		else
		{
			theta=solve(aaBB[,,1], apply(sapply(1:n, get_partial_random_effect_for_1pc, y_mu, B, alpha), 1, sum))
			theta=array(theta, dim=c(n_basis, 1))
			theta=theta/sqrt(sum(theta^2))
		}

# M-step for variance parameters
		D=apply(aa, c(1,2), mean)
		#diag(D)=diag(apply(aa, c(1,2), mean))
		Dinv=solve(D)

		s2=sum(sapply(id_idx, get_s2, y_mu, alpha, alphaVar, Btheta))/n_obs

# Calculate y_random: y minus the random effects
		y_random=sapply(id_idx, get_y_random, y, B, theta, alpha, simplify=F)

# Perform ratotation on theta. If m=1, there is no need to do rotation at all.
		if (m>=2)
		{
			SVD=svd(D)
			u=SVD$u
			for (j in 1:m)
			{
				if (u[j,j] < 0)
					u[,j]=-u[,j]
			}
			theta=theta %*% u
			D[]=0
			diag(D)=SVD$d
		}

# Calculate the difference between the current and last parameter estimate
		parm_diff=max(abs(c(mean_parm-last_mean_parm, last_theta-theta, last_D-D, last_s2-s2)))	
		#if (parm_diff<0.01)
		#print(c(D, s2))

		if (status==T)
			print(parm_diff)
		if (showlik==T)
			print(loglik(formula, grouping, data, m, knots, lambda, mean_parm, theta, D, s2, response_name, time_name))

	}

	list(formula=formula, knots=knots, m=m, lambda=lambda, mean_parm=mean_parm, theta=theta, D=D, s2=s2, alpha=alpha, 
		formula=formula, time_name=time_name,response_name=response_name)
		
}

plot_fpca=function(time, m, knots, theta, ylim, col_list,lty_list, lwd_list, myylab, myxlab)
{
	plot(range(time), ylim, type='n',ylab=myylab, xlab=myxlab, cex=1.5, cex.lab=1.5, cex.axis=1.5 )
	obase=OrthogonalSplineBasis(expand.knots(knots))
	B=evaluate(obase, time)
	for (j in 1:m)
	{
		lines(time, B %*% theta[,j], col=col_list[j],lty=lty_list[j], lwd=lwd_list[j])
	}
}


predict_mean=function(fit, time)
{
	temp_env=list()
	temp_env[[fit$time_name]]=time
	temp_env=as.data.frame(temp_env)
	temp_env=append(temp_env, as.list(fit$mean_parm))

	obase=OrthogonalSplineBasis(expand.knots(fit$knots))
	B=evaluate(obase, time)
	y_pred=eval(as.expression(fit$formula[[3]]), env=temp_env)
	y_pred
}

plot_fpca_with_mean=function(fit, time, ylim,mylab,myxlab)
{
	temp_env=list()
	temp_env[[fit$time_name]]=time
	temp_env=as.data.frame(temp_env)
	temp_env=append(temp_env, as.list(fit$mean_parm))

	obase=OrthogonalSplineBasis(expand.knots(fit$knots))
	B=evaluate(obase, time)
	y_pred=eval(as.expression(fit$formula[[3]]), env=temp_env)

	for (j in 1:ncol(fit$theta))
	{
		plot(range(time), ylim, type='n', xlab=myxlab , ylab=mylab)
		title(paste('PCA function ', j, ' (Percentage of Variability ', format(fit$D[j,j]/(sum(fit$D)+fit$s2)*100, digits=3), ')'))
		lines(time, y_pred)
		lines(time, y_pred+sqrt(fit$D[j,j]) * B %*% fit$theta[,j], lty=2)
		lines(time, y_pred-sqrt(fit$D[j,j]) * B %*% fit$theta[,j], lty=2)
	}
}

plot_fpca_with_mean2=function(fit, time, mu, ylim, prefix=NA)
{

	obase=OrthogonalSplineBasis(expand.knots(fit$knots))
	B=evaluate(obase, time)
	y_pred=mu

	for (j in 1:ncol(fit$theta))
	{
		if (!is.na(prefix))
			pdf(paste(prefix, j, '.pdf', sep=''))
		plot(range(time), ylim, type='n', xlab=fit$time_name, ylab=fit$response_name, cex=1.5, cex.lab=1.5, cex.axis=1.5)
		title(paste('PCA function ', j, ' (Percentage of Variability ', format(fit$D[j,j]/(sum(fit$D)+24*fit$s2)*100, digits=3), ')'), cex=1.5)
		lines(time, y_pred)
		lines(time, y_pred+sqrt(fit$D[j,j]) * B %*% fit$theta[,j], lty=2)
		lines(time, y_pred-sqrt(fit$D[j,j]) * B %*% fit$theta[,j], lty=2)
		if (!is.na(prefix))
			dev.off()
	}
}


compute_penalty=function(theta, pm)
{
	c(theta %*% pm %*% theta)
}


loglik_per_subject=function(i, time, y_mu, B, theta, D, s2)
{
	PI=3.14159
	ni=length(y_mu[[i]])
	BT=B[[i]] %*% theta
	Sigma=BT %*% D %*% t(BT)
	diag(Sigma)=diag(Sigma) + s2

	p1=-ni/2*log(2*PI)
	p2=-1/2*log(det(Sigma))
	p3=-1/2* c( y_mu[[i]] %*% solve(Sigma, y_mu[[i]]) )

	p1+p2+p3
}

loglik=function(formula, grouping, data, m, knots, lambda, mean_parm, theta, D, s2, response_name=NA, time_name=NA)
{
# figure the names
	if(length(formula)==3) 
	# if formula==NA, then length(formula)==1
	# Yes. This is an ugly hack. Hope to find a better solution.
	{
		response_name=formula[[2]]
		mean_parm_name=names(mean_parm)
		time_name=setdiff(all.vars(formula[[3]]), c(response_name, mean_parm_name))
	}
	
# Translate formula
	Y=data[[as.character(response_name)]]
	TIME=data[[as.character(time_name)]]
	SID=data[[as.character(grouping)]]

# Restructure the data
	y=tapply(Y, SID, list)
	time=tapply(TIME, SID, list)
	n=length(unique(SID))
	n_obs=length(Y)


# Pre calculate the bases
	obase=OrthogonalSplineBasis(expand.knots(knots))
	B=list()
	n_basis=dim(obase)[2]
	BB=array(NA, dim=c(n_basis, n_basis, n))
	for (i in 1:n)
	{
		B[[i]]=evaluate(obase, time[[i]])
		BB[,,i]=t(B[[i]]) %*% B[[i]]
	}
	p_m=OuterProdSecondDerivative(obase)
	penalty=sum(apply(theta, 2, compute_penalty, p_m))


	if(length(formula)==3) 
	# if formula==NA, then length(formula)==1
	# Yes. This is an ugly hack. Hope to find a better solution.
	{
		temp_env=append(data, as.list(mean_parm))
		y_mu=eval(as.expression(formula[[3]]), env=temp_env)
	}
	else
	{
		y_mu=0 #R will auto expand.
	}
	y_mu=tapply(Y-y_mu, SID, list)
	llik=sum(sapply(1:n, loglik_per_subject, time, y_mu, B, theta, D, s2))
	c(loglik=llik, penalty=-penalty, l_penalty=-lambda*penalty, combined=llik-lambda*penalty)
}
	

cv_fpca_em=function(formula, data, 
			grouping, 
			m, knots, lambda, 
			mean_parm_start, theta_start, D_start, s2_start, 
			cv_grouping, 
			tol, status=F, showlik=F,
			response_name=NA, time_name=NA)
{
	n.m=length(m)
	n.lambda=length(lambda)
	#cv_group=sort(unique(data[[cv_grouping]]))
	cv_group=sort(unique(data$cv_grouping))
	n.cv_group=length(cv_group)
	cv_score=matrix(NA, nrow=n.m*n.lambda, ncol=n.cv_group)
	cv=data.frame(m=rep(m, each=n.lambda), lambda=rep(lambda, n.m))
	for (j in 1:n.cv_group)
	{
#		training_data=data[data[[cv_grouping]]!=cv_group[j],]
#		test_data=data[data[[cv_grouping]]==cv_group[j],]
		training_data=data[data$cv_grouping!=cv_group[j],]
		test_data=data[data$cv_grouping==cv_group[j],]
		
		for (i in 1:(n.m*n.lambda))
		{
			training_fit=fpca_em(formula, training_data, 
						grouping, m=cv$m[i], knots, lambda=cv$lambda[i], 
						mean_parm_start, theta_start, D_start, s2_start, 
						tol,  status, showlik,
						response_name=response_name, time_name=time_name)

			cv_score[i,j]=loglik(formula, grouping, 
						test_data, m=cv$m[i], knots, lambda=cv$lambda[i],
						training_fit$mean_parm, training_fit$theta, training_fit$D, training_fit$s2,
						response_name=response_name, time_name=time_name)[1]
			#print(cv_score[,j])
		}
	}
  
	cv$score=apply(cv_score, 1, mean)
	cv
}
 

select_m_lambda=function(cv_score, m, lambda,cutoff= 0.001 )
{
  #browser()
  good_score=tapply(cv_score, m, max)
  percent=-(good_score[2:length(good_score)]-good_score[1:(length(good_score)-1)])/good_score[1:(length(good_score)-1)]
  #m_idx=which.min((percent<0.01)*(1:length(percent)))
  
  if (percent[1]<cutoff)
    best_m=1
  else
    if (percent[2]<cutoff)
      best_m=2
  else
    best_m=3
  
  #best_m=as.numeric(names(good_score)[m_idx])
  best_lambda=lambda[m==best_m][which.max(cv_score[m==best_m])]
  c(best_m, best_lambda)
}


