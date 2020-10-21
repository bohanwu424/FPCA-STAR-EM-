library(MASS)
dat = d2; groupings = 'county'; m = m_0; knots = seq(n_dates)/n_dates; lambda = 10000; 
mean_parm_start = mean_parm_0; theta_start = theta_0;efunc = ef; D_start = D_0; s2_start = s2_0; tol = Tol; status=F; showlik=F;
response_name="value"; time_name="dates"; maxiter= 1000

james_em = function(dat, groupings = 'county', m, knots, lambda, 
                 mean_parm_start = 0, theta_start,efunc, D_start, s2_start, tol, status=F, showlik=F,
                 response_name="value", time_name="dates", maxiter= 2000)
  # m is the number of components
{ # Translate formula
  Y=dat[[as.character(response_name)]]
  TIME=dat[[as.character(time_name)]]
  SID=dat[[as.character(groupings)]]
  

  # Restructure the data
  y=tapply(Y, SID, list)
  time=tapply(TIME, SID, list)
  n=length(unique(SID))
  n_obs=length(Y)
  
  
  # Pre calculate the bases
  orders = 4
  obase=OrthogonalSplineBasis(knots  = expand.knots(knots))
  n_basis=dim(obase)[2]
  B=list()
  BB=array(NA, dim=c(n_basis, n_basis, n))
  for (i in 1:n)
  {
    B[[i]]=evaluate(obase, time[[i]])
    length(time[[i]])
    BB[,,i]=t(B[[i]]) %*% B[[i]]
  }
  
  if(any(is.na(theta_start))){
    theta_start = ginv(evaluate(obase, knots))%*% efunc #generalize inverse Theta = S^{-1} %*% F
  }
  # Initialization
  id_idx=1:n
  mean_parm=mean_parm_start
  y_random=y # to begin with, random effects are zero.
  s2=s2_start
  if (any(is.na(theta_start))){
    theta=matrix(rnorm(m*n_basis, mean=0, sd=0.1), nrow=n_basis, ncol=m)
  }else{
    theta=theta_start
  }
  if (any(is.na(D_start))){
    D=diag(1, m)
  }else{
    D=D_start
  }
  Dinv=solve(D)
  
  penalty_matrix_original=OuterProdSecondDerivative(obase)
  
  parm_diff=2*tol
  iter_counter=0
  log_like_old = Inf
  log_like_new = -Inf
  log_like_diff = 2*tol
  
  while (parm_diff > tol && iter_counter<=maxiter)
  {
    iter_counter=iter_counter+1
    
    last_mean_parm=mean_parm
    last_theta=theta
    last_D=D
    last_s2=s2
    formula = "value~"
    # M-step for mean parameters
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
    }else
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
    parm_diff=max(abs(c(mean_parm-last_mean_parm, last_theta-theta, last_D-D, last_s2-s2)),na.rm = TRUE)	
    #if (parm_diff<0.01)
    #print(c(D, s2))
    if (status==T)
      print(parm_diff)
    if (showlik==T)
      print(loglik(formula, groupings, dat, m, knots, lambda, mean_parm, theta, D, s2, response_name, time_name))
    # Calculate the difference between the current and last log-likelihood estimate
    log_like_new = loglik(formula, groupings, dat, m, knots, lambda, mean_parm, theta, D, s2, response_name, time_name)[1]
    log_like_diff = abs(log_like_new - log_like_old)
    log_like_old = log_like_new
    
  }

  output = list(fitted.values = matrix(), coefficients = matrix())
  #fitted values#
  Yy = matrix(nrow = n, ncol = length(knots))
  for(i in seq(n)){
    Yy[i,1:nrow(B[[i]])] = B[[i]]%*%theta %*% alpha[,i]
  }
  output$fitted.values <- as.vector(t(Yy)) %>% na.exclude()
  output$coefficients <- c(length(unique(TIME)), m,
                           as.vector(theta),as.vector(evaluate(obase, knots)%*% theta),
                           diag(D), s2)
  output$iter = iter_counter
    return(output)
}
