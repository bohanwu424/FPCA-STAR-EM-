estimator.combined = function(v = y_f){
  input <<- list(Ly = matrix(nrow = n_counties, ncol = n_dates),Lt = Lt_f,all.Ly = v)
  
  for(i in seq(n_counties)){
    input$Ly[i,input$Lt[[i]]] <<- input$all.Ly[input$Lt[[i]]]
    input$all.Ly <<- input$all.Ly[-input$Lt[[i]]]
  }
  #get non_na values and index for each county
  #Use first 10 FPCs
  FPCA.fit <<- fpca.face(Y = input$Ly,npc = n_fpc)
  #input$all.Ly[fit_em$error_inds]
  
  #Fitted Values and Coefficients
  fitted_f.obs = FPCA.fit$Yhat
  coef_f = FPCA.fit$scores
  
  #Input parameters
  n_dates = ncol(FPCA.fit$Yhat);
  mean_parm_0 = 0; D_0 = diag(1/FPCA.fit$evalues);s2_0 = mean(abs(FPCA.fit$Yhat- FPCA.fit$Y))
  Tol = 10^{-2}; Lam = 10
  
  #Input dataframe
  d2 <<- melt(input$Ly)  %>% 
    rename(dates = Var2, county = Var1) %>%
    group_by(county )  %>% 
    drop_na()%>% mutate(dates = dates/n_dates)%>% 
    arrange(county)
  d2 <<- d2[,c("county","value","dates")]
  
 output <<- james_em(dat = d2, groupings = 'county', m = FPCA.fit$npc, knots = seq(n_dates)/n_dates, lambda = Lam, 
          mean_parm_start = mean_parm_0, theta_0 = FPCA.fit$efunctions, D_start = D_0, s2_start = s2_0, tol = Tol, status=F, showlik=T,
          response_name="value", time_name="dates", maxiter= 1000)
  return(output)
}
  