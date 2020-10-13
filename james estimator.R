james_estimator = function(v = y_f, theta){
  #theta = (n_dates, m, theta_start, efcn,D_0, S2_0)
  m_0 = theta[2]
  n_dates = theta[1]
  theta_0 = matrix(theta[3:(length(theta) - m - m_0*n_dates - 1)], ncol = m_0 )
  ef = matrix(theta[(length(theta) - m - m_0*n_dates):(length(theta) - m - 1)], nrow = n_dates )
  D_0 = diag(theta[(length(theta) - m): (length(theta) - 1)])
  s2_0 = theta[length(theta)]
  mean_parm_0 = 0; 
  Tol = 10^{-2}; Lam = 1
  
  #Input dataframe
  input <<- list(Ly = matrix(nrow = n_counties, ncol = n_dates),Lt = Lt_f,all.Ly = v)
  for(i in seq(n_counties)){
    input$Ly[i,input$Lt[[i]]] <<- input$all.Ly[input$Lt[[i]]]
    input$all.Ly <<- input$all.Ly[-input$Lt[[i]]]
  }
  d2 <<- melt(input$Ly)  %>% 
    rename(dates = Var2, county = Var1) %>%
    group_by(county )  %>% 
    drop_na()%>% mutate(dates = dates/n_dates)%>% 
    arrange(county)
  d2 <<- d2[,c("county","value","dates")]
  
  
  output <<- james_em(dat = d2, groupings = 'county', m = m_0, knots = seq(n_dates)/n_dates, lambda = Lam, 
                      mean_parm_start = mean_parm_0, theta_start = theta_0,efunc = ef, D_start = D_0, s2_start = s2_0, tol = Tol, status=F, showlik=F,
                      response_name="value", time_name="dates", maxiter= 1000)
  
}

#fit = james_estimator(round(output$fitted.values),output$coefficients)
