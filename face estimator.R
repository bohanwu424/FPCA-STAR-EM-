estimator.fpca = function(v = y_f)
{
  input <<- list(Ly = matrix(nrow = n_counties, ncol = n_dates),Lt = Lt_f,all.Ly = v)
  
  for(i in seq(n_counties)){
    input$Ly[i,input$Lt[[i]]] <<- input$all.Ly[input$Lt[[i]]]
    input$all.Ly <<- input$all.Ly[-input$Lt[[i]]]
  }
  #get non_na values and index for each county
  #Use first 10 FPCs
  FPCA.fit <<- fpca.face(Y = input$Ly,npc = n_fpc)

  #Fitted Values and Coefficients
  fitted_f.obs = FPCA.fit$Yhat
  theta_0 = array(dim = (ncol(FPCA.fit$Yhat) + 2) *FPCA.fit$npc)
  coef_f = c(ncol(FPCA.fit$Yhat),FPCA.fit$npc, theta_0, 
             as.vector(FPCA.fit$efunctions),
             as.vector(1/FPCA.fit$evalues), 
             mean(abs(FPCA.fit$Yhat- FPCA.fit$Y)) )
  
  #n_dates,m,theta_0,  D_0,s2_0
  #Output
  output <<- list(fitted.values = matrix(), coefficients = matrix())
  
  for(i in seq(n_counties)){
    output$fitted.values <<- append(output$fitted.values,fitted_f.obs[i, input$Lt[[i]]])
  }
  output$fitted.values <<- output$fitted.values %>% na.exclude()
  output$coefficients <<- coef_f 
  return(output)
}

