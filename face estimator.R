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
  #input$all.Ly[fit_em$error_inds]
  
  #Fitted Values and Coefficients
  fitted_f.obs = FPCA.fit$Yhat
  coef_f = FPCA.fit$scores
  #Output
  output <<- list(fitted.values = matrix(), coefficients = matrix())
  
  for(i in seq(n_counties)){
    output$fitted.values <- append(output$fitted.values,fitted_f.obs[i, input$Lt[[i]]])
    output$coefficients <- append(output$coefficients,coef_f[i,])
  }
  output$coefficients <<- append(output$coefficients,as.vector(FPCA.fit$efunctions))
  output$fitted.values <- output$fitted.values %>% na.exclude()
  output$coefficients <- output$coefficients %>% na.exclude()
  return(output)
}
