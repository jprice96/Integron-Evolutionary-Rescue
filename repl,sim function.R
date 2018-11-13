Q1 = function(x){
  y = as.numeric(quantile(x,0.25))
  return(y)
}
Q3 = function(x){
  y = as.numeric(quantile(x,0.75))
  return(y)
}

repl.sim = function(n, output, simulation, sim.para){
  Timeline = (0:(sim.para$Time * 100))/200
  Storage = data.frame(Time = Timeline)
  
  Values = as.data.frame(matrix(ncol=length(output)))
  names(Values) = output
  reset.values = Values
  
  for(i in 1:n){
    #Run the chosen function with the provided parameters, making sure that the parameters in the sim.para object are correctly named
    test = do.call(simulation, sim.para)
    
    #Create a sliding window between Timeline intervals
    x = 2
    for(j in 2:length(Timeline)){
      Values = reset.values
      x = x - 1
      #If a value is produced within the window being assessed, store it in Values
      while(test[x, "Time"]<Timeline[j]){
        if(test[x, "Time"]>=Timeline[j-1]){
          for(m in 1:length(output)){
            Values[,m] = c(Values[,m], test[x, output[m]])
          }
        }
        x = x + 1
        if(is.na(test[x, "Time"])){break}
      }
      
      #Find the appropriate outputs for the window
      for(p in 1:length(output)){
        Storage[(j-1),(i+1)][p] = do.call(output[p], Values[,m])
      }
    }
  }
  return(Storage)
}





