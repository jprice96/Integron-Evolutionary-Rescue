evol.rescue = function(k, start.pop, repro, death, int.activity, reinsertion, t.max){
  #Values that need to be provided to the function:
  #k, single value giving the number of cassettes in the starting integron
  #start.pop, either a vector giving individual numbers for each potential integron genotype 
  #OR a single value which will be taken as the number of 0...1 individuals and assume all others are zero
  #repro, a vector of length k or k+1 giving the reproductive rates of each integron genotype
  #death, a vector of length k or k+1 giving the death rates of each integron genotype
  #int.activity, a single value giving the rate of integrase activity in non-adaptive integron genotypes
  #reinsertion, probability that a cassette once excised will be reinserted by the integrase
  #t.max, the maximum number of generations to simulate for
  
  ## Create the blank table for the data to be entered in to
  genotypes = NA
  for (i in 1:k){
    array.string = rep(0, i)
    for(j in 1:(i+1)){
      index = ((i+1)*i / 2) - 1
      genotypes[index+j] = paste0(array.string, collapse = "")
      array.string = rep(0,i)
      array.string[((i+1)-j)] = 1
    }
  }
  genotypes = c("empty", genotypes)
  
  pop = matrix(NA, ncol = length(genotypes) + 2, nrow=1)
  pop = as.data.frame(pop)
  names(pop) = c("Time", genotypes, "Total Population")
  
  ## Set up starting population values
  table.size = length(pop[1,])
  pop[1,] = 0
  #Input provided starting population values
  if(length(start.pop) == 1){
    pop[1,(table.size-k)] = start.pop
  }
  if(length(start.pop) == (table.size-2)){
    pop[1,2:(table.size-k)] = start.pop
  }
  pop[1,table.size] = sum(pop[1,2:(table.size-1)])
  if(pop[1,'Total Population'] == 0){
    print("start.pop incorrect number of values")
  }
  
  ## Generate probability modifiers for integrase events
  table.size = ((k+2)*(k+1)/2)
  prob.table = matrix(NA, ncol=table.size, nrow=table.size, dimnames = list(genotypes, genotypes))
  
  start.col = 1
  
  for(i in 1:(k+1)){
    #Find the starting coordinates of the block, block dimensions = i+1 x i+1
    start.row = start.col
    #row = c(1, 1, 2, 4, 7)
    start.col = start.col + (i-1)
    #col = c(1, 2, 4, 7, 11)
    
    # (n-1)/k
    for(j in 1:(i-1)){
      row = start.row + (j-1)
      col = start.col + j
      prob.table[row, col] = (j-1)/(i-1)
    }
    
    # 1/k
    for(j in 1:i){
      row = start.row
      col = start.col + (j-1)
      prob.table[row, col] = 1/(i-1)
    }
    
    # (k-n)/k
    for(j in 1:i){
      row = start.row + (j-1)
      col = start.col + (j-1)
      prob.table[row, col] = ((i-1) - (j-1))/(i-1)
    }
  }
  
  start.row = 1
  
  for(i in 1:(k+1)){
    #Find the starting coordinates of the block, block dimensions = i+1 x i+1
    start.row = start.row + (i-1)
    #row = c(1, 2, 4, 7, 11)
    start.col = start.row
    #col = c(1, 2, 4, 7, 11)
    
    # -(1-n)/k
    if(i>2){
      for(j in 1:(i-2)){
        row = start.row + (j)
        col = start.col + (j+1)
        prob.table[row, col] = (-1)*(1-(j+1))/(i-1)
      }
    }
    
    # (n-(k+1))/(k+1)
    for(j in 1:i){
      row = start.row + (j-1)
      col = start.col + (j-1)
      prob.table[row, col] = (i-j)/(i-1)
    }
    
    # 1/k
    if(i>1){
      for(j in 1:(i-1)){
        row = start.row + (i-1)
        col = start.col + j
        prob.table[row, col] = 1/(i-1)
      }
    }
  }
  
  #Identify which genotypes are valid targets for integrase activity
  #We are assuming that genotypes with '1' in the primary position will not experience integrase activity
  valid.int = 1:length(genotypes)
  primary.pos.genotype = 1
  
  for(i in 1:(k+1)){
    primary.pos.genotype[i] = i*(i+1) / 2
    valid.int[primary.pos.genotype] = NA
  }
  
  valid.int = valid.int[!is.na(valid.int)]
  
  t=1
  tau = NA
  #tau.storage = NA
  
  while(pop[t,"Time"]<t.max && pop[t,"Total Population"]<=start.pop && pop[t, "Total Population"]>0){
    #Data frame that will hold the number of each event that will happen this time interval
    rates = data.frame(Geno=NA, Birth=NA, Dies=NA, Loss=NA, Reinsert=NA)
    
    #Determine birth rate
    rates[1,"Geno"] = genotypes[1]
    rates[1,"Birth"] = pop[t, genotypes[1]] * repro[length(repro)]
    
    for(i in 2:length(genotypes)){
      #Turn genotype into numeric vector for evaluation
      cassette.array = as.numeric(unlist(strsplit(genotypes[i], split="")))
      for(j in 1:length(cassette.array)){
        #Identify which position the '1' is in and assign the appropriate reproductive rate to the current number of individuals
        if(cassette.array[j]==1){
          rates[i,"Geno"] = genotypes[i]
          rates[i,"Birth"] = pop[t, genotypes[i]] * repro[j]
        }
      }
      #If there is no '1' present, and so the value present is still NA, use the kth value in the repro vector (assumed equivalent)
      if(is.na(rates[i,"Birth"])){
        rates[i,"Geno"] = genotypes[i]
        rates[i,"Birth"] = pop[t, genotypes[i]] * repro[k]
      }
    }
    
    #Determine death rate
    rates[1,"Dies"] = pop[t, genotypes[1]] * death[length(death)]
    
    for(i in 2:length(genotypes)){
      #Turn genotype into numeric vector for evaluation
      cassette.array = as.numeric(unlist(strsplit(genotypes[i], split="")))
      for(j in 1:length(cassette.array)){
        #Identify which position the '1' is in and assign the appropriate death rate to the current number of individuals
        if(cassette.array[j]==1){
          rates[i,"Geno"] = genotypes[i]
          rates[i,"Dies"] = pop[t, genotypes[i]] * death[j]
        }
      }
      #If there is no '1' present, and so the value present is still NA, use the kth value in the death vector (assumed equivalent)
      if(is.na(rates[i,"Dies"])){
        rates[i,"Geno"] = genotypes[i]
        rates[i,"Dies"] = pop[t, genotypes[i]] * death[k]
      }
    }
    
    #Determine rate of integrase-related events
    #Firstly, excision of a random cassette followed by its loss due to lack of reintegration
    rates[,"Loss"] = 0
    
    for(i in 1:length(valid.int)){
      rates[valid.int[i], "Loss"] = pop[t, genotypes[valid.int[i]]] * int.activity * (1 - reinsertion)
    }
    
    #Secondly, excision of a random cassette followed by it being successfully reintegrated in the primary position
    rates[,"Reinsert"] = 0
    
    for(i in 1:length(valid.int)){
      rates[valid.int[i], "Reinsert"] = pop[t, genotypes[valid.int[i]]] * int.activity * reinsertion
    }
    
    #Determine 'tau', the time interval based off the rates
    alpha = unlist(as.list(as.matrix(rates[,2:5])))
    alpha.sum = sum(alpha)
    alpha.mean = mean(alpha)
    alpha.var = var(alpha)
    
    tau[1] = (0.01 * alpha.sum)/alpha.mean
    tau[2] = ((0.01^2) * (alpha.sum^2))/alpha.var
    tau = min(tau)
    
    if(tau<(1/alpha.sum)){
      tau = 1/alpha.sum 
    }
    
    #tau.storage[t] = tau
    
    #Now, based off these rates, determine how many of these events will happen using a Poisson distribution
    events = data.frame(Geno=NA, Birth=NA, Dies=NA, Loss=NA, Reinsert=NA)
    for(i in 1:length(genotypes)){
      events[i,"Geno"] = rates[i,"Geno"]
      events[i,"Birth"] = rpois(1, rates[i,"Birth"]*tau)
      events[i,"Dies"] = rpois(1, rates[i,"Dies"]*tau)
      events[i,"Loss"] = rpois(1, rates[i,"Loss"]*tau)
      events[i,"Reinsert"] = rpois(1, rates[i,"Reinsert"]*tau)
    }
    
    #Calculate the result of these events for the next time interval
    #Set the next time point
    pop[t+1, "Time"] = pop[t, "Time"] + tau
    pop[t+1, 2:(length(pop[1,])-1)] = pop[t, 2:(length(pop[1,])-1)]
    
    #First resolve birth events
    for(i in 1:length(genotypes)){
      pop[t+1, genotypes[i]] = pop[t+1, genotypes[i]] + events[i, "Birth"]
    }
    
    #Next resolve death events
    for(i in 1:length(genotypes)){
      pop[t+1, genotypes[i]] = pop[t+1, genotypes[i]] - events[i, "Dies"]
    }
    
    #Next resolve excision with loss events
    for(i in 1:length(genotypes)){
      #If the genotype is a valid genotype for a shuffle event
      if(i %in% valid.int && events[i,"Loss"]>0){
        #Each event potentially results in a different result, so each excision event is run individually
        for(j in 1:events[i,"Loss"]){
          random = runif(n=1, min=0, max=1)
          outcomes = prob.table[,genotypes[i]]
          for(h in 1:length(outcomes)){
            if(random <= sum(outcomes[1:h], na.rm=TRUE) && sum(outcomes[1:h], na.rm=TRUE) <= 1){
              pop[t+1, names(outcomes[h])] = pop[t+1, names(outcomes[h])] + 1
              pop[t+1, genotypes[i]] = pop[t+1, genotypes[i]] - 1
              break
            }
          }
        }
      }
    }
    
    #Finally resolve excision with reintegration events
    for(i in 1:length(genotypes)){
      #If the genotype is a valid genotype for a shuffle event
      if(i %in% valid.int && events[i,"Reinsert"]>0){
        #Each event potentially results in a different result, so each excision event is run individually
        for(j in 1:events[i,"Reinsert"]){
          random = runif(n=1, min=0, max=1)
          outcomes = prob.table[,genotypes[i]]
          for(h in 1:length(outcomes)){
            if((random+1) <= sum(outcomes[1:h], na.rm=TRUE) && sum(outcomes[1:h], na.rm=TRUE) > 1){
              pop[t+1, names(outcomes[h])] = pop[t+1, names(outcomes[h])] + 1
              pop[t+1, genotypes[i]] = pop[t+1, genotypes[i]] - 1
              break
            }
          }
        }
      }
    }
    
    #Correct for cases where more individuals are being removed than existed (e.g. 2 deaths but only 1 individual)
    for(i in 1:length(genotypes)){
      if(pop[t+1, genotypes[i]] < 0){
        pop[t+1, genotypes[i]] = 0
      }
    }
    
    #Calculate the total population size following these events
    pop[t+1, length(pop[1,])] = sum(pop[t+1, 2:(length(pop[1,])-1)])
    
    #Advance to the next time step
    t = t + 1
  }
  return(pop)
}