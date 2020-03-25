      ## physical and temporal distancing COVID-19
      ## example code
      ## CJA Bradshaw
      ## March 2020
      
      ## remove everything
      rm(list = ls())
      
      # set up spatial matrix
      mat.dim <- 20
      spat.layer <- matrix(data=0, nrow=mat.dim, ncol=mat.dim)
      
      popn <- 1000 # population size
      infect.pr <- 0.02
      
      max.move <- round(mat.dim/2, 0) # maximum movement in one coordinate per time interval
      min.move <- 0
      move.prob <- 0.5 # probability of moving in a given time interval
      move.pr.red <- 0.2 # equates to an 80% reduction in movement probability
      move.pr <- move.prob*move.pr.red
      exp.pr <- 0.5 # probability of exposure (decreases with increasing temporal distancing)
        
      # infection spread
      t.steps <- 50
      iter <- 1000
      itdiv <- iter/10
      tot.inf.out <- matrix(data=NA, nrow=iter, ncol=t.steps+1)
      
      for (s in 1:iter) {
        
        # place population randomly into grid
        row.rn <- sample(1:mat.dim, popn, replace=T)
        col.rn <- sample(1:mat.dim, popn, replace=T)
        
        spat.init <- spat.layer
        
        for (p in 1:popn) {
          spat.init[row.rn[p], col.rn[p]] <- spat.init[row.rn[p], col.rn[p]] + 1   
        }
        spat.init
        
        # start infection processes
        init.infect <- 1 # initial number of infected in the matrix
        
        # randomly allocate infected number of people to start matrix
        st.pop.coords <- which(spat.init > 0,  arr.ind=T) # populated cell coordinates
        ran.inf.coords <- sample(1:dim(st.pop.coords)[1], init.infect, replace=T)
        st.inf.coords <- as.matrix(t(st.pop.coords[ran.inf.coords,]))
        st.inf.coords
        
        infect.mat <- spat.layer
        infect.mat[st.inf.coords[1,1], st.inf.coords[1,2]] <- init.infect
        not.infect.mat <- spat.init - infect.mat
        
        infect.n <- rep(0,t.steps+1) # number infected time series
        infect.n[1] <- init.infect
        
        for (t in 2:(t.steps+1)) {
        
          for (i in 1:mat.dim) {
            for (j in 1:mat.dim) {
              
              ## total pop in cell
              n.infected <- infect.mat[i,j]
              n.not.infected <- not.infect.mat[i,j]
              cell.pop <- n.not.infected + n.infected
              prop.infected <- n.infected/cell.pop
              
              if (cell.pop > 0) {
                ## new infections
                pr.new.infect <- ifelse(((exp.pr*infect.pr)^(1-prop.infected)) > 1, 1, ((exp.pr*infect.pr)^(1-prop.infected))) # probability of a new infection
                n.new.infect <- rbinom(n=1, size=n.not.infected, prob=pr.new.infect)
                
                # update infection & non-infection layers
                infect.mat[i,j] <- infect.mat[i,j] + n.new.infect
                not.infect.mat[i,j] <- ifelse(cell.pop - infect.mat[i,j] < 0, 0, cell.pop - infect.mat[i,j])
                
                # recalculate infected & non-infected
                n.infected <- infect.mat[i,j]
                n.not.infected <- not.infect.mat[i,j]
                
                # movement
                # move infected people first
                if (n.infected > 0) {
                  move.row.rn <- sample(min.move:max.move, n.infected, replace=T) # move rows
                  move.col.rn <- sample(min.move:max.move, n.infected, replace=T) # move columns
                  
                  dir.row.rn <- sample(c(-1,1), n.infected, replace=T) # direction rows
                  dir.col.rn <- sample(c(-1,1), n.infected, replace=T) # direction columns
            
                  if (rbinom(1,1,move.pr) == 1)  {
                    ## both directions ok
                    for (n in 1:n.infected) {
                      if (i + move.row.rn[n]*dir.row.rn[n] < mat.dim & i + move.row.rn[n]*dir.row.rn[n] > 0) {
                        if (j + move.col.rn[n]*dir.col.rn[n] < mat.dim & j + move.col.rn[n]*dir.col.rn[n] > 0) {
                          infect.mat[i + move.row.rn[n]*dir.row.rn[n], j + move.col.rn[n]*dir.col.rn[n]] <- infect.mat[i + move.row.rn[n]*dir.row.rn[n], j + move.col.rn[n]*dir.col.rn[n]] + 1
                          infect.mat[i, j] <- ifelse((infect.mat[i, j] - 1) < 0, 0, infect.mat[i, j] - 1)
                        }
                      }
                    }
        
                    # cannot move max rows; try other direction
                    for (n in 1:n.infected) {
                      if (i + move.row.rn[n]*dir.row.rn[n] > mat.dim & i + move.row.rn[n]*dir.row.rn[n] > 0) {
                        if (j + move.col.rn[n]*dir.col.rn[n] < mat.dim & j + move.col.rn[n]*dir.col.rn[n] > 0) {
                          infect.mat[i + move.row.rn[n]*-1*dir.row.rn[n], j + move.col.rn[n]*dir.col.rn[n]] <- infect.mat[i + move.row.rn[n]*-1*dir.row.rn[n], j + move.col.rn[n]*dir.col.rn[n]] + 1
                          infect.mat[i, j] <- ifelse((infect.mat[i, j] - 1) < 0, 0, infect.mat[i, j] - 1)
                        }
                      }
                    }
                    
                    # cannot move min rows; try other direction
                    for (n in 1:n.infected) {
                      if (i + move.row.rn[n]*dir.row.rn[n] < mat.dim & i + move.row.rn[n]*dir.row.rn[n] < 0) {
                        if (j + move.col.rn[n]*dir.col.rn[n] < mat.dim & j + move.col.rn[n]*dir.col.rn[n] > 0) {
                          infect.mat[i + move.row.rn[n]*-1*dir.row.rn[n], j + move.col.rn[n]*dir.col.rn[n]] <- infect.mat[i + move.row.rn[n]*-1*dir.row.rn[n], j + move.col.rn[n]*dir.col.rn[n]] + 1
                          infect.mat[i, j] <- ifelse((infect.mat[i, j] - 1) < 0, 0, infect.mat[i, j] - 1)
                        }
                      }
                  }
      
                  # cannot move max cols; try other direction
                  for (n in 1:n.infected) {
                    if (i + move.row.rn[n]*dir.row.rn[n] < mat.dim & i + move.row.rn[n]*dir.row.rn[n] > 0) {
                      if (j + move.col.rn[n]*dir.col.rn[n] > mat.dim & j + move.col.rn[n]*dir.col.rn[n] > 0) {
                        infect.mat[i + move.row.rn[n]*dir.row.rn[n], j + move.col.rn[n]*-1*dir.col.rn[n]] <- infect.mat[i + move.row.rn[n]*dir.row.rn[n], j + move.col.rn[n]*-1*dir.col.rn[n]] + 1
                        infect.mat[i, j] <- ifelse((infect.mat[i, j] - 1) < 0, 0, infect.mat[i, j] - 1)
                      }
                    }
                  }
      
                  # cannot move min cols; try other direction
                  for (n in 1:n.infected) {
                    if (i + move.row.rn[n]*dir.row.rn[n] < mat.dim & i + move.row.rn[n]*dir.row.rn[n] > 0) {
                      if (j + move.col.rn[n]*dir.col.rn[n] < mat.dim & j + move.col.rn[n]*dir.col.rn[n] < 0) {
                        infect.mat[i + move.row.rn[n]*dir.row.rn[n], j + move.col.rn[n]*-1*dir.col.rn[n]] <- infect.mat[i + move.row.rn[n]*dir.row.rn[n], j + move.col.rn[n]*-1*dir.col.rn[n]] + 1
                        infect.mat[i, j] <- ifelse((infect.mat[i, j] - 1) < 0, 0, infect.mat[i, j] - 1)
                      }
                    }
                  }
                  # recalculate infected & non-infected
                  n.infected <- infect.mat[i,j]
                  n.not.infected <- not.infect.mat[i,j]
                  
                }
              }
              
              # move non-infected people now
              if (n.not.infected > 0) {
                move.row.rn <- sample(min.move:max.move, n.not.infected, replace=T) # move rows
                move.col.rn <- sample(min.move:max.move, n.not.infected, replace=T) # move columns
                
                dir.row.rn <- sample(c(-1,1), n.not.infected, replace=T) # direction rows
                dir.col.rn <- sample(c(-1,1), n.not.infected, replace=T) # direction columns
                
                  ## both directions ok
                  for (n in 1:n.not.infected) {
                      if (rbinom(1,1,move.pr) == 1)  {
                      if (i + move.row.rn[n]*dir.row.rn[n] < mat.dim & i + move.row.rn[n]*dir.row.rn[n] > 0) {
                        if (j + move.col.rn[n]*dir.col.rn[n] < mat.dim & j + move.col.rn[n]*dir.col.rn[n] > 0) {
                          not.infect.mat[i + move.row.rn[n]*dir.row.rn[n], j + move.col.rn[n]*dir.col.rn[n]] <- not.infect.mat[i + move.row.rn[n]*dir.row.rn[n], j + move.col.rn[n]*dir.col.rn[n]] + 1
                          not.infect.mat[i, j] <- ifelse((not.infect.mat[i, j] - 1) < 0, 0, not.infect.mat[i, j] - 1)
                        }
                      }
                    }
                  }
      
                  # cannot move max rows; try other direction
                  for (n in 1:n.not.infected) {
                    if (rbinom(1,1,move.pr) == 1)  {
                      if (i + move.row.rn[n]*dir.row.rn[n] > mat.dim & i + move.row.rn[n]*dir.row.rn[n] > 0) {
                        if (j + move.col.rn[n]*dir.col.rn[n] < mat.dim & j + move.col.rn[n]*dir.col.rn[n] > 0) {
                          not.infect.mat[i + move.row.rn[n]*-1*dir.row.rn[n], j + move.col.rn[n]*dir.col.rn[n]] <- not.infect.mat[i + move.row.rn[n]*-1*dir.row.rn[n], j + move.col.rn[n]*dir.col.rn[n]] + 1
                          not.infect.mat[i, j] <- ifelse((not.infect.mat[i, j] - 1) < 0, 0, not.infect.mat[i, j] - 1)
                        }
                      }
                    }
                  }
       
                  # cannot move min rows; try other direction
                  for (n in 1:n.not.infected) {
                    if (rbinom(1,1,move.pr) == 1)  {
                      if (i + move.row.rn[n]*dir.row.rn[n] < mat.dim & i + move.row.rn[n]*dir.row.rn[n] < 0) {
                        if (j + move.col.rn[n]*dir.col.rn[n] < mat.dim & j + move.col.rn[n]*dir.col.rn[n] > 0) {
                          not.infect.mat[i + move.row.rn[n]*-1*dir.row.rn[n], j + move.col.rn[n]*dir.col.rn[n]] <- not.infect.mat[i + move.row.rn[n]*-1*dir.row.rn[n], j + move.col.rn[n]*dir.col.rn[n]] + 1
                          not.infect.mat[i, j] <- ifelse((not.infect.mat[i, j] - 1) < 0, 0, not.infect.mat[i, j] - 1)
                        }
                      }
                    }
                  }
                  
                  # cannot move max cols; try other direction
                  for (n in 1:n.not.infected) {
                    if (rbinom(1,1,move.pr) == 1)  {
                      if (i + move.row.rn[n]*dir.row.rn[n] < mat.dim & i + move.row.rn[n]*dir.row.rn[n] > 0) {
                        if (j + move.col.rn[n]*dir.col.rn[n] > mat.dim & j + move.col.rn[n]*dir.col.rn[n] > 0) {
                          not.infect.mat[i + move.row.rn[n]*dir.row.rn[n], j + move.col.rn[n]*-1*dir.col.rn[n]] <- not.infect.mat[i + move.row.rn[n]*dir.row.rn[n], j + move.col.rn[n]*-1*dir.col.rn[n]] + 1
                          not.infect.mat[i, j] <- ifelse((not.infect.mat[i, j] - 1) < 0, 0, not.infect.mat[i, j] - 1)
                        }
                      }
                    }
                  }
                  
                  # cannot move min cols; try other direction
                  for (n in 1:n.not.infected) {
                    if (rbinom(1,1,move.pr) == 1)  {
                      if (((i + move.row.rn[n]*dir.row.rn[n]) < mat.dim) & ((i + move.row.rn[n]*dir.row.rn[n]) > 0)) {
                        if (j + move.col.rn[n]*dir.col.rn[n] < mat.dim & j + move.col.rn[n]*dir.col.rn[n] < 0) {
                          not.infect.mat[i + move.row.rn[n]*dir.row.rn[n], j + move.col.rn[n]*-1*dir.col.rn[n]] <- not.infect.mat[i + move.row.rn[n]*dir.row.rn[n], j + move.col.rn[n]*-1*dir.col.rn[n]] + 1
                          not.infect.mat[i, j] <- ifelse((not.infect.mat[i, j] - 1) < 0, 0, not.infect.mat[i, j] - 1)
                        }
                      }
                    }
                  }
                }
              } 
            } # end j loop
          } # end i loop
          
          # sum total infected
          infect.n[t] <- sum(infect.mat)
          #print(t)
        
        } # end t loop
      
        #plot(1:(t.steps+1), infect.n, type="l", xlab="time", ylab="number infected")
        if (s %% itdiv==0) print(s) 
        tot.inf.out[s, ] <- infect.n
        
      } # end s loop
      
      inf.ts.med <- apply(tot.inf.out, MARGIN=2, median, na.rm=T)
      inf.ts.lo <- apply(tot.inf.out, MARGIN=2, quantile, probs=0.025, na.rm=T)
      inf.ts.up <- apply(tot.inf.out, MARGIN=2, quantile, probs=0.975, na.rm=T)
      
      plot(1:(t.steps+1), inf.ts.med, type="l", xlab="time", ylab="number infected")
      lines(1:(t.steps+1), inf.ts.lo, lty=2, col="red")
      lines(1:(t.steps+1), inf.ts.up, lty=2, col="red")
      
      out.dat <- data.frame(1:(t.steps+1), inf.ts.med, inf.ts.up, inf.ts.lo)
      colnames(out.dat) <- c("step", "med", "up", "lo")
      setwd("~/Documents/Papers/Disease/COVID-19")
      write.table(out.dat,file="80pcMoveRedHalfExp.out.out.csv",sep=",", row.names = F, col.names = T)
