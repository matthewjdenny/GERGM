#Loading and pre-processing the migration data
load.data = function(directory){
  #Load the needed libraries
  require(ergm, quietly = TRUE)
  require(statnet.common, quietly = TRUE)
  require(network, quietly = TRUE)
  require(combinat, quietly = TRUE)
  require(fossil, quietly = TRUE)
  require(Rcpp, quietly = TRUE)
  require(RcppArmadillo, quietly = TRUE)
  require(gtools, quietly = TRUE)

  #I: Read in Migration Data
  migr05 <- t(read.csv(paste(directory,"/State_to_State2005.csv", sep = ""),stringsAsFactors=F, row.names=1))
  migr06 <- t(read.csv(paste(directory,"/State_to_State2006.csv", sep = ""), stringsAsFactors=F, row.names=1))
  migr07 <- t(read.csv(paste(directory,"/State_to_State2007.csv", sep = ""),, stringsAsFactors=F, row.names=1))
  migr08 <- t(read.csv(paste(directory,"/State_to_State2008.csv", sep = ""), stringsAsFactors=F, row.names=1))

  #II: Get Unemployment, population, latitutde/longitude, temperature and income data

  pop <- read.csv(paste(directory,"/population.csv", sep = ""),stringsAsFactors = F)
  unemp <- read.csv(paste(directory,"/unemployment.csv", sep = ""),stringsAsFactors = F)
  latlong <- read.csv(paste(directory,"/latlong.csv", sep = ""),stringsAsFactors = F)
  temps <- read.csv(paste(directory,"/jantemp.csv", sep = ""),stringsAsFactors=F)
  income <- read.csv(paste(directory,"/income.csv", sep = ""),stringsAsFactors=F)

  ## Make single array of Migration data

  indpr <- which(toupper(colnames(migr06))=="PUERTO_RICO")
  migr05 <- migr05[-indpr,]
  migr05 <- migr05[,-indpr]
  migr06 <- migr06[-indpr,]
  migr06 <- migr06[,-indpr]
  migr07 <- migr07[-indpr,]
  migr07 <- migr07[,-indpr]
  migr08 <- migr08[-indpr,]
  migr08 <- migr08[,-indpr]

  diag(migr05) <- 0
  diag(migr06) <- 0
  diag(migr07) <- 0
  diag(migr08) <- 0

  # Create a single Array #
  migr <- array(0,dim=c(nrow(migr05),nrow(migr05),4))
  migr[,,1] <- migr05
  migr[,,2] <- migr06
  migr[,,3] <- migr07
  migr[,,4] <- migr08

  #########################################################
  ## Make Population Sender and Reciever ##
  popr06 <- matrix(0,nrow(migr06),nrow(migr06))
  pop06 <- pop[,c(5,13)]
  for(i in 1:nrow(migr06)){
    for(j in 1:nrow(migr06)){
      if(i != j){
        row <- which(toupper(pop06[,1]) == toupper(colnames(migr06))[j])
        popr06[i,j] <- pop06[row,2]
      }
    }
  }

  popr07 <- matrix(0,nrow(migr06),nrow(migr06))
  pop07 <- pop[,c(5,14)]
  for(i in 1:nrow(migr06)){
    for(j in 1:nrow(migr06)){
      if(i != j){
        row <- which(toupper(pop07[,1]) == toupper(colnames(migr06))[j])
        popr07[i,j] <- pop07[row,2]
      }
    }
  }

  popr08 <- matrix(0,nrow(migr06),nrow(migr06))
  pop08 <- pop[,c(5,15)]
  for(i in 1:nrow(migr06)){
    for(j in 1:nrow(migr06)){
      if(i != j){
        row <- which(toupper(pop08[,1]) == toupper(colnames(migr06))[j])
        popr08[i,j] <- pop08[row,2]
      }
    }
  }

  pops06 <- matrix(0,nrow(migr06),nrow(migr06))
  pop06 <- pop[,c(5,13)]
  for(i in 1:nrow(migr06)){
    for(j in 1:nrow(migr06)){
      if(i != j){
        row <- which(toupper(pop06[,1]) == toupper(colnames(migr06))[i])
        pops06[i,j] <- pop06[row,2]
      }
    }
  }

  pops07 <- matrix(0,nrow(migr06),nrow(migr06))
  pop07 <- pop[,c(5,14)]
  for(i in 1:nrow(migr06)){
    for(j in 1:nrow(migr06)){
      if(i != j){
        row <- which(toupper(pop07[,1]) == toupper(colnames(migr06))[i])
        pops07[i,j] <- pop07[row,2]
      }
    }
  }

  pops08 <- matrix(0,nrow(migr06),nrow(migr06))
  pop08 <- pop[,c(5,15)]
  for(i in 1:nrow(migr06)){
    for(j in 1:nrow(migr06)){
      if(i != j){
        row <- which(toupper(pop08[,1]) == toupper(colnames(migr06))[i])
        pops08[i,j] <- pop08[row,2]
      }
    }
  }

  #########################################################
  ## Make Unemployment Sender and Reciever ##

  uner06 <- matrix(0,nrow(migr06),nrow(migr06))
  unes06 <- matrix(0,nrow(migr06),nrow(migr06))
  une06 <- subset(unemp,unemp$Year == 2005)
  for(i in 1:nrow(migr06)){
    for(j in 1:nrow(migr06)){
      if(i != j){
        rowr <- which(toupper(une06$State) == toupper(colnames(migr06))[j])
        rows <- which(toupper(une06$State) == toupper(colnames(migr06))[i])
        uner06[i,j] <- une06[rowr,3]
        unes06[i,j] <- une06[rows,3]
      }
    }
  }

  uner07 <- matrix(0,nrow(migr06),nrow(migr06))
  unes07 <- matrix(0,nrow(migr06),nrow(migr06))
  une07 <- subset(unemp,unemp$Year == 2006)
  for(i in 1:nrow(migr06)){
    for(j in 1:nrow(migr06)){
      if(i != j){
        rowr <- which(toupper(une07$State) == toupper(colnames(migr06))[j])
        rows <- which(toupper(une07$State) == toupper(colnames(migr06))[i])
        uner07[i,j] <- une07[rowr,3]
        unes07[i,j] <- une07[rows,3]
      }
    }
  }


  uner08 <- matrix(0,nrow(migr06),nrow(migr06))
  unes08 <- matrix(0,nrow(migr06),nrow(migr06))
  une08 <- subset(unemp,unemp$Year == 2007)
  for(i in 1:nrow(migr06)){
    for(j in 1:nrow(migr06)){
      if(i != j){
        rowr <- which(toupper(une08$State) == toupper(colnames(migr06))[j])
        rows <- which(toupper(une08$State) == toupper(colnames(migr06))[i])
        uner08[i,j] <- une08[rowr,3]
        unes08[i,j] <- une08[rows,3]
      }
    }
  }
  #########################################################
  ## Make income Sender and Reciever ##

  incr06 <- matrix(0,nrow(migr06),nrow(migr06))
  incs06 <- matrix(0,nrow(migr06),nrow(migr06))
  for(i in 1:nrow(migr06)){
    for(j in 1:nrow(migr06)){
      if(i != j){
        rowr <- which(toupper(income$State) == toupper(colnames(migr06))[j])
        rows <- which(toupper(income$State) == toupper(colnames(migr06))[i])
        incr06[i,j] <- income[rowr,6]
        incs06[i,j] <- income[rows,6]
      }
    }
  }

  incr07 <- matrix(0,nrow(migr06),nrow(migr06))
  incs07 <- matrix(0,nrow(migr06),nrow(migr06))
  for(i in 1:nrow(migr06)){
    for(j in 1:nrow(migr06)){
      if(i != j){
        rowr <- which(toupper(income$State) == toupper(colnames(migr06))[j])
        rows <- which(toupper(income$State) == toupper(colnames(migr06))[i])
        incr07[i,j] <- income[rowr,5]
        incs07[i,j] <- income[rows,5]
      }
    }
  }

  incr08 <- matrix(0,nrow(migr06),nrow(migr06))
  incs08 <- matrix(0,nrow(migr06),nrow(migr06))
  for(i in 1:nrow(migr06)){
    for(j in 1:nrow(migr06)){
      if(i != j){
        rowr <- which(toupper(income$State) == toupper(colnames(migr06))[j])
        rows <- which(toupper(income$State) == toupper(colnames(migr06))[i])
        incr08[i,j] <- income[rowr,4]
        incs08[i,j] <- income[rows,4]
      }
    }
  }

  #########################################################
  ## Make the kilometer arc distance ##

  statedist <- matrix(0,nrow(migr06),nrow(migr06))
  for(i in 1:nrow(migr06)){
    for(j in 1:nrow(migr06)){
      rowi <- which(latlong[,1]==toupper(colnames(migr06)[i]))
      rowj <- which(latlong[,1]==toupper(colnames(migr06)[j]))
      lati <- latlong[rowi,2]
      latj <- latlong[rowj,2]
      longi <- latlong[rowi,3]
      longj <- latlong[rowj,3]
      statedist[i,j] <- deg.dist(longi,lati,longj,latj)
    }
  }

  # Make the 4 nearest neighbor matrix #

  nnmat <- NULL
  for(i in 1:ncol(statedist)){
    disti <- statedist[,i]
    disti[i] <- Inf
    nnmati <- rep(0,length(disti))
    nnmati[order(disti)[1:4]] <- 1
    nnmat <- cbind(nnmat,nnmati)
  }

  #########################################################
  ## Jan Temp Matrices ##

  jtemps <- matrix(0,nrow(migr06),nrow(migr06))
  for(i in 1:nrow(migr06)){
    for(j in 1:nrow(migr06)){
      if(i != j){
        row <- which(toupper(temps[,1]) == toupper(colnames(migr06))[i])
        jtemps[i,j] <- temps[row,2]
      }
    }
  }

  jtempr <- matrix(0,nrow(migr06),nrow(migr06))
  for(i in 1:nrow(migr06)){
    for(j in 1:nrow(migr06)){
      if(i != j){
        row <- which(toupper(temps[,1]) == toupper(colnames(migr06))[j])
        jtempr[i,j] <- temps[row,2]
      }
    }
  }

  #########################################################
  #Construct Z: An array of covariates which parameterize the latent space #

  Z <- array(0,dim=c(nrow(migr05),nrow(migr05),10))
  Z[,,1] <- 1

  Z[,,2] <- pops07-pops06

  Z[,,3] <- popr07-popr06

  Z[,,4] <- unes07-unes06

  Z[,,5] <- uner07-uner06

  Z[,,6] <- incs07-incs06

  Z[,,7] <- incr07-incr06

  Z[,,8] <- jtemps

  Z[,,9] <- jtempr

  Z[,,10] <- statedist


  #Standardize the covariates
  for(i in 2:10){
    Z[,,i] <- (Z[,,i]-mean(c(Z[,,i])))/sd(c(Z[,,i]))
  }
  #return migr07, migr06, and Z
  return(list(Graph = migr07 - migr06, Z = Z))
}
