#=======================================================================
# Julia L. Blanchard's code for a functional group size spectrum
# 15/08/2017
# Contents:
#   1. Model description
#   2. Parameter setting function
#   3. Setting up the model using the parameters function
#   4. Projecting through time function
#   5. Plots and visualisation
#=======================================================================

rm(list=ls())
setwd("/Users/juliab6/Google Drive/Git/FUZE")
list.files()

# install required packages (if not already installed)
list.of.packages <- c("ggplot2", "dplyr", "reshape2")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# ----------------------------------------------------------------------
# 1. Model description
# -----------------------------START------------------------------------

# FUZE -  FUnctional groups within a siZE spectrum

# Models dynamical size spectra in different functional feeding groups as well
#   as incorporating unstructured biomass pools for some groups if needed (eg. detritus, macroalage).
# Originally based on benthic-pelagic size spectrum model but made generic to facilitate collaborative 
#   developments and application to real ecosystems - see separate code for fitting models to data
# Models can be forced with time series of external environmental drivers ( output from physical biogeochemical models)

# Three functions:
#   1. Param - read in model parameters and structure
#   2. Setup - takes the model parameter list and creates the objects used in the projection
#   3. Model - projects forward in time

# ------------------------------END-------------------------------------


# ----------------------------------------------------------------------
# 2. Setting the Parameters of the model
# -----------------------------START------------------------------------

Param <- function(fileGroups="groupparams.csv",fileCoupling="couplingstrengths.csv",ENVinput=NA){
	
  # Read in functional group parameters from file:
  groups = read.csv(fileGroups)
  strengths = read.csv(fileCoupling)
  row.names(strengths) <- strengths[,1]
  #drivers = read.csv(ENVinput)
  
  # Set up parameter list
  param=list( 
            # Model Structure:
              groups = groups,                                # functional group parameters
              theta  = strengths[-1],                         # coupling strengths, or interaction matrix for groups
 
              
            # Grid Structure:
              nss     = length(grep("Spectrum",groups$Type)), # Number of size spectrum groups
              npred   = length(grep("Predator",groups$Type)), # Number of predator groups
              nbio    = length(grep("Pool",groups$Type)),     # Number of unstructured biomass pool groups
              ngrid   = 130,                                  # No. of size classes          
              tmax    = 100,                                  # no. years
              dt      = 0.1,                                  # time step
              isave   = 1,                                    # how often to save results
              w0      = min(groups$w0),                                # minimum size in g of grid
              wMax    = max(groups$wmax),                                #  max size weight g of grid
             
             
            # Resource Spectrum input (power law scaling of number versus body mass, kap is prefactor, lambda is exponent)
             kap = 0.1,
             lambda = -2    
             
            # Temperature dependence - same across all groups for the time being:
              # c1=25.22,                                       # constant used in Jennings et al. 2008 Proc B to standardize metabolism-temperature 
                                                              # # effects for Boltzmann equation. Derived from Simon's fit to Andy Clarke's data
              # E=0.63,                                         # activation energy, eV
              # k=8.62*10^-5                                   # Boltzmann's constant

            # Complexity-size dependent prey vulnerability:
            
            # Environmental Input:
            #  env.drivers =  "drivers"         
          )
  return(param)
}
param <- Param()

# ------------------------------END-------------------------------------


# ----------------------------------------------------------------------
# 3. Creating the objects for the projection from the parameters
# -----------------------------START------------------------------------

Setup <- function(param) {

  # Set up size grid,  note that producer and consumer spectra are combined on this grid
  
  w         <- 10^(seq(from=log10(param$w0),to=log10(param$wMax),length.out=param$ngrid)) # creating sequence from 
  dw        <- diff(w)
  dw[param$ngrid] <- dw[param$ngrid-1] # Set final dw as same as one before
 
  
  # Might keep seperate for future...
  ngridPP <- 30
  
  idxGrid <- (ngridPP+1):param$ngrid # shortcut to index just the dynamic consumer/predator Size Spectra, perhaps not needed here

  nsave  <- floor(param$tmax/(param$dt*param$isave)) # no. of time slots to save
  
  nss <- param$nss
  
  npred <-param$npred
  
  ngrid <- param$ngrid
  
  nbio <- param$nbio

  model =  list(
    
      param = param,

      # Grid parameters
      w = w,
      
      dw = dw,
  
     # Functional Group parameters - for size spectra components only
   
      SearchVol = matrix(NA,nrow=nss,ncol=ngrid),
      
      predkernel = array(NA,dim=c(nss,ngrid,ngrid)),      

      Neff = array(NA,dim=c(nss,ngrid,npred)),
   
      #intrinsic natural mortality
      Z0 =  matrix(NA,nrow=nss,ncol=ngrid),

      # senescence mortality rate to limit large fish from building up in the system
      # same function as in Law et al 2008, with chosen parameters gives similar M2 values as in Hall et al. 2006

      SM = matrix(NA,nrow=nss,ncol=ngrid),
      
      # Arrays to save output 
      # Can add more if necessary

      N = array(0,dim=c(nsave,nss,ngrid)),      # Abundance density spectra
      BP = array(0,dim=c(nsave,nbio)),          # Biomass pools
      Fm = array(0,dim=c(nsave,nss,ngrid)),     # Fishing/harvesting mortality
      # RDI = matrix(0,nrow=nsave,ncol=nspp),   # Recruitment (NEED TO REVISIT THIS IF COEXISTENCE BECOMES AN ISSUE), reporduction off at moment
      M2 = array(0,dim=c(nsave,nss,length(w))), # predation mortality of functional size spectrum groups
      gg = array(0,dim=c(nsave,nss,ngrid)),     # growth rates of functional size spectrum groups
      Yield = matrix(0,nrow=nsave,ncol=(nss + nbio)),    # Total Yields of fisheries/harvesting by group
      Biomass = matrix(0,nrow=nsave,ncol=(nss + nbio)) # Total biomass across sizes by group
     
      )
      
      # Initial population abundance - moved down
      #if (ContinueCalculation == T) model$N[1,,] = initialcommunity$N[dim(initialcommunity$N)[1],,]
      #if (ContinueCalculation == F) model$N[1,,] <- unlist(tapply(w,1:length(w),function(wx,N0,w0,slope0,Wmax) N0 * (wx/w0)^slope0,N0=10^5,w0=w[1],slope0=-2))                                           
     
      # end with(param) 
       
       
      return(model)
		
}

model<-Setup(param)
# ------------------------------END-------------------------------------


# ----------------------------------------------------------------------
# 4. Projecting forward
# -----------------------------START------------------------------------

Project <- function(model,temp.effect=F,complexity.effect = F,ContinueCalculation=FALSE, initialcommunity=NULL) {
	
with (model, {

      
      # Shortcuts
      grp <- param$groups
      dt <- param$dt
      theta <- as.matrix(param$theta)
      
      # max index of array 
      itimemax <- param$tmax/dt 
      
      N <- model$N
      BP <- model$BP
             
       # Matrices for solver
       A <- matrix(0,nrow=nss-1,ncol=ngrid)
       B <- matrix(0,nrow=nss-1,ncol=ngrid)
       S <- matrix(0,nrow=nss-1,ncol=ngrid)       
             
            
	# fixed SearchVol,  predkernal arrays and non-predation mortality terms     

      SearchVol[] <- unlist(tapply(w,1:ngrid,function(wx,A,gamma) A*wx^gamma,gamma=param$groups[grep("Spectrum",param$groups$Type),"gamma"],A=param$groups[grep("Spectrum",param$groups$Type),"A"]))

      # Could improve this...dims are: group x pred sizes x prey sizes
      predkernel[] <- param$groups[grep("Spectrum",param$groups$Type),"beta"]
      predkernel[] <- exp(-0.5*sweep(log(sweep(sweep(predkernel,3,model$w,"*")^-1,2,model$w,"*")),1,param$groups[grep("Spectrum",param$groups$Type),"sigma"],"/")^2)
      
    #  image(predkernel[5,,])
    # predkernel[] <- sweep(predkernel,c(1,2),combn(model$w,1,function(w,w0)w>w0,w0=param$groups[grep("Spectrum",param$groups$Type),"w0"]),"*") # find out the untrues and then multiply
      
    # other mortality   
     Z0[] <- unlist(tapply(w,1:ngrid,function(wx,m0) m0*wx^-0.25,m0=param$groups[grep("Spectrum",param$groups$Type),"mu0"]))

     Z0[] <- sweep(Z0,c(1,2),combn(model$w,1,function(w,w0)w>w0,w0=param$groups[grep("Spectrum",param$groups$Type),"w0"]),"*") # find out the untrues and then multiply by 0
     
    # senescence mortality  
      
     SM[] <- unlist(tapply(w,1:ngrid,function(wx,wmat,zspre,zsexp) zspre*(wx/wmat)^zsexp,wmat=param$groups[grep("Spectrum",param$groups$Type),"wmat"],zspre=param$groups[grep("Spectrum",param$groups$Type),"zspre"],zsexp=param$groups[grep("Spectrum",param$groups$Type),"zsexp"]))

     SM[] <- sweep(SM,c(1,2),combn(model$w,1,function(w,wmat)w>wmat,wmat=param$groups[grep("Spectrum",param$groups$Type),"wmat"]),"*") # find out the untrues and then multiply by 0
     
     # Set up helpful matrices to avoid using sweep statements
     # wmat <- matrix(w,byrow=T,nrow=nss,ncol=ngrid)
     # dwmat <- matrix(dw,byrow=T,nrow=nss,ncol=ngrid)
     
      wmat <- matrix(w,byrow=T,nrow=nss,ncol=length(w))
      dwmat <- matrix(dw,byrow=T,nrow=nss,ncol=length(w))

		  
	   # Initial population abundances 
      
      if (ContinueCalculation == T) {
      	N[1,,] = initialcommunity$N[dim(initialcommunity$N)[1],,]
      	BP[1,] = initialcommunity$BP[dim(initialcommunity$BP)[1],]
      	}
      	
      if (ContinueCalculation == F) {
      	N[1,,] <- unlist(tapply(w,1:length(w),function(wx,N0,w0,slope0,Wmax) N0 * (wx/w0)^slope0,N0=10^5,w0=param$groups[grep("Spectrum",param$groups$Type),"w0"],slope0=-2))                                           
      
      BP[1,] <- 0.01      # arbitrary value !!
	  }
	  
	  
	  
	  # MAIN TIME LOOP
      for (itime in 1:itimemax)
      {	
      	
      
         # Need to stucture this part by each different "Type"
         
                	
         # We have 3 different types of groups: ResourceSpectrum, BiomassPool, PredatorSpectrum and ConsumerSpectrum
             
         	
         # For now the ResourceSpectrum is held constant at intial values (could be semi-chemostat or forced)     
         
         # Set up resource (primary producer) spectrum
          
         N[itime, grep("ResourceSpectrum",param$groups$Type),] <- param$kap*w^(-param$lambda) # the resource carrying capacity - one for each mp and m (130 of them)
         N[itime, grep("ResourceSpectrum",param$groups$Type),w>10] <- 0      #set density of sizes < plankton cutoff size = 10 g
        
         # if (ContinueCalculation==TRUE)
         # model$nPP[1,] <- initialcommunity$nPP[dim(initialcommunity$nPP)[1],]
         # else
         # model$nPP[1,] <- model$NinfPP # plankton numbers - set the initial values at carrying cap

         
         # The PredatorSpectra are size-based predation and use the predkernel for feeding - driving growth and death. They are also predated by the Consumer BiomassPools
             	   	
      	 # This is the numbers of prey in particular group that are exposed to each predator (coupling strengths interaction matrix)
       
         # integrate feeding kernel over all prey, total available food at size for each predators 
         
         Neff <- theta%*%N[itime,,]  # use coupling matrix
 
        
        # Now, with the feeding kernel to determine total biomass of suitably sized prey, across all groups, by predator
 
      	phiprey <- predkernel[predgrp,,]%*%colSums(Neff[,,predgrp]*dwmat*wmat)
  
        # Encountered food      
         	
      	encount <- SearchVol[predgrp,]*phiprey
      	  
      	  
        # Predation mortality - predator size spectra only
         
        # Pred kernel is how much each group (1) by mass (2), predates on other masses (3)
 	      
 	     predrate <- sweep(predkernel[predgrp,,],c(1,2),model$SearchVol*N[itime,predgrp,]*dwmat,"*")
      	  
        #### STOPPED HERE
         
        # Calculate growth
      	 
      	  
      	gg <- K*encount
      	  
      	  # could incorporate feeding level here or handling time for FR Type II, for now it's a Type I FR
         
      	  
      	  # Then apply theta to predation mortality
          
        M2 <- theta[grep("Spectrum",groups$Type),grep("Spectrum",groups$Type)] %*% colSums(aperm(predrate, c(2,1,3)),dims=1)[,model$idxGrid]
   
   
       } # end predgrp loop
          
        otherMort <- Z0 + SM
          
          # Fishing mortality rates
          
          # Fm at time t and size w for each size spectrum group - using minimum and maximum size fished
        Fmort <- sweep(model$selectivity,c(1,3),sweep(model$param$Q,2,model$effort[itime,],"*"),"*")
          
        # Total mortality - could be faster to make Z0 a matrix to make addition simple
        # Need to sum Fmort over gears and add all mrotality terms
        Z = sweep(M2 + rowSums(Fmort,dims=2),1,otherMort,"+")
       
      	  # Iterate each predator size spectrum one time step forward:
          # Set up matrix:
           A[,idx] <- -gg[,idx-1]*dt/dwmat[,idx]
           B[,idx] <- 1 + gg[,idx]*dt/dwmat[,idx] + Z[,idx]*dt
           S[,idx] <- N[,idx]


		  # Invert matrix
 			for (i in 1:nss)
      			for (j in (w0idx[i]+1):ngrid)
        			N[i,j] <- (S[i,j] - A[i,j]*N[i,j-1]) / B[i,j]

       
 
         # NOT DONE:The ConsumerSpectra (e.g.non-size based feeding detritivores,herbivores but that are eaten by PredatorSpectrum compete for shared food and grow in size)   
                
                
         # Biomass pools  - # Need to separate into consumers and producers (later also ectotherms and endotherms as in Yodzis and Innes?)
         
         
         # Consumer pool has mean adult body size and need to feed on size range of prey to get inputs and outputs form other mortality and disease (starvation?)
         
         
         # food intake, single feeding kernal based on fixed (max) size and beta and sigma
         
         Bsearch <- param$group$A*param$group$wMax^param$group$gamma
         
         Bpredkernel <- exp(-0.5*sweep(log(sweep(sweep(predkernel,3,wFull,"*")^-1,2,param$group$wMax,"*")),1,groups[grep("Spectrum",groups$Type),"sigma"],"/")^2)
   
         input<-sum(Bsearch*Bpredkernel)
         
         # fixed mortality rate, no predation
         
         output <-0.2
         
         # NOT DONE: Resource pool has nutrient fluxes in and experiences size-based mortality but not predation. If detritus the flexes in come from egestion.faeces and dying material.
         
                 
         
         
         
                	
     








  # Save results:
    if ((itime %% param$isave)==0)
    {
              isav<-itime/param$isave
              model$gg[isav,,]<- gg
              model$N[isav,,]<-N
              model$Fm[isav,,]<-rowSums(Fmort,dims=2) # sum Fmort over gears
              model$M2[isav,,] <- M2
              # model$BP <- ???????????????
              model$Yield <- rowSums(sweep((model$Fmort * model$Ntemp),3,w*dw,"*"),dims=2)
              model$Biomass <-  rowSums(sweep(model$Ntemp,3,w*dw,"*"),dims=2)         
                          
                          
                          
    }
  }  #end of for time loop

  # Calculate these for easy reference  
  model$Yield[] <- 
  model$Biomass[] <- 







	
# end with(model)
 })	
	
	
	
	
	
	return(model)
	
}

# ------------------------------END-------------------------------------


# ----------------------------------------------------------------------
# 5. Plots and visualisation
# -----------------------------START------------------------------------



# ------------------------------END-------------------------------------