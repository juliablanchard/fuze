# Simple size-structured model of salmon feeding on an unstructured resource (fed)
# JLB 31/10/2019

# Parameters:
ngrid   = 130                                  # No. of size classes          
tmax    = 10                                  # number of days
dt      = 0.1                                  # time step
isave   = 1                                    # how often to save results
wmin    = 10                                   # minimum size in g of grid
wmax    = 4000                                 #  max size weight g of grid
nsave  <- floor(tmax/(dt*isave)) # no. of time slots to save
# feeding
intake_coef = 8
intake_exp = 0.8
max_intake = 30
# intrinsic mortality
mort_coef = 0.01
mort_exp = -1/4
# metabolism
met_coef = 10
met_exp = 3/4
# assimilation efficiency
alpha  = 0.6
# egestion
egest =  1 - alpha

####----optional-------------

# Temperature dependence - same across all groups for the time being:
# c1=25.22                                      # constant used in Jennings et al. 2008 Proc B to standardize metabolism-temperature 
# # effects for Boltzmann equation. Derived from Simon's fit to Andy Clarke's data
# E=0.63                                        # activation energy, eV
# k=8.62*10^-5                                  # Boltzmann's constant

### -----------------------------

# Set up  matrices
  
w         <- 10^(seq(from=log10(wmin),to=log10(wmax),length.out=ngrid)) # creating sequence from 
dw        <- diff(w)
dw[ngrid] <- dw[ngrid-1] # Set final dw as same as one before

#intrinsic natural mortality
Z =  mort_coef*w^mort_exp

# senescence mortality rate to limit large fish from building up in the system
# same function as in Law et al 2008, with chosen parameters gives similar M2 values as in Hall et al. 2006

# Arrays to save output 
# Can add more if necessary

N = array(0,dim=c(nsave,ngrid))     # Abundance at size
gg = array(0,dim=c(nsave,ngrid))    # growth rates of functional size spectrum groups
harvest = matrix(0,nrow=nsave)      # Total number of removals of fisheries/harvesting by group
biomass = matrix(0,nrow=nsave)      # Total biomass across sizes by group in tank
e = array(0,dim=c(nsave,ngrid))     # amount of egested food
detritus = matrix(0,nrow=nsave)     # biomass of detritus

# feed input
rb = array(10,nsave)          # Resource biomass pools (feed)

# initial numbers and sizes at start of experiment

N[,1] = 20

# Matrices for solver
#A <- matrix(0,ncol=ngrid)
#B <- matrix(0,ncol=ngrid)
#S <- matrix(0,ncol=ngrid) 


########---------------------------
# Main loop to iterate over time, N [days]



for (i in 1:(nsave-1)) {

  # growth rate in grams per unit time at size
  gg[i,]<- (alpha*intake_coef*w^intake_exp)*rb[i] - met_coef*w^met_exp
  
  # egested food at size
  egested <-(1-alpha)*intake_coef*w^intake_exp*rb[i]
  
   
  # total detritus input 
  
  detritus[i] <- sum(egested)
  
  ### ---- integration scheme to get number at size in next time step
  # Iterate each predator size spectrum one time step forward:
  # Set up matrix:
  # A[,idx] <- -gg[,idx-1]*dt/dwmat[,idx]
  # B[,idx] <- 1 + gg[,idx]*dt/dwmat[,idx] + Z[,idx]*dt
  # S[,idx] <- N[,idx]
  ### ---------------------------------------
  
  idx = 2:ngrid
  
  N[i+1,idx] = N[i,idx] - (dt/dw[idx])*N[i,idx]*gg[i,idx] + (dt/dw[idx])*N[i,idx-1]*gg[i,idx-1] - dt*N[i,idx]*Z[idx]
  
  
  # current density + growth into size - growth out of size - mortality
 
  
    
}
 

#plot set up to track size spectrum over time
r = rainbow(nsave)
plot(log10(w), log10(N[100,]), type="l", col=r[1], cex=1.2,xlab = "log10(w)", ylab="log10(numbers)")

#plot size spectrum over time
for (i in 3:nsave) points(log10(w), log10(N[i,]), type="l", col=r[i], cex=1.2)


