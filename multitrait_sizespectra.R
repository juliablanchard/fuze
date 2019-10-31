#-----Dynamic size spectra with multiple life-history traits --------------------------------
#
## 17/11/2008 KHA: added plotGrowth() and adjusted plotCommunity()
## 30/05/2008 last modified by JLB
# Program to model multiple size spectra with different life-history traits.
# Originally based on Pope et al. approach but with food dependent growth and other processes
# following approach of Ken Andersen, Martin Pedersen et al.
#
# 1. the life history trait that differs across size spectra is asymptotic body mass (Winf)
# 2. growth is food dependent and based on bioenergetics
# 3. type II functional response allows for feeding level to not exceed a maximum if individuals are satiated under high food conditions
#    - feeding level is determined from von Bertalanffy growth and is a function of body mass
# 4. starvartion mortality is included  - BUT NOT MODELLED AT PRESENT
# 5. feeding kernels are the same across species
# 6. reproduction can be modelled incorporating: !! CURRENTLY HELD CONSTANT!!

                                       #           - bioenergetics
#           - weight-based maturity schedules that differ for each spectra
# 7. dynamics of a "background" (plankton) size spectrum is modelled (using logistic equation)- CURRENTLY HELD CONSTANT
#
# APPLICATION: OPERATIONAL MODEL FOR EU IMAGE PROJECT
# work by Ken Andersen & Julia Blanchard
# Goal is to test response of ecosystem indicators to management scenarios using MSE in FLR
# Output of this code will be used for sampling and management model outside of this code
# This code is first crude attempt at an OM which can be used in preliminary MSE for indicators
#--------------------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Simulate size spectra
#  Input:
#   param  -- parameters which could be set by a call to baseparameters
#   bContinueCalculation -- true if the calculation should continue an old calculation
#                           given by "community".
#  Output:
#   community -- list with all output
# ------------------------------------------------------------------------------
RealCom <- function(param, bContinueCalculation=FALSE, community=NULL) {
    # load matlab libary to use handy logspace and zeros function
    library(matlab)

    attach(param)
    # ------------------------------------------------------------------
    #  Initial setup
    # ------------------------------------------------------------------
    sp = param$species
    
    #
    # Set up grid:
    #
    w = logspace(log10(w0), log10(wMax), ngrid)   #logarithmically spaced vector of mass
    dw = diff(w)
    dw[ngrid] = dw[ngrid-1]
    # wFull is same as above but with plankton spectrum added on
    ngridPP = round(ngrid)*0.3
    wFull= c(logspace(log10(w[1]/(3*mean(beta))),log10(w[1]-dw[1]),ngridPP),w)
    nGridfull<-length(wFull)
    dwFull = diff(wFull)
    dwFull[nGridfull] = dwFull[nGridfull-1]
    idxGrid=(ngridPP+1):nGridfull       # indexing size grid
    #
    # Set up plankton spectrum
    #
    nPP = kap*wFull^(-lambda)
    nPP[wFull>param$wPPcut]=0      #set density of sizes < plankton cutoff size to zero
    ##start with saturated plankton spectrum
    ##rrPP=rPP*wFull^lPP #weight specific plankton growth rate ##
    ##
    # Set up species
    #
    wmat = sp$Wmat # Size at maturation
    Z0 = Z0pre*sp$Winf^param$Z0exp   # background morality
    N = zeros(nspp,length(w))
    psi = zeros(nspp,length(w))
    F = zeros(nspp,length(w))
    IntakeMax = zeros(nspp,ngrid)
    SearchVol = zeros(nspp,ngrid)
    Activity = zeros(nspp,ngrid)
    StdMetab = zeros(nspp,ngrid)

    for (isp in 1:nspp){
        if (bContinueCalculation)
            N[isp,] = community$N[end,isp,]
        else {
            N[isp,] = sp$N0[isp] * (w/param$w0)^param$slope0
            ##make sure zero densities at w>Winf for each species
            N[isp,w>sp$Wmat[isp]]=0               
        }
        # maturity ogives are set such that the expected growth curves will be
        # von Bertalanffy-like:
        tmp=w/sp$Winf[isp]
        psi[isp,]=tmp^(1-n)*1/(1+(w/sp$Wmat[isp])^(-10))  
        psi[isp,w<0.1*sp$Wmat[isp]]=0
        psi[isp,tmp>1]=1

        ## Maximum intake
        IntakeMax[isp,] = sp$h[isp] * w^n 
        ## Search volume
        SearchVol[isp,] = param$gamma[isp] * w^q
        ## Bioenergetics
        Activity[isp,] = param$k[isp] * w
        StdMetab[isp,] = param$ks[isp] * w^p
    }
    ##
    ##Set up internal temporary calculation of feeding kernels
    ##
    predkernel=zeros(nspp,ngrid,length(wFull)) # lookup table for feeding kernel of each
                                        # pred/prey size
                                        # predators in rows,prey in columns
    for (i in 1:nspp)
        for(j in 1:ngrid) {# loop over pred sizes
            predkernel[i,j,] = exp(-0.5*(log(w[j]/(sp$beta[i]*wFull))/sp$sigma[i])^2)
            predkernel[i,j, find(wFull>=w[j])]=0
        } 
    ##
    ## Fishing:
    ##
    for (isp in 1:nspp)
        F[isp,]=F0[isp]*(1/(1+exp((wFstart[isp]-w)/xi*wFstart[isp])))
    #
    # Allocate arrays
    #
    A=zeros(nspp,ngrid)
    B=zeros(nspp,ngrid)
    S=zeros(nspp,ngrid)
    R=zeros(1,nspp)   #CHECK THIS NOT THE SAME AS KENS
    f=zeros(nspp,ngrid)
    #
    #   Initialise and setup arrays to store output
    #
    itimemax=tmax/dt  #max index of time array
    nsave=floor(tmax/(dt*isave)) #num of time slots to save

    M2 = as.single(zeros(nspp+1,length(wFull)))
    Z=zeros(ngrid,1)
    e=(zeros(nspp,ngrid))
    gg=(zeros(nspp,ngrid))
    Nsave=zeros(nsave,nspp,ngrid)
    Ntotsave=(zeros(nsave,ngrid))
    Rsave=(zeros(nsave,nspp))
    M2save=(zeros(nsave,nspp,length(wFull)))
    SSBsave=zeros(nsave,nspp)
    SSBmsave=zeros(nsave,nspp,ngrid)
    Rtotsave=zeros(nsave,1)                    
    fsave=(zeros(nsave,nspp,ngrid))
    gsave=(zeros(nsave,nspp,ngrid))
    MSsave=(zeros(nsave,nspp,ngrid))
    Biomass = zeros(nsave,nspp)
    phiprey = zeros(nspp,ngrid)
    Nppsave = zeros(nsave,length(wFull))
    #
    # Calculate total spectrum (first iteration)
    #
    Ntot=colSums(N)
    idx=2:ngrid
    #-----------------------------------------
    # Main loop
    #-----------------------------------------
    for (itime in 1:itimemax) {
        #
        # Calculate feeding level of all predators
        #
        for (i in 1:nspp) {
            ## Calculate effective specutrum of food:
            Neff = apply((theta[i,] %*% ones(1,ngrid)) * N, 2, sum)
            ## - integrate phi over all prey, to find expression
            ##  for the available food for predators j

            ## Available food:
            phiprey[i,] =
                as.vector((dwFull*wFull*nPP)%*%t(predkernel[i,1:ngrid,])+ # from resource
                          (dw*w*Neff)%*%t(predkernel[i,,idxGrid])) # from fish
            
            ## Encountered food
            encount = SearchVol[i,]*phiprey[i,]

            ## Calculate feeding level:
            f[i,]=encount/(encount + IntakeMax[i,])
        }
        ##
        ## Calculate predation by predator size j
        ##
        M2 = zeros(nspp,length(wFull))
        for (i in 1:nspp) {
            for (jpred in 1:nspp)
                M2[i,] = as.vector(
                  M2[i,] + theta[i,jpred]*
                  (dw*(1-f[jpred,])*SearchVol[jpred,]*
                   N[jpred,]) %*% predkernel[jpred,,])
        }
        ##
        ##  Iterate each species densities
        ##
        for (i in 1:nspp)    {
            #
            # Calculate the somatic growth:
            #
            
            # calculate assimilated food intake
            e[i,]=sp$alpha[i]*f[i,]*IntakeMax[i,]
            # Subtract basal metabolism and activity
            e[i,]=e[i,]-StdMetab[i,]-Activity[i,]
            e[i, e[i,]<0] = 0  ## Do not allow starvation
            #Starvation mortality ###NOT USED AT PRESENT  --need to check through
            ##Ms=zeros(ngrid,1)
            ##ix=e[i,]<0
            ##Ms(ix)=-param.Ms0*e[i,ix]/w[ix]
            ##e[i,ix]=0

            # Calculate the amount of energy allocated to gonads
            SSBm = psi[i,]*e[i,] 
            # subtract from assimilated food and use the rest for growth
            gg[i,] = e[i,]-SSBm
            #
            # Total mortality
            #
            Z = Z0[i] + M2[i,idxGrid] + F[i,]
            #
            # Iterate species on time step forward:
            #
            # Set up matrix:
            A[i,idx] = -gg[i,idx-1]*dt/dw[idx]
            B[i,idx] = 1+gg[i,idx]*dt/dw[idx]+Z[idx]*dt
            S[i,idx] = N[i,idx]

            # Boundary condition upstream end (recruitment)
            # SSB=sum(eRepro*SSBm*N[i,]*dw)
            # R =0.5*rho[i]*SSB/w[1]+R0[i]
            # Set to constant for now:
            B[i,1] = 1+gg[i,1]*dt/dw[1]+Z[1]*dt
            N[i,1]= N[i,1] #(N[i,1]+R*dt/dw[1])/B[i,1] #hold smallest size class constant

            # Invert matrix
            for (j in 2:ngrid)
                N[i,j] = (S[i,j]-A[i,j]*N[i,j-1]) / B[i,j]
            #
            # Save results:
            #
            if (mod(itime,isave)==0){
                isav=itime/isave
                fsave[isav,i,]=f[i,]
                gsave[isav,i,]=gg[i,]
                Nsave[isav,i,]=N[i,]
                Rsave[isav,i]=R[,i]
                Biomass[isav,i]=sum(N[i,]*w*dw) ##check units here
                M2save[isav,i,] = M2[i,]
            }
                                        #Eaten[isav,i]=sum(dw*f[isp,]*IntakeMax*N[isp,])
        }
        #
        # Calculate total community spectrum:
        #
        Ntot=colSums(N)
        # and save it
        if (mod(itime,isave)==0){ 
            isav=itime/isave
            Ntotsave[isav,]=Ntot
            Nppsave[isav,]=nPP
        }
    }
    #
    # Save results in a list:
    #
    community <- list(nSpecies = nspp,
                      w = w,
                      wFull = wFull,
                      Ntot = Ntotsave,
                      f = fsave,
                      N = Nsave,
                      M2 = M2save,
                      Biomass = Biomass,
                      Nresource = Nppsave,
                      g = gsave
                      )
                     
                                         

###------------------- OUTPUT FOR VIRTUAL SURVEY SAMPLING 

#needs to be converted to numbers.m^2
#Nsave has units number.m-3.g-1

#Nout=Nsave[end,,]                    #takes last time step in year 
#Num.out=sweep(Nout,2,dw,"*")         #convert to numbers per m-3
#plot(w,kap*w^(-lambda+1),log="xy")
#points(w,colSums(Num.out)) 
#lm(log(colSums(Num.out)[1:190])~log(w[1:190]))

#Num.out=Num.out*50 #multiply by average depth of North Sea: 50 metres to give numbers per m2

#community<-cbind(expand.grid(lapply(list(winf,w),as.numeric)),as.vector(Num.out), rep(winf.dat[,1],dim(Num.out)[2]), rep(winf.dat[,3],dim(Num.out)[2]), rep(winf.dat[,4],dim(Num.out)[2]),rep(winf.dat[,2],dim(Num.out)[2]))

#names(community)<-c("winf","weight.g","num.m-2","species","a","b","linf.orig")

#community$length.cm<-(community$weight.g/community$a)^(1/community$b)
    detach(param)
    return(community)
}


# ------------------------------------------------------------------------------
# Plot the output of a computation
# Input: community -- output of a computation
# ------------------------------------------------------------------------------
plotCommunity <- function(param, community) {
    attach(param)
    attach(community)
    sp = param$species
    iTend = dim(f)[1]
    xLimit = c(min(w), max(w))
    ## Index used for computing means
    idx = c(round(iTend/2):iTend)
    ## Vector of colors for different species
    colg=rainbow(12,start=0,end=0.8)
    #
    # Set up graphics:
    #
    par(mfrow=c(3,2))
    par(oma=c(2,2,2,2))
    par(mar=c(4,4,4,4))
    #                                    
    # Feeding level at final t
    #
    #meanf = apply(f[idx,param$nspp,], 2, mean)
    #plot(w, meanf,
    #     type="l",log="x",ylim=c(0,1),xlim=xLimit,
    #     xlab="weight (g)",ylab="Feeding level")
    #    points(w, apply(f[idx,nSpecies,], 2, min), type="l",lty=3)
    #points(w, apply(f[idx,nSpecies,], 2, max), type="l",lty=3)

    plot(w,w, type="n",log="x",ylim=c(0,1),xlim=xLimit,
         xlab="weight (g)",ylab="Feeding level")
    for (i in 1:param$nspp) {
        ix = w<=param$species$Winf[i]
        points(w[ix], apply(f[idx,i,ix], 2, mean), type="l", col=colg[i])
    }
    points(w, rep(f0[1],length(w)), type="l", lty=2)
    points(w, rep(fc[1],length(w)), type="l", lty=2)
    ## Mean feeding level:
    ##points(w, apply(f[idx,i,], c(2, mean), type="l", col=colg[i])
    
    #
    # Predation mortality
    #
    alphap = sqrt(2*pi)*param$gamma*sp$sigma*param$kap*
        exp((param$n-1)^2*sp$sigma^2/2)*sp$beta^(param$n-1)
    mu = (1-param$f0[1])*alphap[1]*wFull^(1+param$q-param$lambda)            #equ'n 10. in Andersen & Beyer
    meanM2 = apply(M2[idx,nSpecies,],2,mean)
    plot(wFull, meanM2, type="l",
         log="xy", xlim=xLimit, ylim=c(0.00001,1000),
         xlab="weight (g)",ylab="Predation mortality rate (1/year)")
    for (i in 1:param$nspp) {
        ix = wFull<=param$species$Winf[i]
        points(wFull[ix], apply(M2[idx,i,ix],2,mean), type="l", col=colg[i])
    }
    #points(wFull, apply(M2[idx,nSpecies,],2,min), type="l", lty=3)
    #points(wFull, apply(M2[idx,nSpecies,],2,max), type="l", lty=3)
    points(wFull,mu,type="l",lty=2)
    ##
    ## Community size spectrum
    ##
    plot(w, apply(Ntot[idx,],2,mean),
         type="l", log="xy", lwd=3, xlim=xLimit, ylim=c(1e-3,1e17),
         xlab="weight (g)", ylab="Abundance density (1/g/volume)")
         
    ## Species spectra:
    for(i in 1:nSpecies){
        ix = w<=param$species$Winf[i]
      points(w[ix], apply(N[idx,i,ix],2,mean),type="l",col=colg[i],lwd=2)
    }
    ## Reference spectrum:
    points(w, kap*w^-lambda, type="l", lty=2)
    ## Resource spectrum:
    points(community$wFull, apply(community$Nresource[idx,],2,mean),
           type="l", lty=2, col="green")
    ##
    ## Total biomass of species over time
    ##
    t = dt*c(1:iTend)*isave
    tonnes = 100000
    totBiomass = apply(Biomass, 1, sum)/tonnes
    plot(t, totBiomass,
         type="l",xlab="time (yrs)",ylab="Biomass (t/vol)",
         log="y", lwd=3, ylim=c(min(Biomass/tonnes),max(totBiomass)))
    for(i in 1:nSpecies)
        points(t, Biomass[,i]/tonnes, type="l",col=colg[i],lwd=2)
    ##
    ## Growth curves calculated using mean feeding level:
    ##
    require(odesolve)
    growth <- function(t,w,g) {
        return(list(approx(community$w,g,w)$y))
    }

    t = seq(0,30,by=0.1)
    plot(1,1, xlim=c(0,20),ylim=c(0,1.2),
         xlab="time Winf^(n-1)", ylab="Relative weight (w/Winf)",type="n")
    for (i in 1:nspp) {
        gg = apply(community$g[idx,i,],2,mean) ## Average growth rate
        gr = lsoda(param$w0, t, growth, gg)
        points(gr[,1], gr[,2]/sp$Winf[i],type="l",col=colg[i],cex=0.6,lwd=2) #*sp$Winf[i]^(param$n-1)
    }
   
   #Legend
    
      plot(1, 1, xlim=c(0,20),ylim=c(0,1.2),type="n", lwd=3, bty="n",axes="F",xlab="",ylab="")
            
      legend(0,1, sp$species, col=colg, lty=1,cex = 1,bty="n", ncol=2,lwd=2)
      
   
   
    #
    # Heading
    #
    #mtext(text=paste("gamma=",gamma," kappa=",kap,
    #      " beta=",sp$beta," sigma=",round(sp$sigma,1)," h=",sp$h," F=",F0[1]),outer=TRUE)

    detach(community)
    detach(param)
}
## --------------------------------------------------------
##  Plot growth curves of each species in a separate panel
## --------------------------------------------------------
plotGrowth <- function(param,c) {
    require(odesolve)
    ##
    ## Set up graphics:
    ##
    cols = floor(sqrt(param$nspp))
    par(mfrow=c(1+floor(param$nspp/cols),cols))
    par(mar=c(0.2,2.5,2,0.2))
    #par(mai=c(0.5,0.5,0.5,0.5))
    ##
    ## Growth function:
    ##
    growth <- function(t,w,g) {
        return(list(approx(c$w,g,w)$y))
    }
    ##
    ## Plots:
    ##
    t = seq(0,30,by=0.1)
    ## Index used for computing means
    iTend = dim(c$f)[1]
    idx = c(round(iTend/2):iTend) 
 
    for (i in 1:param$nspp) {
        gg = apply(c$g[idx,i,],2,mean) ## Average growth rate
        gr = lsoda(param$w0, t, growth, gg)
        plot(gr[,1], gr[,2],
             ylim=c(0,1.1*param$species$Winf[i]),
             xlab="age (yrs)", ylab="weight (g)",
             type="l") #, col=colg[i], cex=0.6,lwd=2)
        legend("bottomright", legend=param$species$species[i], bty="n")
        points(t, rep(param$species$Winf[i],length(t)), type="l", lty=2)
    }
}


#--------------------------------------------------------------------

# BASIC PARAMETERS FOR MULTITRAIT SIZE SPECTRA MODEL

# Basic units:
#  weight: grams
#  time: years

#--------------------------------------------------------------------
paramTraitModel <- function(nSpecies=2) {
    # load matlab libary to use handy logspace and zeros function
    library(matlab)
    param = list()
    #
    # Species asymptotic sizes:
    #
    param$nspp = nSpecies
    param$Winf = logspace(log10(10), log10(10^4), param$nspp) #discretised vector of Winf
    param$theta = ones(nspp, nspp) # Interaction matrix
    onez = rep(1,param$nspp)
    #
    # Numerical parameters:
    #
    param$ngrid = 100           # no. of size classes
    param$tmax = 1              # no. years
    param$dt = 0.02             # time step
    param$isave = 1             # how often to save results
    param$w0 = 0.001            # minimum size in grams
    param$wMax = 2*max(param$Winf)    # overall max weight of grid
    #
    # Expected feeding level:
    #
    param$f0est = 0.6          # equilibrium feeding level, for which h-bar was estimated
    #
    # Primary production:
    #
    param$kap=0.005           # Value of kappa (height of plankton) used in Kens code
                        # as a basis for the calculation of gamma (search vol)
                        # Here kappa is an input parameter rather than estimated
    param$wPPcut=min(param$Winf)           #cutoff size of background spectrum
    #
    # Morality:
    #
    param$Z0pre = 0.84          # intrinsic natural mort. prefactor
    param$Z0exp=-1/4            #exponent for intrinsinc nat mort
    #
    # Growth:
    #
    param$alpha=0.6           #assimilation efficiency
    param$n=3/4               #scaling of intake
    param$h=85*onez                #prefactor for intake
    param$k=0*onez                 #activity
    param$p=3/4               #scaling of std metabolism
    param$ks=10*onez               #std metabolism
    #
    # Encounter of food:
    #
    param$q=0.8                 #search volume exponent
    param$beta = 100*onez              #mean PPMR
    param$sigma = 1.5*onez             #width of size pref = stdev of ln(PPMR) in Gaussian function 
    param$gamma=4000*onez            #Prefactor of volumetric search rate
    param$lambda=2+param$q-param$n        # Exponent of resource spectrum
    param$alphae = sqrt(2*pi)*param$gamma*param$sigma*param$beta^(param$lambda-2)*
        exp((param$lambda-2)^2*param$sigma^2/2)
    #
    # Recruitment
    #
    param$alphaMature=0.25      # fraction of Winf to mature
    param$eRepro=0.20           # efficiency of gonad production
    ##
    ## Initial spectrum
    ##
    param$slope0 = -param$n - 0.5 # estimated initial slope
    param$N0 = 10*param$Winf^(-param$lambda-param$slope)
    ##
    # Fishing:
    #
    param$F0=0*ones(param$nspp,1)      # level of fishing mortality
    param$wFstart=0.05*param$Winf      # start of fished sizes
    param$xi=0.01
    #
    # 'Diagnostics'
    #
    param$fc = param$ks[1]/(param$alpha*param$h[1])     #critical feeding level - only enough food eaten to meet standard metabolism
    param$f0 = 1 / (1 + param$h[1]/(param$kap*param$alphae[1]))

    print(c("Background feeding level: ",param$f0))
    print(c("Critical feeding level: ",param$fc))

    return(param)
}
## ----------------------------------------------------
##  Run the trait-based model
## ----------------------------------------------------
runTraitmodel <- function() {
    param<-baseparameters()
    C<-RealCom(param)
    plotCommunity(param,C)
    return(C)
}
#--------------------------------------------------------------------

# BASIC PARAMETERS FOR a community models (i.e. neither traits nor species)

# Basic units:
#  weight: grams
#  time: years

#--------------------------------------------------------------------
paramCommunityModel <- function() {
    param = list()
    #
    # Species asymptotic sizes:
    #
    param$nspp = 1
    param$Winf = 100000
    onez = rep(1,param$nspp)
    #
    # Numerical parameters:
    #
    param$ngrid = 100           # no. of size classes
    param$tmax = 1              # no. years
    param$dt = 0.02             # time step
    param$isave = 1             # how often to save results
    param$w0 = 0.001            # minimum size in grams
    param$wMax = 2*max(param$Winf)    # overall max weight of grid
    #
    # Expected feeding level:
    #
    param$f0est = 0.6          # equilibrium feeding level, for which h-bar was estimated
    #
    # Primary production:
    #
    param$kap=0.005           # Value of kappa (height of plankton) used in Kens code
                        # as a basis for the calculation of gamma (search vol)
                        # Here kappa is an input parameter rather than estimated
    param$wPPcut=param$w0           #cutoff size of background spectrum
    #
    # Morality:
    #
    param$Z0pre = 0.84          # intrinsic natural mort. prefactor
    param$Z0exp=-1/4            #exponent for intrinsinc nat mort
    #
    # Growth:
    #
    param$alpha=0.6           #assimilation efficiency
    param$n=3/4               #scaling of intake
    param$h=85                #prefactor for intake
    param$k=0                 #activity
    param$p=3/4               #scaling of std metabolism
    param$ks=10               #std metabolism
    #
    # Encounter of food:
    #
    param$q=0.8                 #search volume exponent
    param$beta = 100*onez              #mean PPMR
    param$sigma = 1.5*onez             #width of size pref = stdev of ln(PPMR) in Gaussian function 
    param$gamma=4000            #Prefactor of volumetric search rate
    param$lambda=2+param$q-param$n        # Exponent of resource spectrum
    param$alphae = sqrt(2*pi)*param$gamma*param$sigma*param$beta^(param$lambda-2)*
        exp((param$lambda-2)^2*param$sigma^2/2)
    #
    # Recruitment
    #
    param$alphaMature=1      # fraction of Winf to mature
    param$eRepro=0.20           # efficiency of gonad production
    ##
    ## Initial spectrum
    ##
    param$slope0 = -param$lambda # estimated initial slope
    param$N0 = param$kap*param$w0^param$slope0
    #
    # Fishing:
    #
    param$F0=0*ones(param$nspp,1)      # level of fishing mortality
    param$wFstart=0.05*param$Winf      # start of fished sizes
    param$xi=0.01
    #
    # 'Diagnostics'
    #
    param$fc = param$ks/(param$alpha*param$h)     #critical feeding level - only enough food eaten to meet standard metabolism
    param$f0 = 1 / (1 + param$h/(param$kap*param$alphae[1]))

    print(c("Background feeding level: ",param$f0))
    print(c("Critical feeding level: ",param$fc))

    return(param)
}
#--------------------------------------------------------------------

# BASIC PARAMETERS FOR North Sea model

# Basic units:
#  weight: grams
#  time: years

#--------------------------------------------------------------------
paramNorthSeaModel <- function(fileSpecies="nsea_params.txt" ,fileInteraction="interaction_matrix.txt") {
    ## Read in species specific parameter from file:
    
    ## Species params file
    species = read.delim(fileSpecies)   
   
    ## Interaction matrix file
    imat = read.delim(fileInteraction) 
      
    # Set up parameters
    
    param=list(species=species, interaction=imat)
    param$theta = as.matrix(param$interaction[1:12,2:13])

    # Species asymptotic sizes:
    
    param$nspp = dim(species)[1]
    onez = rep(1,param$nspp)
    #
    # Numerical parameters:
    #
    param$ngrid = 100                   # no. of size classes
    param$tmax = 50                      # no. years
    param$dt = 1                     # time step
    param$isave = 1                     # how often to save results
    param$w0 = 0.001                    # minimum size in grams
    param$wMax = 2*max(species$Winf)      # overall max weight of grid
    #
    # Expected feeding level:
    #
    param$f0est = 0.6          # equilibrium feeding level, for which h-bar was estimated
    #
    # Primary production:
    #
    param$kap=1e12           # Value of kappa (height of plankton) used in Kens code
                        # as a basis for the calculation of gamma (search vol)
                        # Here kappa is an input parameter rather than estimated
    param$wPPcut=10           #cutoff size of background spectrum
    #
    # Mortality:
    #
    param$Z0pre = 0.84          # intrinsic natural mort. prefactor
    param$Z0exp=-1/4            #exponent for intrinsinc nat mort
    #
    # Growth:
    #
    param$n=3/4                                       #scaling of intake
    param$k=0*onez                                         #activity
    param$p=3/4                                       #scaling of std metabolism
    param$ks=0*onez                                       #std metabolism
    #
    # Encounter of food:
    #
    param$q=0.8                 #search volume exponent
    param$beta = species$beta   #mean PPMR
    param$sigma = species$sigma #width of size pref = stdev of ln(PPMR) in Gaussian function 
    param$gamma=50*1e-12*onez            #Prefactor of volumetric search rate
    param$lambda=2+param$q-param$n        # Exponent of resource spectrum
    param$alphae = sqrt(2*pi)*param$gamma*species$sigma*species$beta^(param$lambda-2)*
        exp((param$lambda-2)^2*species$sigma^2/2)
    #
    # Recruitment
    #
    param$alphaMature=0.25      # fraction of Winf to mature
    param$eRepro=0.20           # efficiency of gonad production
    ##
    ## Initial spectrum
    ##
    param$slope0 = -param$n - 0.5 # estimated initial slope
    
    ##
    # Fishing:
    #
    param$F0=0*onez      # level of fishing mortality
    param$wFstart=0.05*species$Winf      # start of fished sizes
    param$xi=0.01
    #
    # 'Diagnostics'
    #
    param$fc = param$ks/(species$alpha*species$h)     #critical feeding level - only enough food eaten to meet standard metabolism
    param$f0 = 1 / (1 + species$h/(param$kap*param$alphae[1]))

    print(c("Background feeding level: ",param$f0))
    print(c("Critical feeding level: ",param$fc))

    return(param)
}


paramNorthSeaModelold <- function() {
    sFile = "NorthSeaSpecies.dat"
    param = list()
    ##
    ## Read in species specific parameter from file:
    ##
    param$theta = ones(nspp, nspp)
    #
    # Species asymptotic sizes:
    #
    param$nspp = 10
    param$Winf = logspace(log10(10), log10(10^4), param$nspp) #discretised vector of Winf
    onez = rep(1,param$nspp)
    #
    # Numerical parameters:
    #
    param$ngrid = 100           # no. of size classes
    param$tmax = 5              # no. years
    param$dt = 0.2             # time step
    param$isave = 1             # how often to save results
    param$w0 = 0.001            # minimum size in grams
    param$wMax = 2*max(param$Winf)    # overall max weight of grid
    #
    # Expected feeding level:
    #
    param$f0est = 0.6          # equilibrium feeding level, for which h-bar was estimated
    #
    # Primary production:
    #
    param$kap=0.005           # Value of kappa (height of plankton) used in Kens code
                        # as a basis for the calculation of gamma (search vol)
                        # Here kappa is an input parameter rather than estimated
    param$wPPcut=param$w0           #cutoff size of background spectrum
    #
    # Mortality:
    #
    param$Z0pre = 0.84          # intrinsic natural mort. prefactor
    param$Z0exp=-1/4            #exponent for intrinsinc nat mort
    #
    # Growth:
    #
    param$alpha=0.6           #assimilation efficiency
    param$n=3/4               #scaling of intake
    param$h=85                #prefactor for intake
    param$k=0                 #activity
    param$p=3/4               #scaling of std metabolism
    param$ks=10               #std metabolism
    #
    # Encounter of food:
    #
    param$q=0.8                 #search volume exponent
    param$beta = 100*onez              #mean PPMR
    param$sigma = 1.5*onez             #width of size pref = stdev of ln(PPMR) in Gaussian function 
    param$gamma=4000            #Prefactor of volumetric search rate
    param$lambda=2+param$q-param$n        # Exponent of resource spectrum
    param$alphae = sqrt(2*pi)*param$gamma*param$sigma*param$beta^(param$lambda-2)*
        exp((param$lambda-2)^2*param$sigma^2/2)
    #
    # Recruitment
    #
    param$alphaMature=0.25      # fraction of Winf to mature
    param$eRepro=0.20           # efficiency of gonad production
    ##
    ## Initial spectrum
    ##
    param$slope0 = -param$n - 0.5 # estimated initial slope
    param$N0 = 500*param$Winf^(-param$lambda-param$slope)
    ##
    # Fishing:
    #
    param$F0=0*ones(param$nspp,1)      # level of fishing mortality
    param$wFstart=0.05*param$Winf      # start of fished sizes
    param$xi=0.01
    #
    # 'Diagnostics'
    #
    param$fc = param$ks/(param$alpha*param$h)     #critical feeding level - only enough food eaten to meet standard metabolism
    param$f0 = 1 / (1 + param$h/(param$kap*param$alphae[1]))

    print(c("Background feeding level: ",param$f0))
    print(c("Critical feeding level: ",param$fc))

    return(param)
}
## ----------------------------------------------------
##  Run the trait-based model
## ----------------------------------------------------
runTraitmodel <- function() {
    param<-paramTraitModel()
    Com<-RealCom(param)
    plotCommunity(param,Com)
    return(Com)
}
## ----------------------------------------------------
##  Run the North Sea model
## ----------------------------------------------------
runNorthSeaModel <- function() {
    param<-paramNorthSeaModel()
    Com<-RealCom(param)
    plotCommunity(param,Com)
    return(Com)
}
