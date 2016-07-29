###############################################################################
#-----------------------------------------------------------------------------#
# Model of C. hyperboreus egg vertical distribution in Amundsen Gulf          #
#                                                                             #
# Dufour K. 2015-2016                                                         #
# Maps F.   2015-2016                                                         #
#-----------------------------------------------------------------------------#
###############################################################################

require(deSolve)

require(ReacTran)

require(R.matlab)

require(doParallel)

require(foreach)

#------------------------------#
# Import data                  #
#------------------------------#

#--- Vertical physical profiles
profile <- readMat("Stations_full_5m.mat")

#--- Vertical distribution of M. longa females
Mlonga  <- readMat("Mlonga_CFL_profile.mat")$C6f


#------------------------------#
# Environmental variables      #
#------------------------------#

#--- Reshape data in regular 5m bins

dz  <- 5                             # Layer thickness (m)
zi  <- seq(2.5, 300, by=dz)          # Distance from the center of a layer to the next
Nz  <- length(zi)                    # Number of vertical layers

#--- Physical profiles data start at the third bin! Fill in the gap assuming homogeneous surface layer.

p_prof        <- matrix(0, nrow=Nz, ncol=dim(profile$stationRho)[2]) # density anomaly profile (g m^-3)
p_prof[3:Nz,] <- (profile$stationRho-1e3) * 1e3
for(i in 1:2) {
  p_prof[i,]  <- p_prof[3,]
}
for(i in 1:dim(p_prof)[2]) {
  j <- max(which(is.finite(p_prof[,i])))
  if(j<Nz) {
    p_prof[(j+1):Nz,i] <- p_prof[j,i] 
  }
}

#if(FALSE) { 
K_prof        <- 0*p_prof                # turbulence profile (m^2 h^-1)
K_prof[3:Nz,] <- profile$stationK * 3600
K_prof[2,]  <- K_prof[3,]
for(i in 1:dim(K_prof)[2]) {
  j <- max(which(is.finite(K_prof[,i])))
  if(j<Nz) {
    K_prof[(j+1):Nz,i] <- K_prof[j,i] 
  }
}
#}

T_prof        <- 0*p_prof               # temperature profile (deg C)
T_prof[3:Nz,] <- profile$stationCT
for(i in 1:2) {
  T_prof[i,]  <- T_prof[3,]
}
for(i in 1:dim(T_prof)[2]) {
  j <- max(which(is.finite(T_prof[,i])))
  if(j<Nz) {
    T_prof[(j+1):Nz,i] <- T_prof[j,i] 
  }
}


#------------------------------#
#            MODEL             #
#------------------------------#

#--- Model ODE

degg <- function (t, C, pars)
  
{ with( as.list( c(t,C,pars) ), {
  
  C[C<1e-12] <- 0 # get rid of extremely low values for stability
  
  CC <- matrix(nrow=Nz, ncol=Ne, data=C) # reshape data as a 2D matrix
  
  # Advection rates
  w  <- dt * c(0,W) # Egg velocity profile
  v  <- dt * H      # Egg development profile
  
  # Diffusion rate
  k  <- dt * c(K,0) # Vertical eddy diffusivity
  
  # Reaction terms (biology)
  
  # 1) Egg production
  egg_p     <- matrix(0, nrow=Nz, ncol=Ne) # vertical profile; only contributes to first egg development stage
  egg_p[,1] <- dt/24 * eggi * dnorm(zi, mean=250, sd=15) # hourly egg production rate
  
  # 2) Egg ingestion
  egg_ing <- dt * ig(rowSums(CC), ing)/max(1e-9,rowSums(CC)) * mlon # ingestion by M. longa females
  egg_ing <- matrix(rep(egg_ing,Ne), nrow=Nz, ncol=Ne)            # spread over every egg development stage
  
  # 3) Egg hatching (development)
  #    Use a flux limiter scheme to prevent numerical dispersion between development stages
  #    Vectorize the 2D matrix to speed-up computation
  Ch <- as.vector( rbind( t(CC), rep(0,Nz) ) )        # add an extra empty stage, then ...
  
  vh <- as.vector( rbind( rep(0,Nz), matrix(v, nrow=Ne, ncol=Nz, byrow=T) ) )
  
  dH <- advection.1D(C=Ch, C.up=0, C.down=0, v=vh, dx=1/Ne, adv.method="quick")$dC
  dH <- matrix(dH, nrow=Nz, ncol=Ne+1, byrow=T)[,1:Ne] # ... get rid of the extra empty stage
  
  # Rate of change dC/dt
  
  # Vertical advection-diffusion
  dC <- tran.2D(C=CC, C.x.up=0, C.x.down=0, C.y.up=0, C.y.down=0, v.x=w, v.y=0, D.x=k, D.y=0, dx=dz, dy=1/Ne)$dC
  
  # Add reaction terms
  dC <- dC + dH + egg_p - egg_ing * C
  
  return(list(dC))
} )
}


#------------------------------#
# Model parameters             #
#------------------------------#

#--- Physical functions

# Vertical eddy diffusivity
# Use Bourgault et al. 2011 (doi:10.1029/2011GL047936) 
# function if K vertical profile NOT provided
Kz    <- function (z, K, K0=1.1e-2*3600, delta=-0.17) {
  
  if (is.vector(K)) {
    k <- K
  } else {
    k     <- K+vector("double",length(z))
    
    i0    <- z<10
    k[i0] <- K0*exp(delta*10)
    
    i1    <- z<=46 & z>=10
    k[i1] <- K0*exp(delta*z[i1])
  }
  return(k)
}

#--- Biological functions

# Egg development rate (h^-1)
dv    <- function (temp, dev) { 1/(dev['a']*(temp-dev['alpha'])^dev['b']) / 24 }

# Egg ingestion rate by M. longa females (egg fem^-1 h^-1)
ig    <- function (egg, ing) {
  i <- exp(ing['a']+ing['b']*egg*1e-3) # [egg] m^-3 -> L^-1
  if( sum(ing)==0 ) { i <- 0 }
  return(i)
}

# Egg velocity according to Stokes' law (m h^-1)
Wz    <- function (p_wat, p_egg, d, g, mu) { g*d^2*(p_egg-p_wat)/(mu*18) * 3600 }

#--- Physical parameters

g     <- 9.80665        # gravitational acceleration (m s^2) 
mu    <- 19*1e-3*100    # Dynamic molecular viscosity (g m^-1 s^-1) Visser & Jonasdottir 1999
K     <- 3.4e-6*3600    # Constant mean vertical eddy diffusivity below the surace layer (m^-2 h^-1) Bourgault et al. 2011 doi:10.1029/2011GL047936

#--- Biological parameters

Ne    <- 20             # number of age classes within the egg stage
eggi  <- 3e4            # initial condition from average daily egg production in february-march (d^-1 m^-2) Darnis 2013 
d_egg <- 192*1e-6       # Egg diameter (m) Jung-Madsen et al. 2013 doi:10.4319/lo.2013.58.6.0000
p_egg <- 19.4e3         # Egg density anomaly (g m^-3) Jung-Madsen et al. 2013 doi:10.4319/lo.2013.58.6.0000 range = c(0.6, 26.8)e3

dev   <- c(a=1196, alpha=-12.7, b=-2.05) # Parameters for the egg development Belehradek relationshiop to temperature (d^-1) (Jung-Madsen et al. 2013 doi:10.4319/lo.2013.58.6.0000; their Expt. 1)

ing   <- c(a=-2.965, b=0.054)            # Parameters for the functional response of females M. longa (this study Dufour et al. 2016)

#--- Parameters vector

params <- list( Ne    = Ne,
                Nz    = Nz,
                g     = g,
                mu    = mu,
                d_egg = d_egg,
                p_egg = p_egg,
                K     = K,
                dev   = dev,
                ing   = ing,
                eggi  = eggi
              )


#------------------------------#
# Simulations                  #
#------------------------------#

#--- Temporal discretization

dt <- 12                                # 8-hour timestep
tt <- seq(1, 15*24/dt, by=1)           # 15-day simulation to reach an equilibrium

#--- Initial conditions

Ci     <- matrix(0, nrow=Nz, ncol=Ne)
Ci[,1] <- eggi*dnorm(zi, mean=250, sd=15) # Initial vertical profile of egg concentration

#--- Choose a simulation: physical properties profile x profile of M. longa female density

ip <- rep( seq(1,20), 18)
jp <- matrix( seq(1,18), nrow=20, ncol=18, byrow=T); dim(jp) <- 20*18

#--- parallelized loop

source('egg_profile_func_2D.R')

ncpu <- detectCores()
registerDoParallel(cores=ncpu) # specify number of cores to use

ptime <- system.time({ # this is to check the time of execution
digp <- foreach(i=1:(20*18), .combine='cbind', .packages=c('deSolve','ReacTran')) %dopar% {    
  egg_sim(tt, dt, i, ip, jp, p_prof, K_prof, T_prof, Mlonga, params)
}
})[3]

registerDoSEQ()                # quit the parallel mode

#------------------------------#
# Post-processing              #
#------------------------------#

#--- Exploratory figures

i <- which.max(digp[dim(digp)[1],])

image(tt*dt/24, 1:(20*18), digp)

plot(tt*dt/24, digp[, i], type="l")

hist(digp[dim(digp)[1], ]*100, xlim=c(0,4))

#--- Save results

writeMat('egg_profile_v2_D0.mat', digp=digp, i=i, ip=ip[i], jp=jp[i])
