#-----------------------------------------------------------------------------#
# Setup for executing the model of C. hyperboreus egg vertical distribution   #
# in Amundsen Gulf. Required for parallelization with the "foreach" package.  #
#                                                                             #
# Dufour K. 2016                                                              #
# Maps F.   2016                                                              #
#-----------------------------------------------------------------------------#


egg_sim <- function (tt, dt, i, ip, jp, p_prof, K_prof, T_prof, Mlonga, params) {
  
with( as.list( c(tt, dt, i, ip, jp, p_prof, K_prof, T_prof, Mlonga, params) ), {
  
  #--- Forcing vectors

  # Water density profile
  p_wat <- p_prof[, ip[i]]

  # Water turbulent kinetic energy profile
  K     <- K_prof[, ip[i]]
  
  # Water temperature profile
  temp  <- T_prof[, ip[i]]
  
  # M. longa female abundance profile
  mlon  <- Mlonga[, jp[i]]
  
  # Egg vertical velocity
  W     <- Wz(p_wat, p_egg, d_egg, g, mu)
  
  # Egg hatching rate profile (required since it varies with temperature)
  H     <- dv(temp, dev)

  forcing <- list( dt   = dt,
                   zi   = zi,
                   dz   = dz,
                   Nz   = Nz,
                   Ne   = Ne,
                   eggi = eggi,
                   ing  = ing,
                   mlon = mlon,
                   W    = W,
                   H    = H,
                   K    = K
                  )
  
  #--- Run simulations
 
  # First run WITH Intra-Guild Predation  
  outigp <- ode.2D( y=Ci, times=tt, func=degg, dimens=c(Nz,Ne), parms=forcing, method="ode45" )
  
  # Second run NO IGP
  forcing$ing <- c(a=0, b=0)  
  out0   <- ode.2D( y=Ci, times=tt, func=degg, dimens=c(Nz,Ne), parms=forcing, method="ode45" )
  
  #--- Compute IGP contribution to egg dynamics
  
  x0 <- rowSums( out0[, 2:dim(out0)[2]] )
  x1 <- rowSums( outigp[, 2:dim(outigp)[2]] )
  
  digp <- ( x0 - x1 ) / x0
  
  return(digp)
} )
}
