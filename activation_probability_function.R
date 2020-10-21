# # # # # # # # # # # ## # # # TRANSITION PROBABILITY # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
#input: 
#  energy = activation or desactivation energy
activation_probability <- function(energy, Temp) {
  p <- 1/( 1 + exp(-energy/Temp) )
  return(p)
}
# # # # # # # # # # # ## # # # TRANSITION PROBABILITY # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
# Ejemplo
#activation_probability(activation_energy(H=c(0,0,0,0,0,0,0), J=I, spins=c(0,1,1,0,1,0,0), idspin=2), T=NULL)
# # # # # # # # # # # ## # # # TRANSITION PROBABILITY # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
