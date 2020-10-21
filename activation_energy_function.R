
# # # # # # # # # # # ## # # # TRANSITION PROBABILITY # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
# Funcion para calcular la probabilidad de que un spin Si cambie de estado P(Si = 1) mientras
# que P(Si = 0) = 1 - P(Si = 1)
# Date: 05-12-17
# # # # # # # # # # # ## # # # TRANSITION PROBABILITY # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
#Inputs: 
#  vector de entrada = spins
#  J matriz de acoples = J
#  h vector de magnetizaciones = H
#  spin <- numero de neurona que se desea calcular la probabilidad
#Output: activation energy  

# 11-DIC-17 ACLARACION:
# Si Delta E = E(xi=1) - E(xi=0). Este Delta E es la energia de activacion.
# Si Delta E = E(xi=0) - E(xi=1). Este Delta E es la energia de desactivacion.
activation_energy <- function(H, acoples, spins, idspin) {
  # computo de energia de activacion:
  rows <- length(spins)
  id <- combn(1:rows,2)
  #N_par <- nrow(id) #numero de parametros
  N_par = dim(id)[1] #25-12-17
  #N_units <- nrow(acoples)
  N_units <- dim(acoples)[1] #25-12-17
  E <- H[idspin] #valor de la energia de activacion
  losj <- seq_along(1:N_units)
  losj <- losj[-(idspin)]
  for (j in losj) {
    E <- E + acoples[idspin,j]*spins[j]
  }
  return(E)
}
