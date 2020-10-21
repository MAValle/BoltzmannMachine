# viene de PAPER_MBA by solving the inverse ISING problem (ok)
# # # # # # # # # # # ## # # # METROPOLIS-HASTING # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
#INPUT: 
#   ACOPLES: matriz de acoples J del momento
#   initial_spin (pueden ser todos cero, todos 1, o aleatorios)
#   kmax = maximo numero de iteraciones sweep de MH
#   To = temperatura inicial. Si To = NULL, entonces T = 1. 
#   ks = numero de iteraciones adicionales en temperatura T = 1.
#OUTPUTS:
# matriz de spines
# vector de energias de cada sweep


metropolis.sampling <- function(magnetizaciones, acoples, initial_spin, kmax, ks, To=NULL, alpha) {
  output <- matrix(NA, ncol=length(initial_spin), nrow=kmax+ks)
  energies <- matrix(NA, ncol=2, nrow=kmax+ks)
  k <- 1
  vec <- initial_spin
  if ( is.null(To) ) { To = 1 ; alpha = 1 } 
  Temp <- To
  #sampleo annealing
  while (k <= kmax) {
    spines <- sweep(magnetizaciones=magnetizaciones, acoples = acoples, vector = vec, Temp = Temp)
    output[k, ] <- spines
    energies[k, ] <- c(k, goodness(spins = spines, magnetizaciones=magnetizaciones, acoples = acoples))
    vec <- spines
    Temp <- tempfunc(Temp, alpha = alpha)
    k <- k + 1
  }
  #sampleo con temperatura T = 1.
  if ( ks > 0 ) {
    while (k <= kmax + ks) {
      spines <- sweep(magnetizaciones=magnetizaciones, acoples = acoples, vector = vec, Temp=1)
      output[k, ] <- spines
      energies[k, ] <- c(k, goodness(spins = spines, magnetizaciones=magnetizaciones, acoples = acoples))
      vec <- spines
      #Temp <- tempfunc(Temp, alpha=0)
      k <- k + 1
    }
  }
  
  return(list(output, energies))
}
#ejemplo
#mh <- metropolis.sampling(acoples = I, initial_spin = rep(1,7), kmax=20000, ks=0, To=NULL, alpha=NULL)
# si quiero las ultimas Kcut muestras: muestras <- ss[[1]]; muestras <- muestras[-(1:Kcut), , drop=FALSE]
#ss[[2]] #energias al final de cdada sweep
# # # # # # # # # # # ## # # # METROPOLIS-HASTING # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
