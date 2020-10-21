# viene de PAPER_MBA by solving the inverse ISING problem (ok)
# sweep function
#input: vector de estados del momentos (spines) y matriz de acoples J
sweep <- function(magnetizaciones, acoples, vector, Temp) {
  N_nodos <- length(vector)
  v <- vector
  for (i in seq_along(1:N_nodos)) {
    logico <- as.logical(v)
    #energy_activation <- activation_energy(H=rep(0, N_nodos), acoples=acoples, spins=v, idspin=i)
    energy_activation <- activation_energy(H=magnetizaciones, acoples=acoples, spins=v, idspin=i) #25-12-17
    Pon <- activation_probability(energy_activation, Temp=Temp)
    Poff <- 1 - Pon
    u <- runif(1)
    if (logico[i]) {
      if (Poff > u) { 
        logico[i] <- !logico[i] 
      } 
    } else {
      if ( Pon > u) { 
        logico[i] <- !logico[i] } 
    }
    v <- 1*logico
    #print(v)
  }
  return(v)
}