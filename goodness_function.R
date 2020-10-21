
# viene de PAPER_MBA by solving the inverse ISING problem (ok)
# # # # # # # # # # # ## # # # ENERGY EVALUATION# ## # # # # ## # # # # ## # # # # ##  # # # # ## 
# FUNCION DE ENERGIA (O COSTO) A SER MINIMIZADA
# INPUT: vector de 1XN de spines (N=numero de nodos de la red)
#spins <- vector(mode = "numeric", length = node_numbers)
# OUTPUT: valor de la energia de la red para un spin dado y matiz de acople I dado.
#goodness <- function(spins, acoples) {
#  value = 0
#acoples <- I[upper.tri(I)]
#  rows <- length(spins)
#  id <- which(upper.tri(matrix(, rows, rows)) == TRUE, arr.ind=T)
#  n <- nrow(id)
#  for (i in seq_along(1:n)) {
#    value <- value + acoples[id[i,1], id[i,2]]*spins[id[i,1]]*spins[id[i,2]]
#print(value)
#  }
#  return(value)
#}
goodness <- function(spins, magnetizaciones, acoples) {
  simples <- t(magnetizaciones)%*%(as.matrix(spins)) #25-12-17
  value = 0
  rows <- length(spins)
  colu <- rev(abs(sequence(seq.int(rows - 1)) - rows) + 1)
  fila <- rep.int(seq.int(rows - 1), rev(seq.int(rows - 1)))
  id <- cbind(fila, colu)
  n <- dim(id)[1] #15-12-17
  for (i in seq_along(1:n)) {
    value <- value + acoples[id[i,1], id[i,2]]*spins[id[i,1]]*spins[id[i,2]]
    #print(value)
  }
  value <-  value + simples #25-12-17
  return(value)
}
# # # # # # # # # # # ## # # # ENERGY EVALUATION# ## # # # # ## # # # # ## # # # # ##  # # # # ## 





