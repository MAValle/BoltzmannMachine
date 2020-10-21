# viene de PAPER_MBA by solving the inverse ISING problem (ok)
# # # # # # # # # # # # # # # # # CONDICIONES INCIALES # # # # # # # # # # # # # # # # # # # 
#04-ago-17
# Alternativas de condicionales iniciales de los acoples y/o magnetizaciones - tipe:
# a1 = todos igual a cero.
# a2 = todos igual al vector de medias para h0 y a la mtriz de correlaciones para J
# a3 = las magnetizaciones iniciales h0 son igual a log(pi/(1-pi)) donde pi es la proporcion
#      de spines 1 en la variable i. Los acoples J0 igual a la matriz de correlaciones
# a4 = las magnetizaciones iniciales son iguales a log(pi((1-pi))) donde pi es la prob. que 
#      spin este activada (Si=1). Los acoples son iguales a 0. Se necesita como input la 
#      la matriz de datos 
# a5 = la matriz de acople son valores aleatorios de una gaussiana con media igual a 0 y 
#      desv. estandar igual a desv. y las magnetizaciones todas igual a cero.
# 06 = la matriz de acople son dados de alguna estimacion anterior, al igual que el vector
#       de campos.

#INPUTS:
#  * randomon = obtiene un valor aleatorio de un valor de una gaussiana con media igual al valor 
#  y desv. estandar igual a 0.05 (on, off)
#  * si y sij = vector sigma_i_data y matriz sigma_ij_data de magnetizaciones y correlaciones de la data
#  * df = data original con spines 1 y -1.
#  * N = numero de nodos o variables
# * desv = desviacion estandar en caso que tipe = a5

#OUTPUTS:
#  * ho y J0 : vector y matriz inicial de magnetizacion y acoples.
initial_conditions <- function(df=NULL, N, tipe, si, sij, desv, hg = NULL, Jg = NULL) {
  require(Matrix)
  if (tipe == 'a1') {
    h0 <- t(as.matrix(matrix(0L, nrow=1, ncol=N)))
    J0 <- matrix(0L, nrow=N, ncol=N) 
  } else if (tipe == 'a2') {
    h0 <- si
    J0 <- sij
  } else if (tipe == 'a3') {
    df = 0.5*df+0.5
    p = colSums(df)/nrow(df)
    h0 <- log(p/(1-p))
    J0 <- matrix(0L, nrow=N, ncol=N) 
  } else if (tipe == 'a4'){
    h0 <-matrix(0L, nrow=1, ncol=N)
    for (cl in c(1:N)) {
      s <- sum(df[,cl] == 1)
      p <- s/nrow(df)
      h0[1, cl] <- log(p/(1-p))
    }
    h0 <- t(as.matrix(h0))
    J0 <- matrix(0L, nrow=N, ncol=N) 
  } else if (tipe == 'a5') {
    h0 <- t(as.matrix(matrix(0L, nrow=1, ncol=N)))
    J0 <- matrix(runif(N*N, -1, 1),N )
    ind <- lower.tri(J0) 
    J0[ind] <- t(J0)[ind] 
    diag(J0) <- NA
  }
  else if (tipe == 'a6') {
  h0 <- as.numeric(hg)
  J0 <- Jg
  diag(J0) <- NA
  }
  return(list(h0, J0))
}
# # # # # # # # # # # # # # # # # CONDICIONES INCIALES # # # # # # # # # # # # # # # # # # # 
