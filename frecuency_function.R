# viene de PAPER_MBA by solving the inverse ISING problem (ok)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# FUNCION PARA determinar la frecuencia de los estados de una muestra                   #
# Fecha creacion: 03-dic-17                                                             #
# Fecha modificacion:                                                                   #
# Input: matriz muestra (posiblemente datos simulados de metropolis hasting o datos de compra)                                           #                     
# Output: matriz de todos los estados encontrados en la muestra                         #
# Outputs:  vector con la fecuencia de cada uno de los estados encontrados              #
# Nota: si hay NAs en las salidas es porque No todos los estados estaban presentes      #
# en la muestra                                                                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
frecuency <- function(muestra) {
  if (!is.matrix(muestra)) {
    muestra <- as.matrix(muestra)
  }
  #N = nrow(muestra) #muestra es una matriz con lso estados, simulados de un metropolis hasting
  N = dim(muestra)[1] #25-12-17
  #n = ncol(muestra)
  n = dim(muestra)[2] #15-12-17
  n_states = 2^n
  frqq <- matrix(NA, nrow=n_states, ncol=n+1) #colocamos todos los posibles estados y la frecuencia
  
  it = 1
  matriz_states <- muestra
  #while ( !(is.null(nrow(matriz_states) ) ) | (nrow(matriz_states) != 0)  ) { 
  while ( nrow(matriz_states) != 0) { 
    selected_state <- matriz_states[1,]
    row.is.a.match <- apply(matriz_states, 1, identical, selected_state)
    total.matches <- sum(row.is.a.match)
    frqq[it, ] <- c(selected_state, total.matches) #gradamos estado y su frecuencia
    match.idx <- which(row.is.a.match)
    #ahora borramos todos los selected_state de matriz_state
    matriz_states <- matriz_states[-match.idx, ]
    if ( !(is.matrix(matriz_states)) ) {
      matriz_states <- t(as.matrix(matriz_states))
    }
    #print(nrow(matriz_states))
    it <- it + 1
  }
  frqq  <- frqq[order(-frqq[,n+1]), ]
  states <- frqq[,c(1:n)] #matriz con los estados contados
  counts <- frqq[, n+1] #matriz con el conteo de cada estado
  
  return(list(states, counts))
}
#ejemplo:
#ff <- frecuency(muestra = mh[[1]]) #muestra generada del metropolis hasting
#ff[[1]] #matriz con los estados contados
#ff[[2]] #vector con el conteo de cada estado
# frecuencyV2 corre marginalmente un poco mas rapido que frecuency.
frecuencyV2 <- function(muestra) {
  if (!is.matrix(muestra)) {
    muestra <- as.matrix(muestra)
  }
  # chequeo si muestra viene con columna newid
  if ("newid" %in% colnames(muestra))
  {
    borrar <- which( colnames(muestra)== "newid" )
    muestra <- muestra[, -c(borrar)]
  }
  
  #N = nrow(muestra) #muestra es una matriz con lso estados, simulados de un metropolis hasting
  N = dim(muestra)[1] #25-12-17 numero de filas
  #n = ncol(muestra)
  n = dim(muestra)[2] #15-12-17 numero de variables
  
  # identificamos los estados de muestra
  unique_wb <- as.matrix(unique(muestra[,c(1:ncol(muestra))])) 
  #n_states = 2^n
  n_states <- nrow(unique_wb)
  frqq <- matrix(NA, nrow=n_states, ncol=n+1) #colocamos todos los posibles estados y la frecuencia
  
  # comenzamos:
  matriz_states <- muestra
  it = 1
  while ( nrow(matriz_states) != 0) { 
    selected_state <- matriz_states[1,]
    row.is.a.match <- apply(matriz_states, 1, identical, selected_state)
    total.matches <- sum(row.is.a.match)
    frqq[it, ] <- c(selected_state, total.matches) #gradamos estado y su frecuencia
    match.idx <- which(row.is.a.match)
    #ahora borramos todos los selected_state de matriz_state
    matriz_states <- matriz_states[-match.idx, ]
    if ( !(is.matrix(matriz_states)) ) {
      matriz_states <- t(as.matrix(matriz_states))
    }
    #print(nrow(matriz_states))
    it <- it + 1
  }
  frqq  <- frqq[order(-frqq[,n+1]), ]
  states <- frqq[,c(1:n)] #matriz con los estados contados
  counts <- frqq[, n+1] #matriz con el conteo de cada estado
  
  return(list(states, counts))
}
#Ejemplo
#ff <- frecuencyV2(muestra = mh[[1]]) 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 