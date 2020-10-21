# viene de PAPER_MBA by solving the inverse ISING problem (ok)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# FUNCION PARA determinar KUllback-Leibler tradicional.                                 #
# #https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence                    #
# Fecha creacion: 12-dic-17                                                             #
# Inputs:                                                                               #
# P distribucion real empirica                                                          #
# Q distribucion aproximada (ej de un metropolis hasting)                               #
# Si P es NULL entonces imputs adicionales: P_states y P_frequency                      #
# Output: KUllback-Leibler tradicional index                                            #
#  kullback2 es lo mismo pero se asume que P=NULL                                       #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# note:
# 19-sep-20 : a kullback3 le agrege simplemente el caso en que el KL sea infinito, es decir
#             no existe coincidencia de estados entre la distribucion P (real)
#             con la Q (del sampling).
# 21-sep-20: creamos una nueva manera de calcular la divergencia KL entre P (real - empirico)
#             y la Q (del sampling ) que puede ser mas rapido que las funciones anteriores.

kullback <- function(P=NULL, states=NULL, frequencies=NULL, Q) {
  if ( !is.null(P) ) {
    P_info <- frecuency(muestra = P )
    P_states <- na.omit(P_info[[1]])
    P_frequency <- na.omit(P_info[[2]])
  } else {
    P_states <- states
    P_frequency <- frequencies
  }
  Q_info <- frecuency(muestra = Q )
  Q_states <- na.omit(Q_info[[1]])
  Q_frequency <- na.omit(Q_info[[2]])
  
  # Politica:Todos los estados en Q_states que NO estan en P_states hay que borrarlos.
  #filas <- nrow(P_states )
  filas <- dim(P_states )[1] #25-12-17
  #nodos <- ncol(P_states )
  nodos <- dim(P_states )[2] #25-12-17
  # matriz con: estados de muestra original, frequency de pdf original (P_frequency), frequency de la pdf de sample (Q_frequency)
  M <- matrix(NA, ncol=nodos + 2, nrow=filas)
  for (i in seq_along(1:filas)) {
    vec <- P_states[i,]
    M[i,c(1:nodos)] <- vec
    M[i,nodos+1] <- P_frequency[i]
    id <- which(apply(Q_states, 1, function(x) identical(x[1:nodos], vec )))
    if (length(id) == 0) {
      M[i,nodos+2] <- 0
    } else {
      M[i,nodos+2] <- Q_frequency[id]
    }
  }
  require(entropy)
  # The Kullback-Leibler divergence from Q to P is often denoted DKL(P||Q).
  KL <- KL.empirical(M[,nodos+2], M[,nodos+1])
  
  return(KL)
}
kullback2 <- function(states=NULL, frequencies=NULL, Q) {
  
  P_states <- states
  P_frequency <- frequencies
  
  Q_info <- frecuency(muestra = Q )
  Q_states <- na.omit(Q_info[[1]])
  Q_frequency <- na.omit(Q_info[[2]])
  
  # Politica:Todos los estados en Q_states que NO estan en P_states hay que borrarlos.
  #filas <- nrow(P_states )
  filas <- dim(P_states )[1] #25-12-17
  #nodos <- ncol(P_states )
  nodos <- dim(P_states )[2] #25-12-17
  # matriz con: estados de muestra original, frequency de pdf original (P_frequency), frequency de la pdf de sample (Q_frequency)
  M <- matrix(NA, ncol=nodos + 2, nrow=filas)
  for (i in seq_along(1:filas)) {
    vec <- P_states[i,]
    M[i,c(1:nodos)] <- vec
    M[i,nodos+1] <- P_frequency[i]
    id <- which(apply(Q_states, 1, function(x) identical(x[1:nodos], vec )))
    if (length(id) == 0) {
      M[i,nodos+2] <- 0
    } else {
      M[i,nodos+2] <- Q_frequency[id]
    }
  }
  require(entropy)
  # The Kullback-Leibler divergence from Q to P is often denoted DKL(P||Q).
  KL <- KL.empirical(M[,nodos+2], M[,nodos+1])
  
  return(KL)
}
#Ejemplo: kl <- kullback(P=muestra, Q=mh[[1]])
# ejemplo: kullback(P=NULL, states=ff[[1]], frequencies=ff[[2]], Q=mh[[1]] )
# 22-dc-17
# Vemos que el calulo de KL cuando el numero de nodos es de 20, es extremamemnte largo 
# y por tanto impractico. Por lo tanto creamos kullback3. kullback3, calcula kl pero 
# para un porcentaje acumulado preestablecido (ej. 0.8) de la distribucion de P.
# Inputs : * P (la muestra real) -- opcional
#          * states de ff[1]] que son los estados ya cortados al umbral acumulado preestablecido de P
#          * frequencies  de ff[[2]] que es la frecuencia  de los estados de states de P.
#          * Q muestra en proceso.
kullback3 <- function(P=NULL, states=NULL, frequencies=NULL, Q, umbral=0.8) {
  if ( !is.null(P) ) { #si ya tenemos la frecuencia de P, utilizamos directamente states y frequencies
    P_info <- frecuency(muestra = P )
    P_states <- na.omit(P_info[[1]])
    P_frequency <- na.omit(P_info[[2]])
  } else {
    P_states <- na.omit(states)
    P_frequency <- na.omit(frequencies)
  }
  Q_info <- frecuency(muestra = Q )
  Q_states <- na.omit(Q_info[[1]])
  Q_frequency <- na.omit(Q_info[[2]])
  
  #Cortamos Q_states y Q_frequency al umbral preestablecido
  pacum <- cumsum(Q_frequency/nrow(Q))
  corte <- which.min(abs(pacum-umbral))
  corte <- corte + 1
  Q_states <-  Q_states[-c(corte:nrow(Q_states)), ]
  Q_frequency <- Q_frequency[-c(corte:length(Q_frequency))]
  
  # Politica:Todos los estados en Q_states que NO estan en P_states hay que borrarlos.
  #filas <- nrow(P_states )
  filas <- dim(P_states )[1] #25-12-17
  #nodos <- ncol(P_states )
  nodos <- dim(P_states )[2] #25-12-17
  # matriz con: estados de muestra original, frequency de pdf original (P_frequency), frequency de la pdf de sample (Q_frequency)
  M <- matrix(NA, ncol=nodos + 2, nrow=filas)
  for (i in seq_along(1:filas)) {
    vec <- P_states[i,]
    M[i,c(1:nodos)] <- vec
    M[i,nodos+1] <- P_frequency[i]
    id <- which(apply(Q_states, 1, function(x) identical(x[1:nodos], vec )))
    if (length(id) == 0) {
      M[i,nodos+2] <- 0
    } else {
      M[i,nodos+2] <- Q_frequency[id]
    }
  }
  if (sum(M[,nodos+2]) == 0 ) { # 19-sep-20: en caso que no exista ninguna coincidencia de estados entre Q y P.
    KL <- NA
  } else {
    #require(entropy)
    # The Kullback-Leibler divergence from Q to P is often denoted DKL(P||Q). P es la original, Q es la aproximada
    KL <- KL.empirical(M[,nodos+2], M[,nodos+1])
  }
  
  return(KL)
} 
#ejemplo:
#kl <- kullback3(states=states, frequencies=fr, Q=Q, umbral=1 )
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# 21-sep-20
# creamos una nueva funcion kullback que tiene la siguiente estrategia:
# En vez de tener como base de comparacion los 2^N estados del sistema, hacemos:
# 1. creamos una matriz P_states con todos los estados de la P (data empirica) con su respectiva
#   frecuencia. Por ejemplo, si tenemos N=20, habrian 2^20 estados, pero nosotros observamos
#   100 estado solamente. Entonces creamos una matriz identificando esos 100 estados y su
#   respectiva frecuencia.
# 2. Luego formamos otra matriz identificando los estados de Q (data sampleada de metropolis hasting)
# 3. Luego, agregamos otra columna a P_states que sera la frecuencia de ese estado en Q:
#       Para cada estado de P, buscamos si se encuentra en Q. Si no se encuentra, la frecuencia 
#       de la tercera columna sera cero, si se encuentra, se busca cuantas veces esta ese estado en Q.
# 4. Se calcula el KL basado en el conteo de estados de P y Q.
# Nos percatamos que siempre se toma como refeencia los estados de P, porque se trata de calcular 
# KL(Q||P) cuando difiere Q de P.
# inputs:
# states_P = estados unicos identificados de P (data empirica)
# frequencies_P = frecuencia identificadas de P para cada estado de P.
# sampling_Q = muestra de datos de la cual calculamos la distribucion Q
# umbral = porcentaje acumulado preestablecido (ej. 0.8) de la distribucion de P en caso que sean muchos estados.
# muestra = muestra de datos empiricos en caso que sea necesario
# outputs:
# KL = Kullback-Leibler
# NOTA: SE DEBE TENER INSTALADO PACKAGE philentropy
kullback4 <- function(muestra, states_P, frequencies_P, sampling_Q) {
  if (is.null(states_P) || is.null(frequencies_P)) {
    ff <- frecuencyV2(muestra = muestra ) # #calcula frecuencias de los estados de la muestra 
    states_P <- ff[[1]]
    frequencies_P <- ff[[2]]
    #sum(frequencies_P) 
  }
  ffQ <- frecuencyV2(muestra = sampling_Q )
  states_Q <- ffQ[[1]]
  
  # inicio del algoritmo
  frecuencia_Q <- integer(length=nrow(states_P)) # almacenamos la frecuencia de cada estado de Q, que esta en P
  for (i in 1:nrow(states_P) ) {
    #print(i)
    # buscar el estado states_P[i, ] en states_Q
    state_P <- states_P[i,]
    # buscar cuantas veces se repite state_P  en Q
    matches <- which(apply(states_Q, 1, function(x) return(all(x == state_P   ))))
    frecuencia_Q[i] <- length(matches)
  }
  #sum(frequencies_Q ) # 663

  # Utilizando library(philentropy)
  P <- frequencies_P/sum(frequencies_P, na.rm=TRUE)
  Q <- frecuencia_Q/sum(frecuencia_Q, na.rm=TRUE)
  x <- rbind(P,Q)
  
  value <- KL(x, unit = "log2") 
  return(value)
}
# ejemplo

# ff <- frecuencyV2(muestra = data ) # #calcula frecuencias de los estados de la muestra 
# states_P <- ff[[1]]
# frequencies_P <- ff[[2]]
# last_mh_sampling <- machine_results[[7]] 
# KL <- kullback4(muestra = NULL, states_P = states_P, frequencies_P = frequencies_P, sampling_Q = last_mh_sampling)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 




