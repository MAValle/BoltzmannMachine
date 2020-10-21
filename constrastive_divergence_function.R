# viene de PAPER_MBA by solving the inverse ISING problem (ok)
# # # # # # # # # # # # # # # # # CONTRASTIVE DIVERGENCE # # # # # # # # # # # # # # # # # # # 
#OUTPUTS:  h matriz de Nx1 con los parametros de magnetizacion
#          J matriz NxN con los acoples
# INPUTS: sigma_i_data  medias de los spines de la data empirica
#         sigma_ij_data correlaciones entre i y j de la data empirica
#         nu parametro de aprendizaje 
#         h, J parametros de magnetizacion y acoples en la epoca t-1
contrastive_divergence <- function(sigma_i_data, sigma_ij_data, sigma_i_sample, sigma_ij_sample, nu, h, J) {
  #require(gdata)
  h_old <- h
  J_old <- J
  #contrastive divergence
  h <- h_old + nu*(sigma_i_data - sigma_i_sample)
  J <- J_old + nu*(sigma_ij_data - sigma_ij_sample)
  #trackeo de convergencia de parametros h i J
  diff_hpar <- h - h_old
  diff_Jpar <- J - J_old
  diff_Jpar <- as.matrix(upperTriangle(diff_Jpar, diag = FALSE))
  #mean square of differences
  msd_hpar <- (t(diff_hpar)%*%diff_hpar)/length(diff_hpar) #  suma de las diferencias al cuadrado de hi
  msd_Jpar <- (t(diff_Jpar)%*%diff_Jpar)/length(diff_Jpar) # suma de las diferencias al cuadrado de los Jij
  cat(sprintf("hi square sum diffs %.7f -- Jij square sum diffs: %.7f \n", msd_hpar, msd_Jpar))
  
  
  #trackeo de convergencia de los sigma
  diff_h <- sigma_i_data - sigma_i_sample
  diff_J <- sigma_ij_data - sigma_ij_sample
  diff_J <- as.matrix(upperTriangle(diff_J, diag = FALSE))
  #mean square of differences
  msd_h <- (t(diff_h)%*%diff_h)/length(diff_h) #  suma de las diferencias al cuadrado de hi
  msd_J <- (t(diff_J)%*%diff_J)/length(diff_J) # suma de las diferencias al cuadrado de los Jij
  cat(sprintf("sigmas_i square sum diffs %.7f -- sigmas_ij square sum diffs: %.7f \n", msd_h, msd_J))
  
  
  return(list(h, J, msd_h, msd_J, msd_hpar, msd_Jpar))
}
#El mismo que antes, pero el error en las correlaciones entre el modelo y cada iteracion se calcula distinto.
contrastive_divergence2 <- function(sigma_i_data, sigma_ij_data, sigma_i_sample, sigma_ij_sample, nu, h, J) {
  #sigma_ij_data son los <sisj>+ de la fase posotiva, que en nuestra caso se calcula solo 1 vez a partir de la data de entranmiento
  #sigma_ij_sample son los <sisj>- de la fase negativa que viene de la funcion annealing.
  # require(gdata)
  h_old <- h
  J_old <- J
  #contrastive divergence
  h <- h_old + nu*(sigma_i_data - sigma_i_sample)
  J <- J_old + nu*(sigma_ij_data - sigma_ij_sample)
  #J <- as-matrix(J) #puesto el 031217
  
  #trakeo de <sisj>+ - <sisj>-
  dif <- sigma_ij_data - sigma_ij_sample
  N <- ncol(J)
  N <- N*(N+1)/2 - N
  trk <- matrix(NA, nrow=1, ncol=N)  #trakeo de <sisj>+ - <sisj>-
  idx <- combn(1:ncol(J),2)
  for (i in seq_along(1:N)) {
    trk[ 1, i ] <- dif[ idx[1,i], idx[2,i] ]
  }
  trk <- as.vector(trk)
  
  #trackeo de convergencia de parametros h i J
  diff_hpar <- h - h_old
  diff_Jpar <- J - J_old
  diff_Jpar <- as.matrix(upperTriangle(diff_Jpar, diag = FALSE))
  #mean square of differences
  msd_hpar <- (t(diff_hpar)%*%diff_hpar)/length(diff_hpar) #  suma de las diferencias al cuadrado de hi
  msd_Jpar <- (t(diff_Jpar)%*%diff_Jpar)/length(diff_Jpar) # suma de las diferencias al cuadrado de los Jij
  #cat(sprintf("hi square sum diffs %.7f -- Jij square sum diffs: %.7f \n", msd_hpar, msd_Jpar))
  
  #trackeo de convergencia de los sigma como the absolute correlation difference metric (paper de Broderick: Faster solutions of the inverse pairwise ising problem)
  sigmas_i_diff <- abs(sigma_i_data - sigma_i_sample)
  sigmas_ij_diff <- abs(sigma_ij_data - sigma_ij_sample)
  sigmas_ij_diff <- as.matrix(upperTriangle(sigmas_ij_diff, diag = FALSE))
  #mean 
  mean_sigmas_i_diff <- (mean(sigmas_i_diff))    # log(mean(abs(sigma_i_data - sigma_i_sample))
  mean_sigmas_ij_diff <- (mean(sigmas_ij_diff))  # log(mean(abs(sigma_ij_data - sigma_ij_sample))
  #cat(sprintf("sigmas_i abs means diffs %.7f -- sigmas_ij abs means  diffs: %.7f \n", mean_sigmas_i_diff, mean_sigmas_ij_diff))
  
  return(list(h, J, mean_sigmas_i_diff, mean_sigmas_ij_diff, msd_hpar, msd_Jpar, trk))
}
#contrastive_divergence3 aplica el siguiente criterio para calcular el nuevo J:
#si <sisj>+ > <sisj>- entonces Jij = Jij + nu,
#si <sisj>+ < <sisj>- entonces Jij = Jij - nu,
contrastive_divergence3 <- function(sigma_i_data, sigma_ij_data, sigma_i_sample, sigma_ij_sample, nu, h, J) {
  #sigma_ij_data son los <sisj>+ de la fase posotiva, que en nuestra caso se calcula solo 1 vez a partir de la data de entranmiento
  #sigma_ij_sample son los <sisj>- de la fase negativa que viene de la funcion annealing.
  #require(gdata)
  h_old <- h
  J_old <- J
  #contrastive divergence (segun paper Ackley 1985)
  tempJ <- ((sigma_ij_data > sigma_ij_sample)*1)*2-1
  J <- J_old + nu*tempJ
  tempH <- ((sigma_i_data > sigma_i_sample)*1)*2-1
  h <- h_old + nu*tempH
  #J <- as-matrix(J) #puesto el 031217
  
  #trakeo de <sisj>+ - <sisj>-
  dif <- sigma_ij_data - sigma_ij_sample
  N <- ncol(J)
  N <- N*(N+1)/2 - N
  trk <- matrix(NA, nrow=1, ncol=N)  #trakeo de <sisj>+ - <sisj>-
  idx <- combn(1:ncol(J),2)
  for (i in seq_along(1:N)) {
    trk[ 1, i ] <- dif[ idx[1,i], idx[2,i] ]
  }
  trk <- as.vector(trk)
  
  #trackeo de convergencia de parametros h i J
  diff_hpar <- h - h_old
  diff_Jpar <- J - J_old
  diff_Jpar <- as.matrix(upperTriangle(diff_Jpar, diag = FALSE))
  #mean square of differences
  msd_hpar <- (t(diff_hpar)%*%diff_hpar)/length(diff_hpar) #  suma de las diferencias al cuadrado de hi
  msd_Jpar <- (t(diff_Jpar)%*%diff_Jpar)/length(diff_Jpar) # suma de las diferencias al cuadrado de los Jij
  #cat(sprintf("hi square sum diffs %.7f -- Jij square sum diffs: %.7f \n", msd_hpar, msd_Jpar))
  
  #trackeo de convergencia de los sigma como the absolute correlation difference metric (paper de Broderick: Faster solutions of the inverse pairwise ising problem)
  sigmas_i_diff <- abs(sigma_i_data - sigma_i_sample)
  sigmas_ij_diff <- abs(sigma_ij_data - sigma_ij_sample)
  sigmas_ij_diff <- as.matrix(upperTriangle(sigmas_ij_diff, diag = FALSE))
  #mean 
  mean_sigmas_i_diff <- log(mean(sigmas_i_diff))    # log(mean(abs(sigma_i_data - sigma_i_sample))
  mean_sigmas_ij_diff <- log(mean(sigmas_ij_diff))  # log(mean(abs(sigma_ij_data - sigma_ij_sample))
  #cat(sprintf("sigmas_i abs means diffs %.7f -- sigmas_ij abs means  diffs: %.7f \n", mean_sigmas_i_diff, mean_sigmas_ij_diff))
  
  return(list(h, J, mean_sigmas_i_diff, mean_sigmas_ij_diff, msd_hpar, msd_Jpar, trk))
}
# # # # # # # # # # # # # # # # # CONTRASTIVE DIVERGENCE # # # # # # # # # # # # # # # # # # # 






# 22-sep-20
# creamos contrastive_divergence4 que es el mismo que contrastive_divergence2
# pero que tiene como salida adicional, el trakeo de #trakeo de <si>+ - <si>-
contrastive_divergence4 <- function(sigma_i_data, sigma_ij_data, sigma_i_sample, sigma_ij_sample, nu, h, J) {
  #sigma_ij_data son los <sisj>+ de la fase posotiva, que en nuestra caso se calcula solo 1 vez a partir de la data de entranmiento
  #sigma_ij_sample son los <sisj>- de la fase negativa que viene de la funcion annealing.
  # require(gdata)
  h_old <- h
  J_old <- J
  #contrastive divergence
  h <- h_old + nu*(sigma_i_data - sigma_i_sample)
  J <- J_old + nu*(sigma_ij_data - sigma_ij_sample)
  #J <- as-matrix(J) #puesto el 031217
  
  #trakeo de <si>+ - <si>-
  trk_si <- sigma_i_data - sigma_i_sample
  
  #trakeo de <sisj>+ - <sisj>-
  dif <- sigma_ij_data - sigma_ij_sample
  N <- ncol(J)
  N <- N*(N+1)/2 - N
  trk <- matrix(NA, nrow=1, ncol=N)  #trakeo de <sisj>+ - <sisj>-
  idx <- combn(1:ncol(J),2)
  for (i in seq_along(1:N)) {
    trk[ 1, i ] <- dif[ idx[1,i], idx[2,i] ]
  }
  trk <- as.vector(trk)
  
  #trackeo de convergencia de parametros h i J
  diff_hpar <- h - h_old
  diff_Jpar <- J - J_old
  diff_Jpar <- as.matrix(upperTriangle(diff_Jpar, diag = FALSE))
  #mean square of differences
  msd_hpar <- (t(diff_hpar)%*%diff_hpar)/length(diff_hpar) #  suma de las diferencias al cuadrado de hi
  msd_Jpar <- (t(diff_Jpar)%*%diff_Jpar)/length(diff_Jpar) # suma de las diferencias al cuadrado de los Jij
  #cat(sprintf("hi square sum diffs %.7f -- Jij square sum diffs: %.7f \n", msd_hpar, msd_Jpar))
  
  #trackeo de convergencia de los sigma como the absolute correlation difference metric (paper de Broderick: Faster solutions of the inverse pairwise ising problem)
  sigmas_i_diff <- abs(sigma_i_data - sigma_i_sample)
  sigmas_ij_diff <- abs(sigma_ij_data - sigma_ij_sample)
  sigmas_ij_diff <- as.matrix(upperTriangle(sigmas_ij_diff, diag = FALSE))
  #mean 
  mean_sigmas_i_diff <- (mean(sigmas_i_diff))    # log(mean(abs(sigma_i_data - sigma_i_sample))
  mean_sigmas_ij_diff <- (mean(sigmas_ij_diff))  # log(mean(abs(sigma_ij_data - sigma_ij_sample))
  #cat(sprintf("sigmas_i abs means diffs %.7f -- sigmas_ij abs means  diffs: %.7f \n", mean_sigmas_i_diff, mean_sigmas_ij_diff))
  
  return(list(h, J, mean_sigmas_i_diff, mean_sigmas_ij_diff, msd_hpar, msd_Jpar, trk, trk_si))
}