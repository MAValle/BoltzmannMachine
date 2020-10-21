

# # # # # # # # # # # ## # # # BOLTZMANN-LEARNING # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
# Creation date: 19-may-18    (viene de PAPER_MBA by solving the inverse ISING problem (ok))                                                                     #
# Creamos la funcion bmachine (utilizada en BM_toy_asonam.R) para llevar a cabo                   #
# el Boltzmann-learning datos transaccionales.                                                    #
# Bolstzamnn machine sacado de boltzmann_machine_realdata_V1.R                                    #
# INPUTS:                                                                                         #
# condInitial: puede ser -a2- para la funcion initial_condition                                   #
# muestra: base de datos con columnas (productos) y filas (compras) en formato 0 y 1 matrix       #
# Kmax: Numero maximo de iteraciones de la maquina de Boltzmann                                   #
# parametros aprendizaje                                                                          #
# LearningRate: parametro de aprendizaje 0.9 el original funcionando                              #
# decay: decaimiento del parametro de aprendizaje (en 0.02)                                       #
# N_sa: Numero de muestras para la fase negativa (numero de veces que repetimoes el annealing)    #
# umbral: umbral para el calculo del KL                                                           #
# every: cada cuantas iteraciones medir el KL.                                                    #
# metrop_its: numero iteraciones de metropolis hasting                                            #
# OUTPUTS:                                                                                        #
# just_for_track: trackeo de el ajuste <sisj>+ - <sisj>-)                                         #
# diver : divergencia KL                                                                          #
# cd: diferencias <sisj>+ - <sisj>-                                                               #
# inferences[Kmak, ] : parametros inferidos y los sigmas_{ij} y Hi                                #
# # # # # # # # # # # ## # # # BOLTZMANN-LEARNING # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
bmachine <- function(condInitial, muestra, Kmax, metrop_its, LearningRate, decay, N_sa, umbral, every, hg, Jg) {
  library(progress)
  #require(ggplot2)
  #library(latex2exp)
  #library(ggthemes)
  # ADQUISICION DE ALGUNOS PARAMETROS PREVIOS
  # original: ff <- frecuency(muestra = muestra ) #calcula frecuencias de los estados de la muestra 
  ff <- frecuencyV2( muestra = muestra ) # 21-sep-20
  
  sigma_i_sample <-  as.matrix(colMeans(muestra)) #15-abr-18
  sigmas_sample <- as.matrix(pairproducts(muestra))
  #ic <- initial_conditions(df=NULL, N=ncol(muestra), tipe=condInitial, si=sigma_i_sample, sij=sigmas_sample, desv=0.5) 
  ic <- initial_conditions(df=NULL, N=ncol(muestra), tipe=condInitial, si=sigma_i_sample, sij=sigmas_sample, desv=0.5, hg, Jg) 
  H0 <- ic[[1]] 
  J0 <- ic[[2]]
  numero_mediciones <- Kmax/every
  # generacion de vector senal para medir kl
  if (every == 1) {
    senal <- rep(TRUE, Kmax)
  } else {
    inic <- c(TRUE, rep(FALSE, every-2), TRUE)
    v <- c(rep(FALSE, every-1), TRUE)
    vv <- rep(v, numero_mediciones-1)
    senal <- c(inic, vv)
    rm(inic, v, vv)
  }
  
  # PPRESETEO DE VARIABLES Y MATRICES
  #comenzamos:
  pb <- progress_bar$new(format = "  Progress [:bar] :percent eta: :eta", total = Kmax, clear = FALSE, width= 60) 
  # Annealing network
  initial_spins=rep(1,ncol(sigmas_sample))
  J <- as.matrix(J0)
  H <- as.matrix(H0)
  epoch <- 1
  just_for_track <- matrix(NA, ncol=3, nrow=Kmax) #col1:epoch, col2:log(mean(abs(<sisj>+ - <sisj>-)), col3:suma de las diferencias al cuadrado de los Jij
  numbers_of_nodes = ncol(muestra)
  numbers_of_coupl = (numbers_of_nodes*(numbers_of_nodes+1)/2) - numbers_of_nodes
  #matriz inferences: epoch de la BA / Jij / hi. El orden es J12, J13, J23, J14, J24, J34, J15, J25, J35, etc..
  inferences <- matrix(NA, nrow=Kmax, ncol=numbers_of_nodes+numbers_of_coupl+1) #guardando los parametros h y J de cada epoch de la BM.
  #matriz de las diferencias <sisj>+ - <sisj>- para cada epoch
  cd <-  matrix(NA, nrow=Kmax, ncol=numbers_of_coupl+1)
  #matriz de las diferencias <si>+ - <si>- para cada epoch
  cd2 <-  matrix(NA, nrow=Kmax, ncol=ncol(muestra)+1)
  #matriz de los KUllback-Leibler Divergence
  diver <- matrix(NA, nrow=Kmax, ncol=2) # numero de epoch / KL
  
  # MAQUINA
  while (epoch <= Kmax) {
    # progress bar
    pb$tick()
    Sys.sleep(1 / Kmax)
    # progress bar fin 
    # Iniciamos la fase negativa: network is run free by using metropolis hasting
    sampling <- metropolis.sampling(magnetizaciones=H, acoples = J, initial_spin = initial_spins, kmax=metrop_its, ks=0, To=NULL, alpha=NULL)
    # en caso de que quisiera operar la maquina sin las magnetizaciones.
    #sampling <- metropolis.sampling(magnetizaciones=as.matrix(rep(0, ncol(muestra))), acoples = J, initial_spin = initial_spins, kmax=20000, ks=0, To=NULL, alpha=NULL)
    sigmas <- as.matrix(pairproducts( sampling[[1]] ))
    rescue <- contrastive_divergence4(sigma_i_data = sigma_i_sample, 
                                      sigma_ij_data = sigmas_sample, 
                                      sigma_i_sample = as.matrix( colMeans( sampling [[1]]) ), 
                                      sigma_ij_sample = sigmas,
                                      nu=learning_par(ciclo = epoch, decay = decay, LearningPar = LearningRate),
                                      h=H,
                                      J=J)
    J <- as.matrix(rescue[[2]])
    H <- rescue[[1]]
    put_all <- c(epoch, rescue[[4]],  rescue[[6]]) #epoch / log(mean(abs(<sisj>+ - <sisj>-)) / msd_Jpar
    just_for_track[epoch,] <- put_all
    #Colocamos los acoples Jij y hi en matriz inferences
    inferences[epoch, ] <- c(epoch, upperTriangle(J, diag = FALSE), rescue[[1]])
    #colocamos las diferencias <sisj>+ - <sisj>- de cada epoch
    cd[epoch, ] <- c(epoch, rescue[[7]])
    cd2[epoch, ] <- c(epoch, rescue[[8]])
    #calculamos Jensen-Shannon divergence y la colocamos en matrix diver
    # vamos a ir calculando kl solo cuando sea necesario y no en cada iteracion. Dato de entrada -every- cada cuanto medimos
    if (senal[epoch] == TRUE) {
      kl <- kullback4(muestra = NULL, states_P = ff[[1]], frequencies_P = ff[[2]], sampling_Q = sampling[[1]]) #21-sep-20
      diver[epoch, ] <- c(epoch, kl ) # 21.sep.20
      # original: kl <- kullback3(states=ff[[1]], frequencies=ff[[2]], Q=sampling[[1]], umbral=umbral ) #22-dic-17
      # original: diver[epoch, ] <- c(epoch, log2(kl) )
    }
    epoch <- epoch + 1
  }
  # OUTPUTS
  return(list(just_for_track, diver, cd, H, J, inferences, sampling[[1]], cd2 ))
}
# EJEMPLO:
#system.time(
#  machine_results <- bmachine(condInitial = "a5", muestra = wbasket, Kmax = 100, metrop_its = 20000, LearningRate = 0.8, decay = 0.02, N_sa = 10, umbral = .7, every = 10)
#)
# Recuperacion de la informacion
#just_for_track <- machine_results[[1]]
#diver <- machine_results[[2]]
#cd <- machine_results[[3]]
#H <- machine_results[[4]]
#J <- machine_results[[5]]
#inferences <- machine_results[[6]]
#last_mh_sampling <- machine_results[[7]] # equivale a sampling[[1]] del ultimo metropolis hasting efectuado.
# # # # # # # # # # # ## # # # BOLTZMANN-LEARNING # ## # # # # ## # # # # ## # # # # ##  # # # # ## 


