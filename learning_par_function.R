
# viene de PAPER_MBA by solving the inverse ISING problem (ok)
# # # # # # # # # # # ## # # # LEARNING PARAMETER # ## # # # # ## # # # # ## # # # # ##  # # # # ## 
# FUNCION PARA GENERAR PARAMETRO DE APRENDIZAJE EN LA MAQUINA DE BOLTZMANN
# INPUT: valor de alpha usualmente entre 0.8 y 0.9999
#        temperatura_antes
#        ciclo de iteracion actual
# OUTPUT: una prametro para el ciclo m
learning_par <- function(ciclo, decay, LearningPar) {
  LearningPar <- LearningPar * (1/(1 + decay*ciclo))
  return(LearningPar)
}
# ejemplo
#LearningRate <- learning_par(ciclo = 3, decay = 0.0008, LearningPar = 0.001)
# # # # # # # # # # # ## # # # LEARNING PARAMETER # ## # # # # ## # # # # ## # # # # ##  # # # # ## 

