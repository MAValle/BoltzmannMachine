

# # # # # # # # # # # ## # # # FUNCION PARA GENERAR TEMPERATURA # # # # # # # # # # # ## # # # # #
# INPUT: valor de alpha udualmente entre 0.8 y 0.9999                                            #
#        temperatura_antes                                                                       #
#        ciclo de iteracion actual                                                               #
# OUTPUT: una temperatura para el ciclo k.                                                       #
tempfunc <- function(temperatura_antes, alpha) {
  temp <- temperatura_antes*alpha
  return(temp)
}
# # # # # # # # # # # ## # # # FUNCION PARA GENERAR TEMPERATURA # # # # # # # # # # # ## # # # # #

