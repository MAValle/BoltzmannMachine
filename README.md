# BoltzmannMachine
Set of functions in R to carry out Boltzmann Learning.

Dependences:
* bmachine_function.R 

  -> initial_conditions_function.R
  
  -> learning_par_function.R
  
  -> contrastive_divergence_function.R
  
  -> kullback_leibler_function.R (kullback3 is used)
  
  -> metropolis_sampling_function.R
  
* kullback_leibler_function.R 

  -> frecuency_function.R

* metropolis_sampling_function.R

  -> tempfunc_function.R
  
  -> sweep_function.R
  
  -> goodness_function.R
  
* sweep_function.R

  -> activation_energy_function.R
  
  -> activiation_probability_function.R
  

# Example
The "data" is a dataframe with N columns and M rows. Every can be considered as a sample state of the 
system. Every column can be considered as a spin of the system with state "0" (off) or "1" (on).

Inputs:
* condInitial: see initial_conditions_function.R for options.
* muestra: input sample for training.
* Kmax: Number ib iteratiosn for the Boltzmann Machine.
* metrop_its: number of Metropolis-Hasting iterations.
* LearningRate: See learning_par_function.R 
* decay: decay for learning parameter (see learning_par_function.R) 
* N_sa: Number of samples for the negative phase.
* umbral: Threshold to measure KL divergence (see kullback_leibler_function.R)
* every: cada cuantas iteraciones medir el KL.


	```
{
  "firstName": "John",
  "lastName": "Smith",
  "age": 25
}
```


`machine_results <- bmachine(condInitial = "a2", 
                            muestra = data, 
                            Kmax = 25000, 
                            metrop_its = 1519, 
                            LearningRate = 0.95, 
                            decay = 0.0004, 
                            N_sa = 10,          
                            umbral = 0.8,                 
                            every = 20,                         
                            hg = NULL,                       
                            Jg = NULL)`
                            

Information retrieving:

`just_for_track <- machine_results[[1]]`
`diver <- machine_results[[2]]`
`cd <- machine_results[[3]]`  
`cd2 <- machine_results[[8]]` 
`H <- machine_results[[4]]`
`J <- machine_results[[5]]`


