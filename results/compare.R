library(dplyr)

#setwd("desktop/CN/code/results")
setwd("C:/Users/simon/OneDrive/Desktop/Tesi/results")
#setwd("C:/Users/simon/Downloads/Tesi/results")

#results <- read.csv("test_2.csv", sep = ";")
#results <- read.csv("test_3.csv", sep = ";") 
#results <- read.csv("test_4.csv", sep = ";") #Miglioramento del codice di test_3: aggiunta della local search durante il processo evolutivo
#results <- read.csv("test_5.csv", sep = ";") #Stesso codice di test_4
#results <- read.csv("test_6.csv", sep = ";") #Miglioramento greedy della creazione della popolazione
#results <- read.csv("test_ml_2.csv", sep = ";") #Miglioramento ML della probabilità di mutazione
#results <- read.csv("test_7.csv", sep = ";") #Stesso codice di test_6
#results <- read.csv("results_cot_1.csv", sep = ";") #Test solo con crossover type 1
#results <- read.csv("results_cot_2.csv", sep = ";") #Test solo con crossover type 2
#results <- read.csv("results_cot_vario.csv", sep = ";") #Test con crossover type dinamico
#results <- read.csv("res_ml_004.csv", sep = ";")
results <- read.csv("res_ml_finale.csv", sep = ";")


compare <- read.csv("FeedbackVertexSet.csv", sep = "\t", header = FALSE)

optimum = 0
approx = 0
perc = 0

for (i in 1:nrow(results)) {
  filename = strsplit(results$filename[i], "/")[[1]][3]
  tmp = filter(compare, V1 == filename)
  results$HybridIA[i] <- tmp$V2[1]
  results$diff[i] <- results$best[i] - results$HybridIA[i]
  
  if (results$diff[i] == 0) optimum = optimum + 1
  
  if (results$diff[i] <= 10) approx = approx + 1
  
  results$diffPerc[i] <- results$diff[i] * 100 /results$HybridIA[i]
  if (results$diffPerc[i] <= 5) perc = perc + 1
}

optimumRatio = optimum/nrow(results)*100
approxRatio = approx/nrow(results)*100
percRatio = perc/nrow(results)*100