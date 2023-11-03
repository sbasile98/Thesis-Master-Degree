#Confronto le medie ottenute dall' algoritmo evolutivo con quelle ottenute da Hybrid-IA.
#Tale confronto è rappresentanto da un valore percentuale indicante di quanto le medie ottenute dall'algoritmo si scostino ripetto a quelle di Hybrid-IA

library(dplyr)

setwd("C:/Users/simon/Downloads/Tesi/results")

results <- read.csv("results_cot_1.csv", sep = ";")
compare <- read.csv("FeedbackVertexSet.csv", sep = "\t", header = FALSE)

final = results[2:5]
final = distinct(final)

for (i in 1:nrow(results)) {
  filename = strsplit(results$filename[i], "/")[[1]][3]
  tmp = filter(compare, V1 == filename)
  results$HybridIA[i] <- tmp$V2[1]
}

for (i in 1:nrow(final)) {
  tmp = filter(results, n == final$n[i] & e == final$e[i] & l == final$l[i] & u == final$u[i])
  
  filename = strsplit(tmp$filename[1], "/")[[1]][3]
  final$filename[i] <- filename
  
  avg_1 = sum(tmp$best) / 5
  final$EA[i] <- avg_1
  
  avg_2 = sum(tmp$HybridIA) / 5
  final$HybridIA[i] <- avg_2
  
  final$diff[i] <- (avg_1 - avg_2)/avg_2 * 100
}

final <- final[c(5,1,2,3,4,6,7,8)]

grid = filter(final, grepl("Grid", final$filename))
rand = filter(final, grepl("Rand", final$filename))

write.csv(rand,"rand.csv", row.names = FALSE)
write.csv(grid,"grid.csv", row.names = FALSE)