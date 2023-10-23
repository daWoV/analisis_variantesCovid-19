library(seqinr)
library(stringr)
library(ggplot2)
getwd()

setwd("C:/Users/daeld/OneDrive/Documentos/Programas")

alpha <- read.fasta("sequence_alpha.fasta")
beta <- read.fasta("sequence_beta.fasta")
delta <- read.fasta("sequence_delta.fasta")
gamma <- read.fasta("sequence_gamma.fasta")

alphalen <- length(alpha[[1]])
betalen <- length(beta[[1]])
deltalen <- length(delta[[1]])
gammalen <- length(gamma[[1]])

print("Longitud de las secuencias de las variantes")
print(paste0("Variante Alpha: ", alphalen))
print(paste0("Variante Beta: ", betalen))
print(paste0("Variante Delta: ", deltalen))
print(paste0("Variante Gamma: ", gammalen))

alphadf <- as.data.frame(count(alpha[[1]], 1))
colnames(alphadf) <- c("Nucleotidos", "Alpha")

betadf <- as.data.frame(count(beta[[1]], 1))
colnames(betadf) <- c("Nucleotidos", "Beta")

deltadf <- as.data.frame(count(delta[[1]], 1))
colnames(deltadf) <- c("Nucleotidos", "Delta")

gammadf <- as.data.frame(count(gamma[[1]], 1))
colnames(gammadf) <- c("Nucleotidos", "Gamma")

df1 <- merge(alphadf, betadf)
df2 <- merge(df1, deltadf)
df3 <- merge(df2, gammadf)

cantidad1 <- alphadf[,2]
cantidad2 <- betadf[,2]
cantidad3 <- deltadf[,2]
cantidad4 <- gammadf[,2]

grafica <- data.frame(
  "Nucleótidos" = c("A", "C", "G", "T"),
  "Variantes" = c(rep("Alpha", 4), rep("Beta", 4), rep("Delta", 4), rep("Gamma", 4)),
  "Cantidad" = c(cantidad1, cantidad2, cantidad3, cantidad4)
)

graficar <- print(ggplot(grafica, aes(x = Nucleótidos, y = Cantidad, fill=Variantes)) + geom_bar(position="dodge", stat="identity"))
graficar

#En la gráfica de barras se puede observar la cantidad de nucleótidos que tiene cada una de la variantes que se investigaron del SARS-CoV-2. 
#La variante Alpha está representada con el color rojo, la variante Beta con un color verde, la variante Delta con un color acua y finalmente de color morado la variante Gamma.
#En el eje "x" de la grafica se ubican las bases nitrogenadas (A, C, G, T) y en el eje "y" la cantidad de las mismas. 

porcentajeGC <- function(seq) {
  secuencia <- count(seq[[1]], 1)
  C <- secuencia[[2]]
  G <- secuencia[[3]]
  CG <- C + G
  total <- length(seq[[1]])

  porcentaje <- function(total, CG){
    return(round(100*CG/total, 2))
  }
  print(paste0("GC: ", porcentaje(total, CG), "%"))
}
print("Alpha:")
porcentajeGC(alpha)
print("Beta:")
porcentajeGC(beta)
print("Delta:")
porcentajeGC(delta)
print("Gamma:")
porcentajeGC(gamma)

contrasentido <- function(seq) {
  s <- comp(seq[[1]])
  s <- str_to_upper(paste(s, collapse = ""))
  return(s)
}
print("Alpha:")
print(contrasentido(alpha))
print("Beta:")
print(contrasentido(beta))
print("Delta:")
print(contrasentido(delta))
print("Gamma:")
print(contrasentido(gamma))