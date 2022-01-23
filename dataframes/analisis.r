setwd("/home/usuario/Documentos/GitHub/Conservacion-de-codones-raros/dataframes/")
library(ggplot2)
# Importar datasets
datos_codones <- read.csv("data3_facundo.csv")
datos_bicodones <- read.csv("data6_facundo.csv")

# Relación conservación-referencia
datos_codones$tasaConservacion <- datos_codones$conservacion / datos_codones$referencia

# Lista provisoria, consultar luego si se incluyen los codones subóptimos
codones_raros <- c("TTT", "TTA", "CTT", "CTC", "CTA", "ATT", "ATA", "GTT", "GTC", "GTA", "TCT", 
                   "TCA", "CCT", "ACT", "ACA", "ACG", "GCT", "GCA", "GCG", "TAT", "CAT", "CAA",
                   "AAT", "AAA", "GAA", "TGT", "CGT", "CGA", "CGG", "GCT", "GGG")
vector_logico <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
                   TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, 
                   TRUE, TRUE, TRUE, TRUE, TRUE)
dflogic <- data.frame(codones_raros, vector_logico) # Solución momentánea, buscar una forma más elegante

# Establecer codones raros
# Raro: TRUE, no raro: FALSE
datos_codones$esRaro <- logical(length = 61)

for(i in codones_raros){
  for(j in datos_codones$codon){
    if(i==j){
      datos_codones[datos_codones$codon==j, 5] <- TRUE
    }
  }
}

# Barplot
codones_barplot <- ggplot(datos_codones,
                          aes(   x = codon,
                                 y = tasaConservacion,
                              fill = esRaro)) +
                   geom_bar(stat = "identity")

codones_barplot

# Ahora para bicodones:
datos_bicodones$tasaConservacion <- datos_bicodones$conservacion / datos_bicodones$referencia
datos_bicodones$esRaro <- c(length = 3721) # Default: FALSE
datos_bicodones$primerCodon <- character(length = 3721)
datos_bicodones$segundoCodon <- character(length = 3721)

# Separar los pares de codones en codones individuales, para luego calcular
# el producto de las probabilidades individuales de cada codón y hacer una 
# regresión lineal
iterador = 1
for(i in datos_bicodones$bicodon){
  vector_codones <- strsplit(i, "(?<=.{3})", perl = TRUE)[[1]]
  datos_bicodones[iterador, 6] <- vector_codones[1]
  datos_bicodones[iterador, 7] <- vector_codones[2]
  iterador = iterador + 1
}


# Tasas individuales de cada codon
datos_bicodones$tasaCodon1 <- numeric(length = 3721)
datos_bicodones$tasaCodon2 <- numeric(length = 3721)

for(i in datos_codones$codon){
  for(j in datos_bicodones$primerCodon){
  if(i == j){
    datos_bicodones[datos_bicodones$primerCodon == i, 8] <- datos_codones[datos_codones$codon == i, 4]
    }
  }
}

for(j in datos_codones$codon){
  for(k in datos_bicodones$segundoCodon){
    if(j == k){
      datos_bicodones[datos_bicodones$segundoCodon == j, 9] <- datos_codones[datos_codones$codon == j, 4]    
    }  
  }
}

# Cálculo del producto de tasas
datos_bicodones$productoDeTasas <- datos_bicodones$tasaCodon1 * datos_bicodones$tasaCodon2

# Determinación de pares de codones inhibitorios
datos_bicodones$cod1Bool <- logical(length = 3721)
datos_bicodones$cod2Bool <- logical(length = 3721)

for(i in codones_raros){
  for(j in datos_bicodones$primerCodon){
    if(i == j){
      datos_bicodones[datos_bicodones$primerCodon == i, 11] <- dflogic[1, 2]
    }
  }
}

for(i in codones_raros){
  for(k in datos_bicodones$segundoCodon){
    if(i == k){
      datos_bicodones[datos_bicodones$segundoCodon == i, 12] <- dflogic[1, 2]
    }
  }
}

datos_bicodones$esRaro <- datos_bicodones$cod1Bool * datos_bicodones$cod2Bool # Raro: 1, no raro: 0














