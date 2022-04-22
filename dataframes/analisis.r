setwd("/home/usuario/Documents/GitHub/Conservacion-de-codones-raros/dataframes/")
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
                   geom_bar(stat = "identity") +
                   labs(    x = "Codón", 
                            y = "Tasa de conservación",
                        title = "Conservación de codones individuales") +
                   guides(fill=guide_legend(title="Codón raro:"))

codones_barplot

# Ahora para bicodones:
datos_bicodones$tasaConservacion <- datos_bicodones$conservacion / datos_bicodones$referencia
datos_bicodones$esRaro <- c(length = 3721) 
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
      datos_bicodones[datos_bicodones$primerCodon == j, 11] <- TRUE
    }
  }
}

for(i in codones_raros){
  for(k in datos_bicodones$segundoCodon){
    if(i == k){
      datos_bicodones[datos_bicodones$segundoCodon == k, 12] <- TRUE
    }
  }
}

datos_bicodones$esRaro <- datos_bicodones$cod1Bool * datos_bicodones$cod2Bool # Raro: 1, no raro: 0

modelo <- ggplot(data = datos_bicodones,
                 aes(   x = productoDeTasas,
                        y = tasaConservacion,
                        color = esRaro)) +
          geom_point(stat = "identity")

modelo

bicodones_raros <- subset(datos_bicodones, esRaro == 1)

plot_pares_inhibitorios <- ggplot(data = bicodones_raros,
                                  aes(   x = productoDeTasas,
                                         y = tasaConservacion,
                                     label = bicodon)) +
                          geom_rect(aes(xmin = 0, xmax = 0.05, 
                                        ymin = 0, ymax = 0.05)) +
                           geom_point(color = "blueviolet") +
                           xlim(0, 0.05) +
                           ylim(0, 0.2) +
                           geom_hline(yintercept = 0.05, linetype = "dotdash") +
                           geom_vline(xintercept = 0.05, linetype = "dotdash") +
                           geom_abline(aes(intercept = 0,
                                               slope = 1)) +
                           labs(x = "Producto de tasas de conservación de codones constituyentes",
                                y = "Tasa de conservación real del bicodón raro",
                                title = "Tasa real vs. tasa predicha para bicodones")  +
                           geom_text(aes(label = bicodon), hjust = 0, vjust = 0)
                           
plot_pares_inhibitorios

# Seleccionar aquellos codones cuya tasa de conservación supera la predicha
datos_bicodones$tasaRealSobrePredicha <- datos_bicodones$tasaConservacion / datos_bicodones$productoDeTasas
bicodones_raros <- subset(datos_bicodones, esRaro == 1)
bicodones_raros_conservados <- subset(bicodones_raros, tasaRealSobrePredicha > 25)

plot_pares_inhibitorios_conservados <- ggplot(data = bicodones_raros_conservados,
                                              aes(   x = productoDeTasas,
                                                     y = tasaConservacion,
                                                     label = bicodon)) +
                                       geom_point(color = "blueviolet") +
                                       xlim(0, 0.004) +
                                       ylim(0, 0.2) +
                                       labs(x = "Producto de tasas de conservación de codones constituyentes",
                                            y = "Tasa de conservación real del bicodón raro",
                                            title = "Pares de codones conservados")  +
                                       geom_text(aes(label = bicodon), hjust = 0, vjust = 0)

plot_pares_inhibitorios_conservados

# Ver qué sucede con los bicodones normales
bicodones_normales <- subset(datos_bicodones, esRaro == 0)
bicodones_normales_conservados <- subset(bicodones_normales, tasaRealSobrePredicha > 50)

# Ajuste lineal
modelo_lineal <- lm(tasaConservacion ~ 0 + productoDeTasas, data = datos_bicodones)

summary(modelo_lineal)
#plot(modelo_lineal) # No hay normalidad. Pero ajusta con un r cuadrado de 0.81


# Z-score distribution plots
# Necesito agrupar pares de codones de acuerdo al dipéptido que codifican
# Además necesito agregar pseudocounts a los pares de codones cuya conservación fue cero

codon   <- c("TTT", "TTC", "TTA", "TTG", "CTT", "CTC", "CTA", "CTG", "ATT", "ATC", "ATA", "ATG", "GTT", "GTC", "GTA", "GTG",
             "TCT", "TCC", "TCA", "TCG", "CCT", "CCC", "CCA", "CCG", "ACT", "ACC", "ACA", "ACG", "GCT", "GCC", "GCA", "GCG",
             "TAT", "TAC", "CAT", "CAC", "CAA", "CAG", "AAT", "AAC", "AAA", "AAG", "GAT", "GAC", "GAA", "GAG",
             "TGT", "TGC", "TGG", "CGT", "CGC", "CGA", "CGG", "AGT", "AGC", "AGA", "AGG", "GGT", "GGC", "GGA", "GGG")

peptido <- c("Phe", "Phe", "Leu", "Leu", "Leu", "Leu", "Leu", "Leu", "Ile", "Ile", "Ile", "Met", "Val", "Val", "Val", "Val",
             "Ser", "Ser", "Ser", "Ser", "Pro", "Pro", "Pro", "Pro", "Thr", "Thr", "Thr", "Thr", "Ala", "Ala", "Ala", "Ala",
             "Tyr", "Tyr", "His", "His", "Gln", "Gln", "Asn", "Asn", "Lys", "Lys", "Asp", "Asp", "Glu", "Glu",
             "Cys", "Cys", "Trp", "Arg", "Arg", "Arg", "Arg", "Ser", "Ser", "Arg", "Arg", "Gly", "Gly", "Gly", "Gly")

peptido_segun_codon <- data.frame(codon, peptido, row.names = NULL)
datos_bicodones$primerPeptido <- character(length = 3721)
datos_bicodones$segundoPeptido <- character(length = 3721)

for(i in peptido_segun_codon$codon){
  for(j in datos_bicodones$primerCodon){
    if(i == j){
      datos_bicodones[datos_bicodones$primerCodon == j, 14] <- peptido_segun_codon[peptido_segun_codon$codon == j, 2]
    }
  }
}

for(i in peptido_segun_codon$codon){
  for(k in datos_bicodones$segundoCodon){
    if(i == k){
      datos_bicodones[datos_bicodones$segundoCodon == k, 15] <- peptido_segun_codon[peptido_segun_codon$codon == k, 2]
    }
  }
}

# Algo va mal. Los loops for están insertando el índice del peptido en lugar del péptido



datos_bicodones$dipeptido <- paste(datos_bicodones$primerPeptido, datos_bicodones$segundoPeptido, sep = "-")
datos_bicodones$factorDipeptido <- factor(datos_bicodones$dipeptido)

# Para el Z-score necesito calcular la media de la tasa real sobre predicha
# en cada familia de dipeptidos, ademas del desvio estandar de la familia

datos_bicodones$logTasaNormalizada <- log(datos_bicodones$tasaRealSobrePredicha)

tapply(datos_bicodones$logTasaNormalizada, datos_bicodones$factorDipeptido, mean)

# Grabar dataframes
#write.csv(datos_codones, file = "/home/usuario/Documentos/GitHub/Conservacion-de-codones-raros/dataframes/analisis_codones.csv", row.names = FALSE)
#write.csv(datos_bicodones, file = "/home/usuario/Documentos/GitHub/Conservacion-de-codones-raros/dataframes/analisis_bicodones.csv", row.names = FALSE)











