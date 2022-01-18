# IMPORTANTE: Setear directorio de trabajo!
setwd("~/Documentos/GitHub/Conservacion-de-codones-raros/dataframes")

# 1 - Función que une dataframes en uno solo
combinar_dataframes = function(dataframe_codones, dataframe_bicodones){ 
 
  df1 <- read.csv(dataframe_codones[1])
  df2 <- read.csv(dataframe_bicodones[1])
  data3 = data.frame(codon = character(length = 61), referencia = numeric(length = 61), conservacion = numeric(length = 61))
  data6 = data.frame(bicodon = character(length = 3721), referencia = numeric(length = 3721), conservacion = numeric(length = 3721))
  data3$codon <- df1$codon
  data6$bicodon <- df2$bicodon
  
  directorio <- getwd()
  
  
  for(archivo in dataframe_codones){
    dfc <- read.csv(paste(directorio, archivo, sep = "/", collapse = NULL))
    data3$referencia <- data3$referencia + dfc$referencia
    data3$conservacion <- data3$conservacion + dfc$conservacion
  }

  for(archivo in dataframe_bicodones){
    dfb <- read.csv(paste(directorio, archivo, sep = "/", collapse = NULL))
    data6$referencia <- data6$referencia + dfb$referencia
    data6$conservacion <- data6$referencia + dfb$conservacion
  }
  bicod <- subset(data6, referencia > 0)
  bicod2 <- subset(bicod, conservacion > 0)

  write.csv(data3,"/home/usuario/Documentos/GitHub/Conservacion-de-codones-raros/dataframes/data3_facundo.csv", row.names = FALSE)
  write.csv(data6,"/home/usuario/Documentos/GitHub/Conservacion-de-codones-raros/dataframes/data6_facundo.csv", row.names = FALSE)
  write.csv(bicod2,"/home/usuario/Documentos/GitHub/Conservacion-de-codones-raros/dataframes/data6_infoutil.csv", row.names = FALSE)
}

temp1 = list.files(pattern="historial_codones*")
temp2 = list.files(pattern="historial_bicodones*")

temp1

combinar_dataframes(temp1, temp2)

dataframe_de_prueba = data.frame(codon = character(), referencia = numeric(), conservacion = numeric())

