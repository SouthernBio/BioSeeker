# IMPORTANTE: Setear directorio de trabajo!
setwd("~/Documentos/GitHub/Conservacion-de-codones-raros/dataframes")

# 1 - Función que une dataframes en uno solo
combinar_dataframes = function(dataframe_codones, dataframe_bicodones){
  codones <- as.data.frame(c())
  bicodones <- as.data.frame(c())
  referencia <- as.data.frame(c())
  conservacion <- as.data.frame(c())
  
  df1 <- read.csv(dataframe_codones[1])
  df2 <- read.csv(dataframe_bicodones[1])
  data3 = data.frame(codones, referencia, conservacion)
  data6 = data.frame(bicodones, referencia, conservacion)
  data3$codones <- data3$codones + df1$codones
  data6$bicodones <- data6$bicodones + df2$bicodones

  for(archivo in dataframe_codones){
    read.csv(archivo)
    attach(data3)
    referencia <- referencia + archivo$referencia
    conservacion <- conservacion + archivo$conservacion
    detach(data3)
  }

  for(archivo in dataframe_bicodones){
    read.csv(archivo)
    attach(data6)
    referencia <- referencia + archivo$referencia
    conservacion <- referencia + archivo$conservacion
    detach(data6)
  }
  bicod <- subset(data6, referencia > 0)
  bicod2 <- subset(bicod, conservacion > 0)

  write.csv(data3,"/home/usuario/Escritorio/data3_facundo.csv", row.names = TRUE)
  write.csv(data6,"/home/usuario/Escritorio/data6_facundo.csv", row.names = TRUE)
  write.csv(bicod2,"/home/usuario/Escritorio/data6_infoutil.csv", row.names = TRUE)
}

temp1 = list.files(pattern="historial_codones*")
temp2 = list.files(pattern="historial_bicodones*")

combinar_dataframes(temp1, temp2)
