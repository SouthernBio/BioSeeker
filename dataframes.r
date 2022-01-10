combinar_dataframes = function(dataframe_codones, dataframe_bicodones){
  df1 <- read.csv(dataframe_codones[1])
  df2 <- read.csv(dataframe_bicodones[1])
  data3 = data.frame(codones, referencia, conservacion)
  data6 = data.frame(bicodones, referencia, conservacion)
  data3$codones <- data3$codones + df1$codones
  data6$bicodones <- data6$bicodones + df2$bicodones

  for(archivo in dataframe_codones){
    data3$referencia <- data3$referencia + archivo$referencia
    data3$conservacion <- data3$conservacion + archivo$conservacion
  }

  for(archivo in dataframe_bicodones){
    data6$referencia <- data6$referencia + archivo$referencia
    data6$conservacion <- data6$referencia + archivo$conservacion
  }
  bicod <- subset(data6, referencia > 0)
  bicod2 <- subset(bicod, conservacion > 0)

  write.csv(data3,"/home/usuario/Escritorio/data3_facundo.csv", row.names = TRUE)
  write.csv(data6,"/home/usuario/Escritorio/data6_facundo.csv", row.names = TRUE)
  write.csv(bicod2,"/home/usuario/Escritorio/data6_infoutil.csv", row.names = TRUE)
}


