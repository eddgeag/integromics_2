
source("./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R")
source("./scripts_R/scripts_utiles/scripts_funciones/calculo_medianas.R")
datos <- readRDS("../datos/preprocesado_08_09_23/novoom.rds")

set.seed(126581)

directorio <- "./scripts_R/IP/resultados_univariantes_totales"

if(!dir.exists(directorio)){
  dir.create(directorio)
}

X <- scale(as.matrix(as.data.frame(lapply(datos$totales$ip,as.numeric),
                                   row.names = rownames(datos$totales$ip))))

## search na

filas_y_variables_con_NA <- function(df) {
  # Obtener las filas con NA
  filas_con_NA <- rownames(df)[apply(df, 1, function(x) any(is.na(x)))]
  
  # Obtener las variables con NA
  variables_con_NA <- colnames(df)[apply(df, 2, function(x) any(is.na(x)))]
  
  # Crear una lista con las filas y variables con NA
  resultado <- list(Filas_con_NA = filas_con_NA, Variables_con_NA = variables_con_NA)
  
  return(resultado)
}

matriz_con_NA_a_cero <- function(df) {
  # Reemplazar NA por 0 en el marco de datos
  df_sin_NA <- ifelse(is.na(df), 0, df)
  
  return(df_sin_NA)
}
grupo <- datos$totales$general_data$GROUP
obesidad <- datos$totales$general_data$OBESE

resultado <- analyze_data(as.matrix(X),grupo,obesidad,correccion = 3)
write.csv(resultado$p_valores,file.path(directorio,"p_valores.csv"))
write.csv(resultado$interpretacion,file.path(directorio,"interpretacion.csv"))
write.csv(compute_medians(X,grupo,obesidad),file.path(directorio,"medianas.csv"))


directorio <- "./scripts_R/IP/resultados_univariantes_comunes"

if(!dir.exists(directorio)){
  dir.create(directorio)
}
X <- scale(as.matrix(datos$comunes$ip))
grupo <- datos$comunes$variables_in_bacteria$GROUP
obesidad <- datos$comunes$variables_in_bacteria$OBESE

resultado <- analyze_data(as.matrix(X),grupo,obesidad,correccion = 3)

write.csv(resultado$p_valores,file.path(directorio,"p_valores.csv"))
write.csv(resultado$interpretacion,file.path(directorio,"interpretacion.csv"))
write.csv(compute_medians(X,grupo,obesidad),file.path(directorio,"medianas.csv"))



