library(limma)
library(ALDEx2)
library(MOFA2)
source(
  "./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R"
)
source("./scripts_R/scripts_utiles/scripts_funciones/manova_vanvalen.R")


mofa_componentes_original <-
  function(ncomp,
           semilla,
           directorio_modelo,
           mofa.obj,
           sample_metadata_) {
    samples_metadata(mofa.obj) <- sample_metadata_
    
    if (!dir.exists(directorio_modelo)) {
      dir.create(directorio_modelo)
    }
    data_opts <- get_default_data_options(mofa.obj)
    data_opts$scale_views <- F
    
    
    model_opts <- get_default_model_options(mofa.obj)
    
    model_opts$num_factors <- ncomp
    
    train_opts <- get_default_training_options(mofa.obj)
    train_opts$seed <- semilla
    train_opts$convergence_mode <- "slow"
    
    ### probar estocastico
    
    MOFAobject <- prepare_mofa(
      object = mofa.obj,
      data_options = data_opts,
      model_options = model_opts,
      training_options = train_opts
    )
    
    
    outfile <-
      file.path(directorio_modelo, paste0(semilla, "_", ncomp, "model.hdf5"))
    
    MOFAobject.trained <-
      run_mofa(MOFAobject, outfile, use_basilisk = T)
    return(MOFAobject.trained)
    
  }
source("./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R")
source("./scripts_R/scripts_utiles/scripts_funciones/manova_vanvalen.R")


datos <- readRDS("../datos/preprocesado_08_09_23/novoom-04-10-23.rds")
grupo <- datos$comunes$grupo
obesidad <- datos$comunes$obesidad
sexo <- datos$comunes$variables_in_bacteria$SEX
## metaboloma
x <- datos$comunes$metaboloma[-c(38:40),]
succinato <- colMeans(datos$comunes$metaboloma)
x <- rbind(x,succinato)
rownames(x)[38] <- "Succ"
rownames(x)<- gsub("_mean",replacement = "",x=rownames(x))
x.tmp <- x[-38]
x.tmp <- add_classificacion_metaboloma(x.tmp[-38,])
rownames(x)[-38] <- paste(x.tmp$grupo_meta,rownames(x.tmp),sep="_")

x.t <- scale(t(x))

## microbioma

mat <-  scale(t(mean_aldex(datos$comunes$microbiota$genero)))

## ip

ip <- scale(datos$comunes$ip)
ip.d <- t(ip[,-27])

## clinicas/covariables
clinicas <- datos$comunes$variables_in_bacteria
clinicas.mean <- t((clinicas[, grep("mean", colnames(clinicas))]))
clinicas.other <- as.data.frame(bind_cols(sample=datos$comunes$variables_in_bacteria$Paciente,
                                          (clinicas[, c(6:10, 12, 14)])))
hscrp <- ((datos$comunes$variables_in_bacteria$hsCRP))
clinicas.other$hsCRP <- hscrp
colnames(clinicas.mean) <- datos$comunes$variables_in_bacteria$Paciente


covariables <- clinicas.other



colnames(x.t)[grep("CARBOHIDRATOS_GRASAS_KETONA_GLYCEROL_",colnames(x.t))]<-gsub("CARBOHIDRATOS_GRASAS_KETONA_GLYCEROL_","CGKG_",colnames(x.t)[grep("CARBOHIDRATOS_GRASAS_KETONA_GLYCEROL_",colnames(x.t))])
datos.list <- list(metaboloma=t(x.t),
                   metagenoma=t(mat),
                   ip = ip.d)


lista.datos <- create_mofa_from_matrix(
  datos.list
)

directorio.general <- "./re_analisis_2/modelos_final"

if (!dir.exists(directorio.general)) {
  dir.create(directorio.general)
}


sample_metadata <- (cbind(t(clinicas.mean),covariables))
rownames(sample_metadata) <- colnames(datos.list$metaboloma)
for (ncomp in 4:12) {
  for (semilla in sample(40000, size = 10)) {
    obj <-  mofa_componentes_original(ncomp,
                                      semilla,
                                      directorio.general,
                                      lista.datos,
                                      sample_metadata_ =sample_metadata)
  }
}

#
#
#
#
#


