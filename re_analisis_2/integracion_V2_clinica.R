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
    data_opts$scale_views <- T
    
    
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

datos <-
  readRDS("../datos/preprocesado_08_09_23/novoom-04-10-23.rds")
grupo <- datos$comunes$grupo
obesidad <- datos$comunes$obesidad
sexo <- datos$comunes$variables_in_bacteria$SEX

x <- datos$comunes$metaboloma[-c(38:40), ]
succinato <- colMeans(datos$comunes$metaboloma)
x <- rbind(x, succinato)
rownames(x)[38] <- "Succ"
rownames(x) <- gsub("_mean", replacement = "", x = rownames(x))
x.tmp <- x[-38]
x.tmp <- add_classificacion_metaboloma(x.tmp[-38, ])
rownames(x)[-38] <- paste(x.tmp$grupo_meta, rownames(x.tmp), sep = "_")

tot <- as.data.frame(t(x))
grupo <- datos$comunes$grupo
obesidad <- datos$comunes$obesidad

tot$grupo <- grupo
tot$obesidad <- obesidad
tot.melt <- reshape2::melt(tot)

ggplot(tot.melt, aes(value, color = obesidad)) + geom_density()


ip <- datos$comunes$ip

ip.d <- as.data.frame(scale(ip))
ip.d$grupo <- grupo
ip.d$obesidad <- obesidad

ip.melt <- reshape2::melt(ip.d)
ggplot(ip.melt, aes(value, color = obesidad)) + geom_density()

ip <- as.data.frame(t(scale(ip)))

clinicas <- datos$comunes$variables_in_bacteria
clinicas.mean <- clinicas[, grep("mean", colnames(clinicas))]
clinicas.other <- clinicas[, c(6:10, 12, 14)]

clinicas.tot <-
  as.data.frame(scale(bind_cols(clinicas.mean, clinicas.other)))
rownames(clinicas.tot) <- colnames(ip)


metaboloma.t <- t(scale(t(x)))
metagenoma <- scale((mean_aldex(datos$comunes$microbiota$genero)))
clinicas <- clinicas.tot

bloque_clinico <- t(as.matrix(clinicas[, grep("mean", colnames(clinicas))]))

colnames(bloque_clinico) <- rownames(clinicas)
rownames(bloque_clinico) <- colnames(clinicas[, grep("mean", colnames(clinicas))])

bloque_clinico <- as.data.frame(bloque_clinico)

covariables_clinicas <- clinicas[, -grep("mean", colnames(clinicas))]

samples_metadata_ <- data.frame(sample=rownames(covariables_clinicas),
                                covariables_clinicas)

lista.datos <- create_mofa_from_matrix(
  list(
    metaboloma = as.matrix(metaboloma.t),
    metagenoma = as.matrix(metagenoma),
    ip = as.matrix(ip),
    clinicas = as.matrix(bloque_clinico)
  )
)


directorio.general <- "./re_analisis_2/modelos_clinicos"
if (!dir.exists(directorio.general)) {
  dir.create(directorio.general)
}


for (ncomp in 2:7) {
  for (semilla in sample(40000, size = 10)) {
    obj <-  mofa_componentes_original(ncomp,
                                      semilla,
                                      directorio.general,
                                      lista.datos,
                                      sample_metadata_ =samples_metadata_)
  }
}


modelos <-
  lapply(list.files("./re_analisis_2/modelos_clinicos", full.names = T),
         MOFA2::load_model)

csa <- select_model(modelos, plot = T)
