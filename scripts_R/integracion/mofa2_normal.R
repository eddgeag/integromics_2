library(glmnet)
library(caret)
library(ALDEx2)
library(MOFA2)
library(psych)
library(vegan)
library(mixOmics)
library(caret)
library(glmnet)
source("./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R")
source("./scripts_R/scripts_utiles/scripts_funciones/manova_vanvalen.R")

transformacion <- function(X) {
  return((t(log2(X / colSums(
    X
  )))))
}

mofa_componentes <-
  function(ncomp,
           semilla,
           directorio_modelo,
           variables) {
    idx <- sample(46)
    vars.metaboloma <- variables[grep("mean", variables)]
    vars.metagenoma <- variables[-grep("mean", variables)][-1]
    
    metaboloma <-
      (transformacion(datos$comunes$metaboloma))[idx, vars.metaboloma]
    metagenoma <-
      scale(t(mean_aldex(datos$comunes$microbiota$genero)))[idx, vars.metagenoma]
    mofa.obj <- create_mofa_from_matrix(list(
      metaboloma = t(metaboloma),
      metagenoma = t(metagenoma)
    ))
    samples_metadata(mofa.obj) <-
      data.frame(
        sample = datos$comunes$variables_in_bacteria$Paciente,
        datos$comunes$variables_in_bacteria[, 3:16]
      )
    rownames(samples_metadata(mofa.obj)) <-
      datos$comunes$variables_in_bacteria$Paciente
    
    if (!dir.exists(directorio_modelo)) {
      dir.create(directorio_modelo)
    }
    data_opts <- get_default_data_options(mofa.obj)
    data_opts$scale_views <- T
    
    model_opts <- get_default_model_options(mofa.obj)
    
    model_opts$num_factors <- ncomp
    
    train_opts <- get_default_training_options(mofa.obj)
    train_opts$seed <- semilla
    
    ### probar estocastico
    
    MOFAobject <- prepare_mofa(
      object = mofa.obj,
      data_options = data_opts,
      model_options = model_opts,
      training_options = train_opts
    )
    
    
    outfile <-
      file.path(directorio_modelo, paste0(semilla, "_", ncomp, "model.hdf5"))
    
    MOFAobject.trained <- run_mofa(MOFAobject, outfile)
    
  }

mofa_componentes_original <-
  function(ncomp,
           semilla,
           directorio_modelo,
           metaboloma,
           metagenoma) {
    mofa.obj <- create_mofa_from_matrix(list(
      metaboloma = t(metaboloma),
      metagenoma = t(metagenoma)
    ))
    samples_metadata(mofa.obj) <-
      data.frame(
        sample = datos$comunes$variables_in_bacteria$Paciente,
        datos$comunes$variables_in_bacteria[, 3:16]
      )
    rownames(samples_metadata(mofa.obj)) <-
      datos$comunes$variables_in_bacteria$Paciente
    
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
    
    MOFAobject.trained <- run_mofa(MOFAobject, outfile)
    return(MOFAobject.trained)
    
  }

extraer_variables <- function(modelo.diablo) {
  comps <- 1:ncomps
  
  maxima_contribucion.meta <- lapply(1:ncomps,
                                     function(x)
                                       plotLoadings(
                                         modelo.diablo,
                                         contrib = "max",
                                         method = "median",
                                         block = 1,
                                         comp = x,
                                         plot = F
                                       ))
  
  
  maxima_contribucion.meta.vars <-
    unique(unlist(lapply(maxima_contribucion.meta, function(x)
      rownames(x[x$importance > 0.05, ]))))
  
  
  maxima_contribucion.mic <- lapply(1:ncomps,
                                    function(x)
                                      plotLoadings(
                                        modelo.diablo,
                                        contrib = "max",
                                        method = "median",
                                        block = 2,
                                        comp = x,
                                        plot = F
                                      ))
  
  
  maxima_contribucion.mic.vars <-
    unique(unlist(lapply(maxima_contribucion.mic, function(x)
      rownames(x[x$importance > 0.05, ]))))
  
  variables <-
    unique(c(
      maxima_contribucion.meta.vars,
      maxima_contribucion.mic.vars
    ))
  return(variables)
}

funcion_pca <- function(X, scle, trans) {
  if (!trans) {
    pca <- prcomp(X, scale. = scle)
    
    varianzas <- round(100 * pca$sdev ^ 2 / sum(pca$sdev ^ 2), 2)
    pca.plot <- data.frame(
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 2],
      obesidad = obesidad,
      grupo = grupo,
      sexo = sexo
    )
  } else{
    pca <- prcomp(t(X), scale. = scle)
    
    varianzas <- round(100 * pca$sdev ^ 2 / sum(pca$sdev ^ 2), 2)
    pca.plot <- data.frame(
      PC1 = pca$rotation[, 1],
      PC2 = pca$rotation[, 2],
      obesidad = obesidad,
      grupo = grupo,
      sexo = sexo
    )
  }
  
  p1 <-
    ggplot(pca.plot, aes(PC1, PC2, color = grupo)) + geom_point() + xlab(paste("PC1", varianzas[1])) +
    ylab(paste("PC2", varianzas[2]))
  p2 <-
    ggplot(pca.plot, aes(PC1, PC2, color = grupo)) + geom_point() + xlab(paste("PC1", varianzas[1])) +
    ylab(paste("PC2", varianzas[2])) + facet_grid( ~ obesidad)
  p3 <-
    ggplot(pca.plot, aes(PC1, PC2, color = sexo)) + geom_point() + xlab(paste("PC1", varianzas[1])) +
    ylab(paste("PC2", varianzas[2])) + facet_grid( ~ obesidad)
  p4 <- ggarrange(p1, p2, p3)
  p1 <-
    ggplot(pca.plot, aes(PC1, PC2, color = obesidad)) + xlab(paste("PC1", varianzas[1])) +
    ylab(paste("PC2", varianzas[2])) + geom_point()
  p2 <-
    ggplot(pca.plot, aes(PC1, PC2, color = obesidad)) + xlab(paste("PC1", varianzas[1])) +
    ylab(paste("PC2", varianzas[2])) + geom_point() + facet_grid( ~ grupo)
  p3 <-
    ggplot(pca.plot, aes(PC1, PC2, color = obesidad)) + xlab(paste("PC1", varianzas[1])) +
    ylab(paste("PC2", varianzas[2])) + geom_point() + facet_grid( ~ sexo)
  
  
  p5 <- ggarrange(p1, p2, p3)
  p1 <-
    ggplot(pca.plot, aes(PC1, PC2, color = sexo)) + xlab(paste("PC1", varianzas[1])) +
    ylab(paste("PC2", varianzas[2])) + geom_point()
  p2 <-
    ggplot(pca.plot, aes(PC1, PC2, color = sexo)) + xlab(paste("PC1", varianzas[1])) +
    ylab(paste("PC2", varianzas[2])) + geom_point() + facet_grid( ~ grupo)
  p3 <-
    ggplot(pca.plot, aes(PC1, PC2, color = sexo)) + xlab(paste("PC1", varianzas[1])) +
    ylab(paste("PC2", varianzas[2])) + geom_point() + facet_grid( ~ sexo)
  
  
  p6 <- ggarrange(p1, p2, p3)
  
  return(list(p4, p5, p6))
}


extract_univ <- function(x, correc) {
  univ <- analyze_data(x, grupo, obesidad, correc)
  
  univ.media <- univ$mediana
  univ.p_valores <- univ$p_valores[, -c(1, 3)]
  
  colnames(univ.p_valores) <- colnames(univ.media)
  return(list(media = univ.media,
              pvalores = univ.p_valores))
}


datos  <- readRDS("../datos/preprocesado_05_02_23/novoom.rds")


metaboloma <- voom((datos$comunes$metaboloma))$E
metagenoma <- scale((mean_aldex(datos$comunes$microbiota$genero)))

lista.datos <- list(voom=list(metaboloma=metaboloma,
                              metagenoma=metagenoma))
# lista.datos <-
#   list(
#     scale.scale = list(metaboloma = scale(t(
#       datos$comunes$metaboloma
#     )),
#     metagenoma = scale(t(
#       mean_aldex(datos$comunes$microbiota$genero)
#     ))),
#     trans.scale = list(metaboloma = transformacion(t(
#       datos$comunes$metaboloma
#     )),
#     metagenoma = scale(t(
#       mean_aldex(datos$comunes$microbiota$genero)
#     ))),
#     trans.trans = list(metaboloma = transformacion(t(
#       datos$comunes$metaboloma
#     )),
#     metagenoma = (t(
#       mean_aldex(datos$comunes$microbiota$genero)
#     ))),
#     scale.trans = list(metaboloma = scale(t(
#       datos$comunes$metaboloma
#     )),
#     metagenoma = (t(
#       mean_aldex(datos$comunes$microbiota$genero)
#     )))
#   )

get_pesos <- function(sclae=F,mdl) Reduce(rbind, get_weights(mdl, scale = sclae))
get_factores <- function(sclae=F,mdl) get_factors(mdl, scale = sclae)[[1]]


grupo <- datos$comunes$grupo
obesidad <- datos$comunes$obesidad
sexo <- datos$comunes$variables_in_bacteria$SEX
lista.res <- vector("list",length=length(lista.datos))
for(l in 1:length(lista.datos)){
  
  metaboloma <- lista.datos[[l]]$metaboloma
  metagenoma <- lista.datos[[l]]$metagenoma
  
  if(ncol(metagenoma)==46){
    metagenoma <-t(metagenoma)
  }
  if(ncol(metaboloma)==46){
    metaboloma <- t(metaboloma)
    
  }

  dir.tmp <- paste0("./scripts_R/integracion/tmp_mofa_voom","_iteracion_",l)
  
  modelos <- lapply(2:11, function(x) mofa_componentes_original(x,
                                                                123,
                                                                directorio_modelo = dir.tmp,
                                                                metaboloma,
                                                                metagenoma))
  elbos <- unlist(lapply(modelos,get_elbo))
  w <- which.max(elbos)

  modelo <- modelos[[w]]
  
  factores <- get_factores(F,modelo)
  pesos <- get_pesos(F,modelo)
  R <- factores %*% t(pesos)
  
  funcion_pca(R,trans = F,scle = F)[[1]]
  
  van.ob <- vanValen.test(R,obesidad)
  van.grupo <- vanValen.test(R,grupo)
  van.interaccion <- vanValen.test(R,interaction(grupo,obesidad))
  van.females <- vanValen.test(R[grupo!="Male",],as.factor(as.character(grupo[grupo!="Male"])))
  van.herma <- vanValen.test(R[grupo!="Female",],as.factor(as.character(grupo[grupo!="Female"])))
  van.control <- vanValen.test(R[grupo!="PCOS",],as.factor(as.character(grupo[grupo!="PCOS"])))
  
  van.ob.males <- vanValen.test(R[grupo=="Male",],as.factor(as.character(obesidad[grupo=="Male"])))
  van.ob.females <- vanValen.test(R[grupo=="Female",],as.factor(as.character(obesidad[grupo=="Female"])))
  van.ob.pcos <- vanValen.test(R[grupo=="PCOS",],as.factor(as.character(obesidad[grupo=="PCOS"])))
  
  
  van.noob.males <- vanValen.test(R[grupo!="Male" & obesidad=="No Obese",],as.factor(as.character(grupo[grupo!="Male" & obesidad=="No Obese"])))
  
  van.noob.females <-  vanValen.test(R[grupo!="Female" & obesidad=="No Obese",],as.factor(as.character(grupo[grupo!="Female" & obesidad=="No Obese"])))
  
  van.noob.pcos <-vanValen.test(R[grupo!="PCOS" & obesidad=="No Obese",],as.factor(as.character(grupo[grupo!="PCOS" & obesidad=="No Obese"])))
  
  
  vanvalen.tests <- c(obesidad=van.ob,
                      grupo=van.grupo,
                      interaccion=van.interaccion,
                      van.females=van.females,
                      van.herma=van.herma,
                      van.control = van.control,
                      obesidad.males = van.ob.males,
                      obesidad.females=van.ob.females,
                      obesidad.pcos=van.ob.pcos,
                      noob.males=van.noob.males,
                      noob.females=van.noob.females,
                      noob.pcos=van.noob.pcos)
  
  lista.res[[l]] <- vanvalen.tests
}


csa <- Reduce("rbind",lista.res)

funcion_pca(R,scle = F,trans = F)
