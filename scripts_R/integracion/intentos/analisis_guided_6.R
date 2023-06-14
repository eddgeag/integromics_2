library(glmnet)
library(caret)
library(ALDEx2)
library(MOFA2)
library(psych)
library(vegan)
library(mixOmics)
library(caret)
library(glmnet)


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
           variables) {
    vars.metaboloma <- variables[grep("mean", variables)]
    vars.metagenoma <- variables[-grep("mean", variables)][-1]
    
    metaboloma <-
      (transformacion(datos$comunes$metaboloma))[, vars.metaboloma]
    metagenoma <-
      scale(t(mean_aldex(datos$comunes$microbiota$genero)))[, vars.metagenoma]
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


metaboloma <- scale(t(datos$comunes$metaboloma))
metagenoma <- scale(t(mean_aldex(datos$comunes$microbiota$genero)))
grupo <- datos$comunes$grupo
obesidad <- datos$comunes$obesidad
sexo <- datos$comunes$variables_in_bacteria$SEX

int <- sample(nrow(metaboloma),nrow(metaboloma)*0.7)
X <- cbind(metaboloma,metagenoma)
X.train <- X[int,]
grupo.train <- sexo[int]
X.test <- X[-int,]
grupo.test <- sexo[-int]

set.seed(849
         )
# simulate data

# fit multinomial
ridge_fit <-  cv.glmnet(X.train,
                        grupo.train,
                        family = "binomial",
                        alpha = 0,
                        lambda = seq(0,1,0.0001),nfolds = 10)
best_lambda_ridge <- ridge_fit$lambda.min
ridge_bestfit <- ridge_fit$glmnet.fit
ridge_pred <- as.factor(stats::predict(ridge_bestfit, s = best_lambda_ridge, 
                      newx = X.test,type="class"))

confusionMatrix(ridge_pred,grupo.test)

variables <- (colnames(X))[coef(ridge_fit)@x[-1]>0.01]


metagenoma <- mean_aldex(aldex.clr(datos$comunes$microbiota$genero.abs[,colnames(metagenoma) %in% variables]))
metaboloma <- metaboloma[,colnames(metaboloma) %in% variables]



modelos <- lapply(1:11,function(x) mofa_componentes_original(ncomp = x,
                                                             semilla = 1234,
                                                             directorio_modelo = "./scripts_R/integracion/prueba_n",
                                                             variables = variables))


w <- which.min(unlist(lapply(modelos,get_elbo)))

saveRDS(modelos[[w]],"./modelo_prueba.rds")





