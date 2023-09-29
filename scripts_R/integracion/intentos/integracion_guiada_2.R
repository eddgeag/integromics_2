library(glmnet)
library(caret)
library(ALDEx2)
library(MOFA2)
library(psych)
library(vegan)
library(mixOmics)
library(caret)
library(glmnet)
source("./scripts_R/scripts_utiles/scripts_funciones/otras_funciones_utiles.R")
source(
  "./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R"
)
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
      rownames(x[x$importance > 0.05,]))))
  
  
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
      rownames(x[x$importance > 0.05,]))))
  
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
    ylab(paste("PC2", varianzas[2])) + facet_grid(~ obesidad)
  p3 <-
    ggplot(pca.plot, aes(PC1, PC2, color = sexo)) + geom_point() + xlab(paste("PC1", varianzas[1])) +
    ylab(paste("PC2", varianzas[2])) + facet_grid(~ obesidad)
  p4 <- ggarrange(p1, p2, p3)
  p1 <-
    ggplot(pca.plot, aes(PC1, PC2, color = obesidad)) + xlab(paste("PC1", varianzas[1])) +
    ylab(paste("PC2", varianzas[2])) + geom_point()
  p2 <-
    ggplot(pca.plot, aes(PC1, PC2, color = obesidad)) + xlab(paste("PC1", varianzas[1])) +
    ylab(paste("PC2", varianzas[2])) + geom_point() + facet_grid(~ grupo)
  p3 <-
    ggplot(pca.plot, aes(PC1, PC2, color = obesidad)) + xlab(paste("PC1", varianzas[1])) +
    ylab(paste("PC2", varianzas[2])) + geom_point() + facet_grid(~ sexo)
  
  
  p5 <- ggarrange(p1, p2, p3)
  p1 <-
    ggplot(pca.plot, aes(PC1, PC2, color = sexo)) + xlab(paste("PC1", varianzas[1])) +
    ylab(paste("PC2", varianzas[2])) + geom_point()
  p2 <-
    ggplot(pca.plot, aes(PC1, PC2, color = sexo)) + xlab(paste("PC1", varianzas[1])) +
    ylab(paste("PC2", varianzas[2])) + geom_point() + facet_grid(~ grupo)
  p3 <-
    ggplot(pca.plot, aes(PC1, PC2, color = sexo)) + xlab(paste("PC1", varianzas[1])) +
    ylab(paste("PC2", varianzas[2])) + geom_point() + facet_grid(~ sexo)
  
  
  p6 <- ggarrange(p1, p2, p3)
  
  return(list(p4, p5, p6))
}


extract_univ <- function(x, correc) {
  univ <- analyze_data(x, grupo, obesidad, correc)
  
  univ.media <- univ$mediana
  univ.p_valores <- univ$p_valores[,-c(1, 3)]
  
  colnames(univ.p_valores) <- colnames(univ.media)
  return(list(media = univ.media,
              pvalores = univ.p_valores))
}


datos  <- readRDS("../datos/preprocesado_05_02_23/novoom.rds")


metaboloma <- (transformacion(datos$comunes$metaboloma))
metagenoma <- scale(t(mean_aldex(datos$comunes$microbiota$genero)))
grupo <- datos$comunes$grupo
obesidad <- datos$comunes$obesidad
sexo <- datos$comunes$variables_in_bacteria$SEX


trainControl <- trainControl(method = "repeatedcv",
                             number = 5,
                             # prSummary needs calculated class probs
                             repeats = 10)
O <- (cbind(metaboloma, metagenoma))
set.seed(123456)
idx <- sample(nrow(O), size = nrow(O) * 0.7)
O.train <- O[idx,]
O.test <- O[-idx,]
obesidad.train <- obesidad[idx]
obesidad.test <- obesidad[-idx]
sexo.train <- sexo[idx]
sexo.test <- sexo[-idx]
set.seed(123456)
lasso.train <- cv.glmnet(
  x = O.train,
  y = obesidad.train,
  lambda.grid = seq(0, 1, length.out = 1000),
  type.measure = "class",
  nfolds = 10,
  family = "binomial",
  alpha = 0
)
prediccion <-
  as.factor(stats::predict(
    lasso.train,
    newx = O.test,
    s = "lambda.min",
    type = "class"
  ))

caret::confusionMatrix(prediccion, obesidad.test)$overall[1]
set.seed(123456)
final.model.obesidad <-
  glmnet(
    x = O,
    y = obesidad,
    family = "binomial",
    alpha = 0,
    lambda = lasso.train$lambda.min
  )

coeficientes <- coef(final.model.obesidad)

coeficientes@x

variables.obesidad <-
  rownames(coeficientes)[abs(coeficientes@x) > 0.02]

set.seed(123456)
lasso.train <- cv.glmnet(
  x = O.train,
  y = sexo.train,
  lambda.grid = seq(0, 0.5, length.out = 1000),
  type.measure = "class",
  nfolds = 10,
  family = "binomial"
)

prediccion <-
  as.factor(stats::predict(
    lasso.train,
    newx = O.test,
    s = "lambda.min",
    type = "class"
  ))

caret::confusionMatrix(prediccion, sexo.test)$overall[1]
set.seed(123456)
final.model.sexo <-
  glmnet(
    x = O,
    y = sexo,
    family = "binomial",
    alpha = 1,
    lambda = lasso.train$lambda.min
  )

coeficientes <- coef(final.model.sexo)


variables.sexo <- rownames(coeficientes)[abs(coeficientes@x) > 1]

variables.input <- unique(c(variables.obesidad, variables.sexo))




library(mixOmics) # import the mixOmics library

if (!file.exists("./scripts_R/integracion/diablos.rds")) {
  X.list <-
    list(metaboloma = as.matrix(metaboloma),
         metagenoma = as.matrix(metagenoma))
  Y <- grupo
  # names(Y) <- rownames(metaboloma)
  correlacion <-
    pls(X.list$metaboloma, X.list$metagenoma, ncomp = 1)
  correlacion <- cor(correlacion$variates$X, correlacion$variates$Y)
  
  (diseno <-
      t(matrix(
        c(0, correlacion, correlacion, 0),
        byrow = T,
        ncol = 2
      )))
  colnames(diseno) <- c("metaboloma", "metagenoma")
  rownames(diseno) <- colnames(diseno)
  diablo.grupo <-
    mixOmics::block.plsda(X.list,
                          as.numeric(obesidad),
                          design = diseno,
                          ncomp = 10) # run the method
  
  
  perf.diablo <- perf(
    diablo.grupo,
    validation = "Mfold",
    folds = 5,
    nrepeat = 10,
    # use repeated cross-validation
    progressBar = FALSE,
    auc = TRUE
  ) # include AUC values
  
  
  ncomp <-
    perf.diablo$choice.ncomp$AveragedPredict[2,] # what is the optimal value of components according to perf()
  
  list.keepX <-
    list(
      metaboloma = which(colnames(metaboloma) %in% variables.input),
      metagenoma = which(colnames(metagenoma) %in% variables.input)
    )
  
  
  
  
  tune.splsda.srbct <-
    tune.block.splsda(
      X.list,
      Y,
      ncomp = ncomp,
      # calculate for first 4 components
      validation = 'Mfold',
      folds = 5,
      nrepeat = 10,
      # use repeated cross-validation
      dist = 'max.dist',
      # use max.dist measure
      measure = "BER",
      # use balanced error rate of dist measure
      test.keepX = list.keepX,
      progressBar = T,
      design = diseno
    ) #'
  
  list.keepX <- tune.splsda.srbct$choice.keepX
  
  final.diablo.model = block.splsda(
    X = X.list,
    Y = Y,
    ncomp = ncomp,
    keepX = list.keepX,
    design = diseno
  )
  
} else{
  final.diablo.model <- readRDS("./scripts_R/integracion/diablos.rds")
  
}


obesidad.diablo.metaboloma <-
  final.diablo.model$obesidad$loadings$metaboloma
obesidad.diablo.metagenoma <-
  final.diablo.model$obesidad$loadings$metagenoma

rownames(obesidad.diablo.metaboloma)[apply(obesidad.diablo.metaboloma,
                                           1, function(x)
                                             any(x != 0))]

obesidad.vars <-
  rownames(obesidad.diablo.metagenoma)[apply(obesidad.diablo.metagenoma,
                                             1, function(x)
                                               any(x != 0))]

grupo.diablo.metaboloma <-
  final.diablo.model$grupo$loadings$metaboloma
sexo.diablo.metaboloma <-
  final.diablo.model$sexo$loadings$metaboloma
grupo.diablo.metagenoma <-
  final.diablo.model$grupo$loadings$metagenoma
sexo.diablo.metagenoma <-
  final.diablo.model$sexo$loadings$metagenoma

sexo.vars <-
  rownames(sexo.diablo.metagenoma)[apply(sexo.diablo.metagenoma,
                                         1, function(x)
                                           any(x != 0))]
grupo.vars <-
  rownames(grupo.diablo.metagenoma)[apply(grupo.diablo.metagenoma,
                                          1, function(x)
                                            any(x != 0))]

input_vars <-
  unique(c(sexo.vars, grupo.vars, obesidad.vars, colnames(metaboloma)))



ncomps <- 2:11
directorio_modelo <- "./scripts_R/integracion/mofa2_guided_2"
if (length(list.files(directorio_modelo)) == 0) {
  semillas <- seq(10, 1000, 10)
  for (s in semillas) {
    sapply(ncomps, function(x)
      mofa_componentes_original(
        ncomp = x,
        semilla = s,
        directorio_modelo = directorio_modelo,
        variables = input_vars
      ))
    
  }
} else{
  modelos <-
    lapply(
      list.files("./scripts_R/integracion/mofa2_guided_2/",
                 full.names = T),
      load_model
    )
  
}

# }
w <- which.min(compare_elbo(modelos, return_data = T)$ELBO)






modelo <- modelos[[w]]

pesos <- Reduce(rbind, get_weights(modelo, scale = F))

factores <- get_factors(modelo, scale = F)[[1]]


R <- factores %*% t(pesos)

funcion_pca(R, scle = T, trans = T)


p1 <- plot_variance_explained(modelo, max_r2 = 15)

p2 <- plot_variance_explained(modelo, plot_total = T)[[2]]

p3 <- ggarrange(p1, p2)

p3
dir.res <- "./scripts_R/integracion/resultados_mofa2_guided_2"
if (!dir.exists(dir.res)) {
  dir.create(dir.res)
}

ggsave(plot = p3,
       filename = file.path(dir.res, "explicacion_Varianza.jpeg"))

p4 <- ggcorrplot::ggcorrplot(corr = corr.test(factores)$r)

p4

ggsave(plot = p4,
       filename = file.path(dir.res, "sanity_check.jpeg"))


p5 <- correlate_factors_with_covariates(
  modelo,
  covariates = c(
    "BMI",
    "EDAD",
    "Free_ESTRA",
    "SHBG",
    "Total_ESTR",
    "WC",
    "WHR",
    "hsCRP"
  ),
  plot = "log_pval"
)

ggsave(plot = p5,
       filename = file.path(dir.res, "covaraites_factors.jpeg"))

p5

metagenoma_cor_Fac <-
  corr.test(metagenoma[input_vars, ], factores, method = "spearman")
p.values <- metagenoma_cor_Fac$p
metagenome_important <-
  rownames(metagenoma_cor_Fac$r)[which(apply(p.values, 1, function(x)
    any(x < 0.05)) == T)]

metagenome.cor <- metagenoma_cor_Fac$r[metagenome_important, ]
metagenome.pval <- metagenoma_cor_Fac$p[metagenome_important, ]

jpeg(filename = file.path(dir.res, "corr_metagenome_1_part.jpeg"))
corrplot::corrplot(metagenome.cor[1:20, ], p.mat = metagenome.pval[1:20, ], insig = "label_sig")
dev.off()


metaboloma_cor_Fac <-
  corr.test(metaboloma, factores, method = "spearman")
p.values <- metaboloma_cor_Fac$p
metaboloma_important <-
  rownames(metaboloma_cor_Fac$r)[which(apply(p.values, 1, function(x)
    any(x < 0.05)) == T)]
metabolome.cor <- metaboloma_cor_Fac$r[metaboloma_important, ]
metabolome.pval <- metaboloma_cor_Fac$p[metaboloma_important, ]
jpeg(filename = file.path(dir.res, "corr_metabolome.jpeg"))
corrplot::corrplot(metabolome.cor, p.mat = metabolome.pval, insig = "label_sig")
dev.off()


funcion_pca(R, scle = T, trans = F)

funcion_manova <- function(z) {
  pvals.grupo <-
    sapply(1:1000, function(x)
      summary(aovp(factores[, z] ~ grupo * obesidad, perm = "Exact"))[[1]]$`Pr(Prob)`[1])
  pval.grupo <- mean(pvals.grupo)
  
  pvals.grupo <-
    sapply(1:1000, function(x)
      summary(aovp(factores[, z] ~ grupo * obesidad, perm = "Exact"))[[1]]$`Pr(Prob)`[2])
  pval.obesidad <- mean(pvals.grupo)
  
  pvals.grupo <-
    sapply(1:1000, function(x)
      summary(aovp(factores[, z] ~ grupo * obesidad, perm = "Exact"))[[1]]$`Pr(Prob)`[3])
  pval.interaccion <- mean(pvals.grupo)
  
  return(c(grupo=pval.grupo,
           obesidad=pval.obesidad,
           pval.interaccion = pval.interaccion))
  
}

resultado.manova <- lapply(1:8, funcion_manova)





factores.df <- as.data.frame(factores)
factores.df$grupo <- grupo
factores.df$obesidad <- obesidad
factores.df$sexo <- sexo

ggplot(factores.df,aes(Factor3,Factor2,col=grupo))+geom_point()

R <- factores[,1:3] %*% t(pesos[,1:3])

funcion_pca(R,scle = F,trans = F)


pcx <- prcomp(t(pesos))
plote <- as.data.frame(pcx$x)
plote$label <- as.factor(paste0("PC",rep(1:8)))
ggplot(plote,aes(PC1,PC2,label=label))+geom_point()+geom_text()

csa <- factoextra::fviz_contrib(pcx,"var")
datx <- csa$data
cosa <- rownames(datx)[(100-cumsum(datx$contrib))>80]

R <- factores[,2:3] %*% t(pesos[cosa,2:3])

funcion_pca(R,scle = F,trans = F)[[1]]

#### ahora veremos grupo, obesidad, interaccion

Factores_analisis <- analyze_data(factores,grupo,obesidad,correccion = 1)

jpeg(filename = file.path(dir.res,"association_goi.jpeg"))
corrplot::corrplot(as.matrix(Factores_analisis$p_valores[,1:3]),
                   is.corr = F,
                   p.mat=as.matrix(Factores_analisis$p_valores[,1:3]),
                   insig = "label_sig")

dev.off()


Factores_analisis2 <- extract_univ(factores,correc = 1)

jpeg(filename = file.path(dir.res,"association_goi_2.jpeg"))
corrplot::corrplot(as.matrix(Factores_analisis2$media),
                   is.corr = F,
                   p.mat=as.matrix(Factores_analisis2$pvalores),
                   insig = "label_sig",method = "square")

dev.off()






# ##### Grupo
#
#
# Y <- grupo
# # names(Y) <- rownames(metaboloma)
# correlacion <- pls(X.list$metaboloma,X.list$metagenoma,ncomp = 1)
# correlacion <- cor(correlacion$variates$X,correlacion$variates$Y)
#
# (diseno <- t(matrix(c(0,correlacion,correlacion,0),byrow = T,ncol = 2)))
# colnames(diseno) <- c("metaboloma","metagenoma")
# rownames(diseno) <- colnames(diseno)
# diablo.grupo <- mixOmics::block.plsda(X.list, as.factor(grupo),design = diseno,ncomp = 10) # run the method
#
#
# perf.diablo <- perf(diablo.grupo, validation = "Mfold",
#                     folds = 5, nrepeat = 10, # use repeated cross-validation
#                     progressBar = FALSE, auc = TRUE) # include AUC values
#
#
# ncomp <- perf.diablo$choice.ncomp$AveragedPredict[2,] # what is the optimal value of components according to perf()
#
# pesos.totales <- unique(c(pesos_obesidad,pesos.sexo))
# list.keepX <- list(metaboloma=which(colnames(metaboloma) %in% pesos.totales),
#                    metagenoma=which(colnames(metagenoma) %in% pesos.totales))
#
#
#
#
# tune.splsda.srbct <- tune.block.splsda(X.list, Y, ncomp = ncomp, # calculate for first 4 components
#                                        validation = 'Mfold',
#                                        folds = 5, nrepeat = 10, # use repeated cross-validation
#                                        dist = 'max.dist', # use max.dist measure
#                                        measure = "BER", # use balanced error rate of dist measure
#                                        test.keepX = list.keepX,progressBar = T,
#                                        design = diseno) #
#
# list.keepX <- tune.splsda.srbct$choice.keepX
#
# final.diablo.model.grupo = block.splsda(X = X.list, Y = Y, ncomp = ncomp,
#                                         keepX = list.keepX, design = diseno)
#
#
#
# datos.plot <- as.data.frame((bind_cols(final.diablo.model.grupo$variates$metaboloma[,1:2],final.diablo.model.grupo$variates$metagenoma[,1:2],
#                                        grupo=grupo,obesidad=obesidad,sexo=sexo)))
#
# colnames(datos.plot) <-c("mPC1","mPC2","gPC1","gPC2","grupo","obesidad","sexo")
#
# ggplot(datos.plot,aes(mPC1,gPC1,color=obesidad))+geom_point()
#
# diablo.grupo.metaboloma <- final.diablo.model.grupo$loadings$metaboloma
# diablo.grupo.metagenoma <- final.diablo.model.grupo$loadings$metaboloma
#
# ##### sexo
#
#
# Y <- sexo
# # names(Y) <- rownames(metaboloma)
# correlacion <- pls(X.list$metaboloma,X.list$metagenoma,ncomp = 1)
# correlacion <- cor(correlacion$variates$X,correlacion$variates$Y)
#
# (diseno <- t(matrix(c(0,correlacion,correlacion,0),byrow = T,ncol = 2)))
# colnames(diseno) <- c("metaboloma","metagenoma")
# rownames(diseno) <- colnames(diseno)
# diablo.sexo <- mixOmics::block.plsda(X.list, as.factor(sexo),design = diseno,ncomp = 10) # run the method
#
#
# perf.diablo <- perf(diablo.sexo, validation = "Mfold",
#                     folds = 5, nrepeat = 10, # use repeated cross-validation
#                     progressBar = FALSE, auc = TRUE) # include AUC values
#
#
# ncomp <- perf.diablo$choice.ncomp$AveragedPredict[2,] # what is the optimal value of components according to perf()
#
# pesos.totales <- unique(c(pesos_obesidad,pesos.sexo))
# list.keepX <- list(metaboloma=which(colnames(metaboloma) %in% pesos.totales),
#                    metagenoma=which(colnames(metagenoma) %in% pesos.totales))
#
#
#
#
# tune.splsda.srbct <- tune.block.splsda(X.list, Y, ncomp = ncomp, # calculate for first 4 components
#                                        validation = 'Mfold',
#                                        folds = 5, nrepeat = 10, # use repeated cross-validation
#                                        dist = 'max.dist', # use max.dist measure
#                                        measure = "BER", # use balanced error rate of dist measure
#                                        test.keepX = list.keepX,progressBar = T,
#                                        design = diseno) #
#
# list.keepX <- tune.splsda.srbct$choice.keepX
#
# final.diablo.model.sexo = block.splsda(X = X.list, Y = Y, ncomp = ncomp,
#                                        keepX = list.keepX, design = diseno)
#
#
#
# datos.plot <- as.data.frame((bind_cols(final.diablo.model.grupo$variates$metaboloma[,1:2],final.diablo.model.grupo$variates$metagenoma[,1:2],
#                                        grupo=grupo,obesidad=obesidad,sexo=sexo)))
#
# colnames(datos.plot) <-c("mPC1","mPC2","gPC1","gPC2","grupo","obesidad","sexo")
#
# ggplot(datos.plot,aes(mPC1,gPC2,color=grupo))+geom_point()+facet_grid(~obesidad)
#
# diablo.sexo.metaboloma <- final.diablo.model.sexo$loadings$metaboloma
# diablo.sexo.metagenoma <- final.diablo.model.sexo$loadings$metaboloma
#
# diablos <- list(obesidad=final.diablo.model,
#                 grupo=final.diablo.model.grupo,
#                 sexo=final.diablo.model.sexo)
#
# saveRDS(diablos,"./scripts_R/integracion/diablos.rds")
