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

directorio_modelo <- "./scripts_R/integracion/MOFA2_guided"

if (!dir.exists(directorio_modelo)) {
  dir.create(directorio_modelo)
}



semillas <- 1000
ncomps <- 2:11
matriz.elbo <- matrix(NA, nrow = 20, ncol = length(ncomps))
matriz.dev <- matrix(NA, nrow = 20, ncol = length(ncomps))
matriz.nas <- matrix(NA, nrow = 20, ncol = length(ncomps))
matriz.comps <- matrix(NA, nrow = 20, ncol = length(ncomps))

for (n in 1:length(ncomps)) {
  nombres.metaboloma <- 1:37
  nombres.metagenoma <- 1:213
  for (i in 1:20) {
    dir.tmp <-
      paste0("./scripts_R/integracion/vueltas/vuelta_",
             i,
             "ncomp_",
             n,
             "_")

    modelo <- mofa_componentes_original(
      ncomp = ncomps[n],
      semilla = semillas,
      directorio_modelo = dir.tmp,
      variables = c(colnames(metaboloma)[nombres.metaboloma],
                    colnames(metagenoma)[nombres.metagenoma])
    )


    elbo <-   get_elbo(modelo)

    factores <- get_factors(modelo, as.data.frame = F)[[1]]
    factores.df <- as.matrix(factores)
    pesos <- get_weights(modelo)
    pesos <- rbind(pesos$metaboloma, pesos$metagenoma)
    R <- factores.df %*% t(pesos)
    R <- as.data.frame(bind_cols(R, variable = grupo))
    mdl <- nnet::multinom(variable ~ ., data = R)
    mdl0 <- nnet::multinom(variable ~ 1, data = R)
    ratio.dev.interaccion <- 1 - deviance(mdl) / deviance(mdl0)



    pesos_metaboloma <-
      get_weights(modelo, "metaboloma", scale = T, abs = T)[[1]]
    pesos_metagenoma <-
      get_weights(modelo, "metagenoma", scale = T, abs = T)[[1]]

    nombres.metaboloma <-
      unique(unlist(apply(pesos_metaboloma, 2, function(x)
        which(x > 0.05))))
    nombres.metagenoma <-
      unique(unlist(apply(pesos_metagenoma, 2, function(x)
        which(x > 0.05))))

    matriz.elbo[i, n] <- elbo
    matriz.dev[i, n] <- ratio.dev.interaccion

    p.values <- corr.test(get_factors(modelo, scale = T)[[1]])$p
    matriz.nas[i, n] <- any(p.values[lower.tri(p.values)] < 0.05)
    matriz.comps[i, n] <-
      length(nombres.metaboloma) + length(nombres.metagenoma)

  }
}

csa <- matriz.dev

csa[matriz.nas == T] <- NA

indices <- which(csa==max(csa,na.rm = T),arr.ind=T)
# max(matriz.dev)
# #
# #
modelo <-
  load_model(paste0("./scripts_R/integracion/vueltas/vuelta_",
                    indices[1],
                    "ncomp_",indices[2],"_/1000_",indices[2]+1,"model.hdf5"))


pesos <- Reduce(rbind, get_weights(modelo, scale = F))

nombres.metagenoma <- rownames(pesos)[-grep("mean",rownames(pesos))]
nombres.metaboloma <- rownames(pesos)[grep("mean",rownames(pesos))]
library(xlsx)
metabolome <- read.xlsx("../datos/Integromics_Metabolome.xlsx",sheetIndex = 1)
variables_in_bacteria <- datos$comunes$variables_in_bacteria
sujetos <- variables_in_bacteria$Paciente
grupo <- variables_in_bacteria$GROUP
obesidad <- variables_in_bacteria$OBESE
metabolome.basal.tmp <-
  metabolome[metabolome$Order %in% variables_in_bacteria$Orden, ]
redundant_data_meta <- c("Order", "Sample", "GROUP", "OBESE")
variables_meta <- colnames(metabolome.basal.tmp)
mean_meta <- variables_meta[grep("_[GLP]0$", variables_meta)]
auc_meta <- variables_meta[grep("AU[CG]", variables_meta)]
w <-
  c(
    which(variables_meta %in% redundant_data_meta),
    which(variables_meta %in% mean_meta),
    which(variables_meta %in% auc_meta)
  )
metabolome.basal_mean <- metabolome.basal.tmp[,-w]
rownames(metabolome.basal_mean) <- sujetos

metaboloma.mean <- metabolome.basal_mean

metaboloma.red <- metaboloma.mean[,colnames(metaboloma.mean) %in% nombres.metaboloma]

metaboloma.trans <- scale(metaboloma.red)

metagenoma.raw <- read.xlsx("../datos/integromics_microbiota.xlsx",sheetIndex = 3)
bacteria_phylum.abs <- read.xlsx("../datos/integromics_microbiota.xlsx",sheetIndex = 1)
general_data <- datos$totales$general_data
proces_bacteria <- function(datos, tipo) {
  Order <- bacteria_phylum.abs$Order
  variables_in_bacteria <-
    general_data[general_data$Orden %in% Order , ]
  Sample <- variables_in_bacteria$Paciente
  
  
  if (tipo == "genera") {
    Order <- Order
    
    X <- datos[-1,-1]
    X <- apply(X, 2, as.numeric)
    colnames(X) <- datos[1, 2:ncol(datos)]
    rownames(X) <- Sample
    X.rel <- 100 * X / rowSums(X)
    
    
  } else if (tipo == "phylum") {
    X <- datos[,-c(1:2)]
    colnames(X) <- colnames(datos)[3:ncol(datos)]
    rownames(X) <- Sample
    X.rel <- 100 * X / rowSums(X)
  }
  
  ###
  GROUP <- variables_in_bacteria$GROUP
  SEX <- variables_in_bacteria$SEX
  OBESE <- variables_in_bacteria$OBESE
  
  X <- as.data.frame(X)
  X.rel <- as.data.frame(X.rel)
  
  X$GROUP <- GROUP
  X$SEX <- SEX
  X$OBESE <- OBESE
  
  
  X.rel$GROUP <- GROUP
  X.rel$SEX <- SEX
  X.rel$OBESE <- OBESE
  
  return(list(abs = X, rel = X.rel))
  
}


metagenoma.t <- proces_bacteria(t(metagenoma.raw),"genera")

metagenoma.abs <- metagenoma.t$abs

metagenoma.abs <- metagenoma.abs[,colnames(metagenoma.abs) %in% nombres.metagenoma]

metagenoma.aldex <- aldex.clr(t(metagenoma.abs),conds=grupo,mc.samples = 1000,verbose = T)

metagenoma.a <- as.matrix(mean_aldex(metagenoma.aldex))

funcion_pca(t(metagenoma.a),scle = F,trans = T)[[2]]

metagenoma <- t(metagenoma.a)
metaboloma <- metaboloma.trans

directorio.modelo <- "./scripts_R/integracion/guided"

semillas <- 2:11

modelos <- sapply(semillas,FUN = function(x) mofa_componentes_original(ncomp=x,
                                                            123456,
                                                            directorio_modelo,
                                                            variables=c(colnames(metaboloma),colnames(metagenoma))))


elbos <- unlist(lapply(modelos, get_elbo))

w <- which.max(elbos)

modelo <- modelos[[w]]



pesos <- Reduce(rbind, get_weights(modelo, scale = F))

factores <- get_factors(modelo, scale = F)[[1]]
#
#
R <- factores %*% t(pesos)
#
funcion_pca(R, scle = F, trans = T)





# #
p1 <- plot_variance_explained(modelo, max_r2 = 15)

p2 <- plot_variance_explained(modelo, plot_total = T)[[2]]

p3 <- ggarrange(p1, p2)

p3
dir.res <- "./scripts_R/integracion/resultados_mofa2_guided_4"
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
  corr.test(metagenoma, factores, method = "spearman")
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

resultado.manova <- lapply(1:ncol(factores), funcion_manova)





factores.df <- as.data.frame(factores)
factores.df$grupo <- grupo
factores.df$obesidad <- obesidad
factores.df$sexo <- sexo

ggplot(factores.df,aes(Factor2,Factor1,col=grupo))+geom_point()

R <- factores[,1:2] %*% t(pesos[,1:2])

funcion_pca(R,scle = F,trans = F)


pcx <- prcomp(t(pesos))
plote <- as.data.frame(pcx$x)
plote$label <- as.factor(paste0("PC",rep(1:3)))
ggplot(plote,aes(PC1,PC2,label=label))+geom_point()+geom_text()

csa <- factoextra::fviz_contrib(pcx,"var")
datx <- csa$data
cosa <- rownames(datx)[(100-cumsum(datx$contrib))>80]

R <- factores[,1:3] %*% t(pesos[cosa,1:3])

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

pesos <- Reduce(rbind, get_weights(modelo, scale = F))

factores <- get_factors(modelo, scale = F)[[1]]
R <- factores[,1:2] %*% t(pesos[,1:2])

funcion_pca(R,scle = F,trans = F)


pesos_metaboloma <-
  get_weights(modelo, "metaboloma", scale = T, abs = T)[[1]]
pesos_metagenoma <-
  get_weights(modelo, "metagenoma", scale = T, abs = T)[[1]]

nombres.metaboloma <- names(which(pesos_metaboloma[,1]>0.05))
nombres.metagenoma <- names(which(pesos_metagenoma[,2]>0.05))

R <- factores[,1:2] %*% t(pesos[c(nombres.metaboloma,
                                  nombres.metagenoma),1:2])

funcion_pca(R,scle = F,trans = T)

saveRDS(modelo,"./scripts_R/integracion/MOFA2.rds")


