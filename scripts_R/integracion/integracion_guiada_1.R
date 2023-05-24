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
source("./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R")
transformacion <- function(X) {
  return((t(log2(X / colSums(
    X
  )))))
}

mofa_componentes <- function(ncomp,semilla,directorio_modelo,variables){
  idx <- sample(46)
  vars.metaboloma <- variables[grep("mean",variables)]
  vars.metagenoma <- variables[-grep("mean",variables)][-1]
  
  metaboloma <- (transformacion(datos$comunes$metaboloma))[idx,vars.metaboloma]
  metagenoma <- scale(t(mean_aldex(datos$comunes$microbiota$genero)))[idx,vars.metagenoma]
  mofa.obj <- create_mofa_from_matrix(list(
    metaboloma = t(metaboloma),
    metagenoma = t(metagenoma)
  ))
  samples_metadata(mofa.obj)<- data.frame(sample=datos$comunes$variables_in_bacteria$Paciente,
                                          datos$comunes$variables_in_bacteria[,3:16])
  rownames(samples_metadata(mofa.obj))<-datos$comunes$variables_in_bacteria$Paciente
  
  if (!dir.exists(directorio_modelo)) {
    dir.create(directorio_modelo)
  }
  data_opts <- get_default_data_options(mofa.obj)
  data_opts$scale_views <- T
  
  model_opts <- get_default_model_options(mofa.obj)
  
  model_opts$num_factors <- ncomp
  
  train_opts <- get_default_training_options(mofa.obj)
  train_opts$seed <- 789456
  
  ### probar estocastico
  
  MOFAobject <- prepare_mofa(
    object = mofa.obj,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  
  
  outfile <- file.path(directorio_modelo, paste0(semilla,"_",ncomp,"model.hdf5"))
  
  MOFAobject.trained <- run_mofa(MOFAobject, outfile)
  
}

extraer_variables <- function(modelo.diablo){
  
  comps <- 1:ncomps
  
  maxima_contribucion.meta <- lapply(1:ncomps, 
                                     function(x) plotLoadings(modelo.diablo,contrib = "max",method = "median",block = 1,comp = x,plot = F))
  
  
  maxima_contribucion.meta.vars <- unique(unlist(lapply(maxima_contribucion.meta, function(x) rownames(x[x$importance > 0.05,]))))
  
  
  maxima_contribucion.mic <- lapply(1:ncomps, 
                                    function(x) plotLoadings(modelo.diablo,contrib = "max",method = "median",block = 2,comp = x,plot = F))
  
  
  maxima_contribucion.mic.vars <- unique(unlist(lapply(maxima_contribucion.mic, function(x) rownames(x[x$importance>0.05,]))))
  
  variables <- unique(c(maxima_contribucion.meta.vars,maxima_contribucion.mic.vars))
  return(variables)
}

funcion_pca <- function(X,scle,trans){
  
  
  if(!trans){
    pca <- prcomp(X,scale. = scle)
    
    varianzas <- round(100*pca$sdev^2/sum(pca$sdev^2),2)
    pca.plot <- data.frame(PC1=pca$x[,1],
                           PC2=pca$x[,2],
                           obesidad=obesidad,
                           grupo = grupo,
                           sexo=sexo)
  }else{
    pca <- prcomp(t(X),scale. = scle)
    
    varianzas <- round(100*pca$sdev^2/sum(pca$sdev^2),2)
    pca.plot <- data.frame(PC1=pca$rotation[,1],
                           PC2=pca$rotation[,2],
                           obesidad=obesidad,
                           grupo = grupo,
                           sexo=sexo)
  }
  
  p1 <- ggplot(pca.plot,aes(PC1,PC2,color=grupo))+geom_point()+xlab(paste("PC1",varianzas[1]))+ylab(paste("PC2",varianzas[2]))
  p2 <- ggplot(pca.plot,aes(PC1,PC2,color=grupo))+geom_point()+xlab(paste("PC1",varianzas[1]))+ylab(paste("PC2",varianzas[2]))+facet_grid(~obesidad)
  p3 <- ggplot(pca.plot,aes(PC1,PC2,color=sexo))+geom_point()+xlab(paste("PC1",varianzas[1]))+ylab(paste("PC2",varianzas[2]))+facet_grid(~obesidad)
  p4 <- ggarrange(p1,p2,p3)
  p1 <- ggplot(pca.plot,aes(PC1,PC2,color=obesidad))+xlab(paste("PC1",varianzas[1]))+ylab(paste("PC2",varianzas[2]))+geom_point()
  p2 <- ggplot(pca.plot,aes(PC1,PC2,color=obesidad))+xlab(paste("PC1",varianzas[1]))+ylab(paste("PC2",varianzas[2]))+geom_point()+facet_grid(~grupo)
  p3 <- ggplot(pca.plot,aes(PC1,PC2,color=obesidad))+xlab(paste("PC1",varianzas[1]))+ylab(paste("PC2",varianzas[2]))+geom_point()+facet_grid(~sexo)
  
  
  p5 <- ggarrange(p1,p2,p3)
  p1 <- ggplot(pca.plot,aes(PC1,PC2,color=sexo))+xlab(paste("PC1",varianzas[1]))+ylab(paste("PC2",varianzas[2]))+geom_point()
  p2 <- ggplot(pca.plot,aes(PC1,PC2,color=sexo))+xlab(paste("PC1",varianzas[1]))+ylab(paste("PC2",varianzas[2]))+geom_point()+facet_grid(~grupo)
  p3 <- ggplot(pca.plot,aes(PC1,PC2,color=sexo))+xlab(paste("PC1",varianzas[1]))+ylab(paste("PC2",varianzas[2]))+geom_point()+facet_grid(~sexo)
  
  
  p6<- ggarrange(p1,p2,p3)
  
  return(list(p4,p5,p6))
}


extract_univ <- function(x,correc){
  
  univ<- analyze_data(x,grupo,obesidad,correc)
  
  univ.media <- univ$mediana
  univ.p_valores <- univ$p_valores[,-c(1,3)]
  
  colnames(univ.p_valores) <- colnames(univ.media)
  return(list(media=univ.media,
              pvalores= univ.p_valores))
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
                             repeats=10)
O <- (cbind(metaboloma,metagenoma))
set.seed(123456)
idx <- sample(nrow(O),size = nrow(O)*0.7)
O.train <- O[idx,]
O.test <- O[-idx,]
obesidad.train <- obesidad[idx]
obesidad.test <- obesidad[-idx]
sexo.train <- sexo[idx]
sexo.test <- sexo[-idx]
set.seed(123456)
lasso.train <- cv.glmnet(x=O.train,
                   y=obesidad.train,
                   lambda.grid = seq(0, 1, length.out = 1000),
                   type.measure = "class",
                   nfolds = 10,
                   family="binomial",
                   alpha=0)
prediccion <- as.factor(stats::predict(lasso.train,newx=O.test,s = "lambda.min",type="class"))

caret::confusionMatrix(prediccion,obesidad.test)$overall[1]
set.seed(123456)
final.model.obesidad <- glmnet(x=O,y=obesidad,family = "binomial",alpha = 0,lambda=lasso.train$lambda.min)

coeficientes <- coef(final.model.obesidad)

coeficientes@x

variables.obesidad <- rownames(coeficientes)[abs(coeficientes@x)>0.02]

set.seed(123456)
lasso.train <- cv.glmnet(x=O.train,
                         y=sexo.train,
                         lambda.grid = seq(0, 0.5, length.out = 1000),
                         type.measure = "class",
                         nfolds = 10,
                         family="binomial")

prediccion <- as.factor(stats::predict(lasso.train,newx=O.test,s = "lambda.min",type="class"))

caret::confusionMatrix(prediccion,sexo.test)$overall[1]
set.seed(123456)
final.model.sexo <- glmnet(x=O,y=sexo,family = "binomial",alpha = 1,lambda=lasso.train$lambda.min)

coeficientes <- coef(final.model.sexo)


variables.sexo <- rownames(coeficientes)[abs(coeficientes@x)>1]

variables.input <- unique(c(variables.obesidad,variables.sexo))

semillas <- seq(40,80,10)

directorio_modelo <- "./scripts_R/integracion/modelo_optimal_comps_guided"
if (!dir.exists(directorio_modelo)) {
  dir.create(directorio_modelo)
}

ncomps <- 2:11
# if(length(list.files(directorio_modelo))==0){
  for(s in semillas){
    
    lapply(ncomps, function(x) mofa_componentes(ncomp = x,semilla = s,directorio_modelo = directorio_modelo,variables=variables.input ))
    
  }
# }


modelos <- lapply(list.files(path = directorio_modelo,
                             full.names = T), load_model)

w <- which.min(compare_elbo(modelos,return_data = T)$ELBO)

modelo <- modelos[[w]]

pesos <- Reduce(rbind,get_weights(modelo,scale = T))

factores <- get_factors(modelo,scale = T)[[1]]

p1 <- plot_variance_explained(modelo, max_r2=15)

p2 <- plot_variance_explained(modelo, plot_total = T)[[2]]

p3 <- ggarrange(p1,p2)

p3
dir.res <- "./scripts_R/integracion/resultados_mofa2_guided"
if(!dir.exists(dir.res)){
  dir.create(dir.res)
}

ggsave(plot=p3,filename = file.path(dir.res,"explicacion_Varianza.jpeg"))

p4 <- ggcorrplot::ggcorrplot(corr = corr.test(factores)$r)

p4

ggsave(plot=p4,filename = file.path(dir.res,"sanity_check.jpeg"))


p5 <- correlate_factors_with_covariates(modelo, 
                                        covariates = c("BMI","EDAD","FREE_TEST",
                                                       "Free_ESTRA","RatioFreeTE2",
                                                       "SHBG","TOTAL_TEST","Total_ESTR",
                                                       "WC","WHR","hsCRP"), 
                                        plot="log_pval"
)

ggsave(plot=p5,filename = file.path(dir.res,"covaraites_factors.jpeg"))

p5

metagenoma_cor_Fac <- corr.test(metagenoma,factores,method = "spearman")
p.values <- metagenoma_cor_Fac$p
metagenome_important <- rownames(metagenoma_cor_Fac$r)[which(apply(p.values,1,function(x) any(x<0.05))==T)]

metagenome.cor <- metagenoma_cor_Fac$r[metagenome_important,]
metagenome.pval <- metagenoma_cor_Fac$p[metagenome_important,]

jpeg(filename = file.path(dir.res,"corr_metagenome_1_part.jpeg"))
corrplot::corrplot(metagenome.cor,p.mat=metagenome.pval,insig = "label_sig")
dev.off()


metaboloma_cor_Fac <- corr.test(metaboloma,factores,method = "spearman")
p.values <- metaboloma_cor_Fac$p
metaboloma_important <- rownames(metaboloma_cor_Fac$r)[which(apply(p.values,1,function(x) any(x<0.05))==T)]
metabolome.cor <- metaboloma_cor_Fac$r[metaboloma_important,]
metabolome.pval <- metaboloma_cor_Fac$p[metaboloma_important,]
jpeg(filename = file.path(dir.res,"corr_metabolome.jpeg"))
corrplot::corrplot(metabolome.cor,p.mat=metabolome.pval,insig = "label_sig")
dev.off()


pesos.metagenoma <- pesos[-grep("mean",ignore.case = T,rownames(pesos)),]
p.values.metagenoma <- metagenoma_cor_Fac$p



comparison <- function(x,y){
  
  x <- names(x)[x<0.05]
  y <- names(y)[abs(y)< 0.5]
  comparacion <- intersect(x,y)
  return(comparacion)
}

interseccion <- lapply(1:2, function(i) comparison(p.values.metagenoma[,i],pesos.metagenoma[,i]))

variables_metagenoma <- Reduce(unique,interseccion)



pesos.metaboloma <- pesos[grep("mean",ignore.case = T,rownames(pesos)),]
p.values.metaboloma <- metaboloma_cor_Fac$p
interseccion.boloma <- lapply(1:2, function(i) comparison(p.values.metaboloma[,i],pesos.metaboloma[,i]))

variables.metaboloma <- Reduce(unique,interseccion.boloma)

R <- factores %*% t(pesos)

factores.anal <- analyze_data(R,grupo,obesidad,correccion = 1)

funcion_pca(R,scle = T,trans = F)

csa <- sapply(1:1000,function(x) summary(aovp(factores[,2] ~ grupo*obesidad,perm = "Exact"))[[1]]$`Pr(Prob)`[1])
mean(csa)

# 
# 
# 
# 
# 
# 
# 
