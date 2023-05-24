library(glmnet)
library(caret)
library(ALDEx2)
library(MOFA2)
library(psych)
library(vegan)
library(mixOmics)
source("./scripts_R/scripts_utiles/scripts_funciones/otras_funciones_utiles.R")
source("./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R")
transformacion <- function(X) {
  return((t(log2(X / colSums(
    X
  )))))
}

mofa_componentes <- function(ncomp,semilla,directorio_modelo){
  idx <- sample(46)
  metaboloma <- (transformacion(datos$comunes$metaboloma))[idx,]
  metagenoma <- scale(t(mean_aldex(datos$comunes$microbiota$genero)))[idx,]
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

## Primero vemos cual es el numero optimo de componentes

semillas <- seq(40,80,10)

directorio_modelo <- "./scripts_R/integracion/modelo_optimal_comps"
if (!dir.exists(directorio_modelo)) {
  dir.create(directorio_modelo)
}

ncomps <- 2:11
if(length(list.files(directorio_modelo))==0){
  for(s in semillas){
    
    lapply(ncomps, function(x) mofa_componentes(ncomp = x,semilla = s,directorio_modelo = directorio_modelo ))
    
  }
}


modelos <- lapply(list.files(path = directorio_modelo,
                             full.names = T), load_model)

w <- which.min(compare_elbo(modelos,return_data = T)$ELBO)

modelo <- modelos[[w]]

set.seed(123896)
ncomp<- (get_factors(modelo)[[1]])
semillas <- sample(1000)
directorio.null <- "./scripts_R/integracion/perms"
if(!length(list.files(directorio.null))>2){
  
  lapply(semillas, function(x) mofa_componentes(ncomp = ncomp,semilla = x,directorio_modelo =directorio.null))

}


pesos <- Reduce(rbind,get_weights(modelo,scale = T))

factores <- get_factors(modelo,scale = T)[[1]]

p1 <- plot_variance_explained(modelo, max_r2=15)

p2 <- plot_variance_explained(modelo, plot_total = T)[[2]]

p3 <- ggarrange(p1,p2)

p3
dir.res <- "./scripts_R/integracion/resultados_mofa2"
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

metagenome.cor <- metagenoma_cor_Fac$r[metagenome_important[1:20],]
metagenome.pval <- metagenoma_cor_Fac$p[metagenome_important[1:20],]

jpeg(filename = file.path(dir.res,"corr_metagenome_1_part.jpeg"))
corrplot::corrplot(metagenome.cor,p.mat=metagenome.pval,insig = "label_sig")
dev.off()
metagenome.cor <- metagenoma_cor_Fac$r[metagenome_important[21:length(metagenome_important)],]
metagenome.pval <- metagenoma_cor_Fac$p[metagenome_important[21:length(metagenome_important)],]
jpeg(filename = file.path(dir.res,"corr_metagenome_2_part.jpeg"))
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

interseccion <- lapply(1:7, function(i) comparison(p.values.metagenoma[,i],pesos.metagenoma[,i]))

variables_metagenoma <- Reduce(unique,interseccion)



pesos.metaboloma <- pesos[grep("mean",ignore.case = T,rownames(pesos)),]
p.values.metaboloma <- metaboloma_cor_Fac$p
interseccion.boloma <- lapply(1:7, function(i) comparison(p.values.metaboloma[,i],pesos.metaboloma[,i]))

variables.metaboloma <- Reduce(unique,interseccion.boloma)

R <- factores %*% t(pesos)

factores.anal <- analyze_data(R,grupo,obesidad,correccion = 1)

funcion_pca(R,scle = F,trans = F)

csa <- sapply(1:1000,function(x) summary(aovp(factores[,3] ~ grupo*obesidad,perm = "Exact"))[[1]]$`Pr(Prob)`[3])



### Ahora veremos una aproximacion por reg logistica








# # Comparacion R y O
# 
# reconstruccion <- factores %*% t(pesos)
# original <- cbind(metaboloma,metagenoma)
# 
# clinicas <- datos$comunes$clinicos[,1:8]
# 
# original_clinico <- as.data.frame(cbind(original,clinicas))
# 
# correlacion_par <- sapply(1:ncol(reconstruccion), 
#               function(x) sapply(1:ncol(original),
#                                  function(y) sapply(1:ncol(clinicas),
#                                                     function(z) ppcor::pcor.test(reconstruccion[,x],original[,y],clinicas[,z],method = "pearson")$estimate)))
# correlacion_par.pvalue <- sapply(1:ncol(reconstruccion), 
#                           function(x) sapply(1:ncol(original),
#                                              function(y) sapply(1:ncol(clinicas),
#                                                                 function(z) ppcor::pcor.test(reconstruccion[,x],original[,y],clinicas[,z],method = "pearson")$p.value)))
# cosa <- corr.test(reconstruccion,original)
# 
# p.value <- cosa$p
# p.value[lower.tri(p.value)] <- cosa$p.adj
# csa  <- cosa$p.adj
# 
# rownames(csa)[csa<0.05]
# 
# run_tsne(modelo,perplexity=4)
# 
# csa <- plot_dimred(object = modelo,method = "TSNE",perplexity=15)
# datos2 <- csa$data
# 
# datos2$grupo <- grupo
# datos2$obesidad <- obesidad
# 
# ggplot(datos2,aes(x,y,col=grupo,shape=obesidad))+geom_point()
# # ## correlacion R y O
# # 
# # cor.fo <- corr.test(reconstruccion,original,adjust = "BH")
# # cor.o <- corr.test(original)$r
# # cor.r <- corr.test(reconstruccion)$r
# # 
# # corr_AB <- cor.fo$r
# # corr_AA <- cor.r
# # corr_BB <- cor.r
# # 
# # partial_corr_AB_given_B = -ginv(corr_AA) %*% corr_AB %*% ginv(corr_BB)
# # 
# # La otra aproximacion es hacer mofa por grupo y obesidad. Hacer unbootstrap anova
# # p.values <- cor.fo$p.adj
# # 
# # ### ahora vemos las variables más correlacionadas con los factores
# # 
# # variables.cor1 <- rownames(p.values)[which(apply(p.values,
# #                                                  1,
# #                                                  function(x) any(x<0.05)))]
# # 
# # 
# # library(ppcor)
# # 
# # sapply(1:(ncol(dat)-1), function(x) sapply(1:(ncol(dat)-1), function(y) {
# #   if (x == y) 1
# #   else pcor.test(dat[,x], dat[,y], dat[,ncol(dat)])$estimate
# # }))
# 
# # univ.raw <- extract_univ(original,correc = 2)
# # 
# # univ.r <- extract_univ(reconstruccion,correc = 1)
# # 
# # 
# # X <- as.matrix(univ.r$media)
# # P.r <- as.matrix(univ.r$pvalores)
# # P.r[P.r>0.01]<- NA
# # P.r <- P.r[which(apply(P.r,1,function(x)any(!is.na(x)))),]
# # X <- X[rownames(X) %in% rownames(P.r),]
# # 
# # colnames(X) <- colnames(P.r)
# # 
# # corrplot::corrplot(X,
#                    is.corr = F,
#                    p.mat = P.r,
#                    sig.level = 0.01,
#                    insig = "label_sig",
#                    tl.cex = 0.5,bg = ifelse(X>0,"lightblue","pink"))
# 
# 
# univ.r <- extract_univ(original,correc = 2)
# 
# 
# X <- as.matrix(univ.r$media)
# P <- as.matrix(univ.r$pvalores)
# P <- P[rownames(P) %in% rownames(P.r),]
# X <- X[rownames(X) %in% rownames(P.r),]
# 
# corrplot::corrplot(X,
#                    is.corr = F,
#                    p.mat = P,
#                    sig.level = 0.01,
#                    insig = "label_sig",
#                    tl.cex = 0.5,bg = ifelse(X>0,"lightblue","pink"),
#                    )
# 
# 
# 
# univ.r <- extract_univ(reconstruccion,correc = 1)
# univ.o <- extract_univ(original,correc = 2)
# 
# 
# correlacion_grupos <- corr.test(univ.r$media,univ.o$media,adjust = "none")
# correlacion_grupos.p <- correlacion_grupos$p
# correlacion_grupos.p[lower.tri(correlacion_grupos.p)] <- correlacion_grupos.p[upper.tri(correlacion_grupos.p)]
# correlacion_grupos.r <- correlacion_grupos$r
# 
# 
# corrplot::corrplot(correlacion_grupos.r,
#                    is.corr = F,
#                    p.mat = correlacion_grupos.p,
#                    sig.level = 0.01,
#                    insig = "label_sig",
#                    tl.cex = 0.5,bg = ifelse(X>0,"lightblue","pink"),
# )
# 
# 
# correlacion_grupos <- corr.test(t(univ.r$media),t(univ.o$media),adjust = "none")
# correlacion_grupos.p <- correlacion_grupos$p
# correlacion_grupos.r <- correlacion_grupos$r
# 
# correlacion_grupos.p[(correlacion_grupos.p)>0.05]<- NA
# 
# 
# 
# 
# corrplot::corrplot(correlacion_grupos.r[w,w],
#                    is.corr = F,
#                    p.mat = correlacion_grupos.p[w,w],
#                    sig.level = 0.01,
#                    insig = "label_sig",
#                    tl.cex = 0.5,bg = ifelse(X>0,"lightblue","pink"),
# )
# 
# 
# correlacion.o <- corr.test(factores,original,adjust = "none")
# 
# w <- which(apply(correlacion.o$p,2,function(x) any(x<0.05)))
# 
# O.r <- correlacion.o$r[,w]
# O.p <- correlacion.o$p[,w]
# 
# corrplot::corrplot(O.r,
#                    is.corr = F,
#                    p.mat = O.p,
#                    sig.level = 0.01,
#                    insig = "label_sig",
#                    tl.cex = 0.5,bg = ifelse(X>0,"lightblue","pink"),
# )
# 
# 
# 
# 
# correlacion.o <- corr.test(reconstruccion,original,adjust = "none")
# 
# w <- which(apply(correlacion.o$p,2,function(x) any(x<0.01)))
# 
# O.r <- correlacion.o$r[w,w]
# O.p <- correlacion.o$p[w,w]
# 
# corrplot::corrplot(O.r,
#                    is.corr = F,
#                    p.mat = O.p,
#                    sig.level = 0.01,
#                    insig = "label_sig",
#                    tl.cex = 0.5,bg = ifelse(X>0,"lightblue","pink"),
# )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
