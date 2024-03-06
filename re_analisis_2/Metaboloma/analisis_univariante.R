library(ggplot2)
library(dplyr)
library(limma)
library(psych)

add_classificacion_metaboloma <- function(X) {
  aromatic <- rownames(X)[c(2, 3, 4, 1, 6, 35, 36, 34)]
  other <-
    rownames(X)[c(21, 9, 16, 26, 22, 18, 13, 29, 30, 27)]
  aa <- rownames(X)[c(33, 15, 20, 19, 23)]
  carbo <-
    rownames(X)[!rownames(X) %in% c(aromatic, other, aa)]
  
  X$grupo_meta <-
    ifelse(rownames(X) %in% aromatic, "AROMATICO", NA)
  X$grupo_meta <-
    ifelse(rownames(X) %in% other,
           "OTHER",
           X$grupo_meta)
  X$grupo_meta <-
    ifelse(rownames(X) %in% aa,
           "DERIVADOS AA",
           X$grupo_meta)
  X$grupo_meta <-
    ifelse(
      rownames(X) %in% carbo,
      "CARBOHIDRATOS_GRASAS_KETONA_GLYCEROL",
      X$grupo_meta
    )
  
  return(X)
}

transform.in.t <- function(x) scale(t(x))

dat <- readRDS("../datos/preprocesado_08_09_23/novoom.rds")

tot<- as.data.frame(transform.in.t(dat$totales$metaboloma))


grupo <-dat$totales$grupo
obesidad <- dat$totales$obesidad

tot$grupo <- grupo
tot$obesidad <- obesidad
tot.melt <- reshape2::melt(tot)

ggplot(tot.melt,aes(value,color=obesidad))+geom_density()




tot<- as.data.frame(transform.in.t(dat$totales$metaboloma))

mardia(tot)

mdl2 <- model.matrix(~0+obesidad*grupo)

colnames(mdl2) <- c("NO","Obesidad","PCOS","Male","PCOSObese","MaleObese")
contraste <- c(pcosfemale="PCOS-(Obesidad+NO)",
               malefemale="Male-(Obesidad+NO)",
               pcosmale="PCOS-Male",
               obesidad="Obesidad-NO",
               pcosfamaleno="(PCOS-Obesidad)-NO",
               malefemaleno="(Male-Obesidad)-NO",
               pcosmaleobeno="(PCOS-Obesidad)-(Male-Obesidad)",
               pcosfemaleob = "PCOSObese-Obesidad",
               malefemaleob="MaleObese-Obesidad",
               pcosmaleob = "MaleObese-PCOSObese",
               pcosobese = "PCOSObese",
               MaleObese="MaleObese",
               FemaleObese="Obesidad"
)


contrastes <- makeContrasts(contrasts = contraste,levels = colnames(mdl2))
fit <- lmFit(t(tot),design = mdl2)
fit2 <- contrasts.fit(fit,contrastes)
fit2 <- eBayes(fit2)
csa2 <- topTable(fit2,2)

modelo2 <- coef(fit2,Inf)
modelo2 <- add_classificacion_metaboloma(as.data.frame(modelo2))
colnames(modelo2)[ncol(modelo2)] <- "GRUPO"

colnames(modelo2) <- names(contraste)
colnames(modelo2)[ncol(modelo2)] <- "GRUPO"

aux.fun <- function(toptable,mdl) {
      
  obesidad <- add_classificacion_metaboloma(topTable(fit = fit2,coef = 4,number = Inf))
  mdl.obesidad <- modelo2[obesidad$adj.P.Val<0.05,c("obesidad","GRUPO")]
  mdl.obesidad$signficado <- ifelse(mdl.obesidad$obesidad>0,"Aumenta","Disminuye")
  mdl.obesidad$nombre <- rownames(mdl.obesidad)
  
}
