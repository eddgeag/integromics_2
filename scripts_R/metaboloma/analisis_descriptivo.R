
library(ggplot2)

library(reshape2)
source("./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R")
source("./scripts_R/scripts_utiles/scripts_funciones/calculo_medianas.R")
datos <- readRDS("../datos/preprocesado_08_09_23/novoom.rds")
set.seed(126581)


X <- as.data.frame(t(datos$totales$metaboloma))
X$succinato_mean <- rowMeans(X[,grep("succinate_0",ignore.case = T,x=colnames(X))])
X <- X[,-grep("Succinate",x=colnames(X))]


## oxoisocaproic 

hist(X$Oxoisocaproic_mean)

## leucine
hist(X$Leu_mean)

## Isuoleu
hist(X$Isoleu_mean)

hist(X$Valine_mean)

hist(X$Isobutyric_mean)

hist(X$Oxoisovaleric_mean)

hist(X$OHbutyric_mean)

hist(X$Lactate_mean)

hist(X$Alanine_mean)

hist(X$Acetate_mean)

hist(X$Glycoprot_mean)

hist(X$Acetone_mean)

hist(X$Gluacid_mean)

hist(X$Pyruvate_mean)

hist(X$Pyroglu_mean)

hist(X$Gln_mean)

hist(X$Citrate_mean)

hist(X$Asn_mean)

hist(X$Creatine_mean)

hist(X$Creatinine_mean)

hist(X$Lysine_mean)

hist(X$Ornithine_mean)

hist(X$Choline_mean)

hist(X$Carnitine_mean)

hist(X$Betaine_mean)

hist(X$Glycine_mean)

hist(X$Threonine_mean)

hist(X$Glycerol_mean)

hist(X$Serine_mean)

hist(X$Proline_mean)

hist(X$BGlucose_mean)

hist(X$DGlucose_mean)

hist(X$Methys_mean)

hist(X$Tyrosine_mean)

hist(X$Phe_mean)

hist(X$Tryp_mean)

hist(X$Formate_mean)

hist(log10(X$succinato_mean))


X.scale <- (log2(X/colSums(X)))



grupo <- datos$totales$grupo
obesidad <- datos$totales$obesidad

X.rs <- reshape2::melt(data.frame(X.scale
,grupo=grupo,obesidad=obesidad))

ggplot(X.rs,aes(value,color=grupo))+geom_density()

pcx <- prcomp(X.scale,scale. = T)
plts <- data.frame(pcx$x,grupo=grupo,obesidad=obesidad)
ggplot(plts,aes(PC1,PC2,color=grupo))+geom_point()

ks.test(X$Oxoisocaproic_mean,"")

