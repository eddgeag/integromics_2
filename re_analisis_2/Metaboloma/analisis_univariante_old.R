library(ggplot2)
library(dplyr)


transform.in.t <- function(x) scale(t(x))

dat <- readRDS("../datos/preprocesado_08_09_23/novoom.rds")

tot<- as.data.frame(transform.in.t(dat$totales$metaboloma))

 
grupo <-dat$totales$grupo
obesidad <- dat$totales$obesidad

tot$grupo <- grupo
tot$obesidad <- obesidad
tot.melt <- reshape2::melt(tot)

ggplot(tot.melt,aes(value,color=obesidad))+geom_density()

library(limma)
library(psych)




tot<- as.data.frame(transform.in.t(dat$totales$metaboloma))

mardia(tot)


mdl <- model.matrix(~0+interaction(grupo,obesidad))
colnames(mdl) <- c("FNO","PNO","MNO","FO","PO","MO")
x <- c("PNO+PO-FO-FNO",
       "MNO+MO-FO-FNO",
       "PNO+PO-MNO-MO",
       "PNO-FNO",
       "PNO-MNO",
       "MNO-FNO",
       "PO-FO",
       "PO-MO",
       "MO-FO")
contrastes <- makeContrasts(contrasts = x,levels = colnames(mdl))


fit <- lmFit(t(tot),design = mdl)
fit2 <- contrasts.fit(fit,contrastes)
fit2. <- eBayes(fit2)
pcos.females <- topTable(fit2.,coef = 1,number = Inf)
# males.females <- topTable(fit2,coef = 2,number = Inf)
# pcos.males <- topTable(fit2,coef = 3,number = Inf)
# pcos.females.no <- topTable(fit2,coef = 4,number = Inf)
# pcos.males.no <- topTable(fit2,coef = 5,number = Inf)
# males.females.no <- topTable(fit2,coef = 6,number = Inf)
# pcos.females.o <- topTable(fit2,coef = 7,number = Inf)
# pcos.males.o<- topTable(fit2,coef = 8,number = Inf)
# males.females.o<- topTable(fit2,coef = 9,number = Inf)

modelo1 <- coef(fit2.)

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
fit2_ <- eBayes(fit2)
csa2 <- topTable(fit2_,2)

modelo2 <- coef(fit2_,Inf)

colnames(modelo2) <- names(contraste)

csa3 <- topTable(fit2_,number = Inf,coef=length(contraste)-3)

colnames(csa3)[1:length(contraste)] <- names(contraste)

