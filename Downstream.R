library(ggplot2)
library(psych)
library(limma)
library(car)
library(mixOmics)
library(vegan)
library(ggpubr)
library(MOFA2)
library(dplyr)
library(ALDEx2)
library(PERMANOVA)
library(ggrepel)
library(MKinfer)
library(patchwork)
library(forcats)

# library(dicosar)
source("./scripts_R/scripts_utiles/scripts_funciones/vanvalen2.R")
source("./scripts_R/scripts_utiles/scripts_funciones/manova_vanvalen.R")
source("./scripts_R/scripts_utiles/scripts_funciones/otras_funciones_utiles.R")
datos <- readRDS("../datos/preprocesado_08_09_23/novoom-04-10-23.rds")
# soi <-c("PdGLP6","MdGLP4")
# obesidad <- datos$comunes$obesidad[!rownames(datos$comunes$clinicos) %in% soi]
# grupo <- datos$comunes$grupo[!rownames(datos$comunes$clinicos) %in% soi]
# sexo <- datos$comunes$variables_in_bacteria$SEX[!rownames(datos$comunes$clinicos) %in% soi]
# Función para aproximar la inversa de una matriz utilizando la serie de Neumann
neumann_inverse <- function(A, tol = 1e-6, max_iter = 1000) {
  # Verificar que A es una matriz cuadrada
  if (!is.matrix(A) || nrow(A) != ncol(A)) {
    stop("La matriz A debe ser cuadrada.")
  }
  
  n <- nrow(A)
  I <- diag(n)
  
  # Calcular la norma de A
  norm_A <- norm(A, type = "2")
  
  # Escalar la matriz A si su norma es mayor o igual a 1
  if (norm_A >= 1) {
    alpha <- 1 / (norm_A + 1)
    A_scaled <- alpha * A
  } else {
    alpha <- 1
    A_scaled <- A
  }
  
  # Matriz M para la serie de Neumann
  M <- I - A_scaled
  
  # Iniciar la serie
  X <- I
  term <- I
  for (k in 1:max_iter) {
    print(k)
    term <- term %*% M
    X <- X + term
    # Verificar la convergencia
    if (norm(term, type = "F") < tol) {
      cat("Convergencia alcanzada en", k, "iteraciones.\n")
      break
    }
  }
  
  # Aproximación de la inversa
  inv_A_approx <- alpha * X
  
  # Comprobar si se alcanzó la convergencia
  if (k == max_iter && norm(term, type = "F") >= tol) {
    warning("No se alcanzó la convergencia en el número máximo de iteraciones.")
  }
  
  return(inv_A_approx)
}

mahalanobis_dist <- function(x, y, S_inv) {
  # x y y son los dos vectores
  # S_inv es la matriz inversa de varianzas y covarianzas (ya aproximada con Taylor)
  
  # Diferencia entre los vectores
  diff <- x - y
  
  
  # Calculo de la distancia de Mahalanobis
  dist <- sqrt(t(diff) %*% S_inv %*% diff)
  
  return(dist)
}
obesidad <- datos$comunes$obesidad
grupo <- datos$comunes$grupo
sexo <- datos$comunes$variables_in_bacteria$SEX
fun_plot_limma <- function(X) {
  # tot <- X <- factores.df
  tot <- X
  mdl2 <- model.matrix( ~ 0 + obesidad)
  colnames(mdl2) <- c("No", "O")
  contraste <- c("O-No")
  
  contrastes <-
    makeContrasts(contrasts = contraste, levels = colnames(mdl2))
  fit <- lmFit(t(tot), design = mdl2)
  
  fit2 <- contrasts.fit(fit, contrastes)
  fit2 <- eBayes(fit2)
  top_table <- topTable(fit2, coef = 1, number = Inf)
  logFC <- top_table$logFC ## obtengo las medias
  tstat <- top_table$t
  res_obesidad <- cbind(logFC, tstat, adj.P.Val = top_table$P.Value)
  mdl2 <- model.matrix( ~ 0 + grupo)
  colnames(mdl2) <- c("Female", "PCOS", "Male")
  contraste <- c("Male-Female", "PCOS-Female", "PCOS-Male")
  
  contrastes <-
    makeContrasts(contrasts = contraste, levels = colnames(mdl2))
  fit <- lmFit(t(tot), design = mdl2)
  fit2 <- contrasts.fit(fit, contrastes)
  fit2 <- eBayes(fit2)
  
  res_grupo <- lapply(1:3, function(x)
    topTable(fit2, coef = x, number = Inf))
  res_grupo <- lapply(res_grupo, function(x)
    cbind(x$logFC, x$t, x$P.Value))
  names(res_grupo) <-  c("PF", "PM", "MF")
  interaccion <- interaction(grupo, obesidad)
  mdl2 <- model.matrix( ~ 0 + interaccion)
  colnames(mdl2) <- gsub("interaccion", "", colnames(mdl2))
  colnames(mdl2) <- c("FNO", "PCOSNO", "MALENO", "FO", "PCOSO", "MALEO")
  contraste <- c(
    female = "FO-FNO",
    pcos = "PCOSO-PCOSNO",
    male = "MALEO-MALENO",
    no_sexo_control = "MALENO-FNO",
    no_pcos_female = "PCOSNO-FNO",
    no_pcos_male = "PCOSNO-MALENO",
    o_pcos_female = "PCOSO-FO",
    o_pcos_male = "PCOSO-MALEO",
    o_sexo_control = "MALEO-FO"
    
  )
  
  colnames(mdl2) <- gsub("interaccion", "", colnames(mdl2))
  
  contrastes <-
    makeContrasts(contrasts = contraste, levels = colnames(mdl2))
  fit <- lmFit(t(tot), design = mdl2)
  fit2 <- contrasts.fit(fit, contrastes)
  fit2 <- eBayes(fit2)
  
  res_interaccion <- lapply(1:dim(coef(fit2))[2], function(x)
    topTable(fit2, coef = x, number = Inf))
  res_interaccion <- lapply(res_interaccion, function(x)
    cbind(x$logFC, x$t, x$P.Value))
  names(res_interaccion) <- colnames(coef(fit2))
  
  
  
  logfc_ <- data.frame(
    res_obesidad[, 1],
    sapply(res_grupo, function(x)
      x[, 1]),
    sapply(res_interaccion, function(x)
      x[, 1])
  )
  
  tstats_ <-  data.frame(
    res_obesidad[, 2],
    sapply(res_grupo, function(x)
      x[, 2]),
    sapply(res_interaccion, function(x)
      x[, 2])
  )
  
  
  
  pvals <-  data.frame(
    res_obesidad[, 3],
    sapply(res_grupo, function(x)
      x[, 3]),
    sapply(res_interaccion, function(x)
      x[, 3])
  )
  
  
  rownames(tstats_) <- colnames(tot)
  tstats_$Variable <- rownames(tstats_)
  orden_contrastes <- c(
    "O.No",
    "MF",
    "PF",
    "PM",
    "FO.FNO",
    "PCOSO.PCOSNO",
    "MALEO.MALENO",
    "MALENO.FNO",
    "PCOSNO.FNO",
    "PCOSNO.MALENO",
    "MALEO.FO",
    "PCOSO.FO",
    "PCOSO.MALEO"
  )
  colnames(tstats_) <- c(orden_contrastes, "Variable")
  tstats_$Variable <- factor(tstats_$Variable, levels = colnames(tot))
  
  tstats <- reshape2::melt(tstats_, id.vars = "Variable")
  tstats$variable <- factor(tstats$variable , levels = colnames(tstats_)[-length(colnames(tstats_))])
  # tstats$variable <- factor(tstats$variable,levels = orden_contrastes)
  
  # colnames(csa) <-c("variables","sujetos","categoria","")
  pvalues <- pvals
  rownames(pvalues) <- colnames(tot)
  colnames(pvalues) <- orden_contrastes
  pvalues$Variable <- colnames(tot)
  pvalues$Variable  <- factor(pvalues$Variable, levels = colnames(tot))
  pvalues <- reshape2::melt(pvalues, id.vars = "Variable")
  pvalues$variable <- factor(tstats$variable, levels = colnames(tstats_)[-length(colnames(tstats_))])
  
  # colnames(csa) <-c("variables","sujetos","categoria","")
  pvalues$text <- ifelse(pvalues$value < 0.05, "*", "")
  tstats$text <- ifelse(pvalues$value < 0.05, "*", "")
  logfc <- data.frame(
    res_obesidad[, 1],
    sapply(res_grupo, function(x)
      x[, 1]),
    sapply(res_interaccion, function(x)
      x[, 1])
  )
  rownames(logfc) <- colnames(tot)
  colnames(logfc) <- orden_contrastes
  
  logfc$Variable <- colnames(tot)
  logfc$Variable  <- factor(logfc$Variable, levels = colnames(tot))
  logfc_ <- reshape2::melt(logfc, id.vars = "Variable")
  logfc_$variable <- factor(logfc_$variable, levels = orden_contrastes)
  
  
  
  logfc_$text <- ifelse(pvalues$value < 0.05, "*", "")
  
  
  
  return(list(
    pvals = pvals,
    logFC = logfc_,
    tstats = tstats
  ))
  
}
misuper3 <- function(tot) {
  ### obesidad
  mdl2 <- model.matrix(~ 0 + obesidad)
  colnames(mdl2) <- c("No", "O")
  contraste <- c("O-No")
  
  contrastes <-
    makeContrasts(contrasts = contraste, levels = colnames(mdl2))
  fit <- lmFit(t(tot), design = mdl2)
  
  fit2 <- contrasts.fit(fit, contrastes)
  fit2 <- eBayes(fit2)
  coeficientes <- coef(fit2) ## obtengo las medias
  coeficientes.obesidad <- coeficientes ## medias
  cn <- dim(coeficientes)[2]
  res.pdf <-
    lapply(1:cn, function(x)
      topTable(fit2, coef = x, number = Inf))[[1]]
  res.pdf2 <- res.pdf$P.Value ## obtengo el o p valor
  res.pdf2 <- as.data.frame(res.pdf2)
  colnames(res.pdf2) <- "O.No"
  rownames(res.pdf2) <- rownames(coeficientes)
  obese <- res.pdf2 ## pvalor
  obese.LF <- sapply(cn, function(x)
    topTable(fit2, coef = x, number = Inf)$logFC)
  obese.fit <- apply(tot, 2, function(y)
    auxlm(y, obesidad))
  mdl2 <- model.matrix(~ 0 + grupo)
  colnames(mdl2) <- c("Female", "PCOS", "Male")
  contraste <- c("PCOS-Female", "PCOS-Male", "Male-Female")
  contrastes <-
    makeContrasts(contrasts = contraste, levels = colnames(mdl2))
  fit <- lmFit(t(tot), design = mdl2)
  fit2 <- contrasts.fit(fit, contrastes)
  fit2 <- eBayes(fit2)
  coeficientes <- coef(fit2)
  coefficientes.grupo <- coeficientes ## medias
  cn <- dim(coeficientes)[2]
  res.pdf <-
    lapply(1:cn, function(x)
      as.data.frame(topTable(
        fit2, coef = x, number = Inf
      )))
  res.pdf2 <- Reduce("cbind", lapply(res.pdf, function(x)
    x$P.Value))
  ## grupo pvalor
  colnames(res.pdf2) <- c("PF", "PM", "MF")
  rownames(res.pdf2) <- rownames(coeficientes)
  res.pdf2 <- as.data.frame(res.pdf2)
  grups <- res.pdf2 ## p valor
  grups.fit <- apply(tot, 2, function(y)
    auxlm(y, grupo))
  grups.LF <- sapply(1:cn, function(x)
    topTable(fit2, coef = x, number = Inf)$logFC)
  
  interaccion <- interaction(grupo, obesidad)
  mdl2 <- model.matrix(~ 0 + interaccion)
  
  colnames(mdl2) <- gsub("interaccion", "", colnames(mdl2))
  colnames(mdl2) <- c("FNO", "PCOSNO", "MALENO", "FO", "PCOSO", "MALEO")
  contraste <- c(
    no_pcos_female = "PCOSNO-FNO",
    no_pcos_male = "PCOSNO-MALENO",
    no_sexo_control = "MALENO-FNO",
    o_pcos_female = "PCOSO-FO",
    o_pcos_male = "PCOSO-MALEO",
    o_sexo_control = "MALEO-FO",
    pcos = "PCOSO-PCOSNO",
    male = "MALEO-MALENO",
    female = "FO-FNO"
  )
  
  colnames(mdl2) <- gsub("interaccion", "", colnames(mdl2))
  
  contrastes <-
    makeContrasts(contrasts = contraste, levels = colnames(mdl2))
  fit <- lmFit(t(tot), design = mdl2)
  fit2 <- contrasts.fit(fit, contrastes)
  fit2 <- eBayes(fit2)
  coeficientes <- coef(fit2, Inf)
  coeficientes.interact <- coeficientes ## medias
  cn <- dim(coeficientes)[2]
  res.pdf <-
    lapply(1:cn, function(x)
      as.data.frame(topTable(
        fit2, coef = x, number = Inf
      )))
  res.pdf2 <- Reduce("cbind", lapply(res.pdf, function(x)
    x$adj.P.Val))## pvalor
  rownames(res.pdf2) <- rownames(coeficientes)
  colnames(res.pdf2) <- colnames(coeficientes)
  interact <- res.pdf2
  interact.fit <- apply(tot, 2, function(y)
    auxlm(y, interaccion))
  
  interaccion.LF <- sapply(1:cn, function(x)
    topTable(fit2, coef = x, number = Inf)$logFC)
  
  pvalores <- data.frame(obese, grups, interact)
  
  medias <- data.frame(coeficientes.obesidad,
                       coefficientes.grupo,
                       coeficientes.interact)
  logFC <- data.frame(obese.LF, grups.LF, interaccion.LF)
  colnames(logFC) <- colnames(medias)
  colnames(pvalores) <- colnames(medias)
  medias$variables <- rownames(medias)
  pvalores$variables <- rownames(pvalores)
  rownames(logFC) <- rownames(pvalores)
  logFC$variables <- rownames(pvalores)
  retorno <- list(
    pvalores = pvalores,
    medias = medias,
    gof = list(
      obesidad = obese.fit,
      grupo = grups.fit,
      interaccion = interact.fit
    ),
    LF = logFC
  )
  
}
fun_summary <- function(prueba) {
  p <- ncol(prueba$medias)
  k <- nrow(prueba$medias)
  obese.gof <- as.data.frame(t(prueba$gof$obesidad))
  grupo.gof <- as.data.frame(t(prueba$gof$grupo))
  inter.gof <- as.data.frame(t(prueba$gof$interaccion))
  
  O.No <- data.frame(
    obese.gof,
    mean = prueba$medias$O.No,
    pval = prueba$pvalores$O.No,
    logFC = prueba$LF$O.No
  )
  
  PCOS.Female <- data.frame(
    grupo.gof,
    mean = prueba$medias$PCOS.Female,
    pval = prueba$pvalores$PCOS.Female,
    logFC = prueba$LF$PCOS.Female
  )
  PCOS.Male <- data.frame(
    grupo.gof,
    mean = prueba$medias$PCOS.Male,
    pval = prueba$pvalores$PCOS.Male,
    logFC = prueba$LF$PCOS.Male
  )
  Female.Male <- data.frame(
    grupo.gof,
    mean = prueba$medias$Male.Female,
    pval = prueba$pvalores$Male.Female,
    logFC = prueba$LF$Male.Female
  )
  MALENO.FNO <-  data.frame(
    inter.gof,
    mean = prueba$medias$MALENO.FNO,
    pval = prueba$pvalores$MALENO.FNO,
    logFC = prueba$LF$MALENO.FNO
  )
  PCOSNO.FNO <-  data.frame(
    inter.gof,
    mean = prueba$medias$PCOSNO.FNO,
    pval = prueba$pvalores$PCOSNO.FNO,
    logFC = prueba$LF$PCOSNO.FNO
  )
  PCOSNO.MALENO <-  data.frame(
    inter.gof,
    mean = prueba$medias$PCOSNO.MALENO,
    pval = prueba$pvalores$PCOSNO.MALENO,
    logFC = prueba$LF$PCOSNO.MALENO
  )
  
  MALEO.FO <-  data.frame(
    inter.gof,
    mean = prueba$medias$MALEO.FO,
    pval = prueba$pvalores$MALEO.FO,
    logFC = prueba$LF$MALEO.FO
  )
  PCOSO.FO <-  data.frame(
    inter.gof,
    mean = prueba$medias$PCOSO.FO,
    pval = prueba$pvalores$PCOSO.FO,
    logFC = prueba$LF$PCOSO.FO
  )
  PCOSO.MALEO <-  data.frame(
    inter.gof,
    mean = prueba$medias$PCOSO.MALEO,
    pval = prueba$pvalores$PCOSO.MALEO,
    logFC = prueba$LF$PCOSO.MALEO
  )
  
  PCOSO.PCOSNO <- data.frame(
    inter.gof,
    mean = prueba$medias$PCOSO.PCOSNO,
    pval = prueba$pvalores$PCOSO.PCOSNO,
    logFC = prueba$LF$PCOSO.PCOSNO
  )
  
  FO.FNO <- data.frame(
    inter.gof,
    mean = prueba$medias$FO.FNO,
    pval = prueba$pvalores$FO.FNO,
    logFC = prueba$LF$FO.FNO
  )
  MALEO.MALENO <- data.frame(
    inter.gof,
    mean = prueba$medias$MALEO.MALENO,
    pval = prueba$pvalores$MALEO.MALENO,
    logFC = prueba$LF$MALEO.MALENO
  )
  
  
  retorno <- list(
    O.No = O.No,
    PCOS.Female = PCOS.Female,
    PCOS.Male = PCOS.Male,
    Female.Male = Female.Male,
    PCOSNO.FNO = PCOSNO.FNO,
    PCOSNO.MALENO = PCOSNO.MALENO,
    MALENO.FNO = MALENO.FNO,
    PCOSO.FO = PCOSO.FO,
    PCOSO.MALEO = PCOSO.MALEO,
    MALEO.FO = MALEO.FO,
    PCOSO.PCOSNO = PCOSO.PCOSNO,
    FO.FNO = FO.FNO,
    MALEO.MALENO = MALEO.MALENO
  )
  return(retorno)
  
}
funcion_pca <- function(X, scle, trans) {
  if (!trans) {
    pca <- prcomp(X, scale. = scle)
    
    varianzas <- round(100 * pca$sdev ^ 2 / sum(pca$sdev ^ 2), 2)
    pca.plot <- data.frame(
      PC1 = pca$x[, 1],
      PC2 = pca$x[, 3],
      obesidad = obesidad,
      grupo = grupo,
      sexo = sexo
    )
  } else{
    X <- t(scale(X))
    pca <- prcomp((X), scale. = scle)
    
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
    ggplot(pca.plot, aes(PC1, PC2, color = grupo)) + geom_point() + xlab(paste("PC1", varianzas[1], "%")) +
    ylab(paste("PC2", varianzas[2], "%")) + theme(legend.title = element_blank())
  p2 <-
    ggplot(pca.plot, aes(PC1, PC2, color = grupo)) + geom_point() + xlab(paste("PC1", varianzas[1], "%")) +
    ylab(paste("PC2", varianzas[2], "%")) + facet_grid( ~ obesidad) + theme(legend.title = element_blank())
  p3 <-
    ggplot(pca.plot, aes(PC1, PC2, color = sexo)) + geom_point() + xlab(paste("PC1", varianzas[1], "%")) +
    ylab(paste("PC2", varianzas[2], "%")) + facet_grid( ~ obesidad) + theme(legend.title = element_blank())
  p44 <- factoextra::fviz_eig(pca, "variance")
  p4 <- ggarrange(p1, p2, p3, p44)
  p1 <-
    ggplot(pca.plot, aes(PC1, PC2, color = obesidad)) + xlab(paste("PC1", varianzas[1], "%")) +
    ylab(paste("PC2", varianzas[2], "%")) + geom_point() + theme(legend.title = element_blank())
  p2 <-
    ggplot(pca.plot, aes(PC1, PC2, color = obesidad)) + xlab(paste("PC1", varianzas[1], "%")) +
    ylab(paste("PC2", varianzas[2], "%")) + geom_point() + facet_grid( ~ grupo) +
    theme(legend.title = element_blank())
  p3 <-
    ggplot(pca.plot, aes(PC1, PC2, color = obesidad)) + xlab(paste("PC1", varianzas[1], "%")) +
    ylab(paste("PC2", varianzas[2], "%")) + geom_point() + facet_grid( ~ sexo) +
    theme(legend.title = element_blank())
  
  p44 <- factoextra::fviz_eig(pca, "variance")
  p5 <- ggarrange(p1, p2, p3, p44)
  p1 <-
    ggplot(pca.plot, aes(PC1, PC2, color = sexo)) + xlab(paste("PC1", varianzas[1], "%")) +
    ylab(paste("PC2", varianzas[2], "%")) + geom_point() + theme(legend.title = element_blank())
  p2 <-
    ggplot(pca.plot, aes(PC1, PC2, color = sexo)) + xlab(paste("PC1", varianzas[1], "%")) +
    ylab(paste("PC2", varianzas[2], "%")) + geom_point() + facet_grid( ~ grupo) +
    theme(legend.title = element_blank())
  p3 <-
    ggplot(pca.plot, aes(PC1, PC2, color = sexo)) + xlab(paste("PC1", varianzas[1], "%")) +
    ylab(paste("PC2", varianzas[2], "%")) + geom_point() + facet_grid( ~ sexo) +
    theme(legend.title = element_blank())
  p44 <- factoextra::fviz_eig(pca, "variance")
  
  p6 <- ggarrange(p1, p2, p3, p44)
  
  return(list(p4, p5, p6))
}

auxlm <- function(x, categorica) {
  grups.lm <-  summary(lm(x ~ 0 + categorica))
  grups.r2 <- grups.lm$r.squared
  fstat <- grups.lm$fstatistic[1]
  df1 <- grups.lm$fstatistic[2]
  df2 <- grups.lm$fstatistic[3]
  grups.pval <- pf(fstat, df1, df2, lower.tail = F)
  return(c(
    rsquared = grups.r2,
    fstat = fstat,
    pval = grups.pval
  ))
}

cor_comparison <- function(X1, X2) {
  fisher_z <- function(r) {
    return(0.5 * log((1 + r) / (1 - r)))
  }
  # X1 <- grupo1
  # X2 <- grupo2
  Xtotal <- rbind(X1, X2)
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  
  # Function to compute pairwise differences between correlation matrices using Fisher Z
  # Compute correlation matrices
  R1 <- cor(X1)
  R2 <- cor(X2)
  mtest <- ape::mantel.test(R1, R2)
  names_correlations1 <- colnames(R1) <- rownames(R1) <- colnames(X1)
  names_correlations2 <- colnames(R2) <- rownames(R2) <- colnames(X1)
  
  R1.upper <- R1[upper.tri(R1)]
  R2.upper <- R2[upper.tri(R2)]
  
  names_correlations1 <- outer(colnames(R1), rownames(R1), paste, sep = "-")
  names_correlations1 <- names_correlations1[upper.tri(names_correlations1)]
  
  names_correlations2 <- outer(colnames(R2), rownames(R2), paste, sep = "-")
  names_correlations2 <- names_correlations2[upper.tri(names_correlations2)]
  
  names(R1.upper) <- names_correlations1
  names(R2.upper) <- names_correlations2
  
  sd1 <- (1 / (n1 - 3))
  sd2 <- (1  / (n2 - 3))
  # Compute observed test statistic (sum of Z differences)
  observed_stat <- ((fisher_z(R1.upper) - fisher_z(R2.upper))) / sqrt(sd1 +
                                                                        sd2)
  
  p_value <- 2 * (1 - pnorm(abs(observed_stat)))
  
  
  
  n <- ncol(R1)
  # Crear una matriz vacía
  R_observed <- matrix(0, nrow = n, ncol = n)
  
  # Colocar los elementos del vector en la parte superior triangular
  R_observed[upper.tri(R_observed)] <- observed_stat
  
  
  p_values_matrix_ <- matrix(0, nrow = n, ncol = n)
  
  # Colocar los elementos del vector en la parte superior triangular
  p_values_matrix_[upper.tri(p_values_matrix_)] <- p.adjust(p_value, "BH")
  n <- ncol(R_observed)
  R_symmetric <- R_observed + t(R_observed) + diag(1, n)
  
  p_values_matrix_ <- p_values_matrix_ + t(p_values_matrix_) + diag(1, n)
  
  
  
  
  
  Q_observado = sum(((
    fisher_z(R1.upper) - fisher_z(R2.upper)
  ) ^ 2) - (sd1 + sd2))
  
  
  n_permutations <- 1000
  Q_permutado <- numeric(n_permutations)
  p_valor <- numeric(n_permutations)
  p <- ncol(X1)
  dof <- p
  
  for (i in 1:n_permutations) {
    # Permutar las etiquetas de grupo
    permuted_labels <- sample(c(rep("grupo1", n1), rep("grupo2", n2)))
    grupo1_ <- Xtotal[which(permuted_labels == "grupo1"), ]
    grupo2_ <- Xtotal[which(permuted_labels == "grupo2"), ]
    corObese <- cor(grupo1_)
    corNoObese <- cor(grupo2_)
    R1_i <- fisher_z(corObese[upper.tri(corObese)])
    R2_i <- fisher_z(corNoObese[upper.tri(corObese)])
    # Calcular las correlaciones en varias muestras de bootstrap para ambos conjuntos de dato
    v1 <- (1) / (n1 - 3)
    v2 <- (1) / (n2 - 3)
    Z_statistic_perm <- ((R1_i - R2_i) ^ 2) - (v1 + v2)
    Q_permutado[i] <- (sum(Z_statistic_perm))
    
    p_valor[i] <- pchisq(Q_permutado[i], df = p, lower.tail = F)
    
    
    
    
    # Calcular las matrices de correlación y los estadísticos Z para los datos permutados
    # (Repite los pasos anteriores para obtener Z_statistic_perm)
    
    # Calcular Q para la permutación
    
  }
  (p_tot <- pchisq(Q_observado, df = (p * (p - 1)) / 2, lower.tail = F))
  
  (p_mean <- mean((Q_permutado) >= (Q_observado)))
  (p_median <- median((Q_permutado) >= (Q_observado)))
  hist(Q_permutado)
  chis_stat <- Q_observado
  zA <- fisher_z(R1.upper)
  zB <-  fisher_z(R2.upper)
  n_permutations <- 1000
  dist_mahalanobis.list <- vector("list", length(n_permutations))
  print(paste("permutacion"))
  # Si estás usando bootstrap, puedes calcular el bootstrap de las diferencias
  bootstrap_correlations_diff <- replicate(n_permutations, {
    X_tot <- as.data.frame(rbind(X1, X2))
    X_tot$labels <- c(rep("grupo1", nrow(X1)), rep("grupo2", nrow(X2)))
    X_bootstrap <- X_tot[sample(1:nrow(X_tot), replace = T), ]
    
    
    data_bootstrap_A <- X_bootstrap[X_bootstrap$labels == "grupo1", -ncol(X_bootstrap)]
    data_bootstrap_B <- X_bootstrap[X_bootstrap$labels == "grupo2", -ncol(X_bootstrap)]
    
    cor_A_bootstrap <- cor(data_bootstrap_A)[upper.tri(cor(data_bootstrap_A))]
    cor_B_bootstrap <- cor(data_bootstrap_B)[upper.tri(cor(data_bootstrap_B))]
    
    probA  <- fisher_z(cor_A_bootstrap)
    probB  <- fisher_z(cor_B_bootstrap)
    
    probA - probB
    
    
  })
  bootstrap_correlations_diff <- bootstrap_correlations_diff[, colSums(is.na(bootstrap_correlations_diff)) == 0]
  bootstrap_correlations_diff <- bootstrap_correlations_diff[, colSums(is.infinite(bootstrap_correlations_diff)) == 0]
  print(paste("computacion"))
  
  # cov_matrix_diff <- cov(t(bootstrap_correlations_diff))
  
  cov_matrix_diff <- stats::cov(zA, zB)
  
  pseudoinversa <- MASS::ginv(cov_matrix_diff)
  
  mahalanobis_distance_stat <- (t(zA) %*% (zB)) / pseudoinversa
  
  mahalanobis_pvalue <- pchisq(mahalanobis_distance_stat,
                               df = (p * (p - 1)) / 2,
                               lower.tail = F)
  
  
  
  
  
  return(
    list(
      Q_observado = Q_observado,
      p_val_median = p_median,
      p_val_mean = p_mean,
      p_observado = p_tot,
      q_obeservado = Q_observado,
      chiperm = Q_permutado,
      p_permutado = p_valor,
      mahalanobis_distance_stat = mahalanobis_distance_stat,
      mahalanobis_pvalue = mahalanobis_pvalue,
      p_value_matrix = p_values_matrix_,
      R_symmetric = R_symmetric,
      z_stat = mtest$z.stat,
      p_stat = mtest$p
      
      
    )
  )
  
}
plot_correlation <-  function(R, titulo) {
  # R <- R[,biomarcadores]
  colnames(R) <- colnames(pre)[colnames(pre) %in% colnames(R)]
  # R <- grupo1
  cor_test_result <- psych::corr.test(R)
  
  # Extraer la matriz de correlaciones
  cor_matrix <- as.data.frame(cor_test_result$r)[, rev(1:ncol(R))]  # Matriz de correlaciones
  diag(cor_matrix[, rev(1:ncol(R))]) <- 0
  
  # Extraer el vector de p-valores ajustados (fuera de la diagonal)
  p_adj_vector <- cor_test_result$p.adj  # Vector de p-valores ajustados
  
  # Crear una matriz vacía para los p-valores ajustados
  p_adj_matrix <- matrix(NA,
                         nrow = ncol(cor_matrix),
                         ncol = ncol(cor_matrix))
  
  # Rellenar la matriz vacía con los p-valores ajustados por encima de la diagonal
  p_adj_matrix[upper.tri(p_adj_matrix)] <- p_adj_vector
  
  # Hacer la matriz simétrica (rellenar la parte inferior de la diagonal con los mismos valores)
  p_adj_matrix[lower.tri(p_adj_matrix)] <- t(p_adj_matrix)[lower.tri(p_adj_matrix)]
  
  colnames(p_adj_matrix) <- rownames(p_adj_matrix) <- colnames(R)
  p_adj_matrix <- as.data.frame(p_adj_matrix[, rev(1:ncol(p_adj_matrix))])
  diag(p_adj_matrix) <- 1
  p_adj_matrix$Variable <- colnames(R)
  # Convertir las matrices de correlaciones y p-valores ajustados a formato largo (long format)
  cor_matrix$Variable <- colnames(R)
  cor_matrix$Variable <- factor(cor_matrix$Variable, colnames(R))
  cor_data <- reshape2::melt(cor_matrix, id.vars = "Variable")
  cor_data$Variable <- factor(cor_data$Variable, levels = colnames(R))
  p_data <- melt(p_adj_matrix, id.vars = "Variable")
  cor_data$variable <- factor(cor_data$variable, levels = colnames(R))
  
  # Crear una columna con los asteriscos para las celdas significativas (p.adj < 0.05)
  p_data$signif <- ifelse(p_data$value < 0.05, "*", "")
  
  # Fusionar las matrices de correlación y significancia ajustada
  plot_data <- merge(cor_data, p_data, by = c("Variable", "variable"))
  
  # Renombrar las columnas para mayor claridad
  colnames(plot_data) <- c("Var1", "Var2", "Correlation", "P_value_adj", "Significance")
  
  
  # Fusionar las matrices de correlación y significancia
  # Graficar con geom_tile y marcar las celdas significativas con asteriscos
  p <- ggplot(plot_data, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile() +
    geom_text(aes(label = Significance),
              color = "black",
              size = 2.4) +  # Añadir asteriscos en las celdas significativas
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      space = "Lab",
      name = "Pearson \nCorrelation Coefficient",
      limits = c(-1, 1)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 1,
        hjust = 1,
        size = 10
      ),
      axis.text.y = element_text(size = 10)
    ) +
    labs(x = "", y = "") + ggtitle(titulo)
  
  
  return(p)
  
  
}
extract_features <- function(omic, f, pvals, threshold) {
  prueba <- misuper3(omic)
  omic_vars <- colnames(omic)
  biomarcadores_lista <- fun_summary(prueba)
  selected_biomarkers <-
    lapply(biomarcadores_lista, function(x)
      x[x$pval < pvals & x$pval.value < pvals, ])
  
  varis <-
    unique(unlist(lapply(selected_biomarkers, rownames)))
  pesos_omic <- pesos[omic_vars, f]
  
  pesos_selected <-
    rownames(pesos_omic)[apply(pesos_omic, 2, function(x)
      any(abs(x) > threshold))] ##aqui
  
  
  variables_selected <- intersect(pesos_selected, varis)
  
  res <- lapply(selected_biomarkers, function(x)
    x[rownames(x) %in% variables_selected, ])
  
  return(res)
}
outlier_detection <- function(X, ncomp) {
  pcx <- prcomp(X)
  
  y <- pcx$x[, ncomp]
  mod <- lm(y ~ grupo * obesidad)
  summary(lm(sqrt(abs(residuals(
    mod
  ))) ~ fitted(mod)))
  shapiro.test(residuals(mod))
  hatv <- hatvalues(mod)
  p <- length(mod$coefficients) # k+1
  n <- length(mod$fitted.values)
  leverage.mean <- p / n # (k+1)/n
  leverage_points <- which(hatv > 2 * leverage.mean)
  jack <- rstudent(mod) # jacknife residual
  grlib <- n - p - 1
  outliers <- which(abs(jack) > abs(qt(0.05 / (2 * n), grlib)))
  csa3 <- which(abs(jack) > 2)
  
  return(unique(names(c(
    outliers, leverage_points, csa3
  ))))
}

directorio <- "./resultados/14-10-24_with_outliers"
if (!dir.exists(directorio)) {
  dir.create(directorio, recursive = T)
}

fun_pca_2 <- function(R, biomarcadores, titulo) {
  pcx <- prcomp(R[, biomarcadores])
  
  plotdf <- data.frame(pcx$x,
                       Group = grupo,
                       Obesity = obesidad,
                       sujetos = rownames(R))
  varianzas <- round(100 * (pcx$sdev ^ 2) / sum(pcx$sdev ^ 2), 2)
  interaccion <- interaction(grupo, obesidad)
  levels(interaccion) <- c(
    "Control Women: No Obese",
    "PCOS: No Obese",
    "Men: No Obese",
    "Control Women: Obese",
    "PCOS: Obese",
    "Men: Obese"
  )
  
  plotdf$subjects <- interaccion
  
  
  p4 <- factoextra::fviz_screeplot(pcx) + ylab("% Variance")
  
  
  levels(plotdf$Group) <- c("Control Women", "PCOS", "Men")
  
  
  colores <- ifelse(plotdf$subjects == "Control Women: No Obese", "darkgreen", NA)
  colores <- ifelse(plotdf$subjects == "PCOS: No Obese", "yellow4", colores)
  colores <- ifelse(plotdf$subjects == "Men: No Obese", "red4", colores)
  colores <- ifelse(plotdf$subjects == "Control Women: Obese" ,
                    "lightgreen",
                    colores)
  colores <- ifelse(plotdf$subjects == "PCOS: Obese", "wheat3", colores)
  colores <- ifelse(plotdf$subjects == "Men: Obese", "red", colores)
  
  formas <- ifelse(plotdf$Group == "Control Women", 0, NA)
  formas <- ifelse(plotdf$Group == "PCOS", 1, formas)
  formas <- ifelse(plotdf$Group == "Men", 2, formas)
  
  plotdf$colores <- colores
  plotdf$formas <- formas
  # formas <- c("Control Women" = 0, "PCOS" = 1, "Men" = 2)
  
  # Crear el plot con bordes en negro y relleno de color
  plot1 <- ggplot(plotdf, aes(
    x = PC1,
    y = PC2,
    shape = Group,
    color = subjects,
  )) +
    geom_point(size = 7, aes(shape = Group, color = subjects)) +  # Bordes en negro
    scale_shape_manual(values = c(15, 16, 17)) +  # Definir las formas manualmente
    scale_color_manual(
      values = c(
        "Control Women: No Obese" = "green",
        "PCOS: No Obese" = "yellow",
        "Men: No Obese" = "red",
        "Control Women: Obese" = "seagreen4",
        "PCOS: Obese" = "#EEE685",
        "Men: Obese" = "#8B0000"
      )
    ) + # Definir los colores de relleno manualmente
    xlab(paste0("PC1: ", varianzas[1], "%")) +
    ylab(paste0("PC2: ", varianzas[2], "%"))+theme(text = element_text(size = 20))
  
  # plot1# Gráfico 2 - Solo triángulos (obesidad)
  plot2 <- ggplot(subset(plotdf, Obesity == "Obese"),
                  aes(
                    x = PC1,
                    y = PC2,
                    fill = subjects,
                    shape = Group
                  )) +
    geom_point(size = 7, aes(shape = Group, color = subjects)) +  # Bordes en negro
    scale_shape_manual(values = c(15, 16, 17)) +  # Definir las formas manualmente
    scale_color_manual(
      values = c(
        "Control Women: Obese" = "seagreen4",
        "PCOS: Obese" = "#EEE685",
        "Men: Obese" = "#8B0000"
      )
    ) + # Definir los colores de relleno manualmente
    xlab(paste0("PC1: ", varianzas[1], "%")) +
    ylab(paste0("PC2: ", varianzas[2], "%"))+theme(legend.position = "none",text = element_text(size = 20))
  
  # Gráfico 3 - Solo círculos (no obesidad)
  plot3 <- ggplot(subset(plotdf, Obesity == "No Obese"),
                  aes(
                    x = PC1,
                    y = PC2,
                    fill = subjects,
                    shape = Group
                  )) +
    geom_point(size = 7, aes(shape = Group, color = subjects)) +  # Bordes en negro
    scale_shape_manual(values = c(15, 16, 17)) +  # Definir las formas manualmente
    scale_color_manual(
      values = c(
        "Control Women: No Obese" = "green",
        "PCOS: No Obese" = "yellow",
        "Men: No Obese" = "red"
      )
    ) + # Definir los colores de relleno manualmente
    xlab(paste0("PC1: ", varianzas[1], "%")) +
    ylab(paste0("PC2: ", varianzas[2], "%"))+theme(legend.position = "none",text = element_text(size = 20))
  
  
  
  combined_plot <- (plot2 + plot3) / (plot1 + p4) +
    plot_layout(guides = "collect") +
    plot_annotation(title = titulo)
  combined_plot + theme(legend.key.size = 20,legend.text = 20)
  return(combined_plot)
}

modelos <- lapply(list.files("./modelos/14_10_24_With_OUTLIERS/", full.names = T),
                  load_model)
modelo <- MOFA2::select_model(modelos, plot = T)


p1 <- plot_variance_explained(modelo, max_r2 = 15) + theme(axis.text.x =
                                                             element_text(angle = -45, hjust = 0)) + scale_x_discrete(labels = c(
                                                               "metaboloma" = "Metabolome",
                                                               "metagenoma" = "Microbiome",
                                                               "ip" = "Proteins"
                                                             ))

p2 <- plot_variance_explained(modelo, plot_total = T)[[2]] + theme(axis.text.x =
                                                                     element_text(angle = -45, hjust = 0)) + scale_x_discrete(labels = c(
                                                                       "metaboloma" = "Metabolome",
                                                                       "metagenoma" = "Microbiome",
                                                                       "ip" = "Proteins"
                                                                     ))

p3 <- ggarrange(p1, p2)

p3

directorio_preliminar <- file.path(directorio, "preliminares")
if (!dir.exists(directorio_preliminar)) {
  dir.create(directorio_preliminar, recursive = T)
}
ggsave(
  filename = file.path(directorio_preliminar, "varianza_explicada.jpeg"),
  plot = p3
)

###====factores latentes=====

####====factores univariante=====



factores.df <- as.data.frame(MOFA2::get_expectations(modelo,"Z",as.data.frame = F)[[1]])



res_factores <- fun_plot_limma(factores.df)



p_ <- ggplot(res_factores$logFC,
             aes(Variable, variable, fill = value, label = text)) + geom_tile() + scale_fill_gradient2(high = "red",
                                                                                                       low = "blue",
                                                                                                       mid = "white") + theme(axis.text.x = element_text(
                                                                                                         angle = 90,
                                                                                                         vjust = 0.5,
                                                                                                         hjust = 1
                                                                                                       )) + geom_text(size = 7) + theme(text = element_text(size = 20, face = "bold"))  + xlab("") +
  ylab("") + ggtitle("logFC")
p_
ordern_contrastes <- c(
  "Obese vs Non-obese",
  "Men vs Control Women",
  "PCOS vs Control Women",
  "PCOS vs Men",
  "Control Women (Obese vs Non-obese)",
  "PCOS (Obese vs Non-obese)",
  "Men (Obese vs Non-obese)",
  "Non-obese (Men vs Control Women)",
  "Non-obese (PCOS vs Control Women)",
  "Non-obese (PCOS vs Men)",
  "Obese (Men vs Control Women)",
  "Obese (PCOS vs Control Women)",
  "Obese (PCOS vs Men)"
)

p_ <- p_ + scale_x_discrete(labels = c("Factor 1", "Factor 2", "Factor 3", "Factor 4", "Factor 5")) + scale_y_discrete(labels =
                                                                                                                         ordern_contrastes)
p_

directorio_preliminar <- file.path(directorio, "LatentFactors")
if (!dir.exists(directorio_preliminar)) {
  dir.create(directorio_preliminar, recursive = T)
}
ggsave(
  filename = file.path(directorio_preliminar, "LatentFactors_univariante.jpeg"),
  plot = p_,
  height = 15,
  width = 14
)

####====latent correlated clinical=====


feat_fuera <- c("group", "sample")
idx <- sapply(feat_fuera, function(x)
  grep(x, colnames(modelo@samples_metadata)))

## reorder clinical variables

matrix1 <- modelo@samples_metadata[, -idx]
matrix1 <- matrix1
matrix2 <- as.matrix(factores.df)

corr_result <- corr.test(matrix1, matrix2, adjust = "none", method = "pearson")

# Extract correlations and p-values
correlation_matrix <- t(corr_result$r)     # Get the correlations for matching columns

clinical <- c(
  "BMI",
  "Waist.circumference",
  "Waist.to.Hip.Ratio",
  "Total.Testosterone",
  "Free.Testosterone",
  "Total.Estradiol",
  "Free.Estradiol",
  "SHBG",
  "Glucose",
  "Insulin",
  "HOMA.IR",
  "ISI",
  "Triglycerides",
  "Cholesterol",
  "HDL.Cholesterol",
  "LDL.Cholesterol"
)
correlation_matrix <- correlation_matrix[,clinical]


p_value_matrix <- t(corr_result$p.adj)[,clinical]
# colnames(p_value_matrix) <- clinical
#Get the p-values for matching columns

# Create a data frame for ggplot2
data <- reshape2::melt(correlation_matrix)
data$pvals <- reshape2::melt(p_value_matrix)$value

# Add a column to indicate significance based on p-value < 0.05
data$significant <- ifelse(data$pvals < 0.05, "*", "")

# Create the heatmap
p <- ggplot(data, aes(x = Var1, y = Var2, fill = value)) +  # Heatmap with single row
  geom_tile() +
  geom_text(aes(label = significant),
            color = "black",
            size = 7) +  # Add stars for significance
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1)
  ) +  # Color gradient
  labs(x = "", y = "", fill = "Pearson Correlation\nCoefficient") +
  theme_minimal() + scale_x_discrete(labels = c("Factor 1", "Factor 2", "Factor 3", "Factor 4", "Factor 5")) + theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    text = element_text(size = 20, face = "bold")
  ) + ggtitle("Latent Variables correlated with Covariables")
p
ggsave(
  filename = file.path(
    directorio_preliminar,
    "LatentFactors_correlated_clinical.jpeg"
  ),
  plot = p,
  height = 18,
  width = 15
)
p
##====Feature selction #####
pesos<- Reduce("rbind",get_expectations(modelo,"W"))
factores <- get_expectations(modelo,variable = "Z")[[1]]
R <- factores %*% t(pesos)

### feature selection


p1 <- plot_top_weights(modelo, factors = 1:2, view = c(1))

# ggsave(filename=file.path(directorio_weights_explained_variance,"top_weights_metabo.jpeg"),plot = p1)

p1

p1 <- plot_top_weights(modelo, factors = 2:3, view = 3)

# ggsave(filename=file.path(directorio_weights_explained_variance,"top_wights_proteo.jpeg"),plot = p1)

p1

p1 <- plot_top_weights(modelo, factors = 3:5, view = 2)

# ggsave(filename=file.path(directorio_weights_explained_variance,"top_weights_micro.jpeg"),plot = p1)

p1
### overall threshold

#
metaboloma.r <- factores[, 1] %*% t(pesos[rownames(modelo@data$metaboloma$group1), 1])
proteoma.r <- factores[, 2:3] %*% t(pesos[rownames(modelo@data$ip$group1), 2:3])
micro.r <- factores[, 3:5] %*% t(pesos[rownames(modelo@data$metagenoma$group1), 3:5])

pvals <- 0.05
threshold <- 0.8
res.met <- extract_features(metaboloma.r, 1, pvals = pvals, threshold =
                              threshold)
res.mic <- extract_features(micro.r, 3:5, pvals = pvals, threshold = threshold)
res.prot <- extract_features(proteoma.r, 2:3, pvals = pvals, threshold =
                               threshold)

#

lista.res <- list(metaboloma = res.met,
                  microbioma = res.mic,
                  proteoma = res.prot)

vector.res <- vector("list", length = length(res.met))

for (l in 1:length(lista.res[[1]])) {
  vector.res[[l]] <- rbind(lista.res[[1]][[l]], lista.res[[2]][[l]], lista.res[[3]][[l]])
  
}

names(vector.res) <- names(res.met)
biomarcadores <- unique(unlist(lapply(vector.res, rownames)))
# R <- (cbind(metaboloma.r, micro.r, proteoma.r))
R  <- factores %*% t(pesos)

Rprueba <- R[, biomarcadores]
pre <-
  bind_cols(as.data.frame(scale(t(
    Reduce(
      "rbind",
      list(
        modelo@data$metaboloma[[1]],
        modelo@data$metagenoma[[1]],
        modelo@data$ip[[1]]
      )
    )
  ))))

biomarcadores <- biomarcadores[biomarcadores %in% colnames(pre)]

preprueba <- pre[, biomarcadores] ### matriz original

R_selected <- R[, biomarcadores] ### matriz reconstruida

##========More latent factors=====


pesos <- as.data.frame(Reduce(rbind, get_weights(modelo, scale = T)))
pesos_biomarcadores <- pesos[biomarcadores, ]
pesos_biomarcadores$biomarcador <- rownames(pesos_biomarcadores)
pesos_biomarcadores$biomarcador <- as.factor(pesos_biomarcadores$biomarcador)
pesos_biomarcadores$omic <- ifelse(rownames(pesos_biomarcadores) %in% rownames(modelo@data$metaboloma$group1),"Metabolme",NA)
pesos_biomarcadores$omic <- ifelse(rownames(pesos_biomarcadores) %in% rownames(modelo@data$metagenoma$group1),"Microbiome",pesos_biomarcadores$omic)
pesos_biomarcadores$omic <- ifelse(rownames(pesos_biomarcadores) %in% rownames(modelo@data$ip$group1),"Proteins",pesos_biomarcadores$omic)
pesos_biomarcadores$omic <- as.factor(pesos_biomarcadores$omic)
# Plot
factor1 <- pesos_biomarcadores %>% mutate(biomarcador = fct_reorder(biomarcador, Factor1)) %>% ggplot(aes(x =
                                                                                                            biomarcador, y = Factor1)) +
  geom_segment(aes(
    x = biomarcador,
    xend = biomarcador,
    y = 0,
    yend = Factor1,colour = "black"
  )) +
  geom_point(size = 5,colour = "black") +
  coord_flip() +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Weight") +
  ggtitle("Factor 1") + theme(text = element_text(size = 20)) + ylim(c(-1, 1))
# Combinar la gráfica original con la columna de símbolos

# Mostrar la gráfica final


# Plot
factor2 <- pesos_biomarcadores %>% mutate(biomarcador = fct_reorder(biomarcador, Factor2)) %>% ggplot(aes(x =
                                                                                                            biomarcador, y = Factor2)) +
  geom_segment(aes(
    x = biomarcador,
    xend = biomarcador,
    y = 0,
    yend = Factor2,colour = "black"
  )) +
  geom_point( size = 5) +
  coord_flip() +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Weight") +
  ggtitle("Factor 2") + theme(text = element_text(size = 20)) + ylim(c(-1, 1))
# Combinar la gráfica original con la columna de símbolos




# Plot
factor3 <- pesos_biomarcadores %>% mutate(biomarcador = fct_reorder(biomarcador, Factor3)) %>% ggplot(aes(x =
                                                                                                            biomarcador, y = Factor3)) +
  geom_segment(aes(
    x = biomarcador,
    xend = biomarcador,
    y = 0,
    yend = Factor3,colour = "black"
  )) +
  geom_point(color = "black", size = 5) +
  coord_flip() +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Weight") +
  ggtitle("Factor 3") + theme(text = element_text(size = 20)) + ylim(c(-1, 1))
# Combinar la gr



# Plot
factor4 <- pesos_biomarcadores %>% mutate(biomarcador = fct_reorder(biomarcador, Factor4)) %>% ggplot(aes(x =
                                                                                                            biomarcador, y = Factor4)) +
  geom_segment(aes(
    x = biomarcador,
    xend = biomarcador,
    y = 0,
    yend = Factor4,colour = "black"
  )) +
  geom_point(color = "black", size = 5) +
  coord_flip() +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Weight") +
  ggtitle("Factor 4") + theme(text = element_text(size = 20)) + ylim(c(-1, 1))
# Combinar la gr


# Plot
factor5 <- pesos_biomarcadores %>% mutate(biomarcador = fct_reorder(biomarcador, Factor5)) %>% ggplot(aes(x =
                                                                                                            biomarcador, y = Factor5)) +
  geom_segment(aes(
    x = biomarcador,
    xend = biomarcador,
    y = 0,
    yend = Factor5,colour ="black"
  )) +
  geom_point(color = "black", size = 5) +
  coord_flip() +
  theme(legend.position = "none") +
  xlab("") +
  ylab("Weight") +
  ggtitle("Factor 5") + theme(text = element_text(size = 20)) + ylim(c(-1, 1))
# Combinar la gr
directorio_preliminar <- file.path(directorio, "LatentFactors")
if (!dir.exists(directorio_preliminar)) {
  dir.create(directorio_preliminar, recursive = T)
}
ggsave(
  filename = file.path(directorio_preliminar, "Weights_Factor1.jpeg"),
  plot = factor1,
  height = 18,
  width = 15
)


ggsave(
  filename = file.path(directorio_preliminar, "Weights_Factor2.jpeg"),
  plot = factor2,
  height = 18,
  width = 15
)

ggsave(
  filename = file.path(directorio_preliminar, "Weights_Factor3.jpeg"),
  plot = factor3,
  height = 18,
  width = 15
)

ggsave(
  filename = file.path(directorio_preliminar, "Weights_Factor4.jpeg"),
  plot = factor4,
  height = 18,
  width = 15
)

ggsave(
  filename = file.path(directorio_preliminar, "Weights_Factor5.jpeg"),
  plot = factor5,
  height = 18,
  width = 15
)




ordern_contrastes <- c(
  "Obese vs Non-obese",
  "Men vs Control Women",
  "PCOS vs Control Women",
  "PCOS vs Men",
  "Control Women (Obese vs Non-obese)",
  "PCOS (Obese vs Non-obese)",
  "Men (Obese vs Non-obese)",
  "Non-obese (Men vs Control Women)",
  "Non-obese (PCOS vs Control Women)",
  "Non-obese (PCOS vs Men)",
  "Obese (Men vs Control Women)",
  "Obese (PCOS vs Control Women)",
  "Obese (PCOS vs Men)"
)

###======Univariante reconstruccion====#####
R_biomarcadores <- R[, colnames(pre)[colnames(pre) %in% biomarcadores]]
pre_biomarcadores <- pre[, colnames(pre)[colnames(pre) %in% biomarcadores]]

res_R <- fun_plot_limma(R_biomarcadores)
# res_R$logFC
tstats <- res_R$logFC
p2 <- ggplot(tstats, aes(Variable, variable, fill = value, label = text)) + geom_tile() + scale_fill_gradient2(
  high = "red",
  low = "blue",
  mid = "white",
  limits = c(min(res_R$logFC$value), max(res_R$logFC$value))
) + theme(
  axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1,
    size = 15,
    face = "bold"
  ),
  axis.text.y = element_text(size = 15, face = "bold")
) + geom_text(size = 7)   + xlab("") +
  ylab("") + ggtitle("logFC")


p2 <- p2 + scale_y_discrete(labels =
                              ordern_contrastes)
directorio_preliminar <- file.path(directorio, "Reconstruccion")
if (!dir.exists(directorio_preliminar)) {
  dir.create(directorio_preliminar, recursive = T)
}
directorio_preliminar <- "./resultados/14-10-24_with_outliers/"

ggsave(
  filename = file.path(directorio_preliminar, "Reconstruccion", "Univariante.jpeg"),
  plot = p2,
  width = 18,
  height = 15
)



###=========Bivariante reconstruccion=====

R_obesidad <- R_biomarcadores[obesidad == "Obese", ]
R_NoObesidad <- R_biomarcadores[obesidad == "No Obese", ]
R_females <- R_biomarcadores[grupo == "Female", ]
R_males <- R_biomarcadores[grupo == "Male", ]
R_PCOS <- R_biomarcadores[grupo == "PCOS", ]
R_malenoob <- R_biomarcadores[grupo == "Male" &
                                obesidad == "Obese", ]
R_maleoob <- R_biomarcadores[grupo == "Male" &
                               obesidad == "No Obese", ]
R_PCOSob <- R_biomarcadores[grupo == "PCOS" & obesidad == "Obese", ]
R_PCOSnoob <- R_biomarcadores[grupo == "PCOS" &
                                obesidad == "No Obese", ]
R_femaleOb <- R_biomarcadores[grupo == "Female" &
                                obesidad == "Obese", ]
R_femaleNoob <- R_biomarcadores[grupo == "Female" &
                                  obesidad == "No Obese", ]


grupos_R <- list(
  noob = R_NoObesidad,
  ob = R_obesidad,
  R_females = R_females,
  R_males = R_males,
  R_PCOS = R_PCOS,
  R_malenoob = R_malenoob,
  R_maleoob = R_maleoob,
  R_PCOSob = R_PCOSob,
  R_PCOSnoob = R_PCOSnoob,
  R_femaleOb = R_femaleOb,
  R_femaleNoob = R_femaleNoob,
  all = R[, biomarcadores]
)




nombres <- c(
  "Obese",
  "No Obese",
  "Female",
  "Male",
  "PCOS",
  "Male: No Obese",
  "Male: Obese",
  "PCOS: Obese",
  "PCOS: No Obese",
  "Female: Obese",
  "Female: No Obese",
  "All subjects"
)

plots <-  vector("list", length = length(nombres))
for (gof in 1:length(nombres)) {
  plots[[gof]] <- plot_correlation(grupos_R[[gof]], titulo = nombres[gof])
  
  
}


plots



directorio_preliminar <- file.path(directorio, "Reconstruccion")
if (!dir.exists(directorio_preliminar)) {
  dir.create(directorio_preliminar, recursive = T)
}
ggsave(filename = file.path(directorio_preliminar, "Weights.jpeg"),
       plot = factor1)

for (p in 1:length(plots)) {
  ggsave(
    plot = plots[[p]],
    filename = file.path(
      directorio_preliminar,
      paste0("Correlated_", nombres[p], ".jpeg")
    ),
    height = 10,
    width = 10
  )
  
  
  
}
##======PCA=======
R <- as.matrix(get_expectations(modelo,"Z",as.data.frame = F)[[1]]) %*% t(Reduce("rbind",get_expectations(modelo,"W",as.data.frame=F)))
Original <- fun_pca_2(pre, biomarcadores, "PCA: Original Signal")
reconstructed <- fun_pca_2(R, biomarcadores, "PCA: Reconstructed Signal")

directorio_preliminar <- file.path(directorio, "Reconstruccion", "PCA")
if (!dir.exists(directorio_preliminar)) {
  dir.create(directorio_preliminar, recursive = T)
}
ggsave(
  filename = file.path(directorio_preliminar, "Reconstructed.jpeg"),
  plot = reconstructed,width = 18,height = 18
)


directorio_preliminar <- file.path(directorio, "Original", "PCA")
if (!dir.exists(directorio_preliminar)) {
  dir.create(directorio_preliminar, recursive = T)
}
ggsave(filename = file.path(directorio_preliminar, "Original.jpeg"),
       plot = Original,width = 18,height = 18)


####====Downstream PCA=====
#####====PC1======
######====PC1 grupo======
pcx <- prcomp(R[, biomarcadores], center = T, scale. = T)
interaccion <- interaction(grupo, obesidad)
levels(interaccion) <- c(
  "Control Women: No Obese",
  "PCOS: No Obese",
  "Men: No Obese",
  "Control Women: Obese",
  "PCOS: Obese",
  "Men: Obese"
)
grupo_ <- grupo
levels(grupo_) <- c("Control Women","PCOS","Men")
PCs <- data.frame(pcx$x, Group = grupo_, Obesity = obesidad,Subjects = interaccion)

anova.grupo <- aov(PCs$PC1 ~ 0 + grupo * obesidad)
anova.grupo.sum <- summary(anova.grupo)
titulo <- paste(
  "Group x Obesity ANOVA: F (Group) = ",
  round(anova.grupo.sum[[1]]$`F value`[1], 4),
  "\nP-val: = ",
  round(anova.grupo.sum[[1]]$`Pr(>F)`[1], 4)
)
colores <- c("darkgreen","yellow4","red4","lightgreen","wheat3","red")

p <- ggboxplot(
  PCs,
  x = "Group",
  y = "PC1",
  color = "Subjects",
  add = c("jitter", "mean"),facet.by = "Obesity",palette = colores
)
p
my_comparisons <- list(c("Control Women", "Men"), c("Control Women", "PCOS"), c("Men", "PCOS"))

p <- p + stat_compare_means(comparisons = my_comparisons,
                            method = "t.test",
                            hide.ns = T,ref.group = "Control Women",paired = F) + theme(legend.position = "none") +
  ggtitle(titulo) + xlab("") + coord_cartesian(ylim = c(-20, 20))

p

directorio_preliminar <- file.path(directorio, "Reconstruccion", "PCA")
if (!dir.exists(directorio_preliminar)) {
  dir.create(directorio_preliminar, recursive = T)
}
ggsave(
  filename = file.path(directorio_preliminar, "PC1_F_Group.jpeg"),
  plot = p
)





######====PC1 obesidad======

titulo <- paste(
  "Group x Obesity ANOVA: F (Obesity) = ",
  round(anova.grupo.sum[[1]]$`F value`[2], 4),
  "\nP-val: = ",
  round(anova.grupo.sum[[1]]$`Pr(>F)`[2], 4)
)

p <- ggboxplot(
  PCs,
  x = "Obesity",
  y = "PC1",
  color = "Subjects",
  add = c("jitter", "mean"),facet.by = "Group",palette = colores
)
p
my_comparisons <- list(c("No Obese","Obese"))
p <- p + xlab("") + facet_grid(~ Group) + stat_compare_means(comparisons = my_comparisons,
                                                             method = "t.test",
                                                             hide.ns = T)  + ggtitle(titulo)+ coord_cartesian(ylim = c(min(PCs$PC1),2*max(PCs$PC1)))+theme(legend.position = "none")

p

ggsave(filename = file.path(directorio_preliminar, "PC1_F_Obesity.jpeg"),
       plot = p)


#####===PC2=====

######====PC2 grupo======

anova.grupo <- aov(PCs$PC2 ~ 0 + grupo * obesidad)
anova.grupo.sum <- summary(anova.grupo)
titulo <- paste(
  "Group x Obesity ANOVA: F (Group) = ",
  round(anova.grupo.sum[[1]]$`F value`[1], 4),
  "\nP-val: = ",
  round(anova.grupo.sum[[1]]$`Pr(>F)`[1], 4)
)
colores <- c("darkgreen","yellow4","red4","lightgreen","wheat3","red")

p <- ggboxplot(
  PCs,
  x = "Group",
  y = "PC2",
  color = "Subjects",
  add = c("jitter", "mean"),facet.by = "Obesity",palette = colores
)
p
my_comparisons <- list(c("Control Women", "Men"), c("Control Women", "PCOS"), c("Men", "PCOS"))

p <- p + stat_compare_means(comparisons = my_comparisons,
                            method = "t.test",
                            hide.ns = T,ref.group = "Control Women",paired = F) + theme(legend.position = "none") +
  ggtitle(titulo) + xlab("") + coord_cartesian(ylim = c(-13, 20))

p

directorio_preliminar <- file.path(directorio, "Reconstruccion", "PCA")
if (!dir.exists(directorio_preliminar)) {
  dir.create(directorio_preliminar, recursive = T)
}
ggsave(
  filename = file.path(directorio_preliminar, "PC2_F_Group.jpeg"),
  plot = p
)





######====PC2 obesidad======

titulo <- paste(
  "Group x Obesity ANOVA: F (Obesity) = ",
  round(anova.grupo.sum[[1]]$`F value`[2], 4),
  "\nP-val: = ",
  round(anova.grupo.sum[[1]]$`Pr(>F)`[2], 4)
)

p <- ggboxplot(
  PCs,
  x = "Obesity",
  y = "PC2",
  color = "Subjects",
  add = c("jitter", "mean"),facet.by = "Group",palette = colores
)
p
my_comparisons <- list(c("No Obese","Obese"))
p <- p + xlab("") + facet_grid(~ Group) + stat_compare_means(comparisons = my_comparisons,
                                                             method = "t.test",
                                                             hide.ns = T)  + ggtitle(titulo)+ coord_cartesian(ylim = c(min(PCs$PC2),2*max(PCs$PC2)))+theme(legend.position = "none")

p

ggsave(filename = file.path(directorio_preliminar, "PC2_F_Obesity.jpeg"),
       plot = p)

####===Contribution plots ======


# SE DEBE DE REALIZAR EL ANALISIS CON LA MATRIZ CENTRADA Y ESCALADA LUEGO SE COMPUTAN LOS LOADIGNS CON LA CORRELACION DE LA MATRIZ


### Correlacionamos los biomarcadores reconstruidos, con los scores.



fun_contribution <- function(X,
                             scores,
                             unfactor,
                             referencia,
                             Titulo,
                             componente,colores) {
  R_grupo <- X
  ctrl <- unfactor
  # componente <- 1
  robj <- apply((R_grupo), 2, function(x)
    corr.test(x, scores[, componente]))
  r <- unlist(lapply(robj, function(x)
    x$r))
  # p <- p.adjust(unlist(lapply(robj,function(x) x$p)),"BH")
  # r_thresh <- min(abs(r[which(p<0.01)]))
  toplot <- data.frame(r = r, names = names(r))
  toplot$omic <- ifelse(toplot$names %in% rownames(modelo@data$metaboloma$group1),
                        "Metabolome",
                        NA)
  toplot$omic <- ifelse(
    toplot$names %in% rownames(modelo@data$metagenoma$group1),
    "Microbiome",
    toplot$omic
  )
  toplot$omic <- ifelse(toplot$names %in% rownames(modelo@data$ip$group1),
                        "Proteins",
                        toplot$omic)
  toplot$omic <- as.factor(toplot$omic)
  
  p_contrib <- ggbarplot(
    toplot,
    x = "names",
    y = "r",
    fill = "omic",
    # change fill color by mpg_level
    color = "white",
    # Set bar border colors to white
    palette = "jco",
    # jco journal color palett. see ?ggpar
    sort.val = "desc",
    # Sort the value in descending order
    sort.by.groups = FALSE,
    # Don't sort inside each group
    x.text.angle = 90,
    # Rotate vertically x axis texts
    ylab = "",
    legend.title = "Omic",
    rotate = TRUE,
    ggtheme = theme_minimal(base_size = 20)
  ) + xlab("") + ylab("Pearson Correlation Coefficient")
  
  
  dfgg <- data.frame(scores , GOF = unfactor)
  # componente<-"PC1"
  # comparisons <- list(c("Female","Male"))
  
  componente_ <- paste0("PC", componente)
  p_score <- ggboxplot(
    dfgg,
    x = "GOF",
    y = componente_,
    fill = "GOF",
    palette = colores,ggtheme = theme(text  = element_text(size = 20))
  ) + stat_compare_means(method = "t.test", ref.group = referencia,size=10) +
    xlab("") + ggtitle(Titulo) + theme(legend.position = "none",axis.text.x = element_text(size=20)) + ylim(c(-max(abs(dfgg[, componente])) -
                                                                            2, max(abs(dfgg[, componente])) + 2))
  
  
  p_complete <- ggarrange(p_contrib, p_score)
  
  return(p_complete)
}
### evitamos con el cambio de variables, modificar las originales

grupo_ <- grupo
levels(grupo_) <- c("Control Women","PCOS","Men")
obesidad_ <- obesidad
sexo_ <- sexo
ctrl <- droplevels(grupo_[grupo_ == "Control Women" | grupo_ == "Men"])
females <- droplevels(grupo_[grupo_ == "Control Women" | grupo_ == "PCOS"])
hermas <- droplevels(grupo_[grupo_ == "Men" | grupo_ == "PCOS"])
R_ <- R
R_grupo <- R_[grupo_ %in% ctrl, ]
comparaciones <- list(control = c("Men", "Control Women"))
referencia <- "Control Women"
######====Obesidad=======

p_contrib_obesidad_PC1 <- fun_contribution(
  X = R[, biomarcadores],
  scores = pcx$x,
  unfactor = obesidad,
  referencia = "No Obese",
  Titulo = "Score Plot PC1:\nObese vs No Obese",
  componente = 1,colores = c("blue","lightblue")
)

p_contrib_obesidad_PC2 <- fun_contribution(
  X = R[, biomarcadores],
  scores = pcx$x,
  unfactor = obesidad,
  referencia = "No Obese",
  Titulo = "Score Plot PC2:\nObese vs No Obese",
  componente = 2,colores = c("blue","lightblue")
)
#####=========Female Male======

p_contrib_Female_vs_Male_PC1 <- fun_contribution(
  X = R[grupo != "PCOS", biomarcadores],
  scores = pcx$x[grupo != "PCOS", ],
  unfactor = ctrl,
  referencia = "Control Women",
  Titulo = "Score Plot PC1:\nMen vs Control Women",
  componente = 1,colores = c("blue","lightblue")
)

p_contrib_Female_vs_Male_PC2 <- fun_contribution(
  X = R[grupo != "PCOS", biomarcadores],
  scores = pcx$x[grupo != "PCOS", ],
  unfactor = ctrl,
  referencia = "Control Women",
  Titulo = "Score Plot PC2:\nMen vs Control Women",
  componente = 2,colores = c("blue","lightblue")
)

#####=========Female PCOS======

p_contrib_PCOS_vs_Female_PC1 <- fun_contribution(
  X = R[grupo_ != "Men", biomarcadores],
  scores = pcx$x[grupo_ != "Men", ],
  unfactor = females,
  referencia = "Control Women",
  Titulo = "Score Plot PC1:\nPCOS vs Control Women",
  componente = 1,colores = c("blue","lightblue")
)



p_contrib_PCOS_vs_Female_PC2 <- fun_contribution(
  X = R[grupo_ != "Men", biomarcadores],
  scores = pcx$x[grupo_ != "Men", ],
  unfactor = females,
  referencia = "Control Women",
  Titulo = "Score Plot PC2: PCOS vs Control Women",
  componente = 2,colores = c("blue","lightblue")
)

#####=========Male PCOS======

p_contrib_PCOS_vs_Male_PC1 <- fun_contribution(
  X = R[grupo_ != "Control Women", biomarcadores],
  scores = pcx$x[grupo_ != "Control Women", ],
  unfactor = hermas,
  referencia = "Men",
  Titulo = "Score Plot PC1:\nPCOS vs Men",
  componente = 1,colores = c("blue","lightblue")
)



p_contrib_PCOS_vs_Male_PC2 <- fun_contribution(
  X = R[grupo_ != "Control Women", biomarcadores],
  scores = pcx$x[grupo_ != "Control Women", ],
  unfactor = hermas,
  referencia = "Men",
  Titulo = "Score Plot PC2:\nPCOS vs Men",
  componente = 2,colores = c("blue","lightblue")
)





#####==========No Obese Female Male======


w_ctrl <- (grupo_ != "PCOS") & (obesidad == "No Obese")
w_females <- (grupo_ != "Men") & (obesidad == "No Obese")
w_hermas <- (grupo_ != "Control Women") & (obesidad == "No Obese")

ctrl <- droplevels(grupo_[w_ctrl])
females <-  droplevels(grupo_[w_females])
hermas <-  droplevels(grupo_[w_hermas])


p_contrib_NOOB_Males_vs_Females_PC1 <- fun_contribution(
  X = R[w_ctrl, biomarcadores],
  scores = pcx$x[w_ctrl, ],
  unfactor = ctrl,
  referencia = "Control Women",
  Titulo = "Score Plot PC1:\n(No Obese) Men vs Control Women)",
  componente = 1,colores = c("blue","lightblue")
)

p_contrib_NOOB_Males_vs_Females_PC2 <- fun_contribution(
  X = R[w_ctrl, biomarcadores],
  scores = pcx$x[w_ctrl, ],
  unfactor = ctrl,
  referencia = "Control Women",
  Titulo = "Score Plot PC2:\n(No Obese) Males vs Control Women",
  componente = 2,colores = c("blue","lightblue")
)

#####=========No Obese Female PCOS======

p_contrib_NOOB_PCOS_vs_Females_PC1 <- fun_contribution(
  X = R[w_females, biomarcadores],
  scores = pcx$x[w_females, ],
  unfactor = females,
  referencia = "Control Women",
  Titulo = "Score Plot PC1:\n(No Obese) PCOS vs Control Women",
  componente = 1,colores = c("blue","lightblue")
)



p_contrib_NOOB_PCOS_vs_Females_PC2 <- fun_contribution(
  X = R[w_females, biomarcadores],
  scores = pcx$x[w_females, ],
  unfactor = females,
  referencia = "Control Women",
  Titulo = "Score Plot PC2: (No Obese) PCOS vs Control Women",
  componente = 2,colores = c("blue","lightblue")
)

#####=========No Obese Male PCOS======


p_contrib_NOOB_PCOS_vs_Males_PC1 <- fun_contribution(
  X = R[w_hermas, biomarcadores],
  scores = pcx$x[w_hermas, ],
  unfactor = hermas,
  referencia = "Men",
  Titulo = "Score Plot PC1:\n(No Obese) PCOS vs Men",
  componente = 1,colores = c("blue","lightblue")
)



p_contrib_NOOB_PCOS_vs_Males_PC2 <- fun_contribution(
  X = R[w_hermas, biomarcadores],
  scores = pcx$x[w_hermas, ],
  unfactor = hermas,
  referencia = "Men",
  Titulo = "Score Plot PC2:\n(No Obese) PCOS vs Men",
  componente = 2,colores = c("blue","lightblue")
)


#####==========Obese Female Male======


w_ctrl <- (grupo_ != "PCOS") & (obesidad == "Obese")
w_females <- (grupo_ != "Men") & (obesidad == "Obese")
w_hermas <- (grupo_ != "Control Women") & (obesidad == "Obese")

ctrl <- droplevels(grupo_[w_ctrl])
females <-  droplevels(grupo_[w_females])
hermas <-  droplevels(grupo_[w_hermas])


p_contrib_OB_Males_vs_Females_PC1 <- fun_contribution(
  X = R[w_ctrl, biomarcadores],
  scores = pcx$x[w_ctrl, ],
  unfactor = ctrl,
  referencia = "Control Women",
  Titulo = "Score Plot PC1:\n(Obese) Males vs Control Women)",
  componente = 1,colores = c("blue","lightblue"))

p_contrib_OB_Males_vs_Females_PC2 <- fun_contribution(
  X = R[w_ctrl, biomarcadores],
  scores = pcx$x[w_ctrl, ],
  unfactor = ctrl,
  referencia = "Control Women",
  Titulo = "Score Plot PC2:\n(Obese) Males vs Control Women",
  componente = 2,colores = c("blue","lightblue")
)

#####=========Obese Female PCOS======

p_contrib_OB_PCOS_vs_Females_PC1 <- fun_contribution(
  X = R[w_females, biomarcadores],
  scores = pcx$x[w_females, ],
  unfactor = females,
  referencia = "Control Women",
  Titulo = "Score Plot PC1:\n(Obese) PCOS vs Control Women",
  componente = 1,colores = c("blue","lightblue")
)



p_contrib_OB_PCOS_vs_Females_PC2 <- fun_contribution(
  X = R[w_females, biomarcadores],
  scores = pcx$x[w_females, ],
  unfactor = females,
  referencia = "Control Women",
  Titulo = "Score Plot PC2:\n(Obese) PCOS vs Control Women",
  componente = 2,colores = c("blue","lightblue")
)

#####=========No Obese Male PCOS======

p_contrib_OB_PCOS_vs_Males_PC1 <- fun_contribution(
  X = R[w_hermas, biomarcadores],
  scores = pcx$x[w_hermas, ],
  unfactor = hermas,
  referencia = "Men",
  Titulo = "Score Plot PC1:\n(Obese) PCOS vs Men",
  componente = 1,colores = c("blue","lightblue")
)



p_contrib_OB_PCOS_vs_Males_PC2 <- fun_contribution(
  X = R[w_hermas, biomarcadores],
  scores = pcx$x[w_hermas, ],
  unfactor = hermas,
  referencia = "Men",
  Titulo = "Score Plot PC2:\n(No Obese) PCOS vs Men",
  componente = 2,colores = c("blue","lightblue")
)


#####===Females=====

w_males <- which(grupo_ == "Men")
w_females <- which(grupo_ == "Control Women")
w_pcos <- which(grupo_ == "PCOS")

males <- (obesidad[w_males])
females <-  (obesidad[w_females])
pcos <-  (obesidad[w_pcos])



#####=======Males ====

p_contrib_MALE_NOOB_vs_OB_PC1 <- fun_contribution(
  X = R[w_males, biomarcadores],
  scores = pcx$x[w_males, ],
  unfactor = males,
  referencia = "No Obese",
  Titulo = "Score Plot PC1:\n(Men) Obese vs No Obese",
  componente = 1,colores = c("blue","lightblue")
)
p_contrib_MALE_NOOB_vs_OB_PC1
p_contrib_MALE_NOOB_vs_OB_PC2 <- fun_contribution(
  X = R[w_males, biomarcadores],
  scores = pcx$x[w_males, ],
  unfactor = males,
  referencia = "No Obese",
  Titulo = "Score Plot PC2:\n(Men) Obese vs No Obese",
  componente = 2,colores = c("blue","lightblue")
)
p_contrib_MALE_NOOB_vs_OB_PC2
#####=====PCOS=======
p_contrib_PCOS_NOOB_vs_OB_PC1 <- fun_contribution(
  X = R[w_pcos, biomarcadores],
  scores = pcx$x[w_pcos, ],
  unfactor = pcos,
  referencia = "No Obese",
  Titulo = "Score Plot PC1:\n(PCOS) Obese vs No Obese",
  componente = 1,colores = c("blue","lightblue")
)
p_contrib_PCOS_NOOB_vs_OB_PC1
p_contrib_PCOS_NOOB_vs_OB_PC2 <- fun_contribution(
  X = R[w_pcos, biomarcadores],
  scores = pcx$x[w_pcos, ],
  unfactor = pcos,
  referencia = "No Obese",
  Titulo = "Score Plot PC2:\n(PCOS) Obese vs No Obese",
  componente = 2,colores = c("blue","lightblue")
)
p_contrib_PCOS_NOOB_vs_OB_PC2

#####=====Female=======
p_contrib_Female_NOOB_vs_OB_PC1 <- fun_contribution(
  X = R[w_females, biomarcadores],
  scores = pcx$x[w_females, ],
  unfactor = females,
  referencia = "No Obese",
  Titulo = "Score Plot PC1:\n(Control Women) Obese vs No Obese",
  componente = 1,colores = c("blue","lightblue")
)
p_contrib_Female_NOOB_vs_OB_PC1
p_contrib_Female_NOOB_vs_OB_PC2 <- fun_contribution(
  X = R[w_females, biomarcadores],
  scores = pcx$x[w_females, ],
  unfactor = females,
  referencia = "No Obese",
  Titulo = "Score Plot PC2:\n(Control Women) Obese vs No Obese",
  componente = 2,colores = c("blue","lightblue")
)
p_contrib_Female_NOOB_vs_OB_PC2


plots_contrib <- list(
  p_contrib_obesidad_PC1 = p_contrib_obesidad_PC1,
  p_contrib_obesidad_PC2 = p_contrib_obesidad_PC2,
  p_contrib_Female_vs_Male_PC1 = p_contrib_Female_vs_Male_PC1,
  p_contrib_Female_vs_Male_PC2 = p_contrib_Female_vs_Male_PC2,
  p_contrib_PCOS_vs_Female_PC1 = p_contrib_PCOS_vs_Male_PC1,
  p_contrib_PCOS_vs_Female_PC2 = p_contrib_PCOS_vs_Female_PC2,
  p_contrib_PCOS_vs_Male_PC1 = p_contrib_PCOS_vs_Male_PC1,
  p_contrib_PCOS_vs_Male_PC2 = p_contrib_PCOS_vs_Male_PC2,
  p_contrib_NOOB_Males_vs_Females_PC1 = p_contrib_NOOB_Males_vs_Females_PC1,
  p_contrib_NOOB_Males_vs_Females_PC2 = p_contrib_NOOB_Males_vs_Females_PC2,
  p_contrib_NOOB_PCOS_vs_Females_PC1 = p_contrib_NOOB_PCOS_vs_Females_PC1,
  p_contrib_NOOB_PCOS_vs_Females_PC2 = p_contrib_NOOB_PCOS_vs_Females_PC2,
  p_contrib_NOOB_PCOS_vs_Males_PC1 = p_contrib_NOOB_PCOS_vs_Males_PC1,
  p_contrib_OB_Males_vs_Females_PC1 = p_contrib_OB_Males_vs_Females_PC1,
  p_contrib_OB_Males_vs_Females_PC2 = p_contrib_OB_Males_vs_Females_PC2,
  p_contrib_OB_PCOS_vs_Females_PC1 = p_contrib_OB_PCOS_vs_Females_PC1,
  p_contrib_OB_PCOS_vs_Females_PC2 = p_contrib_OB_PCOS_vs_Females_PC2,
  p_contrib_OB_PCOS_vs_Males_PC1 = p_contrib_OB_PCOS_vs_Males_PC1,
  p_contrib_OB_PCOS_vs_Males_PC2 = p_contrib_OB_PCOS_vs_Males_PC2,
  p_contrib_Female_NOOB_vs_OB_PC1 = p_contrib_Female_NOOB_vs_OB_PC1,
  p_contrib_Female_NOOB_vs_OB_PC2 = p_contrib_Female_NOOB_vs_OB_PC2,
  p_contrib_MALE_NOOB_vs_OB_PC1 = p_contrib_MALE_NOOB_vs_OB_PC1,
  p_contrib_MALE_NOOB_vs_OB_PC2 = p_contrib_MALE_NOOB_vs_OB_PC2,
  p_contrib_PCOS_NOOB_vs_OB_PC1 = p_contrib_PCOS_NOOB_vs_OB_PC1,
  p_contrib_PCOS_NOOB_vs_OB_PC2 = p_contrib_PCOS_NOOB_vs_OB_PC2
)

for(pl in 1:length(plots_contrib)){
  
  folder <- file.path(directorio_preliminar,"Contribuciones")
  if(!dir.exists(folder)){
    dir.create(folder,recursive = T)
  }
  ggsave(filename =file.path(folder,paste0(names(plots_contrib)[pl],".jpeg")),plot = plots_contrib[[pl]],height = 18,width = 18)
  
  
}
###=========outlier detection======

outlier_detection(pre[, biomarcadores], 2)
# outlier_detection(R, 2)


######======Univariante====#####
R_biomarcadores_ <- pre[, colnames(pre)[colnames(pre) %in% biomarcadores]]
res_R <- fun_plot_limma(R_biomarcadores_)
# res_R$logFC

p2 <- ggplot(res_R$logFC, aes(Variable, variable, fill = value, label = text)) + geom_tile() + scale_fill_gradient2(
  high = "red",
  low = "blue",
  mid = "white",
  limits = c(-2, 2)
) + theme(axis.text.x = element_text(
  angle = 90,
  vjust = 0.5,
  hjust = 1,
)) + geom_text() + theme(text = element_text(size = 9, face = "bold"))  + xlab("") +
  ylab("") + ggtitle("logFC")




p2
p2 + scale_y_discrete(labels =
                        ordern_contrastes)



######======Univariante====#####

R_obesidad <- R_biomarcadores_[obesidad == "Obese", ]
R_NoObesidad <- R_biomarcadores_[obesidad == "No Obese", ]
R_females <- R_biomarcadores_[grupo == "Female", ]
R_males <- R_biomarcadores_[grupo == "Male", ]
R_PCOS <- R_biomarcadores_[grupo == "PCOS", ]
R_malenoob <- R_biomarcadores_[grupo == "Male" &
                                 obesidad == "Obese", ]
R_maleoob <- R_biomarcadores_[grupo == "Male" &
                                obesidad == "No Obese", ]
R_PCOSob <- R_biomarcadores_[grupo == "PCOS" &
                               obesidad == "Obese", ]
R_PCOSnoob <- R_biomarcadores_[grupo == "PCOS" &
                                 obesidad == "No Obese", ]
R_femaleOb <- R_biomarcadores_[grupo == "Female" &
                                 obesidad == "Obese", ]
R_femaleNoob <- R_biomarcadores_[grupo == "Female" &
                                   obesidad == "No Obese", ]


grupos_R <- list(
  noob = R_NoObesidad,
  ob = R_obesidad,
  R_females = R_females,
  R_males = R_males,
  R_PCOS = R_PCOS,
  R_malenoob = R_malenoob,
  R_maleoob = R_maleoob,
  R_PCOSob = R_PCOSob,
  R_PCOSnoob = R_PCOSnoob,
  R_femaleOb = R_femaleOb,
  R_femaleNoob = R_femaleNoob,
  all = R[, biomarcadores]
)


plot_correlation <-  function(R, titulo) {
  colnames(R) <- colnames(pre)[colnames(pre) %in% colnames(R)]
  # R <- grupo1
  cor_test_result <- psych::corr.test(R)
  
  # Extraer la matriz de correlaciones
  cor_matrix <- as.data.frame(cor_test_result$r)  # Matriz de correlaciones
  diag(cor_matrix) <- 0
  
  # Extraer el vector de p-valores ajustados (fuera de la diagonal)
  p_adj_vector <- cor_test_result$p.adj  # Vector de p-valores ajustados
  
  # Crear una matriz vacía para los p-valores ajustados
  p_adj_matrix <- matrix(NA,
                         nrow = ncol(cor_matrix),
                         ncol = ncol(cor_matrix))
  
  # Rellenar la matriz vacía con los p-valores ajustados por encima de la diagonal
  p_adj_matrix[upper.tri(p_adj_matrix)] <- p_adj_vector
  
  # Hacer la matriz simétrica (rellenar la parte inferior de la diagonal con los mismos valores)
  p_adj_matrix[lower.tri(p_adj_matrix)] <- t(p_adj_matrix)[lower.tri(p_adj_matrix)]
  
  diag(p_adj_matrix) <- 1
  colnames(p_adj_matrix) <- rownames(p_adj_matrix) <- colnames(R)
  p_adj_matrix <- as.data.frame(p_adj_matrix)
  p_adj_matrix$Variable <- colnames(R)
  # Convertir las matrices de correlaciones y p-valores ajustados a formato largo (long format)
  cor_matrix$Variable <- colnames(R)
  cor_matrix$Variable <- factor(cor_matrix$Variable, colnames(R))
  cor_data <- reshape2::melt(cor_matrix, id.vars = "Variable")
  cor_data$Variable <- factor(cor_data$Variable, levels = colnames(R))
  p_data <- melt(p_adj_matrix, id.vars = "Variable")
  cor_data$variable <- factor(cor_data$variable, levels = colnames(R))
  
  # Crear una columna con los asteriscos para las celdas significativas (p.adj < 0.05)
  p_data$signif <- ifelse(p_data$value < 0.05, "*", "")
  
  # Fusionar las matrices de correlación y significancia ajustada
  plot_data <- merge(cor_data, p_data, by = c("Variable", "variable"))
  
  # Renombrar las columnas para mayor claridad
  colnames(plot_data) <- c("Var1", "Var2", "Correlation", "P_value_adj", "Significance")
  
  
  # Fusionar las matrices de correlación y significancia
  # Graficar con geom_tile y marcar las celdas significativas con asteriscos
  p <- ggplot(plot_data, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile() +
    geom_text(aes(label = Significance),
              color = "black",
              size = 2) +  # Añadir asteriscos en las celdas significativas
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      space = "Lab",
      name = "Pearson \nCorrelation Coefficient",
      limits = c(-1, 1)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 1,
        hjust = 1,
        size = 6
      ),
      axis.text.y = element_text(size = 6)
    ) +
    labs(x = "", y = "") + ggtitle(titulo)
  
  
  return(p)
  
  
}

nombres <- c(
  "Obese",
  "No Obese",
  "Female",
  "Male",
  "PCOS",
  "Male: No Obese",
  "Male: Obese",
  "PCOS: Obese",
  "PCOS: No Obese",
  "Female: Obese",
  "Female: No Obese",
  "total"
)

plots_ <-  vector("list", length = length(nombres))
for (gof in 1:length(nombres)) {
  plots_[[gof]] <- plot_correlation(grupos_R[[gof]], titulo = nombres[gof])
  
  
}



plots_


directorio_preliminar <- file.path(directorio, "Original")
if (!dir.exists(directorio_preliminar)) {
  dir.create(directorio_preliminar, recursive = T)
}

for (p in 1:length(plots_)) {
  ggsave(plot = plots_[[p]],
         filename = file.path(
           directorio_preliminar,
           paste0("Original: ", nombres[p], ".jpeg")
         ))
  
  
  
}


##=====OTROS GRAFICOS====
## Graficos y analisis de factores.
## Correlacion indirecta de los factores a traves de los pesos y las clnicas

factores2 <- as.data.frame(get_expectations(modelo,"Z",as.data.frame = F)[[1]])
factores2$Group <- grupo
interaccion <- interaction(grupo, obesidad)
levels(interaccion) <- c(
  "Control Women: No Obese",
  "PCOS: No Obese",
  "Men: No Obese",
  "Control Women: Obese",
  "PCOS: Obese",
  "Men: Obese"
)
factores2$subjects <-   interaccion
factores2$Obesity <- obesidad

ggscatter(factores2,"Factor1","Factor3",color="subjects",fill = "Obesity",facet.by = "Group")


csa <- rowSums(MOFA2::get_variance_explained(modelo)$r2_per_factor$group1)
plot1 <- ggplot(factores2, aes(
  x = Factor1,
  y = Factor2,
  shape = Group,
  color = subjects,
)) +
  geom_point(size = 7, aes(shape = Group, color = subjects)) +  # Bordes en negro
  scale_shape_manual(values = c(15, 16, 17)) +  # Definir las formas manualmente
  scale_color_manual(
    values = c(
      "Control Women: No Obese" = "green",
      "PCOS: No Obese" = "yellow",
      "Men: No Obese" = "red",
      "Control Women: Obese" = "seagreen4",
      "PCOS: Obese" = "#EEE685",
      "Men: Obese" = "#8B0000"
    )
  ) + # Definir los colores de relleno manualmente
  xlab(paste0("PC1: ", varianzas[1], "%")) +
  ylab(paste0("PC2: ", varianzas[2], "%"))+theme(text = element_text(size = 20))


#
MOFA2::plot_factor(modelo,
                   factors = 1,
                   group_by = grupo,
                   color_by = obesidad)
MOFA2::plot_factor(modelo,
                   factors = 2,
                   group_by = grupo,
                   color_by = obesidad)
MOFA2::plot_factor(modelo,
                   factors = 3,
                   group_by = grupo,
                   color_by = obesidad)
MOFA2::plot_factor(modelo,
                   factors = 4,
                   group_by = grupo,
                   color_by = obesidad)
MOFA2::plot_factor(modelo,
                   factors = 5,
                   group_by = grupo,
                   color_by = obesidad)
#
#

MOFA2::plot_factor(modelo,
                   factors = 1,
                   group_by = obesidad,
                   color_by = grupo)
MOFA2::plot_factor(modelo,
                   factors = 2,
                   group_by = obesidad,
                   color_by = grupo)
MOFA2::plot_factor(modelo,
                   factors = 3,
                   group_by = obesidad,
                   color_by = grupo)
MOFA2::plot_factor(modelo,
                   factors = 4,
                   group_by = obesidad,
                   color_by = grupo)
MOFA2::plot_factor(modelo,
                   factors = 5,
                   group_by = obesidad,
                   color_by = grupo)

