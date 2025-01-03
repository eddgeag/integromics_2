library(ggplot2)
library(psych)
library(car)
library(mixOmics)
library(vegan)
library(ggpubr)
library(MOFA2)
library(dplyr)
library(ALDEx2)
library(ggrepel)
library(patchwork)
library(forcats)
library(Bolstad)
library(tidyr)
library(limma)
library(cowplot)
library(ggsignif)
library(rlang)

datos <- readRDS("../datos/preprocesado_08_09_23/novoom-04-10-23.rds")

obesidad <- datos$comunes$obesidad
grupo <- datos$comunes$grupo
sexo <- datos$comunes$variables_in_bacteria$SEX
###====Establecemos los contrastes=====

## para obesidad

levels(obesidad) <- c("No.Obese", "Obese")

levels(grupo) <- c("Control.Women", "PCOS", "Men")
## contraste obesidad

contraste_obesidad <- c(Obese_vs_NoObese = "Obese-No.Obese")


## contraste para grupo
contrastes_grupo <- c(
  Men_vs_Female = c("Men-Control.Women"),
  PCOS_vs_Female = c("PCOS-Control.Women"),
  PCOS_vs_Male = c("PCOS-Men")
)

## interaccion

interaction_factor <- interaction(grupo, obesidad)
## cambiamos los nombres para que sean legibles
levels(interaction_factor) <- c(
  "Control.Women_No.Obese",
  "PCOS_No.Obese",
  "Male_No.Obese",
  "Control.Women_Obese",
  "PCOS_Obese",
  "Male_Obese"
)
## contrastes para itneraccion
contrastes_interaccion <- c(
  Female_Obese_vs_NoObese = "Control.Women_Obese-Control.Women_No.Obese",
  PCOS_Obese_vs_NoObese = "PCOS_Obese-PCOS_No.Obese",
  Male_Obese_Vs_NoObese = "Male_Obese-Male_No.Obese",
  No_Obese_Male_vs_Female = "Male_No.Obese-Control.Women_No.Obese",
  No_Obese_PCOS_vs_Female = "PCOS_No.Obese-Control.Women_No.Obese",
  No_Obese_PCOS_vs_Male = "PCOS_No.Obese-Male_No.Obese",
  Obese_Male_vs_Female = "Male_Obese-Control.Women_Obese",
  Obese_PCOS_vs_Female = "PCOS_Obese-Control.Women_Obese",
  Obese_PCOS_vs_Male = "PCOS_Obese-Male_Obese"
)

contrastes_lista <- list(
  contr_obesidad = contraste_obesidad,
  contr_grupo = contrastes_grupo,
  contr_interaccion = contrastes_interaccion
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

factores_interes <- list(obesidad = obesidad,
                         grupo = grupo,
                         interaccion = interaction_factor)

###====FUNCIONES======
fun_pca <- function(R,
                    scalar,
                    titulo,
                    biomarcadores,
                    text_size,
                    altura,
                    ancho,
                    dir_to_save,
                    name_to_save,
                    common_axes,
                    scale_labs,
                    legend_size,
                    size_points) {
  pcx <- prcomp(R[, biomarcadores], scale. = scalar)
  
  plotdf <- data.frame(pcx$x,
                       Group = grupo,
                       Obesity = obesidad,
                       sujetos = rownames(R))
  varianzas <- round(100 * (pcx$sdev^2) / sum(pcx$sdev^2), 2)
  
  levels(plotdf$Obesity) <- c("No Obese", "Obese")
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
  
  p4 <- factoextra::fviz_screeplot(pcx) + ylab("% Variance") + ggtitle("") + theme(
    axis.title.y = element_text(size = text_size),
    axis.title.x = element_text(size = text_size),
    axis.text.x = element_text(size = text_size * scale_labs),
    axis.text.y = element_text(size = text_size * scale_labs)
  )
  
  levels(plotdf$Group) <- c("Control Women", "PCOS", "Men")
  
  colores <- ifelse(plotdf$subjects == "Control Women: No Obese", "chartreuse", NA)
  colores <- ifelse(plotdf$subjects == "PCOS: No Obese", "gold33", colores)
  colores <- ifelse(plotdf$subjects == "Men: No Obese", "red", colores)
  colores <- ifelse(plotdf$subjects == "Control Women: Obese" ,
                    "chartreuse4",
                    colores)
  colores <- ifelse(plotdf$subjects == "PCOS: Obese", "darkorange2", colores)
  colores <- ifelse(plotdf$subjects == "Men: Obese", "darkred", colores)
  
  formas <- ifelse(plotdf$Group == "Control Women", 0, NA)
  formas <- ifelse(plotdf$Group == "PCOS", 1, formas)
  formas <- ifelse(plotdf$Group == "Men", 2, formas)
  
  plotdf$colores <- colores
  plotdf$formas <- formas
  
  plot1 <- ggplot(plotdf, aes(
    x = PC1,
    y = PC2,
    shape = Group,
    color = subjects
  )) +
    geom_point(size = size_points) +
    scale_color_manual(
      values = c(
        "Control Women: No Obese" = "chartreuse",
        "PCOS: No Obese" = "gold3",
        "Men: No Obese" = "red",
        "Control Women: Obese" = "chartreuse4",
        "PCOS: Obese" = "darkorange2",
        "Men: Obese" = "darkred"
      )
    ) +
    scale_shape_manual(values = c(
      "Control Women" = 15,
      "PCOS" = 16,
      "Men" = 17
    )) +
    xlab(paste0("PC1: ", varianzas[1], "%")) +
    ylab(paste0("PC2: ", varianzas[2], "%")) + ggtitle("All") +
    theme(
      axis.title.y = element_text(size = text_size),
      axis.title.x = element_text(size = text_size),
      axis.text.x = element_text(size = text_size *  scale_labs),
      axis.text.y = element_text(size = text_size * scale_labs),
      plot.title = element_text(size = text_size * 0.5),
      legend.title = element_blank()
    ) +
    guides(shape = guide_legend(order = 1), color = guide_legend(order = 2))
  
  plot2 <- ggplot(subset(plotdf, Obesity == "Obese"),
                  aes(
                    x = PC1,
                    y = PC2,
                    shape = Group,
                    color = subjects
                  )) +
    geom_point(size = size_points) +
    scale_shape_manual(values = c(15, 16, 17)) +
    scale_color_manual(
      values = c(
        "Control Women: Obese" = "chartreuse4",
        "PCOS: Obese" = "darkorange2",
        "Men: Obese" = "darkred"
      )
    ) +
    xlab("") +
    ylab("") + ggtitle("Obese") +
    theme(
      legend.position = "none",
      axis.title.y = element_text(size = text_size),
      axis.title.x = element_text(size = text_size),
      axis.text.x = element_text(size = text_size * scale_labs),
      axis.text.y = element_text(size = text_size * scale_labs),
      plot.title = element_text(size = text_size * 0.5),
      legend.title = element_blank()
    )
  
  plot3 <- ggplot(
    subset(plotdf, Obesity == "No Obese"),
    aes(
      x = PC1,
      y = PC2,
      shape = Group,
      color = subjects
    )
  ) +
    geom_point(size = size_points) +
    scale_shape_manual(values = c(15, 16, 17)) +
    scale_color_manual(
      values = c(
        "Control Women: No Obese" = "chartreuse",
        "PCOS: No Obese" = "gold3",
        "Men: No Obese" = "red"
      )
    ) +
    xlab("") +
    ylab(paste0("PC2: ", varianzas[2], "%")) + ggtitle("No Obese") +
    theme(
      legend.position = "none",
      axis.title.y = element_text(size = text_size),
      axis.title.x = element_text(size = text_size),
      axis.text.x = element_text(size = text_size * scale_labs),
      axis.text.y = element_text(size = text_size * scale_labs),
      plot.title = element_text(size = text_size * 0.5),
      legend.title = element_blank()
    )
  
  common_xlim <- c(-common_axes, common_axes)
  common_ylim <- c(-common_axes, common_axes)
  
  plot1 <- plot1 + scale_x_continuous(limits = common_xlim) + scale_y_continuous(limits = common_ylim)
  plot2 <- plot2 + scale_x_continuous(limits = common_xlim) + scale_y_continuous(limits = common_ylim)
  plot3 <- plot3 + scale_x_continuous(limits = common_xlim) + scale_y_continuous(limits = common_ylim)
  
  combined_plot <- (plot3 + plot2) / (plot1 + p4) +
    plot_layout(guides = "collect") +
    plot_annotation(title = titulo) &
    theme(
      axis.ticks = element_line(),
      axis.text = element_text(),
      legend.key.size = unit(legend_size, "cm"),
      legend.text = element_text(size = text_size * scale_labs),
      legend.title = element_blank(),
      plot.title = element_text(size =  scale_labs * text_size)
    )
  
  ggsave(
    plot = combined_plot,
    filename = file.path(dir_to_save, paste0(name_to_save, ".jpeg")),
    height = altura,
    width = ancho
  )
  
  return(combined_plot)
}



# Función para combinar makeContrasts y t.test
contrast_t_test <- function(datas, factor_, response, contrast_expr) {
  # Verificar que el factor y la respuesta están en el dataframe
  if (!all(c(factor_, response) %in% colnames(datas))) {
    stop("El factor y la variable de respuesta deben estar en las columnas del dataframe.")
  }
  
  # Extraer los niveles del factor y crear la matriz de diseño
  levels_factor <- levels(datas[[factor_]])
  design <- model.matrix( ~ 0 + datas[[factor_]])  # Diseño sin intercepto
  
  # Asignar nombres a las columnas de la matriz de diseño
  colnames(design) <- levels_factor
  
  # Crear la matriz de contrastes con makeContrasts
  contrasts <- makeContrasts(contrasts = contrast_expr, levels = levels(datas[[factor_]]))
  
  # Verificar que la matriz de contrastes sea válida
  if (nrow(contrasts) != ncol(design)) {
    stop("El contraste especificado no coincide con los niveles del factor.")
  }
  
  # Aplicar el contraste al diseño para obtener grupos
  contrast_values <- design %*% contrasts
  
  # Extraer los grupos contrastados para el t-test
  x <- datas[[response]][contrast_values[, 1] > 0]
  y <- datas[[response]][contrast_values[, 1] < 0]
  
  
  
  # Compute empirical estimates
  m_x <- mean(x)
  m_y <- mean(y)
  v_x <- var(x)
  v_y <- var(y)
  
  # Set prior parameters based on empirical estimates
  m <- c(m_x, m_y)            # Prior means
  n0 <- c(1 / v_x, 1 / v_y)       # Prior precisions (inverse of variances)
  sig.med <- median(c(sd(x), sd(y)))  # Prior median standard deviation
  kappa <- 1                  # Degrees of freedom for prior on sigma
  
  # Perform Bayesian t-test with empirical priors
  t_test_result <- bayes.t.test(
    x,
    y,
    var.equal = F,
    prior = "joint.conj",
    m = m,
    n0 = n0,
    sig.med = sig.med,
    kappa = kappa
  )
  
  
  return(t_test_result)
}

contrasts_bayes <- function(X, initial_factor, contrast_expr) {
  ### Necesitamos que el X este ya con los nombres correctos.
  levels_in_contrasts <-  strsplit2(contrast_expr, "-")
  ## downsampleamos o no el factor.
  factor_of_interest <- factor(initial_factor, levels = levels_in_contrasts)
  ## downsampleamos o no el data freame X
  df_of_interest <- X[initial_factor %in%  factor_of_interest, ]
  ## anadimos el factor al df
  df_of_interest <- as.data.frame(df_of_interest) ## nos aseguramos que sea df
  factor_of_interest <- factor_of_interest[which(!is.na(factor_of_interest))]
  df_of_interest$factor_of_interest <- factor_of_interest
  ## convertimos el nombre del factor de interes
  name_of_contrast <- names(contrast_expr)
  colnames(df_of_interest)[ncol(df_of_interest)] <- name_of_contrast
  
  
  ## creamos una matriz de salida
  
  out_matrix <-  matrix(NA, nrow = ncol(X), ncol = 3) #3 porque son p valor stat y diff
  ## iteramos para cada variable, excepto el factor de interes
  
  for (j in 1:ncol(X)) {
    ## X porque todavia no hemos puesto el factor
    variable_j <- colnames(X)[j]
    contrast_j <- contrast_t_test(
      df_of_interest,
      name_of_contrast,
      response = variable_j,
      contrast_expr = contrast_expr
    )
    
    ## extraer estadistico
    statistic <- contrast_j$statistic
    ## extrar diferencia de posterior
    posterior_Dif <- -diff(contrast_j$estimate)
    ## extraer p valor
    p.value <- contrast_j$p.value
    res_j <- c(t_stat = statistic,
               posterior_diff = posterior_Dif,
               p_value = p.value)
    out_matrix[j, ] <- res_j
    
  }
  
  colnames(out_matrix) <- c("t_stat", "post_dif", "p_val")
  rownames(out_matrix) <- colnames(X)
  out_df <- as.data.frame(out_matrix)
  return(out_df)
  
}



analyze_2by2_t_bayes <- function(df_to_analyze,
                                 factor1,
                                 factor2,
                                 factor_interact,
                                 list_contrasts) {
  contrast_1 <- list_contrasts[[1]]
  contrast_2 <- list_contrasts[[2]]
  contrast_3 <- list_contrasts[[3]]
  list_to_analyze <- list(
    factor1 = list(categoric = factor1, contrast_categoric = contrast_1),
    factor2 = list(categoric = factor2, contrast_categoric = contrast_2),
    factor_interact = list(categoric = factor_interact, contrast_categoric = contrast_3)
  )
  
  ## nota mental
  ## seq_along = 1:length(vector)
  list_output_t_bayes <- vector(mode = "list", length = length(list_to_analyze))
  for (idx_of_interest in seq_along(list_to_analyze)) {
    list_element <- list_to_analyze[[idx_of_interest]]
    factor_i <- list_element$categoric
    expresion_contrastes_i <- list_element$contrast_categoric
    list_idx <- vector(mode = "list",
                       length = length(expresion_contrastes_i))
    for (contraste_idx in seq_along(expresion_contrastes_i)) {
      contraste_j <- expresion_contrastes_i[contraste_idx]
      list_idx[[contraste_idx]] <- contrasts_bayes(
        X = df_to_analyze ,
        initial_factor = factor_i,
        contrast_expr = contraste_j
      )
      
    }
    names(list_idx) <- names(expresion_contrastes_i)
    
    list_output_t_bayes[[idx_of_interest]] <- list_idx
    
    
  }
  
  names(list_output_t_bayes) <- names(list_to_analyze)
  
  
  ## extraer t dif p adj
  num_of_contrasts <- sum(unlist(sapply(list_output_t_bayes, length)))
  t_stat.list <- vector(mode = "list", length = length(list_output_t_bayes))
  post_diffs.list <-  vector(mode = "list", length = length(list_output_t_bayes))
  p_vals.list <-  vector(mode = "list", length = length(list_output_t_bayes))
  for (l in seq_along(list_output_t_bayes)) {
    l_factor <- list_output_t_bayes[[l]]
    tstat.aux <- matrix(NA,
                        nrow = length(l_factor),
                        ncol = ncol(df_to_analyze))
    posts.aux <- matrix(NA,
                        nrow = length(l_factor),
                        ncol = ncol(df_to_analyze))
    pvals.aux <- matrix(NA,
                        nrow = length(l_factor),
                        ncol = ncol(df_to_analyze))
    
    for (con_ in seq_along(l_factor)) {
      con_j <- l_factor[[con_]]
      tstat.aux[con_, ] <- con_j$t_stat
      posts.aux[con_, ] <- con_j$post_dif
      pvals.aux[con_, ] <- con_j$p_val
      
    }
    t_stat.list[[l]] <- as.data.frame(tstat.aux)
    post_diffs.list[[l]] <- as.data.frame(posts.aux)
    p_vals.list[[l]] <- as.data.frame(pvals.aux)
    
  }
  
  
  contrasts_names <- unlist(lapply(list_output_t_bayes, names))
  
  t_stats_bayes <- bind_rows(t_stat.list)
  post_difs <- bind_rows(post_diffs.list)
  pvals_bayes <- bind_rows(p_vals.list)
  
  
  colnames(t_stats_bayes) <- colnames(post_difs) <- colnames(pvals_bayes) <- colnames(df_to_analyze)
  rownames(t_stats_bayes) <- rownames(post_difs) <- rownames(pvals_bayes) <- contrasts_names
  
  t_stats_bayes$contrasts_names <-   factor(rownames(t_stats_bayes), levels =
                                              rownames(t_stats_bayes))
  post_difs$contrasts_names <- factor(rownames(t_stats_bayes), levels =
                                        rownames(t_stats_bayes))
  pvals_bayes$contrasts_names <- factor(rownames(t_stats_bayes), levels =
                                          rownames(t_stats_bayes))
  
  to_Adj <- pvals_bayes[5:nrow(pvals_bayes), -ncol(pvals_bayes)]
  to_Adj <- apply(to_Adj, 2, function(x)
    p.adjust(x, "BH"))
  pvals_Adjs <- pvals_bayes
  pvals_Adjs[5:nrow(pvals_bayes), -ncol(pvals_bayes)] <- to_Adj
  list_to_return <- list(
    to_plot = list(
      tstats = t_stats_bayes,
      post_difs = post_difs,
      pvals_bayes = pvals_bayes,
      pvals_adj = pvals_Adjs
    ),
    to_extract_features = list_output_t_bayes
  )
  
  return(list_to_return)
  
  
}



fun_plot_contrasts <- function(tstats,
                               pvals_bayes,
                               post_difs,
                               threshold,
                               size_text) {
  tstat.melt <- reshape2::melt(tstats)
  pvals.melt <- reshape2::melt(pvals_bayes)
  posts.melt <- reshape2::melt(post_difs)
  
  df.melt_toplot <- data.frame(
    feature = tstat.melt$variable,
    tstat = tstat.melt$value,
    pval = pvals.melt$value,
    stars = ifelse(pvals.melt$value < threshold, "*", ""),
    posts = posts.melt$value,
    contrasts = tstat.melt$contrasts_names
  )
  
  p_ <- ggplot(df.melt_toplot,
               aes(feature, contrasts, fill = tstat, label = stars)) + geom_tile() + scale_fill_gradient2(high = "red",
                                                                                                          low = "blue",
                                                                                                          mid = "white") + theme(axis.text.x = element_text(
                                                                                                            angle = 90,
                                                                                                            vjust = 0.5,
                                                                                                            hjust = 1
                                                                                                          )) + geom_text(size = 5) + theme(text = element_text(
                                                                                                            size = size_text,
                                                                                                            face = "bold",
                                                                                                            color = "black"
                                                                                                          ))  + xlab("") +
    ylab("") + ggtitle(paste0("Empirical Bayesian T test: Pvals ", threshold)) + scale_y_discrete(labels = ordern_contrastes)
  return(p_)
}
feature_selection <- function(omic_features,
                              sujetos,
                              f,
                              threshold,
                              factores,
                              pesos,
                              thresh_pval) {
  ## 1 tomar aquellas proteinas con peso mayor al threshold
  ## factores de la omica
  ## reconstruyo la senal de la omica
  
  omic_r <- factores[, f] %*% t(pesos[omic_features, f])
  ## obtengo los pesos de la omica de interes
  pesos_omic <- pesos[omic_features, f]
  
  ## si solo es un factor
  if (length(f) > 1) {
    ## si es para mas de un factor
    maximum <- apply(pesos_omic, 1, function(x)
      max(abs(x)))
    biomarcadores_omic <- names(maximum[abs(maximum) > threshold])
  } else{
    biomarcadores_omic <- names(pesos_omic[abs(pesos_omic) > threshold])
    
  }
  reconstruct_omic <-  as.data.frame(factores[, f] %*% t(pesos[biomarcadores_omic, f]), row.names = sujetos)
  ## si es para mas de un facto
  
  res <- analyze_2by2_t_bayes(omic_r,
                              obesidad,
                              grupo,
                              interaction_factor,
                              contrastes_lista)
  
  
  pvals <- res$to_plot$pvals_adj
  logic <- (apply(pvals[, -ncol(pvals)], 2, function(x)
    any(x < thresh_pval)))
  significant <- data.frame(biomarcadores = names(which(logic)))
  biomarcadores_omic <- data.frame(biomarcadores = biomarcadores_omic)
  
  biomarcadores_omic_final <- full_join(significant, biomarcadores_omic)
  
  return(biomarcadores_omic_final$biomarcadores)
  
}

feature_selection_threshold <- function(omic_features,
                                        sujetos,
                                        f,
                                        threshold,
                                        factores,
                                        pesos) {
  ## 1 tomar aquellas proteinas con peso mayor al threshold
  ## factores de la omica
  ## reconstruyo la senal de la omica
  
  omic_r <- factores[, f] %*% t(pesos[omic_features, f])
  ## obtengo los pesos de la omica de interes
  pesos_omic <- pesos[omic_features, f]
  
  ## si solo es un factor
  if (length(f) > 1) {
    ## si es para mas de un factor
    maximum <- apply(pesos_omic, 1, function(x)
      max(abs(x)))
    biomarcadores_omic <- names(maximum[abs(maximum) > threshold])
  } else{
    biomarcadores_omic <- names(pesos_omic[abs(pesos_omic) > threshold])
    
  }
  ## si es para mas de un facto
  
  
  
  return(biomarcadores_omic)
  
}

dummy_plot_weight_factor <- function(modelo, biomarcadores, dir_to_save) {
  pesos <- as.data.frame(Reduce(rbind, get_weights(modelo, scale = T)))
  pesos_biomarcadores <- pesos[biomarcadores, ]
  pesos_biomarcadores$biomarcador <- rownames(pesos_biomarcadores)
  pesos_biomarcadores$biomarcador <- as.factor(pesos_biomarcadores$biomarcador)
  pesos_biomarcadores$omic <- ifelse(
    rownames(pesos_biomarcadores) %in% rownames(modelo@data$metaboloma$group1),
    "Metabolome",
    NA
  )
  pesos_biomarcadores$omic <- ifelse(
    rownames(pesos_biomarcadores) %in% rownames(modelo@data$metagenoma$group1),
    "Microbiome",
    pesos_biomarcadores$omic
  )
  pesos_biomarcadores$omic <- ifelse(
    rownames(pesos_biomarcadores) %in% rownames(modelo@data$ip$group1),
    "Proteins",
    pesos_biomarcadores$omic
  )
  pesos_biomarcadores$omic <- as.factor(pesos_biomarcadores$omic)
  # Plot
  factor1 <- pesos_biomarcadores %>% mutate(biomarcador = fct_reorder(biomarcador, Factor1)) %>% ggplot(aes(x =
                                                                                                              biomarcador, y = Factor1, color = omic)) +
    geom_segment(aes(
      x = biomarcador,
      xend = biomarcador,
      y = 0,
      yend = Factor1,
      colour = omic
    )) +
    geom_point(size = 5) +
    coord_flip() +
    theme(legend.position = "none") +
    xlab("") +
    ylab("Weight") +
    ggtitle("Factor 1") + theme(text = element_text(size = 20)) + ylim(c(-1, 1)) +
    scale_colour_manual(values = c("blue", "yellow", "gray"))
  # Combinar la gráfica original con la columna de símbolos
  
  # Mostrar la gráfica final
  
  
  # Plot
  factor2 <- pesos_biomarcadores %>% mutate(biomarcador = fct_reorder(biomarcador, Factor2)) %>% ggplot(aes(x =
                                                                                                              biomarcador, y = Factor2, color = omic)) +
    geom_segment(aes(
      x = biomarcador,
      xend = biomarcador,
      y = 0,
      yend = Factor2,
      colour = omic
    )) +
    geom_point(size = 5) +
    coord_flip() +
    theme(legend.position = "none") +
    xlab("") +
    ylab("Weight") +
    ggtitle("Factor 2") + theme(text = element_text(size = 20)) + ylim(c(-1, 1)) +
    scale_colour_manual(values = c("blue", "yellow", "gray"))
  # Combinar la gráfica original con la columna de símbolos
  
  
  
  
  # Plot
  factor3 <- pesos_biomarcadores %>% mutate(biomarcador = fct_reorder(biomarcador, Factor3)) %>% ggplot(aes(x =
                                                                                                              biomarcador, y = Factor3 , color = omic)) +
    geom_segment(aes(
      x = biomarcador,
      xend = biomarcador,
      y = 0,
      yend = Factor3,
      colour = omic
    )) +
    geom_point(size = 5) +
    coord_flip() +
    theme(legend.position = "none") +
    xlab("") +
    ylab("Weight") +
    ggtitle("Factor 3") + theme(text = element_text(size = 20)) + ylim(c(-1, 1)) +
    scale_colour_manual(values = c("blue", "yellow", "gray"))
  # Combinar la gr
  
  
  
  # Plot
  factor4 <- pesos_biomarcadores %>% mutate(biomarcador = fct_reorder(biomarcador, Factor4)) %>% ggplot(aes(x =
                                                                                                              biomarcador, y = Factor4, color = omic)) +
    geom_segment(aes(
      x = biomarcador,
      xend = biomarcador,
      y = 0,
      yend = Factor4,
      colour = omic
    )) +
    geom_point(size = 5) +
    coord_flip() +
    theme(legend.position = "none") +
    xlab("") +
    ylab("Weight") +
    ggtitle("Factor 4") + theme(text = element_text(size = 20)) + ylim(c(-1, 1)) +
    scale_colour_manual(values = c("blue", "yellow", "gray"))
  # Combinar la gr
  
  
  # Plot
  factor5 <- pesos_biomarcadores %>% mutate(biomarcador = fct_reorder(biomarcador, Factor5)) %>% ggplot(aes(x =
                                                                                                              biomarcador, y = Factor5, color = omic)) +
    geom_segment(aes(
      x = biomarcador,
      xend = biomarcador,
      y = 0,
      yend = Factor5,
      colour = omic
    )) +
    geom_point(size = 5) +
    coord_flip() +
    theme(legend.position = "none") +
    xlab("") +
    ylab("Weight") +
    ggtitle("Factor 5") + theme(text = element_text(size = 20)) + ylim(c(-1, 1)) +
    scale_colour_manual(values = c("blue", "yellow", "gray"))
  # Combinar la gr
  
  ggsave(
    filename = file.path(dir_to_save, "Weights_Factor1.jpeg"),
    plot = factor1,
    height = 18,
    width = 15
  )
  
  
  ggsave(
    filename = file.path(dir_to_save, "Weights_Factor2.jpeg"),
    plot = factor2,
    height = 18,
    width = 15
  )
  
  ggsave(
    filename = file.path(dir_to_save, "Weights_Factor3.jpeg"),
    plot = factor3,
    height = 18,
    width = 15
  )
  
  ggsave(
    filename = file.path(dir_to_save, "Weights_Factor4.jpeg"),
    plot = factor4,
    height = 18,
    width = 15
  )
  
  ggsave(
    filename = file.path(dir_to_save, "Weights_Factor5.jpeg"),
    plot = factor5,
    height = 18,
    width = 15
  )
  
  list_to_return <- list(
    factor1 = factor1,
    factor2 = factor2,
    factor3 = factor3,
    factor4 = factor4,
    factor5 = factor5
  )
  return(list_to_return)
}

wrapper_feature_selection_method_MOFA <- function(modelo, threshold, fs) {
  factores <- get_expectations(modelo, "Z", as.data.frame = F)[[1]]
  pesos <- Reduce("rbind", get_expectations(modelo, "W"))
  omic_metabolome <- rownames(modelo@data$metaboloma$group1)
  omic_ip <- rownames(modelo@data$ip$group1)
  omic_metagenoma <- rownames(modelo@data$metagenoma$group1)
  sujetos <- colnames(modelo@data$metaboloma$group1)
  
  biomarkers_metabolome <- feature_selection_threshold(
    omic_features = omic_metabolome,
    sujetos = sujetos,
    f = fs$metaboloma,
    threshold = threshold,
    factores = factores,
    pesos = pesos
  )
  
  
  biomarkers_ip <- feature_selection_threshold(
    omic_features = omic_ip,
    sujetos = sujetos,
    f = fs$ip,
    threshold = threshold,
    factores = factores,
    pesos = pesos
  )
  
  biomarkers_metagenoma <- feature_selection_threshold(
    omic_features = omic_metagenoma,
    sujetos = sujetos,
    f = fs$metagenoma,
    threshold = threshold,
    factores = factores,
    pesos = pesos
  )
  
  
  biomarcadores <- c(biomarkers_metabolome,
                     biomarkers_ip,
                     biomarkers_metagenoma)
  R <- factores %*% t(pesos)
  
  biomarcadores_mofa <- biomarcadores[match(colnames(R), biomarcadores, nomatch = 0)]
  
  return(biomarcadores_mofa)
  
}
wrapper_feature_selection_method_mon <- function(modelo, threshold, thresh_pval, fs) {
  factores <- get_expectations(modelo, "Z", as.data.frame = F)[[1]]
  pesos <- Reduce("rbind", get_expectations(modelo, "W"))
  omic_metabolome <- rownames(modelo@data$metaboloma$group1)
  omic_ip <- rownames(modelo@data$ip$group1)
  omic_metagenoma <- rownames(modelo@data$metagenoma$group1)
  sujetos <- colnames(modelo@data$metaboloma$group1)
  
  biomarkers_metabolome <- feature_selection(
    omic_features = omic_metabolome,
    sujetos = sujetos,
    f = fs$metaboloma,
    threshold = threshold,
    factores = factores,
    pesos = pesos,
    thresh_pval = thresh_pval
  )
  
  
  biomarkers_ip <- feature_selection(
    omic_features = omic_ip,
    sujetos = sujetos,
    f = fs$ip,
    threshold = threshold,
    factores = factores,
    pesos = pesos,
    thresh_pval = thresh_pval
  )
  
  biomarkers_metagenoma <- feature_selection(
    omic_features = omic_metagenoma,
    sujetos = sujetos,
    f = fs$metagenoma,
    threshold = threshold,
    factores = factores,
    pesos = pesos,
    thresh_pval = thresh_pval
  )
  
  
  biomarcadores <- c(biomarkers_metabolome,
                     biomarkers_ip,
                     biomarkers_metagenoma)
  R <- factores %*% t(pesos)
  
  biomarcadores <- biomarcadores[match(colnames(R), biomarcadores, nomatch = 0)]
  
  return(biomarcadores)
  
}

reconstruccion <- function(modelo) {
  factores <- get_expectations(modelo, "Z", as.data.frame = F)[[1]]
  pesos <- Reduce("rbind", get_expectations(modelo, "W"))
  R <- factores %*% t(pesos)
  return(R)
  
}
plot_contrastes_wrapper <- function(R,
                                    dir_to_save,
                                    name_to_save,
                                    altura,
                                    ancho,
                                    size_text) {
  res_contrastes <- analyze_2by2_t_bayes(
    R,
    factor1 = obesidad,
    factor2 = grupo,
    factor_interact = interaction_factor,
    list_contrasts = contrastes_lista
  )
  
  p1 <- fun_plot_contrasts(
    tstats = res_contrastes$to_plot$tstats,
    pvals_bayes = res_contrastes$to_plot$pvals_adj,
    post_difs = res_contrastes$to_plot$post_difs,
    threshold = 0.05,
    size_text = size_text
  )
  
  ggsave(
    plot = p1,
    filename = file.path(dir_to_save, paste0(name_to_save, ".jpeg")),
    width = ancho,
    height = altura
  )
  return(p1)
  
  
}
wrapper_pca <- function(X,
                        biomarcadores,
                        scale,
                        titulo,
                        dir_to_save,
                        name_to_save) {
  p1 <- fun_pca(X, biomarcadores, titulo = titulo)
  
  ggsave(plot = p1, filename = file.path(dir_to_save, paste0(name_to_save, ".jpeg")))
  return(p1)
}

fun_dummy_plot_clinical_latent <- function(modelo,
                                           dir_to_save,
                                           name_to_save,
                                           ancho ,
                                           altura,
                                           factor_,
                                           text_size) {
  feat_fuera <- c("group", "sample")
  idx <- sapply(feat_fuera, function(x)
    grep(x, colnames(modelo@samples_metadata)))
  
  ## reorder clinical variables
  factores.df <- as.data.frame(MOFA2::get_expectations(modelo, "Z", as.data.frame = F)[[1]])
  
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
  
  ## variables ponemos para graficar bonitas
  clinical_bonito <- c(
    "BMI",
    "Waist circumference",
    "Waist to Hip Ratio",
    "Total Testosterone",
    "Free Testosterone",
    "Total Estradiol",
    "Free Estradiol",
    "SHBG",
    "Glucose",
    "Insulin",
    "HOMA IR",
    "ISI",
    "Triglycerides",
    "Cholesterol",
    "HDL Cholesterol",
    "LDL Cholesterol"
  )
  
  correlation_matrix <- correlation_matrix[, clinical]
  
  
  p_value_matrix <- t(corr_result$p.adj)[, clinical]
  # colnames(p_value_matrix) <- clinical
  #Get the p-values for matching columns
  
  # Create a data frame for ggplot2
  
  
  correlation_matrix <- correlation_matrix[rev(rownames(correlation_matrix)), ]
  p_value_matrix <- p_value_matrix[rev(rownames(p_value_matrix)), ]
  
  data <- reshape2::melt(correlation_matrix)
  data$pvals <- reshape2::melt(p_value_matrix)$value
  
  # Add a column to indicate significance based on p-value < 0.05
  data$significant <- ifelse(data$pvals < 0.05, "*", "")
  
  
  
  # Plot using ggplot2 with circles
  p1 <-  ggplot(data, aes(x = Var2, y = Var1)) +
    geom_tile(color = "grey",
              fill = NA,
              linewidth = 0.5) +  # Cuadrícula
    geom_point(aes(
      fill = value,
      color = value,
      size = abs(value)
    ), shape = 16) + scale_size(range = c(1, factor_), guide = "none") + # Círculos con colores iguales
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      guide = "none"
    ) +  # Escala de relleno sin leyenda
    scale_color_gradient2(low = "blue",
                          mid = "white",
                          high = "red") +  # Escala de borde con leyenda
    geom_text(aes(label = significant),
              color = "black",
              size = 7) +  # Texto en cada celda
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        colour = "black",
        size = text_size
      ),
      axis.title = element_blank(),
      axis.text.y = element_text(colour = "black", size = text_size),
      title = element_text(size = text_size * 2, color = "black"),
      legend.key.size = unit(2, "cm"),
      legend.title = element_text(size = text_size * 1.3),
      legend.text = element_text(size = text_size)
      
      
    ) +
    labs(color = "Pearson\nCorrelation Coefficient")  + ggtitle("Latent Variables correlated with Covariables") + scale_y_discrete(labels = c("Factor 5", "Factor 4", "Factor 3", "Factor 2", "Factor 1")) #tiqueta única para la leyenda
  
  ggsave(
    filename = file.path(dir_to_save, paste0(name_to_save, ".jpeg")),
    plot = p1,
    height = altura,
    width = ancho
  )
  return(p1)
}


fun_plot_dummy_preliminar <- function(modelo,
                                      dir_to_save,
                                      altura,
                                      ancho,
                                      size_var,
                                      size_Axis) {
  data.variance_explained <- get_variance_explained(modelo)
  data.variance_explained <- as.data.frame(data.variance_explained$r2_per_factor$group1)
  data.variance_explained$Factors <- rownames(data.variance_explained)
  # Preparar datos para el heatmap (formato largo)
  data.melt <- reshape2::melt(data.variance_explained)
  
  # Crear el heatmap
  (
    heatmap_plot <- ggplot(data.melt, aes(
      x = variable, y = Factors, fill = value
    )) +
      geom_tile(color = "white") +
      scale_fill_gradient(
        low = "white",
        high = "blue",
        name = "Var. (%)"
      ) +
      geom_text(aes(label = round(value, 2)), size = size_var) +
      labs(title = "Variance Explained by View and Factor", x = NULL, y = NULL) +
      scale_x_discrete(labels = c("Metabolome", "Microbiome", "Proteins")) + theme(
        axis.text.x = element_text(size = size_Axis, colour = "black"),
        axis.text.y = element_text(size = size_Axis, colour = "black"),
        legend.text = element_text(size = size_Axis),
        plot.title = element_text(size = size_Axis +
                                    5, face = "bold"),
        legend.key.size = unit(2, "cm"),
        legend.title = element_text(size =
                                      size_Axis)
      )
  )
  
  # Preparar datos para el barplot
  bar_data <- data.melt %>%
    group_by(variable) %>%
    summarise(Total = sum(value))
  levels(bar_data$variable) <- c("Metabolome", "Microbiome", "Proteins")
  # Crear el barplot
  p2 <- ggbarplot(bar_data,
                  x = "variable",
                  y = "Total",
                  fill = "variable") +
    scale_fill_manual(values = c(
      Metabolome = "blue",
      Microbiome = "yellow",
      Proteins = "grey"
    )) +
    labs(title = "Variance Explained by View", x = "", y = "Var. (%)") +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = size_Axis, colour = "black"),
      plot.title = element_text(size = size_Axis + 5, face = "bold"),
      axis.title.y  = element_text(size = size_Axis),
      axis.text.y = element_text(size = size_Axis)
    )
  
  # Combinar ambos gráficos en una fila
  combined_plot <- plot_grid(heatmap_plot,
                             p2,
                             ncol = 2,
                             rel_widths = c(1, 1))  # Ajustar proporciones)
  
  ggsave(
    plot = combined_plot,
    filename = file.path(dir_to_save, "graficos_preliminares.jpeg"),
    width = ancho,
    height = altura
  )
  
  return(combined_plot)
}



plot_correlation <-  function(R, titulo, text_size) {
  # R <- R[,biomarcadores]
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
  p_data <- reshape2::melt(p_adj_matrix, id.vars = "Variable")
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
              size = 5) +  # Añadir asteriscos en las celdas significativas
    scale_fill_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      space = "Lab",
      name = "Pearson\nCorrelation\nCoefficient",
      limits = c(-1, 1),
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = text_size,
        color = "black"
      ),
      axis.text.y = element_text(size = text_size, color = "black"),
      legend.text = element_text(size = text_size, color = "black"),
      legend.title = element_text(size = text_size * 1.5, color = "black"),
      plot.title = element_text(size = 2 * text_size, color = "black"),
      legend.key.size = unit(1, "cm")
    ) +
    labs(x = "", y = "") + ggtitle(titulo)
  
  
  return(p)
  
  
}

plot_correlation_wrapper <- function(R,
                                     biomarcadores,
                                     dir_to_save,
                                     altura,
                                     anchura,
                                     text_size) {
  R_biomarcadores <- R[, biomarcadores]
  R_obesidad <- R_biomarcadores[obesidad == "Obese", ]
  R_NoObesidad <- R_biomarcadores[obesidad == "No.Obese", ]
  R_females <- R_biomarcadores[grupo == "Control.Women", ]
  R_males <- R_biomarcadores[grupo == "Men", ]
  R_PCOS <- R_biomarcadores[grupo == "PCOS", ]
  R_malenoob <- R_biomarcadores[grupo == "Men" &
                                  obesidad == "Obese", ]
  R_maleoob <- R_biomarcadores[grupo == "Men" &
                                 obesidad == "No.Obese", ]
  R_PCOSob <- R_biomarcadores[grupo == "PCOS" &
                                obesidad == "Obese", ]
  R_PCOSnoob <- R_biomarcadores[grupo == "PCOS" &
                                  obesidad == "No.Obese", ]
  R_femaleOb <- R_biomarcadores[grupo == "Control.Women" &
                                  obesidad == "Obese", ]
  R_femaleNoob <- R_biomarcadores[grupo == "Control.Women" &
                                    obesidad == "No.Obese", ]
  
  
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
    "Control Women",
    "Male",
    "PCOS",
    "Male: No Obese",
    "Male: Obese",
    "PCOS: Obese",
    "PCOS: No Obese",
    "Control Women: Obese",
    "Control Women: No Obese",
    "All subjects"
  )
  
  plots <-  vector("list", length = length(nombres))
  for (gof in 1:length(nombres)) {
    plots[[gof]] <- plot_correlation(grupos_R[[gof]], titulo = nombres[gof], text_size = text_size)
    
    
  }
  
  for (p in 1:length(plots)) {
    ggsave(
      plot = plots[[p]],
      filename = file.path(dir_to_save, paste0("Correlated_", nombres[p], ".jpeg")),
      height = altura,
      width = anchura
    )
    
  }
  
  return(plots)
}


wrapper_aov_PCs <- function(R, biomarcadores, escalado, componente) {
  pcx <- prcomp(R[, biomarcadores], center = T, scale. = escalado)
  sujetos <- colnames(modelo@data$metaboloma$group1)
  PCs <- data.frame(pcx$x,
                    Group = grupo,
                    Obesity = obesidad,
                    Subjects = sujetos)
  model <- lm(PCs[, componente] ~ grupo * obesidad)
  (cte_variance <- summary(lm(sqrt(abs(
    residuals(model)
  )) ~ fitted(model))))
  (normal_Resid <- shapiro.test(resid(model)))
  
  anova_model <- aov(model)
  
  anova_res <- as.data.frame(summary(anova_model)[[1]])
  to_return <- list(
    result = anova_res,
    cte_variance = cte_variance,
    normal_Resid = normal_Resid,
    PCs = PCs
  )
  return(to_return)
}

downstream_pca <- function(R,
                           biomarcadores,
                           componente.,
                           escalado,
                           common_axes,
                           label_y_Noobese,
                           label_y_obese,
                           label_females,
                           label_PCOS,
                           label_men,
                           label_obesity,
                           label_group,
                           dir_to_save,
                           altura,
                           ancho) {
  res_anova_PC <- wrapper_aov_PCs(
    R = R,
    biomarcadores = biomarcadores,
    escalado = escalado,
    componente = componente.
  )
  
  
  F_Value_group <- round(res_anova_PC$result$`F value`[1], 4)
  F_Value_Obesity <- round(res_anova_PC$result$`F value`[2], 4)
  F_Value_interaction <- round(res_anova_PC$result$`F value`[3], 4)
  
  P_value_group <- round(res_anova_PC$result$`Pr(>F)`[1], 4)
  P_value_Obesity <- round(res_anova_PC$result$`Pr(>F)`[2], 4)
  P_value_interaction <- round(res_anova_PC$result$`Pr(>F)`[3], 4)
  
  
  titulo <- paste(
    "Group x Obesity ANOVA: F (Group) = ",
    F_Value_group,
    "\nGroup x Obesity ANOVA: F (Obesity) = ",
    F_Value_Obesity,
    "\nGroup x Obesity ANOVA: F = ",
    F_Value_interaction,
    "\nP-val (Group): = ",
    P_value_group,
    "P-val (Obesity): = ",
    P_value_Obesity,
    "P-val (Interaction): = ",
    P_value_interaction
  )
  
  PCs <- res_anova_PC$PCs
  
  PCs_y <- paste0("PC", componente.)
  
  PCS_obese <- subset(PCs, Obesity == "Obese")
  PCS_No_obese <- subset(PCs, Obesity == "No.Obese")
  
  
  PCS_females <- subset(PCs, Group == "Control.Women")
  PCS_males <- subset(PCs, Group == "Men")
  PCS_PCOS <- subset(PCs, Group == "PCOS")
  
  
  
  
  PCS_obesity.plot <- ggboxplot(
    data = PCs,
    x = "Obesity",
    y = PCs_y,
    color = "Obesity",
    add = c("jitter", "mean")
  ) +
    scale_color_manual(values = c("No.Obese" = "lightblue", "Obese" = "blue")) +
    ylab("") + xlab(PCs_y) +
    theme(legend.position = "none") +
    ylim(-common_axes, common_axes) +
    theme(title = element_text(size = 10), axis.text = element_text(size = 10)) +
    geom_pwc(
      method = "t_test",
      hide.ns = "p",
      label = "p.signif" ,
      y.position = label_obesity,
      symnum.args  = list(
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, Inf),
        symbols = c("****", "***", "**", "*", "\u2020", "ns")
      )
    ) + scale_x_discrete(labels = c("No Obese", "Obese"))
  PCS_obesity.plot <- annotate_figure(PCS_obesity.plot, bottom = text_grob(
    titulo,
    hjust = 0,
    x =
      0.06,
    size = 8
  ))
  
  
  my_comparisons <- list(c("Control.Women", "Men"),
                         c("Control.Women", "PCOS"),
                         c("Men", "PCOS"))
  
  
  (
    PCS_Group.plot <- ggboxplot(
      data = PCs,
      x = "Group",
      y = PCs_y,
      color = "Group",
      add = c("jitter", "mean")
    ) +
      scale_color_manual(
        values = c(
          "Control.Women" = "green2",
          "PCOS" = "orange",
          "Men" = "red3"
        )
      ) +
      ylab("") + xlab(PCs_y) +
      theme(legend.position = "none") +
      ylim(-common_axes, common_axes) +
      theme(title = element_text(size = 10), axis.text = element_text(size = 10)) +
      geom_pwc(
        method = "t_test",
        hide.ns = "p",
        label = "p.signif" ,
        y.position = label_group,
        symnum.args  = list(
          cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, Inf),
          symbols = c("****", "***", "**", "*", "\u2020", "ns")
        )
      ) + scale_x_discrete(labels = c("Control Women", "PCOS", "Men"))
  )
  
  PCS_Group.plot <- annotate_figure(PCS_Group.plot, bottom = text_grob(
    titulo,
    hjust = 0,
    x =
      0.06,
    size = 8
  ))
  
  
  # Crear el gráfico para Obese
  PCS_obese.plot <- ggboxplot(
    data = PCS_obese,
    x = "Group",
    y = PCs_y,
    color = "Group",
    add = c("jitter", "mean")
  ) +
    scale_color_manual(values = c(
      "Control.Women" = "chartreuse4",
      "PCOS" = "darkorange2",
      "Men" = "darkred"
    )) +
    ylab("") + xlab("") +
    theme(legend.position = "none") +
    ylim(-common_axes, common_axes) +
    ggtitle("Obese") +
    theme(title = element_text(size = 10), axis.text = element_text(size = 10)) +
    geom_pwc(
      method = "t_test",
      hide.ns = "p",
      label = "p.signif" ,
      y.position = label_y_obese,
      symnum.args  = list(
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, Inf),
        symbols = c("****", "***", "**", "*", "\u2020", "ns")
      )
    ) + scale_x_discrete(labels = c("Control Women", "PCOS", "Men"))
  
  # Crear el gráfico para No Obese
  PCS_No_obese.plot <- ggboxplot(
    data = PCS_No_obese,
    x = "Group",
    y = PCs_y,
    color = "Group",
    add = c("jitter", "mean")
  ) +
    scale_color_manual(values = c(
      "Control.Women" = "chartreuse",
      "PCOS" = "gold3",
      "Men" = "red"
    )) +
    xlab("") +
    theme(legend.position = "none") +
    ylim(-common_axes, common_axes) +
    ggtitle("No Obese") +
    theme(title = element_text(size = 10), axis.text = element_text(size = 10)) +
    geom_pwc(
      method = "t_test",
      hide.ns = "p",
      label = "p.signif" ,
      y.position = label_y_Noobese,
      symnum.args  = list(
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, Inf),
        symbols = c("****", "***", "**", "*", "\u2020", "ns")
      )
    )  + scale_x_discrete(labels = c("Control Women", "PCOS", "Men"))
  
  
  # Combinar ambos gráficos
  combined_plot_By_Group <- ggarrange(PCS_No_obese.plot, PCS_obese.plot)
  
  combined_plot_By_Group <- annotate_figure(combined_plot_By_Group,
                                            bottom = text_grob(
                                              titulo,
                                              hjust = 0,
                                              x =
                                                0.06,
                                              size = 8
                                            ))
  
  # Crear el gráfico con stat_compare_means
  females_plot <- ggplot(PCS_females, aes(
    x = Obesity,
    y = !!sym(PCs_y),
    color = Obesity
  )) +
    geom_boxplot() +   geom_jitter(width = 0.2, size = 1.5) +  # Añadir jitter
    
    geom_pwc(
      method = "t_test",
      hide.ns = "p",
      label = "p.signif" ,
      y.position = label_females,
      symnum.args  = list(
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, Inf),
        symbols = c("****", "***", "**", "*", "\u2020", "ns")
      )
    )  + scale_color_manual(values = c(
      "No.Obese" = "chartreuse",
      "Obese" = "chartreuse4"
    )) +
    labs(title = "Females", x = "", y = PCs_y) + theme(
      panel.background = element_rect(fill = "white", color = NA),
      # Fondo blanco
      panel.grid = element_blank(),
      # Sin líneas de rejilla
      axis.line = element_line(color = "black"),
      # Ejes en negro
      axis.ticks = element_line(color = "black"),
      # Ticks en negro
      axis.text = element_text(color = "black"),
      # Texto de ejes en negro
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "none"# Fondo de toda la figura en blanco
    ) + ylim(-common_axes, common_axes) + scale_x_discrete(labels = c("No Obese", "Obese"))
  # Mostrar el gráfico
  
  # Mostrar el gráfico
  
  # Crear el gráfico con stat_compare_means
  PCOS_plot <- ggplot(PCS_PCOS, aes(
    x = Obesity,
    y = !!sym(PCs_y),
    color = Obesity
  )) +
    geom_boxplot() + geom_jitter(width = 0.2, size = 1.5) +
    geom_pwc(
      method = "t_test",
      hide.ns = "p",
      label = "p.signif" ,
      y.position = label_PCOS,
      symnum.args  = list(
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, Inf),
        symbols = c("****", "***", "**", "*", "\u2020", "ns")
      )
    )  +
    scale_color_manual(values = c(
      "No.Obese" = "gold3",
      "Obese" = "darkorange2"
    )) + # Colores personalizados) +
    labs(title = "PCOS", x = "", y = "") + theme(
      panel.background = element_rect(fill = "white", color = NA),
      # Fondo blanco
      panel.grid = element_blank(),
      # Sin líneas de rejilla
      axis.line = element_line(color = "black"),
      # Ejes en negro
      axis.ticks = element_line(color = "black"),
      # Ticks en negro
      axis.text = element_text(color = "black"),
      # Texto de ejes en negro
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "none"# Fondo de toda la figura en blanco
    ) + ylim(-common_axes, common_axes) + scale_x_discrete(labels = c("No Obese", "Obese"))
  # Mostrar el gráfico
  
  # Mostrar el gráfico  # Crear el gráfico con stat_compare_means
  MEN_plot <- ggplot(PCS_males, aes(
    x = Obesity,
    y = !!sym(PCs_y),
    color = Obesity
  )) +
    geom_boxplot() + geom_jitter(width = 0.2, size = 1.5) +
    geom_pwc(
      method = "t_test",
      hide.ns = "p",
      label = "p.signif" ,
      y.position = label_men,
      symnum.args  = list(
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, Inf),
        symbols = c("****", "***", "**", "*", "\u2020", "ns")
      )
    )  +
    scale_color_manual(values = c("No.Obese" = "red", "Obese" = "darkred")) + # Colores personalizados) +
    labs(title = "Men", x = "", y = "") + theme(
      panel.background = element_rect(fill = "white", color = NA),
      # Fondo blanco
      panel.grid = element_blank(),
      # Sin líneas de rejilla
      axis.line = element_line(color = "black"),
      # Ejes en negro
      axis.ticks = element_line(color = "black"),
      # Ticks en negro
      axis.text = element_text(color = "black"),
      # Texto de ejes en negro
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "none"# Fondo de toda la figura en blanco
    ) + ylim(-common_axes, common_axes) + scale_x_discrete(labels = c("No Obese", "Obese"))
  # Mostrar el gráfico
  # Mostrar el gráfic
  
  plot_by_obesity <- ggarrange(females_plot,
                               PCOS_plot,
                               MEN_plot,
                               ncol = 3,
                               nrow = 1)
  
  final_plot <- annotate_figure(plot_by_obesity, bottom = text_grob(
    titulo,
    hjust = 0,
    x =
      0.06,
    size = 8
  ))
  
  ggsave(
    plot = combined_plot_By_Group,
    filename = file.path(dir_to_save, paste0("Group_by_Obesity_", PCs_y, ".jpeg")),
    height = altura,
    width = ancho
  )
  ggsave(
    plot = final_plot,
    filename = file.path(dir_to_save, paste0("Obesity_by_Group", PCs_y, ".jpeg")),
    height = altura,
    width = ancho
  )
  ggsave(
    plot = PCS_obesity.plot,
    filename = file.path(dir_to_save, paste0("Obesity_plot", PCs_y, ".jpeg")),
    height = altura,
    width = ancho
  )
  ggsave(
    plot = PCS_Group.plot,
    filename = file.path(dir_to_save, paste0("Group_plot", PCs_y, ".jpeg")),
    height = altura,
    width = ancho
  )
  
  return_plots <- list(group = combined_plot_By_Group, obesity = final_plot)
  return(return_plots)
}

fun_wraper_downstream_PCA <- function(R, biomarcadores, directorios_lista) {
  scaled_dir <- directorios_lista$scaled_dir
  NOT_scaled_dir <- directorios_lista$NOT_scaled_dir
  
  plot_PCA_downstream(
    R,
    escalado = T,
    componente. = 1,
    biomarcadores = biomarcadores,
    dir_to_save = scaled_dir
  )
  
  plot_PCA_downstream(
    R,
    escalado = T,
    componente. = 2,
    biomarcadores = biomarcadores,
    dir_to_save = scaled_dir
  )
  
  plot_PCA_downstream(
    R,
    escalado = F,
    componente. = 1,
    biomarcadores = biomarcadores,
    dir_to_save = NOT_scaled_dir
  )
  
  plot_PCA_downstream(
    R,
    escalado = F,
    componente. = 2,
    biomarcadores = biomarcadores,
    dir_to_save = NOT_scaled_dir
  )
  
  
}

dummy_fun_contrasts <- function(R, grupo, obesidad) {
  grupo_ <- grupo
  # levels(grupo_) <- c("Control Women", "PCOS", "Men")
  obesidad_ <- obesidad
  
  
  
  w_ctrl.group <- (grupo_ != "PCOS")
  w_females.group <- (grupo_ != "Men")
  w_hermas.group <- (grupo_ != "Control.Women")
  
  ctrl.group <- droplevels(grupo_[w_ctrl.group])
  females.group <- droplevels(grupo_[w_females.group])
  hermas.group <- droplevels(grupo_[w_hermas.group])
  
  R.ctrl.group <- R[w_ctrl.group, ]
  R.females.group <- R[w_females.group, ]
  R.hermas.group <- R[w_hermas.group, ]
  
  
  w_ctrl.noob <- (grupo_ != "PCOS") & (obesidad == "No.Obese")
  w_females.noob <- (grupo_ != "Men") & (obesidad == "No.Obese")
  w_hermas.noob <- (grupo_ != "Control.Women") &
    (obesidad == "No.Obese")
  
  
  ctrl.noob <- droplevels(grupo_[w_ctrl.noob])
  females.noob <-  droplevels(grupo_[w_females.noob])
  hermas.noob <-  droplevels(grupo_[w_hermas.noob])
  
  
  R.ctrl.noob <- R[w_ctrl.noob , ]
  R.females.noob <- R[w_females.noob , ]
  R.hermas.noob <- R[w_hermas.noob , ]
  
  w_ctrl.ob <- (grupo_ != "PCOS") & (obesidad == "Obese")
  w_females.ob <- (grupo_ != "Men") & (obesidad == "Obese")
  w_hermas.ob <- (grupo_ != "Control.Women") & (obesidad == "Obese")
  
  ctrl.ob <- droplevels(grupo_[w_ctrl.ob])
  females.ob <-  droplevels(grupo_[w_females.ob])
  hermas.ob <-  droplevels(grupo_[w_hermas.ob])
  
  R.ctrl.ob <- R[w_ctrl.ob, ]
  R.females.ob <- R[w_females.ob, ]
  R.hermas.ob <- R[w_hermas.ob, ]
  
  w_males <- which(grupo_ == "Men")
  w_females <- which(grupo_ == "Control.Women")
  w_pcos <- which(grupo_ == "PCOS")
  
  R.males <- R[w_males, ]
  R.females <- R[w_females, ]
  R.pcos <- R[w_pcos, ]
  males <- (obesidad[w_males])
  females <-  (obesidad[w_females])
  pcos <-  (obesidad[w_pcos])
  
  retorno <- list(
    factores = list(
      control_group = ctrl.group,
      females_group = females.group,
      femalespcos_group = hermas.group,
      control_noob = ctrl.noob,
      females_noob = females.noob,
      hermas_noob = hermas.noob,
      control_ob = ctrl.ob,
      females_ob = females.ob,
      femalespcos_ob = hermas.ob,
      males = males,
      females = females,
      pcos = pcos
    ),
    senales = list(
      control_group = R.ctrl.group,
      females_group = R.females.group,
      femalespcos_group = R.hermas.group,
      control_noob = R.ctrl.noob,
      females_noob = R.females.noob,
      hermas_noob = R.hermas.noob,
      control_ob = R.ctrl.ob,
      females_ob = R.females.ob,
      femalespcos_ob = R.hermas.ob,
      males = R.males,
      females = R.females,
      pcos = R.pcos
    ),
    whiches = list(
      control_group = w_ctrl.group,
      females_group = w_females.group,
      femalespcos_group = w_hermas.group,
      control_noob = w_ctrl.noob,
      females_noob = w_females.noob,
      hermas_noob = w_hermas.noob,
      control_ob = w_ctrl.ob,
      females_ob = w_females.ob,
      femalespcos_ob = w_hermas.ob,
      males = w_males,
      females = w_females,
      pcos = w_pcos
    )
  )
  return(retorno)
  
}




plot_contributions <- function(R,
                               escalado,
                               componente,
                               Titulo,
                               colores,
                               referencia,
                               indice_gof,
                               factor_gof,
                               text.size,
                               common_axes) {
  
  prcomp_obj <- prcomp(R, scale. = escalado)
  
  
  scores <- prcomp_obj$x
  
  sub_scores <- scores[indice_gof, ]
  
  R_sub <- R[indice_gof, ]
  robj <- apply((R_sub), 2, function(x)
    corr.test(x, sub_scores[, componente]))
  r <- unlist(lapply(robj, function(x)
    x$r))
  
  
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
    rotate = TRUE
  ) + xlab("") + ylab("Pearson Correlation Coefficient") + theme(
    axis.text.x = element_text(size = text.size * 1.2),
    title = element_text(size = text.size *
                           1.5),
    axis.text.y = element_text(size = text.size),
    panel.background = element_rect(fill = "white", color = NA),
    # Fondo blanco del panel
    plot.background = element_rect(fill = "white", color = NA),
    # Fondo blanco del gráfico
    panel.grid = element_blank(),
    # Eliminar líneas de la rejilla
    axis.line = element_line(color = "black"),
    # Ejes en negro
    axis.ticks = element_line(color = "black"),
    # Ticks en negro
    axis.text = element_text(color = "black"),
    # Texto de los ejes en negro
    axis.title = element_text(color = "black"),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = text.size * 1.3),
    legend.title = element_blank(),
    axis.title.y = element_text(size = text.size),
    # Títulos de los ejes en negro,
    legend.position = "right"
    
    
  )
  
  
  dfgg <- data.frame(sub_scores , GOF = factor_gof)
  # componente<-"PC1"
  # comparisons <- list(c("Female","Male"))
  
  
  
  componente_ <- paste0("PC", componente)
  
  niveles <- levels(dfgg$GOF)
  nivel_no_ref <- niveles[niveles != referencia]
  print(nivel_no_ref)
  
  my_comparisons = list(c(nivel_no_ref, referencia))
  labels_Aux_ <- c(nivel_no_ref, referencia)
  labels_Aux <- gsub("\\.", " ", labels_Aux_)
  
  labels_box <- setNames(labels_Aux, labels_Aux_)
  colores_ <- setNames(colores, names(colores))

  
  p_score <- ggplot(dfgg, aes(
    x = GOF,
    y = !!sym(componente_),
    color = GOF
  )) +
    geom_boxplot() + geom_jitter(width = 0.2, size = 1.5) +
    geom_pwc(
      method = "t_test",
      hide.ns = "p",
      label = "p.signif" ,
      y.position = 15,
      symnum.args  = list(
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, Inf),
        symbols = c("****", "***", "**", "*", "\u2020", "ns")
      )
    ) +  xlab("") + ggtitle(Titulo) + theme(
    legend.position = "none",
    axis.text.x = element_text(size = text.size * 1.2),
    title = element_text(size = text.size *
                           1.5),
    axis.text.y = element_text(size = text.size),
    panel.background = element_rect(fill = "white", color = NA),
    # Fondo blanco del panel
    plot.background = element_rect(fill = "white", color = NA),
    # Fondo blanco del gráfico
    panel.grid = element_blank(),
    # Eliminar líneas de la rejilla
    axis.line = element_line(color = "black"),
    # Ejes en negro
    axis.ticks = element_line(color = "black"),
    # Ticks en negro
    axis.text = element_text(color = "black"),
    # Texto de los ejes en negro
    axis.title = element_text(color = "black"),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = text.size * 1.3),
    legend.title = element_blank()# Títulos de los ejes en negro
    
  ) + scale_color_manual(values = colores_) + scale_x_discrete(labels = labels_box) + ylim(c(-common_axes,common_axes))
  
  
  
  
  
  p_complete <- ggarrange(p_contrib, p_score)
  
  return(p_complete)
}



wrapper_plot_contribuciones <- function(modelo,
                                        biomarcadores,
                                        escalado,
                                        dir_to_save,
                                        original,
                                        common_axes,
                                        text_size,
                                        altura,
                                        ancho,
                                        grupo = grupo,
                                        obesidad = obesidad) {


  
  #Score Plot PC1:\nObese vs No Obese
  dum_contrastes <- dummy_fun_contrasts(R, grupo, obesidad)
  
  titulos <- vector("list", length = length(names(dum_contrastes$factores)))
  names(titulos) <- names(dum_contrastes$factores)
  titulos$control_group <- list(PC1 = "Score Plot PC1:\nMen vs Control Women", PC2 = "Score Plot PC2:\nMen vs Control Women")
  titulos$females_group <- list(PC1 = "Score Plot PC1:\nPCOS vs Control Women", PC2 = "Score Plot PC2:\nPCOS vs Control Women")
  titulos$femalespcos_group <- list(PC1 = "Score Plot PC1:\nPCOS vs Men", PC2 = "Score Plot PC2:\nPCOS vs Men")
  
  titulos$control_noob <- list(PC1 = "Score Plot PC1:\n(No Obese) Men vs Control Women)", PC2 = "Score Plot PC2:\n(No Obese) Men vs Control Women)")
  titulos$females_noob <- list(PC1 = "Score Plot PC1:\n(No Obese) PCOS vs Control Women", PC2 = "Score Plot PC2:\n(No Obese) PCOS vs Control Women")
  
  titulos$hermas_noob <- list(PC1 = "Score Plot PC1:\n(No Obese) PCOS vs Men", PC2 = "Score Plot PC2:\n(No Obese) PCOS vs Men")
  
  titulos$control_ob <- list(PC1 = "Score Plot PC1:\n(Obese) Men vs Control Women)", PC2 = "Score Plot PC2:\n(Obese) Men vs Control Women)")
  titulos$females_ob <- list(PC1 = "Score Plot PC1:\n(Obese) PCOS vs Control Women", PC2 = "Score Plot PC2:\n(Obese) PCOS vs Control Women")
  
  titulos$femalespcos_ob <- list(PC1 = "Score Plot PC1:\n(Obese) PCOS vs Men", PC2 = "Score Plot PC2:\n(Obese) PCOS vs Men")
  
  titulos$females <- list(PC1 = "Score Plot PC1:\n(Control Women) Obese vs No Obese", PC2 = "Score Plot PC2:\n(Control Women) Obese vs No Obese")
  titulos$males <- list(PC1 = "Score Plot PC1:\n(Men) Obese vs No Obese", PC2 = "Score Plot PC2:\n(Men) Obese vs No Obese")
  titulos$pcos <- list(PC1 = "Score Plot PC1:\n(PCOS) Obese vs No Obese", PC2 = "Score Plot PC2:\n(PCOS) Obese vs No Obese")
  
  if (original == F) {
    R <- reconstruccion(modelo)[, biomarcadores]
    
  } else{
    R <-
      get_senal_original(modelo)[, biomarcadores]
    
  }
  size_barplot <- 10
  
  
  obesidad_PC1 <- plot_contributions(
    R = R,
    escalado = escalado,
    componente = 1,
    Titulo = "Score Plot PC1:\nObese vs No Obese",
    colores =  c("No.Obese" = "lightblue", "Obese" = "blue"),
    referencia = "No.Obese",
    indice_gof =  1:nrow(R),
    factor_gof = obesidad,
    text.size = text_size,
    common_axes = common_axes
  )
  ggsave(
    plot = obesidad_PC1,
    filename = file.path(dir_to_save, "PC1_Obese_vs_No_Obese.jpeg"),
    height = altura,
    width = ancho
  )
  
  obesidad_PC2 <- plot_contributions(
    R = R,
    escalado = escalado,
    componente = 2,
    Titulo = "Score Plot PC2:\nObese vs No Obese",
    colores =  c("No.Obese" = "lightblue", "Obese" = "blue"),
    referencia = "No.Obese",
    indice_gof =  1:nrow(R),
    factor_gof = obesidad,
    text.size = text_size,
    common_axes = common_axes
  )
  ggsave(
    plot = obesidad_PC2,
    filename = file.path(dir_to_save, "PC2_Obese_vs_No_Obese.jpeg"),
    height = altura,
    width = ancho
  )
  # colores <- ifelse(plotdf$subjects == "Control Women: No Obese", "chartreuse", NA)
  # colores <- ifelse(plotdf$subjects == "PCOS: No Obese", "gold3", colores)
  # colores <- ifelse(plotdf$subjects == "Men: No Obese", "red", colores)
  # colores <- ifelse(plotdf$subjects == "Control Women: Obese" ,
  #                   "chartreuse4",
  #                   colores)
  # colores <- ifelse(plotdf$subjects == "PCOS: Obese", "darkorange2", colores)
  # colores <- ifelse(plotdf$subjects == "Men: Obese", "darkred", colores)
  
  
  referencias <- vector("list", length = length(names(dum_contrastes$factores)))
  names(referencias) <- names(dum_contrastes$factores)
  referencias$control_group <- "Control.Women"
  referencias$females_group <- "Control.Women"
  referencias$femalespcos_group <- "Men"
  referencias$control_noob <- "Control.Women"
  referencias$females_noob <- "Control.Women"
  referencias$hermas_noob <- "Men"
  referencias$control_ob <- "Control.Women"
  referencias$females_ob <- "Control.Women"
  referencias$femalespcos_ob <- "Men"
  referencias$females <- "No.Obese"
  referencias$males <- "No.Obese"
  referencias$pcos <- "No.Obese"
  
  plots_PC1s <- vector("list", length = length(names(dum_contrastes$factores)))
  names(plots_PC1s) <- names(dum_contrastes$factores)
  
  colores <- vector("list", length = length(names(dum_contrastes$factores)))
  names(colores) <- names(dum_contrastes$factores)
  
  colores$control_group <- c("Control.Women" = "green2", "Men" = "red3")
  colores$females_group <- c("Control.Women" = "green2",
                             "PCOS" = "orange")
  colores$femalespcos_group <- c("Men" = "red3", "PCOS" = "orange")
  
  colores$control_noob <- c("Control.Women" = "chartreuse",
                            "Men" = "red")
  colores$females_noob <- c("Control.Women" = "chartreuse",
                            "PCOS" = "gold3")
  colores$hermas_noob <- c("Men" = "red", "PCOS" = "gold3")
  
  
  colores$control_ob <- c("Control.Women" = "chartreuse4",
                          "Men" = "darkred")
  colores$females_ob <- c("Control.Women" = "chartreuse4",
                          "PCOS" = "darkorange2")
  colores$femalespcos_ob <- c("Men" = "darkred", "PCOS" = "darkorange2")
  
  colores$males <- c("No.Obese" = "red", "Obese" = "darkred")
  colores$females <-  c("No.Obese" = "chartreuse",
                        "Obese" = "chartreuse4")
  
  colores$pcos <-  c("No.Obese" = "gold3",
                     "Obese" = "darkorange2")
  
  plots_PC1s <- vector("list", length = length(names(dum_contrastes$factores)))
  names(plots_PC1s) <- names(dum_contrastes$factores)
  
  plots_PC2s <- vector("list", length = length(names(dum_contrastes$factores)))
  names(plots_PC2s) <- names(dum_contrastes$factores)
  for (pl in 1:length(names(dum_contrastes$factores))) {
    plot_PC1 <- plot_contributions(
      R = R,
      escalado = escalado,
      componente = 1,
      Titulo = titulos[[pl]]$PC1,
      colores = colores[[pl]],
      referencia = referencias[[pl]],
      indice_gof = dum_contrastes$whiches[[pl]],
      factor_gof = dum_contrastes$factores[[pl]],
      text.size = text_size,
      common_axes = common_axes
    )
    plot_PC2 <- plot_contributions(
      R = R,
      escalado = escalado,
      componente = 2,
      Titulo = titulos[[pl]]$PC2,
      colores = colores[[pl]],
      referencia =  referencias[[pl]],
      indice_gof = dum_contrastes$whiches[[pl]],
      factor_gof = dum_contrastes$factores[[pl]],
      text.size = text_size,
      common_axes = common_axes
    )
    plots_PC1s[[pl]] <- plot_PC1
    plots_PC2s[[pl]] <- plot_PC2
  }
  
  titulos_to_Save <- c(
    "Male_vs_Female",
    "PCOS_vs_Female",
    "PCOS_vs_Male",
    "NoObese_Male_vs_Female",
    "NoObese_PCOS_vs_Female",
    "NoObese_PCOS_vs_Male",
    "Obese_Male_vs_Female",
    "Obese_PCOS_vs_Female",
    "Obese_PCOS_vs_Male",
    "Males_Obese_vs_NoObese",
    "Females_Obese_vs_NoObese",
    "PCOS_Obese_vs_NoObese"
  )
  
  
  dir_to_save_pl_PC1 <- file.path(dir_to_save, "PC1")
  dir_to_save_pl_PC2 <- file.path(dir_to_save, "PC2")
  
  if (!dir.exists(dir_to_save_pl_PC1)) {
    dir.create(dir_to_save_pl_PC1)
    
  }
  if (!dir.exists(dir_to_save_pl_PC2)) {
    dir.create(dir_to_save_pl_PC2)
    
  }
  
  for (pl in seq_along(titulos_to_Save)) {
    ggsave(
      plot = plots_PC1s[[pl]],
      filename = file.path(dir_to_save_pl_PC1, paste0(titulos_to_Save[pl], ".jpeg")),
      height = altura,
      width = ancho
    )
    
    ggsave(
      plot = plots_PC2s[[pl]],
      filename = file.path(dir_to_save_pl_PC2, paste0(titulos_to_Save[pl], ".jpeg")),
      height = altura,
      width = ancho
    )
    
  }
  
  return(list(PC1 = plots_PC1s, PC2 = plots_PC2s))
}

dir_create <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = T)
  }
  
}

get_senal_original <- function(modelo) {
  so <- bind_cols(as.data.frame(scale(t(
    Reduce(
      "rbind",
      list(
        modelo@data$metaboloma[[1]],
        modelo@data$metagenoma[[1]],
        modelo@data$ip[[1]]
      )
    )
  ))))
  return(so)
}


plot_contributions_general <- function(R,
                                       biomarcadores,
                                       escalado,
                                       componente,
                                       text_size,
                                       dir_to_save,
                                       altura,
                                       ancho) {
  R <- R[, biomarcadores]
  prcomp_obj <- prcomp(R, scale. = escalado)
  
  
  scores <- prcomp_obj$x
  
  sub_scores <- scores
  
  R_sub <- R
  robj <- apply((R_sub), 2, function(x)
    corr.test(x, sub_scores[, componente]))
  r <- unlist(lapply(robj, function(x)
    x$r))
  
  
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
  
  
  
  
  
  p_contrib_cor <- ggbarplot(
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
    rotate = TRUE
  ) + xlab("") + ylab("Pearson Correlation Coefficient") + theme(
    axis.text.x = element_text(size = text_size * 1.2),
    title = element_text(size = text_size *
                           1.5),
    axis.text.y = element_text(size = text_size),
    panel.background = element_rect(fill = "white", color = NA),
    # Fondo blanco del panel
    plot.background = element_rect(fill = "white", color = NA),
    # Fondo blanco del gráfico
    panel.grid = element_blank(),
    # Eliminar líneas de la rejilla
    axis.line = element_line(color = "black"),
    # Ejes en negro
    axis.ticks = element_line(color = "black"),
    # Ticks en negro
    axis.text = element_text(color = "black"),
    # Texto de los ejes en negro
    axis.title = element_text(color = "black"),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = text_size * 1.3),
    legend.title = element_blank(),
    axis.title.y = element_text(size = text_size),
    # Títulos de los ejes en negro,
    legend.position = "right"
    
  )+ ggtitle(paste("Principal Component", componente))
  
  contribuciones <- factoextra::fviz_contrib(prcomp_obj, choice = "var", axes = componente)
  dd <- factoextra::facto_summarize(prcomp_obj,
                                    element = "var",
                                    result = "contrib",
                                    axes = componente)
  contrib <- dd$contrib
  names(contrib) <- rownames(dd)
  
  # expected Average contribution
  theo_contrib <- 100 / length(contrib)
  
  datos_Contrib <- contribuciones$data[order(contribuciones$data$contrib, decreasing = T), ]
  colnames(datos_Contrib)[1] <- "names"
  contribuciones_to_plot <- inner_join(toplot, datos_Contrib, by = "names")
  
  contribuciones_to_plot <- contribuciones_to_plot[order(contribuciones_to_plot$contrib, decreasing = T), ]
  
  
  p_contrib <- ggbarplot(
    contribuciones_to_plot,
    x = "names",
    y = "contrib",
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
    rotate = F
  ) + xlab("") + ylab("Contributions (%)") + theme(
    axis.text.x = element_text(size = text_size * 1.2),
    title = element_text(size = text_size *
                           1.5),
    axis.text.y = element_text(size = text_size),
    panel.background = element_rect(fill = "white", color = NA),
    # Fondo blanco del panel
    plot.background = element_rect(fill = "white", color = NA),
    # Fondo blanco del gráfico
    panel.grid = element_blank(),
    # Eliminar líneas de la rejilla
    axis.line = element_line(color = "black"),
    # Ejes en negro
    axis.ticks = element_line(color = "black"),
    # Ticks en negro
    axis.text = element_text(color = "black"),
    # Texto de los ejes en negro
    axis.title = element_text(color = "black"),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = text_size * 1.3),
    legend.title = element_blank(),
    axis.title.y = element_text(size = text_size),
    # Títulos de los ejes en negro,
    legend.position = "right"
    
  ) + ggtitle(paste("Principal Component", componente)) +    geom_hline(yintercept = theo_contrib,
                                                                        linetype = 2,
                                                                        color = "red")
  
  
  
  return_plots <- list(contribuciones = p_contrib,
                       correlations_contrib = p_contrib_cor)
  
  ggsave(
    filename = file.path(
      dir_to_save,
      paste0("Contribution_plot_PC", componente, "_.jpeg")
    ),
    plot = p_contrib,
    height = altura / 2,
    width = ancho*1.5
  )
  
  ggsave(
    filename = file.path(
      dir_to_save,
      paste0("Contribution_correlated_plot_PC", componente, "_.jpeg")
    ),
    plot = p_contrib_cor,
    height = altura,
    width = ancho
  )
  
  return(return_plots)
}
###===Cargamos modelo ====
modelos <- lapply(list.files("./modelos/14_10_24_With_OUTLIERS/", full.names = T),
                  load_model)
modelo <- MOFA2::select_model(modelos, plot = T)
senal_original <- get_senal_original(modelo)


get_clinical <- function(modelo) {
  feat_fuera <- c("group", "sample")
  idx <- sapply(feat_fuera, function(x)
    grep(x, colnames(modelo@samples_metadata)))
  
  ## reorder clinical variables
  factores.df <- as.data.frame(MOFA2::get_expectations(modelo, "Z", as.data.frame = F)[[1]])
  
  matrix1 <- modelo@samples_metadata[, -idx]
  matrix1 <- as.data.frame(matrix1)
  matrix1 <- matrix1[, c(
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
  )]
  return(matrix1)
}

senal_reconstruida <- reconstruccion(modelo)





###===Graficos preliminares====

directorio <- "./resultados/30_1-12-24"



directorio_preliminar <- file.path(directorio, "plots_preliminares")

dir_create(directorio_preliminar)

plot_preliminar <- fun_plot_dummy_preliminar(
  modelo,
  dir_to_save = directorio_preliminar,
  altura = 15,
  ancho = 20,
  size_var = 8,
  size_Axis = 18
)

# Mostrar el gráfico combinado
# Mostrar el gráfico combinado
###====factores latentes=====

directorio_latentes <- file.path(directorio, "latent_variables_plots")
dir_create(directorio_latentes)

####====factores contrastes univariate=====

factores.df <- as.data.frame(MOFA2::get_expectations(modelo, "Z", as.data.frame = F)[[1]])
plot_contrastes_latent_factors <- plot_contrastes_wrapper(
  factores.df,
  dir_to_save = directorio_latentes,
  name_to_save = "contrasts_latent_factors",
  size_text = 12,
  altura = 10,
  ancho = 8
)

# ####====latent correlated clinical=====

plot_clinical_correlated_latent_factors <- fun_dummy_plot_clinical_latent(
  modelo,
  dir_to_save =
    directorio_latentes ,
  name_to_save = "latent_factors_clinical_correlation",
  ancho = 18,
  altura = 8,
  factor_ = 20,
  text_size = 15
  
)


###====Feature selection #####

### factores a tener en cuenta para la extraccion

directorio_feature_selection <- file.path(directorio
                                          , "feature_selection_plots")
dir_create(directorio_feature_selection)
threshold <- 0.6
(p_metaboloma <- plot_weights(modelo, view = 1, factors = 1) + geom_vline(xintercept = c(-threshold, threshold)))

ggsave(
  plot = p_metaboloma,
  filename = file.path(
    directorio_feature_selection,
    paste0("weights_metaboloma", ".jpeg")
  ),
  width = 18,
  height = 12
)

(
  p_metagenoma <- plot_weights(modelo, view = 2, factors = 3:5) + geom_vline(xintercept = c(-threshold, threshold))
)
ggsave(
  plot = p_metagenoma,
  filename = file.path(
    directorio_feature_selection,
    paste0("weights_microbiome", ".jpeg")
  ),
  width = 18,
  height = 12
)

(p_proteins <- plot_weights(modelo, view = 3, factors = 2) + geom_vline(xintercept = c(-threshold, threshold)))
ggsave(
  plot = p_proteins,
  filename = file.path(
    directorio_feature_selection,
    paste0("weights_proteins", ".jpeg")
  ),
  width = 18,
  height = 12
)


fs <- list(metaboloma = 1,
           ip = 2,
           metagenoma = 3:5)

biomarcadores_mon <- wrapper_feature_selection_method_mon(
  modelo = modelo,
  threshold = 0.8,
  thresh_pval = 0.001,
  fs = fs
)
biomarcadores_MOFA <- wrapper_feature_selection_method_MOFA(modelo = modelo,
                                                            threshold = 0.6,
                                                            fs = fs)

### VAMOS a crear 5 subdirectorios de resultados
#### aproximacion Edmond  Reconstrucion y Original
#### aproximacion MOFA Reconstrucion y Original

directorio_mon <- file.path(directorio, "aproximacion_Edmond")
directorio_MOFA <- file.path(directorio, "aproximacion_MOFA")
directorio_original <- file.path(directorio, "original")
directorio_original_mon <- file.path(directorio_original, "aproximacion_Edmond")
directorio_original_MOFA <- file.path(directorio_original, "aproximacion_MOFA")
if (!dir.exists(directorio_mon)) {
  dir.create(directorio_mon)
}
if (!dir.exists(directorio_MOFA)) {
  dir.create(directorio_MOFA)
}
if (!dir.exists(directorio_original)) {
  dir.create(directorio_original)
}
if (!dir.exists(directorio_original_mon)) {
  dir.create(directorio_original_mon)
}
if (!dir.exists(directorio_original_MOFA)) {
  dir.create(directorio_original_MOFA)
}
# ##========More latent factors weights=====


directorio_weights_plots_mon <- file.path(directorio_mon, "weight_latent_Factors_plots")

if (!dir.exists(directorio_weights_plots_mon)) {
  dir.create(directorio_weights_plots_mon)
}

directorio_weights_plots_MOFA <- file.path(directorio_MOFA, "weight_latent_Factors_plots")

if (!dir.exists(directorio_weights_plots_MOFA)) {
  dir.create(directorio_weights_plots_MOFA)
}


modelo_weights_plots_mon <- dummy_plot_weight_factor(modelo = modelo,
                                                     biomarcadores = biomarcadores_mon,
                                                     dir_to_save = directorio_weights_plots_mon)







modelo_weights_plots_MOFA <- dummy_plot_weight_factor(modelo = modelo,
                                                      biomarcadores = biomarcadores_MOFA,
                                                      dir_to_save = directorio_weights_plots_MOFA)
# ###======Univariante reconstruccion====#####

R <- reconstruccion(modelo)
R_MON <- R[, biomarcadores_mon]
R_MOFA <- R[, biomarcadores_MOFA]
O_MON <- senal_original[, biomarcadores_mon]
O_MOFA <- senal_original[, biomarcadores_MOFA]

altura <- 10
ancho <- 16
size_text <- 14
plot_contrastes_wrapper(
  R_MON,
  directorio_mon,
  "Reconstructed Signal (Method 2): Contrasts",
  altura = altura ,
  ancho = ancho,
  size_text = size_text
)
plot_contrastes_wrapper(
  R_MOFA,
  directorio_MOFA,
  "Reconstructed Signal (Method MOFA): Contrasts",
  altura = altura ,
  ancho = ancho,
  size_text = size_text
)
plot_contrastes_wrapper(
  O_MON,
  directorio_original_mon,
  "Reconstructed Signal (Method 2): Contrasts",
  altura = altura ,
  ancho = ancho,
  size_text = size_text
)
plot_contrastes_wrapper(
  O_MOFA,
  directorio_original_MOFA,
  "Reconstructed Signal (Method MOFA): Contrasts",
  altura = altura ,
  ancho = ancho,
  size_text = size_text
)

###=========Bivariante reconstruccion=====

directorio_bivariante_mon <- file.path(directorio_mon, "analisis_bivariante")
directorio_bivariante_MOFA <- file.path(directorio_MOFA, "analisis_bivariante")
directorio_original_mon_bivariante <- file.path(directorio_original_mon, "analisis_bivariante")
directorio_original_MOFA_bivariante <- file.path(directorio_original_MOFA, "analisis_bivariante")

if (!dir.exists(directorio_bivariante_mon)) {
  dir.create(directorio_bivariante_mon)
}
if (!dir.exists(directorio_bivariante_MOFA)) {
  dir.create(directorio_bivariante_MOFA)
}
if (!dir.exists(directorio_original_mon_bivariante)) {
  dir.create(directorio_original_mon_bivariante)
}
if (!dir.exists(directorio_original_MOFA_bivariante)) {
  dir.create(directorio_original_MOFA_bivariante)
}

altura <- 23
anchura <- 23
size_text <- 20
clinical_Vars <- get_clinical(modelo)
colnames.clinical_Vars <- colnames(clinical_Vars)
R <- reconstruccion(modelo)
plot_correlation_wrapper(
  R = cbind(R, clinical_Vars),
  biomarcadores = c(biomarcadores_mon, colnames.clinical_Vars),
  dir_to_save = directorio_bivariante_mon,
  altura = altura,
  anchura = anchura,
  text_size = size_text
)

plot_correlation_wrapper(
  R = cbind(senal_original, clinical_Vars),
  biomarcadores = c(biomarcadores_mon, colnames.clinical_Vars),
  dir_to_save = directorio_original_mon_bivariante,
  altura = altura,
  anchura = anchura,
  text_size = size_text
)
plot_correlation_wrapper(
  R = cbind(R, clinical_Vars),
  biomarcadores = c(biomarcadores_MOFA, colnames.clinical_Vars),
  dir_to_save = directorio_bivariante_MOFA,
  altura = altura,
  anchura = anchura,
  text_size = size_text
)
plot_correlation_wrapper(
  R = cbind(senal_original, clinical_Vars),
  biomarcadores = c(biomarcadores_MOFA, colnames.clinical_Vars),
  dir_to_save = directorio_original_MOFA_bivariante,
  altura = altura,
  anchura = anchura,
  text_size = size_text
)

##======PCA dirs=======
directorio_multivariante_mon <- file.path(directorio_mon, "analisis_multivariante", "Scaled")
directorio_multivariante_mon_Not_Scaled <- file.path(directorio_mon, "analisis_multivariante", "Not_Scaled")

directorio_multivariante_MOFA <- file.path(directorio_MOFA, "analisis_multivariante", "Scaled")
directorio_multivariante_MOFA_Not_Scaled <- file.path(directorio_MOFA, "analisis_multivariante", "Not_Scaled")

directorio_original_mon_multivariante <- file.path(directorio_original_mon, "analisis_multivariante", "Scaled")
directorio_original_mon_multivariante_Not_Scaled  <- file.path(directorio_original_mon,
                                                               "analisis_multivariante",
                                                               "Not_Scaled")

directorio_original_MOFA_multivariante <- file.path(directorio_original_MOFA, "analisis_multivariante", "Scaled")
directorio_original_MOFA_multivariante_Not_Scaled  <- file.path(directorio_original_MOFA,
                                                                "analisis_multivariante",
                                                                "Not_Scaled")




if (!dir.exists(directorio_multivariante_mon)) {
  dir.create(directorio_multivariante_mon, recursive = T)
}
if (!dir.exists(directorio_multivariante_MOFA)) {
  dir.create(directorio_multivariante_MOFA, recursive = T)
}
if (!dir.exists(directorio_original_mon_multivariante)) {
  dir.create(directorio_original_mon_multivariante, recursive = T)
}
if (!dir.exists(directorio_original_MOFA_multivariante)) {
  dir.create(directorio_original_MOFA_multivariante, recursive = T)
}

if (!dir.exists(directorio_multivariante_mon_Not_Scaled)) {
  dir.create(directorio_multivariante_mon_Not_Scaled, recursive = T)
}
if (!dir.exists(directorio_multivariante_MOFA_Not_Scaled)) {
  dir.create(directorio_multivariante_MOFA_Not_Scaled, recursive = T)
}
if (!dir.exists(directorio_original_mon_multivariante_Not_Scaled)) {
  dir.create(directorio_original_mon_multivariante_Not_Scaled,
             recursive = T)
}
if (!dir.exists(directorio_original_MOFA_multivariante_Not_Scaled)) {
  dir.create(directorio_original_MOFA_multivariante_Not_Scaled,
             recursive = T)
}


####======PCA perse=======
params <- list(
  text_size = 17,
  altura = 12,
  ancho = 15,
  common_axes = 15,
  scale_labs = 0.9,
  legend_size = 0.9,
  size_points = 5
)

(
  plot_multivariante_mon_not_scaled <- fun_pca(
    R = reconstruccion(modelo),
    biomarcadores = biomarcadores_mon,
    titulo = "PCA: Not Scaled",
    dir_to_save = directorio_multivariante_mon_Not_Scaled,
    name_to_save = "PCA_NOT_SCALED",
    scalar = F,
    text_size = params$text_size,
    altura = params$altura,
    ancho = params$ancho,
    common_axes = params$common_axes,
    scale_labs = params$scale_labs,
    legend_size = params$legend_size,
    size_points = params$size_points
  )
)

(
  plot_multivariante_mon_scaled <- fun_pca(
    R =  reconstruccion(modelo),
    biomarcadores = biomarcadores_mon,
    titulo = "PCA: Scaled",
    dir_to_save = directorio_multivariante_mon,
    name_to_save = "PCA_SCALED",
    scalar = T,
    text_size = params$text_size,
    altura = params$altura,
    ancho = params$ancho,
    common_axes = params$common_axes,
    scale_labs = params$scale_labs,
    legend_size = params$legend_size,
    size_points = params$size_points
    
  )
)


(
  plot_multivariante_MOFA_not_Scaled <- fun_pca(
    R =  reconstruccion(modelo),
    biomarcadores = biomarcadores_MOFA,
    titulo = "PCA: Not Scaled",
    dir_to_save = directorio_multivariante_MOFA_Not_Scaled,
    name_to_save = "PCA_NOT_SCALED",
    scalar = F,
    text_size = params$text_size,
    altura = params$altura,
    ancho = params$ancho,
    common_axes = params$common_axes,
    scale_labs = params$scale_labs,
    legend_size = params$legend_size,
    size_points = params$size_points
  )
)

(
  plot_multivariante_MOFA_Scaled <- fun_pca(
    R =  reconstruccion(modelo),
    biomarcadores = biomarcadores_MOFA,
    titulo = "PCA: Scaled",
    dir_to_save = directorio_multivariante_MOFA,
    name_to_save = "PCA_SCALED",
    scalar = T,
    text_size = params$text_size,
    altura = params$altura,
    ancho = params$ancho,
    common_axes = params$common_axes,
    scale_labs = params$scale_labs,
    legend_size = params$legend_size,
    size_points = params$size_points
  )
)


(
  plot_multivariante_mon_not_Scaled_original <- fun_pca(
    R = get_senal_original(modelo),
    biomarcadores = biomarcadores_mon,
    titulo = "PCA: Not Scaled",
    dir_to_save = directorio_original_mon_multivariante_Not_Scaled,
    name_to_save = "PCA_NOT_SCALED",
    scalar = F,
    text_size = params$text_size,
    altura = params$altura,
    ancho = params$ancho,
    common_axes = params$common_axes,
    scale_labs = params$scale_labs,
    legend_size = params$legend_size,
    size_points = params$size_points
  )
)

(
  plot_multivariante_mon_Scaled_original <- fun_pca(
    R = get_senal_original(modelo),
    biomarcadores = biomarcadores_mon,
    titulo = "PCA: Scaled",
    dir_to_save = directorio_original_mon_multivariante,
    name_to_save = "PCA_SCALED",
    scalar = T,
    text_size = params$text_size,
    altura = params$altura,
    ancho = params$ancho,
    common_axes = params$common_axes,
    scale_labs = params$scale_labs,
    legend_size = params$legend_size,
    size_points = params$size_points
  )
)

(
  plot_multivariante_MOFA_not_Scaled_original <- fun_pca(
    R = get_senal_original(modelo),
    biomarcadores = biomarcadores_MOFA,
    titulo = "PCA: Not Scaled",
    dir_to_save = directorio_original_MOFA_multivariante_Not_Scaled,
    name_to_save = "PCA_NOT_SCALED",
    scalar = F,
    text_size = params$text_size,
    altura = params$altura,
    ancho = params$ancho,
    common_axes = params$common_axes,
    scale_labs = params$scale_labs,
    legend_size = params$legend_size,
    size_points = params$size_points
  )
)

(
  plot_multivariante_MOFA_Scaled_original <- fun_pca(
    R = get_senal_original(modelo),
    biomarcadores = biomarcadores_MOFA,
    titulo = "PCA: Scaled",
    dir_to_save = directorio_original_MOFA_multivariante,
    name_to_save = "PCA_SCALED",
    scalar = T,
    text_size = params$text_size,
    altura = params$altura,
    ancho = params$ancho,
    common_axes = params$common_axes,
    scale_labs = params$scale_labs,
    legend_size = params$legend_size,
    size_points = params$size_points
  )
)



####====Downstream PCA=====


#####===edmond approx=====

R <- reconstruccion(modelo)
downstream_pca(
  R,
  biomarcadores = biomarcadores_mon,
  componente. = 1,
  escalado = T,
  common_axes = 18,
  label_y_Noobese = c(14, 12, 10),
  label_y_obese = c(12, 8, 5),
  label_females = 13,
  label_PCOS = 13,
  label_men = 12,
  dir_to_save = directorio_multivariante_mon,
  altura = 6,
  ancho = 8,
  label_group = c(14, 13, 10),
  label_obesity = c(14, 14)
)

downstream_pca(
  R,
  biomarcadores = biomarcadores_mon,
  componente. = 2,
  escalado = T,
  common_axes = 18,
  label_y_Noobese = c(17, 13, 13),
  label_y_obese = c(14, 10, 10),
  dir_to_save = directorio_multivariante_mon,
  label_females = 15,
  label_PCOS = 15,
  label_men = 15,
  altura = 6,
  ancho = 8,
  label_group = c(15, 13, 10),
  label_obesity = c(14, 14)
)

downstream_pca(
  R,
  biomarcadores = biomarcadores_mon,
  componente. = 1,
  escalado = F,
  common_axes = 18,
  label_y_Noobese = c(12, 16, 10),
  label_y_obese = c(12, 8, 5),
  label_females = 12,
  label_PCOS = 12,
  label_men = 12,
  dir_to_save = directorio_multivariante_mon_Not_Scaled,
  altura = 6,
  ancho = 8,
  label_group = c(16, 13, 10),
  label_obesity = c(14, 14)
)

downstream_pca(
  R,
  biomarcadores = biomarcadores_mon,
  componente. = 2,
  escalado = F,
  common_axes = 18,
  label_y_Noobese = c(10, 11, 13),
  label_y_obese = c(10, 8, 10),
  dir_to_save = directorio_multivariante_mon_Not_Scaled,
  label_females = 15,
  label_PCOS = 15,
  label_men = 15,
  altura = 6,
  ancho = 8,
  label_group = c(16, 13, 10),
  label_obesity = c(14, 14)
)

#####====MOFA approex =====

downstream_pca(
  R,
  biomarcadores = biomarcadores_MOFA,
  componente. = 1,
  escalado = T,
  common_axes = 18,
  label_y_Noobese = c(10, 15, 8),
  label_y_obese = c(10, 15, 5),
  label_females = 12,
  label_PCOS = 12,
  label_men = 12,
  dir_to_save = directorio_multivariante_MOFA,
  altura = 6,
  ancho = 8,
  label_group = c(16, 13, 10),
  label_obesity = c(14, 14)
)

downstream_pca(
  R,
  biomarcadores = biomarcadores_MOFA,
  componente. = 2,
  escalado = T,
  common_axes = 18,
  label_y_Noobese = c(17, 11, 13),
  label_y_obese = c(14, 8, 10),
  dir_to_save = directorio_multivariante_MOFA,
  label_females = 15,
  label_PCOS = 15,
  label_men = 15,
  altura = 6,
  ancho = 8,
  label_group = c(16, 13, 10),
  label_obesity = c(14, 14)
)

downstream_pca(
  R,
  biomarcadores = biomarcadores_MOFA,
  componente. = 1,
  escalado = F,
  common_axes = 18,
  label_y_Noobese = c(12, 16, 10),
  label_y_obese = c(12, 16, 5),
  label_females = 12,
  label_PCOS = 12,
  label_men = 12,
  dir_to_save = directorio_multivariante_MOFA_Not_Scaled,
  altura = 6,
  ancho = 8,
  label_group = c(16, 13, 10),
  label_obesity = c(14, 14)
)

downstream_pca(
  R,
  biomarcadores = biomarcadores_MOFA,
  componente. = 2,
  escalado = F,
  common_axes = 17,
  label_y_Noobese = c(16, 11, 13),
  label_y_obese = c(14, 8, 10),
  dir_to_save = directorio_multivariante_MOFA_Not_Scaled,
  label_females = 15,
  label_PCOS = 15,
  label_men = 15,
  altura = 6,
  ancho = 8,
  label_group = c(16, 13, 10),
  label_obesity = c(14, 14)
)





#####=== Original edmond approx=====

S <- get_senal_original(modelo)
downstream_pca(
  S,
  biomarcadores = biomarcadores_mon,
  componente. = 1,
  escalado = T,
  common_axes = 18,
  label_y_Noobese = c(12, 17, 10),
  label_y_obese = c(12, 8, 5),
  label_females = 12,
  label_PCOS = 12,
  label_men = 12,
  dir_to_save = directorio_original_mon_multivariante,
  altura = 6,
  ancho = 8,
  label_group = c(16, 13, 10),
  label_obesity = c(14, 14)
)

downstream_pca(
  S,
  biomarcadores = biomarcadores_mon,
  componente. = 2,
  escalado = T,
  common_axes = 18,
  label_y_Noobese = c(14, 9, 13),
  label_y_obese = c(14, 8, 10),
  dir_to_save = directorio_original_mon_multivariante,
  label_females = 15,
  label_PCOS = 15,
  label_men = 15,
  altura = 6,
  ancho = 8,
  label_group = c(16, 13, 10),
  label_obesity = c(14, 14)
)

downstream_pca(
  S,
  biomarcadores = biomarcadores_mon,
  componente. = 1,
  escalado = F,
  common_axes = 18,
  label_y_Noobese = c(12, 16, 10),
  label_y_obese = c(12, 8, 5),
  label_females = 12,
  label_PCOS = 12,
  label_men = 12,
  dir_to_save = directorio_original_mon_multivariante_Not_Scaled,
  altura = 6,
  ancho = 8,
  label_group = c(16, 13, 10),
  label_obesity = c(14, 14)
)

downstream_pca(
  S,
  biomarcadores = biomarcadores_mon,
  componente. = 2,
  escalado = F,
  common_axes = 18,
  label_y_Noobese = c(16, 11, 13),
  label_y_obese = c(14, 8, 10),
  dir_to_save = directorio_original_mon_multivariante_Not_Scaled,
  label_females = 15,
  label_PCOS = 15,
  label_men = 15,
  altura = 6,
  ancho = 8,
  label_group = c(16, 13, 10),
  label_obesity = c(14, 14)
)

#####==== Original MOFA approex =====

downstream_pca(
  S,
  biomarcadores = biomarcadores_MOFA,
  componente. = 1,
  escalado = T,
  common_axes = 18,
  label_y_Noobese = c(8, 16, 10),
  label_y_obese = c(10, 16, 5),
  label_females = 12,
  label_PCOS = 12,
  label_men = 12,
  dir_to_save = directorio_original_MOFA_multivariante,
  altura = 6,
  ancho = 8,
  label_group = c(16, 13, 10),
  label_obesity = c(14, 14)
)

downstream_pca(
  S,
  biomarcadores = biomarcadores_MOFA,
  componente. = 2,
  escalado = T,
  common_axes = 18,
  label_y_Noobese = c(16, 11, 13),
  label_y_obese = c(14, 8, 10),
  dir_to_save = directorio_original_MOFA_multivariante,
  label_females = 15,
  label_PCOS = 15,
  label_men = 15,
  altura = 6,
  ancho = 8,
  label_group = c(16, 13, 10),
  label_obesity = c(14, 14)
)

downstream_pca(
  S,
  biomarcadores = biomarcadores_MOFA,
  componente. = 1,
  escalado = F,
  common_axes = 18,
  label_y_Noobese = c(12, 16, 10),
  label_y_obese = c(8, 12, 5),
  label_females = 12,
  label_PCOS = 12,
  label_men = 12,
  dir_to_save = directorio_original_MOFA_multivariante_Not_Scaled,
  altura = 6,
  ancho = 8,
  label_group = c(16, 13, 10),
  label_obesity = c(14, 14)
)

downstream_pca(
  S,
  biomarcadores = biomarcadores_MOFA,
  componente. = 2,
  escalado = F,
  common_axes = 18,
  label_y_Noobese = c(17, 11, 13),
  label_y_obese = c(14, 8, 10),
  dir_to_save = directorio_original_MOFA_multivariante_Not_Scaled,
  label_females = 15,
  label_PCOS = 15,
  label_men = 15,
  altura = 6,
  ancho = 8,
  label_group = c(16, 13, 10),
  label_obesity = c(14, 14)
)


####==CONTRIBUCIONES====

plots_mon_scaled <- wrapper_plot_contribuciones(
  modelo,
  biomarcadores_mon,
  escalado = T,
  dir_to_save = directorio_multivariante_mon,
  original = F,
  grupo = grupo,
  obesidad = obesidad,
  common_axes = 18,
  text_size = 15,
  altura = 17,
  ancho = 16
  
)

plots_mon_NOT_scaled <- wrapper_plot_contribuciones(
  modelo,
  biomarcadores_mon,
  escalado = F,
  dir_to_save = directorio_multivariante_mon_Not_Scaled,
  original = F,
  grupo = grupo,
  obesidad = obesidad,
  common_axes = 18,
  text_size = 15,
  altura = 17,
  ancho = 16
)


plots_MOFA_scaled <- wrapper_plot_contribuciones(
  modelo,
  biomarcadores_MOFA,
  escalado = T,
  dir_to_save = directorio_multivariante_MOFA,
  original = F,
  grupo = grupo,
  obesidad = obesidad,
  common_axes = 18,
  text_size = 15,
  altura = 17,
  ancho = 16
)


plots_MOFA_NOT_scaled <- wrapper_plot_contribuciones(
  modelo,
  biomarcadores_MOFA,
  escalado = F,
  dir_to_save = directorio_multivariante_MOFA_Not_Scaled,
  original = F,
  grupo = grupo,
  obesidad = obesidad,
  common_axes = 18,
  text_size = 15,
  altura = 17,
  ancho = 16
)

plots_original_mon_scaled <- wrapper_plot_contribuciones(
  modelo,
  biomarcadores_mon,
  escalado = T,
  dir_to_save = directorio_original_mon_multivariante,
  original = T,
  grupo = grupo,
  obesidad = obesidad,
  common_axes = 18,
  text_size = 15,
  altura = 17,
  ancho = 16
)


plots_original_mon_not_Scaled <- wrapper_plot_contribuciones(
  modelo,
  biomarcadores_mon,
  escalado = F,
  dir_to_save = directorio_original_mon_multivariante_Not_Scaled,
  original = T,
  grupo = grupo,
  obesidad = obesidad,
  common_axes = 18,
  text_size = 15,
  altura = 17,
  ancho = 16
)


plots_original_MOFA <- wrapper_plot_contribuciones(
  modelo,
  biomarcadores_mon,
  escalado = T,
  dir_to_save = directorio_original_MOFA_multivariante,
  original = T,
  grupo = grupo,
  obesidad = obesidad,
  common_axes = 18,
  text_size = 15,
  altura = 17,
  ancho = 16
)


plots_original_MOFA_not_Scaled <- wrapper_plot_contribuciones(
  modelo,
  biomarcadores_mon,
  escalado = F,
  dir_to_save = directorio_original_MOFA_multivariante_Not_Scaled,
  original = T,
  grupo = grupo,
  obesidad = obesidad,
  common_axes = 18,
  text_size = 15,
  altura = 17,
  ancho = 16
)
### Correlacionamos los biomarcadores reconstruidos, con los scores, pero
### para tener mejor inciso, correlacionamos aquellos scores correspondinetes a
### al grupo de interes: cor(scores[mujeres,] ,senal[mujeres,])

##====contribuciones generales====
R <- reconstruccion(modelo)
plot_contributions_general(
  R = R,
  biomarcadores = biomarcadores_mon,
  escalado = T,
  componente = 1,
  text_size = 16,
  dir_to_save = directorio_multivariante_mon,
  altura = 18,
  ancho = 15
)
plot_contributions_general(
  R = R,
  biomarcadores = biomarcadores_mon,
  escalado = T,
  componente = 2,
  text_size = 16,
  dir_to_save = directorio_multivariante_mon,
  altura = 18,
  ancho = 15
)

plot_contributions_general(
  R = R,
  biomarcadores = biomarcadores_mon,
  escalado = F,
  componente = 1,
  text_size = 14,
  dir_to_save = directorio_multivariante_mon_Not_Scaled,
  altura = 15,
  ancho = 15
)
plot_contributions_general(
  R = R,
  biomarcadores = biomarcadores_mon,
  escalado = F,
  componente = 2,
  text_size = 14,
  dir_to_save = directorio_multivariante_mon_Not_Scaled,
  altura = 15,
  ancho = 15
)

plot_contributions_general(
  R = R,
  biomarcadores = biomarcadores_MOFA,
  escalado = T,
  componente = 1,
  text_size = 14,
  dir_to_save = directorio_multivariante_MOFA,
  altura = 15,
  ancho = 15
)
plot_contributions_general(
  R = R,
  biomarcadores = biomarcadores_MOFA,
  escalado = T,
  componente = 2,
  text_size = 14,
  dir_to_save = directorio_multivariante_MOFA,
  altura = 15,
  ancho = 15
)

plot_contributions_general(
  R = R,
  biomarcadores = biomarcadores_MOFA,
  escalado = F,
  componente = 1,
  text_size = 14,
  dir_to_save = directorio_multivariante_MOFA_Not_Scaled,
  altura = 15,
  ancho = 15
)
plot_contributions_general(
  R = R,
  biomarcadores = biomarcadores_MOFA,
  escalado = F,
  componente = 2,
  text_size = 14,
  dir_to_save = directorio_multivariante_MOFA_Not_Scaled,
  altura = 15,
  ancho = 15
)

S <- get_senal_original(modelo)
plot_contributions_general(
  R = S,
  biomarcadores = biomarcadores_mon,
  escalado = T,
  componente = 1,
  text_size = 14,
  dir_to_save = directorio_original_mon_multivariante,
  altura = 15,
  ancho = 15
)
plot_contributions_general(
  R = S,
  biomarcadores = biomarcadores_mon,
  escalado = T,
  componente = 2,
  text_size = 14,
  dir_to_save = directorio_original_mon_multivariante,
  altura = 15,
  ancho = 15
)

plot_contributions_general(
  R = S,
  biomarcadores = biomarcadores_mon,
  escalado = F,
  componente = 1,
  text_size = 14,
  dir_to_save = directorio_original_mon_multivariante_Not_Scaled,
  altura = 15,
  ancho = 15
)
plot_contributions_general(
  R = S,
  biomarcadores = biomarcadores_mon,
  escalado = F,
  componente = 2,
  text_size = 14,
  dir_to_save = directorio_original_mon_multivariante_Not_Scaled,
  altura = 15,
  ancho = 15
)

plot_contributions_general(
  R = S,
  biomarcadores = biomarcadores_MOFA,
  escalado = T,
  componente = 1,
  text_size = 14,
  dir_to_save = directorio_original_MOFA_multivariante,
  altura = 15,
  ancho = 15
)
plot_contributions_general(
  R = S,
  biomarcadores = biomarcadores_MOFA,
  escalado = T,
  componente = 2,
  text_size = 14,
  dir_to_save = directorio_original_MOFA_multivariante,
  altura = 15,
  ancho = 15
)

plot_contributions_general(
  R = S,
  biomarcadores = biomarcadores_MOFA,
  escalado = F,
  componente = 1,
  text_size = 14,
  dir_to_save = directorio_original_MOFA_multivariante_Not_Scaled,
  altura = 15,
  ancho = 15
)
plot_contributions_general(
  R = S,
  biomarcadores = biomarcadores_MOFA,
  escalado = F,
  componente = 2,
  text_size = 14,
  dir_to_save = directorio_original_MOFA_multivariante_Not_Scaled,
  altura = 15,
  ancho = 15
)
