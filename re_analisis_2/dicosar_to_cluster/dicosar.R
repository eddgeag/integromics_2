

library(dicosar)


wrapper_dico <- function(R,grupo,obesidad) {
  
  ## obesidad dico
  matrices_obesidad <- split(R,obesidad)
  print("dico obesidad....")
  p.dico_OB <- dicosar(matrices_obesidad$`No Obese`,matrices_obesidad$Obese)$p_global
  ## grupo dico
  matrices_grupo <- split(R,grupo)
  p.dico_FM <- dicosar(matrices_grupo$Female,matrices_grupo$Male)$p_global
  print("dico grupo 1")
  
  p.dico_FP <- dicosar(matrices_grupo$Female,matrices_grupo$PCOS)$p_global
  print("dico grupo 2")
  
  p.dico_PM <- dicosar(matrices_grupo$Female,matrices_grupo$Male)$p_global
  print("dico grupo 3")
  
  
  ## grupo no obesos dico
  grupo_noob_logi  <- obesidad=="No Obese"
  R_grupo_noob <- R[grupo_noob_logi,]
  grupo_noob <- grupo[obesidad=="No Obese"]
  matrices_grupo_noob <- split(R_grupo_noob,grupo_noob)
  p.dico_FMNOOB <- dicosar(matrices_grupo_noob$Female,matrices_grupo_noob$Male)$p_global
  print("dico noob grupo 1")
  
  p.dico_FPNOOB <- dicosar(matrices_grupo_noob$Female,matrices_grupo_noob$PCOS)$p_global
  print("dico noob grupo 2")
  
  p.dico_PMNOOB <- dicosar(matrices_grupo_noob$Female,matrices_grupo_noob$Male)$p_global
  print("dico noob grupo 3")
  
  ## grupo obesos dico
  grupo_ob_logi  <- obesidad=="Obese"
  R_grupo_ob <- R[grupo_ob_logi,]
  grupo_ob <- grupo[obesidad=="Obese"]
  matrices_grupo_ob <- split(R_grupo_ob,grupo_ob)
  p.dico_FMOB <- dicosar(matrices_grupo_ob$Female,matrices_grupo_ob$Male)$p_global
  print("dico ob grupo 1")
  p.dico_FPOB <- dicosar(matrices_grupo_ob$Female,matrices_grupo_ob$PCOS)$p_global
  print("dico ob grupo 2")
  
  p.dico_PMOB <- dicosar(matrices_grupo_ob$PCOS,matrices_grupo_ob$Male)$p_global
  print("dico ob grupo 3")
  
  
  ## por females
  R_ <- R[grupo=="Female",]
  obesidad_ <- obesidad[grupo=="Female"]
  R_split <- split(R_,obesidad_)
  p.dico.f <- dicosar(R_split$Obese,R_split$`No Obese`)$p_global
  print("dico females ")
  
  ## por males
  R_ <- R[grupo=="Male",]
  obesidad_ <- obesidad[grupo=="Male"]
  R_split <- split(R_,obesidad_)
  p.dico.m <- dicosar(R_split$Obese,R_split$`No Obese`)$p_global
  print("dico males")
  
  ## por pcos
  R_ <- R[grupo=="PCOS",]
  obesidad_ <- obesidad[grupo=="PCOS"]
  R_split <- split(R_,obesidad_)
  p.dico.p <- dicosar(R_split$Obese,R_split$`No Obese`)$p_global
  print("dico pcos")
  
  resvan.list <- (
    list(
      obesidad = p.dico_OB,
      sexo = p.dico_FM,
      females = p.dico_FP,
      hermas = p.dico_PM,
      sex.no = p.dico_FMNOOB,
      females.no = p.dico_FPNOOB,
      hermas.no = p.dico_FMNOOB,
      sex.ob  = p.dico_FMOB,
      females.ob = p.dico_FPOB,
      hermas.ob = p.dico.f,
      femalesnoob = p.dico.m,
      malesnoob = p.dico.m,
      pcosnoob = p.dico.p
    )
  )
  

  
  resvan1 <- as.data.frame(Reduce("rbind", resvan.list))
  
  
  return(resvan1)
}

ruta <- "./re_analisis_2/dicosar_data.rds"
dats <- readRDS(ruta)

grupo <- dats$grupo
obesidad <- dats$obesidad
pre <- dats$pre
R <- dats$pre
biomarcadores <- dats$biomarcadores
biomarcadores_vis <- dats$biomarcadores_vis

ret <- wrapper_dico(R[,biomarcadores_vis],grupo = grupo,obesidad = obesidad)


