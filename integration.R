library(limma)
library(ALDEx2)
library(MOFA2)
library(ggrepel)
source(
  "./scripts_R/scripts_utiles/scripts_funciones/analisis_univariante_e_interpretacion.R"
)
source("./scripts_R/scripts_utiles/scripts_funciones/manova_vanvalen.R")


mofa_componentes_original <-
  function(ncomp,
           semilla,
           directorio_modelo,
           mofa.obj,
           sample_metadata_) {
    samples_metadata(mofa.obj) <- sample_metadata_
    
    if (!dir.exists(directorio_modelo)) {
      dir.create(directorio_modelo)
    }
    data_opts <- get_default_data_options(mofa.obj)
    data_opts$scale_views <- F
    
    
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
    
    MOFAobject.trained <-
      run_mofa(MOFAobject, outfile, use_basilisk = T)
    return(MOFAobject.trained)
    
  }

datos__ <- readRDS("../datos/preprocesado_08_09_23/novoom-04-10-23.rds")
datos <- readRDS("../datos/preprocesado_08_09_23-26/novoom.rds")


grupo <- datos$comunes$grupo
obesidad <- datos$comunes$obesidad
sexo <- datos$comunes$variables_in_bacteria$SEX
## metaboloma
x <- datos$comunes$metaboloma[-c(38:40), ]
succinato <- colMeans(datos$comunes$metaboloma)
x <- rbind(x, succinato)
rownames(x)[38] <- "Succ"
rownames(x) <- gsub("_mean", replacement = "", x = rownames(x))
x.tmp <- x[-38]
x.tmp <- add_classificacion_metaboloma(x.tmp[-38, ])
rownames(x)[-38] <- paste(x.tmp$grupo_meta, rownames(x.tmp), sep = "_")

x.t <- scale(t(x))

## microbioma

mat <-  scale(t(mean_aldex(datos$comunes$microbiota$genero)))
fun_pca <- function(x){
  
  
  pcx <- prcomp(x)
  scores <- pcx$x
  dfplot <- data.frame(scores,grupo=grupo,rownames=rownames(x),obesidad=obesidad)
  varianzas <- pcx$sdev^2/sum(pcx$sdev^2)*100
  p <- ggplot(dfplot,aes(PC1,PC2,color=obesidad,label=rownames))+geom_point()+xlab(paste("PC1",round(varianzas[1],2),
                                                                                         "%"))+ylab(paste("PC2",round(varianzas[2],2),
                                                                                                          "%"))+geom_text_repel()
  return(p)
}

## ip

datos$comunes$ip <- as.data.frame(lapply(datos$comunes$ip,as.numeric),row.names = rownames(datos$comunes$ip))
ip <- t(scale(datos$comunes$ip))
# fun_pca(mat,sexo)


# fun_pca(ip,grupo)
# ip.d <- t(ip[,-27])

## covariables

colnames(x.t)[grep("CARBOHIDRATOS_GRASAS_KETONA_GLYCEROL_",colnames(x.t))]<-gsub("CARBOHIDRATOS_GRASAS_KETONA_GLYCEROL_","CGKG_",colnames(x.t)[grep("CARBOHIDRATOS_GRASAS_KETONA_GLYCEROL_",colnames(x.t))])




directorio.general <- "./modelos/14_10_24_With_OUTLIERS-17-04-26"

if (!dir.exists(directorio.general)) {
  dir.create(directorio.general,recursive = T)
}


ISI <- datos$comunes$variables_in_bacteria$SOG_ISI
testosterona <- datos$comunes$variables_in_bacteria[,c("FREE_TEST","TOTAL_TEST")]
clinicas.mean <- datos$comunes$clinicos
covariables <- cbind(ISI,clinicas.mean,testosterona)
covariables <- covariables[,-grep("EDAD",colnames(covariables))]
covariables <- as.data.frame(apply(covariables,2,scale))
rownames(covariables) <- colnames(datos$comunes$metaboloma)
covariables$sample <- rownames(covariables)


fun_pca(covariables[,-ncol(covariables)])

### muestras outliers en clinicas

# soi <- c("PdGLP6", "VOGLP4" ,"MdGLP4" ,"MdGLP7")



# idx <- sapply(soi,function(covariablesx) grep(x,rownames(covariables)))

# covariables
covariables_soi <- as.data.frame((covariables))
## metaboloma
metaboloma_soi <- t(x.t)
## proteoma
proteoma_soi <- ip
## metagenoma
metagenoma_soi <- t(mat)

# covariables
# covariables_soi <- as.data.frame((covariables))
# ## metaboloma
# metaboloma_soi <- t(x.t)
# ## proteoma
# proteoma_soi <- ip
# ## metagenoma
# metagenoma_soi <- t(mat)

datos.list <- list(metaboloma=metaboloma_soi,
                   metagenoma=metagenoma_soi,
                   ip = proteoma_soi)

metabolites <- c(
  "2-oxoisocaproic acid",
  "Leucine",
  "Isoleucine",
  "Valine",
  "Isobutyric acid",
  "2-oxoisovaleric acid",
  "3-hydroxybutyric acid",
  "Lactate",
  "Alanine",
  "Acetate",
  "N-acetyl glycoproteins",
  "Acetone",
  "Glutamic acid",
  "Pyruvate",
  "Pyroglutamic acid",
  "Glutamine",
  "Citrate",
  "Asparagine",
  "Creatinine",
  "Creatine",
  "Lysine",
  "Ornithine",
  "Choline",
  "Carnitine",
  "Betaine",
  "Glycine",
  "Threonine",
  "Glycerol",
  "Serine",
  "Proline",
  "Beta-Glucose",
  "D-Glucose",
  "1-methylhistidine",
  "Tyrosine",
  "Phenylalanine",
  "Tryptophan",
  "Formate",
  "Succinate"
)

# Mostrar el vector
print(metabolites)
lista.datos <- create_mofa_from_matrix(
  datos.list
)
rownames(datos.list$metaboloma) <-metabolites


rownames(datos.list$ip)

biomarkers <- c(
  "LBP",
  "sCD14",
  "GLP-2",
  "ZFP",
  "Adiponectin",
  "Leptin",
  "sLeptinR",
  "Ghrelin",
  "Omentin-1",
  "Vaspin",
  "Lipocalin-2",
  "Adipsin",
  "PAI-1",
  "Chemerin",
  "FGF-21",
  "FGF-23",
  "Galectin-3",
  "Pentraxin-3",
  "IL-18",
  "IL-6",
  "TNF-a",
  "GlycA",
  "GlycB",
  "GlycF",
  "GlycA H/W",
  "GlycB H/W",
  "hsCRP"
)


rownames(datos.list$ip) <- biomarkers

clinical <- c(
  "sample",
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
  "HOMA-IR",
  "ISI",
  "Triglycerides",
  "Cholesterol",
  "HDL-Cholesterol",
  "LDL-Cholesterol"
)



samples_covariables <- covariables_soi[,c("sample","BMI","WC","WHR","TOTAL_TEST","FREE_TEST","Total_ESTR","Free_ESTRA",
                                          "SHBG","GlucosaBasal_mean","InsulinaBasal_mean","HOMAIRmean","ISI",
                                          "TGmean","COLmean","HDLmean","LDLmean")]

colnames(samples_covariables) <- clinical
# Mostrar el vector
lista.datos <- create_mofa_from_matrix(
  datos.list
)


for (ncomp in 4:12) {
  for (semilla in sample(40000, size = 10)) {
    obj <-  mofa_componentes_original(ncomp,
                                      semilla,
                                      directorio.general,
                                      lista.datos,
                                      sample_metadata_ =samples_covariables)
  }
}

# # # #
