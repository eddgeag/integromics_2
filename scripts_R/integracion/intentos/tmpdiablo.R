
#### Entrenamos al DIABLO
### obesidad
set.seed(123456789)
idx <- sample(1:nrow(metaboloma),size = nrow(metaboloma)*0.7)

train.met <- metaboloma[idx,]
train.mic <- metagenoma[idx,]
test.met <- metaboloma[-idx,]
test.mic <- metagenoma[-idx,]
grupo.train <- grupo[idx]
obesidad.train <- obesidad[idx]
grupo.test <- grupo[-idx]
obesidad.test <- obesidad[-idx]

X.list <- list(metaboloma=as.matrix(train.met),metagenoma=as.matrix(train.mic))
correlacion <- pls(X.list$metaboloma,X.list$metagenoma,ncomp = 1)
correlacion <- cor(correlacion$variates$X,correlacion$variates$Y)

(diseno <- t(matrix(c(0,correlacion,correlacion,0),byrow = T,ncol = 2)))


model.preperf <- block.splsda(X=list(metaboloma=train.met,
                                     metaganoma=train.mic),
                              Y=obesidad.train,ncomp = 10,design = diseno)
model.perf <- perf(model.preperf,dist = "all",validation = "Mfold",folds = 3,nrepeat = 15,signif.threshold = 0.05)

ncomps <- 3
tune.diablo <- tune.block.splsda(list(metaboloma=train.met,
                                      metaganoma=train.mic),
                                 Y = obesidad.train,design = diseno,
                                 ncomp = ncomps,
                                 test.keepX = list(metaboloma=which(colnames(metaboloma) %in% names(variables_significativas)),
                                                   metagenoma=which(colnames(metagenoma) %in% names(variables_significativas))),
                                 validation = "Mfold",folds = 3,nrepeat = 15,signif.threshold = 0.05,progressBar = T)


keep.list <- tune.diablo$choice.keepX
modelo.diablo.obesidad <- block.splsda(list(metaboloma=train.met,
                                            metaganoma=train.mic),
                                       Y=obesidad.train,
                                       design = diseno,ncomp = ncomps,keepX = keep.list)


predict.diablo.obesidad <-  stats::predict(modelo.diablo.obesidad,list(metaboloma=test.met,
                                                                       metaganoma=test.mic),type="class")
(confusion.mat.obesidad = caret::confusionMatrix(data = as.factor(predict.diablo.obesidad$class$max.dist$metaboloma[,2]),
                                                 reference=obesidad.test))$overall[1]

### sexo
set.seed(123456789)
idx <- sample(1:nrow(metaboloma),size = nrow(metaboloma)*0.7)

train.met <- metaboloma[idx,]
train.mic <- metagenoma[idx,]
test.met <- metaboloma[-idx,]
test.mic <- metagenoma[-idx,]
sexo.train <- sexo[idx]
sexo.test <- sexo[-idx]

X.list <- list(metaboloma=as.matrix(train.met),metagenoma=as.matrix(train.mic))
correlacion <- pls(X.list$metaboloma,X.list$metagenoma,ncomp = 1)
correlacion <- cor(correlacion$variates$X,correlacion$variates$Y)

(diseno <- t(matrix(c(0,correlacion,correlacion,0),byrow = T,ncol = 2)))


model.preperf <- block.splsda(X=list(metaboloma=train.met,
                                     metaganoma=train.mic),
                              Y=sexo.train,ncomp = 1,design = diseno)
model.perf <- perf(model.preperf,dist = "all",validation = "Mfold",folds = 3,nrepeat = 15,signif.threshold = 0.05)

ncomps <- 3

tune.diablo <- tune.block.splsda(list(metaboloma=train.met,
                                      metaganoma=train.mic),
                                 Y = sexo.train,design = diseno,
                                 ncomp = ncomps,
                                 test.keepX = list(metaboloma=which(colnames(metaboloma) %in% names(variables_significativas)),
                                                   metagenoma=which(colnames(metagenoma) %in% names(variables_significativas))),
                                 validation = "Mfold",folds = 3,nrepeat = 15,signif.threshold = 0.05,progressBar = T)


keep.list <- tune.diablo$choice.keepX
modelo.diablo.sexo <- block.splsda(list(metaboloma=train.met,
                                        metaganoma=train.mic),
                                   Y=sexo.train,
                                   design = diseno,ncomp = ncomps,keepX = keep.list)


predict.diablo.sexo <-  stats::predict(modelo.diablo.sexo,list(metaboloma=test.met,
                                                               metaganoma=test.mic),type="class")
(confusion.mat.sexo = caret::confusionMatrix(data = as.factor(predict.diablo.sexo$class$max.dist$metaboloma[,2]),
                                             reference=sexo.test))$ovarall[1]

### Tenemos que ver las contribuciones

vars.sexo <- extraer_variables(modelo.diablo.sexo)
vars.obesidad <- extraer_variables(modelo.diablo.obesidad)

variables.obesidad.only <- vars.obesidad[!vars.obesidad %in% vars.sexo]
variables.sexo.only <- vars.sexo[!vars.sexo %in% vars.obesidad]

comunes <- intersect(vars.sexo,vars.obesidad)

todas.vars <- unique(c(vars.sexo,vars.obesidad,comunes))

### Ahora realizamos el analisis MOFA2
directorio_mofa_results <- "./scripts_R/integracion/mofa2_Results"
if(!dir.exists(directorio_mofa_results)){
  dir.create(directorio_mofa_results)
}
p1 <- plot_variance_explained(modelo)
p2 <- plot_variance_explained(modelo, x="group", y="factor", plot_total = T)[[2]]
p3 <- ggarrange(p1,p2)

modelo@covariates <-samples_metadata(modelo)


correlate_factors_with_covariates(modelo,
                                  transpose = F,
                                  covariates = colnames(samples_metadata(modelo))[-c(5,6,8,14,16)],
)

factores <- get_factors(modelo,scale = T)[[1]]

grupo <- datos$comunes$grupo
obesidad <- datos$comunes$obesidad
analisis_factores <- analyze_data(factores,grupo,obesidad,correccion = 2)

analisis_factores.pvals <- as.matrix(analisis_factores$p_valores)[,-c(1,3)]
analisis_factores.medias <- as.matrix(analisis_factores$medianas)
colnames(analisis_factores.pvals) <- colnames(analisis_factores.medias)

corrplot::corrplot(analisis_factores.medias,
                   is.corr = F,
                   p.mat = analisis_factores.pvals,
                   sig.level = 0.05,insig = "label_sig")

csa <- plot_weights(modelo,view = 2,factors = 3)

variables_significativas




