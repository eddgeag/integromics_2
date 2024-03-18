




datos <- readRDS("../datos/preprocesado_08_09_23/novoom-04-10-23.rds")


ip <- datos$comunes$ip

for(j in 1:ncol(ip)){
  
  hist(ip[,j],main=colnames(ip)[j])
  
}

ip$ADIPOQ <- log(ip$ADIPOQ)

ip$LEP <- log(ip$LEP)

ip$VASPIN <- log(ip$VASPIN)

ip$hsCRP <-  log(ip$hsCRP)

hist(ip$ADIPOQ)
hist(ip$LEP)
hist(ip$VASPIN)
hist(ip$hsCRP)
