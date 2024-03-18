vanValen.test2 <- function(X, Y) {
  
  X.standard <- scale(X)
  Y.standard <- scale(Y)
  X.dat <- list(X.standard, Y.standard)
  X.dat <- lapply(X.dat, as.matrix)
  medianas <- lapply(X.dat, function(x)
    apply(x, 2, median))
  abs.dev.factor <- vector("list", length = length(X.dat))
  for (j in 1:length(X.dat)) {
    y <- X.dat[[j]]
    mediana_x <- medianas[[j]]
    abs.dev.factor[[j]] <- abs(sweep(y, 2, mediana_x))
    
  }
  
  norma.euclidea <- function(x)    sqrt(sum(x ^ 2))
  distancias1 <- apply(abs.dev.factor[[1]], 1, norma.euclidea)
  distancias2 <- apply(abs.dev.factor[[2]], 1, norma.euclidea)
  aux.factor <- as.factor(c(rep(1, length(distancias1)),
                            rep(2, length(distancias2))))
  (distancias <- reshape2::melt(t(rbind(
    distancias1, distancias2
  ))))
  norm.test <-
    ks.test(residuals(lm(distancias1 ~ distancias2)), "pnorm")$p.value
  vari.test <-
    leveneTest(distancias$value, distancias$Var2)$`Pr(>F)`[1]
  
  if (norm.test > 0.05 && vari.test > 0.05) {
    p.valor <-
      t.test(distancias$value ~ distancias$Var2, var.equal = T)$p.value
    stat <-
      t.test(distancias$value ~ distancias$Var2, var.equal = F)$statistic
    
    
  } else{
    p.valor <-
      t.test(distancias$value ~ distancias$Var2, var.equal = F)$p.value
    stat <-
      t.test(distancias$value ~ distancias$Var2, var.equal = F)$statistic
    
  } 
  res <- c(p=p.valor, t=stat)
  return(res)
  
}

