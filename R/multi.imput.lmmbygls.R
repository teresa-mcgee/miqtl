#' @export
multi.imput.lmmbygls <- function(formula, data=NULL, pheno.id="SUBJECT.NAME",
                                 y=NULL, fit0=NULL, null.formula=NULL,
                                 num.imp, founders=founders, X.probs,
                                 K=NULL,
                                 use.par, fix.par=NULL, model=c("additive", "full"), p.value.method=c("LRT", "ANOVA"), locus.as.fixed=TRUE,
                                 use.lmer, impute.map,
                                 brute=TRUE, seed=1, do.augment,
                                 weights=NULL,
                                 MX0=NULL, My=NULL,
                                 GxT.treatment=NULL, treatment=NULL,
                                 return.allele.effects=FALSE, return.qtl.predictor=FALSE){
  model <- model[1]
  p.value.method <- p.value.method[1]
  eigen.K <- logDetV <- M <- allele.effects <- qtl.predictor <- NULL
  imp.LOD.GxT <- imp.p.value.GxT <- NULL
  imp.allele.effects.GxT <- imp.allele.effects.GxT.delta <- NULL
  GxT.T.design <- GxT.treatment.levels <- NULL
  if(!is.null(GxT.treatment)){
    imp.LOD.GxT <- imp.p.value.GxT <- rep(0, num.imp)
    ## Treatment doesn't vary by imputation, so its dummy design is built
    ## once here rather than rebuilt inside the per-imputation loop below.
    GxT.T.design <- make.treatment.design(treatment)
    GxT.treatment.levels <- levels(as.factor(treatment))
  }
  ## For allele effect plots
  if(return.allele.effects){
    allele.effects <- matrix(NA, nrow=length(founders), ncol=num.imp,
                             dimnames=list(founders, paste0("imp", 1:num.imp)))
    if(!is.null(GxT.treatment)){
      imp.allele.effects.GxT <- array(NA, dim=c(length(founders), length(GxT.treatment.levels), num.imp),
                                      dimnames=list(founders, GxT.treatment.levels, paste0("imp", 1:num.imp)))
      imp.allele.effects.GxT.delta <- array(NA, dim=c(length(founders), length(GxT.treatment.levels) - 1, num.imp),
                                            dimnames=list(founders, GxT.treatment.levels[-1], paste0("imp", 1:num.imp)))
    }
  }
  ## For conditional scans through residuals
  if(return.qtl.predictor){
    if(!is.null(data)){ n <- nrow(data) }
    else{ n <- length(y) }
    qtl.predictor <- matrix(NA, nrow=n, ncol=num.imp)
  }
  if(is.null(fit0)){
    null.formula <- make.null.formula(formula=formula, do.augment=do.augment)
    if(use.lmer){
      fit0 <- lmmbylmer(null.formula, data=data, REML=FALSE, weights=weights)
    }
    else{
      fit0 <- lmmbygls(null.formula, data=data, pheno.id=pheno.id, K=K, use.par=use.par, brute=brute, weights=weights)
      K <- fit0$K
    }
  }
  if(is.null(weights) & !use.lmer){
    eigen.K <- fit0$eigen.K
  }
  if(!is.null(fix.par) & !use.lmer){
    M <- fit0$M
    logDetV <- fit0$logDetV
  }
  full.to.dosages <- straineff.mapping.matrix()
  
  imp.logLik <- imp.h2 <- imp.df <- imp.LOD <- imp.p.value <- rep(0, num.imp)

  set.seed(seed)
  for(i in 1:num.imp){
    if(model == "additive"){
      if(any(X.probs < 0)){
        X.probs[X.probs < 0] <- 0
        X.probs <- t(apply(X.probs, 1, function(x) x/sum(x)))
      }
      X <- run.imputation(diplotype.probs=X.probs, impute.map=impute.map) %*% full.to.dosages
      colnames(X) <- gsub(pattern="/", replacement=".", x=founders, fixed=TRUE)
    }
    if(model == "full"){
      X <- run.imputation(diplotype.probs=X.probs, impute.map=impute.map)
      colnames(X) <- gsub(pattern="/", replacement=".", x=colnames(X.probs), fixed=TRUE)
    }
    keep.col <- 1:ncol(X)
    if(locus.as.fixed){
      max.column <- which.max(colSums(X, na.rm=TRUE))[1]
      keep.col <- keep.col[keep.col != max.column]
      X <- X[,keep.col]
    }

    locus.formula <- make.alt.formula(formula=formula, X=X, do.augment=do.augment)
    if(use.lmer){
      data <- cbind(null.data, X)
      fit1 <- lmmbylmer(formula=locus.formula, data=data, REML=FALSE, weights=weights)
      imp.logLik[i] <- as.numeric(logLik(fit1))
      imp.h2[i] <- NA
      imp.df[i] <- length(fixef(fit1))
      imp.LOD[i] <- log10(exp(imp.logLik[i] - as.numeric(logLik(fit0))))
      imp.p.value[i] <- pchisq(q=-2*(as.numeric(logLik(fit0)) - imp.logLik), df=imp.df - length(fixef(fit0)), lower.tail=FALSE)
      fit1$locus.effect.type <- "fixed"
    }
    else{
      if(locus.as.fixed){
        X.locus <- X # kept for GxT interaction columns, before X is cbound with the null design below
        X <- cbind(fit0$x, X)
        fit1 <- lmmbygls(formula=locus.formula, pheno.id=pheno.id, eigen.K=eigen.K, K=K,
                         y=y, X=X,
                         logDetV=logDetV, M=M,
                         use.par="h2", fix.par=fix.par,
                         MX0=MX0, My=My,
                         brute=brute, weights=weights)
        imp.logLik[i] <- fit1$logLik
        imp.h2[i] <- fit1$h2
        imp.df[i] <- fit1$rank
        imp.LOD[i] <- log10(exp(imp.logLik[i] - fit0$logLik))
        imp.p.value[i] <- get.p.value(fit0=fit0, fit1=fit1, method=p.value.method)
        fit1$locus.effect.type <- "fixed"

        if(return.allele.effects){
          allele.effects[,i] <- get.allele.effects.from.fixef(fit=fit1, founders=founders, allele.in.intercept=founders[max.column])
        }
        if(return.qtl.predictor){
          qtl.predictor[,i] <- regress.out.qtl(fit1=fit1, null.formula=null.formula, alt.formula=locus.formula,
                                               locus.as.fixed=locus.as.fixed)
        }
        ## GxT interaction test for this imputation: null is fit1.treatment
        ## (locus main effect + treatment main effect), not fit1 (no
        ## treatment term) or fit0 -- same marginality rationale as
        ## scan.h2lmm()'s ROP path. Reuses fit0's M/logDetV/MX0/My
        ## throughout, same as fit1.
        if(!is.null(GxT.treatment)){
          X.int <- make.interaction.design(X=X.locus, treatment=treatment)
          fit1.treatment <- lmmbygls(pheno.id=pheno.id, eigen.K=eigen.K, K=K,
                                     y=y, X=cbind(fit1$x, GxT.T.design),
                                     logDetV=logDetV, M=M,
                                     use.par="h2", fix.par=fix.par,
                                     MX0=MX0, My=My,
                                     brute=brute, weights=weights)
          fit1.GxT <- lmmbygls(formula=make.GxT.formula(formula=formula, X=X.locus, GxT.treatment=GxT.treatment, do.augment=do.augment),
                               pheno.id=pheno.id, eigen.K=eigen.K, K=K,
                               y=y, X=cbind(fit1.treatment$x, X.int),
                               logDetV=logDetV, M=M,
                               use.par="h2", fix.par=fix.par,
                               MX0=MX0, My=My,
                               brute=brute, weights=weights)
          imp.LOD.GxT[i] <- log10(exp(fit1.GxT$logLik - fit1.treatment$logLik))
          imp.p.value.GxT[i] <- get.p.value(fit0=fit1.treatment, fit1=fit1.GxT, method=p.value.method)

          if(return.allele.effects){
            imp.allele.effects.GxT[,,i] <- get.allele.effects.from.fixef.GxT(fit.GxT=fit1.GxT, founders=founders,
                                                                              allele.in.intercept=founders[max.column],
                                                                              treatment=treatment)
            imp.allele.effects.GxT.delta[,,i] <- get.interaction.deltas.from.fixef(fit.GxT=fit1.GxT, founders=founders,
                                                                                    allele.in.intercept=founders[max.column],
                                                                                    treatment=treatment)
          }
        }
      }
      else{
        null.formula <- make.null.formula(formula=formula, do.augment=do.augment)
        fit1 <- lmmbygls.random(formula=null.formula, pheno.id=pheno.id,  
                                y=y, X=fit0$x, K=K, eigen.K=eigen.K, Z=X,
                                use.par="h2", null.h2=fix.par,
                                brute=brute, weights=weights)
        imp.logLik[i] <- fit1$REML.logLik
        imp.h2[i] <- fit1$h2
        imp.df[i] <- fit1$rank
        imp.LOD[i] <- log10(exp(imp.logLik[i] - fit0$REML.logLik))
        imp.p.value[i] <- get.p.value(fit0=fit0, fit1=fit1, method=p.value.method)
        fit1$locus.effect.type <- "random"
        
        if(return.allele.effects){
          allele.effects[,i] <- get.allele.effects.from.ranef(fit=fit1, founders=founders)
        }
      }
    }
  }
  return(list(h2=imp.h2,
              LOD=imp.LOD,
              p.value=imp.p.value,
              LOD.GxT=imp.LOD.GxT,
              p.value.GxT=imp.p.value.GxT,
              allele.effects=allele.effects,
              allele.effects.GxT=imp.allele.effects.GxT,
              allele.effects.GxT.delta=imp.allele.effects.GxT.delta,
              qtl.predictor=qtl.predictor,
              locus.effect.type=fit1$locus.effect.type))
}
