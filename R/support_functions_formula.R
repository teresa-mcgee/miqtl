## Formula manipulation functions
make.null.formula <- function(formula, do.augment){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
  this.formula <- as.formula(ifelse(do.augment, paste0(this.formula.string, " + augment.indicator"), this.formula.string))
  return(this.formula)
}
make.alt.formula <- function(formula, X, do.augment){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
  this.formula.string <- ifelse(do.augment, paste0(this.formula.string, " + augment.indicator"), this.formula.string)
  this.formula <- as.formula(paste(this.formula.string, paste(gsub(pattern="/", replacement=".", x=colnames(X), fixed=TRUE), collapse=" + "), sep=" + "))
  return(this.formula)
}
make.snp.null.formula <- function(formula, condition.loci, X.list, model){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
  if(!is.null(condition.loci)){
    for(i in 1:length(condition.loci)){
      if(model == "additive"){
        this.formula.string <- paste(this.formula.string, paste("cond_SNP", i, sep="_"), sep=" + ")
      }
      else if(model == "full"){
        this.formula.string <- paste(this.formula.string, c(paste("cond_SNP", i, "aa", sep="_"), paste("cond_SNP", i, "Aa", sep="_")), sep=" + ")
      }
    }
  }
  return(as.formula(this.formula.string))
}
make.snp.alt.formula <- function(formula, model){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
  if(model == "additive"){
    this.formula <- as.formula(paste(this.formula.string, "SNP", sep=" + "))
  }
  else if(model == "full"){
    this.formula <- as.formula(paste(this.formula.string, "SNP_aa", "SNP_Aa", sep=" + "))
  }
  return(this.formula)
}
remove.whitespace.formula <- function(formula){
  formula.string <- paste0(Reduce(paste, deparse(formula)))
  formula.string <- gsub("[[:space:]]", "", formula.string)
  return(as.formula(formula.string))
}
check.for.lmer.formula <- function(formula){
  formula.string <- paste0(Reduce(paste, deparse(formula)))
  formula.string <- gsub("[[:space:]]", "", formula.string)
  use.lmer <- grepl(pattern="\\([a-zA-Z0-9\\.]+\\|[a-zA-Z0-9\\.]+\\)", x=formula.string, perl=TRUE)
  return(use.lmer)
}
process.random.formula <- function(geno.id){
  random.formula <- as.formula(paste("~", paste(geno.id, "1", sep=" - ")))
  return(random.formula)
}
make.GxT.formula <- function(formula, X, GxT.treatment, do.augment){
  ## Cosmetic only, for output$formula.GxT labeling: whenever a per-locus
  ## fit is built from y/X directly (as every GxT call site is), lmmbygls()'s
  ## formula argument is never actually used to construct the design matrix
  ## (see the `if(!is.null(data) & is.null(y))` branch in lmmbygls()), so
  ## this does not need to exactly reproduce make.interaction.design()'s
  ## column construction -- only to describe it.
  this.formula <- make.alt.formula(formula=formula, X=X, do.augment=do.augment)
  this.formula.string <- Reduce(paste, deparse(this.formula))
  locus.names <- gsub(pattern="/", replacement=".", x=colnames(X), fixed=TRUE)
  interaction.terms <- paste(locus.names, GxT.treatment, sep=":")
  as.formula(paste(this.formula.string, paste(interaction.terms, collapse=" + "), sep=" + "))
}
## Builds the locus x treatment interaction design columns for a GxT test.
## X is the already reference-column-dropped locus dosage matrix used for
## the main-effect fit (as-is; no additional collinearity handling needed
## on the genotype side). `treatment` is dummy-coded with standard
## treatment contrasts (reference level dropped) for factors/character
## vectors with >2 levels; a 2-level factor collapses to a single dummy
## column; a numeric treatment (binary indicator or continuous) is used
## as-is, without model.matrix()'s intercept handling, since there is no
## reference level to drop.
make.interaction.design <- function(X, treatment){
  ## NA treatment values should already be impossible here -- make.processed.data()
  ## folds GxT.treatment into its covariates so NA rows are dropped upstream via
  ## model.frame(). Checked explicitly (rather than left to silently misalign)
  ## because model.matrix() below would otherwise drop NA rows on its own and
  ## desync row counts against X without any other warning.
  if(anyNA(treatment)){
    stop("make.interaction.design(): treatment contains NA values; these should have been removed upstream by make.processed.data()", call.=FALSE)
  }
  if(is.numeric(treatment)){
    T.design <- matrix(treatment, ncol=1, dimnames=list(NULL, "trt"))
  }
  else{
    treatment <- as.factor(treatment)
    T.design <- model.matrix(~treatment)[, -1, drop=FALSE]
    colnames(T.design) <- sub(pattern="^treatment", replacement="", x=colnames(T.design))
  }
  interaction.design <- do.call(cbind, lapply(seq_len(ncol(T.design)), function(j) X * T.design[, j]))
  colnames(interaction.design) <- as.vector(outer(colnames(X), colnames(T.design), paste, sep="."))
  rownames(interaction.design) <- rownames(X)
  return(interaction.design)
}
