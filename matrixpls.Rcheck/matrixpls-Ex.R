pkgname <- "matrixpls"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('matrixpls')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("GSCA")
### * GSCA

flush(stderr()); flush(stdout())

### Name: GSCA
### Title: Generalized structured component analysis (GSCA) weights
### Aliases: GSCA

### ** Examples

if(!require(ASGSCA)){
    print("This example requires the ASGSCA package from Bioconductor")
} else{
# Run the example from ASGSCA package using GSCA estimation

data(GenPhen)
W0 <- matrix(c(rep(1,2),rep(0,8),rep(1,2),rep(0,8),rep(1,3),rep(0,7),rep(1,2)),
             nrow=8,ncol=4)
B0 <- matrix(c(rep(0,8),rep(1,2),rep(0,3),1,rep(0,2)),nrow=4,ncol=4)

# Set seed becayse ASGSCA uses random numbers as starting values 
set.seed(1)

GSCA.res <-GSCA(GenPhen,W0, B0,estim=TRUE,path.test=FALSE, 
                 latent.names=c("Gene1","Gene2",
                                "Clinical pathway 1",
                                "Clinical pathway 2"))

# Setup matrixpls to estimate the same model. Note that ASGSCA places dependent
# variables on columns but matrixpls uses rows for dependent variables

inner <- t(B0)
formative <- t(W0)
reflective <- matrix(0,8,4)

colnames(formative) <- rownames(reflective) <- names(GenPhen)

colnames(inner) <- rownames(inner) <- 
  rownames(formative) <- colnames(reflective) <-
  c("Gene1","Gene2","Clinical pathway 1","Clinical pathway 2")

model <- list(inner = inner, 
              reflective = reflective,
              formative = formative)

# Estimate using alternating least squares

matrixpls.res1 <- matrixpls(cov(GenPhen),  model,
                            outerEstim = outerEstim.gsca,
                            innerEstim = innerEstim.gsca)

# Estimate using direct minimization of the estimation criterion
# Set the convergence criterion to be slightly stricter than normally
# to get indentical results

matrixpls.res2 <- matrixpls(cov(GenPhen),  model,
                            weightFun = weightFun.optim,
                            optimCrit = optimCrit.gsca,
                            control = list(reltol = 1e-12))

# Compare the weights

do.call(cbind,lapply(list(ASGSCA =GSCA.res[["Weight"]],
                          matrixpls_als = t(attr(matrixpls.res1,"W")),
                          matrixpls_optim =t(attr(matrixpls.res2,"W"))),
                     function(W) W[W!=0]))


# Check the criterion function values

optimCrit.gsca(matrixpls.res1)
optimCrit.gsca(matrixpls.res2)

}



cleanEx()
nameEx("cei")
### * cei

flush(stderr()); flush(stdout())

### Name: cei
### Title: Composite Equivalence Indices
### Aliases: cei

### ** Examples

# Load the Tenenhaus et al 2005 model and data from semPLS
library(semPLS)
data(ECSImobi)
data(mobi)

# Reflective and empty formative model
reflective<- ECSImobi$M
formative <- t(reflective)
formative[] <- 0

# Estimation using covariance matrix
model <- list(inner =  t(ECSImobi$D),
              reflective = reflective,
              formative = formative)


S <- cor(mobi)

matrixpls.ModeA <- matrixpls(S, model, innerEstim = innerEstim.centroid)
matrixpls.Fixed <- matrixpls(S, model, weightFun = weightFun.fixed)

cei(matrixpls.ModeA)
cei(matrixpls.ModeA, matrixpls.Fixed)



cleanEx()
nameEx("estimator")
### * estimator

flush(stderr()); flush(stdout())

### Name: estimator
### Title: Parameter estimation of a model matrix
### Aliases: estimator estimator.ols estimator.tsls estimator.plscLoadings
###   estimator.efaLoadings estimator.cfaLoadings estimator.plsreg

### ** Examples

# Run the education example from the book

# Sanchez, G. (2013) PLS Path Modeling with R
# Trowchez Editions. Berkeley, 2013. 
# http://www.gastonsanchez.com/PLS Path Modeling with R.pdf

education <- read.csv("http://www.gastonsanchez.com/education.csv")

Support <- c(0, 0, 0, 0, 0, 0)
Advising <- c(0, 0, 0, 0, 0, 0)
Tutoring <- c(0, 0, 0, 0, 0, 0)
Value <- c(1, 1, 1, 0, 0, 0)
# Omit two paths (compared to the model in the book) to achieve 
# identification of the 2SLS analysis
Satisfaction <- c(0, 0, 1, 1, 0, 0)
Loyalty <- c(0, 0, 0, 0, 1, 0)

inner <- rbind(Support, Advising, Tutoring, Value, Satisfaction, Loyalty)


reflective <- diag(6)[c(rep(1,4),
                        rep(2,4),
                        rep(3,4),
                        rep(4,4),
                        rep(5,3),
                        rep(6,4)),]
formative <- matrix(0, 6, 23)

colnames(inner) <- colnames(reflective) <- rownames(formative) <- rownames(inner)
rownames(reflective) <- colnames(formative) <- colnames(education)[2:24]

education.model <- list(inner = inner,
              reflective = reflective,
              formative = formative)

# Reverse code two variables
education[,c("sup.under","loy.asha")] <- - education[,c("sup.under","loy.asha")]

S <- cor(education[,2:24])

# PLSc with OLS regression

education.out <- matrixpls(S,education.model,
                      disattenuate = TRUE,
                      parametersReflective = estimator.plscLoadings)

# PLSc with 2SLS regresssion

education.out2 <- matrixpls(S,education.model,
                      disattenuate = TRUE,
                      parametersReflective = estimator.plscLoadings,
                      parametersInner = estimator.tsls)


# Disattenuated regression with unit-weighted scales and exploratory factor analysis
# reliability estimates (with unconstrained MINRES estimator)

education.out3 <- matrixpls(S,education.model,
                       disattenuate = TRUE,
                       weightFun = weightFun.fixed,
                       parametersReflective = estimator.efaLoadings)

# Disattenuated GSCA with 2SLS regression after disattenuated based on 
# confirmatory factor analysis reliability estimates


education.out4 <- matrixpls(S,education.model,
                       disattenuate = TRUE,
                       innerEstim = innerEstim.gsca,
                       outerEstim = outerEstim.gsca,
                       parametersInner = estimator.tsls,
                       parametersReflective = estimator.cfaLoadings)


# Compare the results

cbind(PLSc = education.out, PLSc_2sls = education.out2, 
      DR = education.out3, GSCAc = education.out4)

# Compare the reliability estimates

cbind(PLSc = attr(education.out,"Q"), PLSc_2sls = attr(education.out2,"Q"), 
      DR = attr(education.out3,"Q"), GSCAc = attr(education.out4,"Q"))



cleanEx()
nameEx("matrixpls.plspm")
### * matrixpls.plspm

flush(stderr()); flush(stdout())

### Name: matrixpls.plspm
### Title: A plspm compatibility wrapper for matrixpls
### Aliases: matrixpls.plspm

### ** Examples

cores <- getOption("mc.cores")
options(mc.cores=2)

# Run the example from plspm package

# load dataset satisfaction
data(satisfaction)
# inner model matrix
IMAG = c(0,0,0,0,0,0)
EXPE = c(1,0,0,0,0,0)
QUAL = c(0,1,0,0,0,0)
VAL = c(0,1,1,0,0,0)
SAT = c(1,1,1,1,0,0)
LOY = c(1,0,0,0,1,0)
sat_inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
# outer model list
sat_outer = list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)
# vector of modes (reflective indicators)
sat_mod = rep("A", 6)

# apply matrixpls
matrixpls.res <- matrixpls.plspm(satisfaction, sat_inner, sat_outer, sat_mod,
                                 scaled=FALSE, boot.val=FALSE)

print(summary(matrixpls.res))

options(mc.cores=cores)




cleanEx()
nameEx("matrixpls.sempls")
### * matrixpls.sempls

flush(stderr()); flush(stdout())

### Name: matrixpls.sempls
### Title: A semPLS compatibility wrapper for matrixpls
### Aliases: matrixpls.sempls

### ** Examples

if(!require(semPLS)){
    print("This example requires the semPLS package")
} else{
data(ECSImobi)

ecsi.sempls <- sempls(model=ECSImobi, data=mobi, wscheme="pathWeighting")
ecsi <- matrixpls.sempls(model=ECSImobi, data=mobi, wscheme="pathWeighting")

# If RUnit is installed check that the results are identical

if(require(RUnit)){
  checkEquals(ecsi.sempls,ecsi, check.attributes = FALSE)
}

ecsi

## create plots
densityplot(ecsi)
densityplot(ecsi, use="prediction")
densityplot(ecsi, use="residuals")

## Values of 'sempls' objects
names(ecsi)
ecsi$outer_weights
ecsi$outer_loadings
ecsi$path_coefficients
ecsi$total_effects


### using convenience methods to sempls results
## path coefficients
pathCoeff(ecsi)

## total effects
totalEffects(ecsi)

## get loadings and check for discriminant validity
(l <- plsLoadings(ecsi))
# outer loadings
print(l, type="outer", digits=2)
# outer loadings greater than 0.5
print(l,type="outer", cutoff=0.5, digits=2)
# cross loadings greater than 0.5
print(l, type="cross", cutoff=0.5, digits=2)


### R-squared
rSquared(ecsi)

}



cleanEx()
nameEx("matrixpls.sim")
### * matrixpls.sim

flush(stderr()); flush(stdout())

### Name: matrixpls.sim
### Title: Monte Carlo simulations with matrixpls
### Aliases: matrixpls.sim

### ** Examples

if(!require(simsem)){
    print("This example requires the simsem package")
} else{

#
# Runs the second model from
#
# Aguirre-Urreta, M., & Marakas, G. (2013). Partial Least Squares and Models with Formatively 
# Specified Endogenous Constructs: A Cautionary Note. Information Systems Research. 
# doi:10.1287/isre.2013.0493

library(MASS)

X <- diag(15) 
X[upper.tri(X, diag=TRUE)] <- c(
  1.000,
  0.640, 1.000,
  0.640, 0.640, 1.000,
  0.640, 0.640, 0.640, 1.000,
  0.096, 0.096, 0.096, 0.096, 1.000,
  0.096, 0.096, 0.096, 0.096, 0.640, 1.000,
  0.096, 0.096, 0.096, 0.096, 0.640, 0.640, 1.000,
  0.096, 0.096, 0.096, 0.096, 0.640, 0.640, 0.640, 1.000,
  0.115, 0.115, 0.115, 0.115, 0.192, 0.192, 0.192, 0.192, 1.000,
  0.115, 0.115, 0.115, 0.115, 0.192, 0.192, 0.192, 0.192, 0.640, 1.000,
  0.115, 0.115, 0.115, 0.115, 0.192, 0.192, 0.192, 0.192, 0.640, 0.640,
  1.000,
  0.115, 0.115, 0.115, 0.115, 0.192, 0.192, 0.192, 0.192, 0.640, 0.640,
  0.640, 1.000,
  0.000, 0.000, 0.000, 0.000, 0.271, 0.271, 0.271, 0.271, 0.325, 0.325,
  0.325, 0.325, 1.000,
  0.000, 0.000, 0.000, 0.000, 0.271, 0.271, 0.271, 0.271, 0.325, 0.325,
  0.325, 0.325, 0.300, 1.000,
  0.000, 0.000, 0.000, 0.000, 0.271, 0.271, 0.271, 0.271, 0.325, 0.325,
  0.325, 0.325, 0.300, 0.300, 1.000
)
X <- X + t(X) - diag(diag(X)) 

colnames(X) <- rownames(X) <- c(paste("Y",1:12,sep=""),paste("X",1:3,sep=""))

# Print the population covariance matrix X to see that it is correct

X

# The estimated model in Lavaan syntax

analyzeModel <- "
ksi =~ Y1 + Y2 + Y3 + Y4
eta1 <~ X1 + X2 + X3
eta2 =~ Y5 + Y6 + Y7 + Y8
eta3 =~ Y9 + Y10 + Y11 + Y12

eta1 ~ ksi
eta2 ~ eta1
eta3 ~ eta1
"

# Only run 100 replications without bootstrap replications each so that the 
# example runs faster. Generate the data outside simsem because simsem
# does not support drawing samples from population matrix

dataSets <- lapply(1:100, function(x){
  mvrnorm(300,                 # Sample size
          rep(0,15),           # Means
          X)                   # Population covarancematrix
})

Output <- matrixpls.sim(model = analyzeModel, rawData = dataSets, boot.R=FALSE,
                        multicore = FALSE, stopOnError = TRUE)

summary(Output)


}



cleanEx()
nameEx("matrixplsboot")
### * matrixplsboot

flush(stderr()); flush(stdout())

### Name: matrixplsboot
### Title: Bootstrapping of matrixpls function
### Aliases: matrixplsboot matrixpls.boot summary.matrixplsboot

### ** Examples

# Run the customer satisfaction example form plspm

# load dataset satisfaction
data(satisfaction)
# inner model matrix
IMAG = c(0,0,0,0,0,0)
EXPE = c(1,0,0,0,0,0)
QUAL = c(0,1,0,0,0,0)
VAL = c(0,1,1,0,0,0)
SAT = c(1,1,1,1,0,0)
LOY = c(1,0,0,0,1,0)
inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
colnames(inner) <- rownames(inner)

# Reflective model

reflective<- matrix(
  c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1),
  27,6, dimnames = list(colnames(satisfaction)[1:27],colnames(inner)))

# empty formative model

formative <- matrix(0, 6, 27, dimnames = list(colnames(inner),
                                              colnames(satisfaction)[1:27]))

# Only 100 replications to make the example faster when running on single processor. 
# Increase the number of replications (e.g. to 500) to get the BCa intervals

matrixpls.boot.out <- matrixpls.boot(satisfaction[,1:27],
                           model = list(inner = inner,
                                        reflective = reflective,
                                        formative = formative),
                           R = 100)

print(matrixpls.boot.out)
summary(matrixpls.boot.out)



cleanEx()
nameEx("predict.matrixpls")
### * predict.matrixpls

flush(stderr()); flush(stdout())

### Name: predict.matrixpls
### Title: Predict method for matrixpls results
### Aliases: predict.matrixpls

### ** Examples

# Run the customer satisfaction example form plspm

# load dataset satisfaction
data(satisfaction)
# inner model matrix
IMAG = c(0,0,0,0,0,0)
EXPE = c(1,0,0,0,0,0)
QUAL = c(0,1,0,0,0,0)
VAL = c(0,1,1,0,0,0)
SAT = c(1,1,1,1,0,0)
LOY = c(1,0,0,0,1,0)
inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
colnames(inner) <- rownames(inner)

# Reflective model

reflective<- matrix(
  c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1),
  27,6, dimnames = list(colnames(satisfaction)[1:27],colnames(inner)))

# empty formative model

formative <- matrix(0, 6, 27, dimnames = list(colnames(inner),
                                              colnames(satisfaction)[1:27]))

satisfaction.model <- list(inner = inner,
                           reflective = reflective,
                           formative = formative)

# Estimation using covariance matrix


satisfaction.out <- matrixpls(cov(satisfaction[,1:27]),
                           model = satisfaction.model)

print(satisfaction.out)



# Predict indicators using means from the data
predict(satisfaction.out, 
        newData = satisfaction,
        means= sapply(satisfaction, mean))

# Calculate composite scores
predict(satisfaction.out, 
        newData = satisfaction,
        predictionType = "composites")



cleanEx()
nameEx("weightFun")
### * weightFun

flush(stderr()); flush(stdout())

### Name: weightFun
### Title: Indicator weight algoritms
### Aliases: weightFun weightFun.pls weightFun.optim weightFun.fixed
###   weightFun.factor weightFun.principal

### ** Examples

# Run the customer satisfaction examle form plspm

# load dataset satisfaction
data(satisfaction)
# inner model matrix
IMAG = c(0,0,0,0,0,0)
EXPE = c(1,0,0,0,0,0)
QUAL = c(0,1,0,0,0,0)
VAL = c(0,1,1,0,0,0)
SAT = c(1,1,1,1,0,0)
LOY = c(1,0,0,0,1,0)
inner = rbind(IMAG, EXPE, QUAL, VAL, SAT, LOY)
colnames(inner) <- rownames(inner)

# Reflective model
list(1:5, 6:10, 11:15, 16:19, 20:23, 24:27)

reflective<- matrix(
  c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1),
  27,6, dimnames = list(colnames(satisfaction)[1:27],colnames(inner)))

# empty formative model

formative <- matrix(0, 6, 27, dimnames = list(colnames(inner), colnames(satisfaction)[1:27]))

# Estimation using covariance matrix
model <- list(inner = inner,
              reflective = reflective,
              formative = formative)

S <- cov(satisfaction[,1:27])

matrixpls.ModeA <- matrixpls(S, model)
matrixpls.ModeB <- matrixpls(S, model, outerEstim = outerEstim.modeB)
matrixpls.MaxR2 <- matrixpls(S, model, weightFun = weightFun.optim)

# Compare the R2s from the different estimations

R2s <- cbind(r2(matrixpls.ModeA), r2(matrixpls.ModeB), r2(matrixpls.MaxR2))
print(R2s)
apply(R2s,2,mean)

# Optimization against custom function

maximizeSumOfCorrelations <- function(matrixpls.res){
  C <- attr(matrixpls.res,"C")
  model <- attr(matrixpls.res,"model")
  - sum(C[model$inner != 0])
}

matrixpls.MaxCor <- matrixpls(S, model, weightFun = weightFun.optim,
                             optimCrit = maximizeSumOfCorrelations)

# Compare the Mode B and optimized solutions

C <- attr(matrixpls.ModeB,"C")
print(C)
print(sum(C[inner != 0]))
C <- attr(matrixpls.MaxCor,"C")
print(C)
print(sum(C[inner != 0]))



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
