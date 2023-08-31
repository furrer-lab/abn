## ----------------------------------------------------------------------------------------------------------
## R script for reproducing results from paper: "Additive Bayesian Network Modelling with the R Package abn"

## -----------------------------------
## Author: Gilles Kratzer
## Date: 24.06.2020
## contact: gilles.kratzer@gmail.com

## -----------------------------------
## structure of the document:
# - data folder: original data.
# - input folder: model bug input.
# - output folder: Figures (FIG_X.pdf) + Tables (TAB_X.txt) + code output (CODE_X.txt) X is a decription of the object. They are computed by order of appearence in the paper. Addtional intermediate files are stored.
# - Rdat folder: storing intermediate results.

# Note:
# - Figure 9 entitled: Final ABN model with arc width proportional to link strength, needs post processing and can be manually generated with files: CODE_link_strength.txt (for the link strength) and TrimDAG.dot (for the structure).

## ------------------------------------------------------------------------------------------------------------
## loading libraries

## recommended installation
#install.packages('abn')
#install.packages("INLA", repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

library("abn")
library("bnlearn")
library("Rgraphviz")
library("knitr")
library("kableExtra")
library("rjags")
library("ggplot2")
library("reshape2")
library("ggplot2")
library("akima")
library("fields")
library("dplyr")

## ---------------------------------------------------------------------
## ASIA analysis: Motivating example
## ---------------------------------------------------------------------

#load data 
data("asia", package = "bnlearn")
colnames(asia) <- c("Asia", "Smoking", "Tuberculosis", "LungCancer", 
                    "Bronchitis", "Either", "XRay", "Dyspnea")


distrib <- as.list( rep("binomial", 8))
distrib.g <- as.list( rep("gaussian", 8))
names( distrib) <- names( distrib.g) <- names(asia)

## defining the true asia network using adjency matrix
dag.adj<-plotabn(dag.m =~Asia|Tuberculosis+
                   Tuberculosis|Either+
                   Either|XRay:Dyspnea+
                   Smoking|Bronchitis:LungCancer+
                   LungCancer|Either+
                   Bronchitis|Dyspnea,
                 data.dists = distrib.g,
                 edgedir = "cp",
                 plot = FALSE)

mycache <- buildscorecache(data.df = asia, data.dist = distrib, 
                           max.parents = 4)
mp.dag <- mostprobable(score.cache = mycache)

fabn <- fitabn(object = mp.dag, create.graph = TRUE)

# output
pdf('output/FIG_dag_asia_wo_constraints.pdf')
plotabn(dag.m = mp.dag$dag,data.dists = distrib.g,node.fillcolor = "white",edgedir = "pc",fontsize.node = 16)
dev.off()

#comparing DAGs

compareDag(ref = t(dag.adj),
           test = mp.dag$dag)

start_time <- Sys.time()
mycache <- buildscorecache(data.df = asia, data.dist = distrib, 
                           max.parents = 4, dag.retained = ~ LungCancer|Smoking)
mp.dag <- mostprobable(score.cache = mycache)
fabn <- fitabn(object = mp.dag, create.graph = TRUE)
end_time <- Sys.time()

abn.time <- end_time - start_time

#comparison with bnlearn 
start_time <- Sys.time()
fit = bn.fit(hc(asia,whitelist = data.frame(from = c("Smoking"),to = c("LungCancer"))),asia, method = "bayes")
end_time <- Sys.time()

bnlearn.time <- end_time - start_time

# output
pdf('output/FIG_dag_asia_w_constraints.pdf')
plotabn(dag.m = mp.dag$dag,data.dists = distrib.g,node.fillcolor = "white",edgedir = "pc",fontsize.node = 16)
dev.off()

compareDag(ref = t(dag.adj),
           test = mp.dag$dag)

# output: 
capture.output(fabn, file = "output/CODE_fabn_asia.txt")


## ---------------------------------------------------------------------
## Case study: adg
## ---------------------------------------------------------------------

# load data
dt <- readRDS("data/dataOK.RDS")
drop <- which(colnames(dt) %in% c("pneum", "epg5", "worms"))
drop.clustring <- which(colnames(dt) %in% c("pneum", "epg5", "worms")) #used with cluster correction
abndata <- dt[, -drop]

mm <- data.frame(Variable = colnames(abndata),
                 Meaning = c("presence of atrophic rhinitis",
                             "presence of moderate to severe pneumonia",
                             "sex of the pig (1=female, 0=castrated)",
                             "presence of liver damage (parasite-induced white spots)", 
                             "presence of fecal/gastrointestinal nematode eggs at time of slaughter",
                             "count of nematodes in small intestine at time of slaughter",
                             "days elapsed from birth to slaughter (days)",
                             "average daily weight gain (grams)",
                             "farm ID"), Distribution = c("Binomial", "Binomial","Binomial","Binomial","Binomial", "Poisson", "Continuous", "Continuous", "Discrete"))

out <- kable(mm, caption = "\\label{tab:overview}Description of the variables.", row.names = FALSE, digits = 2, align = "ll", "latex", booktabs = TRUE) %>% kable_styling(latex_options = c("hold_position"), full_width = FALSE, position = "left")


## output: Description of the variables.
capture.output(out, file = "output/TABLE_desc_var_case_study.txt")

## output: Descriptive distributions of the variables.
pdf("output/FIG_descr_distr_case_study.pdf")
par(mfrow=c(3,3), mar=c(2,4,1.5,1))
xx <- barplot(table(dt$AR)/341, ylim=c(0,1), main="AR", ylab="proportion")
text(x = xx, y = table(dt$AR)/341, label = table(dt$AR), pos = 3, cex = 0.8, col = "red")
xx <- barplot(table(dt$pneumS)/341, ylim=c(0,1), main="pneumS", ylab="proportion")
text(x = xx, y = table(dt$pneumS)/341, label = table(dt$pneumS), pos = 3, cex = 0.8, col = "red")
xx <- barplot(table(dt$female)/341, ylim=c(0,1), main="female", ylab="proportion")
text(x = xx, y = table(dt$female)/341, label = table(dt$female), pos = 3, cex = 0.8, col = "red")
xx <- barplot(table(dt$livdam)/341, ylim=c(0,1), main="livdam", ylab="proportion")
text(x = xx, y = table(dt$livdam)/341, label = table(dt$livdam), pos = 3, cex = 0.8, col = "red")
xx <- barplot(table(dt$eggs)/341, ylim=c(0,1), main="eggs", ylab="proportion")
text(x = xx, y = table(dt$eggs)/341, label = table(dt$eggs), pos = 3, cex = 0.8, col = "red")
hist(dt$wormCount, main="worms",prob=TRUE,col="grey",border="white", ylim= c(0,0.6))
lines(density(dt$wormCount),lwd=1.5)
hist(dt$age, xlab="", main="age",prob=TRUE,col="grey",border="white")
lines(density(dt$age),lwd=1.5)
hist(dt$adg, xlab="", main="adg",prob=TRUE,col="grey",border="white")
lines(density(dt$adg),lwd=1.5)
barplot(table(dt$farm), main="Farm ID", ylim=c(0,40), ylab="count", col.main = "gray50")
dev.off()

dist <- list(AR = "binomial", pneumS = "binomial", female = "binomial", 
             livdam = "binomial", eggs = "binomial", wormCount = "poisson",
             age = "gaussian", adg = "gaussian")


#data manipulation
abndata <- dt[, -drop]
abndata[,1:5] <- as.data.frame(lapply(abndata[,1:5], factor))
abndata <- abndata[,-ncol(abndata)]

#ban/retain matrices
retain <- matrix(0, ncol(abndata), ncol(abndata))
colnames(retain) <- rownames(retain) <- names(abndata)

banned <- matrix(0, ncol(abndata), ncol(abndata))
colnames(banned) <- rownames(banned) <- names(abndata)

banned[3,-3] <- 1

result <- vector("numeric", 7)

# ban: formula statement
# retain: not constrained
for (i in 1:7) {
  mycache <- buildscorecache(data.df = as.data.frame(abndata), data.dists = dist, 
                             dag.banned = ~female|., dag.retained = NULL, 
                             max.parents = i,method = "bayes")
  mydag <- mostprobable(score.cache = mycache)
  result[i] <- fitabn(object = mydag)$mlik
}


datadir <- tempdir() 

for (i in 1:7) {
  max.par <- i
  mycache <- buildscorecache(data.df = as.data.frame(abndata), data.dists = dist, 
                             dag.banned = banned, dag.retained = retain, 
                             max.parents = max.par)
  mydag <- mostprobable(score.cache = mycache)
  fabn <- fitabn(object = mydag)
  cat(paste("network score for", i, "parents =", fabn$mlik, "\n\n"))
  save(mycache, mydag, fabn, file = paste("Rdat/", "mp_", max.par, ".RData", sep = ""))
}


# get network score for all parent limits

mp.mlik <- c()
for (i in 1:max.par) {
  load(paste("Rdat/","mp_", i,".RData", sep=""))
  mp.mlik <- c(mp.mlik, fabn$mlik)
}

## output: Total network log marginal likelihood as a function of the number of parents
pdf("output/FIG_parents_limit.pdf")
plot(1:max.par, mp.mlik, xlab = "Parent limit", ylab = "Log marginal likelihood", 
     type = "b", col="red", ylim=range(mp.mlik))
abline(v=which(mp.mlik==max(mp.mlik))[1], col="grey", lty=2)
dev.off()



max.par<- which(mp.mlik==max(mp.mlik))[1]
load(paste("Rdat/","mp_", max.par,".RData", sep=""))

# save dot file
tographviz(dag.m = mydag$dag, data.df = abndata, data.dists = dist,
           outfile = paste0("output/DAG_", which(mp.mlik==max(mp.mlik))[1],"p.dot"))

#output: DAG selected using an exact search with a model complexity of four parents.
pdf(file = 'output/FIG_dag_4p.pdf')
plotabn(dag.m = mydag$dag,data.dists = dist,node.fillcolor = "white",edgedir = "pc")

dev.off()
# system specific! e.g. on linux ->
# system("dot -Tpdf -o graph.pdf graph.dot")
# system("evince graph.pdf")


# Fit marginal densities over a fixed grid --> n.grid=1000
# --------------------------------------------------------
marg.f <- fitabn(object = mydag, compute.fixed=TRUE, n.grid=100)

# Extract values 
# --------------
m <- marg.f$marginals[[1]] 
for(i in 2: length(marg.f$marginals)){ 
  m <- c(m, marg.f$marginals[[i]])}

# Bind all the marginals for the same node into a matrix
# ------------------------------------------------------
AR.p <- cbind( m[[ "AR|(Intercept)"]], m[[ "AR|age"]])
pneumS.p <- cbind( m[[ "pneumS|(Intercept)"]], m[[ "pneumS|age"]])
female.p <- cbind( m[[ "female|(Intercept)"]])
livdam.p <- cbind( m[[ "livdam|(Intercept)"]], m[[ "livdam|eggs"]])
eggs.p <- cbind( m[[ "eggs|(Intercept)"]], m[[ "eggs|adg"]])
wormCount.p <- cbind( m[[ "wormCount|(Intercept)"]], m[[ "wormCount|AR"]],
                      m[[ "wormCount|eggs"]], m[[ "wormCount|age"]], m[[ "wormCount|adg"]])
age.p <- cbind( m[[ "age|(Intercept)"]], m[[ "age|female"]])
prec.age.p <- cbind( m[[ "age|precision" ]])
adg.p <- cbind( m[[ "adg|(Intercept)"]], m[["adg|age"]])
prec.adg.p <- cbind( m[[ "adg|precision" ]])

# Save it to a file named PostParams to be read by JAGS
# -----------------------------------------------------
dump(c("AR.p", "pneumS.p", "female.p", "livdam.p", "eggs.p", 
       "wormCount.p", "age.p", "prec.age.p", "adg.p", "prec.adg.p"),
     file="Rdat/PostParams.R")

# set inits
# ---------
init <- list(".RNG.name"="base::Mersenne-Twister", 
             ".RNG.seed"=122)

# load data
# ---------
source("Rdat/PostParams.R")

# run model once
# --------------
jj <- jags.model(file = "input/model8vPois.bug", 
                 data = list('AR.p'=AR.p, 'pneumS.p'=pneumS.p,
                             'female.p'=female.p, 'livdam.p'=livdam.p,
                             'eggs.p'=eggs.p, 'wormCount.p'=wormCount.p , 
                             'age.p'=age.p,'prec.age.p'=prec.age.p, 
                             'adg.p'=adg.p, 'prec.adg.p'=prec.adg.p),
                 inits = init,
                 n.chains = 1, 
                 n.adapt = 50)

# run more iterations
# -------------------
update(jj, 100000)

# set number of observation we want to extract for a dataset 
# ----------------------------------------------------------
n.obs=341

# sample data (same size as original: 341) 
# with a sampling lag (20) to reduce avoid autocorrelation
# --------------------------------------------------------
samp <- coda.samples(jj, c("AR", "pneumS", "female", 
                           "livdam", "eggs", "wormCount", 
                           "age", "prec.age", 
                           "adg", "prec.adg"),
                     n.iter= n.obs*20 , thin =20)


# extract posterior densities and put in a dataframe
# --------------------------------------------------
post.df <- data.frame(AR = unlist(samp[,"AR"]),
                      pneumS = unlist(samp[,"pneumS"]),
                      female = unlist(samp[,"female"]),
                      livdam = unlist(samp[,"livdam"]),
                      eggs = unlist(samp[,"eggs"]),
                      wormCount = unlist(samp[,"wormCount"]),
                      age = unlist(samp[,"age"]),
                      adg = unlist(samp[,"adg"]))


# compare with original data
# --------------------------
df <- abndata

# get centered version of age and adg to compare to bootstrap data
# ----------------------------------------------------------------
df$age.c <- scale(df$age)
df$adg.c <- scale(df$adg)

# Compare distribution of original and simulated data
# ---------------------------------------------------

##output: Example of simulated data (in blue) versus original data (in black)
# ---------------------------------------------------------------------------
pdf("output/FIG_simulated_data_by_bootstrapping.pdf")
par(mfrow=c(4,4), mar=c(2,2,1.5,1))

barplot(table(df$AR)/n.obs, ylim=c(0,1), main="AR - original")
barplot(table(post.df$AR)/n.obs,  ylim=c(0,1), main="AR - simulated", 
        col.main = "blue", border="blue") 

barplot(table(df$pneumS)/n.obs, ylim=c(0,1), main="pneumS - original")
barplot(table(post.df$pneumS)/n.obs,  ylim=c(0,1), main="pneumS - simulated", 
        col.main = "blue", border="blue") 

barplot(table(df$female)/n.obs, ylim=c(0,1), main="female - original")
barplot(table(post.df$female)/n.obs,  ylim=c(0,1), main="female - simulated", 
        col.main = "blue", border="blue") 

barplot(table(df$livdam)/n.obs, ylim = c(0, 1), main = "livdam - original")
barplot(table(post.df$livdam)/n.obs,  ylim = c(0, 1), main = "livdam - simulated", 
        col.main = "blue", border = "blue")

barplot(table(df$eggs)/n.obs, ylim = c(0, 1), main = "eggs - original")
barplot(table(post.df$eggs)/n.obs, ylim = c(0, 1), main = "eggs - simulated", 
        col.main = "blue", border = "blue") 

hist(df$wormCount, xlab = "", main = "wormCount - original",
     prob=TRUE, col = "grey", border = "white", ylim = c(0, 0.6))
lines(density(df$wormCount),lwd=1.5)
hist(post.df$wormCount, xlab = "", main="wormCount - simulated", col.main = "blue", 
     prob=TRUE, col = "grey", border = "white", xlim = c(0, 80), ylim = c(0, 0.6))
lines(density(post.df$wormCount), lwd = 1.5, col = "blue")

hist(df$age.c, xlab="", main="age - original",
     prob=TRUE,col="grey",border="white")
lines(density(df$age.c),lwd=1.5)
hist(post.df$age, xlab="", main="age - simulated", col.main = "blue",
     prob=TRUE,col="grey",border="white")
lines(density(post.df$age),lwd=1.5, col="blue")

hist(df$adg.c, xlab="", main="adg - original",
     prob=TRUE,col="grey",border="white")
lines(density(df$adg.c),lwd=1.5)
hist(post.df$adg, xlab="", main="adg - simulated", col.main = "blue",
     prob=TRUE,col="grey",border="white")
lines(density(post.df$adg),lwd=1.5, col="blue")
dev.off()

vars <- colnames(abndata)

## compute 
# --------

# select nr. bootstrap samples to run
# -----------------------------------
set.seed(123)

# get 5000 random numbers to set different initial values
# -------------------------------------------------------
n <- sample(1:100000, 5000)

# specify max number of parents based on previous search
# ------------------------------------------------------
max.par <- 4

# Simulate data and run ABN on such dataset
# -----------------------------------------

for (i in 1:length(n)) {
  
  print(paste("/n Running simulation", i))
  
  # pick initials
  init <- list(".RNG.name" = "base::Mersenne-Twister", ".RNG.seed" = n[i])
  
  # run model
  jj <- jags.model(file = "input/model8vPois.bug", 
                   data = list('AR.p' = AR.p, 'pneumS.p' = pneumS.p,
                               'female.p' = female.p, 'livdam.p' = livdam.p,
                               'eggs.p' = eggs.p, 'wormCount.p' = wormCount.p , 
                               'age.p' = age.p,'prec.age.p' = prec.age.p, 
                               'adg.p' = adg.p, 'prec.adg.p' = prec.adg.p),
                   inits = init,
                   n.chains = 1, 
                   n.adapt = 5000)
  
  # run more iterations
  update(jj, 100000)
  
  # sample data (same size as original: 341) with a sampling lag (20) to reduce autocorrelation
  samp <- coda.samples(jj, c("AR", "pneumS", "female", 
                             "livdam", "eggs", "wormCount", 
                             "age", "prec.age", 
                             "adg", "prec.adg"),
                       n.iter= n.obs*20 , thin =20)
  
  # build dataframe in the same shape as the original one
  dt.boot <- as.data.frame(as.matrix(samp)) # pay attention at order names
  
  dt.boot<- dt.boot[, vars]
  
  
  # now coerce to factors if need be and set levels
  abndata <- as.data.frame(abndata)
  
  for(j in 1:length(vars)) {
    if(is.factor(abndata[,j])) {
      dt.boot[,j] <- as.factor(dt.boot[,j]);
      levels(dt.boot[,j]) <- levels(abndata[,j]);
    }
  }
  
  # Build a cache of all local computations
  mycache <- buildscorecache(data.df = dt.boot, data.dists = dist, dag.banned = banned,
                             dag.retained = retain, max.parents = max.par)
  
  # Run the EXACT SEARCH
  mp.dag <- mostprobable(score.cache = mycache)
  fabn <- fitabn(object =  mp.dag)
  
  # Save the results obtained
  save(mycache, mp.dag, fabn, dt.boot, file = sprintf('Rdat/BootData/dt.boot.%04d.RData', i))
}

n <- 5000

# load dags and boostrap data
# ------------------------------
dags <- list()
boot <- list()
for (i in 1:n) {
  load(file = sprintf('Rdat/BootData/dt.boot.%04d.RData',i))
  dags[[i]] <- mp.dag$dag
  boot[[i]] <- dt.boot
  rm(mp.dag, dt.boot, fabn, mycache)
} 

saveRDS(dags, file = "Rdat/BootData/BootDAGs5000.RDS")

arcs <- sapply(dags, sum)

# output: Histogram of the distribution of the number of arcs in the bootstrapped searches
# ----------------------------------------------------------------------------------------
pdf("output/FIG_hist_nb_searches.pdf")
barplot(table(arcs))
dev.off()


# Count how many times each arc appear in the bootstrap data
# ----------------------------------------------------------
# function Reduce sums ("+") each element of each maxtrix in the list and store 
# results in a new matrix of same size
alldag <- Reduce("+", dags)  

# express it in percentage
# ------------------------
perdag <- alldag/length(dags)

#output: 
# ------
capture.output(print(round(perdag*100, digits = 0)), file = "output/CODE_percentage_dag.txt")

trim.dag <- (alldag >= (n*0.5))*1


# Send  final pruned DAG to Graphvis for visualization
# -----------------------------------------------------
tographviz(dag.m = trim.dag, data.df = abndata, data.dists = dist,
           outfile = "output/TrimDAG.dot")

#plot daga
pdf('output/FIG_FINAL_DAG.pdf')
plotabn(dag.m = trim.dag,data.dists = dist,node.fillcolor = "white",edgedir = "pc")
dev.off()

# --------------------------------------------------------------------------
## Marginal densities of the model parameter NOT corrected for random effect
# --------------------------------------------------------------------------


marg.f.not.grouped.fitabn <- fitabn(dag.m = trim.dag, 
                                data.df = as.data.frame(abndata), data.dists = dist,
                                compute.fixed = TRUE,n.grid = 1000)

#output: Marginal densities of model parameter corrected for random effect.
# --------------------------------------------------------------------------

pdf("output/FIG_marginal_densities_wo_group_correction.pdf")
par(mfrow=c(5, 4), mar=c(2,2, 1.5, 1))
for(i in 1:length(marg.f.not.grouped.fitabn$marginals)){
  
  # get the marginal for current node, which is a matrix [x, f(x)]
  cur.node <- marg.f.not.grouped.fitabn$marginals[i]
  nom1 <- names(marg.f.not.grouped.fitabn$marginals)[i]
  
  # pick the first value (for models wothout random effects)
  cur.node <- cur.node[[1]]
  for(j in 1:length(cur.node) ) {
    nom2 <- names(cur.node)[j]
    cur.param <- cur.node[[j]]
    plot(cur.param, type = "l", main = paste(nom1, ":", nom2), cex = 0.7)
  }
}
dev.off()


# ------------------------------------------------------------------
## Marginal densities of model parameter corrected for random effect
# ------------------------------------------------------------------

# incorporate grouping factor in the data
#---------------------------------------- 
abndata <- dt[, -drop.clustring]

# set up factors
#---------------
abndata[,c(1:5,9)] <- as.data.frame(lapply(abndata[,c(1:5,9)], factor))
class(abndata) <- "data.frame"


if(FALSE){
# recompute the cache of scores using GLMM
#-----------------------------------------
mycache <- buildscorecache(data.df = as.data.frame(abndata), 
                           data.dists = dist,
                           group.var = "farm",
                           cor.vars = c("AR", "pneumS", "female",
                                        "livdam", "eggs", 
                                        "wormCount", "age", "adg"), 
                           dag.banned = ~female|., 
                           dag.retained = NULL,  
                           max.parents = 4)

# exact search
#-------------
dag.adjusted <- mostprobable(score.cache = mycache)

# recompute the marginal using GLMM
#----------------------------------
marg.f.grouped <- fitabn(object = dag.adjusted, 
                         group.var = "farm",
                         cor.vars = c("AR", "pneumS", "female", 
                                      "eggs", "wormCount", "age", "adg"),
                         compute.fixed = TRUE, 
                         n.grid = 100, 
                         control = list(max.mode.error = 0, 
                                        epsabs.inner = 1e-01, 
                                        max.hessian.error = 5e-01, 
                                        epsabs = 1e-1, 
                                        error.verbose = FALSE, 
                                        hessian.params = c(1E-02, 1E-01), 
                                        factor.brent = 1E+01, 
                                        loggam.inv.scale = 10e-02))


# save computed marginals
# -----------------------
save(file = "Rdat/marg.grouped.RData",marg.f.grouped)

# alrady computed marginals
# -------------------------
load(file = "Rdat/marg.grouped.RData")


par(mfrow=c(7, 5), mar=c(2,2, 1.5, 1))
for(i in 1:length(marg.f.grouped$marginals)){
  
  # get the marginal for current node, which is a matrix [x, f(x)]
  cur.node <- marg.f.grouped$marginals[i]
  nom1 <- names(marg.f.grouped$marginals)[i]
  
  # pick the first value (for models without random effects)
  cur.node <- cur.node[[1]]
  for(j in 1:length(cur.node) ) {
    nom2 <- names(cur.node)[j]
    cur.param <- cur.node[[j]]
    plot(cur.param, type = "l", main = paste(nom1, ":", nom2), cex = 0.7)
  }
}
}

marg.f.grouped.fitabn <- fitabn(dag.m = trim.dag, 
                                data.df = as.data.frame(abndata), data.dists = dist, 
                         group.var = "farm",
                         cor.vars = c("AR", "pneumS", #"female", 
                                      "livdam", "eggs", "wormCount", 
                                      "age", "adg"),
                         compute.fixed = TRUE,n.grid = 1000)

marg.f.grouped <- marg.f.grouped.fitabn
save(file = "Rdat/marg.grouped.fitabn.RData",marg.f.grouped, trim.dag, abndata, dist)

#load(file = "Rdat/marg.grouped.fitabn.RData")

#output: Marginal densities of model parameter corrected for random effect.
pdf("output/FIG_marginal_densities_group_correction.pdf")
par(mfrow=c(7, 4), mar=c(2,2, 1.5, 1))
for(i in 1:length(marg.f.grouped$marginals)){
  
  # get the marginal for current node, which is a matrix [x, f(x)]
  cur.node <- marg.f.grouped$marginals[i]
  nom1 <- names(marg.f.grouped$marginals)[i]
  
  # pick the first value (for models wothout random effects)
  cur.node <- cur.node[[1]]
  for(j in 1:length(cur.node) ) {
    nom2 <- names(cur.node)[j]
    cur.param <- cur.node[[j]]
    plot(cur.param, type = "l", main = paste(nom1, ":", nom2), cex = 0.7)
  }
}
dev.off()

#extract marginals from non corrected model
marg.dens <- marg.f$marginals[[1]]
for (i in 2:length(marg.f$marginals)) {
  marg.dens <- c(marg.dens, marg.f$marginals[[i]])
}

# extract marginals adjusted for grouped data
marg.dens.grouped <- marg.f.grouped$marginals[[1]]
for (i in 2:length(marg.f.grouped$marginals)) {
  marg.dens.grouped <- c(marg.dens.grouped, marg.f.grouped$marginals[[i]])
}

# calculate AUC --> should be ~1
auc <- rep(NA, length(marg.dens.grouped))
names(auc) <- names(marg.dens.grouped)
for(i in 1:length(marg.dens.grouped)) {
  tmp <- spline(marg.dens.grouped[[i]])
  auc[i] <- sum(diff(tmp$x[-length(tmp$x)]) * tmp$y[-(1:2)])
}


mat <- matrix(rep(NA, length(marg.dens.grouped)*3), ncol = 3)
rownames(mat) <- names(marg.dens.grouped)
colnames(mat) <- c("2.5%", "50%", "97.5%")
ignore.me <- union(union(grep("\\(Int", names(marg.dens.grouped)), grep("prec", names(marg.dens.grouped))), grep("group.precision", names(marg.dens.grouped))) #
comment <- rep("", length(marg.dens.grouped))
for (i in 1:length(marg.dens.grouped)) {
  tmp <- marg.dens.grouped[[i]]
  tmp2 <- cumsum(tmp[,2])/sum(tmp[,2])
  mat[i, ] <-c(tmp[which(tmp2 > 0.025)[1],1],## -1 is so use value on the left of the 2.5%
               tmp[which(tmp2 > 0.5)[1],1],
               tmp[which(tmp2 > 0.975)[1],1])
  vec <- mat[i,]
  
  if( !(i%in%ignore.me) && (vec[1]<0 && vec[3]>0)){
    comment[i] <- "not sig. at 5%"
  }
  
  ## truncate for printing
  mat[i,] <- as.numeric(formatC(mat[i,], digits = 3, format = "f"))
}

mar <- data.frame(cbind(mat))
names(mar) <- c("2.5Q", "median", "97.5Q")

mar <- mar[-ignore.me,]
mar[1:7,] <- exp(mar[1:7,])
mar$interpretation <- c(rep("odds ratio", 3),
                        rep("rate ratio", 4),
                        rep("correlation", 2))

#output: Marginals posterior distribution of the parameter estimates
capture.output(mar, file = "output/TABLE_marginal_posterior_densities.txt")


#data manipulation
abndata <- dt[, -drop]
abndata[,1:5] <- as.data.frame(lapply(abndata[,1:5], factor))
abndata <- abndata[,-ncol(abndata)]



## link strength
## -------------

LS <- linkStrength(dag = trim.dag, data.df = as.data.frame(abndata),
                   data.dists = dist, method="ls.pc")
rownames(LS) <- colnames(LS) <- names(abndata)


#output:
capture.output(print(LS, digits = 3), file = "output/CODE_link_strength.txt")


mar$LS <- c(LS["AR", "age"], LS["livdam", "eggs"], LS["eggs", "adg"],
            LS["wormCount", "AR"], LS["wormCount", "eggs"], LS["wormCount", "age"], 
            LS["wormCount", "adg"], LS["age", "female"], LS["adg", "age"])


# output: 
out <- kable(mar, caption = "Marginals posterior distribution of the parameter estimates", row.names = TRUE, digits = 2, align = "rrrrr", "latex", booktabs = TRUE) %>%  kable_styling(latex_options = c("hold_position","basic"), full_width = FALSE, position = "center")

capture.output(out, file = "output/TABLE_marg_post_distr.txt")


##------------------------------------------------------------------------------------
## simulation study
##------------------------------------------------------------------------------------



compareEG <- function(ref, test, node.names = NULL) {
  
  if (any(dim(ref) != dim(test))) {
    stop("The reference or test DAG has not the same size")
  }
  
  n <- dim(ref)[1]
  
  ## unit matrix
  ref[ref != 0] <- 1
  test[test != 0] <- 1
  
  diff.matrix <- ref - (0.5 * test)
  
  diff.matrix.tab <- table(factor(diff.matrix, levels = c(-0.5, 0, 0.5, 1)))
  
  if(sum(ref == 1)==0 | sum(test == 1)==0){
    warning("If the test or reference matrix is an empty matrix some of the estimates are not defined.")
  }
  
  ## output
  out <- list()
  
  out[["TPR"]] <- (as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5]))/(sum(ref == 1))
  out[["FPR"]] <- (as.numeric(diff.matrix.tab[names(diff.matrix.tab) == -0.5]))/(sum(ref == 0))
  out[["Accuracy"]] <- (as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5]) + as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0]))/(dim(ref)[1]^2)
  out[["FDR"]] <- as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 1])/(sum(test == 1))
  out[["G-measure"]] <- sqrt(as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5])/(sum(test == 1)) * (as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5]))/(sum(ref == 1)))
  out[["F1-score"]] <- (2/((1/as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5])/(sum(test == 1))) + (1/(as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5]))/(sum(ref == 1)))))
  out[["PPV"]] <- as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 0.5])/(sum(test == 1))
  out[["FOR"]] <- as.numeric(diff.matrix.tab[names(diff.matrix.tab) == 1])/(sum(test == 1))
  
  # transforming "reverse arc" into -0.5
  for( i in 1:n){
    for( j in i:n){
      if( diff.matrix[i,j]!=0 & diff.matrix[j,i]!=0){
        diff.matrix[i,j] <- -(diff.matrix[i,j] + diff.matrix[j,i])
        diff.matrix[j,i] <- 0
      }
    }
  }
  
  diff.matrix.tab <- table(factor(diff.matrix, levels = c(-0.5, 0, 0.5, 1)))
  
  out[["Hamming-distance"]] <- sum(as.numeric(diff.matrix.tab[names(diff.matrix.tab) %in% c(-0.5, 1)]))
  
  return(out)
}

#distribution list
dist <- list(a="gaussian",b="gaussian",c="gaussian",d="gaussian", e="gaussian", f="gaussian", g="gaussian", h="gaussian", i="gaussian", j="gaussian")

#define parameter matrix
data.param <- matrix(data = c(0,0.2,0.5,0,0,0,0,0.8,0.9,0.5,
                              0,0,0.3,0.1,0,0.8,0.6, 0, 0, 0.6,
                              0,0,0,0,0,0,0,0,0,0.5,
                              0,0,0,0,0,0,0,0,0,0,
                              0.1,0,0,0,0,0,0.8,0.2,0.3,0,
                              0,0,0,0,0,0,0,0,0,0,
                              0,0,0.3,0.8,0,0.8,0, 0, 0, 0.6,
                              0,0,0.3,0.6,0,0.8,0.6, 0, 0, 0.6,
                              0,0,0.3,0.7,0,0.8,0.6, 0, 0, 0.6,
                              0,0,0,0,0,0,0,0,0,0),nrow = 10L,ncol = 10L,byrow = TRUE)

#precision matrix
data.param.var <- matrix(data = 0,nrow = 10L,ncol = 10L)
diag(data.param.var) <- c(10,20,30,40,10,10,20,20,10,10)

#naming
colnames(data.param) <- rownames(data.param) <- colnames(data.param.var) <- rownames(data.param.var) <- names(dist)

grid_abn <- expand.grid(nobs=c(5,10,50,100,500,1000,2000,3000,5000, 7500, 10000), rep = 1:20) 

out <- list()

seed <- 35674

for(i in 1:dim(grid_abn)[1]){
  simulated.data <- simulateAbn(data.dists = dist,n.chains = 1,n.thin = 1,n.iter = grid_abn[i,1],
                                data.param = data.param,data.param.var = data.param.var,
                                data.param.mult = data.param.mult,seed = seed)
  
  seed<-seed+1
  
  mycache.b <- buildscorecache(data.df = simulated.data,data.dists=dist,max.parents = 5,dry.run = FALSE,centre = FALSE)
  dag.b <- mostprobable(score.cache = mycache.b)
  
  t0 <- essentialGraph(dag = (dag.b$dag))
  
  
  mycache.computed.mle <- buildscorecache(data.df = simulated.data,data.dists = dist,max.parents = 5,method = "mle",centre = FALSE)
  dag.mle.mlik <- mostprobable(score.cache = mycache.computed.mle,score = "mlik")
  t1 <- essentialGraph(plotabn(dag.m = dag.mle.mlik$dag,data.dists = dist,plot = FALSE))
  
  dag.mle.aic <- mostprobable(score.cache = mycache.computed.mle,score = "aic")
  t2 <- essentialGraph(plotabn(dag.m = dag.mle.aic$dag,data.dists = dist,plot = FALSE))
  
  dag.mle.bic <- mostprobable(score.cache = mycache.computed.mle,score = "bic")
  t3 <- essentialGraph(plotabn(dag.m = dag.mle.bic$dag,data.dists = dist,plot = FALSE))
  
  dag.mle.mdl<-mostprobable(score.cache = mycache.computed.mle,score = "mdl")
  t4 <- essentialGraph(plotabn(dag.m = dag.mle.mdl$dag,data.dists = dist,plot = FALSE))
  
  abn <- compareEG(ref = essentialGraph(data.param),test = t0)
  mlik <- compareEG(ref = essentialGraph(data.param),test = t1)
  aic <- compareEG(ref = essentialGraph(data.param),test = t2)
  bic <- compareEG(ref = essentialGraph(data.param),test = t3)
  mdl <- compareEG(ref = essentialGraph(data.param),test = t4)
  
  
  out <- rbind(out, c(abn,mlik,aic,bic,mdl))
  #return(out)
}

out <- cbind(out, grid_abn)
df <- as.data.frame(out)
save(df,file = "Rdat/data_sim.rda",compress = TRUE)

## analysis

df.names <- df[,which(names(df) %in% c("TPR","FPR","Accuracy", "nobs", "rep"))]

names(df.names) <- c("abn.TP", "abn.FP","abn.Accuracy","mlik.TP","mlik.FP", "mlik.Accuracy","aic.TP","aic.FP","aic.Accuracy", "bic.TP","bic.FP","bic.Accuracy","mdl.TP","mdl.FP","mdl.Accuracy", "nobs","rep")

df.names[,] <- apply(df.names, 2, function(x) as.numeric(as.character(x)))


#save(df,out,file = "sim.dat")


df.long.scores.Accuracy <- df.names[,c(3,6,9,12,15,16,17)]
df.long.scores.FP<-df.names[,c(2,5,8,11,14,16,17)]
df.long.scores.TP<-df.names[,c(1,4,7,10,13,16,17)]

#df.long.scores.TN<-melt(df.long.scores.TN, id.vars=c("nobs", "rep"))
#df.long.scores.TP<-melt(df.long.scores.TP, id.vars=c("nobs", "rep"))
names(df.long.scores.Accuracy) <- tolower(names(df.long.scores.Accuracy))

df.long.scores.Accuracy <- melt(as.data.frame(df.long.scores.Accuracy), id=c("nobs", "rep"))
df.long.scores.FP<-melt(as.data.frame(df.long.scores.FP), id=c("nobs", "rep"))
df.long.scores.TP<-melt(as.data.frame(df.long.scores.TP), id=c("nobs", "rep"))

#luca
p1<-ggplot(data = df.long.scores.Accuracy)+
  geom_boxplot(aes(x=as.factor(nobs),y=value,fill=variable))+
  ylim(c(0.5,1))+ggtitle("Percentage of arcs retrieved (30% connected BN)")+xlab("")+ylab("")+theme(legend.title = element_blank())

p2<-ggplot(data = df.long.scores.Accuracy)+
  geom_smooth(aes(x=nobs,y=value,fill=variable, color = variable),level=0.95)+scale_x_log10(breaks=c(10,100,1000,10000))+
  ylim(c(0.5,1))+ggtitle("Percentage of arcs retrieved (30% connected BN)")+xlab("")+ylab("")+theme(legend.title = element_blank())

p3<-ggplot(data = df.long.scores.FP)+
  geom_boxplot(aes(x=as.factor(nobs),y=value,fill=variable))+
  ylim(c(0,0.5))+ggtitle("Percentage of false postive")+xlab("")+ylab("")+theme(legend.title = element_blank())

p4<-ggplot(data = df.long.scores.FP)+
  geom_smooth(aes(x=nobs,y=value,fill=variable, color = variable),level=0.95)+scale_x_log10(breaks=c(10,100,1000,10000))+
  ylim(c(0,0.5))+ggtitle("Percentage of false positive")+xlab("")+ylab("")+theme(legend.title = element_blank())

p5<-ggplot(data = df.long.scores.TP)+
  geom_boxplot(aes(x=as.factor(nobs),y=value,fill=variable))+
  ylim(c(0,1))+ggtitle("Percentage of true positive")+xlab("Number of Observations")+ylab("")+theme(legend.title = element_blank())

p6<-ggplot(data = df.long.scores.TP)+
  geom_smooth(aes(x=nobs,y=value,fill=variable, color = variable),level=0.95)+scale_x_log10(breaks=c(10,100,1000,10000))+
  ylim(c(0,1))+ggtitle("Percentage of true positive")+xlab("Number of Observations")+ylab("")+theme(legend.title = element_blank())

##Multiplot function

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library("grid")
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


pdf("output/FIG_boxplot_big_network.pdf")
multiplot(p1,p3,p5)
dev.off()


pdf("output/Fig_loess_network.pdf")
multiplot(p2,p4,p6)
dev.off()




# set.seed
set.seed(567)

out <- list()

grid_abn <- expand.grid(val = seq(0,1,0.1), rep = 1:1000)

# simulate BN and store metrics
#------------------------------

for(i in 1:dim(grid_abn)[1]){
  sim <- simulateDag(node.name = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t","u", "v", "w", "x", "y", "z", "aa", "ab", "ac", "ad", "ae", "af", "ag", "ah", "ai", "aj", "ak", "al", "am", "an"), nn = grid_abn$val[i])
  
  out <- rbind(out, (infoDag(dag = sim[[1]],node.names = names(sim[[2]]))))
}

out <- cbind(out[,-(1)],grid_abn)

# metrics normalization
#----------------------

out[,] = apply(out, 2, function(x) as.numeric(as.character(x)))
out[,1] <- out[,1]/780
out[,2] <- out[,2]/39
out[,3] <- out[,3]/39
out[,4] <- out[,4]/19.5
out[,5] <- out[,5]/19.5

# reshaping
#----------

out.long <- melt(out, id.vars = c("val", "rep"))

# plotting
#---------

p <- ggplot(data = out.long) +
  geom_boxplot(aes(x = as.factor(val), y = value, fill = variable)) +
  ylim(c(0, 1)) + 
  ggtitle("Density metric of DAG with 40 nodes (n=1000)") + 
  xlab("Normalized density of the DAGs") + 
  ylab("Normalized Maximum [%]")

# output
ggsave("output/FIG_normalized_DAG_metrics.pdf", plot = p)

## Comparison between Bayesian and MLE implementations to estimate regression coefficients in an ABN framework. The panels show the maximum Root Mean Squared Error (max RMSE) in function of the network density, the distribution variability and the sample size.
## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


##===========================================================================
## The purpose of this file is to asses efficence of coefficients estimations
## Linked to Paper Unknown structure
##===========================================================================

##===================================================
## Gaussian Network
##===================================================

##test coefficient of dags

##low density network

dist<-as.list(rep("gaussian", 20))
names(dist)<-c("a","b","c","d","e","f","g","h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t")

data.param<-simulateDag(node.name = c("a","b","c","d","e","f","g","h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t"),data.dists = dist,nn = 0.2)[[1]]


data.param.var<-matrix(data = 0,nrow = 20L,ncol = 20L)


colnames(data.param.var)<-rownames(data.param.var)<-names(dist)

#plotabn(dag.m = data.param,data.dists = dist)


lseq <- function(from=0.001, to=1000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

grid_abn<-expand.grid(nobs=c(25,50,100,500,1000,2000,3000,5000), var.var=lseq(), rep = 1:20) 

#

out.bays<-list()
out.mle<-list()

seed<-247913

for(i in 1:dim(grid_abn)[1]){
  
  
  diag(data.param.var)<-rep(grid_abn[i,2],20)
  
  simulated.data<-simulateAbn(data.dists = dist,n.chains = 1,n.thin = 1,n.iter = grid_abn[i,1],
                              data.param = data.param,data.param.var = data.param.var,
                              data.param.mult = data.param.mult,seed = seed)
  
  seed<-seed+1
  
  
  data.bays<-fitabn(dag.m = plotabn(dag.m = data.param,data.dists = dist,plot = FALSE),data.df = simulated.data,data.dists = dist,centre = FALSE)
  data.mle<-fitabn(dag.m = plotabn(dag.m = data.param,data.dists = dist,plot = FALSE),data.df = simulated.data,data.dists = dist,method = "mle")
  out.bays<-rbind(out.bays, c(unlist(data.bays$modes)))
  out.mle<-rbind(out.mle, c(unlist(data.mle$coef)))
}
#colnames(out.mle)<-colnames(out.bays)

##numeric
df.bayes = apply(out.bays, 2, function(x) as.numeric(as.character(x)))
df.mle = apply(out.mle, 2, function(x) as.numeric(as.character(x)))

df.bayes <- df.bayes[, -( grep(paste0("precision") , colnames(df.bayes),perl = TRUE) ) ]
colnames(df.mle)<-colnames(df.bayes)

df.bayes <- df.bayes[, -( grep(paste0("Intercept") , colnames(df.bayes),perl = TRUE) ) ]
df.mle <- df.mle[, -( grep(paste0("Intercept") , colnames(df.mle),perl = TRUE) ) ]

df.bayes <- cbind(grid_abn,df.bayes)
df.mle <- cbind(grid_abn,df.mle)

df_agg.bayes<-aggregate(x = df.bayes,by = list(df.bayes$nobs, df.bayes$var.var), FUN = function(x) sqrt(mean((1-x)^2)))
df_agg.mle<-aggregate(x = df.mle,by = list(df.mle$nobs, df.mle$var.var), FUN = function(x) sqrt(mean((1-x)^2)))

df_agg_mean.bayes<-df_agg.bayes[,-(1:5)]
df_agg_mean.mle<-df_agg.mle[,-(1:5)]

df_agg_mean.bayes<-apply(df_agg_mean.bayes,1,mean)
df_agg_mean.mle<-apply(df_agg_mean.mle,1,mean)

df_var.bayes<-cbind(df_agg_mean.bayes,df_agg.bayes[,(1:2)])
df_var.mle<-cbind(df_agg_mean.mle,df_agg.mle[,(1:2)])

##plotting

s.lowdens.bayes <- interp(x = log10(df_var.bayes$Group.1),y = log10(df_var.bayes$Group.2),z = log10(df_var.bayes$df_agg_mean))
#filled.contour(s.lowdens.bayes)

s.lowdens.mle <- interp(x = log10(df_var.mle$Group.1),y = log10(df_var.mle$Group.2),z = log10(df_var.mle$df_agg_mean))

#filled.contour(s.lowdens.mle)


##test coefficient of dags

##high density network

data.param<-simulateDag(node.name  = c("a","b","c","d","e","f","g","h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t"),data.dists = dist,nn = 0.8)[[1]]


data.param.var<-matrix(data = 0,nrow = 20L,ncol = 20L)


colnames(data.param.var)<-rownames(data.param.var)<-names(dist)

#plotabn(dag.m = data.param,data.dists = dist)


lseq <- function(from=0.01, to=1000, length.out=6) {
  exp(seq(log(from), log(to), length.out = length.out))
}

grid_abn<-expand.grid(nobs=c(25, 50,100,500,1000,2000,3000,5000), var.var=lseq(), rep = 1:20) 

out.bays<-list()
out.mle<-list()

seed<-24791

for(i in 1:dim(grid_abn)[1]){
  
  
  diag(data.param.var)<-rep(grid_abn[i,2],20)
  
  simulated.data<-simulateAbn(data.dists = dist,n.chains = 1,n.thin = 1,n.iter = grid_abn[i,1],
                              data.param = data.param,data.param.var = data.param.var,
                              data.param.mult = data.param.mult,seed = seed)
  
  seed<-seed+1
  
  data.bays<-fitabn(dag.m = plotabn(dag.m = data.param,data.dists = dist,plot = FALSE),data.df = simulated.data,data.dists = dist,centre = FALSE)
  data.mle<-fitabn(dag.m = plotabn(dag.m = data.param,data.dists = dist,plot = FALSE),data.df = simulated.data,data.dists = dist, method = "mle")
  out.bays<-rbind(out.bays, c(unlist(data.bays$modes)))
  out.mle<-rbind(out.mle, c(unlist(data.mle$coef)))
}

##numeric
df.bayes = apply(out.bays, 2, function(x) as.numeric(as.character(x)))
df.mle = apply(out.mle, 2, function(x) as.numeric(as.character(x)))

df.bayes <- df.bayes[, -( grep(paste0("precision") , colnames(df.bayes),perl = TRUE) ) ]
colnames(df.mle)<-colnames(df.bayes)

df.bayes <- df.bayes[, -( grep(paste0("Intercept") , colnames(df.bayes),perl = TRUE) ) ]
df.mle <- df.mle[, -( grep(paste0("Intercept") , colnames(df.mle),perl = TRUE) ) ]

df.bayes <- cbind(grid_abn,df.bayes)
df.mle <- cbind(grid_abn,df.mle)

df_agg.bayes<-aggregate(x = df.bayes,by = list(df.bayes$nobs, df.bayes$var.var), FUN = function(x) sqrt(mean((1-x)^2)))
df_agg.mle<-aggregate(x = df.mle,by = list(df.mle$nobs, df.mle$var.var), FUN = function(x) sqrt(mean((1-x)^2)))

df_agg_mean.bayes<-df_agg.bayes[,-(1:5)]
df_agg_mean.mle<-df_agg.mle[,-(1:5)]

df_agg_mean.bayes<-apply(df_agg_mean.bayes,1,mean)
df_agg_mean.mle<-apply(df_agg_mean.mle,1,mean)

df_var.bayes<-cbind(df_agg_mean.bayes,df_agg.bayes[,(1:2)])
df_var.mle<-cbind(df_agg_mean.mle,df_agg.mle[,(1:2)])


##plotting

s.highdens.bayes <- interp(x = log10(df_var.bayes$Group.1),y = log10(df_var.bayes$Group.2),z = (df_var.bayes$df_agg_mean))

s.highdens.mle <- interp(x = log10(df_var.mle$Group.1),y = log10(df_var.mle$Group.2),z = (df_var.mle$df_agg_mean))
# 
# p1 <- filled.contour(s.lowdens.bayes,plot.title = {par(cex.main=1);title(main = "Low density BN: Bayesian estimation",xlab="Number of observations [log10]",ylab="Coef of Var [log10]")},key.title = {par(cex.main=0.8);title(main="max RMSE [log10]")})
# p2 <- filled.contour(s.lowdens.mle,plot.title = {par(cex.main=1);title(main = "Low density BN: MLE estimation",xlab="Number of observations [log10]",ylab="Coef of Var [log10]")},key.title = {par(cex.main=0.8);title(main="max RMSE [log10]")})
# 
# p3 <- filled.contour(s.highdens.bayes,plot.title = {par(cex.main=1);title(main = "High density BN: Bayesian estimation",xlab="Number of observations [log10]",ylab="Coef of Var [log10]")},key.title = {par(cex.main=0.8);title(main="max RMSE [log10]")})
# p4 <- filled.contour(s.highdens.mle,plot.title = {par(cex.main=1);title(main = "High density BN: MLE estimation",xlab="Number of observations [log10]",ylab="Coef of Var [log10]")},key.title = {par(cex.main=0.8);title(main="max RMSE [log10]")})
# 
# grid.arrange(p1, p2,p3,p4, ncol=4)


pdf(file = "output/Fig_coef.pdf",width =9,height = 9)
par(mfrow=c(2,2),oma = c(0.5,0.5,0.8,0.8)-0.2,mar = c(-3,-3,1,1) + 8)
image.plot(s.lowdens.bayes,main="Low density BN: Bayesian estimation",xlab="Number of observations [log10]",ylab="Coef of Var [log10]",legend.lab = "max RMSE [log10]",legend.line = 4, col = hcl.colors(24, "YlOrRd", rev = TRUE),legend.mar = 15)
image.plot(s.lowdens.mle,main="Low density BN: MLE estimation",xlab="Number of observations [log10]", ylab="Coef of Var [log10]",legend.lab = "max RMSE [log10]",legend.line = 4,col = hcl.colors(12, "YlOrRd", rev = TRUE))

image.plot(s.highdens.bayes,main="High density BN: Bayesian estimation",xlab="Number of observations [log10]",ylab="Coef of Var [log10]",legend.lab = "max RMSE [log10]",legend.line = 4, col = hcl.colors(24, "YlOrRd", rev = TRUE))
image.plot(s.highdens.mle,main="High density BN: MLE estimation",xlab="Number of observations [log10]", ylab="Coef of Var [log10]",legend.lab = "max RMSE [log10]",legend.line = 4, col = hcl.colors(24, "YlOrRd", rev = TRUE))
dev.off()
