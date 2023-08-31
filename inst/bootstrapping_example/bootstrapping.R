##########################################################################
#  R file which peforms 10K bootstrap analysis
#
#  R creates a new seed and script file for use with JAGS in each 
#  bootstrap iteration
#  Order based searches 320 x 32 cpu
##########################################################################

## setwd("")     # possibly set directory
library("abn")
library("coda")

## Bootstrap and MCMC parameters 
nloop <- 2    # number of cores, i.e., number of parallel loops
lloop <- 100  # loop length

update1 <- 100
update2 <- 100
by2 <- 10

load('pigs_post_params.Rds')
dump(list=c("Abscess.p", "EPcat.p", "HS.p", "MS.p", "PC.p", "PDcat.p", 
            "plbinary.p", "PT.p", "Pyaemia.p", "TAIL.p"), 'pigs_post_params.R')



orig.data <- pigs.vienna[,-11]  # get the Pigs data - drop batch variable
max.par <- 3                    # parent limit for original data
start <- seq(1,    lloop*nloop, by=lloop) 
stop <- seq(lloop, lloop*nloop, by=lloop) 

## now have the boot.data in identical format to original to now repeat exact search.
ban <- matrix(0, ncol=dim(orig.data)[2], nrow=dim(orig.data)[2]) 
## the ban matrix must have names set
colnames(ban) <- rownames(ban) <- names(orig.data) 

retain <- ban   # retain nothing
## again must set names
colnames(retain) <- rownames(retain) <- names(orig.data) 

## setup distribution list for each node
mydists <- list(PC="binomial", PT="binomial", MS="binomial", HS="binomial",
  TAIL="binomial", Abscess="binomial", Pyaemia="binomial", EPcat="binomial",
  PDcat="binomial", plbinary="binomial") 


dags <- list() 
## Normally, index is passed. If not, we use value 1
if (!exists('index')) index <- 1

## MASTER LOOP
## each interation creates a bootstrap sample and finds mostprobable model
for(i in start[index]:stop[index]){ 
   # create bootstrap data
   # 1. create parameter file with unique random seed 
   init.file <- paste0("init_", i)  #file to hold jags seed
   cat(paste0("\".RNG.name\" <- \"base::Mersenne-Twister\"\n", 
      "\".RNG.seed\" <- ", i, "\n"), file=init.file) # note i is unique
   
   # 2. create script file with unique output file name
   run.file <- paste0("script_", i)   # file to hold jags seed
   out.file <- paste0("out_", i)
   
## this is needed verbatim     
cat("model in pigs_model.bug
data in  pigs_post_params.R
compile, nchains(1)
parameters in ", init.file, "\n
initialize
update ", update1, "
monitor PC, thin(10)
monitor PT, thin(10)
monitor MS, thin(10)
monitor HS, thin(10)
monitor TAIL, thin(10)
monitor Abscess, thin(10)
monitor Pyaemia, thin(10)
monitor EPcat, thin(10)
monitor PDcat, thin(10)
monitor plbinary, thin(10)
update ", update2, ", by(", by2,")
coda *, stem(\"", out.file, "\")\n", file=run.file, sep="") 
 
   # 3. run the MCMC sampler
   system(paste0("jags ", run.file))  # possibly adapt due to system

   # 4. read in mcmc data and convert to format suitable for mostprobable
   boot.data <- read.coda(paste0(out.file, "chain1.txt"),
                        paste0(out.file, "index.txt")) 
   boot.data <- as.data.frame(boot.data) 
   for(j in 1:dim(orig.data)[2]){
      if(is.factor(orig.data[,j])){
        boot.data[,j] <- as.factor(boot.data[,j]) 
        levels(boot.data[,j]) <- levels(orig.data[,j]) 
      }
   }

   # 5. run the MostProb search on the bootstrap data
   boot1.cache <- buildscorecache(data.df=boot.data, data.dists=mydists, 
                 max.parents=max.par, dag.banned=ban, dag.retained=retain) 
   dags[[i]] <- mostprobable(score.cache=boot1.cache) 
#   unlink(c(init.file,run.file,out.file,paste0(out.file,"chain1.txt"),
#                            paste0(out.file,"index.txt"))) # Tidy up

}

#unlink('pigs_post_params.R')    # more cleaning
save(dags,file=paste0("ResultBootstrap", index, ".RData")) 


# Not to forget:
if (FALSE) {
   sim.dat1 <- read.coda("out_1chain1.txt", "out_1index.txt") 
   sim.dat2 <- read.coda("out_2chain1.txt", "out_2index.txt") 
   PC <- mcmc.list(sim.dat1[,"PC"], sim.dat2[,"PC"]) 
   plot(PC)
   acf(sim.dat1[,"PC"])
}

