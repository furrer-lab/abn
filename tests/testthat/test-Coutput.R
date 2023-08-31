suppressPackageStartupMessages(require(testthat))
suppressPackageStartupMessages(require(abn))

data("adg", package = "abn")
dt <- adg[,c("AR","pneumS","female","livdam","eggs","wormCount","age","adg","farm")] #different order than original data
drop <- which(colnames(dt) %in% c("farm"))
# drop.clustring <- which(colnames(dt) %in% c("pneum", "epg5", "worms")) #used with cluster correction
abndata <- dt[, -drop]
abndata <- dt
dist <- list(AR = "binomial", pneumS = "binomial", female = "binomial",
             livdam = "binomial", eggs = "binomial", wormCount = "poisson",
             age = "gaussian", adg = "gaussian")

retain <- matrix(0, ncol(abndata), ncol(abndata))
colnames(retain) <- rownames(retain) <- names(abndata)

banned <- matrix(0, ncol(abndata), ncol(abndata))
colnames(banned) <- rownames(banned) <- names(abndata)

banned[3,-3] <- 1   # no arc to female


abndata[,c(1:5,9)] <- as.data.frame(lapply(abndata[,c(1:5,9)], factor))
trim.dag <- matrix( c(
 0,      0,      0,      0,    0,         0,   1,   0,
 0,      0,      0,      0,    0,         0,   0,   0,
 0,      0,      0,      0,    0,         0,   0,   0,
 0,      0,      0,      0,    1,         0,   0,   0,
 0,      0,      0,      0,    0,         0,   0,   1,
 1,      0,      0,      0,    1,         0,   1,   1,
 0,      0,      1,      0,    0,         0,   0,   0,
 0,      0,      0,      0,    0,         0,   1,   0), 8,8, byrow=T)
 colnames(trim.dag) <- rownames(trim.dag) <- names(abndata[-9])


if (F) { # the following should not be part of the regular test as they take too long.
marg.f.grouped.fitabn <- fitAbn(dag = trim.dag[1:2,1:2], data.df = as.data.frame(abndata[,c(1,2,9)]),
                                data.dists = dist[c(1,2)],
                                group.var = "farm",
                                cor.vars = c("AR", "pneumS"),
                                compute.fixed = !TRUE, n.grid = 100,
                                control=list(error.verbose=!TRUE))


# time difference between print and non and number of n.grids
system.time(fitAbn(dag = trim.dag[1:2,1:2], data.df = as.data.frame(abndata[,c(1,2,9)]),
  data.dists = dist[c(1,2)], group.var = "farm", cor.vars = c("AR", "pneumS"),
  compute.fixed = TRUE, n.grid = 100, control=list(error.verbose=TRUE)))
system.time(fitAbn(dag = trim.dag[1:2,1:2], data.df = as.data.frame(abndata[,c(1,2,9)]),
  data.dists = dist[c(1,2)], group.var = "farm", cor.vars = c("AR", "pneumS"),
  compute.fixed = TRUE, n.grid = 100, control=list(error.verbose=TRUE, trace=6)))
system.time(fitAbn(dag = trim.dag[1:2,1:2], data.df = as.data.frame(abndata[,c(1,2,9)]),
  data.dists = dist[c(1,2)], group.var = "farm", cor.vars = c("AR", "pneumS"),
  compute.fixed = TRUE, n.grid = 100, control=list(error.verbose=FALSE)))
    # verbose adds a bit, but very negligable...

system.time(fitAbn(dag = trim.dag[1:2,1:2], data.df = as.data.frame(abndata[,c(1,2,9)]),
 data.dists = dist[c(1,2)], group.var = "farm", cor.vars = c("AR", "pneumS"),
  compute.fixed = TRUE, n.grid = 10))
system.time(fitAbn(dag = trim.dag[1:2,1:2], data.df = as.data.frame(abndata[,c(1,2,9)]),
  data.dists = dist[c(1,2)], group.var = "farm", cor.vars = c("AR", "pneumS"),
  compute.fixed = TRUE, n.grid = 10000))
     # only about one second more than with 10.
}

if (F) { # the following should not be part of the regular test as they take too long.

marg.f.grouped.fitabn <- fitAbn(dag = trim.dag, data.df = as.data.frame(abndata),
                                data.dists = dist,
                                group.var = "farm",
                                cor.vars = c("AR", "pneumS", "livdam", "eggs", "wormCount", "age", "adg", "female"),
                                compute.fixed = TRUE, n.grid = 1000)

}

