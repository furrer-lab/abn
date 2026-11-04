rm(list = ls(all.names = TRUE))

bits_to_key <- function(x) {
  as.character(sum(x * 2^seq.int(from = length(x) - 1, to = 0)))
}

string_to_key <- function(x) {
    paste(x, collapse = "")
}

fun.return <- function(x, n) {
  # x: vector or single integer
  # n: number of variables/columns in data.df
  # returns: vector of size n-1 with 1 for all x and zeros otherwise
    v <- rep(0, n - 1)
    v[x] <- 1
    return(v)
}


set.seed(1)
nvars <- 54
max.parents <- 3

####
# From build_score_cache_mle.R
tmp.list <- vector('list', max.parents)
for(i in 1:(max.parents)) {
    tmp.list[[i]] <- t(combn(x = (nvars - 1), m = i, FUN = fun.return, n = nvars, simplify = TRUE))
}
combs <- do.call(rbind, tmp.list)
storage.mode(combs) <- "integer"
cache_nrow <- nrow(combs) + 1L
node.defn <- matrix(data = 0L,
		    nrow = nvars * cache_nrow,
		    ncol = nvars)
for(j in 1:nvars) {
    n_filled <- (j - 1L) * cache_nrow
    node.defn[(n_filled + 2L):(n_filled + cache_nrow), -j] <- combs
}
####

ok_node.defn <- node.defn[duplicated(node.defn)==FALSE,]
# Specify if keys are defined using:
# 0: bits_to_key function
# 1: inline definition of the bits_to_key
# 2: string_to_key function
# 3: inline definition of the string_to_key
key_mode <- "1"
switch (key_mode,
             "0" = {keys <- apply(ok_node.defn, 1, bits_to_key )},
	     "1" = {keys <- apply(ok_node.defn, 1, function(r) as.character(sum(r * 2^seq.int(from = length(r) - 1, to = 0))))},
	     "2" = {keys <- apply(ok_node.defn, 1, string_to_key )},
	     "3" = {keys <- apply(ok_node.defn, 1, function(x)  paste(x, collapse = "") )}
	     )
df <- data.frame(
  key = keys,
  input = apply(ok_node.defn, 1, paste, collapse = "")
)


# e.g. df[apply(df, 2, function(r) r==864691128455135232)==TRUE, ]
nvars
nrow(node.defn)
nrow(ok_node.defn)
length(which(duplicated(df$key)==TRUE))

### Results for key_mode = 0
# nvars = 60, nrow(node.defn)=2056800, nrow(ok_node.defn)=36051, length(which(duplicated(df$key)==TRUE)) = 1540
# nvars = 50, nrow(node.defn)=982500, nrow(ok_node.defn)=20876, length(which(duplicated(df$key)==TRUE))= 0 (collision) 
# nvars = 54, nrow(node.defn)=1342332, nrow(ok_node.defn)=26290, length(which(duplicated(df$key)==TRUE))=53 (collision)
