#' Convert a DAG into graphviz format
#'
#' Given a matrix defining a DAG create a text file suitable for plotting with graphviz.
#'
#' @usage toGraphviz(dag,
#'                   data.df=NULL,
#'                   data.dists=NULL,
#'                   group.var=NULL,
#'                   outfile=NULL,
#'                   directed=TRUE,
#'                   verbose=FALSE)
#'
#' @param dag a matrix defining a DAG.
#' @param data.df a data frame containing the data used for learning the network.
#' @param data.dists a list with named arguments matching the names of the data frame which gives the distribution family for each variable. See \code{\link{fitAbn}} for details.
#' @param group.var only applicable for mixed models and gives the column name in \code{data.df} of the grouping variable (which must be a factor denoting group membership). See \code{fitAbn} for details.
#' @param outfile a character string giving the filename which will contain the graphviz graph.
#' @param directed logical; if TRUE, a directed acyclic graph is produced, otherwise an undirected graph.
#' @param verbose if TRUE more output is printed. If TRUE and 'outfile=NULL' the '.dot' file is printed to console.
#'
#' @details
#' Graphviz (\url{https://www.graphviz.org}) is a visualisation software developed by AT&T and freely available.
#' This function creates a text representation of the DAG, or the undirected graph, so this can be plotted using graphviz.
#' The R package, \code{Rgraphviz} (available through the Bioconductor project \url{https://www.bioconductor.org/}) interfaces R and the working installation of graphviz.
#'
#' Binary nodes will appear as squares, Gaussian as ovals and Poisson nodes as diamonds in the resulting graphviz network diagram.
#' There are many other shapes possible for nodes and numerous other visual enhancements - see online graphviz documentation.
#'
#' Bespoke refinements can be added by editing the raw outfile produced.
#' For full manual editing, particularly of the layout, or adding annotations,
#' one easy solution is to convert a postscript format graph (produced in graphviz using the -Tps switch)
#' into a vector format using a tool such as pstoedit (\url{http://www.pstoedit.net/}),
#' and then edit using a vector drawing tool like xfig.
#' This can then be resaved as postscript or pdf thus retaining full vector quality.

#'
#' @return Nothing is returned, but a file \code{outfile} written.
#'
#' @author Fraser Iain Lewis
#' @author Marta Pittavino
#'
#' @export
#'
#' @examples
#' ## On a typical linux system the following code constructs a nice
#' ## looking pdf file 'graph.pdf'.
#' \dontrun{
#'   ## Subset of a build-in dataset
#'   mydat <- ex0.dag.data[,c("b1","b2","b3","g1","b4","p2","p4")]
#'
#'   ## setup distribution list for each node
#'   mydists <- list(b1="binomial", b2="binomial", b3="binomial",
#'                   g1="gaussian", b4="binomial", p2="poisson",
#'                   p4="poisson")
#'   ## specify DAG model
#'   mydag <- matrix(c(   0,1,0,0,1,0,0, #
#'                        0,0,0,0,0,0,0, #
#'                        0,1,0,0,1,0,0, #
#'                        1,0,0,0,0,0,1, #
#'                        0,0,0,0,0,0,0, #
#'                        0,0,0,1,0,0,0, #
#'                        0,0,0,0,1,0,0  #
#'   ), byrow=TRUE, ncol=7)
#'
#'   colnames(mydag) <- rownames(mydag) <- names(mydat)
#'
#'   ## create file for processing with graphviz
#'   outfile <- paste(tempdir(), "graph.dot", sep="/")
#'   toGraphviz(dag=mydag, data.df=mydat, data.dists=mydists, outfile=outfile)
#'   ## and then process using graphviz tools e.g. on linux
#'   # system(paste( "dot -Tpdf -o graph.pdf", outfile))
#'   # system("evince graph.pdf")
#'
#'   ## Example using data with a group variable  where b1<-b2
#'   mydag <- matrix(c(0,1, 0,0), byrow=TRUE, ncol=2)
#'
#'   colnames(mydag) <- rownames(mydag) <- names(ex3.dag.data[,c(1,2)])
#'   ## specific distributions
#'   mydists <- list(b1="binomial", b2="binomial")
#'
#'   ## create file for processing with graphviz
#'   outfile <- paste0(tempdir(), "/graph.dot")
#'   toGraphviz(dag=mydag, data.df=ex3.dag.data[,c(1,2,14)], data.dists=mydists,
#'              group.var="group",
#'              outfile=outfile, directed=FALSE)
#'   ## and then process using graphviz tools e.g. on linux:
#'   # pdffile <- paste0(tempdir(), "/graph.pdf")"
#'   # system(paste("dot -Tpdf -o ", pdffile, outfile))
#'   # system(paste("evince ", pdffile, " &")   ## or some other viewer
#' }
toGraphviz <- function(dag, data.df=NULL, data.dists=NULL, group.var=NULL, outfile=NULL, directed=TRUE, verbose=FALSE){
  if (inherits(dag, "abnFit")){
    if(any(!is.null(c(data.df, data.dists, group.var)))){
      warning(paste("The argument 'dag' is of class", class(dag), ". The provided arguments 'data.df', data.dists' and 'group.var' are ignored and retrieved from the dag object."))
    }
    group.var <- dag$group.var
    data.df <- dag$abnDag$data.df
    data.dists <- dag$abnDag$data.dists
    dag <- dag$abnDag$dag
  } else if (inherits(dag, "abnLearned")){
    dag <- check.valid.dag(dag = dag$dag, data.df = data.df, group.var = group.var)

    if (any(!(data.dists %in% c("binomial", "multinomial", "poisson", "gaussian")))){
      stop('Unknown data distribution detected. Must be one of "binomial", "multinomial", "poisson", "gaussian".')
    }
  } else {
    if (any(!(data.dists %in% c("binomial", "multinomial", "poisson", "gaussian")))){
      stop('Unknown data distribution detected. Must be one of "binomial", "multinomial", "poisson", "gaussian".')
    }

    if(!is.null(group.var)){## have group variable so just need to rebuild data.df without this
      if(!(is.character(group.var) && (length(group.var)==1))){
        stop("name of group variable is not a character?!")}
      if(!length(which(group.var%in%names(data.df)==TRUE))){
        stop("name of group variable does not match any of those in data.df")}
      group.var.vals <- data.df[,group.var];## get group id data
      data.df <- data.df[,-which(names(data.df)==group.var)];## drop the group variable from original data.frame and overwrite
    }

    # check dag is in a matrix
    if(!is.matrix(dag)){
      stop("The DAG definition 'dag' must be in a matrix")}

    # check data for missing names
    if(is.null(colnames(dag)) || is.null(rownames(dag))){
      stop("'dag' must have both row and column names set")}

    # check dimension
    if(dim(dag)[1]!=dim(data.df)[2] || dim(dag)[2]!=dim(data.df)[2] ){
      stop("'dag' as dimension inconsistent with 'data.df' - if using grouped data you must supply 'group.var' argument");}

    # check binary
    for(i in 1:dim(dag)[1]){for(j in 1:dim(dag)[2]){if(dag[i,j]!=0 && dag[i,j]!=1){stop("'dag' must comprise only 1's or 0's")}}}
  }

  # check outfile
  if(!is.null(outfile)){
    if(file.exists(outfile)){
      warning("'outfile' already exists and will be overwritten.")
      outfile_path <- file.path(outfile)
    } else {
      outfile_path <- file.path(outfile)
    }
  } else if (is.null(outfile)){
    # create temporary file that will be deleted afterwards
    outfile_path <- tempfile(pattern = "model.bugs")
  } else {
    stop("'outfile' must be NULL or an existing path to a file.")
  }


  # create new file
  sink(file = outfile_path)

  ## create header part
  cat(ifelse(directed, "digraph dag {","graph dag {"),"\n\n")
  # Old version: if(directed){ cat("digraph dag {","\n\n",file=outfile,append=FALSE); }
  # Old version: else{ cat("graph dag {","\n\n",file=outfile,append=FALSE);}
  for(i in 1:length(colnames(dag))){
    if(data.dists[[i]]=="binomial"){cat(paste("\"",colnames(dag)[i],"\"[shape=box];\n",sep=""))}
    if(data.dists[[i]]=="gaussian"){cat(paste("\"",colnames(dag)[i],"\"[shape=circle];\n",sep=""))}
    if(data.dists[[i]]=="poisson"){cat(paste("\"",colnames(dag)[i],"\"[shape=oval];\n",sep=""))}
    if(data.dists[[i]]=="multinomial"){cat(paste("\"",colnames(dag)[i],"\"[shape=diamond];\n",sep=""))}
  }
  cat("\n\n\n")

  for(i in colnames(dag)){##for each variable
    children <- which(dag[,i]==1);##get row with children
    if(length(children)>=1){##if have at least one child
      child.nom <- rownames(dag)[children];
      # if(directed) {for(j in child.nom){cat("\"",i,"\"","->","\"",j,"\";","\n",sep="");}}
      #else { for(j in child.nom){cat("\"",i,"\"","--","\"",j,"\";","\n",sep="",file=outfile,append=TRUE);}
      {
        for(j in child.nom){
          cat("\"",i,"\"",ifelse(directed, "->", "--"),"\"",j,"\";","\n",sep="");
        }
      }
    }
  }
  ## footer part
  cat("\n}\n")

  # End of file creation
  sink()

  # print to stdout
  if (is.null(outfile) && verbose){
    tmp <- readLines(outfile_path)
    cat(tmp)
  } else if (verbose){
    message(paste("Printed DAG to:", normalizePath(outfile_path)))
  }
}
