
.wbs.subset.class<-setClass("wbs.subset.class",representation(data="array",betas="vector",alphas="vector",zeta="numeric",M="integer",seed="integer",
                                                  cpts="numeric",affected="list"))

wbs.subset.class<-function(data,betas,alphas,zeta,M,seed,cpts,affected,...)
{
    .wbs.subset.class(data=data,betas=betas,alphas=alphas,zeta=zeta,M=M,seed=seed,cpts=cpts,affected=affected,...)
}



#' Wild Binary Segmentation using SUBSET.
#'
#' @name wbs.subset
#'
#' @description Applies Wild Binary Segmentation using SUBSET to detect changepoints in multivariate data 
#' 
#' @param data An n by m matrix, array or data.frame containg n observations of m variates
#' @param M An integer > 0 indicating the number of random intervals of [0,n] to generate for use with SUBSET
#' @param zeta The threshold value.
#' @param alphas The value of alphas used in SUBET. If not specified alphas is set to \code{rep(2*log(m),m)}. When specified, \code{alphas} should be a vector
#' with length equal to the number of variates.
#' @param betas The value of betas used in SUBET. If not specified betas is set to \code{4*log(1:n)}. When specified, \code{betas} should be a vector of length
#' equal to the number observations or a function. If \code{betas} is a function it should take a single argument (the length of the subinterval) and return a
#' numeric value.
#' 
#' @return An S4 class of type wbs.subset.class
#'
#' 
#' @docType methods
#'
#' @rdname wbs.subset-methods
#'
#' @examples
#' library(subset)
#' x0<-c(rnorm(100,0,1),rnorm(200,10,1))
#' x1<-c(rnorm(200,0,1),rnorm(100,10,1))
#' x2<-c(rnorm(200,0,1),rnorm(100,10,1))
#' X<-matrix(c(x0,x1,x2),300,3)
#' X<-robustscale(X)
#' res<-wbs.subset(X,1000,5)
#' cpt.locations(res)
#'
wbs.subset<-function(data,M,zeta,alphas=NULL,betas=NULL)
{
    # check data
    if(is.data.frame(data))
    {
        data<-data.matrix(data)
    }
    if(!is.array(data))
    {
        stop("data must be a matrix, array, or data frame")
    }
    # check M
    if(!is.numeric(M))
    {
        stop("M must be numeric")
    }
    if(M < 1)
    {
        stop("M must be a integer > 0")
    }
    if(M < 0)
    {
        stop("M must be a positive integer")
    }
    # get seed
    seed<-as.integer(.Machine$integer.max*runif(1))
    # check betas
    if(is.function(betas))
    {
        betas<-unlist(Map(betas,1:dim(data)[1]))
    }
    if(is.null(betas))
    {
        betas<-4*log(seq(1,dim(data)[1],1))
    }
    if(!is.numeric(betas))
    {
        stop("betas must be a numeric vector or a function")
    }
    if(length(betas) != dim(data)[1])
    {
        stop("the length of betas should equal the number of observations")
    }
    # check alphas
    if(is.null(alphas))
    {
        alphas<-rep(2*log(dim(data)[2]),dim(data)[1])
    }
    if(!is.numeric(alphas))
    {
        stop("alphas must be a numeric vector or a function")
    }
    if(length(alphas) != dim(data)[1])
    {
        stop("the length of alphas should equal the number of variates")
    }
    list.data<-Map(function(i) data[,i],1:ncol(data))
    result<-marshall_wbs_subset(list.data,M,zeta,betas,alphas,seed)
    return(wbs.subset.class(data,betas,alphas,zeta,as.integer(M),as.integer(seed),sort(result$breakpoints),result$affected))
}



#' Process changepoint locations
#'
#' @name cpt.locations
#'
#' @description Extracts the changepoint locations and affected variates from the results of \code{\link{wbs.subset}}.
#' 
#' @param object An S4 object produced by \code{\link{wbs.subset}}.
#' 
#' @return A data.frame containing the affected variates and the changepoint locations.
#'
#' 
#' @docType methods
#'
#' @rdname cpt.locations-methods
#'
if(!isGeneric("cpt.locations")) {setGeneric("cpt.locations",function(object,...) {standardGeneric("cpt.locations")})}
setMethod("cpt.locations",signature=list("wbs.subset.class"),
          function(object)
          {
              affected.variates<-sort(unique(unlist(object@affected)))
              selections<-Map(function(variate) Map(function(affected) variate %in% affected,object@affected),affected.variates)
              result<-
                  Reduce(rbind,
                         Map(function(variate,selection)
                         {
                             cpts<-object@cpts[unlist(selection)]
                             return(data.frame("variate"=rep(variate,"location"=length(cpts)),cpts))
                         },
                         affected.variates,selections),
                         data.frame()
                         )
              result$variate<-result$variate+1
              result$cpts<-result$cpts+1
              return(result)
          })
          

#' Summary of results produced by wbs.subset.
#'
#' @name summary
#'
#' @description Summary methods for S4 object returned by \code{wbs.subset}. The output lists the paramaters provided by the user (or the defaults used), the dimensions of the data,
#' the location of each changepoint, and the affected variates.  
#'
#' @docType methods
#'
#' @rdname summary-methods
#'
#' @aliases summary,wbs.subset.class-method
#' 
#' @export
setMethod("summary",signature=list("wbs.subset.class"),
          function(object)
          {
              cat("wbs.subset summary -- ",'\n',sep="")
              cat("number of observations : ",nrow(object@data),'\n',sep="")
              cat("number of variates : ",ncol(object@data),'\n',sep="")
              cat("number of subintervals : ",object@M,'\n',sep="")
              cat("zeta : ",object@zeta,'\n',sep="")
              cpt.info<-cpt.locations(object)
              
              if(dim(cpt.info)[1] != 0)
              {
                  cat("number of change points : ",length(unique(cpt.info$cpts)),'\n',sep="")
                  cat("number of affected variates : ",length(unique(cpt.info$variate)),'\n',sep="")
                  print(cpt.info)
              }
          })

#' Displays S4 objects produced by capa methods.
#'
#' @name show
#'
#' @description Displays S4 object produced by \code{\link{wbs.subset}}. The information produced is the same as that provided by the summary method.
#' The method is used by the S4 system for automatic printing.
#'
#' @docType methods
#'
#' @param object An S4 class produced by \code{\link{wbs.subset}}.
#' 
#' @rdname show-methods
#'
#' @aliases show,capa.class-method
setMethod("show",signature=list("wbs.subset.class"),
          function(object)
          {
              summary(object)
          })






