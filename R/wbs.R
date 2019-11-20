
wbs<-function(data,M,zeta,seed=NULL)
{
    if(is.null(seed))
    {
        seed<-as.integer(.Machine$integer.max*runif(1))
    }
    return(marshall_wbs(data,M,zeta,seed))
}
