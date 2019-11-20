
subset<-function(Y,s,e,betas=NULL,alphas=NULL)
{
    if(is.null(betas))
    {
        betas<-4*log(seq(1,dim(Y)[2],1))
    }
    if(is.null(alphas))
    {
        alphas<-rep(2*log(dim(Y)[1]),dim(Y)[2])
    }
    Y<-Map(function(i) Y[i,],1:nrow(Y))
    return(marshall_subset(Y,s,e,betas,alphas))
}



