
subset<-function(Y,s,e,betas=NULL,alphas=NULL,thresholds=NULL)
{
    if(is.null(betas))
    {
        betas<-4*log(seq(1,dim(Y)[2],1))
    }
    if(is.null(alphas))
    {
        alphas<-rep(2*log(dim(Y)[1]),dim(Y)[2])
    }
    if(is.null(thresholds))
    {
        thresholds<-dim(Y)[2] + sqrt(2.0*dim(Y)[2]*betas)
    }
    Y<-Map(function(i) Y[i,],1:nrow(Y))
    return(marshall_subset(Y,s,e,betas,alphas,thresholds))
}



