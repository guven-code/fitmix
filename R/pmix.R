#' The mixture cumulative distribution
#' @description Computing cumulative distribution function for the well-known mixture models.
#' @param lifespan Vector of samples
#' @param model choice of one of the mixture models; \code{gompertz}, \code{log-logistics}, \code{log-normal}, and \code{weibull}.
#' @param K number of components
#' @param param Vector of weight \eqn{\omega}, shape \eqn{\alpha}, and scale \eqn{\beta} parameters.
#' @return A vector of the same length as lifespan data, given the cdf of the one of the mixture models computed at \code{lifespan}.
#' @export
#' @importFrom stats dlnorm dweibull kmeans median nlminb optim plnorm pweibull qlnorm rmultinom runif rweibull
#' @examples lifespan<-seq(0,30,0.2)
#' K<-2
#' weight<-c(0.5,0.5)
#' alpha<-c(0.5,1)
#' beta<-c(1,0.5)
#' param<-c(weight,alpha,beta)
#' pmix(lifespan, "log-logistic", K, param)

pmix<-function(lifespan,model,K,param){
  if( model!="gompertz" & model!="log-logistic" & model!="log-normal" & model!="weibull"){
    stop ("No model detected!")}
  y<-rep(NA,length(lifespan));

  if(model=="gompertz"){
    cum=function(x,par){a=par[1]; b=par[2]; 1-exp(-(exp(a*x)-1)*b/a)}
  }
  if(model=="log-logistic"){
    cum=function(x,par){a=par[1]; b=par[2]; 1/((x/b)^(-a)+1)}
  }
  if(model=="log-normal"){
    cum=function(x,par){a=par[1]; b=par[2]; plnorm(x,a,sdlog=b)}
  }

  if(model=="weibull"){
    cum=function(x,par){a=par[1]; b=par[2]; pweibull(x,a,scale=b)}
  }
  kappa<-param[1:K];alpha<-param[(K+1):(2*K)];beta<-param[(2*K+1):(3*K)];lambda<-param[(3*K+1):(4*K)]
  for (i in 1:length(lifespan)){
    cdf<-0
    for (j in 1:K){
      cdf<-cdf+kappa[j]*cum(lifespan[i],c(alpha[j],beta[j],lambda[j]))
    }
    y[i]<-cdf
  }
  return(y)
}
