#' The mixture random generation for the well-known models
#' @description Random generation for the well-known mixture models with parameters \code{weigth}, \code{shape} and \code{scale}.
#' @param N Number of inputs for the mixture random generation
#' @param model Choice of one of the mixture models; \code{gompertz}, \code{log-logistics}, \code{log-normal}, and \code{weibull}.
#' @param K Number of components
#' @param param Vector of weight \eqn{\omega}, shape \eqn{\alpha}, and scale \eqn{\beta} parameters.
#' @return Outputs of random generated vector lenght of N from the given mixture model.
#' @export
#' @importFrom stats dlnorm dweibull kmeans median nlminb optim plnorm pweibull qlnorm rmultinom runif rweibull
#' @examples N<-100
#' K<-2
#' weight<-c(0.5,0.5)
#' alpha<-c(0.5,1)
#' beta<-c(1,0.5)
#' param<-c(weight,alpha,beta)
#' rmix(N, "weibull", K, param)

rmix<-function(N,model,K,param){
  if( model!="gompertz" & model!="log-logistic" & model!="log-normal" & model!="weibull"){
    stop ("No model detected!")}

  kappa<-param[1:K];alpha<-param[(K+1):(2*K)];beta<-param[(2*K+1):(3*K)];lambda<-param[(3*K+1):(4*K)]
  size<-apply(as.matrix(rmultinom(N,1,kappa)),1,sum)
  y<-rep(NA,N)
  s<-0

  if(model=="gompertz"){
    for (j in 1:length(kappa)){
      s<-s+size[j]
      if (j==1){h<-0}else{h<-1}
      y[(h*(s-size[j])+1):s]<-log(1-alpha[j]/beta[j]*log(1-runif(size[j])))/alpha[j]
    }
  }
  if(model=="log-logistic"){
    qf=function(par,u){a=par[1]; b=par[2]; b*(1/u-1)^(-1/a)}
    for (j in 1:length(kappa)){
      s<-s+size[j]
      if (j==1){h<-0}else{h<-1}
      y[(h*(s-size[j])+1):s]<-qf(c(alpha[j],beta[j]),runif(size[j]))
    }
  }
  if(model=="log-normal"){
    for (j in 1:length(kappa)){
      s<-s+size[j]
      if (j==1){h<-0}else{h<-1}
      y[(h*(s-size[j])+1):s]<-qlnorm(runif(size[j]),meanlog=alpha[j],sdlog=beta[j])
    }
  }

  if(model=="weibull"){
    for (j in 1:length(kappa)){
      s<-s+size[j]
      if (j==1){h<-0}else{h<-1}
      y[(h*(s-size[j])+1):s]<-rweibull(size[j],alpha[j],scale=beta[j])
    }
  }
  return(y)
}
