#' Fits lifespan data of time units with gompertz, log-logistics, log-normal, and weibull mixture models choice of one.
#' @export
#' @details
#' Estimates parameters of the given mixture model implementing the expectation maximization (EM) algorithm.
#' General form for the cdf of a statistical mixture model is given by
#'  a distribution f is a mixture of \code{K component} distributions of
#' \eqn{f = (f_1, f_2,..f_K)} if \deqn{f(x) = \sum_{k=1}^{K}\lambda_k f_k(x)} with
#' \eqn{\lambda_k > 0}, \eqn{\sum_k \lambda_k = 1}. This equation is a stochastic model, thus
#' it allows to generate new data points; first picks a distribution of choice, with
#' probablities by weight, then generates another observation according to the chosen distribution.
#' In short represenated by,
#' Z \code{~} Mult(\eqn{\lambda_1, \lambda_2,...\lambda_k)} and
#' X|Z \code{~ f_Z}, where \code{Z} is a discrete random variable which component X is drawn from.
#'
#' The families considered for the cdf of Gompertz, Log-normal, Log-logistic, and Weibull.
#' @param lifespan numeric vector of lifespan dataset
#' @param model model name of the one of the well-known model: \code{gompertz},\code{log-logistics},\code{log-normal}, and \code{weibull}.
#' @param K number of well-known model components.
#' @param initial logical true or false
#' @param starts numeric if initial sets to true
#' @return  1.The return has three values; the first value is estimate, measures, and cluster.
#'
#' 2. The second value includes four different measurements of goodness-of-fit tests involving:
#' Akaike Information Criterion \code{(AIC)}, Bayesian Information Criterion \code{(BIC)}, Kolmogorov-Smirnov \code{(KS)}, and log-likelihood \code{(log.likelihood)} statistics.
#'
#' 3. The last value is the output of clustering vector.
#' @importFrom stats dlnorm dweibull kmeans median nlminb optim plnorm pweibull qlnorm rmultinom runif rweibull
#' @references Farewell, V. (1982). The Use of Mixture Models for the Analysis of Survival Data with Long-Term Survivors. \emph{Biometrics, 38(4), 1041-1046}. doi:10.2307/2529885
#' McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley \& Sons, Inc.
#'
#' Essam K. Al-Hussaini, Gannat R. Al-Dayian & Samia A. Adham (2000) On finite mixture of two-component gompertz lifetime model, \emph{Journal of Statistical Computation and Simulation, 67:1, 20-67}, DOI: 10.1080/00949650008812033
#' @examples
#' lifespan<-sample(1000)
#' fitmixEM(lifespan, "weibull", K = 2, initial = FALSE)
fitmixEM<-function(lifespan, model, K, initial=FALSE, starts){
  n<-length(lifespan)
  N<-60000
  cri<-5e-5
  u<-d.single<-d.mix<-cdf0<-pdf0<-clustering<-y<-alpha<-beta<-delta<-gamma<-lambda<-slice<-c()
  alpha.matrix<-beta.matrix<-lambda.matrix<-gamma.matrix<-p.matrix<-matrix(0,ncol=K,nrow=N)
  weight<-alpha<-beta2<-lambda<-delta<-Delta<-Gama<-rep(0,K)
  eu<-eu2<-tau.matrix<-zeta<-s.pdf<-matrix(0,ncol=K,nrow=n)
  clustering<-rep(0,n)
  clust<-kmeans(lifespan,K,50,1,algorithm="Hartigan-Wong")


  if (model=="gompertz"){
    if(initial==FALSE){
      for  (i in 1:K){
        slice[i]<-sum(clust$cluster==i)/n
        y<-sort(lifespan[clust$cluster==i])
        nn<-length(y)
        qp3<-y[floor(.75*nn)]
        m.lifespan<-y[floor(.5*nn)]
        a.cap<-log(log(1-.75)/log(1-.5))/(qp3-m.lifespan)
        b.cap<--a.cap*log(0.5)/(exp(a.cap*m.lifespan)-1)
        zeta1<-c()
        gompertz.log<-function(par){
          zeta1<--nn*log(par[2])-par[1]*sum(y)+par[2]/par[1]*sum(exp(par[1]*y)-1)
          zeta1[zeta1<1e-16]<-1e-16
        }
        p.cap<-suppressWarnings(optim(c(a.cap,b.cap),gompertz.log)$par)
        beta[i]<-p.cap[1]
        alpha[i]<-p.cap[2]
      }
      alpha.matrix[1,]<-alpha
      beta.matrix[1,]<-beta
      p.matrix[1,]<-slice
      for (i in 1:n){
        for (j in 1:K){
          if (clust$cluster[i]==j){zeta[i,j]<-1}
        }
      }
    }
    funcgom<-function(x,par){a=par[1]; b=par[2]; b*exp(a*(x))*exp(-(exp(a*(x))-1)*b/a)}
    if(initial==TRUE){
      p.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
      for (i in 1:K) zeta[,i]<-p.matrix[1,i]*funcgom(lifespan,c(alpha.matrix[1,i],beta.matrix[1,i]))
      d<-cbind(c(1:n),apply(zeta, 1, which.max))
      zeta<-matrix(0,nrow=n,ncol=K)
      zeta[d]<-1
    }
    r<-1
    eps<-1
    while ((eps>cri) && (r<N)){
      for (j in 1:K){
        hgom<-function(x,par){sum(-zeta[,j]*log(funcgom(x,c(par[1],par[2]))))}
        out<-suppressWarnings(nlminb(c(alpha.matrix[r,j],beta.matrix[r,j]),hgom,x=lifespan)$par)
        alpha.matrix[r+1,j]<-out[1]
        beta.matrix[r+1,j]<-out[2]
        s.pdf[,j]<-p.matrix[r,j]*funcgom(lifespan,c(alpha.matrix[r+1,j],beta.matrix[r+1,j]))
      }
      zeta<-matrix(0,ncol=K,nrow=n)
      for (i in 1:n){
        tau.matrix[i,]<-s.pdf[i,]/sum(s.pdf[i,])
        MAX<-tau.matrix[i,1]
        count<-1
        for (t in 1:K){
          if (tau.matrix[i,t]>= MAX){
            MAX<-tau.matrix[i,t]
            count<-t
          }
        }
        zeta[i,count]<-1
      }
      p.matrix[r+1,]<-apply(tau.matrix,2,sum)/n
      eps<-sum(abs(alpha.matrix[r+1,]-alpha.matrix[r,]))
      r<-r+1
    }
    weight<-round(p.matrix[r,],digits=7)
    weight[K]<-1-(sum(weight)-weight[K])
    alpha<-alpha.matrix[r,]
    beta<-beta.matrix[r,]
    pdf0<-dmix(lifespan,"gompertz",K,c(weight,alpha,beta))
    cdf0<-pmix(sort(lifespan),"gompertz",K,c(weight,alpha,beta))
  }
  if (model=="log-logistic"){
    if(initial==FALSE){
      for  (i in 1:K){
        slice[i]<-sum(clust$cluster==i)/n
        y<-sort(lifespan[clust$cluster==i])
        nn<-length(y)
        qp3<-sort(y)[floor(.75*nn)]
        alpha[i]<-log(0.75/(1-0.75))/log(qp3/median(y))
        beta[i]<-median(y)
      }
      alpha.matrix[1,]<-alpha
      beta.matrix[1,]<-beta
      p.matrix[1,]<-slice
      for (i in 1:n){
        for (j in 1:K){
          if (clust$cluster[i]==j){zeta[i,j]<-1}
        }
      }
    }
    funclog<-function(x,par){a=par[1]; b=par[2]; (a*b^(-a)*(x)^(a-1))/(((x/b)^a +1)^2)}
    if(initial==TRUE){
      p.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
      for (i in 1:K) zeta[,i]<-p.matrix[1,i]*funclog(lifespan,c(alpha.matrix[1,i],beta.matrix[1,i]))
      d<-cbind(c(1:n),apply(zeta, 1, which.max))
      zeta<-matrix(0,nrow=n,ncol=K)
      zeta[d]<-1
    }
    r<-1
    eps<-1
    while ((eps>cri) && (r<N)){
      for (j in 1:K){
        hlog<-function(x,par){sum(-zeta[,j]*log(funclog(x,c(par[1],par[2]))))}
        out<-suppressWarnings(nlminb(c(alpha.matrix[r,j],beta.matrix[r,j]),hlog,x=lifespan)$par)
        alpha.matrix[r+1,j]<-out[1]
        beta.matrix[r+1,j]<-out[2]
        s.pdf[,j]<-p.matrix[r,j]*funclog(lifespan,c(alpha.matrix[r+1,j],beta.matrix[r+1,j]))
      }
      zeta<-matrix(0,ncol=K,nrow=n)
      for (i in 1:n){
        tau.matrix[i,]<-s.pdf[i,]/sum(s.pdf[i,])
        MAX<-tau.matrix[i,1]
        count<-1
        for (t in 1:K){
          if (tau.matrix[i,t]>= MAX){
            MAX<-tau.matrix[i,t]
            count<-t
          }
        }
        zeta[i,count]<-1
      }
      p.matrix[r+1,]<-apply(tau.matrix,2,sum)/n
      eps<-sum(abs(alpha.matrix[r+1,]-alpha.matrix[r,]))
      r<-r+1
    }
    weight<-round(p.matrix[r,],digits=7)
    weight[K]<-1-(sum(weight)-weight[K])
    alpha<-alpha.matrix[r,]
    beta<-beta.matrix[r,]
    pdf0<-dmix(lifespan,"log-logistic",K,c(weight,alpha,beta))
    cdf0<-pmix(sort(lifespan),"log-logistic",K,c(weight,alpha,beta))
  }

  if (model=="log-normal"){
    if(initial==FALSE){
      for  (i in 1:K){
        slice[i]<-sum(clust$cluster==i)/n
        y<-sort(lifespan[clust$cluster==i])
        nn<-length(y)
        median.lifespan<-median(y)
        alpha[i]<-log(median.lifespan)
        beta[i]<-sqrt(2*abs(log(mean(y)/median.lifespan)))
      }
      alpha.matrix[1,]<-alpha
      beta.matrix[1,]<-beta
      p.matrix[1,]<-slice
      for (i in 1:n){
        for (j in 1:K){
          if (clust$cluster[i]==j){zeta[i,j]<-1}
        }
      }
    }
    if(initial==TRUE){
      p.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
      for (i in 1:K) zeta[,i]<-p.matrix[1,i]*dlnorm(lifespan,meanlog=alpha.matrix[1,i],sdlog=beta.matrix[1,i])
      d<-cbind(c(1:n),apply(zeta, 1, which.max))
      zeta<-matrix(0,nrow=n,ncol=K)
      zeta[d]<-1
    }
    r<-1
    eps<-1
    while ((eps>cri) && (r<N)){
      for (j in 1:K){
        hlogn<-function(x,par){sum(-zeta[,j]*dlnorm(x,meanlog=par[1],sdlog=par[2],log=TRUE))}
        out<-suppressWarnings(nlminb(c(alpha.matrix[r,j],beta.matrix[r,j]),hlogn,x=lifespan)$par)
        alpha.matrix[r+1,j]<-out[1]
        beta.matrix[r+1,j]<-out[2]
        s.pdf[,j]<-p.matrix[r,j]*dlnorm(lifespan,meanlog=alpha.matrix[r+1,j],sdlog=beta.matrix[r+1,j])
      }
      zeta<-matrix(0,ncol=K,nrow=n)
      for (i in 1:n){
        tau.matrix[i,]<-s.pdf[i,]/sum(s.pdf[i,])
        MAX<-tau.matrix[i,1]
        count<-1
        for (t in 1:K){
          if (tau.matrix[i,t]>= MAX){
            MAX<-tau.matrix[i,t]
            count<-t
          }
        }
        zeta[i,count]<-1
      }
      p.matrix[r+1,]<-apply(tau.matrix,2,sum)/n
      eps<-sum(abs(alpha.matrix[r+1,]-alpha.matrix[r,]))
      r<-r+1
    }
    weight<-round(p.matrix[r,],digits=7)
    weight[K]<-1-(sum(weight)-weight[K])
    alpha<-alpha.matrix[r,]
    beta<-beta.matrix[r,]
    pdf0<-dmix(lifespan,"log-normal",K,c(weight,alpha,beta))
    cdf0<-pmix(sort(lifespan),"log-normal",K,c(weight,alpha,beta))
  }


  if (model=="weibull"){
    if(initial==FALSE){
      for  (i in 1:K){
        slice[i]<-sum(clust$cluster==i)/n
        y<-sort(lifespan[clust$cluster==i])
        nn<-length(y)
        k<-seq(1,nn)
        X<-cbind(rep(1,nn),log(-log(1-k/(nn+1)))+(k*(nn-k+1))/((nn+1)^2*(nn+2))*(log(1-k/(nn+1))+1)/
                   ((1-k/(nn+1))*log(1-k/(nn+1)))^2)
        V<-matrix(0,nn,nn)
        for(ii in 1:nn){
          for(jj in ii:nn){
            V[ii,jj]<-ii/(nn+1-ii)*(log((nn+1-ii)/(nn+1))*log((nn+1-jj)/(nn+1)))^(-1)
          }
        }
        U<-solve(t(X)%*%solve(V)%*%X)%*%t(X)%*%solve(V)%*%log(y)
        beta[i]<-exp(U[1])
        alpha[i]<-1/U[2]
      }
      alpha.matrix[1,]<-alpha
      beta.matrix[1,]<-beta
      p.matrix[1,]<-slice
      for (i in 1:n){
        for (j in 1:K){
          if (clust$cluster[i]==j){zeta[i,j]<-1}
        }
      }
    }
    if(initial==TRUE){
      p.matrix[1,]<-starts[1:K]
      alpha.matrix[1,]<-starts[(K+1):(2*K)]
      beta.matrix[1,]<-starts[(2*K+1):(3*K)]
      for (i in 1:K) zeta[,i]<-p.matrix[1,i]*dweibull(lifespan,shape=alpha.matrix[1,i],scale=beta.matrix[1,i])
      d<-cbind(c(1:n),apply(zeta, 1, which.max))
      zeta<-matrix(0,nrow=n,ncol=K)
      zeta[d]<-1
    }
    r<-1
    eps<-1
    while ((eps>cri) && (r<N)){
      for (j in 1:K){
        hweib<-function(x,par){sum(-zeta[,j]*dweibull(x,shape=par[1],scale=par[2],log=TRUE))}
        out<-suppressWarnings(nlminb(c(alpha.matrix[r,j],beta.matrix[r,j]),hweib,x=lifespan)$par)
        alpha.matrix[r+1,j]<-out[1]
        beta.matrix[r+1,j]<-out[2]
        s.pdf[,j]<-p.matrix[r,j]*dweibull(lifespan,shape=alpha.matrix[r+1,j],scale=beta.matrix[r+1,j])
      }
      zeta<-matrix(0,ncol=K,nrow=n)
      for (i in 1:n){
        tau.matrix[i,]<-s.pdf[i,]/sum(s.pdf[i,])
        MAX<-tau.matrix[i,1]
        count<-1
        for (t in 1:K){
          if (tau.matrix[i,t]>= MAX){
            MAX<-tau.matrix[i,t]
            count<-t
          }
        }
        zeta[i,count]<-1
      }
      p.matrix[r+1,]<-apply(tau.matrix,2,sum)/n
      eps<-sum(abs(alpha.matrix[r+1,]-alpha.matrix[r,]))
      r<-r+1
    }
    weight<-round(p.matrix[r,],digits=7)
    weight[K]<-1-(sum(weight)-weight[K])
    alpha<-alpha.matrix[r,]
    beta<-beta.matrix[r,]
    pdf0<-dmix(lifespan,"weibull",K,c(weight,alpha,beta))
    cdf0<-pmix(sort(lifespan),"weibull",K,c(weight,alpha,beta))
  }

  for(i in 1:n)clustering[i]<-which(zeta[i,]==1)[1]
  for(i in 1:n){
    u[i]<-ifelse(cdf0[i]==1,0.99999999,cdf0[i])
  }
  n.p<-K*3-1
  log.likelihood<-suppressWarnings(sum(log(pdf0)))
  L<-seq(1,n)
  ks.stat<-suppressWarnings(max(L/n-cdf0,cdf0-(L-1)/n))
  AIC<--2*log.likelihood + 2*n.p
  BIC<--2*log.likelihood + n.p*log(n)
  out2<-cbind(AIC, BIC, ks.stat, log.likelihood)
  colnames(out2)<-c("AIC","BIC","KS","log.likelihood")
  if (model=="skew-normal"){
    out1<-cbind(weight,alpha,beta,lambda)
    colnames(out1)<-c("weight","alpha","beta","lambda")
  }else{
    out1<-cbind(weight,alpha,beta)
    colnames(out1)<-c("weight","alpha","beta")
  }
  out3<-clustering
  list("estimate"=out1,"measures"=out2,"cluster"=out3)
}
