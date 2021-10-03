.loss <- function(X,y,beta,phi,omega,kappa){
  ## Define loss objective function for MINB model.

  J <- length(kappa)
  n <- dim(X)[1]
  p <- dim(X)[2]
  mu <- exp(drop(X%*%beta))
  tmp1 <- lgamma( y + 1.0/phi)-lgamma(y+1)-lgamma(1.0/phi)
  tmp1 <- exp(tmp1)

  tmp2 <- (phi*mu/(1.0 + phi*mu))^y

  tmp3 <- (1.0/(1.0 + phi*mu))^(1.0/phi)
  tmp <- tmp1 *tmp2 * tmp3
  tmp[which(tmp == 0)] = 1e-300
  if(is.null(omega)){
    log_pry = log(tmp)
  }else{
    delta <- .Delta(y, kappa)

    location <- apply( delta, 1, function(x) which(x==1))
    if( sum(omega) >= 1 ){
      substract_value <- min(omega) / 2
      omega <- omega - substract_value
    }
    omega_complete <- c(omega, 1.0 - sum(omega))
    omega_prime <- omega_complete[location]*(location!=(J+1))

    log_pry <- log(omega_prime + omega_complete[J+1] * tmp)
  }


  #Compute loss function
  loss <- -1.0 * sum(log_pry)
  return(loss)
}

.bic <- function(X,y,beta,phi,kappa,omega,thresh){

  ## Define function for calculating the value of BIC criterion.

  n <- dim(X)[1]
  p <- dim(X)[2]
  df <- p + sum(omega >= thresh)
  loss1 <- .loss(X,y,beta,phi,omega,kappa)
  bic <- 2*loss1+log(n)*df
  return(bic)
}

.Delta <- function(y,kappa){

  ## Define indicator matrix for MINB model.

  J <- length(kappa)
  n <- length(y)
  delta <- matrix(0, nrow = n, ncol = J+1)
  for(i in 1:length(y))
    {
    delta[i,] <- .delta_ij(y[i],kappa)
  }
  return(delta)
}

.delta_ij <- function(y_i,kappa){

  J <- length(kappa)
  tmp <- which(y_i==kappa)

  delta_ij <- rep(0,J+1)
  if(length(tmp)){delta_ij[tmp]=1}
  else{ delta_ij[J+1] <- 1}
  return(delta_ij)
}

.gammaij <- function(X,y,beta,phi,omega,kappa){

  ## Define function for calculating the hidden variable matrix.

  J <- length(kappa)
  n <- nrow(X)
  delta1 <- .Delta(y,kappa)
  delta1[,J+1] <- rep(1,n)
  prynb <- .pry_nb(X,y,beta,phi)

  tmp1 <- (1-sum(omega))*prynb
  tmp2 <- matrix(rep(omega,n),nrow=n,byrow=TRUE)
  tmp2 <- cbind(tmp2,tmp1)
  tmp3 <- delta1*tmp2

  SUM <- apply(tmp3,1,sum)
  SUM <- matrix(rep(SUM,J+1),nrow=n,byrow =FALSE)
  gammaij <- tmp3/SUM

  return(gammaij)
}

.grbeta<-function(X,y,beta,phi,omega,kappa,gamma_ij){
  # Define the gradient on beta for M-step in EM iterative algorithm

  J <- length(kappa)
  n <- dim(X)[1]
  p <- dim(X)[2]
  mu <- exp(drop(X%*%beta))
  dmu <- mu*X
  delta <- .Delta(y, kappa)
  location <- apply(delta, 1, function(x) which(x==1))
  omega_complete <- c(omega, 1.0 - sum(omega))
  # tmp1 <- rep(0, n)
  # for(i in 1:n){
  #   if(y[i]!= 0){
  #     tmp1[i] <- sum(log(0:(y[i]-1)+1.0/phi))
  #   }
  #   tmp1[i] <- exp(tmp1[i])/gamma(y[i] + 1)
  # }
  tmp1 <- lgamma( y + 1.0/phi)-lgamma(y+1)-lgamma(1.0/phi)
  tmp1 <- exp(tmp1)
  tmp2 <-(phi*mu/(1.0 + phi*mu))^y
  tmp3 <-(1.0/(1.0 + phi*mu))^(1.0/phi)
  tmp <- tmp1 *tmp2 * tmp3
  tmp[which(tmp == 0)] = 1e-300
  pry <- omega_complete[J+1] * tmp
  dtmp2_mu <- y * ((phi*mu/(1.0 + phi*mu))^(y-1)) *
    (phi/((1.0 + phi*mu)^2))

  dtmp3_mu <- (-1.0) * ((1.0/(1.0 + phi*mu))^(1.0/phi-1)) * (1.0/((1.0 + phi*mu)^2))

  dpry_mu <- omega_complete[J+1] * tmp1 *
    (dtmp2_mu * tmp3 + tmp2 * dtmp3_mu)
  dlogpry_mu <- dpry_mu/pry
  gamma_iJ<-gamma_ij[,J+1]
  gradient<-gamma_iJ*dlogpry_mu
  A<-gradient*dmu
  dobj_beta<-apply(A, 2, sum)
  dobj_beta<- -1.0*dobj_beta
  return(dobj_beta)
}

.grphi <- function(X,y,beta,phi,omega,kappa,gamma_ij){
  # Define the gradient on phi for M-step in EM iterative algorithm
  # Args:
  #   X: design matrix. --Matrix
  #   y: response vector --Vector
  #   beta: values of NB process --Vector
  #   phi:   dispersion parameter of NB process -- Scalar
  #   omega:  values of inflated points --Vector
  #   kappa:  inflated values -- Vector
  #   gamma_ij : the conditional expectation of gamma at current step --Matrix
  #
  # Returns:
  #   Gradient of parameter phi --Vector
  J <- length(kappa)
  n <- dim(X)[1]
  p <- dim(X)[2]
  mu <- exp(drop(X%*%beta))
  omega_complete <- c(omega, 1.0 - sum(omega))
  delta <- .Delta(y, kappa)
  location <- apply(delta, 1, function(x) which(x==1))
#  tmp1 <- rep(0, n)
#  for(i in 1:n){
#    if(y[i]!= 0){
#      tmp1[i] <- sum(log(0:(y[i]-1)+1.0/phi))
#    }
#    tmp1[i] <- exp(tmp1[i])/gamma(y[i] + 1)
# }
  tmp1 <- lgamma( y + 1.0/phi)-lgamma(y+1)-lgamma(1.0/phi)
  tmp1 <- exp(tmp1)

  tmp2 <-(phi*mu/(1.0 + phi*mu))^y

  tmp3 <-(1.0/(1.0 + phi*mu))^(1.0/phi)
  tmp <- tmp1 *tmp2 * tmp3
  tmp[which(tmp == 0)] = 1e-300
  pry<-omega_complete[J+1] * tmp
  dtmp1_phi <- tmp1*(digamma(y+1.0/phi)-digamma(1.0/phi))*(-1.0/(phi^2))
  dtmp1_phi[which(dtmp1_phi == -Inf)] = -1e+20
  dtmp2_phi <-  y * ((phi*mu/(1.0 + phi*mu))^(y-1)) * (mu/((1.0 + phi*mu)^2))
  dtmp3_phi <- log(1.0 + phi * mu)-(phi * mu)/(1.0 + phi * mu)
  dtmp3_phi[which(dtmp3_phi == -Inf)] = -1e+20
  dtmp3_phi <- dtmp3_phi * tmp3 / (phi * phi)
  dpry_phi <- omega_complete[J+1] *(dtmp1_phi * tmp2 * tmp3 + tmp1 * dtmp2_phi * tmp3 +
                                      tmp1 * tmp2 * dtmp3_phi)
  dlogpry_phi <- dpry_phi / pry
  dobj_phi <- sum(gamma_ij[,J+1]*dlogpry_phi)
  dobj_phi<- -1.0*dobj_phi
  return(dobj_phi)
}

.init <- function(X, y, nlambda,lambda.min.ratio,lambda=NULL,control){
  #  compute initialized values for the algorithm

  p <- dim(X)[2]
  n <- dim(X)[1]

  kappa0 <- sort(unique(y))
  J <- length(kappa0)
  omega0 <- sapply(kappa0, function(x) sum(y == x))
  omega0 <- omega0 / length(y)
  substract_value <- min(omega0) / 2
  omega0 <- omega0 - substract_value
  weight0<-1.0/abs(omega0)
  names(omega0) <- paste0("omega", kappa0)
  names(weight0) <- paste0("omega", kappa0)

  data <- as.data.frame(cbind(X, y))
  colnames(data) <- c(paste("beta", 0:(p-1), sep=""), "y")
  omega_start <- omega0

  tmp <- glm.nb(y~., data = data[,-1],control=glm.control(maxit=100))
  beta_start <- tmp$coefficients
  phi_start <- 1.0 / tmp$theta

  if(is.null(lambda)){
    if (lambda.min.ratio >= 1)
      stop("lambda.min.ratio should be less than 1")

  for (i in seq(0.1*n,10*n,ifelse(n>=1000,50,10))){
    lambdamax = i
    omega.lambda.max = .iterate(X,y,beta_start,phi_start,omega_start,kappa0,lambdamax,weight0,control)$omega
    if (length(omega.lambda.max)==0){
      break
    }

  }
  lambdamin <- lambdamax * lambda.min.ratio
  loghi <- log(lambdamax)
  loglo <- log(lambdamin)
  logrange <- loghi - loglo
  interval <- -logrange / (nlambda - 1)
  lambda <- exp(seq(loghi, loglo, interval))}
  else{
    if (any(lambda < 0))
      stop("lambdas should be non-negative")
    if (nlambda != length(lambda))
      stop("nlambda should be equal with the number of the lambda")
  }

  list(beta=beta_start,phi=phi_start,omega=omega0,lambda=lambda,weight=weight0,kappa=kappa0)
}

.iterate <- function(X,y,beta,phi,omega,kappa,lambda,weight,control){
  # Define the EM iterative algorithm for estimating parameters in stage one

  alpha = 100
  ui = rbind(c(1),c(-1))
  ci = c(1e-20,-100)
  i=1
  while(alpha > control$toler && i< control$maxit){
    tmpbeta = beta
    tmpphi = phi
    tmpomega = omega
    gamma_ij = .gammaij(X,y,beta,phi,omega,kappa)
    beta <- optim(beta, .obj, .grbeta, method=control$method,
                X=X,y=y,phi=phi,omega=omega,kappa=kappa,gamma_ij=gamma_ij)$par
    phi <- constrOptim(phi, .obj, .grphi,ui,ci,method=control$method,
                     X=X,y=y,beta=beta,omega=omega,kappa=kappa,gamma_ij=gamma_ij)$par
    tmp <- .omega.iterate(X,y,beta,phi,omega,kappa,lambda,weight,gamma_ij)
    omega <- tmp$omega

    if(!length(omega)){
      omega = kappa = weight = NULL
      fit.nb <- glm.nb(y~., data = as.data.frame(cbind(X,y))[,-1],control=glm.control(maxit=100))
      beta <- fit.nb$coefficients
      phi <- 1.0 / fit.nb$theta
      break
    }
    if(sum(omega)>=1){
      substract_value <- min(omega) / 2
      omega <- omega - substract_value
    }
    kappa<-tmp$kappa
    weight<-tmp$weight
    name <-names(omega)
    tmpomega2<- rep(0,length(tmpomega))
    names(tmpomega2)<- names(tmpomega)
    tmpomega2[name]<- omega
    norm1 <- sum(abs(beta-tmpbeta))/sqrt(sum(tmpbeta^2))
    norm2 <- abs(phi-tmpphi)/abs(phi)
    norm3 <- sum(abs(tmpomega2-tmpomega))
    alpha <- max(norm1, norm2, norm3)

    i<-i+1
    if(i == control$maxit)
      warning("fitting results did not converge -- increase 'maxit'? ")
  }
  list(omega=omega,kappa=kappa,weight=weight,beta=round(beta,6),phi=phi)
}

.iterate2 <- function(X,y,beta,phi,omega,kappa,control){
  # Define the EM iterative algorithm for estimating parameters without penalization


  alpha = 100
  ui = rbind(c(1),c(-1))
  ci = c(1e-20,-100)
  i=1
  while(alpha > control$toler && i< control$maxit){
    tmpbeta=beta
    tmpphi=phi
    tmpomega=omega
    gamma_ij=.gammaij(X,y,beta,phi,omega,kappa)
    omega<-.omega_em(X,y,beta,phi,omega,kappa,gamma_ij)
    beta <- optim(beta, .obj, .grbeta, method=control$method,
                  X=X,y=y,phi=phi,omega=omega,kappa=kappa,gamma_ij=gamma_ij)$par
    phi <- constrOptim(phi, .obj, .grphi,ui,ci,method=control$method,
                       X=X,y=y,beta=beta,omega=omega,kappa=kappa,gamma_ij=gamma_ij)$par
    norm1 <- sum(abs(beta-tmpbeta))/sqrt(sum(tmpbeta^2))
    norm2 <- abs(phi-tmpphi)/abs(phi)
    norm3 <- sum(abs(omega-tmpomega))
    alpha <- max(norm1, norm2, norm3)

    i<-i+1
    if(i == control$maxit)
      warning("fitting results did not converge -- increase 'maxit'?")
  }
  list(omega=omega,kappa=kappa,beta=round(beta,6),phi=phi)
}

.obj<-function(X,y,beta,phi,omega,kappa,gamma_ij){

  ## Define complete-data objective function for EM algorithm.


  J <- length(kappa)
  prynb <- .pry_nb(X,y,beta,phi)
  tmp1 <- log((1-sum(omega))*prynb)
  obj1 <- sum(gamma_ij[,J+1]*tmp1)
  obj1 <- -1.0*obj1
  return(obj1)
}

.omega_em <- function(X,y,beta,phi,omega,kappa,gamma_ij){
  # Define the function for optimizing omega at stage-two in M-step.


  name <- names(omega)
  gamma_hat <- apply(gamma_ij,2,sum)
  omega <- gamma_hat/sum(gamma_hat)
  omega <- omega[-length(omega)]
  names(omega) <- name
  return(omega)
}

.omega.optimize<-function(X,y,beta,phi,omega,kappa,lambda,weight,gamma_ij){
  # Define the function for optimizing omega.


  J <- length(kappa)
  n <- dim(X)[1]
  p <- dim(X)[2]
  omega_complete <- c(omega, 1.0 - sum(omega))
  gamma_hat <- apply(gamma_ij,2,sum)
  omega.updated <- c(lambda*weight*omega,0)-gamma_hat
  omega.updated[which(omega.updated*omega.updated[length(omega.updated)]<0)] = 0
  return(omega.updated/sum(omega.updated))
}

.omega.iterate<-function(X,y,beta,phi,omega,kappa,lambda,weight,gamma_ij)
{
  # Define the iterative function for the iterative solution to updated omega at M-step

  alpha=100
  i=1
  while(alpha>=0.01 && i<50)
  {
    tmpomega <- omega
    tmpomega <-c(tmpomega,1-sum(tmpomega))
    omega <-.omega.optimize(X,y,beta,phi,omega,kappa,lambda,weight,gamma_ij)
    alpha <- sum(abs(omega-tmpomega))/sqrt(sum(tmpomega^2))
    omega <- omega[-length(omega)]
    i=i+1
  }

  tmp.omega.extract <- .omega.extract(omega,weight)

  list(omega  = tmp.omega.extract$omega,
       kappa  = tmp.omega.extract$kappa,
       weight = tmp.omega.extract$weight)
}

.omega.extract <- function(omega,weight)
{
  # Filter and extract  updated inflated points

  b<-which(omega>0)
  name <- names(b)
  kappa <- gsub('omega','',name)
  kappa <- as.numeric(kappa)
  omega <- omega[b]
  weight <- weight[b]
  list(omega=omega,kappa=kappa,weight=weight)
}
.pry_nb<-function(X,y,beta,phi)
{
  ## Define likelihood function for NB process.

  J <- length(kappa)
  n <- dim(X)[1]
  p <- dim(X)[2]
  mu <- exp(drop(X%*%beta))
  # tmp1 <- rep(0, n)
  # for(i in 1:n){
  #   if(y[i]!= 0){
  #     tmp1[i] <- sum(log(0:(y[i]-1)+1.0/phi))
  #   }
  #   tmp1[i] <- exp(tmp1[i])/gamma(y[i] + 1)
  # }
  tmp1 <- lgamma( y + 1.0/phi)-lgamma(y+1)-lgamma(1.0/phi)
  tmp1 <- exp(tmp1)

  tmp2 <-(phi*mu/(1.0 + phi*mu))^y

  tmp3 <-(1.0/(1.0 + phi*mu))^(1.0/phi)
  prynb <- tmp1 * tmp2 * tmp3
  prynb[which(prynb == 0)] = 1e-300
  return(prynb)
}

.stageone<-function(X,y,nlambda,
                   lambda.min.ratio,lambda,
                   control){
  # Define function to find the first-stage estimator by applying the iterative


  ini <- .init(X,y,nlambda,lambda.min.ratio,lambda,control)
  p  <- dim(X)[2]
  if(!is.null(control$start)){
    kappa0 <- control$start$kappa
    omega <- control$start$omega
    beta <- control$start$beta
    phi <- control$start$phi
    weight = 1.0/omega
    if (length(kappa0) != length(omega))
      stop("the number of 'kappa' should be equal with the number of 'omega' ")
    if (length(beta) != p)
      stop("error with the length of starting 'beta'")
    if (is.null(phi))
      stop("value of 'phi' is not supplied ")
  }else{
  kappa0 <- ini$kappa
  beta <- ini$beta
  phi <- ini$phi
  omega<- ini$omega
  weight <- ini$weight}

  J <- length(kappa0)
  lambda_sequence <- ini$lambda


  bics <- rep(0,nlambda)
  betaomegaphis <- matrix(0,nrow=p+J+1,ncol=nlambda)
  rownames(betaomegaphis) <- c(paste0("beta", 0:(p-1)), paste0("omega", kappa0),"phi")
  pb <- txtProgressBar(style=3)
  for(i in 1:nlambda){
    lambda <- lambda_sequence[i]
    tmp <- .iterate(X,y,beta,phi,omega,kappa0,lambda,weight,control)
    beta_hat <- tmp$beta
    phi_hat <- tmp$phi
    omega_hat <- tmp$omega
    kappa <- tmp$kappa
    name <- names(omega_hat)
    betaomegaphis[1:p,i] <- beta_hat
    betaomegaphis[name,i] <- omega_hat
    betaomegaphis[p+J+1,i] <- phi_hat
    bics[i] <- .bic(X,y,beta_hat,phi_hat,kappa,omega_hat,control$thresh)
    setTxtProgressBar(pb, i/nlambda)
  }

  close(pb)

  a<-which.min(bics)
  args_stage1<-betaomegaphis[,a]

  list(betaomegaphis = betaomegaphis, a=a,bics = bics,
       args=args_stage1,lambda=lambda_sequence)
}


.minb.fit <- function(X, y, nlambda,
                     lambda.min.ratio, lambda,
                     control){

  # Define function to find the final estimator by applying the iterative
  # EM algorithm


  n <- dim(X)[1]
  X <- cbind(rep(1, n),X)
  p <- dim(X)[2]
  stage_1 <- .stageone(X,y,nlambda,lambda.min.ratio,lambda,control)
  betaomegaphis = stage_1$betaomegaphis
  bics = stage_1$bics
  lambda_sequence =stage_1$lambda
  args <- stage_1$args
  beta <- args[1:p]
  tmp <- args[-c(1:p)]
  omega <- tmp[-length(tmp)]
  phi <- tmp[length(tmp)]
  kappa <-as.numeric(gsub('omega','',names(which(omega>0))))
  omega <- omega[which(omega>0)]
  if(!length(omega)){
    omega = kappa = NULL
  }else{
  result <- .iterate2(X,y,beta,phi,omega,kappa,control)

  omega <- result$omega*(result$omega >= control$thresh)
  kappa <- as.numeric(gsub('omega','',names(which( omega > 0))))
  omega <- omega[which(omega > 0)]
  }

  rval <- list(omega = omega,
               coefficients = result$beta,
               phi = result$phi,
               kappa = kappa,
               betaomegaphis = betaomegaphis,
               bics = bics,
               lambda_sequence = lambda_sequence)
  return(rval)
}
.bic_p <- function(X,y,beta,kappa,omega,thresh){
  ## Define function for calculating the value of BIC criterion for mip.
  n <- dim(X)[1]
  p <- dim(X)[2]
  df <- p + sum(omega >= thresh)
  loss1 <- .loss_p(X,y,beta,omega,kappa)
  bic <- 2*loss1+log(n)*df
}

.loss_p<-function(X,y,beta,omega,kappa){
  ## Define function for calculating the log-likelihood  for mip
  J <- length(kappa)
  n <- dim(X)[1]
  p <- dim(X)[2]
  mu <- exp(drop(X%*%beta))
  delta <- .Delta(y, kappa)

  tmp1 <- 1.0 / gamma(y + 1)
  tmp2 <- mu^y
  tmp3 <- exp(-1.0 * mu)
  if(is.null(omega)){
    log_pry = log(tmp1 *tmp2 *tmp3)
  }else{
    location <- apply(delta, 1, function(x) which(x==1))
    omega_complete <- c(omega, 1.0 - sum(omega))
    omega_prime <- omega_complete[location]*(location!=(J+1))
    log_pry <- log(omega_prime + omega_complete[J+1] * tmp1 *tmp2 *tmp3)}

  loss <- -1.0 * sum(log_pry)
  return(loss)
}

.init_p<-function(X,y,nlambda,lambda.min.ratio,lambda=NULL,control){
  #  compute initialized values for mip algorithm
  p <- dim(X)[2]
  n <- dim(X)[1]
  kappa0 <- sort(unique(y))
  J <- length(kappa0)
  omega0 <- sapply(kappa0, function(x) sum(y == x))
  omega0 <- omega0 / length(y)
  substract_value <- min(omega0) / 2
  omega0 <- omega0 - substract_value
  weight0<-1.0/abs(omega0)
  names(omega0) <- paste0("omega", kappa0)
  names(weight0) <- paste0("omega", kappa0)
  data <- as.data.frame(cbind(X, y))
  colnames(data) <- c(paste("beta", 0:(p-1), sep=""), "y")

  omega_start <- omega0

  tmp <- glm(y~., data = data[,-1],family = poisson())
  beta_start <- tmp$coefficients

  if(is.null(lambda)){
    if (lambda.min.ratio >= 1)
      stop("lambda.min.ratio should be less than 1")

    for (i in seq(0.1*n,10*n,ifelse(n>=1000,50,10))){
      lambdamax = i
      omega.lambda.max = .iterate_p(X,y,beta_start,omega_start,kappa0,lambdamax,weight0,control)$omega
      if (length(omega.lambda.max)==0){
        break
      }
    }
    lambdamin <- lambdamax * lambda.min.ratio
    loghi <- log(lambdamax)
    loglo <- log(lambdamin)
    logrange <- loghi - loglo
    interval <- -logrange / (nlambda - 1)
    lambda <- exp(seq(loghi, loglo, interval))}
  else{
    if (any(lambda < 0))
      stop("lambdas should be non-negative")
    if (nlambda != length(lambda))
      stop("nlambda should be equal with the number of the lambda")
  }

  list(beta=beta_start,omega=omega0,lambda=lambda,weight=weight0,kappa=kappa0)
}

.omega.iterate.p<-function(X,y,beta,omega,kappa,lambda,weight,gamma_ij)
{
  alpha=100
  i=1
  while(alpha>=0.01 && i< 50)
  {
    tmpomega<-omega
    tmpomega<-c(tmpomega,1-sum(tmpomega))
    omega <- .omega.optimize.p(X,y,beta,omega,kappa,lambda,weight,gamma_ij)

    alpha<- sum(abs(omega-tmpomega))/sqrt(sum(tmpomega^2))
    omega<-omega[-length(omega)]
    i=i+1
  }
  tmp<-.omega.extract(omega,weight)
  omega<-tmp$omega
  kappa<-tmp$kappa
  weight<-tmp$weight
  list(omega=omega,kappa=kappa,weight=weight)
}

.omega.optimize.p<-function(X,y,beta,omega,kappa,lambda,weight,gamma_ij){
  J <- length(kappa)
  n <- dim(X)[1]
  p <- dim(X)[2]
  omega_complete <- c(omega, 1.0 - sum(omega))
  gamma_hat<-apply(gamma_ij,2,sum)
  omega.updated <- c(lambda*weight*omega,0)- gamma_hat
  omega.updated[which(omega.updated*omega.updated[length(omega.updated)]<0)] = 0
  return(omega.updated/sum(omega.updated))
}

.gammaij_p<-function(X,y,beta,omega,kappa){
  J<-length(kappa)
  n<-nrow(X)
  delta1<-.Delta(y,kappa)
  delta1[,J+1]<-rep(1,n)
  pry<-.pry_p(X,y,beta)
  tmp1<-(1-sum(omega))*pry
  tmp2<-matrix(rep(omega,n),nrow=n,byrow=TRUE)
  tmp2<-cbind(tmp2,tmp1)
  tmp3<-delta1*tmp2
  SUM<-apply(tmp3,1,sum)
  SUM<-matrix(rep(SUM,J+1),nrow=n,byrow =FALSE)
  gammaij<-tmp3/SUM
  return(gammaij)
}
.pry_p<-function(X,y,beta){
  J <- length(kappa)
  n <- dim(X)[1]
  p <- dim(X)[2]
  mu <- exp(drop(X%*%beta))
  tmp1 <- 1.0 / gamma(y + 1)
  tmp2 <- mu^y
  tmp3 <- exp(-1.0 * mu)
  pry<-tmp1*tmp2*tmp3
  return(pry)
}
.grbeta_p<-function(X,y,beta,omega,kappa,gamma_ij){
  J <- length(kappa)
  n <- dim(X)[1]
  p <- dim(X)[2]
  mu <- exp(drop(X%*%beta))
  dmu<-mu*X
  delta<-.Delta(y, kappa)
  location <- apply(delta, 1, function(x) which(x==1))
  omega_complete <- c(omega, 1.0 - sum(omega))
  tmp1 <- 1.0 / gamma(y + 1)
  tmp2 <- mu^y
  tmp3 <- exp(-1.0 * mu)
  pry <- omega_complete[J+1] * tmp1 *tmp2 *tmp3
  dtmp2_mu <- y * (mu^(y-1))
  dtmp3_mu <- -1.0 * exp(-1.0 * mu)
  dpry_mu <- omega_complete[J+1] * tmp1 *
    (dtmp2_mu * tmp3 + tmp2 * dtmp3_mu)
  dlogpry_mu <- dpry_mu/pry
  gamma_iJ<-gamma_ij[,J+1]
  gradient<-gamma_iJ*dlogpry_mu
  A<-gradient*dmu
  dobj_beta<-apply(A, 2, sum)
  dobj_beta<- -1.0*dobj_beta
  return(dobj_beta)
}



.iterate_p <- function(X,y,beta,omega,kappa,lambda,weight,control){
  alpha=100
  i=1
  while(alpha > control$toler && i< control$maxit){
    tmpbeta=beta
    tmpomega=omega
    gamma_ij=.gammaij_p(X,y,beta,omega,kappa)
    beta<-optim(beta, .obj_p, .grbeta_p,method="BFGS",
                X=X,y=y,omega=omega,kappa=kappa,gamma_ij=gamma_ij)$par
    tmp<-.omega.iterate.p(X,y,beta,omega,kappa,lambda,weight,gamma_ij)
    omega<-tmp$omega
    if(!length(omega)){
      omega = kappa = weight = NULL
      fit.poisson <- glm(y~., data = as.data.frame(cbind(X,y))[,-1],control=glm.control(maxit=100))
      beta <- fit.poisson$coefficients
      break
    }
    if(sum(omega)>=1){
      substract_value <- min(omega) / 2
      omega <- omega - substract_value
    }
    kappa<-tmp$kappa
    weight<-tmp$weight
    name<-names(omega)
    tmpomega2<-rep(0,length(tmpomega))
    names(tmpomega2)<-names(tmpomega)
    tmpomega2[name] <- omega
    norm1 <- sum(abs(beta-tmpbeta))/sqrt(sum(tmpbeta^2))
    norm2 <- sum(abs(tmpomega2-tmpomega))/sqrt(sum(tmpomega^2))
    alpha <- max(norm1, norm2)
    i <- i+1
    if(i == control$maxit)
      warning("fitting results did not converge -- increase 'maxit'?")
  }
  list(omega=omega,kappa=kappa,weight=weight,beta=round(beta,6))
}

.stageone_p<-function(X,y,
                     nlambda,lambda.min.ratio,lambda,
                     control){
  # Define function to find the first-stage estimator by applying the iterative
  # EM algorithm

  ini <- .init_p(X,y,nlambda,lambda.min.ratio,lambda,control)
  p <- dim(X)[2]
  if(!is.null(control$start)){
    kappa0 <- control$start$kappa
    omega <- control$start$omega
    beta <- control$start$beta
    weight = 1.0/omega
    if (length(kappa0) != length(omega))
      stop("the number of 'kappa' should be equal with the number of 'omega' ")
    if (length(beta) != p)
      stop("error with the length of starting 'beta'")
  }else{
    kappa0 <- ini$kappa
    beta <- ini$beta
    omega<-ini$omega
    weight <- ini$weight}

  lambda_sequence <- ini$lambda
  J <- length(kappa0)

  bics<-rep(0,nlambda)
  betaomegas <- matrix(0,nrow=p+J,ncol=nlambda)
  rownames(betaomegas) <- c(paste0("beta", 0:(p-1)), paste0("omega", kappa0))

  pb <- txtProgressBar(style=3)

  for(i in 1:nlambda){
    lambda <- lambda_sequence[i]
    tmp <- .iterate_p(X,y,beta,omega,kappa0,lambda,weight,control)
    beta_hat <- tmp$beta
    omega_hat <- tmp$omega
    kappa <- tmp$kappa
    name <- names(omega_hat)
    betaomegas[1:p,i] <- beta_hat
    betaomegas[name,i] <- omega_hat
    bics[i] <- .bic_p(X,y,beta_hat,kappa,omega_hat,control$thresh)
    setTxtProgressBar(pb, i/nlambda)
  }

  close(pb)

  a <- which.min(bics)
  args_stage1<-betaomegas[,a]

  list(betaomegas = betaomegas, a=a,bics = bics,
       args=args_stage1,lambda=lambda_sequence)
}

.iterate2_p<-function(X,y,beta,omega,kappa,control){

  alpha=100
  i=1
  while(alpha > control$toler && i< control$maxit){
    tmpbeta=beta
    tmpomega=omega
    gamma_ij=.gammaij_p(X,y,beta,omega,kappa)
    omega<-.omega_em_p(X,y,beta,omega,kappa,gamma_ij)
    beta<-optim(beta, .obj_p, .grbeta_p,method=control$method,
                X=X,y=y,omega=omega,kappa=kappa,gamma_ij=gamma_ij)$par
    norm1 <- sum(abs(beta-tmpbeta))/sqrt(sum(tmpbeta^2))
    norm3 <- sum(abs(omega-tmpomega))/sqrt(sum(tmpomega^2))
    alpha <- max(norm1,norm3)
    i<-i+1
    if(i == control$maxit)
      warning("fitting results did not converge -- increase 'maxit'?" )
  }
  list(omega=omega,kappa=kappa,beta=round(beta,6))
}

.omega_em_p<-function(X,y,beta,omega,kappa,gamma_ij){
  name<-names(omega)
  gamma_hat<-apply(gamma_ij,2,sum)
  omega<-gamma_hat/sum(gamma_hat)
  omega<-omega[-length(omega)]
  names(omega)<-name
  return(omega)
}

.obj_p<-function(X,y,beta,omega,kappa,gamma_ij){
  J <- length(kappa)
  pry <- .pry_p(X,y,beta)
  tmp1 <- log((1-sum(omega))*pry)
  obj <- sum(gamma_ij[,J+1]*tmp1)
  obj <- -1.0*obj
  return(obj)
}

.mip.fit <- function(X, y, nlambda,
                    lambda.min.ratio, lambda,
                    control){

  # Define function to find the final estimator by applying the iterative
  n <- dim(X)[1]
  X <- cbind(rep(1, n),X)
  p <- dim(X)[2]
  stage_1 <- .stageone_p(X,y,nlambda,lambda.min.ratio,lambda,control)
  betaomegas = stage_1$betaomegas
  bics = stage_1$bics
  lambda_sequence =stage_1$lambda
  args <- stage_1$args

  beta <- args[1:p]
  omega <- args[-c(1:p)]
  kappa <-as.numeric(gsub('omega','',names(which(omega>0))))
  omega <- omega[which(omega>0)]

  if(!length(omega)){
    omega = kappa = NULL
  }else{
    result <- .iterate2_p(X,y,beta,omega,kappa,control)

    omega <- result$omega*(result$omega >= control$thresh)
    kappa <- as.numeric(gsub('omega','',names(which( omega > 0))))
    omega <- omega[which(omega > 0)]
  }

  rval <- list(omega = omega,
               coefficients = result$beta,
               kappa = kappa,
               betaomegas = betaomegas,
               bics = bics,
               lambda_sequence = lambda_sequence)
  return(rval)
}

# main function
minb <- function(X, y, dist = 'nb',nlambda = 20,
                 lambda.min.ratio = ifelse(dim(X)[1] < dim(X)[2], 0.01,0.0001 ),
                 lambda = NULL,
                 control = minb.control(),
                 modelfit = FALSE){
  # Define function to find the final estimator by applying the iterative
  # EM algorithm
  n = dim(X)[1]
  # if (!requireNamespace("MASS", quietly = TRUE)) {
  #   stop("Package MASS needed for glm method to work. Please install it.",
  #        call. = FALSE)}else{
  #         require(MASS)
  #        }
  if(! is.character(dist))
    dist = as.character(dist)

  if (!dist %in% c('nb','poisson') ) {
    print(dist)
    stop("'dist' not recognized ")
  }

  rval <- switch (dist,
                 nb = .minb.fit(X, y, nlambda, lambda.min.ratio, lambda,control),
                 poisson = .mip.fit(X, y, nlambda, lambda.min.ratio, lambda,control) )

  X <- cbind(rep(1,n),X)
  if(!is.null(colnames(X))){
      names(rval$coefficients ) = c('Intercept',colnames(X)[-1])
      if(dist == 'nb') rownames(rval$betaomegaphis)[1:dim(X)[2]]= c('Intercept',colnames(X)[-1])
      if(dist == 'poisson') rownames(rval$betaomegas)[1:dim(X)[2]] = c('Intercept',colnames(X)[-1])
  }
  mu.fit <- drop(exp(X %*% rval$coefficients ))

  if(!is.null(rval$omega)){
  rval$fitted.values = (1-sum(rval$omega))*mu.fit +
    c(rval$omega%*%rval$kappa)}else{
      rval$fitted.values = mu.fit }
  rval$residuals = y - rval$fitted.values

  if (dist == 'nb' && rval$phi <= control$thresh){
    class(rval) = 'Mip'
    if(is.null(rval$omega)){
      class(rval) = 'Poisson'
      rval$loglik = -1.0 * .loss_p(X, y, rval$coefficients ,rval$omega,rval$kappa)
      }
    rval$loglik = -1.0 * .loss_p(X, y, rval$coefficients ,rval$omega,rval$kappa)
  }else if(is.null(rval$omega)){
      class(rval)  = switch(dist,nb='Negavtive binomial',poisson='Poisson')
  }else{
    class(rval) = dist
  rval$loglik <- switch (dist,
                         nb = -1.0 * .loss(X, y, rval$coefficients , rval$phi, rval$omega,rval$kappa),
                         poisson = -1.0 * .loss_p(X, y, rval$coefficients, rval$omega,rval$kappa) )
  }

  if(modelfit){
    rval$modelfit = list(coefficients = switch(dist, nb = rval$betaomegaphis,
                                               poisson = rval$betaomegas),
                         lambda = rval$lambda_sequence,
                         bics = rval$bics)
  }
  rval$betaomegaphis= rval$betaomegas =rval$bics = rval$lambda_sequence = NULL

  return(rval)
}

minb.control <- function(method = 'BFGS',
                         toler = 1e-03,
                         maxit = 500,
                         thresh = 1e-02,
                         start = NULL,
                         ...){
  rval = list(method = method, toler = toler, maxit = maxit,
              thresh = thresh, start = start)
  rval = c(rval,list(...))
  if (!rval$method %in% c('BFGS','L-BFGS-B','CG') ){
    print(rval$method)
    stop("'method' not recognized, must be modified within 'BFGS','L-BFGS-B',and 'CG' ")
    }
  if (!is.numeric(toler) || toler <= 0)
      stop("value of 'toler' must be > 0")
  if (!is.numeric(thresh) || thresh <= 0)
      stop("value of 'thresh' must be > 0")
  if (!is.numeric(maxit) || maxit <= 0)
      stop("maximum number of iterations must be non-negative integer")

  return(rval)

}
