# =============================================
# Algorithm: Direct Lagrangian Approach
# =============================================

library(quantreg)

####### function for estimating the mean outcome
mean_est<-function(beta, x, y, a, prob){
  d <- 1 / (1 + exp(x %*% beta))
  c <- a * d  +(1-a) * (1-d)
  wts <- d * 1/prob + (1-d) * (1/(1-prob))
  val <- mean(c * wts * y)
  return(val)
}

####### function for estimating the quantile outcome
quant_est <- function(beta, x, y, a, prob, tau){
  d <- 1 / (1 + exp(x %*% beta))
  c <- a * d + (1-a) * (1-d)
  wts <- d * 1/prob + (1-d) * (1/(1-prob))
  model <- rq(y ~ 1, weights = c*wts, tau=tau)
  return(coefficients(model)[1])
}


####### function for solving the inner maximization problem
solve_beta <- function(x, y, a, prob, tau, lambda){
  d <- dim(x)[2]
  beta <- rep(0, d)
  beta0 <- beta
  iter <- 0
  tol <- 1e-6
  err <- 10000
  lambda0 <- 0
  # Gradient ascent loop
  while ((err > tol) && (iter < 10000)){
    grad <- -quant_grad(beta, x, y, a, prob, tau) + mean_grad(beta, x, y, a, prob) * lambda
    lambda1 <- (1 + sqrt(1 + 4 * lambda0)) / 2
    gamma <- (1 - lambda0) / lambda1
    bint = beta - 0.03*grad
    beta <- (1 - gamma) * bint + gamma * beta0
    iter <- iter + 1
    lambda0 <- lambda1
    beta0 <- bint
    dim(beta) = NULL
  }
  return(beta)
}


quant_grad <- function(beta, x, y, a, prob, tau) {
  n <- dim(x)[1]
  g <- dim(x)[2]
  xb <- x %*% beta
  dt <- 1 / (1 + exp(xb))
  tmp <- exp(xb) / ((1 + exp(xb))^2)
  wts <- dt * 1/prob + (1 - dt) * (1/(1 - prob))
  ct <- a * dt + (1 - a) * (1 - dt)
  cb <- ct * wts

  q_index <- quant_index(cb, y, tau)
  
  # Order
  ord <- order(y)
  c_before <- ord[1:q_index]
  c_after <- ord[(q_index + 1):n]
  grdCq1 <- (2 * a[c_before] * tmp[c_before] - tmp[c_before]) %*% 
    x[c_before, , drop = FALSE]
  grdCq2 <- (2 * a[c_after] * tmp[c_after] - tmp[c_after]) %*% 
    x[c_after, , drop = FALSE]
  y_diff <- y[ord[q_index + 1]] - y[ord[q_index]]
  grdq <- y_diff * (-grdCq1 * tau + grdCq2 * (1 - tau)) * runif(g)

  return(grdq)
}

quant_index <- function(c, y, tau) {
  ord <- order(y)
  c_order <- c[ord]
  n <- length(y)
  cum_fwd <- cumsum(c_order)  
  cum_rev <- sum(c_order) - cumsum(c_order)  
  
  tmp <- tau * cum_fwd > (1 - tau) * cum_rev
  qind <- min(which(tmp == 1))
  
  return(qind)
}

mean_grad <- function(beta, x, y, a, prob) {
  d <- length(beta)
  m0 <- mean_est(beta, x, y, a, prob)
  grd <- rep(0, d)
  
  for (j in 1:d) {
    tmp <- rep(0, d)
    tmp[j] <- 0.01 * (2 * rbinom(1, size = 1, prob = 0.5) - 1)
    grd[j] <- -(mean_est(beta + tmp, x, y, a, prob) - m0) / tmp[j]
  }
  
  return(grd)
}


####### function for estimating proposed Safe-IDR
qmestimate <- function(x,y,a,prob,tau,mcon){
  lambda1 = 0
  lambda2 = 40
  
  beta1  = solve_beta(x,y,a,prob,tau,lambda1)
  beta10 = beta1
  val1m = mean_est(beta1,x,y,a,prob)
  val1q = quant_est(beta1,x,y,a,prob,tau)
  val1  = -val1q + lambda1*(-val1m+mcon)
  
  beta2  = solve_beta(x,y,a,prob,tau,200)  
  beta20 = beta2
  val2m = mean_est(beta2,x,y,a,prob)
  val2q = quant_est(beta2,x,y,a,prob,tau)
  if (val2m < mcon)
  {
    return(9999999)
  }
  val2  = -val2q + lambda2*(-val2m+mcon)
  
  gr  = (1+sqrt(5))/2
  
  lambda3 = lambda2-(lambda2-lambda1)/gr
  lambda4 = lambda1+(lambda2-lambda1)/gr
  
  while (abs(lambda3-lambda4)>0.05){
    # print(abs(lambda3-lambda4))
    beta3 = solve_beta(x,y,a,prob,tau,lambda3)
    val3m = mean_est(beta3,x,y,a,prob)
    val3q = quant_est(beta3,x,y,a,prob,tau)
    val3  = -val3q + lambda3*(-val3m+mcon)
    
    
    tmp = gr_search(beta1,beta2,x,y,a,prob,tau,lambda3,mcon)
    if (tmp$val<val3){
      val3 = tmp$val
      beta3 = tmp$beta
    }
    
    beta4 = solve_beta(x,y,a,prob,tau,lambda4)
    val4m = mean_est(beta4,x,y,a,prob)
    val4q = quant_est(beta4,x,y,a,prob,tau)
    val4  = -val4q + lambda4*(-val4m+mcon)
    
    tmp = gr_search(beta1,beta2,x,y,a,prob,tau,lambda4,mcon)
    if (tmp$val<val4){
      val4 = tmp$val
      beta4 = tmp$beta
    }
    
    if (val3>val4)
    {
      lambda2 = lambda4
      beta2   = beta4
      lambda4 = lambda3
      beta4   = beta3
      lambda3 = lambda2-(lambda2-lambda1)/gr
    } else{
      lambda1 = lambda3
      beta1   = beta3
      lambda3 = lambda4
      beta3   = beta4
      lambda4 = lambda1+(lambda2-lambda1)/gr
    }
  }
  
  if (val3<val4)
  {
    beta  = beta3
    hatQ  = mean_est(beta3,x,y,a,prob)
    hatQ2 = quant_est(beta3,x,y,a,prob,tau)
    dualv = hatQ2 + lambda3*(hatQ - mcon)
  } else{
    beta  = beta4
    hatQ  = mean_est(beta4,x,y,a,prob)
    hatQ2 = quant_est(beta4,x,y,a,prob,tau)
    dualv = hatQ2 + lambda4*(hatQ - mcon)
  }
  
  {
    beta1 = beta10
    beta2 = beta20
    while (sum(abs(beta1-beta2))>0.01){
      beta3 = (beta1+beta2)/2
      val3m = mean_est(beta3,x,y,a,prob)
      if (val3m < mcon){
        beta1 = beta3
      } else {
        beta2 = beta3
      }
    }
    val2m = mean_est(beta2,x,y,a,prob)
    val2q = quant_est(beta2,x,y,a,prob,tau)
  }
  
  if (val2q > hatQ2){
    beta  = beta2
    hatQ  =  val2m
    hatQ2 = val2q
  } else if (hatQ<mcon){
    beta  = beta2
    hatQ  = val2m
    hatQ2 = val2q
  }
  
  dgap = abs(dualv-hatQ2)
  
  
  summary <- list(
    beta = beta,
    mval = hatQ,
    qval = hatQ2,
    dval = dualv,
    dgap = dgap
  )
  return(summary)
}


gr_search <- function(beta1,beta2,x,y,a,prob,tau,lambda,mcon){
  gr <- 1.618
  lambda1 <- 0
  lambda2 <- 1
  val1    <- -quant_est(beta1,x,y,a,prob,tau) + lambda*(-mean_est(beta1,x,y,a,prob)+mcon)
  val2    <- -quant_est(beta2,x,y,a,prob,tau) + lambda*(-mean_est(beta2,x,y,a,prob)+mcon)
  lambda3 <- lambda2-(lambda2-lambda1)/gr
  lambda4 <- lambda1+(lambda2-lambda1)/gr
  while (abs(lambda3-lambda4)>0.01){
    beta3 <- beta1 + lambda3*(beta2-beta1) 
    val3 <- -quant_est(beta3,x,y,a,prob,tau) + lambda*(-mean_est(beta3,x,y,a,prob)+mcon)
    beta4 <- beta1 + lambda4*(beta2-beta1)
    val4 <- -quant_est(beta4,x,y,a,prob,tau) + lambda*(-mean_est(beta4,x,y,a,prob)+mcon)
    if (val3<val4)
    {
      lambda2 <- lambda4
      lambda4 <- lambda3
      lambda3 <- lambda2-(lambda2-lambda1)/gr
    } else{
      lambda1 <- lambda3
      lambda3 <- lambda4
      lambda4 <- lambda1+(lambda2-lambda1)/gr
    }
    
  }
  beta <- beta3
  val <- val3
  
  tmp <- c()
  tmp$beta <- beta
  tmp$val <- val
  return(tmp)
}


