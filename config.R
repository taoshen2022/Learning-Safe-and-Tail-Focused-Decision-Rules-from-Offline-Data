b0 <- 1
bx1 <- 1
bx2 <- -1
bx3_f <- function(x) { x^3 }
bx4_f <- function(x) { exp(x) }

ba <- 3  
bax1 <- -5  
bax2 <- 2  
bax3 <- -3  
bax4 <- 1  

b0_scale <- 1
ba_scale <- 1
bax1_scale <- 1
bax2_scale <- 1
bax3_scale <- 1
bax4_scale <- 1


gen_x <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  x3 <- runif(n)
  x4 <- runif(n)
  x <- cbind(1, x1, x2, x3, x4)
  return(x)
}


gen_error <- function(n) {
  error <- rnorm(n, sd = 1)
  return(error)
}


gen_tp <- function(x) {
  gamma0 <- -0.5
  gamma1 <- 0.5
  tp <- exp(gamma0 + gamma1 * (x[, 2] + x[, 3] + x[, 4] + x[, 5])) / 
    (1 + exp(gamma0 + gamma1 * (x[, 2] + x[, 3] + x[, 4] + x[, 5])))
  return(tp)
}


gen_y <- function(x, a, error) {
  y <- b0 + bx1 * x[, 2] + bx2 * x[, 3] + bx3_f(x[, 4]) + bx4_f(x[, 5]) + 
    ba * a + bax1 * a * x[, 2] + bax2 * a * x[, 3] + bax3 * a * x[, 4] + bax4 * a * x[, 5] + 
    (b0_scale + ba_scale * a + bax1_scale * a * x[, 2] + bax2_scale * a * x[, 3] + 
       bax3_scale * a * x[, 4] + bax4_scale * a * x[, 5]) * error
  return(y)
}


# Generate example data
generate_data <- function(n){
  x <- gen_x(n)
  error <- gen_error(n)
  tp <- gen_tp(x)
  a <- rbinom(n, 1, tp)
  y <- gen_y(x, a, error)
  
  logit <- glm(a ~ x[, 2] + x[, 3] + x[, 4] + x[, 5], family = binomial, epsilon = 1e-14)
  ph <- as.vector(logit$fit)
  
  list(x = x, y = y, a = a, ph = ph, error = error)
}