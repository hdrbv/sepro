#' EM Function
#' 
#' 
#' This function realized EM algorithm (Expectation-Maximization algorithm)
#' for data clustering. In this case, for mixture models.
#' 
#' In initial step (Step_0) explorer must determine the initial gaussian mixture model 
#' parameters. Explorer can take any random parameters or (using surface analysis and some assumptions) 
#' take anothers.
#' 
#' In the first step (called E_Step) with function E_Step explorer calculate 
#' the posterior probabilities (which named "posterior.df") for each 
#' item of initial dataset.
#' 
#' In the second step (called M_Step) with function M_Step explorer
#' update component parameters (using likelihood function).
#' 
#' After that, explorer repeat M_Step (100 times or until difference between current value of 
#' logarithm of likelihood function and previous value of logarithm of likelihood function
#' will be less than 10^(-6) 
#' in other words: logarithm of likelihood function didn’t change much). 
#' 
#' In the end explorer has updated probabilities for each item and
#' updated parameters for each distribution (explorer should look at E.Step and M.Step)
#' 
#' @param x0 - input data (vector)
#' @param k - amount of clusters (mixture components) 
#' @return Parameters for each distribution 
#' @author hdrbv

EM <- function(x0, k){  #Beginning of EM function
  #Function for E-Step
  E_Step <- function(x, mu.vector, sd.vector, w.vector) {
    comp1.prod <- dnorm(x, mu.vector[1], sd.vector[1]) * w.vector[1]
    comp2.prod <- dnorm(x, mu.vector[2], sd.vector[2]) * w.vector[2]
    sum.of.comps <- comp1.prod + comp2.prod
    comp1.post <- comp1.prod / sum.of.comps
    comp2.post <- comp2.prod / sum.of.comps
    sum.of.comps.ln <- log(sum.of.comps, base = exp(1))
    sum.of.comps.ln.sum <- sum(sum.of.comps.ln)
    list("loglik" = sum.of.comps.ln.sum,
         "posterior.df" = cbind(comp1.post, comp2.post))
  }
  #Function for M-Step
  M_Step <- function(x, posterior.df) {
    comp1.n <- sum(posterior.df[, 1])
    comp2.n <- sum(posterior.df[, 2])
    comp1.mu <- 1/comp1.n * sum(posterior.df[, 1] * x)
    comp2.mu <- 1/comp2.n * sum(posterior.df[, 2] * x)
    comp1.var <- sum(posterior.df[, 1] * (x - comp1.mu)^2) * 1/comp1.n
    comp2.var <- sum(posterior.df[, 2] * (x - comp2.mu)^2) * 1/comp2.n
    comp1.w <- comp1.n / length(x)
    comp2.w <- comp2.n / length(x)
    list("mu" = c(comp1.mu, comp2.mu),
         "var" = c(comp1.var, comp2.var),
         "w" = c(comp1.w, comp2.w))
  }
  for (i in 1:100) {
    if (i == 1) {
      # Step_0. Initialization
      E.Step <- E_Step(x0, c(0,1), c(1,2), c(1/k, 1/k))
      M.Step <- M_Step(x0, E.Step[["posterior.df"]])
      cur.loglik <- E.Step[["loglik"]]
      loglik.vector <- E.Step[["loglik"]]
    } else {
      # Repeat E and M steps till convergence
      E.Step <- E_Step(x0, M.Step[["mu"]], sqrt(M.Step[["var"]]), M.Step[["w"]])
      M.Step <- M_Step(x0, E.Step[["posterior.df"]])
      loglik.vector <- c(loglik.vector, E.Step[["loglik"]])
      loglik.diff <- abs((cur.loglik - E.Step[["loglik"]]))
      if(loglik.diff < 1e-6) {
        break
      } else {
        cur.loglik <- E.Step[["loglik"]]
      }
    }
  }
  print(loglik.vector)
  print(E.Step)
  print(M.Step)
} #Ending of EM function
