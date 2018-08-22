## Generalized Extreme Value Binomial Regression Link Function
## (c) Daniel S. Mork 2018


## This function allows for model building where bernoulli response data has a
## very low/high probability of occurance. It is recommended to compare
## predictive values for a test set of data to best optimize the xi value of the
## gev link function. A graphical simulation of this function under different
## values of xi is available here: https://www.desmos.com/calculator/mvqzlylab1

## NOTES:
## - Currently this link will work with glm and mgcv::bam for GAMs, but needs
##   additional functions to work with mgcv::gam.
## - Using rev = T may need to have 'etastart' specified in glm


## REF:
##	- Wang, X., & Dey, D. K. (2010). Generalized extreme value regression for binary response data:
##	    An application to B2B electronic payments system adoption. The Annals of Applied Statistics, 4(4), 2000-2023.
##	- Calabrese, R., & Osmetti, S. A. (2013). Modelling small and medium enterprise loan defaults as rare events:
##	    the generalized extreme value regression model. Journal of Applied Statistics, 40(6), 1172-1188.
##	- https://stackoverflow.com/questions/15931403/modify-glm-function-to-adopt-user-specified-link-function-in-r
##	- https://stackoverflow.com/questions/39793085/custom-link-function-works-for-glm-but-not-mgcv-gam
##  - https://stats.stackexchange.com/questions/46523/how-to-simulate-artificial-data-for-logistic-regression/46525
##	- https://gist.github.com/frbl/1411ee1df13154bd22092a7894503eb2#L22
##	- https://github.com/cran/mgcv/blob/master/R/bam.r
##	- https://github.com/cran/mgcv/blob/master/R/efam.r
##	- R help files for make.link, binomial

## USE:
## - Binomial link function based on GEV cdf. Parameter xi modifies shape
##    of distribution.
## - xi = 0 => Gumbel (also cloglog link)
## - rev will reflect the Gumbel distribution over x=0 (i.e. reverses skew)
## - xi < 0: Use when small amt of 1's compared to 0's
## - xi > 0: Use when large amt of 1's compared to 0's
## - Compare parameterizations of with deviance or AIC
## - glm: add mustart = rep(mean(y), length(y)) {where y is binomial response}
gev <- function (xi = 0, rev = F) {

  # Link function: Map cumulative probability (mu) to quantile value
  linkfun <- function(mu) {
    if (xi == 0) # Gumbel
      ifelse(!rev & validmu(mu),
             -log(-log(mu)),
             log(-log(1 - mu)))

    else # GEV
       ((-log(mu))^(-xi) - 1) / xi
  }

  # Inverse link function: Map quantile value to cumulative probability (cdf)
  linkinv <- function(eta) {

    if (xi == 0) { # Gumbel
      ifelse(!rev & valideta(eta),
        pmax(pmin(exp(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps),
        pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps))

    } else { # GEV
      ifelse(valideta(eta),
             pmax(pmin(exp(-(1 + xi*eta)^(-1/xi)), 1 - .Machine$double.eps), .Machine$double.eps),
             ifelse(xi > 0, .Machine$double.eps, 1 - .Machine$double.eps))

    }
  }

  # Derivative, d.mu/d.eta, of inverse link (pdf)
  mu.eta <- function(eta) {
    eta <- pmax(pmin(eta, 700), -700)
    if (xi == 0) { # Gumbel
      if (!rev)
        pmax(exp(-eta-exp(-eta)), .Machine$double.eps)
      else
        pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)

    } else { # GEV
      ifelse(valideta(eta),
        exp(-(1 + xi*eta)^(-1/xi)) * (1 + xi*eta)^(-(1 + 1/xi)), 0)
    }
  }

  # For GEV dist, there is a domain restriction, must specify valid eta values
  valideta <- function(eta) {
    ifelse(1 + xi*eta > 0, TRUE, FALSE)
  }

  # Specify valid probability values
  validmu <- function(mu) {
    ifelse(mu > 0 & mu < 1, TRUE, FALSE)
  }

  # Need to set inital eta/mu for calculation with mgcv::bam,
  # used 0 and corresponding mu value, as these are always valid
  initialize <- expression({
    etastart <- 0
    mustart <- linkinv(etastart)
  })

  link <- paste("gev(xi = ", xi, ")", sep = "")
  structure(list(linkfun = linkfun, linkinv = linkinv,
                 mu.eta = mu.eta, valideta = valideta, validmu = validmu,
                 initialize = initialize, name = link),
            class = "link-glm")
}

## Example
set.seed(1)
n = 10000
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
z <- x1 - 2*x2 + x3 - 6 + .05 * rnorm(n) # Shift center, add small noise
fam <- gev(-0.5)    # Specify model
p <- fam$linkinv(z) # Probability value from GEV link
plot(sort(p))       #
y <- rbinom(n, 1, p)# Draw bernoulli r.v. with prob. p
table(y)/n

# Fit model with correct GEV parameterization
dat <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)
glm(y ~ x1 + x2 + x3, dat, family = binomial(link = fam),
            mustart = rep(mean(y), length(y)))
# Compare results to logit, probit, cauchit, and gumbel links (logit, cloglog)
glm(y ~ x1 + x2 + x3, dat, family = binomial(link = logit))
glm(y ~ x1 + x2 + x3, dat, family = binomial(link = probit))
glm(y ~ x1 + x2 + x3, dat, family = binomial(link = cauchit))
glm(y ~ x1 + x2 + x3, dat, family = binomial(link = cloglog)) # gev(0)

# Fit to GEV link with different shape parameters, choose best model with AIC
xi.val <- seq(-1, 0, 0.1)
xi.val[which.min(sapply(xi.val, function(xi) glm(y = x1 + x2, dat,
                                   family = binomial(link = gev(xi)),
                                   mustart = rep(mean(y), length(y)))$aic))]
