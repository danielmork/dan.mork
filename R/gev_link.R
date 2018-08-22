#' Generalized Extreme Value Binomial Regression Link Function
#' @author Daniel S. Mork (c)2018
#'
#'
#' @description This function allows for model building where the binomial response
#' data is very unbalanced. It is recommended to compare
#' predictive values for a test set of data to best optimize the xi value of the
#' gev link function. A graphical simulation of this function under different
#' values of xi is available here: \url{https://www.desmos.com/calculator/mvqzlylab1}.
#' \cr\cr
#' The GEV binomial link function is based on the GEV cdf. Parameter
#' xi modifies the shape of the distribution, with xi = 0 being equivalent to
#' the \code{cloglog} link available in \code{R}. The link function often throws
#' errors due to its domain restriction and works better for values of xi < 0 and
#' when there are a large amount of 0's compared to 1's in the response. Setting
#' \code{rev=T} will reflect the \code{cloglog} link, reversing the skew of the
#' distribution.
#' \cr\cr
#' This function currently works in \code{glm} and \code{mgcv::bam} for GAMs, but does
#' not contain all the function needed to work in \code{mgcv::gam}.
#'
#'
#' @param xi Parameter in GEV distribution (negative values produce fewer errors)
#' @param rev If TRUE, reverses skew when xi = 0
#'
#' @return An object of class \code{link-glm}, which has components:\cr
#' \code{linkfun}  Link function \code{function(mu)} (inverse GEV cdf)\cr
#' \code{linkinv}  Inverse link function \code{function(eta)} (GEV cdf)\cr
#' \code{mu.eta}   Derivative of linkinv w.r.t eta \code{function(eta)}\cr
#' \code{valideta} Specifies domain of linkinv \code{function(eta)}\cr
#' \code{validmu} Specifies domain of linkfun, mu in (0,1) \code{function(mu)}\cr
#' \code{initialize} Specifies beginning values of mu and eta \code{expression}
#'
#' @export
#' @references
#' Wang, X., & Dey, D. K. (2010). Generalized extreme value regression for binary response data:\cr
#' 	    An application to B2B electronic payments system adoption. The Annals of Applied Statistics, 4(4), 2000-2023.\cr
#' 	Calabrese, R., & Osmetti, S. A. (2013). Modelling small and medium enterprise loan defaults as rare events:\cr
#' 	    the generalized extreme value regression model. Journal of Applied Statistics, 40(6), 1172-1188.\cr
#' 	\url{https://stackoverflow.com/questions/15931403/modify-glm-function-to-adopt-user-specified-link-function-in-r}\cr
#' 	\url{https://stackoverflow.com/questions/39793085/custom-link-function-works-for-glm-but-not-mgcv-gam}\cr
#' 	\url{https://stats.stackexchange.com/questions/46523/how-to-simulate-artificial-data-for-logistic-regression/46525}\cr
#' 	\url{https://gist.github.com/frbl/1411ee1df13154bd22092a7894503eb2#L22}\cr
#'  R help files for make.link, binomial
#'
#' @examples
#' set.seed(1)
#' n = 10000
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' x3 <- rnorm(n)
#' z <- x1 - 2*x2 + x3 - 6 + .05 * rnorm(n) # Shift center, add small noise
#' fam <- gev(-0.5)    # Specify model
#' p <- fam$linkinv(z) # Probability value from GEV link
#' plot(sort(p))
#' y <- rbinom(n, 1, p)# Draw bernoulli r.v. with prob. p
#' table(y)/n
#'
#' # Fit model with correct GEV parameterization
#' dat <- data.frame(y = y, x1 = x1, x2 = x2, x3 = x3)
#' glm(y ~ x1 + x2 + x3, dat, family = binomial(link = fam),
#'      mustart = rep(mean(y), length(y)))
#'
#' # Compare results to logit, probit, cauchit, and gumbel links (logit, cloglog)
#' glm(y ~ x1 + x2 + x3, dat, family = binomial(link = logit))
#' glm(y ~ x1 + x2 + x3, dat, family = binomial(link = probit))
#' glm(y ~ x1 + x2 + x3, dat, family = binomial(link = cauchit))
#' glm(y ~ x1 + x2 + x3, dat, family = binomial(link = cloglog)) # gev(0)
#'
#' # Fit to GEV link with different shape parameters, choose best model with AIC
#' xi.val <- seq(-1, 0, 0.1)
#' xi.val[which.min(sapply(xi.val, function(xi) glm(y = x1 + x2, dat,
#'                                                  family = binomial(link = gev(xi)),
#'                                                  mustart = rep(mean(y), length(y)))$aic))]
#'
#'
#'

link.gev <- function (xi = 0, rev = F) {

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
