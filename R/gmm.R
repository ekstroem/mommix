#' Fit 2-component Gaussian mixture model
#' 
#' Fit a 2-component Gaussian mixture model.
#' 
#' @param data a data frame or tibble
#' @param formula a formula model formula for the relationship
#' @param B a positive integer giving the number of bootstrap samples used to obtain the estimates of the variance associated with the parameter estimates.
#' @param se.method The method to compute the standard errors. Defaults to "bootstrap"
#' but "wlbs" (weighted likelihood bootstrap) is probably a more robust bet.
#' @param mu A starting value for the mu component. Defaults to NA (which means that the starting value is calculated from the average value of the outcomes) 
#' @param mufixed A boolean stating whether the mu component should be fixed at the given value. Defaults to FALSE.
#' @param maxit The maximum number of iterations used for computing the estimates from the gaussian mixture model. Defaults to 200.
#' @param ... Additional arguments passed on to the fitting process.
#' @return Returns a list with the following elements: ...
#' @author Claus Ekstrom \email{ekstrom@@sund.ku.dk} and Christian Pipper \email{pipper@@sund.ku.dk}
#' @references Unpublished manuscript by authors and Newton and Raftery (1994): "Approximate Bayesian Inference with the Weighted Likelihood Boostrap". JRSS-
#' @keywords manip
#' @examples
#'
#' p <- 0.7  #
#' N <- 100
#' x <- rnorm(N)
#' x2 <- rnorm(N)
#' y <- rnorm(N, mean=x)
#' y[1:(N*p)] <- rnorm(N*p, mean=0, sd=.25)
#' 
#' DF <- data.frame(x, x2, y)
#'
#' gaussian_mixture_model(DF, y ~ x + x2)
#' 
#' @export
gaussian_mixture_model <- function(data, formula, B=200, se.method=c("bootstrap", "wlbs"), mu=NA, mufixed=FALSE, maxit=200, ...) {

    ## Figure out the default method
    se.method <- match.arg(se.method)
    
    ## We run a fit for a lm model to obtain
    ## starting values and to extract the variables
    ## from the relevant environment
    call <- match.call(expand.dots = TRUE)
    origcall <- call    
    lmcall <- origcall
    lmcall[[1]] <- as.name("lm")
    lmcall$maxit <- lmcall$tol <- lmcall$B <- lmcall$se.method <- lmcall$mu <- lmcall$mufixed <- NULL
    lmFit <- eval(lmcall, parent.frame())

    args <- list(...)
#    maxit <- ifelse("maxit" %in% names(args), args$maxit, 200)
    tol <- ifelse("tol" %in% names(args), args$tol, 1e-07)

    ## Jump a hoop by getting the evaluated formulas from lm
    Y <- model.response(model.frame(lmFit), "numeric")
    X <- model.matrix(lmFit)

    ## Should probablky update the CPP function s.t. starting parameters could be given
    fitresult <- mgrwc(y=Y, X=X, maxit=maxit, tol=tol, weight=rep(1,length(Y)), mu=mu, mufixed=mufixed)

    ## Now compute bootstrap SE
    if (se.method == "wlbs") {
        bootresult <- sapply(1:B, function(x) {
            indexsample <- sample(1:length(Y), replace=TRUE)
            mgrwc(y=Y[indexsample], X=X[indexsample,], maxit=maxit, tol=tol, weight=rep(1,length(Y)), alpha=1-(fitresult$alpha+0.05)/1.1, mu=fitresult$mu)
        } )
    } else {
        bootresult <- sapply(1:B, function(x) {
            weight <- rexp(length(Y))
            weight <- length(Y)*weight/sum(weight)
            mgrwc(y=Y, X=X, maxit=maxit, tol=tol, weight=weight, alpha=1-(fitresult$alpha+0.05)/1.1, mu=fitresult$mu)
        } )
    }

    se.alpha <- sd(unlist(bootresult[7,]))

    Bestimates <- matrix(unlist(bootresult[3,]), ncol=B)
    Bmean <- rowMeans(Bestimates)
    tmpB <- Bestimates - Bmean
    Bcovariance <- matrix(B/(B-1)*rowMeans(apply(tmpB, 2, tcrossprod)), ncol=length(fitresult$coefficients))

    RES <- 
    list(coefficients=setNames(as.numeric(fitresult$coefficients), colnames(X)),
         alpha=fitresult$alpha,
         mu2=fitresult$mu,
         covariance=Bcovariance,
         se.beta=sqrt(diag(Bcovariance)),
         se.alpha=se.alpha,
         convergence=c(iterations=fitresult$iterationsused),
         method="gaussian",
         call=call)
    
    class(RES) <- "mommix"

    RES
}
