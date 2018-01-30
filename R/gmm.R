#' Fit 2-component Gaussian mixture model
#' 
#' Fit a 2-component Gaussian mixture model.
#' 
#' @param data a data frame or tibble
#' @param formula a formula model formula for the relationship
#' @param ... Additional arguments passed on to the fitting process.
#' @return Returns a list with the following elements: ...
#' @author Claus Ekstrom \email{ekstrom@@sund.ku.dk} and Christian Pipper \email{pipper@@sund.ku.dk}
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
gaussian_mixture_model <- function(data, formula, ...) {

    ## We run a fit for a lm model to obtain
    ## starting values and to extract the variables
    ## from the relevant environment
    call <- match.call(expand.dots = TRUE)
    origcall <- call    
    lmcall <- origcall
    lmcall[[1]] <- as.name("lm")
    lmcall$maxit <- lmcall$tol <- NULL
    lmFit <- eval(lmcall, parent.frame())

    args <- list(...)
    maxit <- ifelse("maxit" %in% names(args), args$maxit, 200)
    tol <- ifelse("tol" %in% names(args), args$tol, 1e-07)

    ## Jump a hoop by getting the evaluated formulas from lm
    Y <- model.response(model.frame(lmFit), "numeric")
    X <- model.matrix(lmFit)

    ## Fix ... saa parameters can be passed

    fitresult <- mgrwc (y=Y, X=X, maxit=maxit, tol=tol)

    list(coefficients=setNames(as.numeric(fitresult$coefficients), colnames(X)),
         alpha=fitresult$alpha,
         mu2=fitresult$mu,
         covariance=diag(length(fitresult$coefficients)),
         convergence=c(iterations=fitresult$iterationsused),
         call=call)

}
