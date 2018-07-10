#' Mixture modeling using moments
#' 
#' Fit a mixture model using moments
#' 
#' @param data a data frame or tibble
#' @param formula a formula model formula for the relationship
#' @param weights a string that defines the weights used for the estimation using weighted regression. Possibilities are "equal" (the default) for unweighted regression estimation and "square" for squared linear predictor estimation.
#' @param maxiter integer
#' @return Returns a list with the following elements:
#' 
#' * coefficients The estimates of the parameters for 
#' * alpha
#' * mu
#' * covariance
#' * call
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
#' moment_mixture_model(DF, y ~ x + x2)
#' 
#' @export
moment_mixture_model <- function(data, formula, weights=c("equal", "square"), maxiter=5) {

    weightmethod <- match.arg(weights)
    call <- match.call(expand.dots = TRUE)

    ## Ensure that the formula includes an intercept
    iform <- update.formula(formula, . ~ . +1)

    ## Initial main model fit
    mainmodel <- lm(formula=iform, data=data)

    ## Extract the linear predictor for all but the intercept (scaled by the mixture probability)
    lambda1 <- coef(mainmodel)[-1]                               # Remove intercept
    eta <- predict(mainmodel) - coef(mainmodel)["(Intercept)"]   #

    ## Some scaling would probably by good to prevent numerical problems with squaring !    
    Y <- model.frame(mainmodel)[,1]
    X <- model.frame(mainmodel)[,-1]
    ## scalingfactor <- sqrt(max(abs(c(Y, eta))))
    scalingfactor <- 1
    Y2 <- (Y/scalingfactor)^2        ## We scale it to ease the computations

    ## Define the weights
    weight.lin <- rep(1, length(eta))
    weight.sqr <- rep(1, length(eta))
    ## Estimate with weights
    if (weightmethod=="square") {
        weight.lin <- 1/(1+eta^2)
        weight.sqr <- 1/(1+eta^4)
    }

    ## Trying to stabilize something
    weight.lin <- weight.lin/sum(weight.lin)
    weight.sqr <- weight.sqr/sum(weight.sqr)

    ## weight.lin needs to be part of the data data frame
    data$weight.lin <- weight.lin
    
    ## Set maxiter to 1 if weight is equal
    if (weightmethod=="equal") {    
        maxiter <- 1
    }
    ## Iterate between
    for (i in seq(maxiter)) {
        mainmodel <- lm(formula=iform, data=data, weights=weight.lin)
        eta <- predict(mainmodel) - coef(mainmodel)["(Intercept)"]   #
        lambda1 <- coef(mainmodel)[-1]                               # Remove intercept

        ## Update the weights based on the weighted regression
        if (weightmethod=="square") {
            weight.lin <- 1/(1+eta^2)
            weight.sqr <- 1/(1+eta^4)
            weight.lin <- weight.lin/sum(weight.lin)
            data$weight.lin <- weight.lin
            weight.sqr <- weight.sqr/sum(weight.sqr)
        }
                        
        ## Old school with OLS
        sqmodel <- lm(Y2 ~ I(eta/scalingfactor) + I((eta/scalingfactor)^2), weights=weight.sqr)
        lambda2 <- unname(coef(sqmodel)[2])
        lambda3 <- unname(coef(sqmodel)[3]) 

        ## Compute the parts that are used for a's and b's
        ## Better to do this in steps as we need the derivatives as well
        ## Note the weights have been standardized above
        zbar1 <- sum(weight.sqr*eta)
        zbar2 <- sum(weight.sqr*eta^2)
        zbar3 <- sum(weight.sqr*eta^3)
        zbar4 <- sum(weight.sqr*eta^4)
#        dweights <- (-4)*(eta^3)*x*weights*weights
        
        
        ## Derivatives needed for variance
##        dzbar1 <- mean(x*weights+dweights*eta)/wbar-zbar1*mean(dweights)/(wbar)
        
        
        a1n <- sum(weight.sqr*eta^4) - (sum(weight.sqr*eta^2))^2
        a2n <- sum(weight.sqr*eta^3) - (sum(weight.sqr*eta^2))*(sum(weight.sqr*eta))
        a3n <- sum(weight.sqr*eta^2) - (sum(weight.sqr*eta))^2
        b1n <- sum(weight.sqr*Y2*eta)   - (sum(weight.sqr*eta))*(sum(weight.sqr*Y2))
        b2n <- sum(weight.sqr*Y2*eta^2) - (sum(weight.sqr*eta^2))*(sum(weight.sqr*Y2))
        
        ## Estimate p
        ## Should do a check for coefficient not too close to zero
        phat <- unname(1/lambda3)
        beta <- lambda1*lambda3
    }

    ## Now estimate the covariance
    ## Compute the influence functions including epsilon, and derivatives

 ####   epsiloni <- diag(solve(var(X))) * (X-colMeans(X))*resid(mainmodel)

    # Derivatives used for deriviation of a's and b's
#    dzbar1 <- 

    ##ksi = (lambda3*diag(length(lambda1)))

    ## Compute the derivatives of a's and b's
#    a1n<-zbar4-zbar2^2
#    da1n<-dzbar4-2*zbar2*dzbar2
#    a2n<-zbar3-zbar1*zbar2
#    da2n<-dzbar3-dzbar1*zbar2-zbar1*dzbar2
#    a3n<-zbar2-zbar1^2
#    da3n<-dzbar2-2*zbar1*dzbar1
#    b1n<-y2z1bar-y2bar*zbar1
#    db1n<-dy2z1bar-dy2bar*zbar1-y2bar*dzbar1
#    b2n<-y2z2bar-y2bar*zbar2
#    db2n<-dy2z2bar-dy2bar*zbar2-y2bar*dzbar2
#    g1n<-a1n*b1n-a2n*b2n   ## Taeller for lambda2                                                                                             
#    dg1n<-da1n*b1n+a1n*db1n-da2n*b2n-a2n*db2n
#    g2n<-a3n*b2n-a2n*b1n   ## Taeller for lambda3                                                                                             
#    dg2n<-da3n*b2n+a3n*db2n-da2n*b1n-a2n*db1n
#    g3n<-a1n*a3n-a2n*a2n   ## Naevner                                                                                                         
#    dg3n<-da1n*a3n+a1n*da3n-2*a2n*da2n
    
    

    ## Fix boundary problems. Note that this crude fix my introduce some discontinuities in the results
    ## This will ensure that we just return the results from a model with one component.
    if (phat>1) {        
        phat <- 1
    }
    if (phat<0) {
        phat <- 0
    }
        
    RET <- list(coefficients = c(`(Intercept)`=lambda2/2, ifelse(is.infinite(beta), rep(NA, length(beta)), beta)),
                alpha=phat,
##                mu1 = lambda2/2,
                mu2 = ifelse(phat==1, NA, (coef(mainmodel)[1] - lambda2/(2*lambda3))/(1-phat)),
                covariance=diag(length(beta)),
                method="moment",
                call=call
                )
    class(RET) <- "mommix"
    RET
}
