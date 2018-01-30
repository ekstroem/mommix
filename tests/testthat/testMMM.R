context("test results for MMM")

set.seed(1)
p <- 0.7                                                                                                            
N <- 100                                                                                                                                
x <- rnorm(N)                                                                                                                           
x2 <- rnorm(N)                                                                                                                          
y <- rnorm(N, mean=x)                                                                                                                   
y[1:(N*p)] <- rnorm(N*p, mean=0, sd=.25)                                                                                                
DF <- data.frame(x, x2, y)

res <- moment_mixture_model(DF, y ~ x + x2)
res2 <- moment_mixture_model(DF, y ~ x + x2, weight="square")

test_that("computations are consistant", {
    expect_equivalent(coef(res), c(0.1654199, 0.8233733, -0.2667224), tolerance=1e-4)
    expect_equivalent(coef(res2), c(0.1495042, 0.8259178, -0.2469152), tolerance=1e-4)
    expect_equivalent(res$alpha,0.3041905, tolerance=1e-4)
    expect_equivalent(res2$alpha,0.3162248, tolerance=1e-4)    
} )
