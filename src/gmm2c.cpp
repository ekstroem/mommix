# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

//' Gaussian mixture of regression models with contamination
//' 
//' Fast computation of simple regression slopes for each predictor represented by a column in a matrix
//'
//' Missing values (NA, Inf, NaN) are completely disregarded and pairwise complete cases are used for the analysis.
//' 
//' @param y A vector of outcomes.
//' @param X A design matrix of regressor variables. Must have the same number of rows as the length of y.
//' @param weight A vector of weights for each observation that are used in the estimation.
//' @param maxit the maximum number of iterations to use for the EM algorithm
//' @param tol the tolerance used for determining convergence
//' @param alpha the starting probability for an observation originating from the contamination distribution. Must be strictly between 0 and 1
//' @param mu the starting value for the mean of the contamination distribution. If set to NA (the default) then the mean of the y vector is used.
//' @return A list with the variables: N, K, coefficients, mu, sigma1, sigma2, alpha, iterationsused, and groupprob which contains
//' @author Claus Ekstrom <ekstrom@@sund.ku.dk>
//' @examples
//' x <- rnorm(1000)
//' y <- rnorm(1000, mean=c(x[1:700], rep(0, 300)), sd=rep(c(1,2), times=c(700,300)))
//' mgrwc(y, cbind(rep(1, 1000), x), weight=rep(1, 1000))
//'
//' @export
// [[Rcpp::export]]
List mgrwc (const arma::colvec y,
            const arma::mat X,
	    const arma::colvec weight,
            const int maxit = 400,
	    const double tol = 1e-07,
	    double alpha = 0.3,   // Probability of contamination distribution. Should be in (0,1)
	    double mu=NA_REAL
            ) {
    // inputs
    const int N = y.n_rows ;
    const int K = X.n_cols ;
    int it;
    // containers

    double nll = 0;
    double oldnll = 0;

    // Use the mean of the y's as a starding point for the distribution of mu
    if (! arma::is_finite(mu)) {
      mu = mean(y);
    }
    
    double sigma1 = 0;       // log Variance of contamination group
    double sigma2 = 0;       // log Variance of regression group
    //    double sumzi = 0;       // Unweighted version
    double sumziw = 0;
    arma::colvec beta(K) ;   // Set regression parameters    
    beta = arma::solve(X, y);  // Use full data regression estimates as starting point

    arma::mat Xb = X * beta ;

    arma::colvec resid = y - Xb;
    sigma2 = .5*log(arma::as_scalar(arma::dot(resid, resid)/(N-K)));
    sigma1 = .5*log(arma::sum(square(y - arma::mean(y))) / (N-1));
    
    arma::colvec zi(N) ;
    //    arma::colvec sqrtzi(N) ;    // Unweighted version
    arma::colvec sqrtziw(N) ;

    double weightsum = arma::sum(weight);
      
    // algorithm
    for (it = 0 ; it < maxit ; it++) {      
      // First find classification probabilities
      for (int n = 0 ; n < N ; n++) {
	zi(n) = 1 / (1 + exp(log(1-alpha) - log(alpha) + R::dnorm(y(n), Xb(n,0), exp(sigma2), true) - R::dnorm(y(n), mu, exp(sigma1), true)));
      }
      
      // Add weights here to zi for weighted likelihood
      
      // sqrtzi = sqrt(1-zi);  // Unweighted version
      sqrtziw = sqrt((1-zi)%weight);

      //      sumzi = arma::sum(zi);  // Unweighted version
      sumziw = arma::sum(zi%weight);

      // Check if we are close to the border
      if (alpha < 0.01 || alpha > .99) {
	beta = arma::solve(X, y);
	sigma2 = .5*log(arma::as_scalar(arma::dot(resid, resid)/(N-K)));
	mu = NA_REAL;
	sigma1 = NA_REAL;
	alpha=1.0;
	break;
      }

      // alpha = arma::mean(zi);  // Probability of contamination
      alpha = arma::sum(zi%weight)/weightsum;  // Probability of contamination
      // mu = arma::sum(zi%y) / arma::sum(zi);   // SAVE
      mu = arma::sum(zi%weight%y) / arma::sum(zi%weight);  // Mean of contamination distribution
      // sigma1 = .5*log(arma::sum( zi%square(y-mu)) / (sumzi));
      sigma1 = .5*log(arma::sum( weight%zi%square(y-mu)) / (sumziw));

      // Hvad sker der her?
      // beta = arma::solve(sqrtzi%X.each_col(), sqrtzi%y);
      beta = arma::solve(sqrtziw%X.each_col(), sqrtziw%y);

      // Update prediction 
      Xb = X * beta ;
      resid = y - Xb;
      
      // sigma2 = .5*log (arma::sum(arma::square(sqrtzi%resid))/(N-sumzi));
      sigma2 = .5*log (arma::sum(arma::square(sqrtziw%resid))/(N-sumziw));

      // TO DO
      // If alpha close to border then remove the contamination component
      // If sigma1/2 becomes small then remove the contamination component

      oldnll = nll;

      // Evaluate the log likelihood
      nll = 0;
      for (int n = 0 ; n < N ; n++) {
        nll -=  log(alpha * R::dnorm(y(n), mu, exp(sigma1), false) + (1-alpha)*R::dnorm(y(n), Xb(n,0), exp(sigma2), false));
      }
      
      // Exit if the change in log likelihood is below the tolerance
      if (fabs(nll - oldnll) < tol)
      	break;
    }

    // Now we are done with the main estimation

    // returns
    List ret ;
    ret["N"] = N ;
    ret["K"] = K ;
    ret["coefficients"] = beta ;
    ret["mu"] = mu ;
    ret["sigma1"] = exp(sigma1) ;
    ret["sigma2"] = exp(sigma2) ;
    ret["alpha"] = 1-alpha ;
    ret["iterationsused"] = it ;
    ret["groupprob"] = zi;
    return(ret) ;
}
