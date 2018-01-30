#' Functions for mixture estimation based on moments
#'
#' Functions for mixture estimation based on moments.
#'
#' \tabular{ll}{ Package: \tab mommix\cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2012-03-29\cr License: \tab GPL-2\cr } % ~~ An overview of
#' how to use the package, including the most important ~~ % ~~ functions ~~
#'
#' @name mommix
#' @aliases mommix-package mommix
#' @docType package
#' @useDynLib mommix, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @author Claus Thorn Ekstrom \email{ekstrom@@sund.ku.dk}\cr Maintainer: Claus Thorn Ekstrom
#' \email{ekstrom@@sund.ku.dk}
#' @references Some paper
#' @import utils stats graphics
#' @keywords package
NULL

#' Portuguese wine quality 
#'
#' The datasets contains red and white variants of the Portuguese "Vinho Verde" wine.
#' These datasets can be viewed as classification or regression tasks. The classes are ordered and not balanced
#' (e.g. there are munch more normal wines than excellent or poor ones). Outlier detection algorithms could be used
#' to detect the few excellent or poor wines. Also, we are not sure if all input variables are relevant.
#'
#' Due to privacy and logistic issues, only physicochemical (inputs) and sensory (the output) variables
#' are available (e.g. there is no data about grape types, wine brand, wine selling price, etc.). 
#' 
#' @name wine
#' @docType data
#' @format A data frame with 6497 observations on the following 13 variables.
#' \describe{ \item{fixed.acidity}{a numeric vector}
#' \item{volatile.acidity}{a numeric vector}
#' \item{citric.acid}{a numeric vector}
#' \item{residual.sugar}{a numeric vector}
#' \item{chlorides}{a numeric vector}
#' \item{free.sulfur.dioxide}{a numeric vector}
#' \item{total.sulfur.dioxide}{a numeric vector}
#' \item{density}{a numeric vector}
#' \item{pH}{a numeric vector}
#' \item{sulphates}{a numeric vector}
#' \item{alcohol}{a numeric vector}
#' \item{quality}{a numeric vector}
#' \item{colour}{a character vector} }
#' @references P. Cortez, A. Cerdeira, F. Almeida, T. Matos and J. Reis. Modeling wine preferences by data mining from physicochemical properties. In Decision Support Systems, Elsevier, 47(4):547-553, 2009. 
#' @source Data was obtained from the UCI Machine Learning Repository, https://archive.ics.uci.edu/ml/datasets/wine+quality
#' @keywords datasets
#' @examples
#' 
#' data(wine)
#' 
NULL
