#' CBP model fit
#' 
#' Fit a mixture model using moments. Only one component prototype
#' 
#' @param y a data frame or tibble
#' @param x asdasd
#' @param probud probability of group 1
#' @param alpha mu1
#' @param beta effect
#' @return Returns a list with the following elements:.
#' @author Christian Pipper \email{pipper@@sund.ku.dk}
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
#' cbpcode(y, x)
#' 
#' @export
cbpcode <- function(y, x, probud=0.5, alpha=1, beta=1) {

    probdud <- 1-probud
    truebeta <- beta

    Nsim<-1
    Nsamp <- length(y)
    
# probdud<-0.5
#alpha<-1
#beta<-1
# sigma.signal<-1
#sigma.duds<-1
# varx<-1
est.gamma1<-rep(0,Nsim)
est.betabest<-rep(0,Nsim)
var.betabest<-rep(0,Nsim)
est.gamma2<-rep(0,Nsim)
est.gamma3<-rep(0,Nsim)
est.beta<-rep(0,Nsim)
var.gamma1<-rep(0,Nsim)
var.gamma2<-rep(0,Nsim)
var.gamma3<-rep(0,Nsim)
var.beta<-rep(0,Nsim)


i <- 1
    ## First model
    fit1ini<-lm(y~x)


  
#  est.gamma1alt[i]<-coef(fit1ini)[2]
    ztildeini<-predict(fit1ini)-coef(fit1ini)[1]
    mygivenx<-predict(fit1ini)
    vx<-var(x)
    epsiloninii<-(vx^-1)*(x-mean(x))*(y-mygivenx)   ### <- hvorfor mean(x) og var?
#    print(epsiloninii)
#    print((vx^-1)*(x-mean(x))*(y-mygivenx))
    #  var.gamma1alt[i]<-mean(epsiloninii*epsiloninii)/Nsamp
    


  ## ztildeini er eta

  
  weights1<-1/(1+ztildeini^2)
  fit1<-lm(y~x,weights=weights1)
  est.gamma1[i]<-coef(fit1)[2]
  ztilde<-predict(fit1)-coef(fit1)[1]
  
  weights<-1/(1+ztilde^4)
  fit2<-lm(I(y^2)~ztilde+I(ztilde^2),weights=weights)
  est.gamma2[i]<-coef(fit2)[2]    ## 
    est.gamma3[i]<-coef(fit2)[3]    ##

    ## Added CE
    if (is.na(est.gamma3[i]))
        est.gamma3[i] <- 0

    est.beta[i]<-est.gamma3[i]*est.gamma1[i]


## Good until here
  
  dweights1<--2*x*ztildeini*weights1*weights1
  xbar<-sum(weights1*x)/sum(weights1)
  x2bar<-sum(weights1*x*x)/sum(weights1)
  ybar<-sum(weights1*y)/sum(weights1)
  yxbar<-sum(weights1*y*x)/sum(weights1)
  dxbar<-(sum(dweights1*x)*sum(weights1)-sum(dweights1)*sum(weights1*x))/(sum(weights1)*sum(weights1))
  dx2bar<-(sum(dweights1*x*x)*sum(weights1)-sum(dweights1)*sum(weights1*x*x))/(sum(weights1)*sum(weights1))
  dybar<-(sum(dweights1*y)*sum(weights1)-sum(dweights1)*sum(weights1*y))/(sum(weights1)*sum(weights1))
  dyxbar<-(sum(dweights1*y*x)*sum(weights1)-sum(dweights1)*sum(weights1*y*x))/(sum(weights1)*sum(weights1))
  h1n<-yxbar-xbar*ybar
  dh1n<-dyxbar-dybar*xbar-ybar*dxbar
  h2n<-x2bar-xbar*xbar
  dh2n<-dx2bar-2*dxbar*xbar
  dlambda1n<-(dh1n*h2n-h1n*dh2n)/(h2n*h2n)
  vbar<-mean(weights1)  
  epsiloni<-dlambda1n*epsiloninii+(1/(h2n*vbar))*(x-xbar)*(y-predict(fit1))*weights1
  wbar<-mean(weights)
  dweights<-(-4)*(ztilde^3)*x*weights*weights
  
  zbar1<-sum(ztilde*weights)/sum(weights)
  dzbar1<-mean(x*weights+dweights*ztilde)/wbar-zbar1*mean(dweights)/(wbar)
  zbar2<-sum(ztilde*ztilde*weights)/sum(weights)
  dzbar2<-mean(2*x*ztilde*weights+ztilde*ztilde*dweights)/wbar-zbar2*mean(dweights)/(wbar)
  zbar3<-sum(ztilde*ztilde*ztilde*weights)/sum(weights)
  dzbar3<-mean(3*x*ztilde*ztilde*weights+ztilde*ztilde*ztilde*dweights)/wbar-zbar3*mean(dweights)/(wbar)
  zbar4<-sum(ztilde*ztilde*ztilde*ztilde*weights)/sum(weights)
  dzbar4<-mean(4*x*ztilde*ztilde*ztilde*weights+ztilde*ztilde*ztilde*ztilde*dweights)/wbar-zbar4*mean(dweights)/(wbar)

  y2z1bar<-sum(weights*ztilde*y*y)/sum(weights)
  dy2z1bar<-mean(x*y*y*weights+ztilde*y*y*dweights)/wbar-y2z1bar*mean(dweights)/(wbar)
  y2z2bar<-sum(weights*ztilde*ztilde*y*y)/sum(weights)
  dy2z2bar<-mean(2*x*ztilde*y*y*weights+ztilde*ztilde*y*y*dweights)/wbar-y2z2bar*mean(dweights)/(wbar)
  y2bar<-sum(weights*y*y)/sum(weights)
  dy2bar<-mean(y*y*dweights)/wbar-y2bar*mean(dweights)/(wbar)

  
  a1n<-zbar4-zbar2^2
  da1n<-dzbar4-2*zbar2*dzbar2
  a2n<-zbar3-zbar1*zbar2
  da2n<-dzbar3-dzbar1*zbar2-zbar1*dzbar2
  a3n<-zbar2-zbar1^2
  da3n<-dzbar2-2*zbar1*dzbar1
  b1n<-y2z1bar-y2bar*zbar1
  db1n<-dy2z1bar-dy2bar*zbar1-y2bar*dzbar1
  b2n<-y2z2bar-y2bar*zbar2
  db2n<-dy2z2bar-dy2bar*zbar2-y2bar*dzbar2
  g1n<-a1n*b1n-a2n*b2n   ## Tæller for lambda2
  dg1n<-da1n*b1n+a1n*db1n-da2n*b2n-a2n*db2n
  g2n<-a3n*b2n-a2n*b1n   ## Tæller for lambda3
  dg2n<-da3n*b2n+a3n*db2n-da2n*b1n-a2n*db1n
  g3n<-a1n*a3n-a2n*a2n   ## Nævner
  dg3n<-da1n*a3n+a1n*da3n-2*a2n*da2n


  
  dlambda2n<-(dg1n*g3n-g1n*dg3n)/(g3n*g3n)
  dlambda3n<-(dg2n*g3n-g2n*dg3n)/(g3n*g3n)
  f1i<-weights*(ztilde-zbar1)*(y*y-predict(fit2))
  f2i<-weights*(ztilde^2-zbar2)*(y*y-predict(fit2))
  infl.lambda2<-dlambda2n*epsiloni+((g3n*mean(weights))^-1)*(a1n*f1i-a2n*f2i)
  infl.lambda3<-dlambda3n*epsiloni+((g3n*mean(weights))^-1)*(a3n*f2i-a2n*f1i)
  var.gamma1[i]<-(1/Nsamp)*mean(epsiloni*epsiloni)
  var.gamma2[i]<-(1/Nsamp)*mean(infl.lambda2*infl.lambda2)
  var.gamma3[i]<-(1/Nsamp)*mean(infl.lambda3*infl.lambda3)
  infl.beta<-est.gamma3[i]*epsiloni+est.gamma1[i]*infl.lambda3  ## Ksi - næsten?
  var.beta[i]<-(1/Nsamp)*mean(infl.beta*infl.beta)


    
#true gamma1
gamma1true<-beta*probdud
#gamma1true

#empirical performance
#mean(est.gamma1)
#sd(est.gamma1)
#mean(sqrt(var.gamma1))
#95% CP
#lower<-est.gamma1-1.96*sqrt(var.gamma1)
#upper<-est.gamma1+1.96*sqrt(var.gamma1)
#1-(length(lower[lower>gamma1true])+length(upper[upper<gamma1true]))/Nsim

#true gamma2
gamma2true<-2*alpha
#gamma2true

#empirical performance
#mean(est.gamma2)
#sd(est.gamma2)
#mean(sqrt(var.gamma2))
#95% CP
#lower<-est.gamma2-1.96*sqrt(var.gamma2)
#upper<-est.gamma2+1.96*sqrt(var.gamma2)
#1-(length(lower[lower>gamma2true])+length(upper[upper<gamma2true]))/Nsim


#true gamma3
gamma3true<-1/probdud
gamma3true

#empirical performance
#mean(est.gamma3)
#sd(est.gamma3)
#mean(sqrt(var.gamma3))
#95% CP
#lower<-est.gamma3-1.96*sqrt(var.gamma3)
#upper<-est.gamma3+1.96*sqrt(var.gamma3)
#1-(length(lower[lower>gamma3true])+length(upper[upper<gamma3true]))/Nsim


#true beta
#beta

#empirical performance
#mean(est.beta)
#sd(est.beta)
#mean(sqrt(var.beta))
#relative efficiency
#var(est.beta)/((1/Nsamp)*(1/probdud)*(1/varx)*sigma.signal*sigma.signal)



#95% CP
#lower<-est.beta-1.96*sqrt(var.beta)
#upper<-est.beta+1.96*sqrt(var.beta)
#1-(length(lower[lower>beta])+length(upper[upper<beta]))/Nsim

#par(mfrow=c(1,2))
#plot(fit1ini,which=1:2)
#plot(y~x,col=p+1)
#abline(0,0,lwd=2,col=1)
#abline(alpha,beta,lwd=2,col=2)
#abline(lm(y~x),col="blue",lwd=2)#

#plot(est.beta~est.betabest,xlim=c(0,2),ylim=c(0,2))
#abline(h=1,lwd=2)
                                        #abline(v=1,lwd=2)
    list(beta=est.beta[i],
         gamma1=est.gamma1[i],
         gamma2=est.gamma2[i],
         gamma3=est.gamma3[i],
         segamma1=sqrt(var.gamma1[i]),
         segamma2=sqrt(var.gamma2[i]),
         segamma3=sqrt(var.gamma3[i]),
         sebeta=sqrt(var.beta[i]),
         pinv=est.gamma3[i],
         sepinv=sqrt(var.gamma3[i]),
         gamma1true=gamma1true,
         gamma2true=gamma2true,
         gamma3true=gamma3true,
         betatrue=beta,
         pinvtrue=probud
         )
    
}
