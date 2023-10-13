## Aysu Ismayilova s2295782

##Code to efficiently fit p_spline smoothers to x and y data with GCV smoothing parameter using selection basis expansions and penalized regression
#k is chosen generously, so that the approximation can accurately capture a wide range of function shapes and
#complexities, but that means that there is a danger of over-fitting noisy data (i.e. fitting the ϵi component of y not just the f(xi) component). To avoid this we can impose a smoothing penalty to encourage βj s corresponding to
#neighbouring bj (x)s to vary smoothly from one to the next.
#Computing GCV for many trial values of λ would be rather expensive if we directly use fomrulas for GCV but  it turns out that if some matrix decompositions are done before we start searching for
#the optimal λ, then the search can be made very efficient.
##We use  inversion of diagonal matrix  I + λΛ in our computations which  takes only O(k)
##operations for each new λ value, as opposed to the O(k^3) needed by the original expressions(without calculating upfront QR Decomposition and eigen decomposition)



#The function pspline that returns a list containing coefficient, smoothing parameter mu.hat, residual variance, knots vector,gcv, r_squared,
#residual standard deviation,the B-spline and the order of difference to use in the penalty
#Input: the original data x and y obtained from the code above using library(MASS), k---- the number of basis functions to use,
#logsp ---ends of the interval over which we should search for the smoothing parameter in log scale 
#bord---the B-spline order.In our function we use bord =3 that corresponds to cubic 
#pord--- the order of difference that we use in penalty. In our function we use penalty =2
#ngrid ---the number of smoothing parameter values to try. (we should use even spacing of values on the log scale)


pspline <- function(x,y,k, logsp, bord, pord,ngrid){
 
  
  
  lsp = seq(logsp[1], logsp[2], length =ngrid) #sequence of 100 smoothing parameters in log scale
  dk <- diff(range(x))/(k-bord)
  knots <- seq(min(x)-dk*bord, by = dk, length = k+bord+1) #Evaluating the design matrix for the B-splines defined by knots at the values in x.
  X <- splines::splineDesign(knots,x,ord = bord+1, outer.o = TRUE)
  
  
  
  
  


  
  mu.hat <- list() #initialise fitted values (mu hat)
  b.hat <- list() #initialise penalised least squares (beta hat)
  edf <- gcv  <- res_var<- lsp*0 ## initialise effective degree of freedom, generalized cross-validation,residual variance,log smoothing parameter 
  for (i in 1:length(lsp)){ ##for loop over all log smoothing parameters 
    
    fr <- fit.pspline(y,x, X,XX,k,knots,gcv, DD, D, exp(lsp[i]),qrx, lambda, U, bord, pord) ## fit
    gcv[i] <- fr$gcv 
    edf[i] <- fr$edf
    mu.hat[[i]] <- fr$mu.hat 
  
    res_var[i] <- fr$res_var 
  }
  i.opt <- max(which(gcv==min(gcv))) ##optimal index i can be obtained by finding minimum gcv value
                                    #An effective way to choose λ is to find the value of λ minimizing the generalized cross validation criterion
  

  fit.pspline(y,x, X,XX,k,knots,gcv,DD, D,exp(lsp[i.opt]), qrx, lambda, U, bord, pord) 
}

##function pspline takes as an input vector y,matrix X, cross product of X (XX), coefficient k,cross product of the matrix D,
##matrix D, smoothing parameter, QR factorisation of the matrix X, eigenvalues, eigenvectors,he B-spline order and the order of difference that we use in penalty
fit.pspline <- function(y,x,X,XX,k,knots, gcv,DD,D,  sp, qrx, lambda, U, bord, pord){
  D <- diff(diag(k), differences =pord)##matrix D which is a k-2 by k matrix of zeroes, except for Di,i = Di,i+2 = 1 and Di,i+1 = -2,for i=1,....,k-2 
  DD = crossprod(D)##cross product of D (D^T*D)
  
  XX <- crossprod(X)##cross product of X (X^T*X)
  qrx <- qr(X) ##QR factorisation of X
  Q <- qr.Q(qrx)## The matrix Q
  R <- qr.R(qrx)## The matrix R
  
  n <- length(y) ##number of y_i data
 

  
  eigen_decomposition = eigen(solve(t(R))%*%t(D)%*%D%*%(solve(R))) ## eigen decomposition 
  U =eigen_decomposition$vectors ## eigenvectors
 
  lambda = diag(eigen_decomposition$values)## eigenvalues
  b.hat = backsolve(qr.R(qrx),U%*%solve(diag(k)+lambda*sp)%*%t(U)%*% qr.qty(qrx,y)[1:k]) ##we can find beta hat by using this formula that involves QR decomposition of X and then eigen decompostion
  
  mu.hat = X%*%b.hat ##fitted values 
  kappa <- sum(diag(solve(diag(k)+sp*lambda))) ##effective degrees of freedom
  res_var <- sum((y-mu.hat)^2)/(n-kappa) ##residual variance
  V <- solve(XX+sp*DD)*res_var ##covariance matrix V
  

  
  gcv <- res_var/(n-kappa) #finding the generalised cross validation
  
  y_bar = mean(y) #mean value of y
  vec = c() #initialise vector 
  
  for (s in 1:n){     #for loop from 1 to n
    a= (y[s]-y_bar)^2 #finding (y_i-y_bar)^2
    vec <- c(vec, a) #appending the value to the empty vector vec
    
    
  }
  sum_a = sum(vec) #finding sum of all elements
  
  #rsquared is calculated using the following formula
  r_squared = 1-((n-1)*res_var)/sum_a 
  residual_std_dev = sqrt(res_var)
  d <- list(X=X,V = V,x= x,y=y,knots = knots,D = D, coef = k,mu.hat = mu.hat,b.hat = b.hat, gcv = gcv, r_squared = r_squared,edf = kappa, res_var = res_var,residual_std_dev =residual_std_dev, bord = bord, pord = pord)
  class(d) <- 'pspline' #assigning the class pspline for the list d 
  invisible(d)
}



m <- pspline(x, y, k = 5, logsp =c(-5,5),bord = 3, pord =2, ngrid = 100)



#creating function print.pspline which is a method function that reports some details of the model fit m

print.pspline <- function(m){ 
  cat('Order', m$bord,  'p_spline with order', m$pord, 'penalty',"\n")
  cat('Effective degrees of freedom',m$edf, " ")
  cat('Coefficients:',m$coef, "\n")
  cat('residual std dev:',m$residual_std_dev," " )
  cat('r-squared:',m$r_squared, " ")
  cat('GCV',m$gcv, "\n")
  li <- list(m$gcv, m$edf, m$r_squared)
  invisible(li) ##the function silently returns list consisting of generalised cross validation, effective degrees of freedom
                ##and r_squared error
}


print.pspline(m) 




#creating function print.pspline which takes as an input model fit m, new data x and makes predictions for new x values within the range of 
##of the original data x and standard errors 
##if se = TRUE the function returns 2 item list with predicted values in fit item and corresponding standard errors in se item
##if se = FALSE the function returns vector of predictions
predict.pspline<- function(m,x, se){
  k = m$coef ##extracting the coefficient from the list 
  bord = m$bord ##extracting the B-spline order from the list 
  
  Xp <- splines::splineDesign(m$knots,x,ord = bord+1, outer.o = TRUE) ##Creating new matrix Xp for the new data with the original knot positions
  if(se == FALSE){
    
    y_pred = Xp%*%m$b.hat ##vector of predictions
    return(y_pred)
  }
  if(se == TRUE){
    
    y_pred = Xp%*%m$b.hat ##vector of predictions
    se <- rowSums(Xp*(Xp%*%m$V))^0.5 ##vector of standard errors
    
    list_pred_se = list(m = y_pred, se= se) ##creating a 2 item list of predicted values and standard errors
    
    return(list_pred_se)
    
  }
}

predict.pspline(m, newx, se = TRUE)


#creating plot.pspline function that takes a model fit m and produces 3 plots: 
#1) plot of the original x,y data with the estimated smooth function along with approximate 95% credible intervals for the smooth
#2) plot of the model residuals against fitted values
#3) qqplot of the residuals
plot.pspline <- function(m){
  standard_error <-rowSums(m$X*(m$X%*%m$V))^0.5 #formula for calculating the standard error. Here V is the covariance matrix
  lower_b = m$mu.hat-1.96*(standard_error) #lower bound of the confidence interval is obtained by using f_estimated - 1.96*standard errors
  upper_b = m$mu.hat +1.96*(standard_error) #upper bound of the confidence interval is obtained by using f_estimated + 1.96*standard errors
  
  
  plot(x, y, xlab="data x", ylab="data y", main="Plot of the original x and y data")
  lines(x, lty = 1,m$mu.hat) #adding lines to the plot
  lines(x, lower_b, lty =2, col = 2) #adding dashed lines to the plot that corressponds to the lower bound of the confidence interval
  lines(x, upper_b, lty = 2, col = 2 )#adding dashed lines to the plot that corressponds to the upper bound of the confidence interval
  legend(x = "topleft", box.col = "brown",
        box.lwd = 2 , title="Plot of the original x and
                          y data along with 
                          95% confidence intervals", 
         legend=c("upper bound", "smoothing function", "lower bound"), 
        lty = c(2,1,2), col = c("red","black", "red"), text.font = 1, cex = 0.7)
  
  res <-  c() #creating an empty vector for residuals
  
  for (elem in 1:length(y)){#loop over the n which is length of the original data y
    eps <-  (y[elem]-m$mu.hat[elem]) #formula for finding residuals
    res <- c(res, eps) #append residuals to the empty list
  } 
  plot(m$mu.hat, res, xlab="Fitted values", ylab = 'Residuals')
  
  qqnorm(res, main = "Normal Q-Q Plot of Residuals",
         xlab = "x", ylab = "y")
  list1 = list(ll = lower_b, ul = upper_b, x =x) 
  invisible(list1) #silently returns a list containing lower bound, upper bound and vector x
}
plot.pspline(m)
