distrHsDiscr <- function(y2, eta2, sigma2, sigma2.st, nu, nu.st, margin2, naive = FALSE, y2m){

p2 <- derp2.dersigma.st <- derp2.dereta2 <- der2p2.dereta2eta2 <- der2p2.dersigma2.st2 <- der2p2.dereta2dersigma2.st <- indx <- 1

sigma <- sigma2


der2pdf2.dereta2dernu.st    = 1
der2pdf2.sigma2.st2dernu.st = 1
derpdf2.dernu.st            = 1
der2pdf2.dernu.st2          = 1
derp2.nu.st                 = 1
der2p2.dernu.st2            = 1
der2p2.dereta2dernu.st      = 1
der2p2.dersigma2.stdernu.st = 1

cont1par <- c("PO","ZTP")
cont2par <- c("NBI","NBII","NBIa","NBIIa","PIG","PO","ZTP","DGP","DGPII")
cont3par <- c("DEL","SICHEL")

# library(Deriv); library(numDeriv)
# expr <- expression(  )
# Simplify( D(D(expr, "mu2"),"sigma2") )
# func0 <- function(mu2){   }
# grad(func0 , mu2)
#func <- function(x) pNBI(y2, mu = x[1], sigma = sqrt(x[2])) 
#hessian(func, c(mu2,sigma2))
############################################################################
# remember that eta2 will have to disappear if we change default link on mu2
# this only applies to cases in which mu2 must be positive
# otherwise things are fine
############################################################################

#######################################################################
# Define generic numerical derivative functions
#######################################################################



if(margin2 %in% c("PIG","NBI","NBII","NBIa","NBIIa")){

derpdf2.dermu2FUNC2p <- function(func, y2, mu2, sigma) numgh(func, mu2)
derpdf2.sigma2FUNC2p <- function(func, y2, mu2, sigma) numgh(func, sigma)  

der2pdf2.mu2dersigma2FUNC2p <- function(func, y2, mu2, sigma) numch(func, mu2, sigma)

}



#######################################################################

if(margin2 == "PO"){

mu2 <- exp(eta2)
dermu2.dereta2 <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2)

dersigma2.dersigma2.st  <- 0  
dersigma2.dersigma2.st2 <- 0


if(max(y2) > 170){

prec <- pmax(53, getPrec(mu2), getPrec(y2))
        
mu2 <- mpfr(mu2, prec)
y2  <- mpfr( y2, prec)        
        
}        
        
        
# exp(y2*log(mu2) - mu2 - log(gamma(y2 + 1)))
# should be more stable but looks the same
# from a few experiments
        
        
pdf2 <-  as.numeric( (exp(-mu2)*mu2^y2)/factorial(y2) )   

derpdf2.dermu2FUNCpo <- function(y2, mu2) exp(-mu2) * (mu2^(y2 - 1) * y2 - mu2^y2)/factorial(y2) 
derpdf2.dermu2       <- as.numeric( derpdf2.dermu2FUNCpo(y2, mu2) )
    
der2pdf2.dermu2FUNCpo <- function(y2, mu2) exp(-mu2) * (mu2^y2 + y2 * (mu2^(y2 - 2) * (y2 - 1) - 2 * mu2^(y2 - 1)))/factorial(y2)  
der2pdf2.dermu2       <- as.numeric( der2pdf2.dermu2FUNCpo(y2, mu2) )
    
derpdf2.sigma2        <- 0
der2pdf2.dersigma22   <- 0
der2pdf2.mu2dersigma2 <- 0



if(naive == FALSE){   # needs y2m 
 
p2  <- pPO(as.numeric(y2), mu = as.numeric(mu2)) 
 
mu2 <- c(mu2)

derp2.dermu2           <- rowSums( matrix(as.numeric(derpdf2.dermu2FUNCpo(y2m, mu2)), dim(y2m)[1],dim(y2m)[2]), na.rm = TRUE ) 
derp2.dersigma2        <- 0
der2p2.dermu22         <- rowSums( matrix(as.numeric(der2pdf2.dermu2FUNCpo(y2m, mu2)),dim(y2m)[1],dim(y2m)[2]), na.rm = TRUE ) 
der2p2.dersigma22      <- 0
der2p2.derdermu2sigma2 <- 0


                      
    
                   }

}



if(margin2 == "ZTP"){

mu2 <- exp(eta2)
dermu2.dereta2 <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2) 

dersigma2.dersigma2.st  <- 0
dersigma2.dersigma2.st2 <- 0

if(max(y2) > 170){

prec <- pmax(53, getPrec(mu2), getPrec(y2))
        
mu2 <- mpfr(mu2, prec)
y2  <- mpfr( y2, prec) 

}


pdf2FUNCztp <- function(y2, mu2) mu2^y2/(exp(mu2)-1)*1/factorial(y2)  
pdf2        <- as.numeric( pdf2FUNCztp(y2, mu2) )  

derpdf2.dermu2FUNCztp <- function(y2, mu2) (mu2^(y2 - 1) * y2 - mu2^y2 * exp(mu2)/(exp(mu2) - 1))/(factorial(y2) * (exp(mu2) - 1))  
derpdf2.dermu2        <- as.numeric( derpdf2.dermu2FUNCztp(y2, mu2) )
    
der2pdf2.dermu2FUNCztp <- function(y2, mu2) (y2 * (mu2^(y2 - 2) * (y2 - 1) - mu2^(y2 - 1) * exp(mu2)/(exp(mu2) - 
    1)) - exp(mu2) * (mu2^(y2 - 1) * y2 + mu2^y2 - 2 * (mu2^y2 * 
    exp(mu2)/(exp(mu2) - 1)))/(exp(mu2) - 1))/(factorial(y2) * (exp(mu2) - 
    1))   
der2pdf2.dermu2        <- as.numeric( der2pdf2.dermu2FUNCztp(y2, mu2) ) 
    
derpdf2.sigma2        <- 0
der2pdf2.dersigma22   <- 0
der2pdf2.mu2dersigma2 <- 0


if(naive == FALSE){   #needs y2m

mu2 <- c(mu2)

p2  <- rowSums( matrix(as.numeric(pdf2FUNCztp(y2m, mu2)),dim(y2m)[1],dim(y2m)[2]), na.rm = TRUE )

derp2.dermu2           <- rowSums( matrix(as.numeric(derpdf2.dermu2FUNCztp(y2m, mu2)), dim(y2m)[1],dim(y2m)[2]), na.rm = TRUE )
derp2.dersigma2        <- 0
der2p2.dermu22         <- rowSums( matrix(as.numeric(der2pdf2.dermu2FUNCztp(y2m, mu2)),dim(y2m)[1],dim(y2m)[2]), na.rm = TRUE )
der2p2.dersigma22      <- 0
der2p2.derdermu2sigma2 <- 0
                      
    
                   }

}




####


# 11 oct 2018, decided to use all numerical for now...

if(margin2 == "NBI"){ # all Numerical - does not need y2m

sigma <- sigma2 <- ifelse(sigma2 < 4.151334e-06, 4.151334e-06, sigma2) # related to gamma function
sigma2.st <- log(sigma2) 

mu2 <- exp(eta2)
dermu2.dereta2 <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2) 

dersigma2.dersigma2.st  <- exp(sigma2.st)  
dersigma2.dersigma2.st2 <- exp(sigma2.st)  


pdf2 <- dNBI(y2, mu = mu2, sigma = sigma)  

derpdf2.dermu2F <- derpdf2.dermu2FUNC2p(function(mu2) dNBI(y2, mu = mu2, sigma = sigma), y2, mu2, sigma) 
derpdf2.dermu2  <- derpdf2.dermu2F$fi 
der2pdf2.dermu2 <- derpdf2.dermu2F$se  


derpdf2.sigma2F     <- derpdf2.sigma2FUNC2p(function(sigma) dNBI(y2, mu = mu2, sigma = sigma), y2, mu2, sigma) 
derpdf2.sigma2      <- derpdf2.sigma2F$fi      
der2pdf2.dersigma22 <- derpdf2.sigma2F$se    


der2pdf2.mu2dersigma2 <- der2pdf2.mu2dersigma2FUNC2p(function(mu2, sigma) dNBI(y2, mu = mu2, sigma = sigma), y2, mu2, sigma)


if(naive == FALSE){   
 
p2  <- pNBI(y2, mu = mu2, sigma = sigma)  
 
derp2.dermu2F  <- derpdf2.dermu2FUNC2p(function(mu2) pNBI(y2, mu = mu2, sigma = sigma), y2, mu2, sigma)
derp2.dermu2   <- derp2.dermu2F$fi
der2p2.dermu22 <- derp2.dermu2F$se

derp2.dersigma2F <- derpdf2.sigma2FUNC2p(function(sigma) pNBI(y2, mu = mu2, sigma = sigma), y2, mu2, sigma) 
derp2.dersigma2  <- derp2.dersigma2F$fi 
der2p2.dersigma22<- derp2.dersigma2F$se 

der2p2.derdermu2sigma2 <- der2pdf2.mu2dersigma2FUNC2p(function(mu2, sigma) pNBI(y2, mu = mu2, sigma = sigma), y2, mu2, sigma) 

  
                   }

}





#######





if(margin2 == "NBII"){ # all Numerical - does not need y2m

sigma <- sigma2 <- ifelse(sigma2 < 4.151334e-06, 4.151334e-06, sigma2) # related to gamma function
sigma2.st <- log(sigma2) 

mu2 <- exp(eta2)
dermu2.dereta2 <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2) 

dersigma2.dersigma2.st  <- exp(sigma2.st)  
dersigma2.dersigma2.st2 <- exp(sigma2.st)  


pdf2 <- dNBII(y2, mu = mu2, sigma = sigma)  

derpdf2.dermu2F <- derpdf2.dermu2FUNC2p(function(mu2) dNBII(y2, mu = mu2, sigma = sigma), y2, mu2, sigma) 
derpdf2.dermu2  <- derpdf2.dermu2F$fi 
der2pdf2.dermu2 <- derpdf2.dermu2F$se  


derpdf2.sigma2F     <- derpdf2.sigma2FUNC2p(function(sigma) dNBII(y2, mu = mu2, sigma = sigma), y2, mu2, sigma) 
derpdf2.sigma2      <- derpdf2.sigma2F$fi      
der2pdf2.dersigma22 <- derpdf2.sigma2F$se    


der2pdf2.mu2dersigma2 <- der2pdf2.mu2dersigma2FUNC2p(function(mu2, sigma) dNBII(y2, mu = mu2, sigma = sigma), y2, mu2, sigma)


if(naive == FALSE){   
 
p2  <- pNBII(y2, mu = mu2, sigma = sigma)  
 
derp2.dermu2F  <- derpdf2.dermu2FUNC2p(function(mu2) pNBII(y2, mu = mu2, sigma = sigma), y2, mu2, sigma)
derp2.dermu2   <- derp2.dermu2F$fi
der2p2.dermu22 <- derp2.dermu2F$se

derp2.dersigma2F <- derpdf2.sigma2FUNC2p(function(sigma) pNBII(y2, mu = mu2, sigma = sigma), y2, mu2, sigma) 
derp2.dersigma2  <- derp2.dersigma2F$fi 
der2p2.dersigma22<- derp2.dersigma2F$se 

der2p2.derdermu2sigma2 <- der2pdf2.mu2dersigma2FUNC2p(function(mu2, sigma) pNBII(y2, mu = mu2, sigma = sigma), y2, mu2, sigma) 

  
                   }

}





#############################

if(margin2 == "PIG"){ # K all numerical as well
                      # sigma <- ifelse(sigma>precision, sigma, precision)
                      # sigma <- ifelse(sigma<10^7, sigma, 10^7)
                      # tolerances for derivatives may be changed to get better
                      # performace, difficult to find a general rule here
    
    

mu2 <- exp(eta2)
dermu2.dereta2 <- exp(eta2)
der2mu2.dereta2eta2 <- exp(eta2) 

dersigma2.dersigma2.st  <- exp(sigma2.st)  
dersigma2.dersigma2.st2 <- exp(sigma2.st)  


pdf2 <- dPIG(y2, mu = mu2, sigma = sigma)    

derpdf2.dermu2F <- derpdf2.dermu2FUNC2p(function(mu2) dPIG(y2, mu = mu2, sigma = sigma), y2, mu2, sigma) 
derpdf2.dermu2  <- derpdf2.dermu2F$fi 
der2pdf2.dermu2 <- derpdf2.dermu2F$se     
    
 
derpdf2.sigma2F     <- derpdf2.sigma2FUNC2p(function(sigma) dPIG(y2, mu = mu2, sigma = sigma), y2, mu2, sigma) 
derpdf2.sigma2      <- derpdf2.sigma2F$fi      
der2pdf2.dersigma22 <- derpdf2.sigma2F$se 
   
der2pdf2.mu2dersigma2 <- der2pdf2.mu2dersigma2FUNC2p(function(mu2, sigma) dPIG(y2, mu = mu2, sigma = sigma), y2, mu2, sigma)


if(naive == FALSE){   
 
p2  <- pPIG(y2, mu = mu2, sigma = sigma) 
 
derp2.dermu2F  <- derpdf2.dermu2FUNC2p(function(mu2) pPIG(y2, mu = mu2, sigma = sigma), y2, mu2, sigma)
derp2.dermu2   <- derp2.dermu2F$fi
der2p2.dermu22 <- derp2.dermu2F$se
    
derp2.dersigma2F <- derpdf2.sigma2FUNC2p(function(sigma) pPIG(y2, mu = mu2, sigma = sigma), y2, mu2, sigma) 
derp2.dersigma2  <- derp2.dersigma2F$fi 
der2p2.dersigma22<- derp2.dersigma2F$se 
   
der2p2.derdermu2sigma2 <- der2pdf2.mu2dersigma2FUNC2p(function(mu2, sigma) pPIG(y2, mu = mu2, sigma = sigma), y2, mu2, sigma) 
  

                   }

}



if(margin2 %in% c("DGP","DGPII") ){




if(margin2 == "DGP"){
    
                   mu2  <- eta2 # exi
dermu2.dereta2          <- 1
der2mu2.dereta2eta2     <- 0 


}


if(margin2 == "DGPII"){

                   mu2  <- eta2^2 # exi
dermu2.dereta2          <- 2*eta2
der2mu2.dereta2eta2     <- 2 

}


dersigma2.dersigma2.st  <- exp(sigma2.st)  # sigma
dersigma2.dersigma2.st2 <- exp(sigma2.st)  


indx1 <- as.numeric( ((1 + mu2*y2/sigma)     > 0) == FALSE ) 
indx2 <- as.numeric( ((1 + mu2*(y2+1)/sigma) > 0) == FALSE ) # this is not needed
indx  <- rowSums(cbind(indx1, indx2))


pdf2FUNC <- function(y2, mu2, sigma) suppressWarnings(   (1 + mu2*y2/sigma)^(-1/mu2) - (1 + mu2*(1+y2)/sigma)^(-1/mu2)    )  


derpdf2.dermu2FUNC        <- function(y2, mu2, sigma) suppressWarnings(  (((1 + y2)/(1 + mu2 * (1 + y2)/sigma)^(1 + 1/mu2) - y2/(1 + mu2 * 
    y2/sigma)^(1 + 1/mu2))/sigma + (log1p(mu2 * y2/sigma)/(1 + 
    mu2 * y2/sigma)^(1/mu2) - log1p(mu2 * (1 + y2)/sigma)/(1 + 
    mu2 * (1 + y2)/sigma)^(1/mu2))/mu2)/mu2   )


derpdf2.sigma2FUNC        <- function(y2, mu2, sigma)  suppressWarnings(  -(((1 + y2)/(1 + mu2 * (1 + y2)/sigma)^(1 + 1/mu2) - y2/(1 + 
    mu2 * y2/sigma)^(1 + 1/mu2))/sigma^2)  )

       
der2pdf2.dermu2FUNC       <- function(y2, mu2, sigma) suppressWarnings( (((log1p(mu2 * y2/sigma) * (log1p(mu2 * y2/sigma)/(mu2 * (1 + 
    mu2 * y2/sigma)^(1/mu2)) - y2/(sigma * (1 + mu2 * y2/sigma)^(1 + 
    1/mu2))) - log1p(mu2 * (1 + y2)/sigma) * (log1p(mu2 * (1 + 
    y2)/sigma)/(mu2 * (1 + mu2 * (1 + y2)/sigma)^(1/mu2)) - (1 + 
    y2)/(sigma * (1 + mu2 * (1 + y2)/sigma)^(1 + 1/mu2))))/mu2 + 
    (y2/(sigma * (1 + mu2 * y2/sigma)) - 2 * (log1p(mu2 * y2/sigma)/mu2))/(1 + 
        mu2 * y2/sigma)^(1/mu2) - ((1 + y2)/(sigma * (1 + mu2 * 
    (1 + y2)/sigma)) - 2 * (log1p(mu2 * (1 + y2)/sigma)/mu2))/(1 + 
    mu2 * (1 + y2)/sigma)^(1/mu2))/mu2 + (y2 * ((1/(1 + mu2 * 
    y2/sigma)^(1 + 1/mu2) - log1p(mu2 * y2/sigma)/(mu2 * (1 + 
    mu2 * y2/sigma)^(1 + 1/mu2)))/mu2 + y2 * (1 + 1/mu2)/(sigma * 
    (1 + mu2 * y2/sigma)^(1/mu2 + 2))) - ((1 + 1/mu2) * (1 + 
    y2)/(sigma * (1 + mu2 * (1 + y2)/sigma)^(1/mu2 + 2)) + (1/(1 + 
    mu2 * (1 + y2)/sigma)^(1 + 1/mu2) - log1p(mu2 * (1 + y2)/sigma)/(mu2 * 
    (1 + mu2 * (1 + y2)/sigma)^(1 + 1/mu2)))/mu2) * (1 + y2))/sigma)/mu2  )

    
der2pdf2.dersigma22FUNC   <- function(y2, mu2, sigma) suppressWarnings(  (y2 * (mu2 * y2 * (1 + 1/mu2)/(sigma * (1 + mu2 * y2/sigma)^(1/mu2 + 
    2)) - 2/(1 + mu2 * y2/sigma)^(1 + 1/mu2)) - (1 + y2) * (mu2 * 
    (1 + 1/mu2) * (1 + y2)/(sigma * (1 + mu2 * (1 + y2)/sigma)^(1/mu2 + 
    2)) - 2/(1 + mu2 * (1 + y2)/sigma)^(1 + 1/mu2)))/sigma^3  )

    
der2pdf2.mu2dersigma2FUNC <- function(y2, mu2, sigma) suppressWarnings(  -(((1 + y2) * (log1p(mu2 * (1 + y2)/sigma)/(mu2^2 * (1 + mu2 * 
    (1 + y2)/sigma)^(1 + 1/mu2)) - (1 + 1/mu2) * (1 + y2)/(sigma * 
    (1 + mu2 * (1 + y2)/sigma)^(1/mu2 + 2))) - y2 * (log1p(mu2 * 
    y2/sigma)/(mu2^2 * (1 + mu2 * y2/sigma)^(1 + 1/mu2)) - y2 * 
    (1 + 1/mu2)/(sigma * (1 + mu2 * y2/sigma)^(1/mu2 + 2))))/sigma^2)    )




pdf2                  <- as.numeric( pdf2FUNC(y2, mu2, sigma) )  
derpdf2.dermu2        <- derpdf2.dermu2FUNC(y2, mu2, sigma) 
derpdf2.sigma2        <- derpdf2.sigma2FUNC(y2, mu2, sigma) 
der2pdf2.dermu2       <- der2pdf2.dermu2FUNC(y2, mu2, sigma) 
der2pdf2.dersigma22   <- der2pdf2.dersigma22FUNC(y2, mu2, sigma) 
der2pdf2.mu2dersigma2 <- der2pdf2.mu2dersigma2FUNC(y2, mu2, sigma)



pdf2                  <- ifelse( indx == 0, pdf2, 1) # 'cause in output log(1) = 0 hence it will not add anything to the lik
derpdf2.dermu2        <- ifelse( indx == 0, derpdf2.dermu2, 0) 
derpdf2.sigma2        <- ifelse( indx == 0, derpdf2.sigma2, 0) 
der2pdf2.dermu2       <- ifelse( indx == 0, der2pdf2.dermu2, 0) 
der2pdf2.dersigma22   <- ifelse( indx == 0, der2pdf2.dersigma22, 0) 
der2pdf2.mu2dersigma2 <- ifelse( indx == 0, der2pdf2.mu2dersigma2, 0) 








if(naive == FALSE){   # needs y2m
 
ly2 <- length(y2)
if(length(sigma) == 1) sigma <- c(rep(sigma, ly2))

mu2    <- c(mu2)
sigma <- c(sigma)

derp2.dermu2           <- rowSums( derpdf2.dermu2FUNC(        y2m, mu2, sigma ) , na.rm = TRUE )
derp2.dersigma2        <- rowSums( derpdf2.sigma2FUNC(        y2m, mu2, sigma ) , na.rm = TRUE )
der2p2.dermu22         <- rowSums( der2pdf2.dermu2FUNC(       y2m, mu2, sigma ) , na.rm = TRUE )
der2p2.dersigma22      <- rowSums( der2pdf2.dersigma22FUNC(   y2m, mu2, sigma ) , na.rm = TRUE )
der2p2.derdermu2sigma2 <- rowSums( der2pdf2.mu2dersigma2FUNC( y2m, mu2, sigma ) , na.rm = TRUE )


                   }

}






##########################################################




if(margin2 %in% c(cont2par,cont3par)){ 

derpdf2.dereta2              <- derpdf2.dermu2*dermu2.dereta2       
der2pdf2.dereta2             <- der2pdf2.dermu2* dermu2.dereta2^2 + derpdf2.dermu2*der2mu2.dereta2eta2        


der2pdf2.dereta2dersigma2    <- der2pdf2.mu2dersigma2* dermu2.dereta2
der2pdf2.dereta2dersigma2.st <- der2pdf2.dereta2dersigma2 *  dersigma2.dersigma2.st

  
derpdf2.dersigma2.st         <- derpdf2.sigma2 * dersigma2.dersigma2.st   
der2pdf2.dersigma2.st2       <- der2pdf2.dersigma22 * dersigma2.dersigma2.st^2 + derpdf2.sigma2  * dersigma2.dersigma2.st2     

                  

if(naive == FALSE){  



derp2.dereta2                <- derp2.dermu2*dermu2.dereta2
der2p2.dereta2eta2           <- der2p2.dermu22*dermu2.dereta2^2+derp2.dermu2*der2mu2.dereta2eta2      


der2p2.dereta2dersigma2      <- der2p2.derdermu2sigma2* dermu2.dereta2    
der2p2.dereta2dersigma2.st   <- der2p2.dereta2dersigma2 *  dersigma2.dersigma2.st  

derp2.dersigma.st            <- derp2.dersigma2 *  dersigma2.dersigma2.st 
der2p2.dersigma2.st2         <- der2p2.dersigma22 * dersigma2.dersigma2.st^2 + derp2.dersigma2 * dersigma2.dersigma2.st2





                 }
                 
                 
                 

}###############



epsilon <- sqrt(.Machine$double.eps)
max.p   <- 0.9999999

  pdf2 <- ifelse(pdf2 < epsilon, epsilon, pdf2 )

  p2   <- mm(p2) 



ifef <- function(dv){

epsilon <- sqrt(.Machine$double.eps)
dv <- ifelse(is.na(dv), epsilon, dv ) 
dv <- ifelse(dv == Inf ,  8.218407e+20, dv )
dv <- ifelse(dv == -Inf ,  -8.218407e+20, dv )
dv

}


# for safety

pdf2                         = ifef(pdf2)
p2                           = ifef(p2) 
derpdf2.dereta2              = ifef(derpdf2.dereta2) 
derpdf2.dersigma2.st         = ifef(derpdf2.dersigma2.st) 
derp2.dersigma.st            = ifef(derp2.dersigma.st)
derp2.dereta2                = ifef(derp2.dereta2)
der2p2.dereta2eta2           = ifef(der2p2.dereta2eta2) 
der2pdf2.dereta2             = ifef(der2pdf2.dereta2)
der2p2.dersigma2.st2         = ifef(der2p2.dersigma2.st2) 
der2pdf2.dersigma2.st2       = ifef(der2pdf2.dersigma2.st2)
der2p2.dereta2dersigma2.st   = ifef(der2p2.dereta2dersigma2.st)  
der2pdf2.dereta2dersigma2.st = ifef(der2pdf2.dereta2dersigma2.st)
der2pdf2.dereta2dernu.st     = ifef(der2pdf2.dereta2dernu.st)   
der2pdf2.sigma2.st2dernu.st  = ifef(der2pdf2.sigma2.st2dernu.st)
derpdf2.dernu.st             = ifef(derpdf2.dernu.st)           
der2pdf2.dernu.st2           = ifef(der2pdf2.dernu.st2)         
derp2.nu.st                  = ifef(derp2.nu.st)                
der2p2.dernu.st2             = ifef(der2p2.dernu.st2)           
der2p2.dereta2dernu.st       = ifef(der2p2.dereta2dernu.st)     
der2p2.dersigma2.stdernu.st  = ifef(der2p2.dersigma2.stdernu.st)




list(pdf2                         = pdf2,
     p2                           = p2, 
     derpdf2.dereta2              = derpdf2.dereta2, 
     derpdf2.dersigma2.st         = derpdf2.dersigma2.st, 
     derp2.dersigma.st            = derp2.dersigma.st,
     derp2.dereta2                = derp2.dereta2,
     der2p2.dereta2eta2           = der2p2.dereta2eta2, 
     der2pdf2.dereta2             = der2pdf2.dereta2,
     der2p2.dersigma2.st2         = der2p2.dersigma2.st2, 
     der2pdf2.dersigma2.st2       = der2pdf2.dersigma2.st2,
     der2p2.dereta2dersigma2.st   = der2p2.dereta2dersigma2.st,            
     der2pdf2.dereta2dersigma2.st = der2pdf2.dereta2dersigma2.st,
     der2pdf2.dereta2dernu.st     = der2pdf2.dereta2dernu.st,   
     der2pdf2.sigma2.st2dernu.st  = der2pdf2.sigma2.st2dernu.st,
     derpdf2.dernu.st             = derpdf2.dernu.st,           
     der2pdf2.dernu.st2           = der2pdf2.dernu.st2,         
     derp2.nu.st                  = derp2.nu.st,                
     der2p2.dernu.st2             = der2p2.dernu.st2,           
     der2p2.dereta2dernu.st       = der2p2.dereta2dernu.st,     
     der2p2.dersigma2.stdernu.st  = der2p2.dersigma2.stdernu.st, 
     indx = indx == 0)     


}




    