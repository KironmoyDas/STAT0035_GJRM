gjrm <- function(formula, data = list(), weights = NULL, subset = NULL,  
                             BivD = "N", margins, Model, dof = 3, ordinal = FALSE,
                             surv = FALSE, cens1 = NULL, cens2 = NULL, cens3 = NULL, dep.cens = FALSE,  
                             gamlssfit = FALSE, fp = FALSE, infl.fac = 1, 
                             rinit = 1, rmax = 100, iterlimsp = 50, tolsp = 1e-07,
                             gc.l = FALSE, parscale, extra.regI = "t", k1.tvc = 0, k2.tvc = 0, 
                             knots = NULL, penCor = "unpen",
                             sp.penCor = 3, Chol = FALSE, gamma = 1, w.alasso = NULL,
                             drop.unused.levels = TRUE, ind.ord = FALSE,
                             min.dn = 1e-40, min.pr = 1e-16, max.pr = 0.999999){
  
  
  bl <- c("probit", "logit", "cloglog")  
  

  #print(( !(margins[1] %in% bl) || surv == TRUE ) && ordinal == FALSE)
  if(  ( !(margins[1] %in% bl) || surv == TRUE ) && ordinal == FALSE ){
    
  ##########################################################################################################################
  # preamble
  ##########################################################################################################################  
  
  robust <- FALSE; t.c = 3
  sp <- qu.mag <- y1.y2 <- y1.cy2 <- cy1.y2 <- cy1.cy2 <- cy <- cy1 <- gamlss1 <- gamlss2 <- gam1 <- gam2 <- y1m <- y2m <- indexTeq1 <- indexTeq2 <- NULL  
  i.rho <- log.sig2.2 <- log.nu.2 <- log.nu.1 <- log.sig2.1 <- dof.st <- NULL
  end <- X3.d2 <- X4.d2 <- X5.d2 <- X6.d2 <- X7.d2 <- X8.d2 <- l.sp3 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- l.sp8 <- 0
  sp1 <- sp2 <- NULL
  sp3 <- gp3 <- gam3 <- X3 <- sp4 <- gp4 <- gam4 <- X4 <- sp5 <- gp5 <- gam5 <- X5 <- NULL    
  sp6 <- gp6 <- gam6 <- X6 <- sp7 <- gp7 <- gam7 <- X7 <- sp8 <- gp8 <- gam8 <- X8 <- NULL   
  c11 <- c10 <- c01 <- c00 <- NA
  
  Sl.sf <- NULL
  sp.method <- "perf"
  
  Xd1 <- Xd2 <- mono.sm.pos1 <- mono.sm.pos2 <- mono.sm.pos <- NULL
  surv.flex <- FALSE
  
  Deq1 <- pos.pbeq1 <- Deq2 <- pos.pbeq2 <- list()
  
  ###################################
  
  BivD2 <- c("C0C90","C0C270","C180C90","C180C270",
             "J0J90","J0J270","J180J90","J180J270",
             "G0G90","G0G270","G180G90","G180G270")
             
  opc  <- c("N","C0","C90","C180","C270","J0","J90","J180","J270","G0","G90","G180","G270","F","AMH","FGM","T","PL","HO")
  scc  <- c("C0", "C180", "J0", "J180", "G0", "G180", BivD2)
  sccn <- c("C90", "C270", "J90", "J270", "G90", "G270")
  m2   <- c("N","GU","rGU","LO","LN","WEI","iG","GA","BE","FISK","GP","GPII","GPo")
  m3   <- c("DAGUM","SM","TW")
  m1d  <- c("PO", "ZTP")
  m2d  <- c("NBI", "NBII","PIG","DGP","DGPII")
  m3d  <- c("DEL","SICHEL")
  
  ct  <- data.frame( c(opc), c(1:14,55,56,57,60,61) )
  cta <- data.frame( c(opc), c(1,3,23,13,33,6,26,16,36,4,24,14,34,5,55,56,2,60,61) )     
  
  #print(!(BivD %in% BivD2))
  if(!(BivD %in% BivD2)){
    
  nC  <-  ct[which( ct[,1]==BivD),2]
  nCa <- cta[which(cta[,1]==BivD),2]     
    
  }
  
 #######################################################################################  
 
  l.flist <- length(formula)

  form.check(formula, l.flist) 
  cl <- match.call()       
  mf <- match.call(expand.dots = FALSE)
            

  pred.varEXPANDED <- function(formula, l.flist, gaml = FALSE, triv = FALSE, informative = "no"){
    
    ig <- interpret.gam(formula)
    v3 <- v2 <- NULL
    
    #print(gaml == FALSE && triv == FALSE)
    if(gaml == FALSE && triv == FALSE){    
      
      or1 <- as.character(formula[[1]][2])
      or2 <- as.character(formula[[2]][2])
      
      #print(l.flist == 5)
      if( l.flist == 5 ){  
        v1 <- all.vars(as.formula(formula[[1]]))[1]
        v1 <- c(v1, ig[[1]]$pred.names)
        v2 <- all.vars(as.formula(formula[[2]]))[1]
        v2 <- c(v2, ig[[2]]$pred.names)
        v3 <- ig[[3]]$pred.names 
        v4 <- ig[[4]]$pred.names
        v5 <- ig[[5]]$pred.names 
        pred.n <- union(v1,c(v2,v3,v4,v5,or1,or2))
      }   
      
    }

    list(v1 = v1, v2 = v2, v3 = v3, pred.n = pred.n)      
    
  }
  
  pred.varR <- pred.varEXPANDED(formula, l.flist)


  v1     <- pred.varR$v1  
  v2     <- pred.varR$v2
  pred.n <- pred.varR$pred.n  

  fake.formula <- paste(v1[1], "~", paste(pred.n, collapse = " + ")) 
  environment(fake.formula) <- environment(formula[[1]])
  mf$formula <- fake.formula 
  
  mf$min.dn <- mf$min.pr <- mf$max.pr <- mf$dep.cens <- mf$ordinal <- mf$Model <- mf$knots <- mf$k1.tvc <- mf$k2.tvc <- mf$surv <- mf$BivD <- mf$margins <- mf$fp <- mf$dof <- mf$infl.fac <- mf$rinit <- mf$rmax <- mf$iterlimsp <- mf$tolsp <- mf$gc.l <- mf$parscale <- mf$extra.regI <- mf$gamlssfit <- NULL                           
  mf$drop.unused.levels <- drop.unused.levels 
  mf[[1]] <- as.name("model.frame")
  data <- eval(mf, parent.frame())
 
  n <- dim(data)[1]
        
  #print(!("(weights)" %in% names(data)))
  if(!("(weights)" %in% names(data))) {weights <- rep(1,dim(data)[1]) 
                        data$weights <- weights
                        names(data)[length(names(data))] <- "(weights)"}
  
  #print(!("(cens1)" %in% names(data)))
  if(!("(cens1)" %in% names(data))) {cens1 <- rep(0,dim(data)[1]) 
                        data$cens1 <- cens1
                        names(data)[length(names(data))] <- "(cens1)"}
  
  #print(!("(cens2)" %in% names(data)))
  if(!("(cens2)" %in% names(data))) {cens2 <- rep(0,dim(data)[1]) 
                        data$cens2 <- cens2
                        names(data)[length(names(data))] <- "(cens2)"}
  
  #print(!("(cens3)" %in% names(data)))
  if(!("(cens3)" %in% names(data))) {cens3 <- rep(0,dim(data)[1]) 
                        data$cens3 <- cens3
                        names(data)[length(names(data))] <- "(cens3)"}                          
                        
  M <- list(m1d = m1d, m2 = m2, m2d = m2d, m3 = m3, m3d = m3d, BivD = BivD, bl = bl, 
            robust = robust, opc = opc, extra.regI = extra.regI, margins = margins, BivD2 = BivD2, dof = dof,
            surv = surv, c1 = cens1, c2 = cens2, c3 = cens3, dep.cens = dep.cens) 
 
  M$K1 <- NULL
  
  formula.eq1 <- formula[[1]]
  formula.eq2 <- formula[[2]] 
    
 ##############################################################  
 # Equation 1
 ############################################################## 

  form.eq12EXPANDED <- function(formula.eq1, data, v1, margins, m1d, m2d, copSS = FALSE, inde = NULL){
    
    y1m <- f.eq1 <- NULL
    
    formula.eq1r <- formula.eq1   
    y1 <- y1.test <- data[, v1[1]]
    
    #print( margins %in% c("N","LO","GU","rGU","GAi","TW") )
    if( margins %in% c("N","LO","GU","rGU","GAi","TW") )   formula.eq1 <- update(formula.eq1, (. + mean(.))/2 ~ . ) 
    
    f.eq1LI <- temp.respV ~ urcfcphmwicu # specific to surv model with L and I
    
    list(f.eq1LI = f.eq1LI, formula.eq1 = formula.eq1, formula.eq1r = formula.eq1r, y1 = y1, y1.test = y1.test, y1m = y1m, f.eq1 = f.eq1)
    
  }
  

 form.eq12R <- form.eq12EXPANDED(formula.eq1, data, v1, margins[1], m1d, m2d)   
 
 formula.eq1  <- form.eq12R$formula.eq1
 formula.eq1r <- form.eq12R$formula.eq1r
 y1           <- form.eq12R$y1
 y1.test      <- form.eq12R$y1.test 
 y1m          <- form.eq12R$y1m

 #print(surv == FALSE)
 if(surv == FALSE)  gam1 <- eval(substitute(gam(formula.eq1, gamma=infl.fac, weights=weights, data=data, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights)))

 gam1$formula <- formula.eq1r  
 lsgam1 <- length(gam1$smooth)
 
 y1 <- y1.test
 
 attr(data,"terms") <- NULL ## to make it work when using log(y1) for instance, this will have to be checked if we need it or not ##
 
 #print(!(surv == TRUE && margins[1] %in% bl))
 if(!(surv == TRUE && margins[1] %in% bl)){
 
     names(gam1$model)[1] <- as.character(formula.eq1r[2])
     X1 <- predict(gam1, type = "lpmatrix")
     l.sp1 <- length(gam1$sp)
     sp1 <- gam1$sp
 }
 
 gp1 <- gam1$nsdf 
 X1.d2 <- dim(X1)[2]
 
 ##############################################################
 # Equation 2 
 ############################################################## 
 
 form.eq12EXPANDED2 <- function(formula.eq1, data, v1, margins, m1d, m2d, copSS = FALSE, inde = NULL){
   
   y1m <- f.eq1 <- NULL
   
   formula.eq1r <- formula.eq1   
   y1 <- y1.test <- data[, v1[1]]
   
   
   if( margins %in% c("N","LO","GU","rGU","GAi","TW") )   formula.eq1 <- update(formula.eq1, (. + mean(.))/2 ~ . ) 
   
   f.eq1LI <- temp.respV ~ urcfcphmwicu # specific to surv model with L and I
   
   list(f.eq1LI = f.eq1LI, formula.eq1 = formula.eq1, formula.eq1r = formula.eq1r, y1 = y1, y1.test = y1.test, y1m = y1m, f.eq1 = f.eq1)
   
 }
 
 

 form.eq12R <- form.eq12EXPANDED2(formula.eq2, data, v2, margins[2], m1d, m2d)   
 
 formula.eq2  <- form.eq12R$formula.eq1
 formula.eq2r <- form.eq12R$formula.eq1r
 y2           <- form.eq12R$y1
 y2.test      <- form.eq12R$y1.test 
 y2m          <- form.eq12R$y1m

 #print(surv == FALSE)
 if(surv == FALSE)  gam2 <- eval(substitute(gam(formula.eq2, gamma=infl.fac, weights=weights, data=data, knots = knots, drop.unused.levels = drop.unused.levels),list(weights=weights)))
 
 gam2$formula <- formula.eq2r  
 lsgam2 <- length(gam2$smooth)
 
 y2 <- y2.test
 
 attr(data,"terms") <- NULL ## to make it work when using log(y1) for instance, this will have to be checked if we need it or not ##
 
 #print(!(surv == TRUE && margins[2] %in% bl))
 if( !(surv == TRUE && margins[2] %in% bl) ){
 
     names(gam2$model)[1] <- as.character(formula.eq2r[2])
     X2 <- predict(gam2, type = "lpmatrix")
     l.sp2 <- length(gam2$sp)
     sp2 <- gam2$sp
 }
 
 gp2 <- gam2$nsdf 
 X2.d2 <- dim(X2)[2]
  
#################################################################
# Starting value for dependence parameter (and dof for T if used)
#################################################################

 res1 <- residuals(gam1)
 res2 <- residuals(gam2)
 
 ass.s <- cor(res1, res2, method = "kendall")
 ass.s <- sign(ass.s)*ifelse(abs(ass.s) > 0.9, 0.9, abs(ass.s))
 
 
 ass.dpEXPANDED <- function(ass.s, BivD, scc, sccn, nCa){
   
   eps <- sqrt(.Machine$double.eps) # this looks fine here and it is not that dangerous 
   
   #print(!(BivD %in% c("AMH","FGM","PL","HO")))
   if(!(BivD %in% c("AMH","FGM","PL","HO"))) i.rho <- BiCopTau2Par(family = nCa, tau = ass.s)
   
   #print(BivD %in% c("N","AMH","FGM","T"))
   if(BivD %in% c("N","AMH","FGM","T"))         i.rho <- atanh( i.rho )
   
   names(i.rho) <- "theta.star"   
   
   i.rho
   
 }
 
 
 i.rho <- ass.dpEXPANDED(ass.s, BivD, scc, sccn, nCa)

 dof.st <- log(dof - 2) 
 names(dof.st) <- "dof.star"   
                           
##############################################################
# Other starting values + overall
##############################################################
 
#print(!(margins[1] %in% c(m1d,bl)) )
 if( !(margins[1] %in% c(m1d,bl)) ){
   
   
   startsnEXPANDED <- function(margins, y1){
     
     log.nu.1 <- NULL      
     
     #print(!(margins %in% c("GO")))
     if( !(margins %in% c("GO")) ) par.est <- try( resp.check(y1, margin = margins, plots = FALSE, print.par = TRUE, i.f = TRUE), silent = TRUE)
     
     #print(!(margins %in% c("GO")))
     if( !(margins %in% c("GO")) ){
       
       log.sig2.1 <- par.est[2]
       
     }
     
     list(log.sig2.1 = log.sig2.1, log.nu.1 = log.nu.1)        
     
     
   }
   


 start.snR <- startsnEXPANDED(margins[1], y1)
    
 log.sig2.1 <- start.snR$log.sig2.1; names(log.sig2.1) <- "sigma1.star"

 }


#print(!(margins[2] %in% c(m1d,bl)))
 if( !(margins[2] %in% c(m1d,bl)) ){
   
   
   startsnEXPANDED2 <- function(margins, y1){
     
     log.nu.1 <- NULL      

     #print( !(margins %in% c("GO")))
     if( !(margins %in% c("GO")) ) par.est <- try( resp.check(y1, margin = margins, plots = FALSE, print.par = TRUE, i.f = TRUE), silent = TRUE)
     
     
     
     #print(!(margins %in% c("GO")))
     if( !(margins %in% c("GO")) ){
       
       
       log.sig2.1 <- par.est[2]
       
       if( margins %in% c("DAGUM","SM") ){}                                                                                    
       
     }
     
     list(log.sig2.1 = log.sig2.1, log.nu.1 = log.nu.1)        
     
   }
   

   

 start.snR <- startsnEXPANDED2(margins[2], y2)
    
 log.sig2.2 <- start.snR$log.sig2.1; names(log.sig2.2) <- "sigma2.star"


 }

 vo <- list(gam1 = gam1, gam2 = gam2, i.rho = i.rho, log.sig2.2 = log.sig2.2, log.nu.2 = log.nu.2, log.nu.1 = log.nu.1, log.sig2.1 = log.sig2.1, dof.st = dof.st, n = n, drop.unused.levels = drop.unused.levels )

 
 
 overall.svEXPANDED <- function(margins, M, vo, type = "copR", c.gam2 = NULL){
   

   #print(type == "copR")
   if(type == "copR"){
     
     BivD <- M$BivD
     
       #print(margins[1] %in% c(M$m2,M$m2d) && margins[2] %in% c(M$m2,M$m2d))
       if( margins[1] %in% c(M$m2,M$m2d) && margins[2] %in% c(M$m2,M$m2d) ) start.v <- c(vo$gam1$coefficients,vo$gam2$coefficients, vo$log.sig2.1, vo$log.sig2.2,                            vo$i.rho)
     
   }
   
   start.v
   
 }
 
 start.v <- overall.svEXPANDED(margins, M, vo)
   			
##############################################################  
# starting values for case of predictors on all parameters
############################################################## 
 
 
 overall.svGEXPANDED <- function(formula, data, ngc, margins, M, vo, gam1, gam2, type = "copR", inde = NULL, c.gam2 = NULL, gam3 = NULL, knots = NULL){
   
   X3 = X4 = X5 = X6 = X7 = X8 = X3.d2 = X4.d2 = X5.d2 = X6.d2 = X7.d2 = X8.d2 = NULL
   gp3 = gp4 = gp5 = gp6 = gp7 = gp8 = NULL
   gam4 = gam5 = gam6 = gam7 = gam8 = NULL
   l.sp3 = l.sp4 = l.sp5 = l.sp6 = l.sp7 = l.sp8 = 0  
   sp3 = sp4 = sp5 = sp6 = sp7 = sp8 = NULL  
   X3s = X4s = X5s = NULL
   Sl.sf2 <- Sl.sf3 <- NULL
   
   
   #print(type != "triv")
   if(type != "triv") gam3 <- NULL    
   

   #print(type == "copR")
   if(type == "copR"){   
     
     BivD <- M$BivD
     

     
       
       #print(margins[1] %in% c(M$m2,M$m2d) && margins[2] %in% c(M$m2,M$m2d))
       if(margins[1] %in% c(M$m2,M$m2d) && margins[2] %in% c(M$m2,M$m2d)){
         
         formula.eq3 <- formula[[3]] 
         formula.eq4 <- formula[[4]] 
         formula.eq5 <- formula[[5]]     
         nad2.1 <- "sigma2.1" 
         nad2.2 <- "sigma2.2" 
         nad3   <- "theta"     
         formula.eq3 <- as.formula( paste(nad2.1,"~",formula.eq3[2],sep="") ) 
         formula.eq4 <- as.formula( paste(nad2.2,"~",formula.eq4[2],sep="") ) 
         formula.eq5 <- as.formula( paste(  nad3,"~",formula.eq5[2],sep="") )   
         
         set.seed(1)
         sigma2.1 <- rnorm(vo$n, vo$log.sig2.1, 0.001) 
         sigma2.2 <- rnorm(vo$n, vo$log.sig2.2, 0.001) 
         theta    <- rnorm(vo$n, vo$i.rho, 0.001)      
         rm(list=".Random.seed", envir=globalenv()) 
         
         gam3 <- gam(formula.eq3, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels) 
         gam4 <- gam(formula.eq4, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)   
         gam5 <- gam(formula.eq5, data = data, gamma = ngc, knots = knots, drop.unused.levels = vo$drop.unused.levels)     
         l.sp3 <- length(gam3$sp)   
         l.sp4 <- length(gam4$sp)    
         l.sp5 <- length(gam5$sp)  
         
         
         #print(l.sp3 != 0)
         if(l.sp3 != 0){
           ngc <- 2
           while( any(round(summary(gam3)$edf, 1) > 1) ) {gam3 <- gam(formula.eq3, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
         } 
         
         
         #print(l.sp4 != 0)
         if(l.sp4 != 0){
           ngc <- 2
           while( any(round(summary(gam4)$edf, 1) > 1) ) {gam4 <- gam(formula.eq4, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
         }                    
         
         
         #print(l.sp5 != 0)
         if(l.sp5 != 0){
           ngc <- 2
           while( any(round(summary(gam5)$edf, 1) > 1) ) {gam5 <- gam(formula.eq5, data = data, gamma = ngc + 1, knots = knots, drop.unused.levels = vo$drop.unused.levels); ngc <- ngc + 1; if(ngc > 5) break}  
         }                    
         
         
         
         X3 <- model.matrix(gam3)
         X3.d2 <- dim(X3)[2]
         X4 <- model.matrix(gam4)
         X4.d2 <- dim(X4)[2]  
         X5 <- model.matrix(gam5)
         X5.d2 <- dim(X5)[2]     
         
         #print(l.sp3 != 0)
         if(l.sp3 != 0) sp3 <- gam3$sp 
         environment(gam3$formula) <- environment(gam2$formula)
         gp3 <- gam3$nsdf 
         
         
         #print(l.sp4 != 0)
         if(l.sp4 != 0) sp4 <- gam4$sp 
         environment(gam4$formula) <- environment(gam2$formula)
         gp4 <- gam4$nsdf     
         
         
         #print(l.sp5 != 0)
         if(l.sp5 != 0) sp5 <- gam5$sp 
         environment(gam5$formula) <- environment(gam2$formula)
         gp5 <- gam5$nsdf      
         
         start.v  <- c(gam1$coefficients, gam2$coefficients, gam3$coefficients, gam4$coefficients, gam5$coefficients )
         
       }   
     
   }
   
   L <- list(start.v = start.v,
             X3 = X3, X4 = X4, X5 = X5, X6 = X6, X7 = X7, X8 = X8, X3.d2 = X3.d2, X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2,
             X7.d2 =	X7.d2, X8.d2 =	X8.d2, gp3 = gp3, gp4 = gp4, gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8,
             gam3 = gam3, gam4 = gam4, gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8,
             l.sp3 = l.sp3, l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, l.sp8 = l.sp8,
             sp3 = sp3, sp4 = sp4, sp5 = sp5, sp6 = sp6, sp7 = sp7, sp8 = sp8, X3s = X3s, X4s = X4s, X5s = X5s,  Sl.sf2 =  Sl.sf2,  Sl.sf3= Sl.sf3)
   
   
   L
   
 }
 
 
 
    #print(l.flist > 2)
 if(l.flist > 2){
    
    overall.svGR <- overall.svGEXPANDED(formula, data, ngc = 2, margins, M, vo, gam1, gam2, knots = knots)
    
    start.v = overall.svGR$start.v 
    X3 = overall.svGR$X3; X4 = overall.svGR$X4; X5 = overall.svGR$X5
    X6 = overall.svGR$X6; X7 = overall.svGR$X7; X8 = overall.svGR$X8
    X3.d2 = overall.svGR$X3.d2; X4.d2 = overall.svGR$X4.d2; X5.d2 = overall.svGR$X5.d2
    X6.d2 = overall.svGR$X6.d2; X7.d2 = overall.svGR$X7.d2; X8.d2 = overall.svGR$X8.d2
    gp3 = overall.svGR$gp3; gp4 = overall.svGR$gp4; gp5 = overall.svGR$gp5
    gp6 = overall.svGR$gp6; gp7 = overall.svGR$gp7; gp8 = overall.svGR$gp8
    gam3 = overall.svGR$gam3; gam4 = overall.svGR$gam4; gam5 = overall.svGR$gam5
    gam6 = overall.svGR$gam6; gam7 = overall.svGR$gam7; gam8 = overall.svGR$gam8
    l.sp3 = overall.svGR$l.sp3; l.sp4 = overall.svGR$l.sp4; l.sp5 = overall.svGR$l.sp5
    l.sp6 = overall.svGR$l.sp6; l.sp7 = overall.svGR$l.sp7; l.sp8 = overall.svGR$l.sp8
    sp3 = overall.svGR$sp3; sp4 = overall.svGR$sp4; sp5 = overall.svGR$sp5
    sp6 = overall.svGR$sp6; sp7 = overall.svGR$sp7; sp8 = overall.svGR$sp8
    
    }
    

##########################################################
# SPs and penalties
##########################################################
  

 GAM <- list(gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, 
            gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8)   


#print((l.sp1!=0 || l.sp2!=0 || l.sp3!=0 || l.sp4!=0 || l.sp5!=0 || l.sp6!=0 || l.sp7!=0 || l.sp8!=0) && fp==FALSE )
 if( (l.sp1!=0 || l.sp2!=0 || l.sp3!=0 || l.sp4!=0 || l.sp5!=0 || l.sp6!=0 || l.sp7!=0 || l.sp8!=0) && fp==FALSE ){ 

 L.GAM <- list(l.gam1 = length(gam1$coefficients), l.gam2 = length(gam2$coefficients), l.gam3 = length(gam3$coefficients), l.gam4 = length(gam4$coefficients),
              l.gam5 = length(gam5$coefficients), l.gam6 = length(gam6$coefficients), l.gam7 = length(gam7$coefficients), l.gam8 = length(gam8$coefficients))

 L.SP <- list(l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, l.sp4 = l.sp4, 
             l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, l.sp8 = l.sp8)

                 sp <- c(sp1, sp2, sp3, sp4, sp5, sp6, sp7, sp8)
                 qu.mag <- S.m(GAM, L.SP, L.GAM)                             
                                                        }
  

##########################################################
# general lists
##########################################################

#print(missing(parscale))
 if(missing(parscale)) parscale <- 1   

  respvec <- respvec2 <- respvec3 <- list(y1 = y1, y2 = y2,
                                          y1.y2 = NULL, y1.cy2 = NULL, 
                                          cy1.y2 = NULL, cy1.cy2 = NULL, 
                                          cy1 = NULL, cy = NULL, univ = 0)
 
  my.env <- new.env()
  my.env$signind <- 1 # this is for mixed copulae

  lsgam3 <- length(gam3$smooth)
  lsgam4 <- length(gam4$smooth)
  lsgam5 <- length(gam5$smooth)
  lsgam6 <- length(gam6$smooth)
  lsgam7 <- length(gam7$smooth)
  lsgam8 <- length(gam8$smooth)




 #my.env      <- new.env()
 my.env$k1   <- k1.tvc
 my.env$k2   <- k2.tvc


 VC <- list(lsgam1 = lsgam1, indexTeq1 = indexTeq1, indexTeq2 = indexTeq2, 
             lsgam2 = lsgam2, Deq1 = Deq1, pos.pbeq1 = pos.pbeq1, Deq2 = Deq2, pos.pbeq2 = pos.pbeq2,
             lsgam3 = lsgam3, robust = FALSE, sp.fixed = NULL,
             lsgam4 = lsgam4, Sl.sf = Sl.sf, sp.method = sp.method,
             lsgam5 = lsgam5, K1 = NULL,
             lsgam6 = lsgam6, 
             lsgam7 = lsgam7,
             lsgam8 = lsgam8,
             X1 = X1, 
             X2 = X2, 
             X3 = X3,
             X4 = X4, 
             X5 = X5, 
             X6 = X6,  
             X7 = X7, 
             X8 = X8,
             X1.d2 = X1.d2, 
             X2.d2 = X2.d2,
             X3.d2 = X3.d2,
             X4.d2 = X4.d2,
             X5.d2 = X5.d2,
             X6.d2 = X6.d2,
             X7.d2 = X7.d2,
             X8.d2 = X8.d2,
             gp1 = gp1, 
             gp2 = gp2,
             gp3 = gp3,
             gp4 = gp4, 
             gp5 = gp5,
             gp6 = gp6, 
             gp7 = gp7,
             gp8 = gp8,
             l.sp1 = l.sp1, 
             l.sp2 = l.sp2,
             l.sp3 = l.sp3,
             l.sp4 = l.sp4, 
             l.sp5 = l.sp5,
             l.sp6 = l.sp6, 
             l.sp7 = l.sp7, 
             l.sp8 = l.sp8, my.env = my.env,
             infl.fac = infl.fac,
             weights = weights,
             fp = fp, 
             gamlssfit = gamlssfit,
             hess = NULL,
             Model = "CC", univ.gamls = FALSE,
             end = end,
             BivD = BivD, nCa = nCa,
             nC = nC, gc.l = gc.l, 
             n = n, extra.regI = extra.regI,
             parscale = parscale, margins = margins,
             Cont = "YES", ccss = "no", m2 = m2, m3 = m3, 
             m1d = m1d, m2d = m2d, m3d = m3d, 
             bl = bl, triv = FALSE,
             y1m = y1m, y2m = y2m, 
             tc = t.c,
             i.rho = i.rho, dof = dof,
             dof.st = dof.st, BivD2 = BivD2, cta = cta, ct = ct,
             zerov = -10,
             c11 = c11,
             c10 = c10,
             c01 = c01,
             c00 = c00, surv = surv,
             Xd1 = Xd1, Xd2 = Xd2,
             mono.sm.pos1 = mono.sm.pos1, mono.sm.pos2 = mono.sm.pos2, 
             surv.flex = surv.flex,
             mono.sm.pos = mono.sm.pos, gp2.inf = NULL,
             informative = "no",
             zero.tol = 1e-02,
             min.dn = min.dn, min.pr = min.pr, max.pr = max.pr) # original n only needed in SemiParBIV.fit
  
  
         
             
  ##########################################################################################################################
  ##########################################################################################################################
  # GAMLSS fit
  ##########################################################################################################################
  ##########################################################################################################################

 form.gamlEXPANDED <- function(formula, l.flist, M, type = "copR"){
   
   formula.gamlss1 <- formula.gamlss2 <- NULL
   
   if(type == "copR"){
     
     #print(l.flist > 2)
     if(l.flist > 2){##
       
       #print(M$margins[1] %in% c(M$m2,M$m2d) && M$margins[2] %in% c(M$m2,M$m2d))
       if(M$margins[1] %in% c(M$m2,M$m2d) && M$margins[2] %in% c(M$m2,M$m2d)){
         
         formula.gamlss1 <- list(formula[[1]],formula[[3]])
         formula.gamlss2 <- list(formula[[2]],formula[[4]])      
         
       }   
       
       
     } ##
     
     
   }
   
   
   list(formula.gamlss1 = formula.gamlss1, formula.gamlss2 = formula.gamlss2)
   
 }
 

#print(gamlssfit == TRUE)
 if(gamlssfit == TRUE){ 

  form.gamlR <- form.gamlEXPANDED(formula, l.flist, M)

  surv1 <- surv2 <- surv
  
  
  

  
  
  

  gamlss1 <- eval(substitute(gamlss(form.gamlR$formula.gamlss1, data = data, weights = weights, subset = subset,  
                   margin = margins[1], surv = surv1, cens = cens1, infl.fac = infl.fac, 
                   rinit = rinit, rmax = rmax, iterlimsp = iterlimsp, tolsp = tolsp,
                   gc.l = gc.l, parscale = 1, extra.regI = extra.regI, k.tvc = k1.tvc, drop.unused.levels = drop.unused.levels), list(weights=weights,cens1=cens1)))

  gamlss2 <- eval(substitute(gamlss(form.gamlR$formula.gamlss2, data = data, weights = weights, subset = subset,  
                   margin = margins[2], surv = surv2, cens = cens2, infl.fac = infl.fac, 
                   rinit = rinit, rmax = rmax, iterlimsp = iterlimsp, tolsp = tolsp,
                   gc.l = gc.l, parscale = 1, extra.regI = extra.regI, k.tvc = k2.tvc, drop.unused.levels = drop.unused.levels), list(weights=weights,cens2=cens2)))   
                      
  # updated starting values   

  SP <- list(sp1 = sp1, sp2 = sp2, sp3 = sp3, sp4 = sp4, sp5 = sp5, sp6 = sp6, sp7 = sp7, sp8 = sp8)
  gamls.upsvR <- gamls.upsv(gamlss1, gamlss2, margins, M, l.flist, nstv = names(start.v), VC, GAM, SP)
  sp <- gamls.upsvR$sp
  start.v <- gamls.upsvR$start.v 

}

  ##########################################################################################################################
  ##########################################################################################################################
  
  func.opt <- func.OPT(margins, M)  
  
  SemiParFit <- SemiParBIV.fit(func.opt = func.opt, start.v = start.v, 
                         rinit = rinit, rmax = rmax, iterlim = 100, iterlimsp = iterlimsp, tolsp = tolsp,
                         respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag) 
    
  ##########################################################################################################################
  # post estimation
  ##########################################################################################################################

  SemiParFit.p <- copulaReg.fit.post(SemiParFit = SemiParFit, VC = VC, GAM)
 
  y1.m <- y1 
  y2.m <- y2

  SemiParFit <- SemiParFit.p$SemiParFit

  ##########################################################################################################################

 cov.c(SemiParFit)

  ##########################################################################################################################
  
 gam1$call$data <- gam2$call$data <- gam3$call$data <- gam4$call$data <- gam5$call$data <- gam6$call$data <- gam7$call$data <- gam8$call$data <- cl$data 
  
  # for all.terms
  ##########################################################################################################################


 L <- list(fit = SemiParFit$fit, dataset = NULL, n = n, gamlss1 = gamlss1, gamlss2 = gamlss2, formula = formula, robust = FALSE,     
          edf11 = SemiParFit.p$edf11, surv = surv, 
          gam1 = gam1, gam2 = gam2, gam3 = gam3, gam4 = gam4, gam5 = gam5, gam6 = gam6, gam7 = gam7, gam8 = gam8,  
          coefficients = SemiParFit$fit$argument, coef.t = SemiParFit.p$coef.t, 
          iterlimsp = iterlimsp,
          weights = weights, cens1 = cens1, cens2 = cens2, cens3 = cens3,
          sp = SemiParFit.p$sp, iter.sp = SemiParFit$iter.sp, 
          l.sp1 = l.sp1, l.sp2 = l.sp2, l.sp3 = l.sp3, 
          l.sp4 = l.sp4, l.sp5 = l.sp5, l.sp6 = l.sp6, l.sp7 = l.sp7, l.sp8 = l.sp8, bl = bl,
          fp = fp,  
          iter.if = SemiParFit$iter.if, iter.inner = SemiParFit$iter.inner,
          theta = SemiParFit.p$theta, 
          theta.a = SemiParFit.p$theta.a,  
          sigma21 = SemiParFit.p$sigma21, sigma22 = SemiParFit.p$sigma22, 
          sigma21.a = SemiParFit.p$sigma21.a, sigma22.a = SemiParFit.p$sigma22.a,
          sigma1 = SemiParFit.p$sigma21, sigma2 = SemiParFit.p$sigma22, 
          sigma1.a = SemiParFit.p$sigma21.a, sigma2.a = SemiParFit.p$sigma22.a,                    
          nu1 = SemiParFit.p$nu1, nu2 = SemiParFit.p$nu2, 
          nu1.a = SemiParFit.p$nu1.a, nu2.a = SemiParFit.p$nu2.a,
          dof.a = SemiParFit.p$dof.a, dof = SemiParFit.p$dof,
          X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, X6 = X6, X7 = X7, X8 = X8,
          X1.d2 = X1.d2, X2.d2 = X2.d2, X3.d2 = X3.d2, 
          X4.d2 = X4.d2, X5.d2 = X5.d2, X6.d2 = X6.d2, X7.d2 = X7.d2, X8.d2 = X8.d2,            
          He = SemiParFit.p$He, HeSh = SemiParFit.p$HeSh, Vb = SemiParFit.p$Vb, Ve = SemiParFit.p$Ve, 
          F = SemiParFit.p$F, F1 = SemiParFit.p$F1,  
          t.edf = SemiParFit.p$t.edf, edf = SemiParFit.p$edf, 
          edf1 = SemiParFit.p$edf1, edf2 = SemiParFit.p$edf2, edf3 = SemiParFit.p$edf3,
          edf4 = SemiParFit.p$edf4, edf5 = SemiParFit.p$edf5, edf6 = SemiParFit.p$edf6, edf7 = SemiParFit.p$edf7,
          edf8 = SemiParFit.p$edf8,
          edf1.1 = SemiParFit.p$edf1.1, edf1.2 = SemiParFit.p$edf1.2, edf1.3 = SemiParFit.p$edf1.3,
          edf1.4 = SemiParFit.p$edf1.4, edf1.5 = SemiParFit.p$edf1.5, edf1.6 = SemiParFit.p$edf1.6, edf1.7 = SemiParFit.p$edf1.7, 
          edf1.8 = SemiParFit.p$edf1.8, 
          R = SemiParFit.p$R,
          bs.mgfit = SemiParFit$bs.mgfit, conv.sp = SemiParFit$conv.sp, 
          wor.c = SemiParFit$wor.c,  
          eta1 = SemiParFit$fit$eta1, eta2 = SemiParFit$fit$eta2, 
          etad=SemiParFit$fit$etad, etas1 = SemiParFit$fit$etas1, etas2 = SemiParFit$fit$etas2,
          y1 = y1.m, y2 = y2.m, 
          BivD = BivD, margins = margins,   
          logLik = SemiParFit.p$logLik,
          nC = nC, 
          respvec = respvec, hess = TRUE,
          qu.mag = qu.mag, 
          gp1 = gp1, gp2 = gp2, gp3 = gp3, gp4 = gp4, gp5 = gp5, gp6 = gp6, gp7 = gp7, gp8 = gp8, 
          VC = VC, magpp = SemiParFit$magpp,
          gamlssfit = gamlssfit, Cont = "YES",
          tau = SemiParFit.p$tau, tau.a = SemiParFit.p$tau.a, l.flist = l.flist, v1 = v1, v2 = v2, triv = FALSE, univar.gamlss = FALSE,
          BivD2 = BivD2, call = cl, surv = surv, surv.flex = surv.flex,
          Vb.t = SemiParFit.p$Vb.t, coef.t = SemiParFit.p$coef.t)
  

class(L) <- c("gjrm","SemiParBIV")


}

L


}

