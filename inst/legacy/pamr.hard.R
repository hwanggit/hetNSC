pamr.train.hard <-
function(data, gene.subset=NULL, sample.subset=NULL,
         threshold = NULL, n.threshold = 30,
        scale.sd = TRUE, threshold.scale = NULL, se.scale = NULL, offset.percent = 50, hetero=NULL,
         prior = NULL,  remove.zeros = TRUE, sign.contrast="both", ngroup.survival=2)

{
  
# modified aug 2003 to add survival analysis facilities
#
# possibilities for different outcomes in data object "data":
#         y= class variable => classification problem
#         proby= matrix of class probabilities => "soft classification"
#              (not currently used by Excel interface)
#        survival time, censoring status present => survival analysis problem
#
# here is how the two  problem types  are passed to nsc:
#
#       class: y is class label, proby, prob.ytest not used
#       surv.km: proby are soft class labels, computed from kaplan-meier
  
        this.call <- match.call()

 if(!is.null(data$y)){problem.type <- "class"}
  if(!is.null(data$survival.time)) {problem.type <- "surv.km"}

        
          if(!is.null(data$proby) & !is.null(data$y)) {
           stop("Can't have both proby and y present in data object")
         }
        
          if(!is.null(data$y) & !is.null(data$survival.time)) {
    stop("Can't have both class label y and survival.time present in data object")
  }
        
        if(!is.null(data$y) & !is.null(data$censoring.status)) {
           stop("Can't have both class label y and censoring status present in data object")
         }
         if(!is.null(data$survival.time) & is.null(data$censoring.status)) {
           stop("Survival time specified but censoring status missing")
         }
          if(is.null(data$survival.time) & !is.null(data$censoring.status)) {
           stop("Censoring status specified but survival time missing")
         }
  
        
        y <- data$y
        proby <- data$proby
        ytest <- NULL
        xtest <- NULL
        prob.ytest <- NULL
       
        
        cutoffs.survival <- NULL
        
#       estimate class probabilities via Kaplan-Meier
#        use Cox score cutoff of 2.0 or 20th largest score, whichever is smaller
#       this ensures at least 20 genes are used for clustering
        
        if(!is.null(data$survival.time)){
          junk <- pamr.surv.to.class2(data$survival.time, data$censoring.status,
             n.class=ngroup.survival)
                 
            proby <- junk$prob
            cutoffs.survival <- junk$cutoffs
        }

        
# ytemp is just used for computation of the prior

 if(!is.null(y)){ ytemp <- y}
 if(is.null(y) & !is.null(proby)){ ytemp <- apply(proby,1, which.is.max)}
 if(is.null(sample.subset)){sample.subset <-1:ncol(data$x)}
 if(is.null(gene.subset)){gene.subset <-1:nrow(data$x)}

# for survival analysis, make default prior the equal prior

  if(is.null(prior) & !is.null(data$survival.time) ){
         prior <- rep(1/ngroup.survival, ngroup.survival)
}

        if(is.null(prior) & is.null(data$survival.time) )
          {prior <- table(ytemp[sample.subset])/length(ytemp[sample.subset])
           prior <- prior[prior!=0]
        }
      
    if(!is.null(y)){y <-  factor(y[sample.subset])}
        
    if(!is.null(proby)){
        proby <- proby[sample.subset,]}
        junk <- nsc.hard(data$x[gene.subset, sample.subset], y=y, proby=proby,
 xtest=xtest, ytest=ytest, prob.ytest=prob.ytest,
          offset.percent=offset.percent,  threshold = threshold, hetero=hetero,
        n.threshold = n.threshold,  scale.sd= scale.sd,
                    threshold.scale=threshold.scale,
           se.scale= se.scale, prior=prior, remove.zeros=remove.zeros,
            sign.contrast=sign.contrast, problem.type=problem.type)

        junk$call <- this.call
        junk$gene.subset <- gene.subset
        junk$sample.subset <- sample.subset
        junk$survival.time <- data$survival.time
        junk$censoring.status <- data$censoring.status
       junk$cutoffs.survival <- cutoffs.survival
        junk$ngroup.survival <- ngroup.survival
        junk$problem.type <- problem.type
        class(junk)="pamrtrained"
        junk
}


nsc.hard <-
  
  function(x, y = NULL, xtest = NULL, proby = NULL, ytest = NULL, prob.ytest = 
        NULL, threshold = NULL, n.threshold = 30, hetero = NULL, scale.sd = 
        TRUE, threshold.scale = NULL, se.scale = NULL, offset.percent = 50, 
        prior = table(y)/length(y), remove.zeros = TRUE, sign.contrast = "both",
           problem.type=c("class", "surv.km","surv.latent"))
{

# modified aug 2003 to add survival analysis facilities
#
#         problem.type can be "class", "surv.km", or "surv.latent"         
#
#         y= class variable => classification problem (problem.type="class")
#         proby= matrix of class probabilities => "soft classification"
#               from Kaplan-Meier estimate, for survival analysis
#         in this case, nsc computes probability-weighted centroids and
#             training error 

        this.call <- match.call()

        argy <- ytest
        if(is.null(ytest)) {
                argy <- y
        }
        if(!is.null(y) & !is.null(proby) & problem.type!="surv.latent") {
                stop("Can't specify both y and proby")
        }
        if(!is.null(ytest) & !is.null(prob.ytest)) {
                stop("Can't specify both ytest and prob.ytest")
        }
        if(is.null(y)) {
                y <- apply(proby, 1, which.is.max)
        }
        n.class <- table(y)
        if(min(n.class) == 1) {
                cat("Warning: a class contains only 1 sample")
        }
        if(is.null(xtest)) {
                xtest <- x
                ytest <- y
                prob.ytest <- proby
        }
        norm.cent <- NULL
        if(!is.null(hetero)) {
                norm.cent <- apply(x[, y == hetero], 1, mean)
                x <- abs(t(scale(t(x), center = norm.cent, scale = FALSE)))
                if(!missing(xtest)) {
                        xtest <- abs(t(scale(t(xtest), center = norm.cent, 
                                scale = FALSE)))
                }
        }
        n <- sum(n.class)
        ntest <- ncol(xtest)
        K <- length(prior)
        p <- nrow(x)
        if(is.null(proby)) {
                Y <- model.matrix( ~ factor(y) - 1, data = list(y = y))
        }
        if(!is.null(proby)) {
                Y <- proby
        }
        dimnames(Y) <- list(NULL, names(n.class))

        centroids <- scale(x %*% Y, FALSE, n.class)
        sd <- rep(1, p)
        if(scale.sd) {
                xdif <- x - centroids %*% t(Y)
                sd <- (xdif^2) %*% rep(1/(n - K), n)
                sd <- drop(sqrt(sd))
                offset <- quantile(sd, offset.percent/100)
                sd <- sd + offset
        }
        centroid.overall <- drop(x %*% rep(1/n, n))
        if(is.null(threshold.scale)) {
                threshold.scale <- rep(1, K)
                names(threshold.scale) <- names(n.class)
        }
### Now make an adjustment for the sample sizes in the "t" ratios

        if(is.null(se.scale))
                se.scale <- sqrt(1/n.class - 1/n)
        delta <- (centroids - centroid.overall)/sd
        delta <- scale(delta, FALSE, threshold.scale * se.scale)
        if(sign.contrast == "positive") {
                delta <- delta * (delta > 0)
        }
        if(sign.contrast == "negative") {
                delta <- delta * (delta < 0)
        }
#allows differential shrinkage
        if(!is.null(threshold)) {
                n.threshold <- length(threshold)
        }
        else {
                #threshold <- seq((sqrt(2*log(length(delta)))-1)*sd(c(delta)), max(abs(delta)), length = n.threshold)
                #threshold <- quantile(abs(delta), probs = seq(0, 1, length = n.threshold))
				threshold <- seq(0, max(abs(delta)), length = n.threshold)
        }
        nonzero <- seq(n.threshold)
        errors <- threshold
        yhat <- as.list(seq(n.threshold))
        prob <- array(0, c(ntest, K, n.threshold))
        for(ii in 1:n.threshold) {
                cat(ii)
                delta.shrunk <- hard.shrink(delta, threshold[ii])
                delta.shrunk <- scale(delta.shrunk, FALSE, 1/(threshold.scale * 
                        se.scale))
                nonzero[ii] <- attr(delta.shrunk, "nonzero")
                posid <- drop(abs(delta.shrunk) %*% rep(1, K)) > 0
                dd <- diag.disc((xtest - centroid.overall)/sd, delta.shrunk, 
                        prior, weight = posid)
                yhat[[ii]] <- softmax(dd)
                dd <- safe.exp(dd)
                prob[,  , ii] <- dd/drop(dd %*% rep(1, K))
                if(!is.null(ytest)) {
                        errors[ii] <- sum(yhat[[ii]] != ytest)
                }
                if(!is.null(prob.ytest)) {

# use of temp below is to ensure that Yhat doesn;t drop a column
#  when no predictions are made to that class
                        temp <- c(yhat[[ii]], names(table(y)))
                        Yhat <- model.matrix( ~ factor(temp) - 1, data = list(y
                                 = temp))
                        Yhat <- Yhat[1:length(yhat[[ii]]),  ]
                     
                        errors[ii] <- length(yhat[[ii]]) - sum(Yhat * prob.ytest)
                }
               
        }
        thresh.names <- format(round(threshold, 3))
        names(yhat) <- thresh.names
        attr(yhat, "row.names") <- paste(seq(ntest))
        class(yhat) <- "data.frame"
        if(remove.zeros)
                n.threshold <- match(0, nonzero, n.threshold)
        dimnames(prob) <- list(paste(seq(ntest)), names(n.class), thresh.names)
        object <- list(y = argy, proby = prob.ytest, yhat = yhat[, seq(
                n.threshold)], prob = prob[,  , seq(n.threshold)], centroids = 
                centroids, centroid.overall = centroid.overall, sd = sd, 
                threshold = threshold[seq(n.threshold)], nonzero = nonzero[seq(
                n.threshold)], delta=delta, threshold.scale = threshold.scale, se.scale = 
                se.scale, scale.sd=scale.sd, call = this.call, hetero = hetero, norm.cent = 
                norm.cent, prior = prior, offset = offset, sign.contrast = 
                sign.contrast)
        if(!is.null(ytest) | !is.null(prob.ytest))
                object$errors <- errors[seq(n.threshold)]
        class(object) <- "nsc"
        object
}


safe.exp=function(x){
 xx=sign(x)*pmin(abs(x),500)
 return(exp(xx))
}


# The x argument of softmax is a matrix of n_t by K,
# where n_t is the number of observations to be classfied. 
# Each row of x is the vector of the discriminant function ( \delta_1(x^*), ...,  \delta_K(x^*)  )'
# The softmax function assigns each observation to the class which has the smallest discriminant function
#
#  The article equation [7] assign the observation to the class which has the 
#        smallest discriminant function
softmax <-function(x, gap = FALSE) {
  d <- dim(x)
  maxdist <- x[, 1]
  pclass <- rep(1, d[1])
  for(i in seq(2, d[2])) {
    l <- x[, i] > maxdist
    pclass[l] <- i
    maxdist[l] <- x[l, i]
  }  #above loop finds the maximum value in each row and record which column achieved the maximum
     # This loop can be replaced with the following two lines
     # pclass =apply(x, 1, which.max)
	 # maxdist = apply(x, 1, max)

  dd <- dimnames(x)[[2]]
  if(gap) {
    x <- abs(maxdist - x)
    x[cbind(seq(d[1]), pclass)] <- drop(x %*% rep(1, d[2]))
    gaps <- do.call("pmin", data.frame(x))  # this is same as apply(x, 1, pmin)
  }
  pclass <- if(is.null(dd) || !length(dd))
    pclass
  else
    factor(pclass, levels = seq(d[2]), labels = dd)
  if(gap)
    list(class = pclass, gaps = gaps)
  else
    pclass
}


hard.shrink <-function(delta, threshold) {
  dif <- abs(delta) - threshold
  delta <- delta * (dif > 0)
  nonzero <- sum(drop((dif > 0) %*% rep(1, ncol(delta))) > 0)
  attr(delta, "nonzero") <- nonzero
  delta
}


pamr.cv.hard <-

function(fit, data, nfold = NULL, folds = NULL ,...)
{
        x <- data$x[fit$gene.subset, fit$sample.subset]

        if( !is.null(data$y) & !is.null(data$proby)){
           stop("Must have exactly one of y and  proby  present in the data object")
         }
        
        y <- NULL
        proby <- NULL
        
        if(!is.null(fit$y)){
           y<-  factor(fit$y[fit$sample.subset])
         }
        
        if(!is.null(fit$proby)){
           proby<-  fit$proby[fit$sample.subset,]
         }
        
        this.call <- match.call()
        
# three possibilities, 
# problem.type= class: y are class labels, proby=NULL
#               surv.km: y=NULL, proby are soft class probabilities from KM
#               surv.latent: y are latent class labels,
#                   proby are soft class probabilities from KM
# note; problem type is in fit$problem.type
        
        junk <- nsccv.hard(x, y=y, proby=proby, object = fit, nfold=nfold, folds=folds, 
survival.time=data$survival.time, censoring.status = data$censoring.status, 
ngroup.survival=fit$ngroup.survival, problem.type=fit$problem.type, ...)

        junk$call <- this.call
        
        junk$sample.subset <- fit$sample.subset
        class(junk)="pamrcved"
        junk
}

nsccv.hard <- function(x, y=NULL, proby=NULL, nfold = min(table(y)), folds = NULL, threshold =
        NULL, threshold.scale = NULL, survival.time=NULL, censoring.status=NULL, ngroup.survival=NULL,prior, object, ...)
{
        this.call <- match.call()

        argy <- y
        
#         if( !is.null(y) & !is.null(proby)){
#           stop("Must have at most one of y and  proby  present in the data object")
#         }

        if(is.null(y)){ y <- as.factor(apply(proby,1,which.is.max))}
        
        n <- length(y)

if(is.null(nfold) & is.null(survival.time)) {nfold <- min(table(y))}
if(is.null(nfold) & !is.null(survival.time)) {nfold <- 10}


 if(is.null(survival.time)){
        if(is.null(folds)) {
                folds <-balanced.folds(y)
        }
       }


        if(!is.null(survival.time)){
        if(is.null(folds)) {
                folds <- split(sample(1:n), rep(1:nfold, length = n))
        }
        }
         
nfold<- length(folds)

        if(missing(prior)) {
                if(missing(object))
                        prior <- table(y)/n
                else prior <- object$prior
        }
    
        if(missing(threshold)) {
                if(missing(object))
                        stop("Must either supply threshold argument, or an nsc object"
                                )
                else {
                        threshold <- object$threshold
                        threshold.scale <- object$threshold.scale
                        se.scale <- object$se.scale
                }
        }
       
        n.threshold <- length(threshold)        ### Set up the data structures
        yhat <- rep(list(y), n.threshold)
        names(yhat) <- paste(seq(n.threshold))
        yhat <- data.frame(yhat)
        n.class <- table(y)
        prob <- array(1, c(n, length(n.class), n.threshold))
        size <- double(n.threshold)
        hetero <-object$hetero
        cv.objects=vector("list",nfold)
        for(ii in 1:nfold) {
                cat("Fold", ii, ":")
                a <- nsc.hard(x[,  - folds[[ii]]], y=argy[ - folds[[ii]]], x[, folds[[ii
                        ]], drop = FALSE], proby=proby[-folds[[ii]],],
                         threshold = threshold, threshold.scale
                         = threshold.scale, se.scale = se.scale, prior = prior,
                          hetero=hetero,
                        ..., remove.zeros = FALSE)
                size <- size + a$nonzero
                prob[folds[[ii]],  ,  ] <- a$prob
                yhat[folds[[ii]],  ] <- a$yhat
                cat("\n")
        cv.objects[[ii]]=a
        }
        if(missing(object))
                size <- round(size/nfold)
        else size <- object$nonzero
        error <- rep(NA, n.threshold)
        loglik <- error
        pvalue.survival <- error
        
        pvalue.survival.func <- function(group, survival.time, censoring.status,ngroup.survival){
            temp <- coxph(Surv(survival.time, censoring.status)~as.factor(group))
            loglik <- 2*(temp$loglik[2]-temp$loglik[1])
            return(1-pchisq(loglik, ngroup.survival-1))
          }
        
        if(!is.null(proby)){proby.temp <-proby}
        else if(!is.null(survival.time)){proby.temp <- pamr.surv.to.class2(survival.time,
                                       censoring.status, n.class=ngroup.survival)$prob
                                       }
        
        for(i in 1:n.threshold) {
      
                if(is.null(survival.time) & is.null(proby)){error[i] <- sum(yhat[, i] != y)/n}
                if(!is.null(survival.time)){
                    
                    temp <- c(yhat[,i],names(table(y)))
                    Yhat <- model.matrix( ~ factor(temp) - 1,
                                       data = list(y = temp))
                     Yhat <- Yhat[1:length(yhat[[ii]]),]
                     error[i] <- (length(yhat[,i])-sum(Yhat*proby.temp))/n
                  }
            
                
                if(is.null(survival.time)){
                  loglik[i] <- sum(log(prob[,  , i][cbind(seq(1, n), unclass(y))]))/                        n}
                
                if(!is.null(survival.time)){
                  pvalue.survival[i]<- pvalue.survival.func(yhat[,i], survival.time,censoring.status, ngroup.survival)
                }
        }

obj<- list(threshold=threshold, error=error, loglik=loglik,size=size, yhat=yhat,y=y,prob=prob,folds=folds, cv.objects=cv.objects, pvalue.survival=pvalue.survival,
                call = this.call)
        class(obj) <- "nsccv"
        obj
}

 balanced.folds <- function(y, nfolds = min(min(table(y)), 10)) {
   totals <- table(y)
   fmax <- max(totals)
   nfolds <- min(nfolds, fmax)     
   nfolds= max(nfolds, 2)
         # makes no sense to have more folds than the max class size
   folds <- as.list(seq(nfolds))
   yids <- split(seq(y), y) 
         # nice we to get the ids in a list, split by class
###Make a big matrix, with enough rows to get in all the folds per class
   bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
   for(i in seq(totals)) {
cat(i)
     if(length(yids[[i]])>1){bigmat[seq(totals[i]), i] <- sample(yids[[i]])}
     if(length(yids[[i]])==1){bigmat[seq(totals[i]), i] <- yids[[i]]}

   }
   smallmat <- matrix(bigmat, nrow = nfolds)# reshape the matrix
### Now do a clever sort to mix up the NAs
   smallmat <- permute.rows(t(smallmat))   ### Now a clever unlisting
         # the "clever" unlist doesn't work when there are no NAs
         #       apply(smallmat, 2, function(x)
         #        x[!is.na(x)])
   res <-vector("list", nfolds)
   for(j in 1:nfolds) {
     jj <- !is.na(smallmat[, j])
     res[[j]] <- smallmat[jj, j]
   }
   return(res)
 }


diag.disc <-function(x, centroids, prior, weight) {
### Computes the class discriminant functions assuming scaled x and centroids
# The returned value is -0.5*diag.disc -\sum_{i=1}^p x_i^{*2}/((s_i+s_0)^2)   
# But the subtracted part is indept. of class. So it is not necessary to calculate it. 
# The function diag.disc returns a n_t by K matrix, where n_t is the number of columns of x.  
  if(! missing(weight)) {
    posid <- (weight > 0)
    if(any(posid)) {
      weight <- sqrt(weight[posid])
      centroids <- centroids[posid,  , drop = FALSE] * weight
      x <- x[posid,  , drop = FALSE] * weight
    }
    else {
      mat <- outer(rep(1, ncol(x)), log(prior), "*")
      dimnames(mat) <- list(NULL, dimnames(centroids)[[2]])
      return(mat)
    }
  }
  dd <- t(x) %*% centroids
  dd0 <- drop(rep(1, nrow(centroids)) %*% (centroids^2))/2 - log(prior)
  names(dd0) <- NULL
  scale(dd, dd0, FALSE)
}


permute.rows <-function(x)
{
        dd <- dim(x)
        n <- dd[1]
        p <- dd[2]
        mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
        matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}

pamr.predict.hard <-  function(fit, newx, threshold, type = c("class", "posterior", "centroid", "nonzero"), 
                          prior = fit$prior,  threshold.scale = fit$threshold.scale) {
  norm.cen <- fit$norm.cen
  if(!is.null(norm.cen)) {
    newx <- abs(t(scale(t(newx), center = norm.cen, scale = FALSE)))
  }
  type <- match.arg(type)
  sd <- fit$sd
  centroid.overall <- fit$centroid.overall
  centroids <- fit$centroids
  se.scale <- fit$se.scale
  delta <- scale((centroids - centroid.overall)/sd, FALSE, threshold.scale * 
                 se.scale)

  if(fit$sign.contrast=="positive"){delta <- delta*(delta>0)}
  if(fit$sign.contrast=="negative"){delta <- delta*(delta<0)}


  delta.shrunk <- scale(hard.shrink(delta, threshold), FALSE,
                        1/(  threshold.scale * se.scale))
  posid <- drop(abs(delta.shrunk) %*% rep(1, length(prior))) > 0
                
  if(!match(type, c("centroid", "nonzero"), FALSE))
    dd <- diag.disc((newx - centroid.overall)/sd, delta.shrunk, 
                    prior, posid)
  switch(type,
         class = softmax(dd),
         posterior = {
           dd <- safe.exp(dd)
           dd/drop(dd %*% rep(1, length(prior)))
         }
         ,
         centroid = centroid.overall + delta.shrunk * sd,
         nonzero = {
           nz <- drop(abs(delta.shrunk) %*% rep(1, ncol(centroids)
                                                )) > 0
           seq(nz)[nz]
         }
         )
}


#############################################################################
my.which.min=function(ddd){ length(ddd)+1-which.min(rev(ddd))}

# perform cross-validation deeply until the threshold parameter is determined completely. 
my.pamr.Scv.hard=function(train.dat){
khan.train <- pamr.train.hard(train.dat)
khan.results<- pamr.cv.hard(khan.train,train.dat ,nfold=10)
result.matrix= cbind(khan.results$threshold, khan.results$size, khan.results$error)

loca=my.which.min(khan.results$error)

nonGenDif=khan.results$size[max(1,loca-1)]-khan.results$size[min(loca+1,nrow(result.matrix))]
if(nonGenDif>1){
q1= seq( khan.results$threshold[max(1,loca-1)], khan.results$threshold[min(loca+1,nrow(result.matrix))], length=min(30, nonGenDif)  )

condition=T
while (condition){
 
   TR=pamr.train.hard(train.dat, threshold=q1)
   CV=pamr.cv.hard(TR,train.dat,nfold=10)
    Smtrx=cbind(CV$threshold,CV$size,CV$error)
  
 result.matrix=rbind(result.matrix, Smtrx)
 where= my.which.min(Smtrx[,3])

condition1=(Smtrx[max(1,where-1), 2] - Smtrx[min(where+1,nrow(Smtrx)), 2] )
q1= seq(Smtrx[max(1,where-1), 1],Smtrx[min(where+1,nrow(Smtrx)), 1], length=min(30, condition1) )
 
    condition=(condition1>1 )& (nrow(Smtrx)>2)
}}

#list(finalresult=result.matrix[my.which.min(result.matrix[,3]  ),  ], process=result.matrix)
result.matrix[my.which.min(result.matrix[,3]  ),  ]

}