
library(MASS)
library(pamr)

code_dir = '~/Mohammad/Threshold/both_sides_deep_search/'
#code_dir='C:/Users/everyday/Dropbox/my.notes.on.student.advising/Mohammad_heteros/revision/further_code/sahtout_heter_deep_search_code/'
source(paste0(code_dir, 'pamr.hard.R'))
source(paste0(code_dir, "pamr.soft.R"))
source(paste0(code_dir, "pamr.order.R"))
source(paste0(code_dir, "pamr.soft.hetero.R"))
source(paste0(code_dir, "pamr.hard.hetero.R"))
source(paste0(code_dir, "pamr.order.hetero.R"))


# Class2 first ngenes have high correlations (AR1) with nonzero mean signals
# Class1 all genes are independent with mean 0. 

# ntr1, ntr2: class1 and class2 training sample sizes 
# nte1, nte2: class1 and class2 test sample sizes 
# ngene: number of informative genes
# p: total number of genes
# s: vector of standard deviations of genes in the same class. 
#    Different genes could have different variance but the same gene's variance is same 
#    for all classes. 
generate_sample = function(ntr1, ntr2, nte1, nte2,ngene = 2, p=10, s=rep(3, p), a=100){
  N1 =ntr1+nte1
  N2= ntr2 +nte2 

  
  mu2_nonzero = a/10*seq(ngene)
  rho=0.95
  cors = rho^(0:(ngene-1))   #correlation between successive informative genes

  # define covariance matrix
  Sigma = matrix(sqrt(cors)*s[1:ngene], ncol=1)%*% matrix(sqrt(cors)*s[1:ngene], nrow=1)

  #Generate ngene random samples that are correlated with covariance matrix Sigma
  sim_data1_informative <- mvrnorm(n = N1, mu = rep(0, ngene), Sigma = Sigma)
  sim_data2_informative <- mvrnorm(n = N2, mu = mu2_nonzero, Sigma = Sigma)
  
  x1_other = scale(matrix(rnorm((p-ngene)*N1, mean=0, sd=1) , nrow=N1, ncol=p-ngene), center=F, scale=1/s[(ngene+1):p] )
  x2_other = scale(matrix(rnorm(N2*(p-ngene), 0, 1), nrow=N2, ncol=p-ngene), center=F, scale=1/s[(ngene+1):p])
  
  x1 = cbind(sim_data1_informative, x1_other)
  x2 = cbind(sim_data2_informative, x2_other)
  

  Xtr = rbind(x1[1:ntr1,], x2[1:ntr2,])
  ytr = c(rep('A', ntr1), rep('B', ntr2))

  Xte = rbind(x1[-(1:ntr1),], x2[-(1:ntr2),])
  yte = c(rep('A', nte1), rep('B', nte2))
  rm(x1, x2, Sigma, sim_data1_informative,sim_data2_informative, x1_other, x2_other)
  gc()
return(list(Xtr=Xtr, ytr=ytr, Xte=Xte, yte=yte) )
}



#Core metrics functions
# =============================================================================

# ----  Per-class recall (sensitivity) from a confusion matrix ------------
# cm: square matrix, rows = true class, cols = predicted class
per_class_recall <- function(cm) {
  diag(cm) / pmax(rowSums(cm), 1L)
}

# ----  Balanced accuracy = macro-averaged recall -------------------------
balanced_acc <- function(cm) mean(per_class_recall(cm))

#Accumulate a confusion matrix from one run's predictions ----------
make_cm <- function(true_labels, pred_labels, classes) {
  cm <- matrix(0L, nrow = length(classes), ncol = length(classes),
               dimnames = list(true = classes, pred = classes))
  for (i in seq_along(true_labels)) {
    t <- as.character(true_labels[i]); p <- as.character(pred_labels[i])
    if (t %in% classes && p %in% classes)
      cm[t, p] <- cm[t, p] + 1L
  }
  cm
}




## a good setting to keep when x1 directly us s in rnorm rather than scale
#   ntr1=20, ntr2=30, nte1=50, nte2=50, ngene = 20, p=10000, s=1:p, a=100
# But that is not homo 
  
# for one dataset, run N_RUNS times by shuffling the order of training samples
# x: training x of dim pxNtr
# y: training y vector
# xte1: test x of dim px Nte
# yte1: test y vector

run_one_dataset = function( N_RUNS, x, y, xte1, yte1, whatdata, OUT_DIR){
  ALL_METHODS <- c(
    "STh3", "HTh3", "OTh3",           # proposed heteroscedastic NSC
    "STh",  "HTh",  "OTh")
  
  cl  <- sort(unique(y))
  K   <- length(cl);  n = nrow(x)
  
  # Storage: list[method] -> list of K x K CMs (one per run)
  cm_list <- setNames(
    lapply(ALL_METHODS, function(m) vector("list", N_RUNS)),
    ALL_METHODS
  )
  
  ave_num_genes = list()
  error = list()
  
  for (method in c("Hard", "Soft", "Order")){
    method_sub = c(NULL, NULL)
    n_genes_homo = numeric(N_RUNS)
    n_genes_heter = numeric(N_RUNS)
    errors_homo = numeric(N_RUNS)
    errors_heter = numeric(N_RUNS)
    
  
    
    for(itr in 1: N_RUNS){
      set.seed(1000+itr*10)
      
      randomorder=sample(length(y))
      x=x[,randomorder];y=y[randomorder]
      train.dat=list(x=x,y=y)
      
      if (method=="Soft"){
        method_sub = c('STh', 'STh3')
        pick=my.pamr.Scv.old(train.dat)
        thr=pick[1]
        TR=pamr.train.old(train.dat, threshold=thr)
        pred=pamr.predict(TR, xte1 , threshold=thr)
        
        pick.heter=my.pamr.Scv.heter(train.dat)
        thr.heter=pick.heter[1]
        TR.heter=pamr.train.heter(train.dat, threshold=thr.heter)
        pred.heter=pamr.predict.heter(TR.heter, xte1 , threshold=thr.heter)
        
        
      }
      
      if (method=="Order"){
        method_sub = c('OTh', 'OTh3')
        pick=my.pamr.Scv.order(train.dat)
        thr=pick[1]
        TR=pamr.train.order(train.dat, threshold=thr)
        pred=pamr.predict.order(TR, xte1 , threshold=thr)
        
        pick.heter=my.pamr.Scv.O.heter(train.dat)
        thr.heter=pick.heter[1]
        TR.heter=pamr.train.O.heter(train.dat, threshold=thr.heter)
        pred.heter=pamr.predict.O.heter(TR.heter, xte1 , threshold=thr.heter)
      }
      
      if (method=="Hard"){
        method_sub = c('HTh', 'HTh3')
        pick=my.pamr.Scv.hard(train.dat)
        thr=pick[1]
        TR=pamr.train.hard(train.dat, threshold=thr)
        pred=pamr.predict.hard(TR, xte1 , threshold=thr)
        
        pick.heter=my.pamr.Scv.H.heter(train.dat)
        thr.heter=pick.heter[1]
        TR.heter=pamr.train.H.heter(train.dat, threshold=thr.heter)
        pred.heter=pamr.predict.H.heter(TR.heter, xte1 , threshold=thr.heter)
      }
      
      cm_list[[method_sub[1] ]][[itr]] <- make_cm(yte1, as.character(pred), cl)
      cm_list[[method_sub[2] ]][[itr]] <- make_cm(yte1, as.character(pred.heter), cl)
      
      n_genes_homo[itr] = TR$nonzero
      n_genes_heter[itr] = TR.heter$nonzero
      errors_homo[itr] = mean(as.character(pred) != yte1) 
      errors_heter[itr] = mean(as.character(pred.heter) != yte1) 
      
      #collection= c(thr,TR$nonzero,TR$errors, sum(as.character(pred) != yte1) )
      #col=rbind(col,collection)
      
      #collection.heter= c(thr.heter,TR.heter$nonzero,TR.heter$errors, sum(as.character(pred.heter) != yte1) )
      #col.heter=rbind(col.heter,collection.heter)
      #cat(paste("LL=", LL,sep=""), collection,  "\n", file=paste(method,  whatdata,"txt", sep="."), append=T  )
      
    }
    a = switch(method, 'Hard'='HTh', 'Soft'='STh', 'Order'='OTh')
    b = switch(method, 'Hard'='HTh3', 'Soft'='STh3', 'Order'='OTh3')
    ave_num_genes[[a]] = mean(n_genes_homo)
    ave_num_genes[[b]] = mean(n_genes_heter)
    error[[a]] = mean(errors_homo)
    error[[b]] = mean(errors_heter)
    }
  rho2 = compute_rho2(x, y)
  res =  list(cm_list = cm_list, classes = cl, n = n, K = K, ave_num_genes=ave_num_genes, error = error, rho2= rho2)
  rds_path <- file.path(OUT_DIR, paste0("simres_Ntr",dim(x)[2],'_Nte',dim(xte1)[2],'_a', whatdata, ".rds"))
  saveRDS(res, rds_path)
  res
  }


# =============================================================================
# compute_rho2.R
#
# Computes the empirical relative heterogeneity rho_k^2 for each class k,
# as defined in the manuscript (equation after eq. 10):
#
#   rho_k^2 = (sigma_k^2 / sigma_bar^2 - 1)^2
#
# where sigma_k^2  is the true (unknown) class-k variance, estimated here
#              by s_ik^2 from the TRAINING data (unbiased, ddof = 1)
#   and sigma_bar^2 = (1/n) * sum_k n_k * sigma_k^2  (sample-size-weighted
#              average variance across classes), also estimated from training.
#
# Because the training sample sizes n_k may be very small (e.g. n_k = 3 for
# Basal and ERBB2 in the Breast dataset), the within-class variance estimates
# s_ik^2 from the training data are noisy.  An optional argument allows the
# CLASS VARIANCES to be estimated from a SEPARATE (test/held-out) dataset
# instead, while the CLASS WEIGHTS n_k still come from the training partition.
# This is the approach used in the manuscript for the Breast dataset.
#
# DATA FORMAT
# -----------
# X_train, X_test : numeric matrices with GENES as rows and SAMPLES as columns.
# y_train, y_test : character (or factor) vectors of class labels, one per sample.
#                   Length must equal ncol(X_train) / ncol(X_test) respectively.
#
# RETURNS
# -------
# A data frame with one row per class containing:
#   class        : class label
#   n_k          : training sample size for that class
#   sigma2_hat   : estimated class variance (median across genes)
#   sigma_bar2   : weighted-average variance (scalar, same for all rows)
#   rho2_hat     : median empirical rho_k^2 across genes  [MAIN OUTPUT]
#   threshold    : simple crossover threshold 2/(n_k - 1)  (w_k = 1 bound)
#   crossover_met: logical, TRUE when rho2_hat > threshold
# =============================================================================


#' Compute empirical rho_k^2 for each class
#'
#' @param X_train  Gene x sample matrix of training expression values.
#' @param y_train  Class labels for training samples (length = ncol(X_train)).
#' @param X_test   Gene x sample matrix used to ESTIMATE class variances.
#'                 Defaults to X_train.  Supply a separate test/held-out set
#'                 when training class sizes are too small for stable estimates.
#' @param y_test   Class labels for X_test (length = ncol(X_test)).
#'                 Defaults to y_train.
#' @param summary_fn Function applied across genes to summarise rho_k^2.
#'                   Default is median.  Use mean for a different summary.
#'
#' @return A data frame (one row per class) with columns described above.
compute_rho2 <- function(X_train,
                         y_train,
                         X_test      = X_train,
                         y_test      = y_train,
                         summary_fn  = median) {
  
  # ---- Input checks ---------------------------------------------------------
  stopifnot(is.matrix(X_train), is.matrix(X_test))
  stopifnot(length(y_train) == ncol(X_train))
  stopifnot(length(y_test)  == ncol(X_test))
  stopifnot(nrow(X_train)   == nrow(X_test))
  
  y_train <- as.character(y_train)
  y_test  <- as.character(y_test)
  classes <- sort(unique(y_train))
  p       <- nrow(X_train)
  
  # ---- Step 1: per-class training sample sizes n_k --------------------------
  n_k <- setNames(
    sapply(classes, function(k) sum(y_train == k)),
    classes
  )
  n   <- sum(n_k)   # total training n
  
  # ---- Step 2: per-class variance estimates (from X_test) -------------------
  # s2[k, ] is a length-p vector of per-gene within-class sample variances
  # for class k, estimated from the test (or training) data.
  s2 <- sapply(classes, function(k) {
    Xk <- X_test[, y_test == k, drop = FALSE]
    if (ncol(Xk) < 2L) {
      warning(sprintf(
        "Class '%s' has only %d sample(s) in the variance-estimation dataset; ",
        k, ncol(Xk),
        "variance estimates will be NA."
      ))
      return(rep(NA_real_, p))
    }
    apply(Xk, 1L, var)          # var() uses ddof = 1 by default in R
  })
  # s2 is now p x K (genes x classes)
  
  # ---- Step 3: sigma_bar^2 = (1/n) * sum_k n_k * sigma_k^2 (per gene) ------
  #   sigma_bar2_vec[i] = (1/n) * sum_k n_k * s2[i, k]
  sigma_bar2_vec <- as.numeric(s2 %*% n_k) / n   # length p
  
  # ---- Step 4: rho_k^2 = (sigma_k^2 / sigma_bar^2 - 1)^2  (per gene) ------
  rho2_mat <- sweep(s2, 1L, sigma_bar2_vec, "/") - 1L
  rho2_mat <- rho2_mat^2    # p x K matrix
  
  # ---- Step 5: summarise across genes for each class ------------------------
  rho2_summary    <- apply(rho2_mat, 2L, summary_fn, na.rm = TRUE)
  sigma2_summary  <- apply(s2,       2L, summary_fn, na.rm = TRUE)
  sigma_bar2_val  <- summary_fn(sigma_bar2_vec, na.rm = TRUE)
  
  # ---- Step 6: crossover threshold 2 / (n_k - 1)  (simple bound, w_k = 1) -
  threshold <- ifelse(n_k > 1L, 2 / (n_k - 1L), Inf)
  
  # ---- Assemble output ------------------------------------------------------
  out <- data.frame(
    class         = classes,
    n_k           = n_k[classes],
    sigma2_hat    = round(sigma2_summary[classes],  4L),
    sigma_bar2    = round(sigma_bar2_val,            4L),
    rho2_hat      = round(rho2_summary[classes],     4L),
    threshold     = round(threshold[classes],         4L),
    crossover_met = rho2_summary[classes] > threshold[classes],
    row.names     = NULL,
    stringsAsFactors = FALSE
  )
  
  out
}


N_RUNS = 100
OUT_DIR     <- "~/Mohammad/Threshold/sahtout_homo_with_corr_simu/result"    
#OUT_DIR = getwd() 
#if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
#if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

# This setting does not satisfies 
#  - Reviewer1 comemnt 4, varying sample sizes while fixing the 
#      effect size and variance structure
#  - Reviewer1 comment 7, homoscedastic settings where the equal-
#      variance assumption holds




m=20
p=10000; s=rep(rep(1:m), p/m)

As = 9:23
res=list()
for (a in As){
     set.seed(iter*8900)
     dat2  = generate_sample(ntr1=15, ntr2=20, nte1=50, nte2=50, ngene = 20, p=p, s=s, a=a)
     res2=run_one_dataset( N_RUNS, x=t(dat2$Xtr), y=dat2$ytr, xte1=t(dat2$Xte), yte1=dat2$yte, 
            whatdata=paste0(a, '_iter', iter), OUT_DIR) 
     res[[a-8]] =res2
}




result = rbind(data.frame(res[[1]]$error), 
    data.frame(res[[2]]$error), data.frame(res[[3]]$error), 
    data.frame(res[[4]]$error), data.frame(res[[5]]$error), 
    data.frame(res[[6]]$error), data.frame(res[[7]]$error), 
    data.frame(res[[8]]$error), data.frame(res[[9]]$error), 
    data.frame(res[[10]]$error), data.frame(res[[11]]$error), 
    data.frame(res[[12]]$error),data.frame(res[[13]]$error),
    data.frame(res[[14]]$error),data.frame(res[[15]]$error))   
result$a =As

write.csv(result, file =file.path(OUT_DIR, paste0('iter', iter, '_sahtout_sim_homo_correlated.csv')) )


# =============================================================================
# SECTION 3: Summary tables for balanced accuracy and average recall sensitivity
# =============================================================================


# ----  Per-class recall (sensitivity) from a confusion matrix ------------
# cm: square matrix, rows = true class, cols = predicted class
per_class_recall <- function(cm) {
  diag(cm) / pmax(rowSums(cm), 1L)
}

# ---- Balanced accuracy = macro-averaged recall -------------------------
balanced_acc <- function(cm) mean(per_class_recall(cm))


# ---- 3.1  Mean CM over runs (average count per cell) ------------------------
mean_cm <- function(cm_runs) {
  valid <- Filter(Negate(is.null), cm_runs)
  if (!length(valid)) return(NULL)
  Reduce("+", valid) / length(valid)
}

# ---- 3.2  Build per-class recall table: one row per (dataset x class),
#           one column per method -------------------------------------------
build_recall_table <- function(all_res3ults) {
  rows <- list()
  for (res3 in all_res3ults) {
    #ds      <- res3$dsname
    classes <- res3$classes
    for (k in classes) {
      row <- data.frame(#Dataset = ds, 
        Class = k, stringsAsFactors = FALSE)
      for (meth in names(res3$cm_list)) {
        mcm <- mean_cm(res3$cm_list[[meth]])
        if (is.null(mcm)) {
          row[[paste0("recall_", meth)]] <- NA_real_
        } else {
          row[[paste0("recall_", meth)]] <- round(per_class_recall(mcm)[k], 3)
        }
      }
      rows <- c(rows, list(row))
    }
  }
  do.call(rbind, rows)
}

# ---- 3.3  Balanced accuracy table: one row per dataset, one column per method
build_balanced_acc_table <- function(all_res3ults) {
  do.call(rbind, lapply(all_res3ults, function(res3) {
    row <- data.frame( 
                      K = res3$K, n = res3$n,
                      stringsAsFactors = FALSE)
    for (meth in names(res3$cm_list)) {
      mcm <- mean_cm(res3$cm_list[[meth]])
      row[[paste0("BalAcc_", meth)]] <-
        if (is.null(mcm)) NA_real_ else round(balanced_acc(mcm) * 100, 2)
    }
    row
  }))
}



cat("\nBuilding summary tables ...\n")
recall_tab  <- build_recall_table(res)
recall_tab$a = kronecker( As, c(1,1))
recall_tab$iter =iter
bal_acc_tab <- build_balanced_acc_table(res)
bal_acc_tab$a = As
bal_acc_tab$iter = iter

write.csv(recall_tab,  file.path(OUT_DIR, paste0("homo_sahtout_iter", iter,"_corr_class_sensitivity.csv")),
          row.names = FALSE)
write.csv(bal_acc_tab, file.path(OUT_DIR, paste0("homo_sahtout_iter", iter,"corr_balanced_accuracy.csv")),
          row.names = FALSE)

cat("\n=== BALANCED ACCURACY (%) BY METHOD AND DATASET ===\n")
print(bal_acc_tab, row.names = FALSE)



