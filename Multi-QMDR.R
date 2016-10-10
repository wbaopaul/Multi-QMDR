##=======================================================================##
## **** this program is written by Wenbao Yu to implement Multi-QMDR 
##=======================================================================##


## list all functions 

## calculate training score for a given ids for train samples and
## snp.dat, n by k, -- the k-th column correspond to the k-th snps
## snp.dat could have missing values and test score was calculated for each model
qmdr <- function(train.ids, test.ids, snp.dat, S, phes, test.type = 'ht2', method){
  snp.dat = as.matrix(snp.dat)
  
  ## remove missing values
  fids = which(complete.cases(snp.dat))
  snp.dat = as.matrix(snp.dat[fids, ])
  
  test.ids = intersect(test.ids, fids)
  train.ids = intersect(train.ids, fids)
  
  #sbar = mean(S[train.ids])
  sbar = mean(S)
  k = ncol(snp.dat) 
  
  ## split data into cells
  tlist = vector('list', k)
  for(i in 1:k) tlist[[i]] = snp.dat[, i]
  cells = split(data.frame(cbind(fids, snp.dat)), tlist)
  
  ## delete NULL cells
  obs.cell = sapply(cells, function(x) nrow(x))
  cell.null = which(obs.cell == 0)
  if(length(cell.null) > 0 ) cells = cells[-cell.null]
  
  ## get trainid in each cell
  cells.trainid = lapply(cells, function(x) return(intersect(x[, 1], train.ids)))
  
  cells.num = length(cells)
  #cells.label = rep(0, cells.num)
  high.all  = NULL
  for(i in 1:cells.num){
    temp.ids = cells.trainid[[i]]
    if(length(temp.ids) == 0) next
    if (mean(S[temp.ids]) >= sbar){
      high.all = c(high.all, cells[[i]][, 1])
    }
  }
  
  if(test.type == 't'){
    train.stat = cal_tstat(train.ids, high.all, S)
    test.stat = cal_tstat(test.ids, high.all, S)
  }else{
    train.stat = cal_ht2(train.ids, high.all, phes)
    test.stat = cal_ht2(test.ids, high.all, phes)
  }
  return(list('train.stat' = train.stat, 'test.stat' = test.stat ))
}


## calculating t score 
cal_tstat <- function(ids, high.all, S){  
  high.ids = intersect(ids, high.all)
  low.ids = setdiff(ids, high.ids)
  
  s1 = S[high.ids]
  
  s2 = S[low.ids]
  
  if(length(high.ids) == 0 || length(low.ids) == 0) return(0)
  stat = t.test(s1, y = s2, var.equal = TRUE)$statistic
  
  return(abs(stat))
}


## calculating ht2 score 
cal_ht2 <- function(ids, high.all, phes){  

  high.ids = intersect(ids, high.all)
  low.ids = setdiff(ids, high.ids)
  d = ncol(phes)
  
  s1 = phes[high.ids, ]
  
  s2 = phes[low.ids, ]
  
  if(length(high.ids) == 0 || length(low.ids) == 0) return(0)
  
  s1 = as.matrix(s1)
  s2 = as.matrix(s2)
  
  if(ncol(s1) == 1) s1 = t(s1)
  if(ncol(s2) == 1) s2 = t(s2)
  
           
  stat = dire_ht2(s1, s2, phes)$fstat                   ## another version
  ## stat is scaled that it follows a F distribution with degree d, n-1-d under the null
  return(stat)
}


## calculate HT2 directly
dire_ht2 <- function(X, Y, phes){
  
  # number of observations for two group:
  l1 <- nrow(X)
  l2 <- nrow(Y)
  d = ncol(X)
  
  # Sample mean vectors for the each group:
  m1 <- apply(X, 2, mean)
  m2 <- apply(Y, 2, mean)
  
  # "pooled" sample covariance matrix:
  poolS <- ((l1 - 1) * cov(X) + 
              (l2 - 1) * cov(Y))/(l1 + l2 - 2)
  
  if(any(is.na(poolS)) || abs(det(poolS)) < 0.00001) poolS <- cov(phes)
  
  # Hotelling T^2, the F-statistic, and the P-value:
  T2 <- ((l1 * l2)/(l1 + l2))* (t(m1 - m2) %*% solve(poolS) %*% (m1 - m2) ) 
  
  Fstat <- ((l1 + l2- d - 1) * T2)/((l1 + l2 - 2) * d)
  # pvalue <- pf(Fstat, d, l1 + l2 - d - 1, lower.tail = FALSE)
  
  return(list("stat" = round(T2, 4), 
              "fstat" = round(Fstat, 4) ))
  
}

## choose the best k-way interaction by cross validation
## return the snp combination and its test score
## snp.combs includes a possible k-way snps combinations in each column
## using cv consistency or testing score to decide best model
tune_kmodel <- function(folds = 10, snp.all, S, phes, 
                        test.type = 'ht2', sele.type = 'cvc', snp.combs, method){
  ns = ncol(snp.combs)
  n = length(S)
  
  ## split the whole data into folds
  cvlen = floor(n/folds)
  cc = 1:n
  test.stats = train.stats = rep(0, ns)
  temptest.stats = matrix(0, ns, folds)
  best.comb = rep(0, folds)
  
  ## select best model(i.e. snp combination)
  for(i in 1:folds){
    testid = ((i - 1) * cvlen + 1) : (i * cvlen)
    trainid = cc[-testid]
    
    for(j in 1:ns){
      temp.result = qmdr(trainid, testid, snp.all[, snp.combs[, j]], S, phes, test.type, method)
      train.stats[j] = temp.result$train.stat
      temptest.stats[j, i] = temp.result$test.stat
    }
    # which snp pair has best training stat for each trainind set
    best.comb[i] = j0 = which.max(train.stats)      
  }
  
  test.stats = rowMeans(temptest.stats, na.rm = TRUE)  ## average testing stats for all snp pairs
  
  if(sele.type == 'cvc'){
    ta = table(best.comb)
    cvc = max(ta)          ## the largest cvc
    best.pair = as.numeric(names(which(ta == cvc)))[1]  ## the pair gets largest cvc
  }
  if(sele.type == 'score'){
    best.pair = which.max(test.stats)  ## the pair gives larges test score
    cvc = length(which(best.comb == best.pair))
  }
  
  sele.score = test.stats[best.pair]
  
  ## sele.score -- corresponding to the test score of the final selected model
  ## test.stats -- record test.scores for all possible k-way model (snp interactions)
  
  ## Save the test score of the best model in each cv
  scores.cv = temptest.stats[best.pair, ]
  
  return(list('cvc' = cvc, 'score' = sele.score, 'scores.cv' = scores.cv,
               'best.pair' = best.pair, 'test.stats' = test.stats))
  
}


## calculate empirical distribution under the null using all pairs
## this empirical distribution is used to calculate pvalues for all snp
emp_null_perm <- function(K, S, snp.all, phes, B = 1000, folds = 10, 
                       test.type = 'ht2', method){
  set.seed(1234)
  
  n = nrow(snp.all)
  p = ncol(snp.all)
  
  
  snp.combs <- combn(p, K)  ## all possible combinatory pairs
  ns  = ncol(snp.combs)
  
  test.stats = matrix(0, ns, folds)
  cvlen = floor(n / folds)
  cc = 1:n 
  
  ## defube high/low by train and calculate test score and average them
  ## only permute phenotype S and phes
  run = 0
  stats = NULL
  repeat{
    run = run + 1
    perm.id = sample(1:n, n)
    perm.phes = phes[perm.id, ]
    perm.S = S[perm.id]
    for(i in 1:(ns)){
      for(j in 1:folds){
        testid = ((j - 1) * cvlen + 1) : (j * cvlen)
        trainid = cc[-testid]
        temp.result = qmdr(trainid, testid, snp.all, perm.S, perm.phes, test.type, method)
        test.stats[i, j] = temp.result$test.stat
      }  
    }
    stats = rbind(stats, test.stats)
    if(run == B) break
  }
  
  emp_stats_null = rowMeans(stats, na.rm = TRUE)
  return(emp_stats_null)
}


## function to do k-way interaction using Multi-QMDR 
## output including the pvalue for the best model
multi_qmdr <- function(phes, snp.mat, method = 'FPC', K = 2, nperm = 1000, test.type = 't', 
                       sele.type = 'cvc', kfolds = 10, covrt = NULL){
 
 # adjust covariant's effect for each phenotype
  if(!is.null(covrt)) {
    fun <- function(y){
      resid = lm(y ~ covrt)$residuals
      return(resid)
    }
    phes = apply(phes, 2, fun)
  }

  SS = summary_score(phes, method) 
  
  set.seed(1)
  n = nrow(phes) 
  p = ncol(snp.mat)
  snp.combs <- combn(p, K)  ## all possible combinatory pairs
  ns  = ncol(snp.combs)
  test.stats = rep(0L, ns)
  
  aa = sample(1:n, n)  ## shuffle samples
 
  result = tune_kmodel(10, snp.mat[aa, ], SS[aa], phes[aa, ], test.type, 
                       sele.type, snp.combs, method)
  
  model.cons = result$cvc
  model.sele = result$best.pair
  model.score = result$score
  test.stats = result$test.stats
  scores.cv = result$scores.cv
  
  best.ksnps = snp.combs[, model.sele]
  
  perm.pv = NULL
  if(nperm > 0){
    nperm = max(1, ceiling(nperm/ns))
    emp_stats_null = emp_null_perm(K, SS, snp.all = snp.mat, 
                                   phes, B = nperm, folds = kfolds, test.type, method)
    perm.pv = mean(ifelse(emp_stats_null > model.score, 1, 0))
  }
  
  
  #perm.pv = permute_pv(tscore = model.score, loci = best.ksnps, SS, snp.all = snp.mat, 
  #                     phes, B = nperm, folds = 10, test.type, method)
  
  
  
  
  res = list('best_ksnps' = best.ksnps, 'cvc' = model.cons, 'score' = model.score, 'pv' = perm.pv)
  return(res)
  
}



## return weighted summation of pc scores
wpc <- function(pys){
  pys = scale(pys, center = TRUE, scale = FALSE)
  d = ncol(pys) 
  covm = cov(pys)
  eigens = eigen(covm)
  lbd = eigens$values
  P = eigens$vectors
  #wP = P * matrix(rep(sqrt(lbd), d), d, byrow=T)
  wP = pys %*% P
  wP = t(apply(wP, 1, function(x) x * sqrt(lbd)))
  AS = rowSums(wP)  ## weighted summation of pcs
  return(AS)
}

wspc <- function(pys){
  pys = scale(pys, center = TRUE, scale = FALSE)
  d = ncol(pys) 
  covm = cov(pys)
  eigens = eigen(covm)
  lbd = eigens$values
  P = eigens$vectors
  #wP = P * matrix(rep(sqrt(lbd), d), d, byrow=T)
  #AS = rowSums((pys %*% wP)^2)  ## weighted square summation of pcs
  wP = (pys %*% P)^2
  wP = t(apply(wP, 1, function(x) x * lbd))
  AS = rowSums(wP)  ## weighted summation of pcs
  return(AS)
}

## return fist PC score
fpc <- function(pys){
  pys = scale(pys, center = TRUE, scale = FALSE)
  covm = cov(pys)
  eigens = eigen(covm)
  P = eigens$vectors
  AS = pys %*% P[, 1]  ## first PC scores
  return(AS)
}


## aggregating phenotypes
summary_score <- function(phes, method){
  d = ncol(phes)
  if(method %in% 1:d){
      SS = phes[, method]
  }
  
  if(method == 'WPC'){
    SS = wpc(phes)
  }
  if(method == 'WSPC'){
    SS = wspc(phes)
  }
  if(method == 'FPC'){
    SS = fpc(phes) 
  }
  return(SS)
}




