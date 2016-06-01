# Multi-QMDR
Detect gene-gene interaction/epistasis by using multivariate phenotypes and MDR

## Notice
All source codes were listed in file "Multi-QMDR.R" for implemeting Multi-QMDR 

## Citation
Yu, Wenbao, M-S. Kwon, and T. Park. "Multivariate Quantitative Multifactor Dimensionality Reduction for Detecting Gene-Gene Interactions." Human heredity 79.3-4 (2015): 168-181.

## Example 
try run_example.R to get a quick start, and an example data was provided too


## Usage
In R:

```
source("Multi-QMDR.R")

multi_qmdr(phes, snp.mat, method = 'FPC', K = 2, nperm = 0, test.type = 'ht2', 
           sele.type = 'cvc', kfolds = 10, covrt = NULL)
```

* inputs: 
  * phes      ---- phenotypes, n times d matrix
  * snp.mat   ---- snp matrix, n times p matrix
   * method   ---- one of 'FPC', 'WPC', 'WSPC'; if you just use single phenotype (qmdr),
                  specify method=i, which will use qmdr and the ith column of phes as phenotype; defalt 'FPC'
   * K        ---- K-way interactions, default 2
   * nperm    ---- permutation times for calculating pvalue for the best model (0 if pvalue if not needed; default)
   * test.type ---- test statistics, could be 'ht2' or 't', corresponding to hotelling t2 test and t test; default 'ht2'
   * sele.type ---- the way to tune the best model, 'cvc' or 'score'; default 'cvc'
   * kfolds    ---- k-fold cross validation; default 10
   * covrt     ---- the covariate matrix; default NULL (no covariates)

* output: a list with elements as follows
    *  best_ksnps ---- the snp ids for the best model
    *  cvc        ---- the cvc number of the best model (out of 10)
    *  score      ---- the test statstics for the best model
    * pvalue      ---- the corresponding empirical pvalue for the best model
      
