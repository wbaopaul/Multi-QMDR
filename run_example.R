
## run Multi-QMDR on example data


source("Multi-QMDR.R")

## load example data
load(file = 'example_data.Rdata')


multi_qmdr(phes, snp.mat[, 1:5], method = 'FPC', K = 2, nperm = 0, 
           test.type = 'ht2', sele.type = 'cvc', kfolds = 10, covrt = covrt)