library(Rcpp)
sourceCpp('co_two.cpp')
sourceCpp('dist_co.cpp')
library(fastDummies)
library(sparcl)
library(cluster)
library(VarSelLCM)
library(flexclust)
library(pgnorm)




################ Simulation Data Functions ###################
### There are K=2,3,4 three cases.                         ###
### Generate normal and p-generalized normal-polynomial distribution.
### n1 is the sample size of the first cluster.            ###
### n2 is the sample size of the other cluster.            ###
### p is the number of the variables.                      ###



# K=2
# p-generalized normal-polynomial distribution.
generatesimd <- function(n1, n2, p){
  #conti_mu <- c(rep(1,5),rep(0,195))
  #contivar <- as.data.frame(rbind(mvrnorm(n = 50, mu = conti_mu, Sigma = diag(nrow = 200)), mvrnorm(n = 50, mu = rep(0,200), Sigma = diag(nrow = 200))))
  conti_signal <- NULL
  for (i in 1:5) {
    conti_signal <- cbind(conti_signal,
                          c(rpgnorm(n = n1,p = 0.7785,mean = 1.5,sigma = 1),
                            rpgnorm(n = n2,p = 0.7785,mean = 0,sigma = 1)))
  }
  for (i in 1:5) {
    conti_signal <- cbind(conti_signal,
                          c(rpgnorm(n = n1,p = 0.7785,mean = 0,sigma = 1),
                            rpgnorm(n = n2,p = 0.7785,mean = 1.5,sigma = 1)))
  }
  conti_noise <- NULL
  for (i in 1:(p-10)) {
    conti_noise <- cbind(conti_noise,
                         c(rpgnorm(n = (n1+n2),p = 0.7785,mean = 0,sigma = 1)))
  }
  contivar <- data.frame(cbind(conti_signal, conti_noise))
  
  cate_v1 <- c(sample(1:2,size = n1,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)))
  cate_v2 <- c(sample(1:2,size = n1,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)))
  cate_v3 <- c(sample(1:3,size = n1,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3))
  cate_v4 <- c(sample(1:3,size = n1,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3))
  cate_v5 <- c(sample(1:4,size = n1,replace = T,prob = c(0.4,0.4,0.1,0.1)),
               sample(1:4,size = n2,replace = T,prob = c(0.25,0.25,0.25,0.25)))
  
  cate_v6 <- c(sample(1:2,size = n1,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)))
  cate_v7 <- c(sample(1:2,size = n1,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)))
  cate_v8 <- c(sample(1:3,size = n1,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3))
  cate_v9 <- c(sample(1:3,size = n1,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3))
  cate_v10 <- c(sample(1:4,size = n1,replace = T,prob = c(0.4,0.4,0.1,0.1)),
                sample(1:4,size = n2,replace = T,prob = c(0.25,0.25,0.25,0.25)))
  
  
  catevar <- data.frame(cate_v1,cate_v2,cate_v3,cate_v4,cate_v5,
                        cate_v6,cate_v7,cate_v8,cate_v9,cate_v10)
  
  cate_sig <- c(rep(1,length(unique(cate_v1))),
                rep(2,length(unique(cate_v2))),
                rep(3,length(unique(cate_v3))),
                rep(4,length(unique(cate_v4))),
                rep(5,length(unique(cate_v5))),
                rep(6,length(unique(cate_v6))),
                rep(7,length(unique(cate_v7))),
                rep(8,length(unique(cate_v8))),
                rep(9,length(unique(cate_v9))),
                rep(10,length(unique(cate_v10))))
  
  
  catevardum <- dummy_cols(catevar,select_columns = c('cate_v1','cate_v2','cate_v3','cate_v4','cate_v5',
                                                      'cate_v6','cate_v7','cate_v8','cate_v9','cate_v10'),
                           remove_first_dummy = F, remove_selected_columns = T)
  
  
  catenoi <- NULL
  cate_noi <- NULL
  for (i in 1:(p-10)) {
    temp_level <- sample(2:4,1)
    catenoi <- cbind(catenoi, factor(sample(1:(temp_level),(n1+n2),replace = T)))
    cate_noi <- c(cate_noi, rep(i, temp_level))
  }
  catenoi <- as.data.frame(catenoi)
  catenoidum <- dummy_cols(catenoi,select_columns = colnames(catenoi),
                           remove_first_dummy = F, remove_selected_columns = T)
  
  simd <- cbind(contivar,catevardum,catenoidum)
  simgroup <- c(1:p,cate_sig+p,cate_noi+p+10)
  
  conti_df <- data.frame(cbind(conti_signal,conti_noise,catenoi,
                               cate_v1,cate_v2,cate_v3,cate_v4,cate_v5,
                               cate_v6,cate_v7,cate_v8,cate_v9,cate_v10))
  conti_df[,(p+1):(2*p)] <- lapply(conti_df[,(p+1):(2*p)], factor)
  
  return(list(simd=simd,
              simgroup=simgroup,
              conti_df=conti_df))
}


# K=2
# normal distribution
generatesimd <- function(n1, n2, p){
  conti_mu1 <- c(rep(1.5,5),rep(0,(p-5)))
  conti_mu2 <- c(rep(0,5),rep(1.5,5),rep(0,(p-10)))
  Sigma1 <- genPositiveDefMat(
    dim = p,
    covMethod = 'onion', 
    rangeVar = c(1.5,2.5))
  Sigma2 <- genPositiveDefMat(
    dim = p,
    covMethod = 'onion', 
    rangeVar = c(1.5,2.5))
  contivar <- as.data.frame(rbind(mvrnorm(n = n1, mu = conti_mu1, Sigma = Sigma1$Sigma), 
                                  mvrnorm(n = n2, mu = conti_mu2, Sigma = Sigma2$Sigma)))
  
  
  
  cate_v1 <- c(sample(1:2,size = n1,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)))
  cate_v2 <- c(sample(1:2,size = n1,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)))
  cate_v3 <- c(sample(1:3,size = n1,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3))
  cate_v4 <- c(sample(1:3,size = n1,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3))
  cate_v5 <- c(sample(1:4,size = n1,replace = T,prob = c(0.4,0.4,0.1,0.1)),
               sample(1:4,size = n2,replace = T,prob = c(0.25,0.25,0.25,0.25)))
  
  cate_v6 <- c(sample(1:2,size = n1,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)))
  cate_v7 <- c(sample(1:2,size = n1,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)))
  cate_v8 <- c(sample(1:3,size = n1,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3))
  cate_v9 <- c(sample(1:3,size = n1,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3))
  cate_v10 <- c(sample(1:4,size = n1,replace = T,prob = c(0.4,0.4,0.1,0.1)),
                sample(1:4,size = n2,replace = T,prob = c(0.25,0.25,0.25,0.25)))
  
  
  catevar <- data.frame(cate_v1,cate_v2,cate_v3,cate_v4,cate_v5,
                        cate_v6,cate_v7,cate_v8,cate_v9,cate_v10)
  
  cate_sig <- c(rep(1,length(unique(cate_v1))),
                rep(2,length(unique(cate_v2))),
                rep(3,length(unique(cate_v3))),
                rep(4,length(unique(cate_v4))),
                rep(5,length(unique(cate_v5))),
                rep(6,length(unique(cate_v6))),
                rep(7,length(unique(cate_v7))),
                rep(8,length(unique(cate_v8))),
                rep(9,length(unique(cate_v9))),
                rep(10,length(unique(cate_v10))))
  
  
  catevardum <- dummy_cols(catevar,select_columns = c('cate_v1','cate_v2','cate_v3','cate_v4','cate_v5',
                                                      'cate_v6','cate_v7','cate_v8','cate_v9','cate_v10'),
                           remove_first_dummy = F, remove_selected_columns = T)
  
  
  catenoi <- NULL
  cate_noi <- NULL
  for (i in 1:(p-10)) {
    temp_level <- sample(2:4,1)
    catenoi <- cbind(catenoi, factor(sample(1:(temp_level),n1+n2,replace = T)))
    cate_noi <- c(cate_noi, rep(i, temp_level))
  }
  catenoi <- as.data.frame(catenoi)
  catenoidum <- dummy_cols(catenoi,select_columns = colnames(catenoi),
                           remove_first_dummy = F, remove_selected_columns = T)
  
  simd <- cbind(contivar,catevardum,catenoidum)
  simgroup <- c(1:p,cate_sig+p,cate_noi+p+10)
  
  conti_df <- data.frame(cbind(contivar,catenoi,
                               cate_v1,cate_v2,cate_v3,cate_v4,cate_v5,
                               cate_v6,cate_v7,cate_v8,cate_v9,cate_v10))
  conti_df[,(p+1):(2*p)] <- lapply(conti_df[,(p+1):(2*p)], factor)
  
  return(list(simd=simd,
              simgroup=simgroup,
              conti_df=conti_df))
}


# K=3
# p-generalized normal-polynomial distribution.
generatesimd <- function(n1, n2, p){
  #conti_mu <- c(rep(1,5),rep(0,195))
  #contivar <- as.data.frame(rbind(mvrnorm(n = 50, mu = conti_mu, Sigma = diag(nrow = 200)), mvrnorm(n = 50, mu = rep(0,200), Sigma = diag(nrow = 200))))
  conti_signal <- NULL
  for (i in 1:5) {
    conti_signal <- cbind(conti_signal,
                          c(rpgnorm(n = n1,p = 0.7785,mean = 1.5,sigma = 1),
                            rpgnorm(n = n2,p = 0.7785,mean = 0,sigma = 1),
                            rpgnorm(n = n2,p = 0.7785,mean = -1.5,sigma = 1)))
  }
  for (i in 1:5) {
    conti_signal <- cbind(conti_signal,
                          c(rpgnorm(n = n1,p = 0.7785,mean = 0,sigma = 1),
                            rpgnorm(n = n2,p = 0.7785,mean = 1.5,sigma = 1),
                            rpgnorm(n = n2,p = 0.7785,mean = 0,sigma = 1)))
  }
  conti_noise <- NULL
  for (i in 1:(p-10)) {
    conti_noise <- cbind(conti_noise,
                         c(rpgnorm(n = n1+2*n2,p = 0.7785,mean = 0,sigma = 1)))
  }
  contivar <- data.frame(cbind(conti_signal, conti_noise))
  
  cate_v1 <- c(sample(1:2,size = n1,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.2,0.8)))
  cate_v2 <- c(sample(1:2,size = n1,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.2,0.8)))
  cate_v3 <- c(sample(1:3,size = n1,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.1,0.1,0.8)))
  cate_v4 <- c(sample(1:3,size = n1,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.1,0.1,0.8)))
  cate_v5 <- c(sample(1:4,size = n1,replace = T,prob = c(0.4,0.4,0.1,0.1)),
               sample(1:4,size = n2,replace = T,prob = c(0.25,0.25,0.25,0.25)),
               sample(1:4,size = n2,replace = T,prob = c(0.1,0.1,0.4,0.4)))
  
  cate_v6 <- c(sample(1:2,size = n1,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)))
  cate_v7 <- c(sample(1:2,size = n1,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)))
  cate_v8 <- c(sample(1:3,size = n1,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3))
  cate_v9 <- c(sample(1:3,size = n1,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3))
  cate_v10 <- c(sample(1:4,size = n1,replace = T,prob = c(0.25,0.25,0.25,0.25)),
                sample(1:4,size = n2,replace = T,prob = c(0.4,0.4,0.1,0.1)),
                sample(1:4,size = n2,replace = T,prob = c(0.25,0.25,0.25,0.25)))
  
  
  catevar <- data.frame(cate_v1,cate_v2,cate_v3,cate_v4,cate_v5,
                        cate_v6,cate_v7,cate_v8,cate_v9,cate_v10)
  
  cate_sig <- c(rep(1,length(unique(cate_v1))),
                rep(2,length(unique(cate_v2))),
                rep(3,length(unique(cate_v3))),
                rep(4,length(unique(cate_v4))),
                rep(5,length(unique(cate_v5))),
                rep(6,length(unique(cate_v6))),
                rep(7,length(unique(cate_v7))),
                rep(8,length(unique(cate_v8))),
                rep(9,length(unique(cate_v9))),
                rep(10,length(unique(cate_v10))))
  
  
  catevardum <- dummy_cols(catevar,select_columns = c('cate_v1','cate_v2','cate_v3','cate_v4','cate_v5',
                                                      'cate_v6','cate_v7','cate_v8','cate_v9','cate_v10'),
                           remove_first_dummy = F, remove_selected_columns = T)
  
  
  catenoi <- NULL
  cate_noi <- NULL
  for (i in 1:(p-10)) {
    temp_level <- sample(2:4,1)
    catenoi <- cbind(catenoi, factor(sample(1:(temp_level),n1+2*n2,replace = T)))
    cate_noi <- c(cate_noi, rep(i, temp_level))
  }
  catenoi <- as.data.frame(catenoi)
  catenoidum <- dummy_cols(catenoi,select_columns = colnames(catenoi),
                           remove_first_dummy = F, remove_selected_columns = T)
  
  simd <- cbind(contivar,catevardum,catenoidum)
  simgroup <- c(1:p,cate_sig+p,cate_noi+p+10)
  
  conti_df <- data.frame(cbind(conti_signal,conti_noise,catenoi,
                               cate_v1,cate_v2,cate_v3,cate_v4,cate_v5,
                               cate_v6,cate_v7,cate_v8,cate_v9,cate_v10))
  conti_df[,(p+1):(2*p)] <- lapply(conti_df[,(p+1):(2*p)], factor)
  
  return(list(simd=simd,
              simgroup=simgroup,
              conti_df=conti_df))
}


# K=3
# normal distribution
generatesimd <- function(n1, n2, p){
  conti_mu1 <- c(rep(1.5,5),rep(0,p-5))
  conti_mu2 <- c(rep(0,5),rep(1.5,5),rep(0,p-10))
  conti_mu3 <- c(rep(-1.5,5),rep(0,p-5))
  Sigma1 <- genPositiveDefMat(
    dim = p,
    covMethod = 'onion', 
    rangeVar = c(1.5,2.5))
  Sigma2 <- genPositiveDefMat(
    dim = p,
    covMethod = 'onion', 
    rangeVar = c(1.5,2.5))
  Sigma3 <- genPositiveDefMat(
    dim = p,
    covMethod = 'onion', 
    rangeVar = c(1.5,2.5))
  contivar <- as.data.frame(rbind(mvrnorm(n = n1, mu = conti_mu1, Sigma = Sigma1$Sigma), 
                                  mvrnorm(n = n2, mu = conti_mu2, Sigma = Sigma2$Sigma),
                                  mvrnorm(n = n2, mu = conti_mu3, Sigma = Sigma3$Sigma)))
  
  
  cate_v1 <- c(sample(1:2,size = n1,replace = T,prob = c(0.9,0.1)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.1,0.9)))
  cate_v2 <- c(sample(1:2,size = n1,replace = T,prob = c(0.9,0.1)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.1,0.9)))
  cate_v3 <- c(sample(1:3,size = n1,replace = T,prob = c(0.9,0.05,0.05)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.05,0.05,0.9)))
  cate_v4 <- c(sample(1:3,size = n1,replace = T,prob = c(0.9,0.05,0.05)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.05,0.05,0.9)))
  cate_v5 <- c(sample(1:4,size = n1,replace = T,prob = c(0.45,0.45,0.05,0.05)),
               sample(1:4,size = n2,replace = T,prob = c(0.25,0.25,0.25,0.25)),
               sample(1:4,size = n2,replace = T,prob = c(0.05,0.05,0.45,0.45)))
  
  cate_v6 <- c(sample(1:2,size = n1,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.9,0.1)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)))
  cate_v7 <- c(sample(1:2,size = n1,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.9,0.1)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)))
  cate_v8 <- c(sample(1:3,size = n1,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.9,0.05,0.05)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3))
  cate_v9 <- c(sample(1:3,size = n1,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.9,0.05,0.05)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3))
  cate_v10 <- c(sample(1:4,size = n1,replace = T,prob = c(0.25,0.25,0.25,0.25)),
                sample(1:4,size = n2,replace = T,prob = c(0.45,0.45,0.05,0.05)),
                sample(1:4,size = n2,replace = T,prob = c(0.25,0.25,0.25,0.25)))
  
  
  catevar <- data.frame(cate_v1,cate_v2,cate_v3,cate_v4,cate_v5,
                        cate_v6,cate_v7,cate_v8,cate_v9,cate_v10)
  
  cate_sig <- c(rep(1,length(unique(cate_v1))),
                rep(2,length(unique(cate_v2))),
                rep(3,length(unique(cate_v3))),
                rep(4,length(unique(cate_v4))),
                rep(5,length(unique(cate_v5))),
                rep(6,length(unique(cate_v6))),
                rep(7,length(unique(cate_v7))),
                rep(8,length(unique(cate_v8))),
                rep(9,length(unique(cate_v9))),
                rep(10,length(unique(cate_v10))))
  
  
  catevardum <- dummy_cols(catevar,select_columns = c('cate_v1','cate_v2','cate_v3','cate_v4','cate_v5',
                                                      'cate_v6','cate_v7','cate_v8','cate_v9','cate_v10'),
                           remove_first_dummy = F, remove_selected_columns = T)
  
  
  catenoi <- NULL
  cate_noi <- NULL
  for (i in 1:(p-10)) {
    temp_level <- sample(2:4,1)
    catenoi <- cbind(catenoi, factor(sample(1:(temp_level),n1+2*n2,replace = T)))
    cate_noi <- c(cate_noi, rep(i, temp_level))
  }
  catenoi <- as.data.frame(catenoi)
  catenoidum <- dummy_cols(catenoi,select_columns = colnames(catenoi),
                           remove_first_dummy = F, remove_selected_columns = T)
  
  simd <- cbind(contivar,catevardum,catenoidum)
  simgroup <- c(1:p,cate_sig+p,cate_noi+p+10)
  
  conti_df <- data.frame(cbind(contivar,catenoi,
                               cate_v1,cate_v2,cate_v3,cate_v4,cate_v5,
                               cate_v6,cate_v7,cate_v8,cate_v9,cate_v10))
  conti_df[,(p+1):(2*p)] <- lapply(conti_df[,(p+1):(2*p)], factor)
  
  return(list(simd=simd,
              simgroup=simgroup,
              conti_df=conti_df))
}


# K=4
# p-generalized normal-polynomial distribution.
generatesimd <- function(n1, n2, p){
  #conti_mu <- c(rep(1,5),rep(0,195))
  #contivar <- as.data.frame(rbind(mvrnorm(n = 50, mu = conti_mu, Sigma = diag(nrow = 200)), mvrnorm(n = 50, mu = rep(0,200), Sigma = diag(nrow = 200))))
  conti_signal <- NULL
  for (i in 1:5) {
    conti_signal <- cbind(conti_signal,
                          c(rpgnorm(n = n1,p = 0.7785,mean = 1.5,sigma = 1),
                            rpgnorm(n = n2,p = 0.7785,mean = 0,sigma = 1),
                            rpgnorm(n = n2,p = 0.7785,mean = -1.5,sigma = 1),
                            rpgnorm(n = n2,p = 0.7785,mean = 0,sigma = 1)))
  }
  for (i in 1:5) {
    conti_signal <- cbind(conti_signal,
                          c(rpgnorm(n = n1,p = 0.7785,mean = 0,sigma = 1),
                            rpgnorm(n = n2,p = 0.7785,mean = 1.5,sigma = 1),
                            rpgnorm(n = n2,p = 0.7785,mean = 0,sigma = 1),
                            rpgnorm(n = n2,p = 0.7785,mean = -1.5,sigma = 1)))
  }
  conti_noise <- NULL
  for (i in 1:(p-10)) {
    conti_noise <- cbind(conti_noise,
                         c(rpgnorm(n = n1+3*n2,p = 0.7785,mean = 0,sigma = 1)))
  }
  contivar <- data.frame(cbind(conti_signal, conti_noise))
  
  cate_v1 <- c(sample(1:2,size = n1,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.2,0.8)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)))
  cate_v2 <- c(sample(1:2,size = n1,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.2,0.8)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)))
  cate_v3 <- c(sample(1:3,size = n1,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.1,0.1,0.8)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3))
  cate_v4 <- c(sample(1:3,size = n1,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.1,0.1,0.8)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3))
  cate_v5 <- c(sample(1:4,size = n1,replace = T,prob = c(0.4,0.4,0.1,0.1)),
               sample(1:4,size = n2,replace = T,prob = c(0.25,0.25,0.25,0.25)),
               sample(1:4,size = n2,replace = T,prob = c(0.1,0.1,0.4,0.4)),
               sample(1:4,size = n2,replace = T,prob = c(0.25,0.25,0.25,0.25)))
  
  cate_v6 <- c(sample(1:2,size = n1,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.2,0.8)))
  cate_v7 <- c(sample(1:2,size = n1,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.2,0.8)))
  cate_v8 <- c(sample(1:3,size = n1,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.1,0.1,0.8)))
  cate_v9 <- c(sample(1:3,size = n1,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.1,0.1,0.8)))
  cate_v10 <- c(sample(1:4,size = n1,replace = T,prob = c(0.25,0.25,0.25,0.25)),
                sample(1:4,size = n2,replace = T,prob = c(0.4,0.4,0.1,0.1)),
                sample(1:4,size = n2,replace = T,prob = c(0.25,0.25,0.25,0.25)),
                sample(1:4,size = n2,replace = T,prob = c(0.1,0.1,0.4,0.4)))
  
  
  catevar <- data.frame(cate_v1,cate_v2,cate_v3,cate_v4,cate_v5,
                        cate_v6,cate_v7,cate_v8,cate_v9,cate_v10)
  
  cate_sig <- c(rep(1,length(unique(cate_v1))),
                rep(2,length(unique(cate_v2))),
                rep(3,length(unique(cate_v3))),
                rep(4,length(unique(cate_v4))),
                rep(5,length(unique(cate_v5))),
                rep(6,length(unique(cate_v6))),
                rep(7,length(unique(cate_v7))),
                rep(8,length(unique(cate_v8))),
                rep(9,length(unique(cate_v9))),
                rep(10,length(unique(cate_v10))))
  
  
  catevardum <- dummy_cols(catevar,select_columns = c('cate_v1','cate_v2','cate_v3','cate_v4','cate_v5',
                                                      'cate_v6','cate_v7','cate_v8','cate_v9','cate_v10'),
                           remove_first_dummy = F, remove_selected_columns = T)
  
  
  catenoi <- NULL
  cate_noi <- NULL
  for (i in 1:(p-10)) {
    temp_level <- sample(2:4,1)
    catenoi <- cbind(catenoi, factor(sample(1:(temp_level),n1+3*n2,replace = T)))
    cate_noi <- c(cate_noi, rep(i, temp_level))
  }
  catenoi <- as.data.frame(catenoi)
  catenoidum <- dummy_cols(catenoi,select_columns = colnames(catenoi),
                           remove_first_dummy = F, remove_selected_columns = T)
  
  simd <- cbind(contivar,catevardum,catenoidum)
  simgroup <- c(1:p,cate_sig+p,cate_noi+p+10)
  
  conti_df <- data.frame(cbind(conti_signal,conti_noise,catenoi,
                               cate_v1,cate_v2,cate_v3,cate_v4,cate_v5,
                               cate_v6,cate_v7,cate_v8,cate_v9,cate_v10))
  conti_df[,(p+1):(2*p)] <- lapply(conti_df[,(p+1):(2*p)], factor)
  
  return(list(simd=simd,
              simgroup=simgroup,
              conti_df=conti_df))
}



# K=4
# normal distribution
generatesimd <- function(n1, n2, p){
  conti_mu1 <- c(rep(1.5,5),rep(0,p-5))
  conti_mu2 <- c(rep(0,5),rep(1.5,5),rep(0,p-10))
  conti_mu3 <- c(rep(-1.5,5),rep(0,p-5))
  conti_mu4 <- c(rep(0,5),rep(-1.5,5),rep(0,p-10))
  Sigma1 <- genPositiveDefMat(
    dim = p,
    covMethod = 'onion', 
    rangeVar = c(1.5,2.5))
  Sigma2 <- genPositiveDefMat(
    dim = p,
    covMethod = 'onion', 
    rangeVar = c(1.5,2.5))
  Sigma3 <- genPositiveDefMat(
    dim = p,
    covMethod = 'onion', 
    rangeVar = c(1.5,2.5))
  Sigma4 <- genPositiveDefMat(
    dim = p,
    covMethod = 'onion', 
    rangeVar = c(1.5,2.5))
  contivar <- as.data.frame(rbind(mvrnorm(n = n1, mu = conti_mu1, Sigma = Sigma1$Sigma), 
                                  mvrnorm(n = n2, mu = conti_mu2, Sigma = Sigma2$Sigma),
                                  mvrnorm(n = n2, mu = conti_mu3, Sigma = Sigma3$Sigma),
                                  mvrnorm(n = n2, mu = conti_mu4, Sigma = Sigma4$Sigma)))
  
  cate_v1 <- c(sample(1:2,size = n1,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.2,0.8)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)))
  cate_v2 <- c(sample(1:2,size = n1,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.2,0.8)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)))
  cate_v3 <- c(sample(1:3,size = n1,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.1,0.1,0.8)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3))
  cate_v4 <- c(sample(1:3,size = n1,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.1,0.1,0.8)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3))
  cate_v5 <- c(sample(1:4,size = n1,replace = T,prob = c(0.4,0.4,0.1,0.1)),
               sample(1:4,size = n2,replace = T,prob = c(0.25,0.25,0.25,0.25)),
               sample(1:4,size = n2,replace = T,prob = c(0.1,0.1,0.4,0.4)),
               sample(1:4,size = n2,replace = T,prob = c(0.25,0.25,0.25,0.25)))
  
  cate_v6 <- c(sample(1:2,size = n1,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.2,0.8)))
  cate_v7 <- c(sample(1:2,size = n1,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.8,0.2)),
               sample(1:2,size = n2,replace = T,prob = c(0.5,0.5)),
               sample(1:2,size = n2,replace = T,prob = c(0.2,0.8)))
  cate_v8 <- c(sample(1:3,size = n1,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.1,0.1,0.8)))
  cate_v9 <- c(sample(1:3,size = n1,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.8,0.1,0.1)),
               sample(1:3,size = n2,replace = T,prob = c(1,1,1)/3),
               sample(1:3,size = n2,replace = T,prob = c(0.1,0.1,0.8)))
  cate_v10 <- c(sample(1:4,size = n1,replace = T,prob = c(0.25,0.25,0.25,0.25)),
                sample(1:4,size = n2,replace = T,prob = c(0.4,0.4,0.1,0.1)),
                sample(1:4,size = n2,replace = T,prob = c(0.25,0.25,0.25,0.25)),
                sample(1:4,size = n2,replace = T,prob = c(0.1,0.1,0.4,0.4)))
  
  
  catevar <- data.frame(cate_v1,cate_v2,cate_v3,cate_v4,cate_v5,
                        cate_v6,cate_v7,cate_v8,cate_v9,cate_v10)
  
  cate_sig <- c(rep(1,length(unique(cate_v1))),
                rep(2,length(unique(cate_v2))),
                rep(3,length(unique(cate_v3))),
                rep(4,length(unique(cate_v4))),
                rep(5,length(unique(cate_v5))),
                rep(6,length(unique(cate_v6))),
                rep(7,length(unique(cate_v7))),
                rep(8,length(unique(cate_v8))),
                rep(9,length(unique(cate_v9))),
                rep(10,length(unique(cate_v10))))
  
  
  catevardum <- dummy_cols(catevar,select_columns = c('cate_v1','cate_v2','cate_v3','cate_v4','cate_v5',
                                                      'cate_v6','cate_v7','cate_v8','cate_v9','cate_v10'),
                           remove_first_dummy = F, remove_selected_columns = T)
  
  
  catenoi <- NULL
  cate_noi <- NULL
  for (i in 1:(p-10)) {
    temp_level <- sample(2:4,1)
    catenoi <- cbind(catenoi, factor(sample(1:(temp_level),n1+3*n2,replace = T)))
    cate_noi <- c(cate_noi, rep(i, temp_level))
  }
  catenoi <- as.data.frame(catenoi)
  catenoidum <- dummy_cols(catenoi,select_columns = colnames(catenoi),
                           remove_first_dummy = F, remove_selected_columns = T)
  
  simd <- cbind(contivar,catevardum,catenoidum)
  simgroup <- c(1:p,cate_sig+p,cate_noi+p+10)
  
  conti_df <- data.frame(cbind(conti_signal,conti_noise,catenoi,
                               cate_v1,cate_v2,cate_v3,cate_v4,cate_v5,
                               cate_v6,cate_v7,cate_v8,cate_v9,cate_v10))
  conti_df[,(p+1):(2*p)] <- lapply(conti_df[,(p+1):(2*p)], factor)
  
  return(list(simd=simd,
              simgroup=simgroup,
              conti_df=conti_df))
}



#############################################################
#####                                                   #####
#####                                                   #####
#####                    Main Functions                 #####
#####                                                   #####
#####                                                   #####
#############################################################

################## Preparatory functions ####################
library(Rcpp)
library(fastDummies)
library(sparcl)
library(cluster)
library(VarSelLCM)
library(flexclust)
library(pgnorm)
library(clusterGeneration)

########### Calculation of co-occurence distance ############
sourceCpp('co_two.cpp')
sourceCpp('dist_co.cpp')

# Calculation of the adjusted DBI for categorical variables #
GetcateDBI <- function(x, Cs, vg){
    ##### Input: #####
    #   x: Data with categorical variables.
    #   Cs: The clustering labels.
    #   vg: The indicator of the relationships among dummy variables.
    ##### Output: #####
    #   dbi_cate_v: The adjusted DBI for each categorical variable.
    #   dist_cate: The distance matrix for each samples on categorical variables.
    #   t_v: The detail co-occurence distance for each level.
    #   mdbi_cate: The categorical part of mDBI.
  dist_x <- dist_co(as.matrix(x))
  dist_m <- dist_x$d
  v_dbi_d <- NULL
  m_dbi_d <- NULL
  for (v in unique(vg)) {
    v_dcc <- matrix(0,nrow = max(Cs),ncol = max(Cs))
    which_v <- which(vg==v)
    dist_v <- as.matrix(dist(as.matrix(x[,which_v]),method = 'manhattan'))
    for (k in 1:(max(Cs)-1)) {
      if(sum(Cs==k)==1){
        temp_ki <- which(Cs==k)
        temp_sigmak <- 0
      } else{
        temp_ki <- which(Cs==k)
        temp_sigmak <- sum(dist_v[temp_ki,temp_ki])/length(temp_ki)/(length(temp_ki)-1)
      }
      for (l in (k+1):max(Cs)) {
        if(sum(Cs==l)==1){
          temp_li <- which(Cs==l)
          temp_sigmal <- 0
        } else{
          temp_li <- which(Cs==l)
          temp_sigmal <- sum(dist_v[temp_li,temp_li])/length(temp_li)/(length(temp_li)-1)
        }
        if((temp_sigmak+temp_sigmal)==0){
          temp_dcc <- sum(dist_v[temp_ki,temp_li])
          v_dcc[k,l] <- v_dcc[l,k] <- temp_dcc
        } else{
          temp_dcc <- sum(dist_v[temp_ki,temp_li])/length(temp_ki)/(length(temp_li))
          v_dcc[k,l] <- v_dcc[l,k] <- temp_dcc/(temp_sigmak+temp_sigmal)
        }
        
        #temp_dbi <- temp_dcc/(temp_sigmak+temp_sigmal)
        #if(temp_dbi>v_dcc){
        #v_dcc <- temp_dbi
        #}
      }
    }
    #v_dbi_d <- c(v_dbi_d, max(v_dcc))
    v_dbi_d <- c(v_dbi_d, mean(apply(v_dcc,1,max)))
    #m_dbi_d <- c(m_dbi_d, mean(apply(1/v_dcc,1,min)))
    m_dbi_d <- c(m_dbi_d, 1/(max(v_dcc)))
    #m_dbi_d <- c(m_dbi_d, (sum(1/(v_dcc+diag(1,nrow=nrow(v_dcc))))-nrow(v_dcc))/(nrow(v_dcc)*(nrow(v_dcc)-1)))
    #v_dbi_d <- c(v_dbi_d, mean(apply((v_dcc+diag(max(v_dcc),nrow = nrow(v_dcc))),1,min)))
  }
  
  dist_per_v <- aggregate(dist_x$v,by=list(vg),FUN=max)[,2]
  #dbi_v <- dist_per_v*v_dbi_d
  dbi_v <- v_dbi_d
  result <- list(dbi_cate_v = dbi_v,
                 dist_cate = dist_m,
                 t_v = dist_x$v,
                 mdbi_cate = m_dbi_d)
  return(result)
}


# Calculation of the adjusted DBI for continuous variables #
GetcontDBI <- function(x,Cs){
  
  ##### Input: #####
  #   x: Data with continuous variables.
  #   Cs: The clustering labels.
  ##### Output: #####
  #   dbi_cont_v: The adjusted DBI for each continuous variable.
  #   dist_cont: The distance matrix for each samples on continuous variables.
  #   mdbi_cont: The continuous part of mDBI.

  
  v_dbi_d <- NULL
  m_dbi_d <- NULL
  for (v in 1:ncol(x)) {
    v_dcc <- matrix(0,nrow = max(Cs),ncol = max(Cs))
    dist_v <- as.matrix(dist(x[,v]))
    for(k in 1:(max(Cs)-1)){
      if(sum(Cs==k)==1){
        temp_ki <- which(Cs==k)
        temp_sigmak <- 0
      } else{
        temp_ki <- which(Cs==k)
        temp_sigmak <- sum(dist_v[temp_ki,temp_ki])/length(temp_ki)/(length(temp_ki)-1)
      }
      for (l in (k+1):max(Cs)) {
        if(sum(Cs==l)==1){
          temp_li <- which(Cs==l)
          temp_sigmal <- 0
        } else{
          temp_li <- which(Cs==l)
          temp_sigmal <- sum(dist_v[temp_li,temp_li])/length(temp_li)/(length(temp_li)-1)
        }
        if((temp_sigmak+temp_sigmal)==0){
          temp_dcc <- sum(dist_v[temp_ki,temp_li])
          v_dcc[k,l] <- v_dcc[l,k] <- temp_dcc
        } else{
          temp_dcc <- sum(dist_v[temp_ki,temp_li])/length(temp_ki)/(length(temp_li))
          v_dcc[k,l] <- v_dcc[l,k] <- temp_dcc/(temp_sigmak+temp_sigmal)
        }
        #temp_dbi <- temp_dcc/(temp_sigmak+temp_sigmal)
        #if(temp_dbi>v_dcc){
        #v_dcc <- temp_dbi
        #}
      }
    }
    #v_dbi_d <- c(v_dbi_d, max(v_dcc))
    v_dbi_d <- c(v_dbi_d, mean(apply(v_dcc,1,max)))
    #m_dbi_d <- c(m_dbi_d, mean(apply(1/v_dcc,1,min)))
    m_dbi_d <- c(m_dbi_d, 1/(max(v_dcc)))
    #m_dbi_d <- c(m_dbi_d, (sum(1/(v_dcc+diag(1,nrow=nrow(v_dcc))))-nrow(v_dcc))/(nrow(v_dcc)*(nrow(v_dcc)-1)))
    #v_dbi_d <- c(v_dbi_d, mean(apply((v_dcc+diag(max(v_dcc),nrow = nrow(v_dcc))),1,min)))
  }
  result <- list(dbi_cont_v = v_dbi_d,
                 dist_cont = as.matrix(dist(x)),
                 mdbi_cont = m_dbi_d)
  return(result)
}

############### Calculate CER of results ####################
cal_CER <- function(predict, real){
  P <- matrix(0,nrow = length(predict), ncol = length(predict))
  for (i in 1:(length(predict)-1)) {
    for (j in (i+1):length(predict)) {
      if(predict[i]==predict[j]) P[i,j] <- 1
    }
  }
  Q <- matrix(0,nrow = length(real), ncol = length(real))
  for (i in 1:(length(real)-1)) {
    for (j in (i+1):length(real)) {
      if(real[i]==real[j]) Q[i,j] <- 1
    }
  }
  CER <- sum(abs(P-Q))*2/(length(predict)*(length(predict)-1))
  return(CER)
}



############# Simple calculation functions ##################

BinarySearch <- function(argu,sumabs){
  if(l2n(argu)==0 || sum(abs(argu/l2n(argu)))<=sumabs) return(0)
  lam1 <- 0
  lam2 <- max(abs(argu))-1e-5
  iter <- 1
  while(iter<=15 && (lam2-lam1)>(1e-4)){
    su <- soft(argu,(lam1+lam2)/2)
    if(sum(abs(su/l2n(su)))<sumabs){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    iter <- iter+1
  }
  return((lam1+lam2)/2)
}

soft <- function(x,d){
  return(sign(x)*pmax(0, abs(x)-d))
}

l2n <- function(vec){
  return(sqrt(sum(vec^2)))
}


############### Update function for weights #################
UpdateWs_DBI <- function(x, y, vg, Cs, w, v){
  ##### Input: #####
  #   x: Data with continuous variables.
  #   y: Data with categorical variables.
  #   Cs: The clustering labels.
  #   vg: The indicator of the relationships among dummy variables.
  #   w: The weights for categorical variables.
  #   v: The weights for continuous variables.
  ##### Output: #####
  #   ws: The weights of all variables. 
  #   vs: The weights of continuous variables. 
  #   ws_cate: The The weights of categorical variables. 
  #   Cs: The clustering results.
  #   t_v: 
  
  cont <- GetcontDBI(x, Cs)
  cate <- GetcateDBI(y, Cs, vg)
  DBI.perfeaturegroup <- cate$dbi_cate_v
  DBI.perfeaturegroup <- DBI.perfeaturegroup - 0.5
  DBI.perfeaturegroup[DBI.perfeaturegroup<0] <- 0
  lam <- BinarySearch(DBI.perfeaturegroup, w)
  ws.unscaled <- soft(DBI.perfeaturegroup,lam)
  ws <- ws.unscaled/l2n(ws.unscaled)
  ws_cate <- ws
  DBI.perfeature <- cont$dbi_cont_v
  DBI.perfeature <- DBI.perfeature - 0.5
  DBI.perfeature[DBI.perfeature<0] <- 0
  lam <- BinarySearch(DBI.perfeature, v)
  vs.unscaled <- soft(DBI.perfeature, lam)
  vs <- vs.unscaled/l2n(vs.unscaled)
  ws_cont <- vs
  rep_time <- table(vg)
  ws <- c(vs,rep(ws,rep_time))
  return(list(ws=ws,vs=vs,
              ws_cate=ws_cate,ws_cont=ws_cont,
              t_v = cate$t_v))
}



########### Update function for cluster lables ##############
UpdateCs_DBI <- function(x, K, vs, Cs, t_v){
  ##### Input: #####
  #   x: Data with all variables.
  #   K: The number of clusters.
  #   Cs: The clustering labels.
  #   vs: The weights of all variables. 
  #   t_v: The detail co-occurence distance for each level.
  ##### Output: #####
  #   The updated clustering results.
  if(sum(vs!=0)==1){
    x <- x[,vs!=0]
    z <- x*vs[vs!=0]
    nrowz <- length(z)
    mus <- NULL
    if(!is.null(Cs)){
      for(k in unique(Cs)){
        if(sum(Cs==k)>1) mus <- c(mus, mean(z[Cs==k]))
        if(sum(Cs==k)==1) mus <- c(mus, z[Cs==k])
      }
    }
    if(is.null(mus)){
      km <- kmeans(as.matrix(z), centers=K, nstart=10)
    } else {
      distmat <- as.matrix(dist(c(z, mus)))[1:nrowz, (nrowz+1):(nrowz+K)]
      nearest <- apply(distmat, 1, which.min)
      if(length(unique(nearest))==K){
        km <- kmeans(as.matrix(z), K)
      } else {
        km <- kmeans(as.matrix(z), centers=K, nstart=10)
      }
    }
  } else{
    x <- x[,vs!=0]
    #z <- x#sweep(x, 2, sqrt(vs[vs!=0]), "*")
    z <- sweep(x, 2, sqrt(t_v[vs!=0]), "*")
    nrowz <- nrow(z)
    mus <- NULL
    if(!is.null(Cs)){
      for(k in unique(Cs)){
        if(sum(Cs==k)>1) mus <- rbind(mus, apply(z[Cs==k,],2,mean))
        if(sum(Cs==k)==1) mus <- rbind(mus, z[Cs==k,])
      }
    }
    if(is.null(mus)){
      km <- kmeans(z, centers=K, nstart=10)
    } else {
      distmat <- as.matrix(dist(rbind(z, mus)))[1:nrowz, (nrowz+1):(nrowz+K)]
      nearest <- apply(distmat, 1, which.min)
      if(length(unique(nearest))==K){
        km <- kmeans(z, centers=mus)
      } else {
        km <- kmeans(z, centers=K, nstart=10)
      }
    }
  }
  
  return(km$cluster)
}




#################### The proposed method. ######################
###### x is the dataset containing continuous variables. #######
###### y is the dataset containing categorical variables. ######
###### K is the number of clusters. ############################
###### wbounds is the penalty parameter for continuous. ########
###### vbounds is the penalty parameter for categorical. #######
###### group indicates the relationships among dummy variables.#

`COC` <-
  function(x, y, K=NULL, wbounds=NULL, vbounds=NULL, group,
           nstart=20, silent=FALSE, maxiter=6, centers=NULL){
    
    ##### Input: #####
    #   x: Data with continuous variables.
    #   y: Data with categorical variables.
    #   K: The number of clusters.
    #   wbounds: The penalty parameter controlling continuous variables.
    #   vbounds: The penalty parameter controlling categorical variables.
    #   group: The indicator of the relationships among dummy variables.
    ##### Output: #####
    #   vs: The weights of continuous and categorical (dummy) variables. 
    #   Cs: The clustering results.
    #   DBI.perfeature: The adjusted DBI for each continuous variable.
    #   DBI.perfeaturegroup: The adjusted DBI for each categorical variable.
    #   mDBI: The mDBI critrion of each iteration.
    
    
    
    ### Check the input
    if(is.null(K) && is.null(centers)) stop("Must provide either K or centers.")
    if(!is.null(K) && !is.null(centers)){
      if(nrow(centers)!=K) stop("If K and centers both are provided, then nrow(centers) must equal K!!!")
      if(nrow(centers)==K) K <- NULL
    }
    if(!is.null(centers) && ncol(centers)!=(ncol(x)+ncol(y))) stop("If centers is provided, then ncol(centers) must equal ncol(x).")
    if(is.null(wbounds)) wbounds <- seq(1.1, sqrt(max(group)), len=20)
    if(is.null(vbounds)) vbounds <- (wbounds+1)/2
    if(min(wbounds)<=1) stop("wbounds should be greater than 1")
    wbounds <- c(wbounds) # In case wbounds is a single number, turn it into a vector
    out <- list()
    
    total_data <- cbind(x,y)
    t_v <- rep(1,ncol(x)+ncol(y))
    
    if(!is.null(K)) Cs <- kmeans(total_data, centers=K, nstart=nstart)$cluster
    if(is.null(K)) Cs <- kmeans(total_data, centers=centers)$cluster
    for(i in 1:length(wbounds)){
      if(length(wbounds)>1 && !silent) cat(i,fill=FALSE)
      vs <- rep(1/sqrt(ncol(x)+ncol(y)), (ncol(x)+ncol(y))) # Start with equal weights on each feature
      vs.old <- rnorm(ncol(x)+ncol(y))
      store.bcss.vs <- NULL
      store.mdbi <- NULL
      niter <- 0
      while((sum(abs(vs-vs.old))/sum(abs(vs.old)))>1e-4 && niter<maxiter){
        if(!silent) cat(niter, fill=FALSE)
        niter <- niter+1
        vs.old <- vs
        if(!is.null(K)){
          if(niter>1) Cs <- UpdateCs_DBI(total_data, K, vs, Cs, t_v) # if niter=1, no need to update!!
        } else {
          if(niter>1) Cs <- UpdateCs_DBI(total_data, nrow(centers), vs, Cs) # if niter=1, no need to update!!
        }
        #vs <- UpdateWs_DBI(x, Cs, group, category, wbounds[i], vbounds[i])$vs
        ws_dbi <- UpdateWs_DBI(x, y, group, Cs, wbounds[i], vbounds[i])
        t_v <- c(rep(1,ncol(x)),ws_dbi$t_v)
        vs <- ws_dbi$ws
        ws_cate <- rep(ws_dbi$ws_cate,table(group))
        ws_cont <- ws_dbi$ws_cont
        vg <- rep(1:length(unique(group[ws_cate>0])),table(group[ws_cate>0]))
        mdbi_temp <- (sum(GetcontDBI(x, Cs)$mdbi_cont*ws_cont)/sum(ws_cont)
                      +sum(rep(GetcateDBI(y, Cs, group)$mdbi_cate,table(group))*ws_cate)/sum(ws_cate))
        store.mdbi <- c(store.mdbi, mdbi_temp)
        store.bcss.vs <- c(store.bcss.vs, 
                           mdbi_temp+((length(unique(Cs)))^(1/(mdbi_temp*(max(1-0.1*(length(unique(Cs))-1),0.6))))))
        final_dbi <- c(GetcontDBI(x[,ws_cont>0], Cs)$dbi_cont_v,rep(GetcateDBI(y[,ws_cate>0], Cs, vg)$dbi_cate_v,table(vg)))
        #store.bcss.vs <- c(store.bcss.vs, sum(1/final_dbi))
      }
      out[[i]] <- list(vs=vs, Cs=Cs,
                       DBI.perfeaturegroup=GetcateDBI(y, Cs, group)$dbi_cate_v,
                       DBI.perfeature=GetcontDBI(x, Cs)$dbi_cont_v,
                       crit=store.bcss.vs, wbound=wbounds[i],
                       mDBI=store.mdbi)
    }
    if(!silent) cat(fill=TRUE)
    #  if(length(wbounds)==1){
    #    out <- out[[1]]
    #    class(out) <- "kmeanssparse"
    #    return(out)
    #  }
    #  class(out) <- "multikmeanssparse"
    #  return(out)
    class(out) <- "KMeansSparseCluster"
    return(out)
  }




############### Update function for weights #################
UpdateWs_SC <- function(x, y, vg, Cs, w, v){
  cont <- GetcontSC(x, Cs)
  cate <- GetcateSC(y, Cs, vg)
  DBI.perfeaturegroup <- cate$sc_cate_v
  DBI.perfeaturegroup <- DBI.perfeaturegroup + 0.5
  #DBI.perfeaturegroup[DBI.perfeaturegroup<0] <- 0
  lam <- BinarySearch(DBI.perfeaturegroup, w)
  ws.unscaled <- soft(DBI.perfeaturegroup,lam)
  ws <- ws.unscaled/l2n(ws.unscaled)
  ws_cate <- ws
  DBI.perfeature <- cont$sc_cont_v
  DBI.perfeature <- DBI.perfeature + 0.5
  #DBI.perfeature[DBI.perfeature<0] <- 0
  lam <- BinarySearch(DBI.perfeature, v)
  vs.unscaled <- soft(DBI.perfeature, lam)
  vs <- vs.unscaled/l2n(vs.unscaled)
  ws_cont <- vs
  rep_time <- table(vg)
  ws <- c(vs,rep(ws,rep_time))
  return(list(ws=ws,vs=vs,
              ws_cate=ws_cate,ws_cont=ws_cont,
              t_v = cate$t_v))
}


########### Update function for cluster lables ##############
UpdateCs_SC <- function(x, K, vs, Cs, t_v){
  if(sum(vs!=0)==1){
    x <- x[,vs!=0]
    z <- x*vs[vs!=0]
    nrowz <- length(z)
    mus <- NULL
    if(!is.null(Cs)){
      for(k in unique(Cs)){
        if(sum(Cs==k)>1) mus <- c(mus, mean(z[Cs==k]))
        if(sum(Cs==k)==1) mus <- c(mus, z[Cs==k])
      }
    }
    if(is.null(mus)){
      km <- kmeans(as.matrix(z), centers=K, nstart=10)
    } else {
      distmat <- as.matrix(dist(c(z, mus)))[1:nrowz, (nrowz+1):(nrowz+K)]
      nearest <- apply(distmat, 1, which.min)
      if(length(unique(nearest))==K){
        km <- kmeans(as.matrix(z), K)
      } else {
        km <- kmeans(as.matrix(z), centers=K, nstart=10)
      }
    }
  } else{
    x <- x[,vs!=0]
    #z <- x#sweep(x, 2, sqrt(vs[vs!=0]), "*")
    z <- sweep(x, 2, sqrt(t_v[vs!=0]), "*")
    nrowz <- nrow(z)
    mus <- NULL
    if(!is.null(Cs)){
      for(k in unique(Cs)){
        if(sum(Cs==k)>1) mus <- rbind(mus, apply(z[Cs==k,],2,mean))
        if(sum(Cs==k)==1) mus <- rbind(mus, z[Cs==k,])
      }
    }
    if(is.null(mus)){
      km <- kmeans(z, centers=K, nstart=10)
    } else {
      distmat <- as.matrix(dist(rbind(z, mus)))[1:nrowz, (nrowz+1):(nrowz+K)]
      nearest <- apply(distmat, 1, which.min)
      if(length(unique(nearest))==K){
        km <- kmeans(z, centers=mus)
      } else {
        km <- kmeans(z, centers=K, nstart=10)
      }
    }
  }
  
  return(km$cluster)
}


# Calculation of the silhouette score for categorical variables #
GetcateSC <- function(x, Cs, vg){
  #dbi_c <- matrix(0,nrow = max(Cs),ncol = max(Cs))
  dist_x <- dist_co(as.matrix(x))
  dist_m <- dist_x$d
  v_dbi_d <- NULL
  sc_v <- NULL
  for (v in unique(vg)) {
    v_dcc <- matrix(0,nrow = max(Cs),ncol = max(Cs))
    which_v <- which(vg==v)
    dist_v <- as.matrix(dist(as.matrix(x[,which_v]),method = 'manhattan'))
    temp_SC <- silhouette(Cs,dist_v)
    sc_v <- c(sc_v, mean(temp_SC[,3]))
    for (k in 1:(max(Cs)-1)) {
      if(sum(Cs==k)==1){
        temp_ki <- which(Cs==k)
        temp_sigmak <- 0
      } else{
        temp_ki <- which(Cs==k)
        temp_sigmak <- sum(dist_v[temp_ki,temp_ki])/length(temp_ki)/(length(temp_ki)-1)
      }
      for (l in (k+1):max(Cs)) {
        if(sum(Cs==l)==1){
          temp_li <- which(Cs==l)
          temp_sigmal <- 0
        } else{
          temp_li <- which(Cs==l)
          temp_sigmal <- sum(dist_v[temp_li,temp_li])/length(temp_li)/(length(temp_li)-1)
        }
        if((temp_sigmak+temp_sigmal)==0){
          temp_dcc <- sum(dist_v[temp_ki,temp_li])
          v_dcc[k,l] <- v_dcc[l,k] <- temp_dcc
        } else{
          temp_dcc <- sum(dist_v[temp_ki,temp_li])/length(temp_ki)/(length(temp_li))
          v_dcc[k,l] <- v_dcc[l,k] <- temp_dcc/(temp_sigmak+temp_sigmal)
        }
        
        #temp_dbi <- temp_dcc/(temp_sigmak+temp_sigmal)
        #if(temp_dbi>v_dcc){
        #v_dcc <- temp_dbi
        #}
      }
    }
    #v_dbi_d <- c(v_dbi_d, max(v_dcc))
    v_dbi_d <- c(v_dbi_d, mean(apply(v_dcc,1,max)))
    #v_dbi_d <- c(v_dbi_d, mean(apply((v_dcc+diag(max(v_dcc),nrow = nrow(v_dcc))),1,min)))
  }
  
  dist_per_v <- aggregate(dist_x$v,by=list(vg),FUN=max)[,2]
  #dbi_v <- dist_per_v*v_dbi_d
  dbi_v <- v_dbi_d
  result <- list(dbi_cate_v = dbi_v,
                 dist_cate = dist_m,
                 t_v = dist_x$v,
                 sc_cate_v = sc_v)
  return(result)
}

# Calculation of the silhouette score for continuous variables #
GetcontSC <- function(x,Cs){
  v_dbi_d <- NULL
  sc_v <- NULL
  for (v in 1:NCOL(x)) {
    v_dcc <- matrix(0,nrow = max(Cs),ncol = max(Cs))
    dist_v <- as.matrix(dist(x[,v]))
    temp_SC <- silhouette(Cs,dist_v)
    sc_v <- c(sc_v, mean(temp_SC[,3]))
    for(k in 1:(max(Cs)-1)){
      if(sum(Cs==k)==1){
        temp_ki <- which(Cs==k)
        temp_sigmak <- 0
      } else{
        temp_ki <- which(Cs==k)
        temp_sigmak <- sum(dist_v[temp_ki,temp_ki])/length(temp_ki)/(length(temp_ki)-1)
      }
      for (l in (k+1):max(Cs)) {
        if(sum(Cs==l)==1){
          temp_li <- which(Cs==l)
          temp_sigmal <- 0
        } else{
          temp_li <- which(Cs==l)
          temp_sigmal <- sum(dist_v[temp_li,temp_li])/length(temp_li)/(length(temp_li)-1)
        }
        if((temp_sigmak+temp_sigmal)==0){
          temp_dcc <- sum(dist_v[temp_ki,temp_li])
          v_dcc[k,l] <- v_dcc[l,k] <- temp_dcc
        } else{
          temp_dcc <- sum(dist_v[temp_ki,temp_li])/length(temp_ki)/(length(temp_li))
          v_dcc[k,l] <- v_dcc[l,k] <- temp_dcc/(temp_sigmak+temp_sigmal)
        }
        #temp_dbi <- temp_dcc/(temp_sigmak+temp_sigmal)
        #if(temp_dbi>v_dcc){
        #v_dcc <- temp_dbi
        #}
      }
    }
    #v_dbi_d <- c(v_dbi_d, max(v_dcc))
    v_dbi_d <- c(v_dbi_d, mean(apply(v_dcc,1,max)))
    #v_dbi_d <- c(v_dbi_d, mean(apply((v_dcc+diag(max(v_dcc),nrow = nrow(v_dcc))),1,min)))
  }
  result <- list(dbi_cont_v = v_dbi_d,
                 dist_cont = as.matrix(dist(x)),
                 sc_cont_v = sc_v)
  return(result)
}


######## The proposed method with silhouette score. ############
###### x is the dataset containing continuous variables. #######
###### y is the dataset containing categorical variables. ######
###### K is the number of clusters. ############################
###### wbounds is the penalty parameter for continuous. ########
###### vbounds is the penalty parameter for categorical. #######
###### group indicates the relationships among dummy variables.#


`COSC` <-
  function(x, y, K=NULL, wbounds=NULL, vbounds=NULL, group,
           nstart=20, silent=FALSE, maxiter=6, centers=NULL){
    ###检查函数输入是否有错误
    if(is.null(K) && is.null(centers)) stop("Must provide either K or centers.")
    if(!is.null(K) && !is.null(centers)){
      if(nrow(centers)!=K) stop("If K and centers both are provided, then nrow(centers) must equal K!!!")
      if(nrow(centers)==K) K <- NULL
    }
    if(!is.null(centers) && ncol(centers)!=(ncol(x)+ncol(y))) stop("If centers is provided, then ncol(centers) must equal ncol(x).")
    if(is.null(wbounds)) wbounds <- seq(1.1, sqrt(max(group)), len=20)
    if(is.null(vbounds)) vbounds <- (wbounds+1)/2
    if(min(wbounds)<=1) stop("wbounds should be greater than 1")
    wbounds <- c(wbounds) # In case wbounds is a single number, turn it into a vector
    out <- list()
    
    total_data <- cbind(x,y)
    t_v <- rep(1,ncol(x)+ncol(y))
    
    if(!is.null(K)) Cs <- kmeans(total_data, centers=K, nstart=nstart)$cluster
    if(is.null(K)) Cs <- kmeans(total_data, centers=centers)$cluster
    for(i in 1:length(wbounds)){
      if(length(wbounds)>1 && !silent) cat(i,fill=FALSE)
      vs <- rep(1/sqrt(ncol(x)+ncol(y)), (ncol(x)+ncol(y))) # Start with equal weights on each feature
      vs.old <- rnorm(ncol(x)+ncol(y))
      store.bcss.vs <- NULL
      niter <- 0
      while((sum(abs(vs-vs.old))/sum(abs(vs.old)))>1e-4 && niter<maxiter){
        if(!silent) cat(niter, fill=FALSE)
        niter <- niter+1
        vs.old <- vs
        if(!is.null(K)){
          if(niter>1) Cs <- UpdateCs_DBI(total_data, K, vs, Cs, t_v) # if niter=1, no need to update!!
        } else {
          if(niter>1) Cs <- UpdateCs_DBI(total_data, nrow(centers), vs, Cs) # if niter=1, no need to update!!
        }
        #vs <- UpdateWs_DBI(x, Cs, group, category, wbounds[i], vbounds[i])$vs
        ws_dbi <- UpdateWs_SC(x, y, group, Cs, wbounds[i], vbounds[i])
        t_v <- c(rep(1,ncol(x)),ws_dbi$t_v)
        vs <- ws_dbi$ws
        ws_cate <- rep(ws_dbi$ws_cate,table(group))
        ws_cont <- ws_dbi$ws_cont
        vg <- rep(1:length(unique(group[ws_cate>0])),table(group[ws_cate>0]))
        final_dbi <- c(GetcontSC(x[,ws_cont>0], Cs)$dbi_cont_v,rep(GetcateSC(y[,ws_cate>0], Cs, vg)$dbi_cate_v,table(vg)))
        store.bcss.vs <- c(store.bcss.vs, mean(1/final_dbi))
      }
      out[[i]] <- list(vs=vs, Cs=Cs,
                       DBI.perfeaturegroup=GetcateDBI(y, Cs, group)$dbi_cate_v,
                       DBI.perfeature=GetcontDBI(x, Cs)$dbi_cont_v,
                       crit=store.bcss.vs, wbound=wbounds[i])
    }
    if(!silent) cat(fill=TRUE)
    #  if(length(wbounds)==1){
    #    out <- out[[1]]
    #    class(out) <- "kmeanssparse"
    #    return(out)
    #  }
    #  class(out) <- "multikmeanssparse"
    #  return(out)
    class(out) <- "KMeansSparseCluster"
    return(out)
  }



#### Sparse Alternate Sum Clustering Method (SAS) ####

require(sparcl)
require(mclust)
require(phyclust)

withinss = function(X,G,K){
  if(is.vector(X) && is.atomic(X)){
    group = list()
    wcss = 0
    for(j in 1:K){
      group[[j]] = which(G == j)
      l = length(group[[j]])
      if(l>1){
        wcss = wcss + var(X[group[[j]]])*(l-1)
      }
    }
  } else if (is.matrix(X)){
    group = list()
    wcss = 0
    for(j in 1:K){
      group[[j]] = which(G == j)
      l = length(group[[j]])
      if(l>1){
        wcss = wcss + sum(apply(X[group[[j]],],2,var)*(l-1))
      }   
    }
  } else {
    cat("X is niether a vector nor a matrix! \n")
    return(NULL)
  } 
  
  return(wcss)
}

Alternate= function(X, k,tot, initial_set, s, itermax, threshold){
  p = dim(X)[2]
  set0 = initial_set
  set1 = c(0)
  iternum = 0
  while(iternum<= itermax && length(setdiff(set1,set0)) + length(setdiff(set0,set1)) > threshold ){
    clustering = kmeans(X[,set0],iter.max = 20, centers=k,algorithm = "Hartigan-Wong",trace = 0,nstart=2)
    result = clustering$cluster
    
    wcss = apply(X,2,withinss,G = result, K = k)
    iternum = iternum + 1
    
    set1 = set0
    set0 = which(rank((tot-wcss)/tot,ties.method = "random") > p-s)
  }
  out = list(final_set = set0, iternum = iternum, result = result, betweenss = clustering$betweenss)
  return(out)
}

#compute within-cluster distance by clustering feature by feature, select S of size s based on this
hill_climb_GSS = function(X,k,nperms=20,itermax,threshold,tolerance){
  X = scale(X)
  n = dim(X)[1]
  p = dim(X)[2]
  tot = apply(X,2,var)*(n-1)
  permx <- list()
  for(i in 1:nperms){
    permx[[i]] <- matrix(NA, nrow=n, ncol=p)
    for(j in 1:p) permx[[i]][,j] <- sample(X[,j])
  }
  wcss = rep(0,p)

  for(j in 1:p){
    clustering = kmeans(X[,j],iter.max = 10, centers = 2,algorithm = "Hartigan-Wong", trace = 0)
    wcss[j] = clustering$tot.withinss
  }
  rank0 = rank((tot-wcss)/tot,ties.method = "random")
  golden.ratio = 2/(sqrt(5) +1)
  iteration = 0
  upper.bound = p
  lower.bound = 1
  p1 = floor(upper.bound - golden.ratio*(upper.bound-lower.bound))
  p2 = floor(lower.bound + golden.ratio*(upper.bound-lower.bound))
  #evaluate the gap statistics using p1 and p2
  initial_set = which(rank0 > p-p1)
  out1 = Alternate(X, k,tot, initial_set, p1, itermax, threshold)
  permtots = rep(0,nperms)
  for(t in 1:nperms){
    permresult = kmeans(permx[[t]][,out1$final_set], iter.max = 20, centers=k, algorithm = "Hartigan-Wong", trace = 0)
    permtots[t] <- permresult$betweenss
  }
  gap1 = log(out1$betweenss) - mean(log(permtots))
  initial_set = which(rank0 > p-p2)
  out2 = Alternate(X, k,tot, initial_set, p2, itermax, threshold)
  permtots = rep(0,nperms)
  for(t in 1:nperms){
    permresult = kmeans(permx[[t]][,out2$final_set], iter.max = 20, centers=k, algorithm = "Hartigan-Wong", trace = 0)
    permtots[t] <- permresult$betweenss
  }
  gap2 = log(out2$betweenss) - mean(log(permtots))  
  
  while(abs(upper.bound - lower.bound) > tolerance)
  {
    iteration = iteration + 1
    if(gap2 < gap1) # then the maximum is to the left of x2
    {
      upper.bound = p2
      p2 = p1
      gap2 = gap1
      p1 = floor(upper.bound - golden.ratio*(upper.bound - lower.bound))
      #evaluate gaps for p1
      initial_set = which(rank0 > p-p1)
      out1 = Alternate(X, k,tot, initial_set, p1, itermax, threshold)
      
      permtots = rep(0,nperms)
      for(t in 1:nperms){
        permresult = kmeans(permx[[t]][,out1$final_set], iter.max = 20, centers=k, algorithm = "Hartigan-Wong", trace = 0)
        permtots[t] <- permresult$betweenss
      }
      gap1 = log(out1$betweenss) - mean(log(permtots))
    } else {
      # the minimum is to the right of x1
      lower.bound = p1
      p1 = p2
      gap1 = gap2
      p2 = floor(lower.bound + golden.ratio * (upper.bound - lower.bound))
      #evaluate gaps for p2
      initial_set = which(rank0 > p-p2)
      out2 = Alternate(X, k,tot, initial_set, p2, itermax, threshold)
      
      permtots = rep(0,nperms)
      for(t in 1:nperms){
        permresult = kmeans(permx[[t]][,out2$final_set], iter.max = 20, centers=k, algorithm = "Hartigan-Wong", trace = 0)
        permtots[t] <- permresult$betweenss
      }
      gap2 = log(out2$betweenss) - mean(log(permtots))
    }   
  }
  s = floor((lower.bound + upper.bound)/2)
  initial_set = which(rank0 > s)
  out = Alternate(X, k,tot, initial_set, s, 2*itermax, threshold)
  output = list(final_set = out$final_set, iternum = iteration, result = out$result, s = s)
  return(output)
}