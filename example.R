source('/functions.R')



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




#### Simulation for K = 4 and p-generalized normal-polynomial distribution.
n1 <- 100
n2 <- 100
p <- 100


CER <- NULL
RI <- NULL
TPR_conti <- NULL
TNR_conti <- NULL
TPR_cate <- NULL
TNR_cate <- NULL
total_var <- NULL



for (i in 1:50) {
  startdata <- generatesimd(n1 = n1, n2 = n2, p = p)
  simd <- startdata$simd
  simgroup <- startdata$simgroup
  x <- simd[,1:p]
  y <- simd[,-(1:p)]
  simgroup <- startdata$simgroup[-(1:p)]-p
  simresult <- COC(x,y,K = 4,group = simgroup,
                   wbounds = c(2.5),vbounds = 2.5)
  CER <- c(CER, cal_CER(simresult[[1]]$Cs, c(rep(1,n1),rep(2,n2),rep(3,n2),rep(4,n2))))
  RI <- c(RI, randIndex(simresult[[1]]$Cs, c(rep(1,n1),rep(2,n2),rep(3,n2),rep(4,n2))))
  TPR_conti <- c(TPR_conti, sum(simresult[[1]]$vs[1:10]>0)/10)
  TNR_conti <- c(TNR_conti, 1-sum(simresult[[1]]$vs[11:p]>0)/(p-10))
  TPR_cate <- c(TPR_cate, sum(simresult[[1]]$vs[(p+1):(p+28)]>0)/28)
  TNR_cate <- c(TNR_cate, 1-sum(simresult[[1]]$vs[(p+29):(p+ncol(y))]>0)/(ncol(y)-14))
  total_var <- c(total_var, sum(simresult[[1]]$vs>0))
}





