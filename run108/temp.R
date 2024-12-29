versionid = 93
setwd(paste0("/n/home07/rbhaduri/projects/Markov blanket for manifold supported covariates/run", versionid))
result = readRDS("result files/Result-model=softmax_log_normal,n=100,m=100,m_1=10,alpha=0.1,CRT.replications=1500,simulation.iter=1,snr=0.01,p_filter_tolerance=0.9,seed=1")
result[[1]][[1]]$d_CRT_with_r_squared$runtime

length(c(0.01,0.1,0.25,0.5,0.75, 1:15))
length(c(0.01,0.1,0.25,0.5,0.75,1,3,5,7,9,11,13,15))

20*250

13*3000 / 6
