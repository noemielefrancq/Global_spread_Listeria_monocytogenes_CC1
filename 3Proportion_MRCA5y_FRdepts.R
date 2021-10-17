#######################################################
## Figure 4b: Proportion of MRCA <5y, across distance in France
#######################################################
### Author: Noemie Lefrancq
### Date creation: 03/02/2020
### Last modification: 17/10/2021
#######################################################

#####################################################################
## Load relevant datasets and matrices
#####################################################################
## Metadata
data.seq = readRDS('metadata_and_matrices/data.seq.rds')
## Time
time_mat = readRDS('metadata_and_matrices/time_mat.rds')
## Geography
geo_mat_depart_france = readRDS('metadata_and_matrices/geo_mat_depart_france.rds')
geo_mat_km = readRDS('metadata_and_matrices/geo_mat_km.rds')
## Genetic distances 
sim.mats = readRDS('metadata_and_matrices/sim.mats.rds')

#####################################################################
## Parameters to compute the proportions
#####################################################################
## Number of trees to consider (same as the length of sim.mats)
nsim = 100
## Number of bootstrap event to perform, of each tree
nboot = 1 ## Number of bootstrap iterations

## MRCA to consider 
MRCA = 5

## List to store the results
boot_out_total = NULL

############################################################
## Proportion within department
############################################################
boot_out_total$'within_FR_depts' = rep(NA, nboot*nsim)
for(ii in 1:nsim){
  print(paste0('nsim: ', ii, '/', nsim))
  for (i in 1:nboot){
    # data
    a = which(data.seq$Country == "FR" & data.seq$Year>=2010)
    b = which(data.seq[which(data.seq$Country == "FR"),]$Year>=2010)
    geo_mat = geo_mat_depart_france[b,b]
    time_mat2 = time_mat[a,a]<=1
    MRCA_mat2 = sim.mats[,,ii][a,a]
    nseq = length(a)
    
    # Bootstrap
    tmp = sample(nseq, nseq, replace = T)
    geo_mat = geo_mat[tmp,tmp]
    time_mat2 = time_mat2[tmp,tmp]
    MRCA_mat2 = MRCA_mat2[tmp,tmp]
    
    # Compute proportion
    a = sum((MRCA_mat2<=MRCA) * time_mat2 *geo_mat, na.rm = T)
    b = sum(time_mat2 *geo_mat, na.rm = T)
    
    rr.out = (a/b)
    
    boot_out_total$'within_FR_depts'[(ii-1)*nboot + i] = rr.out
  }
}
m_same_dept = quantile(boot_out_total$'within_FR_depts', probs = 0.5, na.rm =T)
ci_same_dept = quantile(boot_out_total$'within_FR_depts', probs = c(0.025, 0.975), na.rm =T)

############################################################
## Mean proportion, across all pairs of sequences in France
############################################################
boot_out_total$'meanFR' = rep(NA, nboot*nsim)
for(ii in 1:nsim){
  print(paste0('nsim: ', ii, '/', nsim))
  for (i in 1:nboot){
    # data
    a = which(data.seq$Country == "FR" & data.seq$Year>=2010)
    b = which(data.seq[which(data.seq$Country == "FR"),]$Year>=2010)
    time_mat2 = time_mat[a,a]<=1
    MRCA_mat2 = sim.mats[,,ii][a,a]
    nseq = length(a)
    
    # Bootstrap
    tmp = sample(nseq, nseq, replace = T)
    time_mat2 = time_mat2[tmp,tmp]
    MRCA_mat2 = MRCA_mat2[tmp,tmp]
    
    # Compute proportion
    a = sum((MRCA_mat2<=MRCA) * time_mat2, na.rm = T)
    b = sum(time_mat2, na.rm = T)
    
    rr.out = (a/b)
    
    boot_out_total$'meanFR'[(ii-1)*nboot + i] = rr.out
  }
}
m_mean = quantile(boot_out_total$'meanFR', probs = 0.5, na.rm =T)
ci_mean = quantile(boot_out_total$'meanFR', probs = c(0.025, 0.975), na.rm =T)

############################################################
## Proportion between french departments separated by >500 km
############################################################
boot_out_total$'prop>500km' = rep(NA, nboot*nsim)
for(ii in 1:nsim){
  print(paste0('nsim: ', ii, '/', nsim))
  for (i in 1:nboot){
    # data
    a = which(data.seq$Country == "FR" & data.seq$Year>=2010)
    b = which(data.seq[which(data.seq$Country == "FR"),]$Year>=2010)
    geo_mat = (1-geo_mat_depart_france[b,b])*(geo_mat_km[a,a]>=500)
    time_mat2 = time_mat[a,a]<=1
    MRCA_mat2 = sim.mats[,,ii][a,a]
    nseq = length(a)
    
    # Bootstrap
    tmp = sample(nseq, nseq, replace = T)
    geo_mat = geo_mat[tmp,tmp]
    time_mat2 = time_mat2[tmp,tmp]
    MRCA_mat2 = MRCA_mat2[tmp,tmp]
    
    # Compute proportion
    a = sum((MRCA_mat2<=MRCA) * time_mat2 *geo_mat, na.rm = T)
    b = sum(time_mat2 *geo_mat, na.rm = T)
    
    rr.out = (a/b)
    
    boot_out_total$'prop>500km'[(ii-1)*nboot + i] = rr.out
  }
}
m_diff_dept_500_fr = quantile(boot_out_total$'prop>500km', probs = 0.5, na.rm =T)
ci_diff_dept_500_fr = quantile(boot_out_total$'prop>500km', probs = c(0.025, 0.975), na.rm =T)

############################################################
## Proportion between french departments separated by xx km
############################################################
proportion_bootstrap_dist <- function(x, y, geo_mat, time_mat2, MRCA_mat2){
  geo_mat.tmp = geo_mat[x,x]
  time_mat.tmp2 = time_mat2[x,x]
  MRCA_mat.tmp2 = MRCA_mat2[x,x]
  
  tmp = geo_mat.tmp * time_mat.tmp2
  tmp[which(tmp == 0)] = NA
  
  a1 = cumsum(hist((MRCA_mat.tmp2<=MRCA) * tmp, breaks = c(0,Pmax,1E10), plot = F)$counts)
  a2 = cumsum(hist((MRCA_mat.tmp2<=MRCA) * tmp, breaks = c(0,Pmin,1E10), plot = F)$counts)
  a = a1 - a2
  
  b1 = cumsum(hist(tmp, breaks = c(0,Pmax,1E10), plot = F)$counts)
  b2 = cumsum(hist(tmp, breaks = c(0,Pmin,1E10), plot = F)$counts)
  b = b1 - b2
  
  rr.out = (a/b)
  
  return(rr.out[-length(rr.out)])
}

## Windows of distances to use
Pmax = seq(50, 500, 25) ## max distance
windows = 50 ## window
Pmin = Pmax - windows ## min distance
Pmin[which(Pmin<0)]<-0 
pmid<-(Pmin+Pmax)/2 ## mid-point (to plot)

## Data to use
a = which(data.seq$Country == "FR" & data.seq$Year>=2010)
b = which(data.seq[which(data.seq$Country == "FR"),]$Year>=2010)
boot_out_total$'prop_0_500km' = matrix(NA, length(pmid), nboot*nsim)

for(ii in 1:nsim){
  print(paste0('nsim: ', ii, '/', nsim))
  for (i in (1:nboot)){
    nseq = length(a)
    geo_mat = (1-geo_mat_depart_france[b,b])*geo_mat_km[a,a]
    time_mat2 = time_mat[a,a]<=1
    MRCA_mat2 = sim.mats[,,ii][a,a]
    nseq = length(a)
    
    ##Bootstrap
    tmp = sample(nseq, nseq, replace = T)
    rr.out= proportion_bootstrap_dist(tmp, tmp_ref, geo_mat, time_mat2, MRCA_mat2)
    boot_out_total$'prop_0_500km'[,(ii-1)*nboot + i] = rr.out
  }
}
boot.ci_diff_dept_france = apply(boot_out_total$'prop_0_500km', 1, quantile, probs = c(0.025,0.975), na.rm = T)
boot.ci.m2_diff_dept_france = apply(boot_out_total$'prop_0_500km', 1, quantile, probs = c(0.5), na.rm = T)
############################################################

#####################################################################
## Write results
#####################################################################
boot_out_total$pmid = pmid
boot_out_total$MRCA = MRCA

## Write results
saveRDS(boot_out_total, 'results/Proportion_MRCA5y_FRdepts.rds')
#####################################################################
