#######################################################
## Supplementary FigureS13 : Relative risk within departments, in France
#######################################################
### Author: Noemie Lefrancq
### Date creation: 03/02/2020
### Last modification: 17/10/2021
#######################################################

## Main function
ratio_bootstrap_dist_discrete_auto <- function(x, y, geo_mat, time_mat2, MRCA_mat2, geo_mat_ref, time_mat2_ref, MRCA_mat2_ref){
  ## This function computes the relative risk, for a given 
  ##     -location matrix (geo_mat), 
  ##     -time matrix (time_mat2), 
  ##     -genetic matrix (MRCA_mat2), 
  ## compared to a reference:
  ##     -location matrix (geo_mat_ref), 
  ##     -time matrix (time_mat2_ref), 
  ##     -genetic matrix (MRCA_mat2_ref),
  ## for given bootstrap vectors x and y.
  
  geo_mat.tmp = geo_mat[x,x]
  time_mat.tmp2 = time_mat2[x,x]
  MRCA_mat.tmp2 = MRCA_mat2[x,x]
  geo_mat_ref.tmp = geo_mat_ref[y,y]
  time_mat2_ref.tmp2 = time_mat2_ref[y,y]
  MRCA_mat2_ref.tmp2 = MRCA_mat2_ref[y,y]
  
  tmp = MRCA_mat.tmp2 * time_mat.tmp2
  tmp[which(tmp == 0)] = NA
  
  tmp2 = time_mat.tmp2
  tmp2[which(tmp2 == 0)] = NA
  
  tmp_ref = MRCA_mat2_ref.tmp2 * time_mat2_ref.tmp2
  tmp_ref[which(tmp_ref == 0)] = NA
  
  tmp2_ref = time_mat2_ref.tmp2
  tmp2_ref[which(tmp2_ref == 0)] = NA
  
  a1 = cumsum(hist(tmp*geo_mat.tmp, breaks = c(0,Pmax,1E10), plot = F)$counts)
  a2 = cumsum(hist(tmp*geo_mat.tmp, breaks = c(0,Pmin,1E10), plot = F)$counts)
  a = a1 - a2
  
  b1 = cumsum(hist(tmp_ref * geo_mat_ref.tmp, breaks = c(0,Pmax,1E10), plot = F)$counts)
  b2 = cumsum(hist(tmp_ref * geo_mat_ref.tmp, breaks = c(0,Pmin,1E10), plot = F)$counts)
  b = b1 - b2
  
  c = sum(tmp2*geo_mat.tmp,na.rm=T)
  
  d = sum(tmp2_ref*geo_mat_ref.tmp,na.rm=T)
  
  rr.out = (a/c)/(b/d) 
  
  return(rr.out[-length(rr.out)])
}

#####################################################################
## Load relevant datasets and matrices
#####################################################################
## Metadata
data.seq = readRDS('metadata_and_matrices/data.seq.rds')
## Time
time_mat = readRDS('metadata_and_matrices/time_mat.rds')
## Geography
geo_mat_country = readRDS('metadata_and_matrices/geo_mat_country.rds')
geo_mat_depart_france = readRDS('metadata_and_matrices/geo_mat_depart_france.rds')
## Genetic distances 
sim.mats = readRDS('metadata_and_matrices/sim.mats.rds')

#####################################################################
## Parameters for the computation of the relative risks
#####################################################################
## Number of trees to consider (same as the length of sim.mats)
nsim = 100
## Number of bootstrap event to perform, of each tree
nboot = 10

## MRCA windows on which to compute the relative risk
Pmax = c(5,20,50,1000) ## max windows
Pmin = c(0,5,20,50) ## min windows
pmid<-(Pmin+Pmax)/2 ## mid-point
int = c(0,5,20,50,1000)
l = length(Pmax) ## number of intervals

## Set the boot matrix to save the results
rr = matrix(NA,l, nboot*nsim)

#####################################################################
## Compute relative risk
#####################################################################
for(ii in 1:nsim){
  print(ii)
  ## Choose isolates
  a = c(which(data.seq$Country == 'FR' & data.seq$Year>=2010 & data.seq$Sublineage!= 'SL150' & data.seq$Sublineage!= 'SL404'))
  b = which(data.seq[which(data.seq$Country == 'FR'),]$Year>=2010 & data.seq[which(data.seq$Country == 'FR'),]$Sublineage!= 'SL150' &
              data.seq[which(data.seq$Country == 'FR'),]$Sublineage!= 'SL404')
  ## Geography matrix: same department in France
  geo_mat = geo_mat_depart_france[b,b]
  time_mat2 = time_mat[a,a]<=2
  MRCA_mat = sim.mats[,,ii]
  MRCA_mat2 = MRCA_mat[a,a]
  nseq = length(b)
  
  ## Reference: Betwwen departments
  ref = c(which(data.seq$Year>=2010 & data.seq$Sublineage!= 'SL150' & data.seq$Sublineage!= 'SL404'))
  geo_mat_ref = 1-geo_mat_depart_france[b,b]
  geo_mat_ref[which(geo_mat_ref == 0)] = NA
  time_mat2_ref = time_mat2
  MRCA_mat2_ref = MRCA_mat2
  nseq_ref = nseq
  
  ##Bootstrap to create the ci
  for (i in (1:nboot)){
    tmp = sample(nseq, nseq, replace = T)
    tmp_ref = sample(nseq_ref, nseq_ref, replace = T)
    rr.out = ratio_bootstrap_dist_discrete_auto(tmp, tmp_ref, geo_mat, time_mat2, MRCA_mat2, geo_mat_ref, time_mat2_ref, MRCA_mat2_ref)
    rr[,((ii-1)*nboot + i)] = rr.out
  }
}
#####################################################################

#####################################################################
## Write results
#####################################################################
res = list('rr' = rr,
           'int' = int,
           'nboot' = nboot,
           'nsim' = nsim)
## Write results
saveRDS(res, 'results/Relative_risk_within_between_depts_France.rds')
#####################################################################


