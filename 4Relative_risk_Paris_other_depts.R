#######################################################
## Figure 4C: Relative risk wothin Paris / oeght departments
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
data.seq$dept = as.numeric(unlist(lapply(as.character(data.seq$Region), function(x)str_split(x, pattern = 'dpt')[[1]][2]))) ## extract department ID
## French department IDs
depts = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", 
          "16", "17", "18", "19", "21", "22", "23", "24", "25", "26", "27", "28", "29", "2A", "2B",
          "30", "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", "41", "42", "43", "44",
          "45", "46", "47", "48", "49", "50", "51", "52", "53", "54", "55", "56", "57", "58", "59",
          "60", "61", "62", "63", "64", "65", "66", "67", "68", "69", "70", "71", "72", "73", "74",
          "75", "76", "77", "78", "79", "80", "81", "82", "83", "84", "85", "86", "87", "88", "89",
          "90", "91", "92", "93", "94", "95")
dept_paris = c('75')
## Time
time_mat = readRDS('metadata_and_matrices/time_mat.rds')
## Geography
geo_mat_depart_france = readRDS('metadata_and_matrices/geo_mat_depart_france.rds')
geo_mat_depart_Paris_france = readRDS('metadata_and_matrices/geo_mat_depart_Paris_france.rds')
geo_mat_km = readRDS('metadata_and_matrices/geo_mat_km.rds')
## Genetic distances 
sim.mats = readRDS('metadata_and_matrices/sim.mats.rds')

#####################################################################
## Parameters for the computation of the relative risks
#####################################################################
## Number of trees to consider (same as the length of sim.mats)
nsim = 100
## Number of bootstrap event to perform, of each tree
nboot = 10

## MRCA to consider 
MRCA = 5

## Vector designing which dept is parisian
which_depts_is_paris = rep(NA, length(depts))
for (i in 1:length(depts)){
  if (!is.na(match(depts[i], dept_paris))) which_depts_is_paris[i] = 'paris'
  if (is.na(match(depts[i], dept_paris))) which_depts_is_paris[i] = 'others'
}

## List to save the results
boot_out = NULL
#####################################################################

#####################################################################
## Compute the relative risks
#####################################################################
##### Within Paris
boot_out$'within_Paris' = rep(NA, nboot*nsim)

for(ii in 1:nsim){
  print(paste0('nsim: ', ii, '/', nsim))
  for (i in 1:nboot){
    # data
    c = depts[which(which_depts_is_paris == 'paris')]
    d = match(data.seq$dept, c)
    e = match(data.seq[which(data.seq$Country == "FR"),]$dept, c)
    a = which(data.seq$Country == "FR" & data.seq$Year>=2000 & is.na(d) == F)
    b = which(data.seq[which(data.seq$Country == "FR"),]$Year>=2000 & is.na(e) == F)
    
    geo_mat = geo_mat_depart_Paris_france[b,b]
    time_mat2 = time_mat[a,a]<=2
    MRCA_mat2 = sim.mats[,,ii][a,a]
    nseq = length(a)
    
    # Bootstrap
    tmp = sample(nseq, nseq, replace = T)
    geo_mat = geo_mat[tmp,tmp]
    time_mat2 = time_mat2[tmp,tmp]
    MRCA_mat2 = MRCA_mat2[tmp,tmp]
    
    # reference
    ref = which(data.seq$Country == "FR" & data.seq$Year>=2000)
    b = which(data.seq[which(data.seq$Country == "FR"),]$Year>=2000)
    geo_mat_ref = (1-geo_mat_depart_Paris_france[b,b])
    time_mat2_ref = time_mat[ref,ref]<=2
    MRCA_mat2_ref = sim.mats[,,ii][ref,ref]
    
    # bootstrap reference
    nref = length(ref)
    tmp = sample(nref, nref, replace = T)
    geo_mat_ref = geo_mat_ref[tmp, tmp]
    time_mat2_ref = time_mat2_ref[tmp,tmp]
    MRCA_mat2_ref = MRCA_mat2_ref[tmp,tmp]
    
    # Compute rr
    a = sum((MRCA_mat2<=MRCA) * time_mat2 *geo_mat, na.rm = T)
    b = sum(time_mat2 *geo_mat, na.rm = T)
    
    c = sum((MRCA_mat2_ref<=MRCA) * time_mat2_ref *(geo_mat_ref), na.rm = T)
    d = sum(time_mat2_ref *(geo_mat_ref), na.rm = T)
    
    rr.out = (a/b)/(c/d) 
    
    if  (rr.out == 0 | is.na(rr.out)){
      rr.out = ((a+1/b+1))/((c+1/d+1)) 
    }
    
    boot_out$'within_Paris'[(ii-1)*nboot + i] = rr.out
  }
}

##### Within other departments
boot_out$'within_other_depts' = rep(NA, nboot*nsim)

for(ii in 1:nsim){
  print(paste0('nsim: ', ii, '/', nsim))
  for (i in 1:nboot){
    # data
    c = depts[which(which_depts_is_paris != 'paris')]
    d = match(data.seq$dept, c)
    e = match(data.seq[which(data.seq$Country == "FR"),]$dept, c)
    a = which(data.seq$Country == "FR" & data.seq$Year>=2000 & is.na(d) == F)
    b = which(data.seq[which(data.seq$Country == "FR"),]$Year>=2000 & is.na(e) == F)
    
    geo_mat = geo_mat_depart_france[b,b]
    time_mat2 = time_mat[a,a]<=2
    MRCA_mat2 = sim.mats[,,ii][a,a]
    nseq = length(a)
    
    # Bootstrap
    tmp = sample(nseq, nseq, replace = T)
    geo_mat = geo_mat[tmp,tmp]
    time_mat2 = time_mat2[tmp,tmp]
    MRCA_mat2 = MRCA_mat2[tmp,tmp]
    
    # reference
    ref = which(data.seq$Country == "FR" & data.seq$Year>=2010)
    b = which(data.seq[which(data.seq$Country == "FR"),]$Year>=2010)
    geo_mat_ref = (1-geo_mat_depart_france[b,b])
    time_mat2_ref = time_mat[ref,ref]<=2
    MRCA_mat2_ref = sim.mats[,,ii][ref,ref]
    
    # bootstrap reference
    nref = length(ref)
    tmp = sample(nref, nref, replace = T)
    geo_mat_ref = geo_mat_ref[tmp, tmp]
    time_mat2_ref = time_mat2_ref[tmp,tmp]
    MRCA_mat2_ref = MRCA_mat2_ref[tmp,tmp]
    
    # Compute relative risk
    a = sum((MRCA_mat2<=MRCA) * time_mat2 *geo_mat, na.rm = T)
    b = sum(time_mat2 *geo_mat, na.rm = T)
    
    c = sum((MRCA_mat2_ref<=MRCA) * time_mat2_ref *(geo_mat_ref), na.rm = T)
    d = sum(time_mat2_ref *(geo_mat_ref), na.rm = T)
    
    rr.out = (a/b)/(c/d) 
    
    boot_out$'within_other_depts'[(ii-1)*nboot + i] = rr.out
  }
}

#####################################################################
## Write results
#####################################################################
boot_out$'nboot' = nboot
boot_out$'nsim' = nsim

## Write results
saveRDS(boot_out, 'results/Relative_risk_Paris_other_depts.rds')
#####################################################################
