#######################################################
## Data preparation for post-hoc analysis
#######################################################
### Author: Noemie Lefrancq
### Date creation: 03/02/2020
### Last modification: 24/10/2021
#######################################################

## Load packages
library(ape)
library(stringr)
library(REdaS)

##############################
## Input trees
##############################
print('loading trees...')
list_trees = list.files(path = 'trees', pattern = "*_g1972tdrmean.prunnedotbks.nwk", all.files = FALSE,
                        full.names = FALSE, recursive = FALSE,
                        ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
trees = lapply(list_trees, function(x) read.tree(paste0('trees/', x)))
print('trees loaded')
ntrees<-length(trees)

## Select trees for bootstrap
nsim = 100 # number of trees to use to bootstrap
index_selection = sample(1:ntrees,nsim)
trees_selected = trees[index_selection]
trees = trees_selected
remove(trees_selected)

##############################
## Extract metadata for each isolate (data.seq)
##############################
master.seq.names = trees[[1]]$tip.label ## Isolates isolates 
data.seq = data.frame('tip_name' = master.seq.names) ## Isolates isolates 
data.seq$Id = unlist(lapply(master.seq.names, function(x)str_split(x, '_')[[1]][1])) ## Extract ID 
data.seq$Source = as.factor(unlist(lapply(master.seq.names, function(x)str_split(x, '_')[[1]][2]))) ## Extract Source
data.seq$Year = as.numeric(unlist(lapply(master.seq.names, function(x)str_split(x, '_')[[1]][3]))) ## Extract Year of collection
data.seq$Continent =  as.factor(unlist(lapply(master.seq.names, function(x)str_split(x, '_')[[1]][4]))) ## Extract continent of collection
data.seq$Country =  as.factor(unlist(lapply(master.seq.names, function(x)str_split(x, '_')[[1]][5]))) ## Extract country of collection
data.seq$Sublineage =  as.factor(unlist(lapply(master.seq.names, function(x)str_split(x, '_')[[1]][6]))) ## Extract isolate's sublineage
data.seq$Region =  as.factor(unlist(lapply(master.seq.names, function(x)str_split(x, '_')[[1]][7]))) ## Extract sub-national region of collection

## Only keep sequences isolated in humans
data.seq = data.seq[which(data.seq$Source == 'H'),]
master.seq.names = data.seq$tip_name

##Add coordinates of each French department
loc_data = read.csv("metadata_and_matrices/coordinates_FRdpts.txt", sep = '\t')
data.seq$DPT_long =  rep(NA, nrow(data.seq))
data.seq$DPT_lat =  rep(NA, nrow(data.seq))
loc_data$code = as.factor(paste0("dpt", loc_data$dpt))
for (i in levels(loc_data$code)){
  data.seq$DPT_long[which(data.seq$Region == i)] = loc_data[which(loc_data$code == i),]$longitude
  data.seq$DPT_lat[which(data.seq$Region == i)] = loc_data[which(loc_data$code == i),]$latitude
}

##Add coordinates of each country
loc_data = read.csv("metadata_and_matrices/countrycoordinates.txt", sep = '\t')
data.seq$Country_long =  rep(NA, nrow(data.seq))
data.seq$Country_lat =  rep(NA, nrow(data.seq))
for (i in levels(data.seq$Country)){
  data.seq$Country_long[which(data.seq$Country == i)] = loc_data[which(loc_data$code == i),]$longitude
  data.seq$Country_lat[which(data.seq$Country == i)] = loc_data[which(loc_data$code == i),]$latitude
}

##Combine both coordinates
data.seq$Combined_long = data.seq$Country_long
data.seq$Combined_lat = data.seq$Country_lat
data.seq$Combined_long[which(is.na(data.seq$DPT_long) == F)] = data.seq$DPT_long[which(is.na(data.seq$DPT_long) == F)]
data.seq$Combined_lat[which(is.na(data.seq$DPT_lat) == F)] = data.seq$DPT_lat[which(is.na(data.seq$DPT_lat) == F)]

##################################################
## Construct matrix: time between infections
##################################################
time_mat = abs(outer(data.seq$Year,data.seq$Year,"-")) ## compute the sampling time difference between all isolates 
# colnames(time_mat) = master.seq.names
# rownames(time_mat) = master.seq.names
diag(time_mat)<-NA ## set the diagonal to NAs.
##################################################
##################################################
## Construct matrix: geography (same/different)
##################################################
## Generates matrix of 0s and 1s
## 0 = same location
## 1 = different location 
## locations considered: departments in France, Paris department in France, countries, continents

# Matrix Department in France
a = length(which(data.seq$Country=='FR'))
data.seq_tmp = data.seq[which(data.seq$Country=='FR'),]
geo_mat_depart_france<-matrix(0,a,a)
for (i in levels(data.seq_tmp$Region)){
  n = which(data.seq_tmp$Region == i)
  geo_mat_depart_france[n,n] = 1
}
n = which(data.seq_tmp$Region == "nd")
geo_mat_depart_france[n,n] = NA
# colnames(geo_mat_depart_france) = master.seq.names[which(data.seq$Country=='FR')]
# rownames(geo_mat_depart_france) = master.seq.names[which(data.seq$Country=='FR')]
diag(geo_mat_depart_france)<-NA

# Matrix Department PARIS in France
a = length(which(data.seq$Country=='FR'))
data.seq_tmp = data.seq[which(data.seq$Country=='FR'),]
geo_mat_depart_Paris_france<-matrix(NA,a,a)
for (i in 'dpt75'){
  n = which(data.seq_tmp$Region == i)
  geo_mat_depart_Paris_france[n,] = 0
  geo_mat_depart_Paris_france[,n] = 0
  geo_mat_depart_Paris_france[n,n] = 1
}
n = which(data.seq_tmp$Region == "nd")
geo_mat_depart_Paris_france[n,n] = NA
# colnames(geo_mat_depart_Paris_france) = master.seq.names[which(data.seq$Country=='FR')]
# rownames(geo_mat_depart_Paris_france) = master.seq.names[which(data.seq$Country=='FR')]
diag(geo_mat_depart_Paris_france)<-NA

#Matrix countries
geo_mat_country<-matrix(0,length(data.seq$Country),length(data.seq$Country))
for (i in levels(data.seq$Country)){
  n = which(data.seq$Country == i)
  geo_mat_country[n,n] = 1
}
# colnames(geo_mat_country) = master.seq.names
# rownames(geo_mat_country) = master.seq.names
diag(geo_mat_country)<-NA

#Matrix continent
geo_mat_continent<-matrix(0,length(data.seq$Continent),length(data.seq$Continent))
for (i in levels(data.seq$Continent)){
  n = which(data.seq$Continent == i)
  geo_mat_continent[n,n] = 1
}
# colnames(geo_mat_continent) = master.seq.names
# rownames(geo_mat_continent) = master.seq.names
diag(geo_mat_continent)<-NA

##################################################
## Construct matrix: distance between isolates
##################################################
#Matrix of distances between isolates, as precise as possible
## each cell corresponds to the distance between isolate i and j
get_dist = function(x1,x2,y1,y2){
  ## Function to compute the distance between 2 points of longitude/latitude coordinates (x1,y1) and (x2,y2)
  R=6371000 ## radius of Earth
  x1 = deg2rad(x1)
  x2 = deg2rad(x2)
  y1 = deg2rad(y1)
  y2 = deg2rad(y2)
  delta.x = x2 - x1
  delta.y = y2 - y1
  a = sin(delta.x/2.0)**2+ cos(x1)*cos(x2)*sin(delta.y/2.0)**2
  c=2*atan2(sqrt(a),sqrt(1-a))
  
  self.meters=R*c 
  km=self.meters/1000.0              # output distance in kilometers
  
  return(km)
}
geo_mat_km<-matrix(0,length(master.seq.names),length(master.seq.names))
a = match(master.seq.names, data.seq$tip_name)
for (i in 1:length(master.seq.names)){
  geo_mat_km[i,] = get_dist(rep(data.seq$Combined_lat[i], length(master.seq.names)), data.seq$Combined_lat, 
                            rep(data.seq$Combined_long[i], length(master.seq.names)), data.seq$Combined_long)
}
diag(geo_mat_km) = NA
# colnames(geo_mat_km) = master.seq.names
# row.names(geo_mat_km) = master.seq.names

#Matrix of distances between isolates, at the country scale: only take the coordinates of the centroid into account
## each cell corresponds to the distance between isolate i and j
geo_mat_km_centroids<-matrix(0,length(master.seq.names),length(master.seq.names))
a = match(master.seq.names, data.seq$tip_name)
for (i in 1:length(master.seq.names)){
  geo_mat_km_centroids[i,] = get_dist(rep(data.seq$Country_lat[i], length(master.seq.names)), data.seq$Country_lat, 
                                      rep(data.seq$Country_long[i], length(master.seq.names)), data.seq$Country_long)
}
diag(geo_mat_km_centroids) = NA
# colnames(geo_mat_km_centroids) = master.seq.names
# row.names(geo_mat_km_centroids) = master.seq.names
geo_mat_km_centroids[which(geo_mat_km_centroids == 0)] = NA ## isolates from the same country are not considered in the matrix: NA
######################################################################

######################################################################
## Construct matrix: genetic distance between isolates, for each tree
######################################################################
########### Store distance matrices into a list of matrices ##########
sim.mats<-array(NaN,c(length(master.seq.names),length(master.seq.names),nsim))

for (ii in 1:nsim){
  dist.mat<-cophenetic.phylo(trees[[ii]]) ## Compute cophenetic distances between all pairs of isolates, for the tree ii
  seq.names<-row.names(dist.mat) ## extract tips names
  a<-sapply(seq.names, function(x)strsplit(x,split="_")[[1]][[1]]) ## extract ID of the tips
  b<-match(as.character(data.seq$Id), a) ## match tree ii IDs with the master dataset
  b<-b[which(is.na(b)==F)] ## remove any tip in tree ii that is not present in the master dataset (eg non-human isolates)
  
  gene_mat = matrix(dist.mat,nrow(dist.mat),nrow(dist.mat))[b,b] ## use the matching vector to reorder the matrix, to have comparable matrices across all trees
  # colnames(gene_mat) = master.seq.names ## rename columns
  # rownames(gene_mat) = master.seq.names ## rename columns
  diag(gene_mat)<-NA ## make sure diagonal is NAs
  
  MRCA_mat = (gene_mat - time_mat)/2 # make the genetic matrix independent of the difference in the isolation dates 
  MRCA_mat[which(MRCA_mat < 0)] = 0 ## sanity check, do not allow distance <0, make them 0
  
  sim.mats[,,ii]<- MRCA_mat ## add this matrix to the main list of matrices 
  
  remove(dist.mat)
  remove(gene_mat)
  remove(MRCA_mat)
  
  print(paste0('Generation of the dist matrix ', ii, '/', nsim))
}
remove(trees)
######################################################################

######################################################################
## Save dataset and matrices
######################################################################
## Metadata - some columns are removed for confidentiality reasons - please email the authors for any inquiry 
data.seq = data.seq[,-c(1,3,8,9,10,11,12,13,14)]
saveRDS(data.seq, 'metadata_and_matrices/data.seq.rds')

## Time
saveRDS(time_mat, 'metadata_and_matrices/time_mat.rds')

## Geography
saveRDS(geo_mat_depart_france, 'metadata_and_matrices/geo_mat_depart_france.rds')
saveRDS(geo_mat_depart_Paris_france, 'metadata_and_matrices/geo_mat_depart_Paris_france.rds')
saveRDS(geo_mat_country, 'metadata_and_matrices/geo_mat_country.rds')
saveRDS(geo_mat_continent, 'metadata_and_matrices/geo_mat_continent.rds')
saveRDS(geo_mat_km, 'metadata_and_matrices/geo_mat_km.rds')
saveRDS(geo_mat_km_centroids, 'metadata_and_matrices/geo_mat_km_centroids.rds')

## Genetic distances 
for(i in 1:nsim){
  saveRDS(sim.mats[,,i], paste0('metadata_and_matrices/sim.mats_',i ,'.rds'))
}
######################################################################
