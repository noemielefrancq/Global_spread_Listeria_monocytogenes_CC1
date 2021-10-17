## Plots all figures
#######################################################
### Author: Noemie Lefrancq
### Date creation: 03/02/2020
### Last modification: 17/10/2021
#######################################################

################################################################################################
## Figure 4a
################################################################################################
res = readRDS('results/Relative_risk_across_scales_and_MRCA.rds')
res = readRDS('E:/Stages/Project_listeria/Multinomial_analysis/Figures/results_listeria_rr_19052020.rds')
int = res$int
nboot = res$nboot
nsim = res$nsim

## Plot 
rr[which(is.infinite(rr))] = NA
grDevices::windows(height=9.5, width =3)
par(mfrow = c(l,1), mar = c(2, 4, 1, 0.5), oma = c(7, 2, 0, 0))
n_steps = 4
for (i in (1:l)){
  boot.ci = apply(rr[c(0:(n_steps-1)*(l)+i),(1:(nboot*nsim))], 1, quantile, probs = c(0.025,0.975), na.rm = T)
  boot.ci.m2 = apply(rr[c(0:(n_steps-1)*(l)+i),(1:(nboot*nsim))], 1, quantile, probs = c(0.5), na.rm = T)
  print(boot.ci)
  print(boot.ci.m2)
  
  boot.ci[which(boot.ci > 1E7)] = 1E7
  boot.ci[which(boot.ci < 0.01)] = 0.01
  boot.ci.m2[which(boot.ci.m2 < 0.01)] = 0.01
  
  plot(1:(n_steps), boot.ci.m2, type="p", pch=16,cex=1.4, xlab="", xlim = c(0.5,n_steps+0.5),
       ylim = c(0.01,500), log='y',
       ylab="", main = paste0(int[i], " <MRCA< ", int[i+1], " years"), xaxt="n", cex.axis= 1.1, yaxt="n")
  axis(2, at = c(0.01, 0.1, 1, 10, 100, 1000), label =  c('<0.01', 0.1, 1, 10, 100 ,1000), las = 2, cex.axis= 1.1)
  if(i == l){
    axis(1, at = 1:(n_steps), labels = c("Wihtin countries",
                                         "Between countries \n <1000km",
                                         "Between countries \n >1000km (ref)",
                                         "Between continents"), cex.axis=1, las =2)
  }
  else{
    axis(1, at = 1:(n_steps), labels = c("", "", "", ""), cex.axis=1, las =2)
  }
  abline(h = 1, col = 'red')
  arrows(1:(n_steps), boot.ci[1,], 1:(n_steps), boot.ci[2,], length=0.05, angle=90, code=3)
}
mtext(paste0("Relative risk"), outer = T, cex = 1.3, side = 2)
################################################################################################

################################################################################################
## Figure 4b
################################################################################################
boot_out_total = readRDS('results/Proportion_MRCA5y_FRdepts.rds')
pmid = boot_out_total$pmid

## Compute mean and confidence intervals, within French departments
m_same_dept = quantile(boot_out_total$'within_FR_depts', probs = 0.5, na.rm =T)
ci_same_dept = quantile(boot_out_total$'within_FR_depts', probs = c(0.025, 0.975), na.rm =T)

## Compute mean and confidence intervals, mean proportion
m_mean = quantile(boot_out_total$'meanFR', probs = 0.5, na.rm =T)
ci_mean = quantile(boot_out_total$'meanFR', probs = c(0.025, 0.975), na.rm =T)

## Compute mean and confidence intervals, between French departments >500km apart
m_diff_dept_500_fr = quantile(boot_out_total$'prop>500km', probs = 0.5, na.rm =T)
ci_diff_dept_500_fr = quantile(boot_out_total$'prop>500km', probs = c(0.025, 0.975), na.rm =T)

## Compute mean and confidence intervals, between French departments xxkm apart
boot.ci_diff_dept_france = apply(boot_out_total$'prop_0_500km', 1, quantile, probs = c(0.025,0.975), na.rm = T)
boot.ci.m2_diff_dept_france = apply(boot_out_total$'prop_0_500km', 1, quantile, probs = c(0.5), na.rm = T)

## Plot
grDevices::windows(width = 6, height = 6)
plot(NULL, ylim = c(0, 0.15), xlim=c(-10, 550), xaxt = 'n',
     xlab = "Distance in km", yaxt = 'n', ylab = '',
     main = '0<MRCA<5y', pch = 16, lwd = 1.5)
## line for the mean
abline(h=m_mean, col = 'forestgreen', lwd = 2)
polygon(c(seq(-30, 580, 1),rev(seq(-30, 580, 1))), 
        c(rep(ci_mean[1], length(seq(-30, 580, 1))), rev(rep(ci_mean[2], length(seq(-30, 580, 1))))), col=rgb(0.13,0.55,0.13,alpha=0.2), border=NA)
## Continuous distance
points(pmid, boot.ci.m2_diff_dept_france, pch = 16, lwd = 1.5)
arrows(pmid,boot.ci_diff_dept_france[1,],pmid, boot.ci_diff_dept_france[2,], length=0.05, angle=90, code=3, lwd = 1.5)
# same departments france
points(0, m_same_dept, pch = 16)
arrows(0, ci_same_dept[1],0, ci_same_dept[2], length=0.05, angle=90, code=3, lwd = 1.5)
# different department >500 in France
points(550, m_diff_dept_500_fr, pch = 16)
arrows(550,ci_diff_dept_500_fr[1],550, ci_diff_dept_500_fr[2], length=0.05, angle=90, code=3, lwd = 1.5)
# axis
axis(side = 2, las = 2, at = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4,1, 10, 100), 
     labels = c(0, 0.05, 0.1,  0.15, 0.2,0.25 ,0.3, 0.4, 1, 10, 100))
axis(side = 1, las = 1, at = c(0, seq(50, 450, 50), 550), labels = c('same dpt',  seq(50, 450, 50), '>500'))
################################################################################################

################################################################################################
## Figure 4c
################################################################################################
boot_out = readRDS('results/Relative_risk_Paris_other_depts.rds')
nsim = boot_out$nsim
nboot = boot_out$nboot

## Compute mean and confidence intervals - Paris
m_same_dept_paris = quantile(boot_out$within_Paris, probs = 0.5, na.rm =T)
ci_same_dept_paris = quantile(boot_out$within_Paris, probs = c(0.05, 0.95), na.rm =T)

## Compute mean and confidence intervals - other depts
m_same_dept_others = quantile(boot_out$within_other_depts, probs = 0.5, na.rm =T)
ci_same_dept_others = quantile(boot_out$within_other_depts, probs = c(0.05, 0.95), na.rm =T)

## Plot
grDevices::windows(width = 7, height = 7)
par(mar = c(2,5,2,2))
plot(NULL, ylim = c(0.01, 150), xlim=c(-1, 5), xaxt = 'n',
     xlab = "", yaxt = 'n', ylab = 'Relative risk', log = 'y',
     main = 'Relative risk of having 0<MRCA<5y', pch = 16, lwd = 1.5)
abline(h = 1, col = 'red', lwd = 2)
axis(side = 2, las = 2, at = c(0.01, 0.1, 1, 10, 100), labels = c("<0.01", "0.1", "1", "10", "100"), cex.axis = 1)
axis(side = 1, las = 1, at = c(0, 1, 3, 4), labels = c('Within \n Paris', 'Ref: \n Paris -\n other depts', 'Within \n other \n depts', 'Ref: \n Between \n all depts'), cex.axis = 1)

# # dot same departments France
points(0, m_same_dept_paris, pch = 16)
arrows(0, ci_same_dept_paris[1],0, ci_same_dept_paris[2], length=0.05, angle=90, code=3, lwd = 1.5)

# dot reference, different department >500 in France
points(1, 1, pch = 15)

# # dot same departments France
points(3, m_same_dept_others, pch = 16)
arrows(3, ci_same_dept_others[1],3, ci_same_dept_others[2], length=0.05, angle=90, code=3, lwd = 1.5)

# dot reference, different department >500 in France
points(4, 1, pch = 15)
################################################################################################

################################################################################################
## Supplementary Figure S13 : Relative risk within departments, in France
################################################################################################
res = readRDS('results/Relative_risk_within_between_depts_France.rds')
rr = res$rr
rr[which(rr == 0)] = NA ## Remove 0s, as they are likely due to lack of data, for some bootstrap iterations
nboot = res$nboot
nsim = res$nsim

## Compute confidence intervals
boot.ci = apply(rr[,(1:(nboot*nsim))], 1, quantile, probs = c(0.025,0.975), na.rm = T)
## Compute mean
boot.ci.m = apply(rr[,(1:(nboot*nsim))], 1, quantile, probs = c(0.5), na.rm = T)
## Print estimates
print(boot.ci)
print(boot.ci.m)

## To plot on log scale, make sure all values are >0
boot.ci[which(boot.ci < 0.01)] = 0.01
boot.ci.m[which(boot.ci.m < 0.01)] = 0.01

## Plot
grDevices::windows(width = 5, height = 5)
par(mar = c(4,4,4,1), oma = c(0,0,0,0))
plot(1:4, boot.ci.m, type="p", pch=16,cex=1.4, xlab="MRCA window (years)", xlim = c(0.5,4+0.5),
     ylim = c(0.1,100), log='y', 
     main = "Relative risk of coming from the same department, \n when having a MRCA within a specific window of time",
     ylab = "Relative risk", xaxt="n", cex.axis= 1.1, yaxt="n", cex.main = 0.9)
axis(2, at = c(0.01, 0.1, 1, 10, 100, 1000), label =  c('<0.01', 0.1, 1, 10, 100 ,1000), las = 2, cex.axis= 1.1) 
axis(1, at = 1:4, labels = c("0-5", "5-20", "20-50", "50+"), cex.axis=1, las =1)
abline(h = 1, col = 'red')
arrows(1:4, boot.ci[1,], 1:4, boot.ci[2,], length=0.05, angle=90, code=3)
################################################################################################ 
