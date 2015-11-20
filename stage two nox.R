#Stange Two model for NOx 
# Load the necessary packages 
library(RPostgreSQL) 
library(rgdal)   
library(spdep)    
library(maptools)
library(rgeos)
library(BayesXsrc)
library(R2BayesX)
library(BayesX)
library(ggplot2)
# Establish the database connection 
drv=dbDriver("PostgreSQL")
con=dbConnect(drv,dbname="la_birth",user="postgres", password="*",host='localhost', port=9434) 
# Load the spatial polygon frame objects from the table, la_tract_map 
myploylayer=readOGR(dsn="PG:host=localhost user=postgres dbname=la_birth password=*  port=9434", layer="tracts_nox_effects") 
# Get the data of attributes from the spatial polygon objects 
dataSet=attr(myploylayer, "data") 
# Convert geographical objects of class "SpatialPolygons" to objects of class "bnd" from use in BayesX 
# The gid of spatial objects also regarded as the unique region identifier used in spatial modeling 
mybnd=sp2bnd(myploylayer,regionNames=as.character(dataSet$gid)) 
# Filter out those with NULL for the effect;  sp_nox_ap_mean: health effects of NOx output from Stage One 
dset=dataSet[!is.na(dataSet$sp_nox_ap_mean),] 
# Formula for NOx spatial modeling, only the predictors with p-values<0.05 selected according to exploratory test of BayesX  
formulaStr="sp_nox_ap_mean~sx(gid,bs='mrf',map=mybnd) + sx(dist0) + sx(p_per_tra) + sx(p_wbo_tra) +
            sx(p_trtime_30m_sm)+ sx(race_black) + sx(med_h_income) + sx(lu_0_ele) + sx(lu_h_ind) + 
            sx(p_clean_energy)"  

bmodel=bayesx(as.formula(formulaStr),method="MCMC",iterations=30000,burnin=10000, data = dset,chains=9)  
# bayesx_logfile(bmodel)

final_fitted=fitted(bmodel) 
## Total diagnostic test for MCMC convergence  
gdiag_total=GRstats(bmodel)
gdiag_total 

#Set of the predictors for influence of NOx effect across Census tracts 
selcovs=c("dist0","p_per_tra", "p_wbo_tra","p_trtime_30m_sm","race_black",
          "med_h_income","lu_0_ele","lu_h_ind","p_clean_energy")

#Get the summary such as DIC from the fitted models  
aa=(summary(bmodel)[[1]])$model.fit
dic=aa$DIC

#Iterate with every predictor to get their respective non-linear association plots, 
#    the change in effects between the 1st and 3rd quantile  and other statistics 
for( i in c(1:length(selcovs))){ # i=1 
  acov=selcovs[i]
  sxcov=paste("sx(",acov,")",sep="") 
#Extractthe the posterior (fitted) model term partial effects. 
  f=fitted(bmodel[[1]],term=sxcov) 
  covvalues=dset[,acov]
  xlim=range(covvalues,na.rm=TRUE)
  qun_vv=quantile(covvalues)   
  low_mean=mean(final_fitted[[1]][which(covvalues<=qun_vv[2])])  
  high_mean=mean(final_fitted[[1]][which(covvalues>=qun_vv[4])])  
#Obtain the change in effect from the 1st and 3rd quantile of the predictor 
  per=(high_mean-low_mean)/low_mean
#Set the plot file path to save the non-linear association plot between the predictor (acov) and the effect 
  pltfl=paste0("/tmp/",acov,"_effect.png")
#Conduct the Gelman and Rubin's convergence diagnostic of the model parameter, acov   
  gdiag=GRstats(bmodel,term=sxcov)
#Make the non-linear association plot and save it into the plot file (png format)   
  png(filename=pltfl,width=500,height=500) 
  par(mar=c(2,4.3,1,1))    
  ylab = expression (paste("Influence on","  ",NO[x], " effect (g/ppb)"))
  plot2d(f,residuals=FALSE,rug=FALSE,jitter=FALSE,xlim=xlim,xlab="",ylab=ylab)
  dev.off()
#Store the output into the target data frame, covs_eff_statistics 
  sts=data.frame(cov=acov,sample=length(covvalues),lowmean=low_mean,highmean=high_mean,
                 mpsrf=gdiag$mpsrf,coef_chg=high_mean-low_mean,perchange=per)
  if(i==1){
    covs_eff_statistics=sts
  }else{
    covs_eff_statistics=rbind(covs_eff_statistics,sts)  
  }
} 
#Save the effect file into a csv file (nox_inf_effect.csv) 
write.csv(covs_eff_statistics,file="nox_inf_effect.csv",row.names=F)  

#Select the first model for posterior estimates 
abmodel=bmodel[[1]]
#Calculate the variance explained from the complete model selected 
obs_v=dset$sp_nox_ap_mean; poster_v=fitted(abmodel)
var_explained=cov(obs_v,poster_v)^2/(var(obs_v)*var(poster_v))

#Revise the formula with removal of spatial effect 
formulaStr_nosp="sp_nox_ap_mean~ sx(dist0) + sx(p_per_tra) + sx(p_wbo_tra) +
            sx(p_trtime_30m_sm)+ sx(race_black) + sx(med_h_income) + sx(lu_0_ele) + sx(lu_h_ind) + 
            sx(p_clean_energy) "  
#Fitted the model with removal of spatial effect 
bmodel_nosp=bayesx(as.formula(formulaStr_nosp),method="MCMC",iterations=30000,burnin=10000, data = dset,chains=9)  
# bayesx_logfile(bmodel)
#Calculate the variance explained from the model with the removal of spatial effect  
abmodel_nosp=bmodel_nosp[[1]]
obs_v=dset$sp_nox_ap_mean; poster_v=fitted(abmodel_nosp)
var_explained_nosp=cov(obs_v,poster_v)^2/(var(obs_v)*var(poster_v))

#Obtain the variance explained by spatial effects 
var_explained_bysp=var_explained-var_explained_nosp  

#Obtain the posterior estimation 
poster_estimates=data.frame(tract=dset$tract,observed=dset$sp_nox_ap_mean,post_est=final_fitted[[1]]) 

#Obtain the posterior probabilities of effects below zero
#In BayesX, we used classical Gaussian regression model specified with homoscedastic variances, i.e. the responses
#yi are Gaussian with expected value μ = ηi depending on covariates and variance σ2 not depending on covariates.
poster_estimates$std=sqrt(abmodel$variance[1,"Mean"])  
poster_estimates$prom0=pnorm(0,mean=poster_estimates$post_est,sd=poster_estimates$std) 

#Output the posterior file for mapping in ArcGIS 
write.csv(poster_estimates,file="poster_estimates_nox.csv",row.names=F) 


 

