#Stange One model for NO2 
# Load the necessary libraries:  
library(RPostgreSQL)  # the library for database link   
library(mgcv)        # the library of generalized additive models 
library(rjags)         # the rjags library for link to JAGS (required for installation) 
library(R2jags)   

# Define the driver and connector of PostgreSQL with PostGIS 
drv=dbDriver("PostgreSQL") 
con=dbConnect(drv,dbname="la_birth",user="postgres", password="*", host='localhost', port=9434)
# Load the census tracts table of LA from the database into the variable, tracts 
sql=paste("select * from la_tracts ")
tracts = dbGetQuery(con, statement=sql)
# Load the points of all the subject locations with their covariates into the variable, all_subs 
sql=paste("select *  from la_birth_sublocations ",sep="") 
all_subs = dbGetQuery (con, statement=sql) 
# Remove the records of missing values for the covariates selected 
# Definition of covariates: (1) brrawid_integ: unique id of the subject location (integer type); (2) sp_no2: NO2 exposure (weekly average) for the entire pregnancy period (sp_nox: NOx exposure (weekly average) for the entire pregnancy period); (3) age_mother: age of the pregnant mother; (4) ndvi100: NDVI with the 100 m buffer radius; (5) birth_weight: term birth weight (unit: g); (6) gest_len: length of gestation; (7) ethnicity; (8) edu_lev: educational level; (9) parity; (10) p_hea_c: primary health care; (11) int_gen: infant's gender.  
all_subs2=na.omit(all_subs[,c("brrawid_integ","sp_no2","age_mother","ndvi100",
                              "birth_weight","gest_len","edu_lev","ethnicity","parity", "p_hea_c")])
all_subs=all_subs[all_subs$brrawid_integ %in% all_subs2$brrawid_integ,]   
# Normalize the continuous covariates (age_mother and ndvi100) with 0 mean and 1 variance:  
all_subs$age_mother_1=(all_subs$age_mother-mean(all_subs$age_mother))/sd(all_subs$age_mother) 
all_subs$ndvi100_1=(all_subs$ndvi100-mean(all_subs$ndvi100))/sd(all_subs$ndvi100) 
all_subs$gest_len_1=(all_subs$gest_len-mean(all_subs$gest_len))/sd(all_subs$gest_len) 
# Initial modeling by generalized additive models (GAM) to obtain the mean and variance for the interceptions, as the priors for the interception in the JAGS modeling within each tract  
gm1=gam(birth_weight~sp_no2+factor(edu_lev)+factor(ethnicity)+factor(parity)+ factor(p_hea_c)+s(age_mother_1)+s(ndvi100_1)+s(gest_len_1), data = all_subs)   
# Extract the statistics (mean and precision) of the intercept from the GAM modeling for the whole sample, as priors for JAGS modeling 
sm=summary(gm1)  
int_mean=sm$p.coeff["(Intercept)"] 
int_sd=sm$se["(Intercept)"] 
int_precision=1/(int_sd*int_sd)

# Define the data frame (tracts_res) to store the result of Bayesian inference within each tract with its index, itr.  
itr=0 
tracts_res=data.frame()  
# Define the priors for health effects according to the summary (Table S1 and Fig. S3) 
prior_eff_mean=-1.12 
prior_eff_precision=1/82.9 
# Set up the file to store the Bayesian model in JAGS 
jagfl="/tmp/md1.jags" 
burn.in =10000
# Major loop module to conduct Bayesian inference within each Census tract and store the results to the data frame, tracts_res 
for(i in 1:nrow(tracts)){ 
  print(i) 
  atract=tracts[i,"tract"] 
  sql=paste("select t.*  from la_birth_sublocations t inner join (select * from la_tracts where tract='",atract,"')b on ST_Within(t.geom,b.geom) ",sep="")
  tr_subs =dbGetQuery(con, statement=sql)
  tr_subs2=na.omit(tr_subs[,c("brrawid_integ","sp_no2","age_mother","ndvi100",
                              "birth_weight","gest_len","edu_lev","ethnicity","parity", "p_hea_c")])
  tr_subs=tr_subs[tr_subs$brrawid_integ %in% tr_subs2$brrawid_integ,]   
  # Extract the normalized covariates for tr_subs from all_subs   
  match_index=match(tr_subs$tr_subs, all_subs$brrawid_integ)  
  tr_subs$age_mother_1= all_subs[match_index, "age_mother_1"]
  tr_subs$ndvi100_1= all_subs[match_index, " ndvi100_1"]
  tr_subs$ gest_len_1= all_subs[match_index, " gest_len_1"]  
  # Delete the existing JAGS model file 
  unlink(jagfl) 
  # Generate the model file and store into the path, jagfl 
  jd=try(jagam(sp_no2~factor(edu_lev)+factor(ethnicity)+factor(parity)+factor(p_hea_c)
               +s(age_mother_1)+s(ndvi100_1)+s(gest_len_1),data=tr_subs,file=jagfl,
               sp.prior="gamma",diagonalize=TRUE),silent=TRUE)  
  if(class(jd)=="try-error"){
    outStr=paste0("Modeling fails for census tract ",atract,"!")  
    cat(outStr,file="/tmp/failed_tracts.txt", append=TRUE)  
    next 
  } 
  # Revise the model file with the priors for the intercept and health effects included 
  al_lines = readLines(jagfl)
  line_no=grep("b\\[i\\] \\~ dnorm\\(0,0.001\\)", al_lines)
  al_lines=al_lines[-line_no] 
  al_lines=append(al_lines, paste0("  b[1] ~ dnorm(",int_mean,",",int_precision," )"), 
                  after=(line_no-1))  # including the priors for the interceptions 
  al_lines=append(al_lines, paste0("b[2] ~ dnorm(",prior_eff_mean,",",prior_eff_precision,")"), after=(line_no))   # including the priors for the health effects 
  writeLines(al_lines,jagfl) # file.show(jagfl) 
  jm =jags.model(jagfl,data=jd$jags.data,inits=jd$jags.ini,n.chains=9)
  params=c("b","rho","scale") # parameters of interest 
  # Sets a trace monitor, updates the model, and coerces the output to a single mcmc. 
  samps =coda.samples(jm, params, n.iter = 20000,thin=1)  
  # Summarize the inference output 
  asum=summary(window(samps, start =burn.in))  
  # Convergence diagnostic 
  gdiag=gelman.diag(samps) 
  
  if(itr==0){
    itr=1
  }else{
    itr=itr+1   
  }
  tracts_res[itr,"tract"]=atract 
  tracts_res[itr,"i"]=i  
  tracts_res[itr,"Mean_int_mean"]=asum$statistics["b[1]","Mean"] 
  tracts_res[itr,"Mean_int_sd"]=asum$statistics["b[1]","SD"] 
  tracts_res[itr,"Mean_eff_mean"]=asum$statistics["b[2]","Mean"] 
  tracts_res[itr,"Mean_eff_sd"]=asum$statistics["b[2]","SD"] 
}
# Output the results into the file 
write.csv(tracts_res,file="tracts_no2_effects.csv",row.names=FALSE) 
# Then, the output file, tracts.csv was used to create the table, la_tract_map with the tract id and geometry fields from the table, la_tracts.    
