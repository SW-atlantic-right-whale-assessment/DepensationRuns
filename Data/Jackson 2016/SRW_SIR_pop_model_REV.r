# R Script for Bayesian Assessment Model of Southern Right Whales- originally designed for humpback whales using the equations in Zerbini et al. (2011) Journal of Cetacean Research and Management 3: 131-144. All citations and use of this script please cite Zerbini et al. 2011 as the source of this code.

## EXAMPLE DATA AND MAIN FUNCTION CALL ############
# Male captures and pecaptures for fitting into the model, as used in Carroll et al. (2013) Ecol Appl 23:1677-1690.
NZmale.captures<-data.frame(Year=c(1995,1996,1997,1998,2006,2007,2008,2009),captures=c(29,21,31,51,53,60,48,76))
male.recaptures=c(0,4,3,2,0,0,0,0,0,0,2,2,0,0,1,0,0,0,0,6,0,1,0,0,0,0,0,0,0,3,2,0,0,0,0,0,0,5,1,2,0,0,0,0,0,0,4,8,0,0,0,0,0,0,0,9)
NZ.malerecaptures<-matrix(male.recaptures,ncol=8,byrow=TRUE)
Year=c(1995,1996,1997,1998,2006,2007,2008,2009)
NZ.malerecaptures<-data.frame(NZ.malerecaptures)
rownames(NZ.malerecaptures)<-Year[-8]
colnames(NZ.malerecaptures)<-Year

# Female captures and pecaptures for fitting into the model, as used in Carroll et al. (2013) Ecol Appl 23:1677-1690.
NZ.femcaptures<-data.frame(Year=c(1995,1996,1997,1998,2006,2007,2008,2009),captures=c(28,20,19,46,50,86,92,103))
fem.recaptures=c(0,0,2,2,0,3,0,2,0,0,1,2,0,1,1,1,0,0,0,0,0,0,1,3,0,0,0,0,2,1,8,3,0,0,0,0,0,5,4,3,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,3)
NZ.femrecaptures<-matrix(fem.recaptures,ncol=8,byrow=TRUE)
NZ.femrecaptures<-data.frame(NZ.femrecaptures)
rownames(NZ.femrecaptures)<-Year[-8]
colnames(NZ.femrecaptures)<-Year

#Absolute Abundance (used for tests, not used in the final models)
Abs.Abundance.2009=data.frame(Year=c(2009), N.obs=c(2148), CV.obs=c(0.20))
Abs.Abundance.2006=data.frame(Year=c(2006), N.obs=c(1820), CV.obs=c(0.20))
Abs.Abundance.1998=data.frame(Year=c(1998), N.obs=c(1220), CV.obs=c(0.20))
Abs.Abundance.1995=data.frame(Year=c(1995), N.obs=c(1066), CV.obs=c(0.20))

#Relative Abundance (used for tests, not used in the final models)
Rel.Abundance=data.frame(Index=c(1,1,1,1), Year=c(1995, 1998, 2006,2009), IA.obs=c(533,610,910,1074), CV.IA.obs=c(0.20,0.20,0.20,0.20))
Count.Data=data.frame(Index=1)

#Catch Series, from Carroll et al. (2014) PLoS One 9:e93789
# **** Catch_inputs.csv is provided within this package. A full path to this file is required on the line below in order to read it in ****
Catches=read.csv(".../Catch_inputs.csv")
Catch.names=c("Year","CoastNZlow","CoastNZhigh","CoastEA","AmericanBaylow","AmericanBayhigh","AmericanoffNZ","AmericanoffNZSE","AmericanoffEA","AmericanoffEASE","French","SovietEA","SovietNZ")
Catch.data=data.frame(Catches)
colnames(Catch.data)<-Catch.names

# Smearing of US whaling catches through years following departure from Table 3, Carroll et al. (2014) PLoS One 9:e93789
smearing=c(0.007,0.238,0.274,0.262,0.193,0.026)
# Struck but lost rates from coastal and offshore whaling, from Table 8, Carroll et al. (2014) PLoS One 9:e93789
coastal.SL=data.frame(LRF=c(1.27),LRFSE=c(0.05))
offshore.SL=data.frame(LRF=c(1.45),LRFSE=c(0.054))

debug(HUMPBACK.SIR)

#Example calls
# scenario where relative abundance measures are fitted into the model rather than mark recaptures, with no Nfloor, for high case southwest Pacific catches
Popgrowth_SRW_lowNZ=HUMPBACK.SIR(file.name="Popgrowth_SRW_lowNZ", n.samples=NULL, n.resamples=50, prior.r_max=c("uniform", 0, 0.12), r_max.bound=c(0, 0.12), prior.N.obs=c("uniform", 500, 20000), prior.add.CV=c("uniform", 0, 1, FALSE), prior.z=c(NA, 2.39, NA), q.prior.IA=c("uniform", 0, 1, FALSE), q.prior.Count=c("uniform", 0, 1, FALSE), Klim=c(1, 500000), target.Yr=2009, num.haplotypes=12, hap.conversion=0,tolerance.for.bisection=0.0001, output.Yrs=c(2009,2015,2020), abs.abundance.key=TRUE, abs.abundance=Abs.Abundance.2009, rel.abundance=Rel.Abundance, rel.abundance.key=TRUE, MRC.key=c(FALSE, survival=0.84), captures=NZ.captures, recaptures=NZ.recaptures, count.data=Count.Data, count.data.key=FALSE, growth.rate.obs=c(0.074, 0.033, FALSE), growth.rate.Yrs=c(1995, 1996, 1997, 1998), catch.data=Catch.data, serieslow=c(FALSE,"NZEA"),smear=smearing, coast.SL=coastal.SL, offsh.SL=offshore.SL,Threshold=3e-12, Print=0)

## scenario using female captures, Nfloor=36, for low case New Zealand only catches
test.fem_SRW_lowNZ=HUMPBACK.SIR(file.name="test", n.samples=NULL, n.resamples=50, prior.r_max=c("uniform", 0, 0.12), r_max.bound=c(0, 0.12), prior.N.obs=c("uniform", 500, 20000), prior.add.CV=c("uniform", 0, 1, FALSE), prior.z=c(NA, 2.39, NA), q.prior.IA=c("uniform", 0, 1, FALSE), q.prior.Count=c("uniform", 0, 1, FALSE), Klim=c(1, 500000), target.Yr=2009, num.haplotypes=12, hap.conversion=3,tolerance.for.bisection=0.0001, output.Yrs=c(2009,2015,2020), abs.abundance.key=FALSE, abs.abundance=Abs.Abundance.2009, rel.abundance=Rel.Abundance, rel.abundance.key=FALSE, MRC.key=c(TRUE, survival=0.97), captures=NZ.femcaptures, recaptures=NZ.femrecaptures, count.data=Count.Data, count.data.key=FALSE, growth.rate.obs=c(0.074, 0.033, FALSE), growth.rate.Yrs=c(1995, 1996, 1997, 1998), catch.data=Catch.data, serieslow=c(TRUE,"NZ"),smear=smearing, coast.SL=coastal.SL, offsh.SL=offshore.SL,Threshold=3e-4, Print=0)



############ HOUSEKEEPING FUNCTION #############
# HUMPBACK SIR controls the sampling from the priors, the bisection and likelihoods and the output functions
################################################

HUMPBACK.SIR=function(file.name="NULL", n.samples=1000, n.resamples=1000, prior.r_max=c("uniform", 0, 0.12), r_max.bound=c(0, 0.12), prior.N.obs=c("uniform", 0, 50000), prior.add.CV=c("uniform", 0, 1, TRUE), prior.z=c(NA, 2.39, NA), q.prior.IA=c("uniform", 0, 1, FALSE), q.prior.Count=c("uniform", 0, 1, FALSE), Klim=c(1, 500000), target.Yr=2015, num.haplotypes=12, hap.conversion=3,tolerance.for.bisection=0.001, output.Yrs=c(2008), abs.abundance.key=TRUE, abs.abundance=Abs.Abundance, rel.abundance=Rel.Abundance, rel.abundance.key=FALSE, MRC.key=c(FALSE, 1.00),captures=Captures,recaptures=Recaptures, count.data=NULL, count.data.key=FALSE, growth.rate.obs=c(0.074, 0.033, TRUE), growth.rate.Yrs=c(1995, 1996, 1997, 1998), catch.data=Catch.data, serieslow=c(TRUE,"NZ"),smear=NULL, coast.SL=NULL, offsh.SL=NULL,Threshold=1e100, Print=0)
  
  # CALL ELEMENTS
  # file name = name of a file to identified the files exported by the function
  # n.samples = number of samples for the Rubin SIR: NOT USED
  # n.resamples = number of resamples to compute the marginal posterior distributions
  # prior.r_max = prior for r_max. The first element identifies the sampling distribution, the send identifies the lower bound (for uniform distribution) or the mean (for a normal distribution), and the third element corresponds to the upper bound (uniform distribution) or the standard error (normal distribution). The default is Uniform with bounds c(0, 0.12)
  # r_max.bound = bounds for the r_max prior. Default is c(0, 0.12)
  # prior.N.obs = prior distribution for a recent abundance estimate. Elements are equivalent to the prior on r_max
  # prior.CV.add = prior for additional CV if applicable. It is defined by four elements: (1) the sampling distribution, (2) lower bound or mean, (3) upper bound or SE, (4) boolean variable to specify whether CV add is used or not in the likelihood. Default is a uniform distribution bounded by (0,1) and FALSE (=CV.add not used)
  # prior.z = prior on the shape parameter. NOT USED, assumed z=2.39 (max productivity at K=0.6)
  # q.prior.IA = prior on q for indices of abundance. Definition of elements is similar to the prior on CV.add. If the fourth element = FALSE, an analytical solution for q is used (as in Zerbini et al. 2011)
# q.prior.Count = similar to q.prior.IA, but for count data. NOT CURRENTLY USED
# K.lim = bounds for K when preforming the bisection method of Punt and Butterworth (1995). Defined by two elements, the lower and upper bounds. Default is (1, 500000)
# target.Yr = year of the target population estimate for the bisection method. Default is 2008
# num.haplotypes = number of haplotypes to compute lower bound population constraint, Nfloor (from Jackson et al., 2006; 2008 and IWC, 2007)
#hap.conversion = multiplier for Nfloor- in this study it is 3 or 4 to reflect num.haplotypes * 3 or * 4 respectively.
# tolerance.for.bisection = tolerance value for performing the bisection method
# output.Yrs = year for outputing the predicted abundance estimates. Default is 2008, but multiple years can be specified. For example, if outputs for 2005 and 2008 are needed, output.Yrs = c(2005, 2008)
# abs.abundance = R object containing year, estimate of absolute abundance, and CV (see example) NOT USED
# rel.abundance = R object containing years, estimates of relative abudnance and CVs (see example) NOT USED
# rel.abundance.key = key to speficy if relative abundance data are used in the likelihood. Default is TRUE
# MRC.key= key to specify if mark recapture data are used in the likelihood, and survival rates if MRC is true
# number of captures in the MRC dataset
# number of recaptures in the MRC dataset
# count.data= R object containing years, estimates of counts and effort. NOT USED
# count.data.key = key to specify in count data are used. Default is FALSE. NOT USED
# growth.rate.obs = observed growth rate (1st element) and standard error (2nd element) as in Zerbni et al. (2011). If third element is FALSE, the growth rate is not included in the likelihood
# growth.rate.Yrs = Years for which the growth.rate.obs were computed (as in Zerbini et al., 2011)
# catch.data = R object containing the years and catches (see example)
# Threshold = threshold for the McCallister et al. (1994) SIR. This is data-specific. Default is 1e-100.
# serieslow = [1] specifies whether low or high catch scenario is used, TRUE is low, FALSE is high, [2] is either "NZ" for NZ only catches, "NZEA" for NZ plus EA catches or "EA" for EA catches only
# Print = key to print various pieces of information as the code runs. Default is 0 (nothing is printed)

{
  require(utils)
  begin.time=Sys.time()
  
  ################################  
  #Assigning variables
  ################################  
  n.samples=n.samples #not used as I am not doing the Rubin (1988) SIR anymore.
  target.Yr=target.Yr
  start.Yr=catch.data$Year[1] #the first year of the projection is set as the first year in the catch series
  end.Yr=max(tail(catch.data$Year,1), max(abs.abundance$Year), max(rel.abundance$Year)) #the last year of the projection is set as the last year in the catch or abundance series, whichever is most recent
  bisection.Yrs=target.Yr-start.Yr+1 #setting the target year for the bisection method
  projection.Yrs=end.Yr-start.Yr+1 #setting the years to project
  
  z=prior.z[2]
  num.IA=max(rel.abundance$Index) #determining the number of Indices of Abundance available
  num.Count=max(count.data$Index) #determining the number of Count Data sets available
  rel.abundance$Sigma=sqrt(log(1+rel.abundance$CV.IA.obs^2)) #computing the value of sigma as in Zerbini et al. 2011
  if(Print==1){print(paste("Rel.Abundance$Sigma =", rel.abundance$Sigma))}
  #count.data$Sigma=sqrt(log(1+count.data$CV.IA.obs^2)) #computing the value of sigma for the count data as in Zerbini et al. (2011)
  abs.abundance$Sigma=sqrt(log(1+abs.abundance$CV.obs^2)) #computing the value of sigma as in Zerbini et al. 2011
  MVP=hap.conversion*num.haplotypes #computing the minimum viable population, if num.haplotypes=0, assumes no MVP
  
  # checking whether to apply low 'NZ' or 'NZEA' case for coastal catch series
  if (serieslow[1]==TRUE) 
  {
    if (serieslow[2]=="NZ")
    {
      coastal.catch=catch.data$CoastNZlow + catch.data$AmericanBaylow
    }
    if (serieslow[2]=="NZEA")
    {
      coastal.catch=catch.data$CoastNZlow + catch.data$AmericanBaylow + catch.data$CoastEA
    }
    if (serieslow[2]=="EA")
    {
      coastal.catch=catch.data$CoastEA 
    }
    if(Print==1){print(paste("Low catch series"))}
  }	
  #		print(paste("CoastNZlow=",sum(catch.data$CoastNZlow)))
  #		print(paste("CoastAmericanBaylow=",sum(catch.data$AmericanBaylow)))
  # applies high 'NZ' or 'NZEA' case for coastal catch series
  else 
  { 
    if (serieslow[2]=="NZ")
    {
      coastal.catch=catch.data$CoastNZhigh + catch.data$AmericanBayhigh
    }
    if (serieslow[2]=="NZEA")
    {
      coastal.catch=catch.data$CoastNZhigh + catch.data$AmericanBayhigh + catch.data$CoastEA
    }		
    if (serieslow[2]=="EA")
    {
      coastal.catch=catch.data$CoastEA 
    }
    if(Print==1){print(paste("High catch series"))}
  }
  offshore.AmericanNZ=catch.data$AmericanoffNZ
  offshore.AmericanNZSE=catch.data$AmericanoffNZSE
  offshore.AmericanEA=catch.data$AmericanoffEA
  offshore.AmericanEASE=catch.data$AmericanoffEASE
  Soviet.catchNZ=catch.data$SovietNZ
  Soviet.catchEA=catch.data$SovietEA
  French.catch=catch.data$French #this French catch data is already smeared (equally between t+1, t+2, t+3), see Carroll et al. 2014 PLoS One 9:e93789
  AmericanBaymean=(catch.data$AmericanBaylow+catch.data$AmericanBayhigh)/2 # average of low & high American bay whaling catches, used to convert French total catches into 'bay' and 'offshore' components
  
  tol=tolerance.for.bisection #assigning the value for tolerance in computing K using the uniroot function (bisection method)
  i=0 #start the loop
  draw=1 #keep track of number of draws
  Cumulative.Likelihood=0
  
  #Creating output vectors
  #-------------------------------------
  names=c("r_max", "K", "sample.N.obs", "add.CV", "Nmin", "YearMin","Catches","violate_MVP", paste("N", output.Yrs, sep=""), paste("ROI_IA", unique(rel.abundance$Index), sep=""), paste("q_IA", unique(rel.abundance$Index), sep=""), paste("ROI_Count", unique(count.data$Index), sep=""), paste("q_Count", unique(count.data$Index), sep=""), "NLL.IAs", "NLL.Count", "NLL.N", "NLL.GR", "NLL.MRC","NLL", "Likelihood", "Max_Dep", paste("status", output.Yrs, sep=""), "draw", "save")
  
  samples.output=matrix(0, nrow=1, ncol=length(names))
  resamples.output=matrix(0, nrow=1, ncol=length(names))
  resamples.trajectories=matrix(NA, nrow=1, ncol=projection.Yrs)
  resamples.catches=matrix(NA,nrow=1,ncol=projection.Yrs)
  final.trajectory=matrix(NA, nrow=projection.Yrs, ncol=6)
  Year=seq(start.Yr, end.Yr, by=1)
  
  #Initiating the SIR loop
  
  while(i<n.resamples)
  {
    
    #Sampling from Priors
    #-------------------------------
    save=FALSE #variable to indicate whether a specific draw is kept
    
    #Sampling for r_max
    sample.r_max=0.2 #setting sample.r_max outside of the bound
    while(sample.r_max<r_max.bound[1] | sample.r_max>r_max.bound[2])
    {
      sample.r_max=SAMPLE.PRIOR(name=prior.r_max[1], Val.1=prior.r_max[2], Val.2=prior.r_max[3]) #prior on r_max, keep if within boundaries
    }
    
    #sampling from the N.obs prior- this is kept very broad and spans all realistic estimate, so as to be minimally informative
    sample.N.obs=SAMPLE.PRIOR(name=prior.N.obs[1], Val.1=prior.N.obs[2], Val.2=prior.N.obs[3]) #prior on N_obs
    
    if(prior.add.CV[4]==TRUE) #prior on additional CV, not used in this model
    {
      sample.add.CV=SAMPLE.PRIOR(name=prior.add.CV[1], Val.1=prior.add.CV[2], Val.2=prior.add.CV[3])
    }
    else
    {
      sample.add.CV=0
    }
    # sampling from struck and lost priors to form the unique catch series
    offshore.AmericanNZ2=rep(0,times=length(Year))
    offshore.AmericanNZ3=rep(0,times=length(Year))
    offshore.AmericanEA2=rep(0,times=length(Year))
    offshore.AmericanEA3=rep(0,times=length(Year))
    offshore.catch2=rep(0,times=length(Year))
    coastal.catch2=coastal.catch
    French.catch2=rep(0,times=length(Year))
    # this picks a unique American catch series for NZ and EA, given the SEs on each annual estimate.
    for (k in 1:length(Year))
    {
      if  (offshore.AmericanNZ[k]>0)   
      {
        AmericanNZ.corr = rnorm(1,offshore.AmericanNZ[k],offshore.AmericanNZSE[k])
        if (AmericanNZ.corr>0)
        {	
          offshore.AmericanNZ2[k]=AmericanNZ.corr
        } 			
      }
      if  (offshore.AmericanEA[k]>0)   
      {
        AmericanEA.corr = rnorm(1,offshore.AmericanEA[k],offshore.AmericanEASE[k])
        if (AmericanEA.corr>0)
        {	
          offshore.AmericanEA2[k]=AmericanEA.corr
        } 			
      }
    }
    # this smears the American NZ catches over the years and proportions specified by the 'smear' vector
    for (ii in 1:length(Year))
    {	
      if (offshore.AmericanNZ2[ii]>0)
      {
        for (j in 1:length(smear)) 
        {
          offshore.AmericanNZ3[(ii+j)-1]=(offshore.AmericanNZ2[ii]*smear[j])+offshore.AmericanNZ3[(ii+j)-1]
        }
      }			    
    }
    if(Print==1){print(paste("American Offshore NZ catch =", sum(offshore.AmericanNZ3)))}
    
    # this smears the American EA catches over the years and proportions specified by the 'smear' vector 
    for (jj in 1:length(Year))
    {	
      if (offshore.AmericanEA2[jj]>0)
      {
        for (ss in 1:length(smear)) 
        {
          offshore.AmericanEA3[(jj+ss)-1]=(offshore.AmericanEA2[jj]*smear[ss])+offshore.AmericanEA3[(jj+ss)-1]
        }
      }			    
    }
    if (serieslow[2]=="NZEA")
    {
      if(Print==1){print(paste("American Offshore EA catch =", sum(offshore.AmericanEA3)))}
    }		    
    
    # this breaks down the French catch series into component types of whaling, based on relative proportion estimated for the American campaign
    France.offEA=(French.catch*offshore.AmericanEA3)/(AmericanBaymean+offshore.AmericanEA3+offshore.AmericanNZ3)
    
    #this gets rid of any NaN in the series
    for (zz in 1:length(Year))
    {
      if (!is.finite(France.offEA[zz]))
      {
        France.offEA[zz]=0
      }
    }
    if ((serieslow[2]=="NZEA")||(serieslow[2]=="EA"))
    {
      if(Print==1){print(paste("French Offshore EA catch =", sum(France.offEA)))}    
    }
    
    # this derives the French catch off NZ and breaks down into relative bay and offshore components, as above.	
    France.offNZ=(French.catch-France.offEA)*(offshore.AmericanNZ3/(offshore.AmericanNZ3+AmericanBaymean))
    for (zz in 1:length(Year))
    {
      if (!is.finite(France.offNZ[zz]))
      {
        France.offNZ[zz]=0
      }
    }
    France.bayNZ=(French.catch-France.offEA)*(AmericanBaymean/(offshore.AmericanNZ3+AmericanBaymean))
    for (zz in 1:length(Year))
    {
      if (!is.finite(France.bayNZ[zz]))
      {
        France.bayNZ[zz]=0
      }
    }
    if(Print==1){print(paste("French Bay NZ catch =", sum(France.bayNZ)))}	
    
    # this sums coastal and offshore catch components for each 'NZ', 'EA' and 'NZEA' catch scenario	 
    if (serieslow[2]=="NZ")
    {	    
      coastal.catch2=coastal.catch2+France.bayNZ
      offshore.catch=France.offNZ+offshore.AmericanNZ3
    }	
    if (serieslow[2]=="NZEA")
    {
      coastal.catch2=coastal.catch2+France.bayNZ
      offshore.catch=France.offNZ+offshore.AmericanNZ3+France.offEA+offshore.AmericanEA3
    }	
    if (serieslow[2]=="EA")
    {
      offshore.catch=France.offEA+offshore.AmericanEA3
    }	
    if(Print==1)
    {print(paste("offshore no SL=",sum(offshore.catch)," coastal no SL=",sum(coastal.catch2)," French coastal no SL=",sum(France.bayNZ)))
    }
    
    #this adds the struck but lost (SL) correction to offshore catch (American and French offshore whaling in NZ region) after checking the SL is >0  
    for (xv in 1:length(Year))
    {
      if (offshore.catch[xv]>0)
      {
        repeat 
        {
          offshore.corr=rnorm(1,offsh.SL[[1]],offsh.SL[[2]])
          if (offshore.corr>0) {break}
        }
        offshore.catch[xv]=offshore.corr*offshore.catch[xv]
      } 			
    }
    
    #this adds the SL correction to coastal catch (coastal NZ landings plus American bay whaling)
    for (k in 1:length(Year))
    {
      if  (coastal.catch2[k]>0)
      {	
        repeat 
        {  				
          Coast.corr =rnorm(1,coast.SL[[1]],coast.SL[[2]])
          if (Coast.corr>0) {break}
        }
        coastal.catch2[k]= Coast.corr*coastal.catch2[k]  
      }
    }
    
    # this builds the total catch series for the run
    catches=c()
    totcatch=0
    for (m in 1:length(Year))
    {
      #		catches[m]=coastal.catch2[m]+offshore.catch2[m]+offshore.American3[m]+Soviet.catch[m]+French.catch2[m]
      if (serieslow[2]=="NZ")
      {
        catches[m]=coastal.catch2[m]+offshore.catch[m]+Soviet.catchNZ[m]
      }
      if (serieslow[2]=="NZEA")
      {
        catches[m]=coastal.catch2[m]+offshore.catch[m]+ Soviet.catchNZ[m]+Soviet.catchEA[m]
      }
      if (serieslow[2]=="EA")
      {
        catches[m]=coastal.catch2[m]+offshore.catch[m]+ Soviet.catchEA[m]
      }
      totcatch=totcatch+catches[m]
    }
    if(Print==1){print(paste("total coastal=",sum(coastal.catch2)," total offshore=",sum(offshore.catch)))}
    
    #Sampling from q priors if q.prior is TRUE (not used in this model)
    if(q.prior.IA[4]==TRUE) #priors on q for indices of abundance
    {  
      q.sample.IA=rep(NA, num.IA)  
      
      for (i in 1:num.IA)
      {
        q.sample.IA[i]=SAMPLE.PRIOR(name=q.prior.IA[1], Val.1=q.prior.IA[2], Val.2=q.prior.IA[3])
      }
    }  
    else
    {
      q.sample.IA=rep(-9999, length(unique(rel.abundance$Index)))
    }  
    
    if(q.prior.Count[4]==TRUE) #priors on q for count data
    {  
      q.sample.Count=rep(NA, 1:num.Count)
      
      for (i in num.Count)
      {
        q.sample.Count[i]=SAMPLE.PRIOR(name=q.prior.Count[1], Val.1=q.prior.Count[2], Val.2=q.prior.Count[3])
      }
    }
    else
    {
      q.sample.Count=rep(-9999, length(unique(count.data$Index)))
    }
    if(Print==1){print(paste("got here"))}   
    
    # this retreives K (or N[0]) using the logistic model and applying the bisection/uniroot method
    
    sample.K=LOGISTIC.BISECTION.K(K.low=Klim[1], K.high=Klim[2], r_max=sample.r_max, z=z, num.Yrs=bisection.Yrs, start.Yr=start.Yr, target.Pop=sample.N.obs, catches=catches, MVP=MVP, tol=tol) #K_prior computed with the bisection/uniroot method
    if(Print==1){print(paste("sample K", sample.K))} 
    
    #Computing the predicted abundances with the samples from the priors
    #----------------------------------------
    Pred.N=GENERALIZED.LOGISTIC(r_max=sample.r_max, K=sample.K, N1=sample.K, z=z, start.Yr=start.Yr, num.Yrs=projection.Yrs, catches=catches, MVP=MVP) 
    if(Print==1){print(paste("sample.r_max=", sample.r_max, "sample.N.obs=", sample.N.obs, "sample.K=", sample.K, "Pred.N.target=", Pred.N$Pred.N[bisection.Yrs]))}
    
    
    #Computing the predicted ROI for the IAs and Count data, if applicable (not used for this model)
    #----------------------------------------    
    #For IAs
    if(rel.abundance.key==TRUE)
    {
      Pred.ROI.IA=COMPUTING.ROI(data=rel.abundance, Pred.N=Pred.N$Pred.N, start.Yr=start.Yr)
      if(Print==1){print(paste("Pred.ROI.IA=",Pred.ROI.IA))}
    }
    else
    {
      Pred.ROI.IA=rep(0, num.IA)
    }
    
    #For Count Data
    if(count.data.key==TRUE)
    {
      Pred.ROI.Count=COMPUTING.ROI(data=count.data, Pred.N=Pred.N$Pred.N, start.Yr=start.Yr)  
    }  
    else
    {
      Pred.ROI.Count=rep(0, num.Count) 
    }
    
    
    
    #Calculate Analytical Qs if rel.abundance.key is TRUE (not used in this model)
    #---------------------------------------------------------
    if(rel.abundance.key==TRUE)
      if(q.prior.IA[4]==FALSE)
      {
        q.sample.IA=CALC.ANALYTIC.Q(rel.abundance, Pred.N$Pred.N, start.Yr, sample.add.CV, num.IA)
      }
    else
    {
      q.sample.IA=q.sample.IA
    }   
    if(Print==1){print(paste("q.sampleIA=",q.sample.IA))}
    
    #browser()
    
    #Calculate Analytical Qs if count.data.key is TRUE (NOT USED YET - AZerbini, Feb 2013)
    #--------------------------------------------------------
    if(rel.abundance.key==TRUE)
      if(q.prior.Count[4]==TRUE)
      {
        q.sample.Count=CALC.ANALYTIC.Q(count.data, Pred.N$Pred.N, start.Yr, sample.add.CV, num.Count)
      }
    else
    {
      q.sample.Count=q.sample.Count
    }   
    
    if(Print==1){print(paste("q.IAs=", q.sample.IA)); print(paste("q.Count=", q.sample.Count))}
    
    
    #Compute the likelihoods
    #--------------------------------   
    # (1) relative indices (if rel.abundance.key is TRUE)
    if(rel.abundance.key==TRUE)
    {
      lnlike.IAs=LNLIKE.IAs(rel.abundance, Pred.N$Pred.N, start.Yr, q.sample.IA, sample.add.CV, num.IA, log=TRUE)
    }
    else
    {
      lnlike.IAs=0
    }
    if(Print==1){print(paste("lnlike.IAs=", lnlike.IAs))}
    
    # (2) count data (if count.data.key is TRUE, not used in this model)
    if(count.data.key==TRUE)
    {
      lnlike.Count=LNLIKE.IAs(count.data, Pred.N$Pred.N, start.Yr, q.sample.Count, sample.add.CV, num.Count, log=TRUE)
    }
    else
    {
      lnlike.Count=0
    }
    if(Print==1){print(paste("lnlike.Count=", lnlike.Count))}
    
    # (3) absolute abundance
    if(abs.abundance.key==TRUE)
    { 
      lnlike.Ns=LNLIKE.Ns(abs.abundance, Pred.N$Pred.N, start.Yr, sample.add.CV, log=TRUE)
    }
    else
    {
      lnlike.Ns=c(0,0)
    }
    if(Print==1){print(paste("lnlike.Ns=", lnlike.Ns))}
    
    # (4) growth rate if applicable
    if(growth.rate.obs[3]==TRUE)
    {
      Pred.GR=PRED.GROWTH.RATE(growth.rate.Yrs=growth.rate.Yrs, Pred.N=Pred.N$Pred.N, start.Yr=start.Yr)
      lnlike.GR=LNLIKE.GR(Obs.GR=growth.rate.obs[1], Pred.GR=Pred.GR, GR.SD.Obs=growth.rate.obs[2])  
    }  
    else
    {
      lnlike.GR=0
    }    
    if(Print==1){print(paste("lnlike.GR=", lnlike.GR))}
    
    #(5) Mark recapture data if applicable (used in this model)
    if(MRC.key[1]==TRUE)
    {
      lnlike.MRC=LNLIKE.MR(Pred.N$Pred.N,MRC.key[2],captures,recaptures,start.Yr)
    }
    else
    {
      lnlike.MRC=0
    }
    if(Print==1){print(paste("lnlike.MRC=", lnlike.MRC))}
    
    # sums likelihood components   
    LL=lnlike.IAs[[1]]+lnlike.Count[[1]]+lnlike.Ns[[1]]+lnlike.GR[[1]]+lnlike.MRC#these use the likelihoods in Zerbini et al. (2011)
    Likelihood=exp(-LL)
    if(Print==1){print(paste("NLL=", LL, "Likelihood=", Likelihood))}
    
    # checks whether Nfloor has been violated and reports   
    if(Pred.N$Violate.MVP==TRUE)
    {
      Likelihood=0
      print(paste("MVP violated on draw", draw))
    }
    
    # measures cumulative likelihood for resampling
    Cumulative.Likelihood=Cumulative.Likelihood+Likelihood
    print(paste("Sample=", i, "Likelihood=", Likelihood))
    
    #Resampling of draws that take cumulative likelihood over the threshhold value (see Text S1 for threshold values used in each model)
    # accrues resample (i.e. posterior) trajectory, catch and output parameter files
    if(Pred.N$Violate.MVP==FALSE)
    {
      while(Cumulative.Likelihood>Threshold)
      {
        print(paste("sample=", i, " & draw=", draw))
        if(Print==1){print(paste("draw=", draw, "Likelihood=", Likelihood, "Cumulative=", Cumulative.Likelihood))}
        save=TRUE
        Cumulative.Likelihood=Cumulative.Likelihood-Threshold
        resamples.trajectories=rbind(resamples.trajectories, Pred.N$Pred.N)
        resamples.catches=rbind(resamples.catches,catches)
        resamples.output=rbind(resamples.output, c(sample.r_max, sample.K, sample.N.obs, sample.add.CV, Pred.N$Min.Pop, Pred.N$Min.Yr[[1]], totcatch, Pred.N$Violate.MVP, c(Pred.N$Pred.N[output.Yrs-start.Yr+1]), Pred.ROI.IA, q.sample.IA, Pred.ROI.Count, q.sample.Count, lnlike.IAs[[1]], lnlike.Count[[1]], lnlike.Ns[[1]], lnlike.GR[[1]], lnlike.MRC,LL, Likelihood, Pred.N$Min.Pop/sample.K, c(Pred.N$Pred.N[output.Yrs-start.Yr+1]/sample.K), draw, save))
        i=i+1
      }
    }  
    # accrues sample (prior) output parameter files
    samples.output=rbind(samples.output, c(sample.r_max, sample.K, sample.N.obs, sample.add.CV, Pred.N$Min.Pop, Pred.N$Min.Yr[[1]], totcatch,Pred.N$Violate.MVP, c(Pred.N$Pred.N[output.Yrs-start.Yr+1]), Pred.ROI.IA, q.sample.IA, Pred.ROI.Count, q.sample.Count, lnlike.IAs[[1]], lnlike.Count[[1]], lnlike.Ns[[1]], lnlike.GR[[1]], lnlike.MRC,LL, Likelihood, Pred.N$Min.Pop/sample.K, c(Pred.N$Pred.N[output.Yrs-start.Yr+1]/sample.K), draw, save))
    # iterates draw until the given number of resamples (n.resamples) is reached
    draw=draw+1
  }  
  
  # constructs output summary files
  
  samples.output=as.data.frame(samples.output)
  names(samples.output)=names
  samples.output=samples.output[-1,]
  #samples.output.summary=SUMMARY.SIR(x=samples.output, scenario=file.name)
  write.csv(samples.output, paste(file.name,"_","samples.output.csv", sep=""))
  
  resamples.output=as.data.frame(resamples.output)
  names(resamples.output)=names
  resamples.output=resamples.output[-1,]
  #resamples.output.summary=SUMMARY.SIR(x=resamples.output, scenario=file.name)
  write.csv(resamples.output, paste(file.name,"_","resamples.output.csv", sep=""))
  
  resamples.trajectories=data.frame(resamples.trajectories)
  names(resamples.trajectories)=seq(start.Yr, end.Yr, 1)
  resamples.trajectories=resamples.trajectories[-1,]
  write.csv(resamples.trajectories, paste(file.name,"_","resample.trajectories.csv", sep=""))
  
  resamples.catches=data.frame(resamples.catches)
  names(resamples.catches)=seq(start.Yr, end.Yr, 1)
  resamples.catches=resamples.catches[-1,]
  write.csv(resamples.catches, paste(file.name,"_","resample.catches.csv", sep="",row.names=NULL))
  
  final.trajectory[,1]=sapply(resamples.trajectories, mean)
  final.trajectory[,2]=sapply(resamples.trajectories, median)
  final.trajectory[,3]=sapply(resamples.trajectories, quantile, probs=c(0.025))
  final.trajectory[,4]=sapply(resamples.trajectories, quantile, probs=c(0.975))
  final.trajectory[,5]=sapply(resamples.trajectories, quantile, probs=c(0.05))
  final.trajectory[,6]=sapply(resamples.trajectories, quantile, probs=c(0.95))
  final.trajectory=data.frame(final.trajectory)
  names(final.trajectory)= c("mean", "median", "PI.2.5%", "PI.97.5%", "PI.5%", "PI.95%")
  final.trajectory=data.frame(Year, final.trajectory)
  
  # calculates output ratio of resamples per sample, listed to screen below
  resamples.per.samples=dim(samples.output)[1]/dim(resamples.output)[1]
  
  end.time=Sys.time()
  print(paste("Time to Compute=", (end.time-begin.time), " minutes"))
  
  return(list(call=call, file.name=file.name, Date.Time=Sys.time(), Time.to.compute.in.minutes=paste(end.time-begin.time), Threshold=Threshold, Ratio.Resamples.per.Sample=paste("1 resample",":", resamples.per.samples, "samples"),  resamples.output=resamples.output,final.trajectory=final.trajectory, inputs=list(draws=draw, n.resamples=n.resamples, prior.r_max=prior.r_max, prior.N.obs=prior.N.obs, target.Yr=target.Yr, MVP=paste("num.haplotypes=", num.haplotypes, "MVP=", 3*num.haplotypes), tolerance=tolerance.for.bisection, output.Years=output.Yrs))) 
}

#### END OF FUNCTION ###################################
# -------------------------

#################################################
# FUNCTION TO SAMPLE FROM PRIORS
# This function samples from 4 distributions given two input values (Val.1 and Val.2). If the distribution
# is uniform or log-uniform, then Val.1 and Val.2 correspond to the lower and upper bounds. If normal
# or log-normal, Val.1 corresponds to the mean and Val.2 to the SD.
# -----------------------------------------------
SAMPLE.PRIOR=function(name=NA, Val.1=NA, Val.2=NA)
{
  #uniform prior
  if(name=="uniform"){return(runif(1,as.numeric(Val.1), as.numeric(Val.2)))}
  #log-uniform prior
  else if(name=="log-uniform"){return(exp(runif(1, as.numeric(Val.1), as.numeric(Val.2))))}
  #normal prior
  else if(name=="normal"){return(rnorm(1, as.numeric(Val.1), as.numeric(Val.2)))}
  #log-normal prior
  else if(name=="log-normal"){return(rlnorm(1, as.numeric(Val.1), as.numeric(Val.2)))}
}

# END OF FUNCTION
# -------------------------


###################################################
# GENERALISED LOGISTIC MODEL  
# This function does the population projection using a Pella-Tomlison population dynamics model
###################################################
GENERALIZED.LOGISTIC=function(r_max, K, N1, z, start.Yr, num.Yrs, catches, MVP=0)
{
  Violate.Min.Viable.Pop=FALSE #variable to indicate whether min population is reached
  Pred.N=rep(NA, num.Yrs) #create a vector to hold the model predicted population size
  Pred.N[1]=N1 #The first year in the vector above is N1
  
  #Project the population
  for(t in 1:(num.Yrs-1))
  {
    Pred.N[t+1]=max(Pred.N[t]+r_max*Pred.N[t]*(1-(Pred.N[t]/K)^z)-catches[t],1)
  }
  
  Min.Pop=min(Pred.N) #Compute the Nmin
  Min.Yr=which(Pred.N==Min.Pop)+start.Yr-1 #Compute the year at which Nmin occurred
  if(Min.Pop<MVP) {Violate.Min.Viable.Pop=TRUE} #Determine whether Nmin is below Min Viable Population
  
  return(list(Min.Pop=Min.Pop, Min.Yr=Min.Yr[[1]], Violate.MVP=Violate.Min.Viable.Pop, Pred.N=Pred.N))
}
#END FUNCTION
#---------------------------------

#######################################################
#THiS FUNCTION COMPUTES THE PREDICTED GROWTH RATE IF SUCH INFORMATION IS AVAILABLE
#FROM AN INDEPENDENT ESTIMATE (see Zerbini et al., 2011) AND WILL NOT BE 
#COMPUTED FROM INPUT DATA IN THIS MODEL
#######################################################
PRED.GROWTH.RATE=function(growth.rate.Yrs, Pred.N, start.Yr=start.Yr)
{
  GR.Yrs=growth.rate.Yrs-start.Yr+1 #computing the growth rate years
  Pred.N.GR=Pred.N[GR.Yrs]
  
  Pred.GR=(log(Pred.N.GR[length(Pred.N.GR)])-log(Pred.N.GR[1]))/(length(Pred.N.GR)-1)
  
  return(Pred.GR)
}
# END OF FUNCTION
#----------------------------------------

########################################  
#THIS FUNCTION COMPUTES THE PREDICTED RATE OF INCREASE FOR A SET OF SPECIFIED
#YEARS FOR COMPARISON WITH TRENDS ESTIMATED SEPARATELY WITH ANY OF THE INDICES OF
#ABUNDANCE OR COUNT DATA
########################################
COMPUTING.ROI=function(data=data, Pred.N=Pred.N, start.Yr=NULL)
{
  num.indices=max(data$Index)
  Pred.ROI=rep(NA, num.indices)
  
  for(i in 1:num.indices)
  {
    index.ini.year=(head(subset(data, Index==i)$Year,1)-start.Yr)
    index.final.year=(tail(subset(data, Index==i)$Year, 1)-start.Yr)
    elapsed.years=index.final.year-index.ini.year
    
    Pred.ROI[i]=exp((log(Pred.N[index.final.year])-log(Pred.N[index.ini.year]))/(elapsed.years))-1
    
  }
  return(Pred.ROI=Pred.ROI)
}
#END OF FUNCTION
#--------------------------


####################################################################  
# FUNCTION TO CALCULATE A TARGET K FOR THE BISECTION METHOD
###################################################################  

#example call: TARGET.K(r_max, K, N1, z, start.Yr=start.Yr, num.Yrs=bisection.Yrs, target.Pop=target.Pop, catches=catches, MVP=MVP)
#  TARGET.K(r_max=sample.r_max, K, N1, z, start.Yr=start.Yr, num.Yrs=bisection.Yrs, target.Pop=sample.N.obs, catches=catches, MVP=MVP)

TARGET.K = function(r_max, K, N1, z, num.Yrs, start.Yr, target.Pop, catches=catches, MVP=0)
{
  Pred.N=GENERALIZED.LOGISTIC(r_max=r_max, K=K, N1=K, z=z, start.Yr=start.Yr, num.Yrs=num.Yrs, catches=catches, MVP=MVP)
  diff=Pred.N$Pred.N[num.Yrs]-target.Pop
  #plot(seq(1901, 1900+num.Yrs, 1), Pred.N$Pred.N, ylim=c(0,50000))
  #print(paste("K=", K, "Target.Pop=", target.Pop, "N-Target=", Pred.N$Pred.N[num.Yrs], "Total catch=",sum(catches)))
  return(diff=diff)  
  #print(paste("Diff =",diff))
}
#END FUNCTION
#--------------------------------- 


#####################################################################
# THIS FUNCTION DOES THE LOGISTIC BISECTION
#####################################################################
#ex call: LOGISTIC.BISECTION.K(K.low=1, K.high=100000, r_max=r_max, z=z, num.Yrs=bisection.Yrs, start.Yr=start.Yr, target.Pop=target.Pop, catches=catches, MVP=MVP, tol=0.001)
LOGISTIC.BISECTION.K=function(K.low, K.high, r_max, z, num.Yrs, start.Yr, target.Pop, catches, MVP, tol=0.001)
{
  Kmin=uniroot(TARGET.K, tol=tol, c(K.low, K.high), r_max=r_max, z=z, num.Yrs=num.Yrs, start.Yr=start.Yr, target.Pop=target.Pop, catches=catches, MVP=MVP)
  #print(paste("Kmin=", Kmin))
  return(Kmin$root)
}
#END OF FUNCTION
#------------------------

######################################################
# THIS FUNCTION COMPUTES THE ANALYTICAL ESTIMATES OF Q
######################################################
CALC.ANALYTIC.Q=function(rel.Abundance, Pred.N, start.Yr, add.CV=0, num.IA)
{
  analytic.Q=rep(NA, num.IA) #vector to store the q values
  
  for(i in 1:num.IA)
  {
    IA=subset(rel.Abundance, Index==i) #subseting across each index of abundance
    IA.yrs=IA$Year-start.Yr+1 #years for which IAs are available
    IA$Sigma=sqrt(log(1+IA$CV.IA.obs^2)) #computing the value of sigma as in Zerbini et al. 2011
    qNumerator=sum((log(IA$IA.obs/Pred.N[IA.yrs]))/(IA$Sigma*IA$Sigma+add.CV*add.CV)) #numerator of the analytic q estimator (Zerbini et al., 2011 - eq. (3))
    qDenominator=sum(1/(IA$Sigma*IA$Sigma)) #denominator of the analytic q estimator (Zerbini et al., 2011 - eq. (3))
    analytic.Q[i]=exp(qNumerator/qDenominator) #estimate of q
  }
  return(analytic.Q=analytic.Q)  
}
# END OF FUNCTION 
#-------------------------------------------------------  

################################################################################
# THIS FUNCTION COMPUTES THE LN LIKELIHOOD OF THE INDICES OF ABUNDANCE
################################################################################
LNLIKE.IAs=function(Rel.Abundance, Pred.N, start.Yr, q.values, add.CV, num.IA, log=TRUE)
{
  loglike.IA1=0
  loglike.IA2=0
  
  for(i in 1:num.IA)
  {
    IA=subset(Rel.Abundance, Index==i) #subseting across each index of abundance
    IA.yrs=IA$Year-start.Yr+1 #years for which IAs are available
    loglike.IA1=loglike.IA1+((sum(log(IA$Sigma)+log(IA$IA.obs)+0.5*((((log(q.values[i]*Pred.N[IA.yrs])-log(IA$IA.obs))^2)/(IA$Sigma*IA$Sigma)))))) #this is the likelihood from Zerbini et al. 2011 (eq. 5) 
    
    loglike.IA2=loglike.IA2 + CALC.LNLIKE(Obs.N=IA$IA.obs, Pred.N=(q.values[i]*Pred.N[IA.yrs]), CV= sqrt(IA$Sigma*IA$Sigma + add.CV*add.CV), log=log) #this is the log-normal distribution from R (using function dnorm)
    
  }
  return(list(loglike.IA1=loglike.IA1, loglike.IA2=loglike.IA2))
  
}
# END OF FUNCTION
################################################################################

################################################################################
# THIS FUNCTION COMPUTES THE LN LIKELIHOOD OF THE ABSOLUTE ABUNDANCE
################################################################################
LNLIKE.Ns=function(Obs.N, Pred.N, start.Yr, add.CV, log=TRUE)
{
  loglike.Ns1=0
  loglike.Ns2=0
  
  N.yrs=Obs.N$Year-start.Yr+1 #years for which Ns are available
  loglike.Ns1=loglike.Ns1+((sum(log(Obs.N$Sigma)+log(Obs.N$N.obs)+0.5*((((log(Pred.N[N.yrs])-log(Obs.N$N.obs))^2)/(Obs.N$Sigma*Obs.N$Sigma+add.CV*add.CV)))))) #this is the likelihood from Zerbini et al. 2011 (eq. 4) 
  loglike.Ns2=loglike.Ns2 + CALC.LNLIKE(Obs.N=Obs.N$N.obs, Pred.N=(Pred.N[N.yrs]), CV= sqrt(Obs.N$Sigma*Obs.N$Sigma + add.CV*add.CV), log=log) #this is the log-normal distribution from R (using function dnorm)
  
  return(list(loglike.Ns1=loglike.Ns1, loglike.Ns2=loglike.Ns2))
  
}

# END OF FUNCTION
################################################################################

################################################################################
# THIS FUNCTION COMPUTES THE LN LIKELIHOOD OF THE GROWTH RATE
################################################################################
LNLIKE.GR=function(Obs.GR, Pred.GR, GR.SD.Obs)
{
  loglike.GR1=0
  loglike.GR2=0
  
  loglike.GR1=loglike.GR1+(((log(GR.SD.Obs)+0.5*(((Pred.GR-Obs.GR)/GR.SD.Obs)^2)))) #this is the likelihood from Zerbini et al. 2011 (eq. 6) 
  
  loglike.GR2=loglike.GR2 + CALC.LNLIKE(Obs.N=Obs.GR, Pred.N=Pred.GR, CV=GR.SD.Obs, log=FALSE)
  
  return(list(loglike.GR1=loglike.GR1, loglike.GR2=loglike.GR2))
  
}
# END OF FUNCTION
################################################################################

################################################################################
# THIS FUNCTION COMPUTES THE LOG LIKELIHOOD 
################################################################################
CALC.LNLIKE=function(Obs.N, Pred.N, CV, log=F)
{
  return(sum(dnorm(x=log(Obs.N), mean=log(Pred.N), sd=CV, log=F)))
}
# END OF FUNCTION
################################################################################
################################################################################
# THIS FUNCTION COMPUTES THE LN LIKELIHOOD OF A MARK RECAPTURE TREND (following Johnston and Butterworth 2008 Paper SC/60/SH37 presented to the IWC Scientific Committee. Available from www.iwc.int.)
################################################################################
LNLIKE.MR=function(Pred.N,survival, n.captures,n.recaptures,start.Yr) 
{
  modelprob.capture <- c()
  N.yrs=n.captures$Year-start.Yr+1 #years for which Ns are available
  years=length(N.yrs)
  modelprob.capture <- n.captures$captures/(Pred.N[N.yrs]/2)
  modelprob.recapture <- array(0,dim=c(years,years))	
  for (k in 1:(years-1)) 
  {
    for (l in (k+1):years) 
    {
      modelprob.recapture[k,l] <- modelprob.capture[k] * modelprob.capture[l] * (Pred.N[N.yrs[k]]/2) * exp(-(1-survival)*(N.yrs[l]-N.yrs[k]))
    }
  }
  total.logs <- 0
  for (m in 1:(years-1)) 
  {
    for (n in (m+1):years) 
    {
      total.logs <- total.logs + (((-1*n.recaptures[m,n])*log(modelprob.recapture[m,n]))+modelprob.recapture[m,n])
    }
  }
  return(total.logs)
}
# END OF FUNCTION
################################################################################