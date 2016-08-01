##### PAH Data Summaries

##Load Packages
library(NADA)
library(dplyr)
library(boot)
library(stringr)
library(EnvStats)
library(readr)
library(RODBC)
library(ggplot2)
library(scales)
library(grid)
library(tidyr)
#=================
##Define Terms
#=================
Percent_Cutoff   = 20
Boot_Repeats = 2000
repeats    = 2000
seed       = 27

#source("X:/Programs/Air_Quality_Programs/Air Focus Areas/EPA Community Scale Monitoring of PAHs Project/Data Analysis/Data Processing/PAH Data processing.R")


##Pull in the allPAHconcentrations.csv. This file is created with this script: PAH Data processing.R and contains gas, particle, total and passive concentrations with appropriate censoring.
pahconcdata <- read.csv("X:/Programs/Air_Quality_Programs/Air Focus Areas/EPA Community Scale Monitoring of PAHs Project/Data Analysis/Data Processing/allPAHconcentrations.csv", stringsAsFactors=FALSE)

##=======================================================================================
##Adjust date formatting and take out the in CAS numbers
##Pull in only columns that are needed for these summaries and rename columns
##=======================================================================================
setwd("X:/Programs/Air_Quality_Programs/Air Focus Areas/EPA Community Scale Monitoring of PAHs Project/Data Analysis/Annual and Seasonal Statistics and Summaries")
pahconcdata<-pahconcdata[, c("Sampler_Type","Result", "MDL", "Start_Run_Date","Analyte", "Units","Year", "MPCA_Site_ID", "CAS", "Season", "Study_Year", "Primary_Occurence_Code", "location_if_collocated", "Study_Location", "Censored", "Season_StudyYear")]
pahconcdata$CAS <- gsub("a", "", pahconcdata$CAS)

pahconcdata$Start_Run_Date<-as.character(pahconcdata$Start_Run_Date)
pahconcdata$Start_Run_Date <- as.Date(pahconcdata$Start_Run_Date, format = "%Y-%m-%d")
colnames(pahconcdata)<-c("Sampler_Type", "result","MDL", "date","Analyte","unit","Year", "siteid", "CAS", "Season", "Study_Year", "poc", "location_if_collocated", "Sampling_Location", "censored", "Season_StudyYear")

##======================================================
#######Average Same Site POCs
##======================================================
pahconcdata<-group_by(pahconcdata, Sampler_Type, date, Analyte, unit, Year, siteid, CAS, Season, Study_Year, censored, Season_StudyYear) %>% summarise(result=mean(result, na.rm=T), MDL=mean(MDL, na.rm=T)) %>% ungroup()

##==================================================================
#######Take out Year and Dates prior to time grouping and averaging
##==================================================================
pahconcdata<- pahconcdata[, c("Sampler_Type", "Analyte", "unit", "siteid", "CAS", "Season", "Study_Year", "censored", "Season_StudyYear", "MDL", "result")]

##==================================================================================
#######Create Time Based Summary grouping variables: entire study, study year, season, season_study year
##==================================================================================
entirestudy <- mutate(pahconcdata, Grouping_Factor = paste0(Sampler_Type, "_", CAS, "_", siteid), Time_Summary = "Whole Study", Season_StudyYear = "All Seasons and Study Years", Season="All Seasons", Study_Year="All Years")

season <-mutate(pahconcdata, Grouping_Factor = paste0(Sampler_Type, "_", CAS, "_", siteid, "_", Season), Time_Summary = "Season", Season_StudyYear = paste0(Season, "All Years"), Study_Year = "All Years")

seasonstudyyear <-mutate(pahconcdata, Grouping_Factor = paste0(Sampler_Type, "_", CAS, "_", siteid, "_", Season_StudyYear), Time_Summary = "Season Study Year")

studyyear <-mutate(pahconcdata, Grouping_Factor = paste0(Sampler_Type, "_", CAS, "_", siteid, "_", Study_Year), Time_Summary = "Study Year", Season_StudyYear = paste0("All Seasons", Study_Year), Season="All Seasons")

##=====================================================================
#######R Bind the Data Together with the 4 Tiem Based Grouping Factors
##=====================================================================

pahconcdata <- rbind(entirestudy, season, seasonstudyyear, studyyear)
pahconcdata$Grouping_Factor <- as.character(pahconcdata$Grouping_Factor)

##=================================================
#Create a copy of the data
##=================================================
pahconcdata2 <- pahconcdata

#==============================================================================
#Calculate total and percent above detection limit and censored values
#==============================================================================
pahconcdata2$censored <- as.logical(pahconcdata2$censored)
pahconcdata2$CensoredDummy <- ifelse(pahconcdata2$censored==TRUE,1,0)
pahconcdata2$ResultZero <- ifelse(pahconcdata2$CensoredDummy==1, 0, pahconcdata2$result)
detectandcensor <- group_by(pahconcdata2, Grouping_Factor) %>% 
  mutate(Count_Sampled  = n(), 
         Count_Censored = sum(CensoredDummy, na.rm=T),
         Count_AboveDetection   = Count_Sampled - Count_Censored,
         Pct_AboveDetection = round(100 * Count_AboveDetection/Count_Sampled),
         Pct_Censored = round(100*Count_Censored/Count_Sampled),
         uniQ= length(unique(result[!censored])))

##====================================================================================
#Create two datasets with different rules for calculating basic statistical summaries
##====================================================================================
passivesnozeroes <- filter(detectandcensor, Sampler_Type %in% "Passive", Count_AboveDetection>=1)
activesnozeroes <- filter(detectandcensor, !Sampler_Type %in% "Passive", Pct_AboveDetection>=20)

activeszeroes <- filter(detectandcensor, !Sampler_Type %in% "Passive", Pct_AboveDetection<20)
passiveszeroes <- filter(detectandcensor, Sampler_Type %in% "Passive", Count_AboveDetection<1)

#=====================================================================================
##calculate means, medians, max, min, stddev, second highs, and full dataset percentiles for ACTIVES
#=====================================================================================
activesnozeroessummary <- group_by(activesnozeroes, Grouping_Factor) %>% mutate(Mean=mean(result, na.rm=T), Median=median(result, na.rm=T), Minimum=min(result, na.rm=T), Maximum=max(result, na.rm=T), StandardDeviation=sd(result, na.rm=T), Percentile95=quantile(result, probs=0.95, names=F, na.rm=T), Percentile5=quantile(result, probs=0.05, names=F, na.rm=T), secondhigh = sort(result, TRUE)[2])
  
activeszeroessummary <- group_by(activeszeroes, Grouping_Factor) %>% mutate(Mean=mean(ResultZero, na.rm=T), Median=median(ResultZero, na.rm=T), Minimum=min(ResultZero, na.rm=T), Maximum=max(ResultZero, na.rm=T), StandardDeviation=sd(ResultZero, na.rm=T), Percentile95=quantile(ResultZero, probs=0.95, names=F, na.rm=T), Percentile5=quantile(ResultZero, probs=0.05, names=F, na.rm=T), secondhigh = sort(ResultZero, TRUE)[2])

actives <- rbind(activesnozeroessummary, activeszeroessummary)
#write.csv(actives, "actives.csv")
actives$ResultZero <- NULL

#=======================================================================================
##calculate means, medians, max, min, stddev, second highs, and full dataset percentiles for PASSIVES
#=======================================================================================
passivesnozeroessummary <- group_by(passivesnozeroes, Grouping_Factor) %>% mutate(Mean=mean(result, na.rm=T), Median=median(result, na.rm=T), Minimum=min(result, na.rm=T), Maximum=max(result, na.rm=T), StandardDeviation=sd(result, na.rm=T), Percentile95=quantile(result, probs=0.95, names=F, na.rm=T), Percentile5=quantile(result, probs=0.05, names=F, na.rm=T), secondhigh = sort(result, TRUE)[2])

passiveszeroessummary <- group_by(passiveszeroes, Grouping_Factor) %>% mutate(Mean=mean(ResultZero, na.rm=T), Median=median(ResultZero, na.rm=T), Minimum=min(ResultZero, na.rm=T), Maximum=max(ResultZero, na.rm=T), StandardDeviation=sd(ResultZero, na.rm=T), Percentile95=quantile(ResultZero, probs=0.95, names=F, na.rm=T), Percentile5=quantile(ResultZero, probs=0.05, names=F, na.rm=T), secondhigh = sort(ResultZero, TRUE)[2])

passives <- rbind(passivesnozeroessummary, passiveszeroessummary)
#write.csv(passives, "passives.csv")
passives$ResultZero <- NULL

pahconcsummary <- rbind(passives, actives)

#================================================
##merge the summaries together
##Calculate coefficient of variation and range
#================================================
pahconcsummary$siteid <- as.character(pahconcsummary$siteid)
pahconcsummary <- unique(pahconcsummary)

pahconcsummary <- ungroup(pahconcsummary)
pahconcsummary<-mutate(pahconcsummary, "Coefficient of Variation" = pahconcsummary$StandardDeviation/pahconcsummary$Mean, "Range" = pahconcsummary$Maximum - pahconcsummary$Minimum)


##========================================================================
#Filter for ROS Mean and UCL Calculations
#ROS UCLs cannot be calcualted without at least 3 unique uncensored values
##========================================================================
ND_filter <- detectandcensor$uniQ < 3
countcens <- detectandcensor[!ND_filter, ]
rm(ND_filter)


#=================================================================================
##Calculate Kaplan Meier stats (no rank in censored results): The Kaplan-Meier method is a nonparametric technique for calculating the (cumulative) probability distribution and for estimating means, sums, and variances with censored data. Originally, the Kaplan-Meier approach was developed for right-censored survival data. More recently, the method was reformulated for left-censored environmental measurements (e.g., nondetects). USEPA's Unified Guidance also recommends the Kaplan-Meier method for use as an intermediate step in calculating parametric prediction limits, control charts, and confidence limits for censored data sets. In this latter application, the Kaplan-Meier estimate of the mean and variance is substituted for the sample mean and variance in the appropriate parametric formula.
#================================================================================

KMstatspahstudy<-cenfit(as.numeric(countcens$result),as.logical(countcens$censored), as.factor(countcens$Grouping_Factor),conf.int=0.95)
write.csv(show(KMstatspahstudy), file="KMstatspahstudy.csv")
KMstatspahstudy2<-read_csv("KMstatspahstudy.csv")
colnames(KMstatspahstudy2)<-c("Grouping_Factor",  "KM_n",  "KM_n.cen",  "KM_median",  "KM_mean",	"KM_sd")
KMstatspahstudy2$Grouping_Factor<-str_sub(KMstatspahstudy2$Grouping_Factor, start=30L)
KMstatspahstudy2$Grouping_Factor <- gsub("Factor)=", "", KMstatspahstudy2$Grouping_Factor)
#==============================================================================
#   Grab ROS means using EnvStats 
##http://www.itrcweb.org/gsmc-1/Content/GW%20Stats/5%20Methods%20in%20indiv%20Topics/5%207%20Nondetects.htm
#Robust regression on order statistics (ROS) is a semi-parametric method that can be used to estimate means and other statistics with censored data. Unlike Kaplan-Meier, ROS internally assumes that the underlying population is approximately normal or lognormal. However, the assumption is directly applied to only the censored measurements and not to the full data set (hence the term 'semi-parametric'). In particular, ROS plots the detected values on a probability plot (with a regular or log-transformed axis) and calculates a linear regression line in order to approximate the parameters of the underlying (assumed) distribution. This fitted distribution is then utilized to generate imputed estimates for each of the censored measurements, which are then combined with the known (detected) values to summary statistics of interest (e.g., mean, variance). The method is labeled 'robust' because the detected measurements are used 'as is' to make estimates, rather than simply using the fitted distributional parameters from the probability plot.
#enormCensored will remove missing values, NA and NaN without a na.rm statement.
#==============================================================================

#==============================================================================
# Suppress repeated warnings
options(warn = -1)
#==============================================================================

countcens$result <- signif(countcens$result, 4)
pahros <- group_by(countcens, Grouping_Factor) %>% 
  mutate(ROSMean= ifelse(Count_Censored[1] < 1, mean(result, na.rm=T), enormCensored(result, as.logical(censored), method="impute.w.qq.reg", ci=F, lb.impute=min(MDL)*.75)$parameters[[1]]), TYPEOFROSMEAN = ifelse(Count_Censored[1] < 1, "ArithmeticMean", "ROSMEAN"))

######Bootstrapping UCL for PAHs considering Non-detects#######
######################################
# 
#   Bootstrap ROS UCL95
#
######################################
#source(bootPAHs.R)
# Define boot_ROS75
boot_ROS75 <- function(group      = Grouping_Factor,     
                       results    = result, 
                       censored   = censored,
                       mdl        = min(MDL),
                       repeats    = 2000,
                       seed       = 27){
  
  options(warn=-1)
  #Set seed for random number generator
  set.seed(seed)
  
  print(group)
  # Set censored as logical (True/False)
  censored <- as.logical(censored)
  
  # Get length of results
  N <- length(censored)
  
  # Define function to record ROS mean of re-sampled table
  getROS <- function(n=N){
    random.rows <- sample(1:n, replace=T)
    #print(results[random.rows][!censored[random.rows]])
    while ((sum(censored[random.rows], na.rm=T) > n-2) | (length(unique(results[random.rows][!censored[random.rows]])) < 3)) {
      random.rows <- sample(1:n, replace=T) 
      }
    
    if (sum(censored[random.rows]) < 1) return(mean(results[random.rows], na.rm=T))
    
    enormCensored(results[random.rows], censored[random.rows], method="impute.w.qq.reg", ci=F, lb.impute=mdl*.75)$parameters[[1]]
  }
  # Use lapply to repeat EnvStats ROS mean 
  bootedMeans <- lapply(1:repeats, FUN=function(x) getROS())
  boot_UCL <- sort(unlist(bootedMeans))[repeats*.95 +1]
  options(warn=1)
  
  return(boot_UCL[[1]]) 
}
#================
##Grab Values
#UCL95 is booted ROS if Count Censored is more than 1, if less it is a boot strapped AR mean
#================
set.seed(seed)

pahros$result <- signif(pahros$result, 4)
pahros <- group_by(pahros, Grouping_Factor) %>% mutate(UCL95=if (Count_Censored[1] < 1) sort(replicate(Boot_Repeats, mean(sample(result, replace=T), na.rm=T)))[Boot_Repeats*.95+1]
      else boot_ROS75(Grouping_Factor[1], result, censored, min(MDL), Boot_Repeats, seed=seed))

#saveRDS(pahros, "sample_data.rdata")

#pahros <- readRDS("sample_data.rdata")

pahros <- group_by(pahros, Grouping_Factor) %>% mutate(SubDL_UCL95= sort(replicate(Boot_Repeats, mean(sample(result, replace=T), na.rm=T)))[Boot_Repeats*.95+1])


nocounts <- names(pahros) %in% c("Count_Censored", "Count_Sampled", "Count_NonCensored", "Pct_AboveDetection", "Pct_Censored", "CensoredDummy", "censored", "result", "resultcen", "MDL") 
pahros <- pahros[, !nocounts]
pahros <- unique(pahros)

#========================================================
#Add Kaplan Meier Statistics to data frame for comparison
#========================================================

addKM<-merge(pahconcsummary, KMstatspahstudy2, by="Grouping_Factor", all.x=T, sort=FALSE)

bighugesummarystudy<-merge(addKM, pahros, by=c("Grouping_Factor", "Sampler_Type", "Analyte", "CAS", "unit", "siteid", "Season", "Study_Year", "Season_StudyYear", "Time_Summary", "uniQ", "Count_AboveDetection"), all.x=T, all.y=T)    
addKM <- unique(addKM)

##===============================================================================
##Provide a column to identify which mean and UCL calculation is representated
##Create a "Reported Mean" field in addition to the arithmetic mean field for comparison
##===============================================================================
bighugesummarystudy <- mutate(bighugesummarystudy, ReportedMean_Type = ifelse(Count_Sampled==1, "Only One Sample in Group", ifelse(Count_Censored[1] < 1, "Arithmetic Mean", "ROS Mean")))

bighugesummarystudy <- mutate(bighugesummarystudy, UCL95Type = ifelse(is.na(UCL95), "Too Few Data for UCl Calculation", ifelse(Count_Censored[1] < 1, "UCL95_Students-T", "Booted ROS UCL95")))

bighugesummarystudy <- ungroup(bighugesummarystudy)
bighugesummarystudy <- mutate(bighugesummarystudy, ReportedMean = ifelse(Count_Censored < 1, Mean, ROSMean))

bighugesummarystudy <- bighugesummarystudy[, c("Grouping_Factor", "Sampler_Type", "Analyte", "CAS", "unit", "siteid", "Season", "Study_Year", "Season_StudyYear", "Time_Summary", "uniQ", "Count_AboveDetection", "Mean", "Median", "Minimum", "Maximum", "StandardDeviation", "Percentile95", "Percentile5", "secondhigh", "Count_Sampled", "Count_Censored", "Pct_AboveDetection", "Pct_Censored", "Coefficient of Variation", "Range", "ROSMean", "UCL95", "UCL95Type", "ReportedMean_Type", "ReportedMean")]

bighugesummarystudy <- unique(bighugesummarystudy)

#==========================================
##add site and cancer potency  information
##=========================================
db="X:/Programs/Air_Quality_Programs/Air Focus Areas/EPA Community Scale Monitoring of PAHs Project/Data Management/MPCA_MDH_DNRE Data/PAHs_Air_EPA_Grant.accdb"
bigdb=odbcConnectAccess2007(db)
pahinfo=sqlFetch(bigdb, "PAH Information")
pahinfo=data.frame(pahinfo, stringsAsFactors=F)
siteinfo=sqlFetch(bigdb, "Site Information")
siteinfo=data.frame(siteinfo, stringsAsFactors=F)
##Close connection with database
odbcClose(bigdb)

pahinfo <- pahinfo[!pahinfo$Chemical.Name %in% c("Benzo[a]pyrene-d12", "Pyrene-d10", "Fluorene-d10", "Fluoranthene-d10"), c("Chemical.Name", "CAS", "Molecular.Weight", "PAH.Unit.Risk", "MDH.RPF", "Phenyl.Rings", "Acronym", "PAH.16.or.7", "Nitro..Methyl..Thio.PAH")]
pahinfo$Chemical.Name <- as.character(pahinfo$Chemical.Name)
pahinfo$CAS <- as.character(pahinfo$CAS)

siteinfo <- siteinfo[, c("AQS.Code", "MPCA.Codes", "Site.Name", "City", "County", "Zipcode", "lat", "long", "Census.Tract", "Census.Block", "Location.Setting", "sensitive.populations", "potential.sources")]
siteinfo <- unique(siteinfo)

bighugesummarystudy <- merge(bighugesummarystudy, pahinfo, by.x = c("Analyte", "CAS"), by.y=c("Chemical.Name", "CAS"), all=FALSE)

bighugesummarystudy <- merge(bighugesummarystudy, siteinfo, by.x="siteid", by.y = "MPCA.Codes", all=TRUE)

write_csv(bighugesummarystudy, "PAHs Time Based Summaries.csv")
