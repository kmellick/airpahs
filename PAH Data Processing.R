##Find Libraries
library(dplyr)
library(RODBC)
library(stringr)
library(readr)
#library(rmarkdown)

##Set working Directory
setwd("X:/Programs/Air_Quality_Programs/Air Focus Areas/EPA Community Scale Monitoring of PAHs Project/Data Analysis/Data Processing") 

###################################################################
##Pull In Data. Data are from an Access Table Named All Results. This query relates the concentration sheets provided by MDH, MPCA field sheets, a Site Information Table and a time summarizing sheet. 
###################################################################
conn = "X:/Programs/Air_Quality_Programs/Air Focus Areas/EPA Community Scale Monitoring of PAHs Project/Data Management/MPCA_MDH_DNRE Data/PAHs_Air_EPA_Grant.accdb"
conn = odbcConnectAccess2007(conn) ##Connect to the ACCESS database
alldata = sqlFetch(conn, "All Results", stringsAsFactors=FALSE) # read a table
qcdata = sqlFetch(conn, "QC Results", stringsAsFactors=FALSE)
siteinfo = sqlFetch(conn, "Site Information", stringsAsFactors=FALSE)
pahinfo = sqlFetch(conn, "PAH Information", stringsAsFactors=FALSE)
odbcClose(conn) # close the connection to the file

alldata <- unique(alldata)
alldata=as.data.frame(alldata, stringsAsFactors=FALSE)

names(alldata) <- str_replace_all(str_trim(names(alldata)), " ", "_")

qcdata <- as.data.frame(qcdata)
qcdata$ID <- NULL
qcdata <- filter(qcdata, QCCode %in% c("LCS", "LCS Dup"))

##Some analyses were run twice to confirm QC results. These QC results are averaged below. Then the QC table is used to create a validation table based on analyte, batch and sampler_type.
qcdata <- group_by(qcdata, LabNumber, Batch, QCCode, Analysis, Analyte, Units, Dilution, MDL, MRL, Recovery, UCL, LCL, RPD, 'RPD Limit', Analyzed, Qualifier, 'Qualifier Text') %>% summarise(Result=mean(Result, na.rm=T)) %>% ungroup()



#####################################################
##Laboratory Control Standards are acceptable if within 30 - 150% recovery. They are rejected if they are below 10% recover, and are considered "estimated" if they are between 10% and 30$ recovery. No samples were above 150% recovery.The Recovery column is from the qcdata table, not the % Recovery column from the concentration table.
#####################################################
qcdata_batch <- mutate(qcdata, LCS_Validation = ifelse(Recovery>=10 & Recovery<30, "Estimated", ifelse(Recovery<10, "Reject", "Valid"))) %>% ungroup()

qcdata_batch$LCS_Validation_Reason <- ifelse(qcdata_batch$LCS_Validation=="Reject", "LCS Recovery less than 10%", ifelse(qcdata_batch$LCS_Validation=="Estimated", "LCS Recovery between 10% and 30%", "LCS Recovery within Acceptance Criteria"))

#####################################################
##Laboratory Control Standard Duplicates are accetable if relative percent difference is above 50%within 30 - 150% recovery. They are rejected if they are below 10% recover, and are considered "estimated" if they are between 10% and 30$ recovery. No samples were above 150% recovery.
#####################################################
qcdata_batch <- mutate(qcdata_batch, RPD_Validation = ifelse(RPD>100, "Reject", ifelse(RPD>50 & RPD<=100, "Estimated", "Valid")))

qcdata_batch$RPD_Validation_Reason <- ifelse(qcdata_batch$RPD_Validation=="Reject", "Relative Percent Different of LCS Duplicate is Above 100%", ifelse(qcdata_batch$RPD_Validation=="Estimated", "Relative Percent Different of LCS Duplicate is Above 50% and Below 100%", "Relative Percent Difference of LCS Duplicate is within Acceptance Criteria"))

qcdata_batch <- qcdata_batch[!is.na(qcdata_batch$RPD_Validation), c("Batch", "Analysis", "Analyte", "LCS_Validation", "LCS_Validation_Reason", "RPD_Validation", "RPD_Validation_Reason")]
qcdata_batch <- unique(qcdata_batch)

#####################################################
##Some concentrations have 2 qualifiers, so this line concatenates the two or more qualifiers and comma separates them.
#####################################################
alldata <- group_by(alldata, MDH_Lab_ID_Code, Sampler_Type, MPCA_Site_ID, Units, CAS, Start_Run_Date, Study_Location, passive_active, Primary_Occurence_Code, location_if_collocated, Year, Season, Flags_Text, Study_Year, Analyte, Batch) %>% summarise(Result=mean(Result, na.rm=T), MDL = mean(MDL, na.rm=T), Percent_recover = mean('%_Recovery', na.rm=T),Qualifier= paste(Qualifier, collapse=",")) %>% ungroup()

#####################################################
##Will qualify samples where the surrogate recovery was below the acceptable level, matrix affects or degradation is possible in these samples. Deuterated surrogate recoveries are acceptable if they are within 30 - 150% recovery. The Benzo[a]pyrene surrogate was the only aurrogate measured below acceptable limits.
#####################################################
alldata$Surrogate_Validation <- ifelse(grepl("(S2)|(S3)", alldata$Qualifier), "Estimated", "Valid")

alldata$Surrogate_Validation_Reason <- ifelse(alldata$Surrogate_Validation=="Estimated", "Deuterated Surrogate Recovery was less than 30%", "Deuterated Surrogate Recovery was within acceptable limits")


###################################################################
##Pull In QC Data. Join by Analyte, Sampler_Type, Batch
##Create Validation Column and a Validation Reason
###################################################################
data_validate <- merge(alldata, qcdata_batch, by.x=c("Analyte", "Batch", "Sampler_Type"), by.y=c("Analyte", "Batch", "Analysis"), all=T)

#####################################################
##Not all Batches have LCS runs. This allows grouping by Batch month as well as Batch
#####################################################
data_validate$Batch_Month <- substr(data_validate$Batch, 1, 3)

data_validate_noNA <- filter(data_validate, !is.na(LCS_Validation) | !is.na(RPD_Validation))
data_validate_NA <- filter(data_validate, is.na(LCS_Validation) | is.na(RPD_Validation))

data_validate_NA <- mutate(data_validate_NA, LCS_Validation = "No information to Validate", LCS_Validation_Reason= "No information to Validate", RPD_Validation= "No information to Validate", RPD_Validation_Reason= "No information to Validate")
  
data_validate_final <- rbind(data_validate_NA, data_validate_noNA)
  
data_validate_final <- mutate(data_validate_final, Data_Validation_ID=paste0(Surrogate_Validation, "_", LCS_Validation, "_", RPD_Validation), Data_Validation_Reason=paste0(Surrogate_Validation_Reason, "_", LCS_Validation_Reason, "_", RPD_Validation_Reason))

data_validate_final <- filter(data_validate_final, !grepl("Reject", Data_Validation_ID))

data_validate_final <- data_validate_final[, c("Analyte", "Sampler_Type", "MDH_Lab_ID_Code", "MPCA_Site_ID", "Units", "CAS", "Start_Run_Date", "Study_Location", "passive_active", "Primary_Occurence_Code", "location_if_collocated", "Year", "Season", "Flags_Text", "Study_Year", "Result", "MDL", "Qualifier", "Data_Validation_ID", "Data_Validation_Reason")]

alldata <- data_validate_final
alldata <- unique(alldata)
#write.csv(alldata, "alldata.csv")
#####################################################
##Remove lab control samples from the data set in the following line
#####################################################
alldata=subset(alldata,!(alldata$Analyte%in%c("Benzo[a]pyrene-d12","Fluoranthene-d10","Fluorene-d10","Pyrene-d10", "Benzo[e]pyrene-d12", "Benzo[g,h,i]perylene-d12")))

#############################################################
##DBTSulfone does not have an MDL for some of the non-reported values. The MDL is constant for passive samples and active samples.In order to make the data set consistent, these missing MDLs will be replaced with the maximum MDL within sample type and analyte. This only applies to DBTS in active samples. This will be cHecked once the data set is final.
#############################################################
minimumDBTS <- filter(alldata, Analyte=="Dibenzothiophene sulfone", Result>0, !is.na(Result), Sampler_Type %in% "PAHs Air, HiVol XAD")
minimumDBTS <- max(minimumDBTS$MDL)
alldata$MDL <- ifelse(is.na(alldata$MDL), minimumDBTS, alldata$MDL)

#############################################################
#Convert NA's to zero, add ND to qualifier string, add Censored field****
#############################################################
alldata$Qualifier=ifelse(is.na(alldata$Result),paste("ND",alldata$Qualifier,sep=","),paste(alldata$Qualifier))
alldata$Result=ifelse(is.na(alldata$Result),0,alldata$Result)

#############################################################
##First we calculate the Sampling Rates of the passive monitors compared to the active monitors. Sampling rates have the units m^3/day. They are a calibration factor for estimating air concentrations for the passive samplers.
##Pull Out Passives Only with no Field Blanks
#############################################################
alldata$passive_active <- as.character(alldata$passive_active)
passives2 <- filter(alldata, passive_active %in% "passive", MPCA_Site_ID %in% c("962", "963", "3051")) 
passives2$CenDummy=ifelse(passives2$Result<passives2$MDL | passives2$Result==0,1,0)
passives <- filter(alldata, passive_active %in% "passive") 

##Calculate a Season-Sampling Location-Analyte Specific Mean and include in all data
passivecalcs<- group_by(passives2, CAS, Study_Location, Season, Study_Year, MPCA_Site_ID) %>% summarise(passive_cal_mean=mean(Result, na.rm=T), CenDummy = sum(CenDummy), Count=n()) %>% ungroup()

##Pull IN ALL ACTIVES
##Remove the data for the two events that caused high total values (mainly particulate levels). These were the Cedar Riverside apartment fire and Mille Lacs Band of Ojibwe Pop Wow.
TOTALactives<- filter(alldata, !Sampler_Type %in% "PAHs Air, Passive XAD")

TOTALactives <- filter(TOTALactives, !(as.character(Start_Run_Date)=="2013-12-30" & as.character(TOTALactives$MPCA_Site_ID)=="963"))
TOTALactives <- filter(TOTALactives, !(as.character(Start_Run_Date)=="2014-08-15" & as.character(TOTALactives$MPCA_Site_ID)=="3051"))

TOTALactives$CenDummy=ifelse(TOTALactives$Result<TOTALactives$MDL | TOTALactives$Result==0,1,0)

##Sum gas and particle 
TOTALactives<- group_by(TOTALactives, Analyte, MDH_Lab_ID_Code, MPCA_Site_ID, Units, CAS, Start_Run_Date, Study_Location, Study_Year, Primary_Occurence_Code, location_if_collocated, Year, Season) %>% summarize(Result=sum(Result, na.rm=T), MDL = max(MDL, na.rm=T), Qualifier = paste0(Qualifier, collapse=","), CenDummy = ifelse(sum(CenDummy)==2|sum(Result, na.rm=T)==0, 1, 0), Count=n()) %>% ungroup()

##Eliminate samples without paired gas and particle
TOTALactives <- filter(TOTALactives, Count!=1)

##Convert Counts to 1 to align with the CenDummy Counts
TOTALactives$Count <- 1

##Average POCS
TOTALactives<- group_by(TOTALactives, Analyte, MPCA_Site_ID, Units, CAS, Start_Run_Date, Study_Location, Study_Year, Year, Season) %>% summarise(Result=mean(Result, na.rm=T), MDL = mean(MDL, na.rm=T), CenDummy = sum(CenDummy), Count=sum(Count)) %>% ungroup()

##Calculate the mean of the total actives by Season, Analyte, Study_Year, and Sampling Location
TOTALactives<- group_by(TOTALactives, CAS, Study_Location, Season, Study_Year, MPCA_Site_ID) %>% summarise(totalactive_cal_mean=mean(Result, na.rm=T), CenDummy = sum(CenDummy), Count=sum(Count)) %>% ungroup()

##Pull gasactive and passive means together into one data frame
passivecalcs <- merge(passivecalcs, TOTALactives, by=c("CAS", "Study_Location", "Study_Year" , "Season", "MPCA_Site_ID"), all=T)

##This is a place holder for the time when passives analysis is behind active reporting
passivecalcs$passive_cal_mean=ifelse(is.na(passivecalcs$passive_cal_mean),0,passivecalcs$passive_cal_mean)
passivecalcs$totalactive_cal_mean=ifelse(is.na(passivecalcs$totalactive_cal_mean),0,passivecalcs$totalactive_cal_mean)

##Calculate the sampling rate by Season, CAS, and Sampling Location
passivecalcs <- group_by(passivecalcs,  CAS, Study_Location, Season, Study_Year, MPCA_Site_ID) %>% mutate(samplingrate = ifelse(passive_cal_mean>0 & totalactive_cal_mean>0, passive_cal_mean/(totalactive_cal_mean*90), 0), CenDummy = sum(CenDummy.x, CenDummy.y), Count = sum(Count.x, Count.y), CenDummy.x=NULL, CenDummy.y=NULL, Count.x=NULL, Count.y=NULL) %>% ungroup()

##Average by Study Locations, taking out individual sites
passivecalcs <- group_by(passivecalcs,  CAS, Study_Location, Study_Year, Season) %>% summarize(samplingrate = mean(samplingrate, na.rm=T), percensored = (sum(CenDummy, na.rm=T)/sum(Count, na.rm=T))) %>% ungroup()

## Merge with all the other passive results
passives <- merge(passives, passivecalcs, by=c("CAS", "Study_Location", "Study_Year", "Season"), all=T) %>% ungroup()

passives$mass_loading_ng <- passives$Result

passives <- mutate(passives, Result = ifelse(samplingrate > 0, Result/(samplingrate*90), 0.00), MDL =  ifelse(samplingrate > 0, MDL/(samplingrate*90), 0.00), mass_loading_ng=Result) 

#check this equation
passives <- filter(passives, !is.na(Analyte))
#write.csv(passives, file="samplingrates.csv", row.names=TRUE) #for checking calculations

##Changing to air concentration units
passives$Units <- "ng/m3"

##Next create a set of data where gas phase active results are summed with particulate phase active results to calculated a total PAH concentration
actives<- filter(alldata, passive_active %in% "active") 


##Perform Censor Test on Data prior to summation for totals
passives$CenTest=ifelse(passives$Result<passives$MDL | passives$Result==0,"Yes","No")
actives$CenTest=ifelse(actives$Result<actives$MDL | actives$Result==0,"Yes","No")

##Sum gas and particle results
totals<- group_by(actives, MDH_Lab_ID_Code, MPCA_Site_ID, Units, CAS, Start_Run_Date, Study_Location, Study_Location, passive_active, Primary_Occurence_Code, location_if_collocated, Year, Season, Study_Year, Analyte) %>% summarise(Result=sum(Result, na.rm=T), Sampler_Type="Total", MDL=min(MDL, na.rm=T), Qualifier=paste(Qualifier, collapse=","), Flags_Text=paste(Flags_Text, collapse=","), CenTest=paste(CenTest, collapse=","), Data_Validation_ID = paste(Data_Validation_ID, collapse=","), Count = n(), Data_Validation_Reason = paste(Data_Validation_Reason, collapse=","))  

totals <- filter(totals, Count!=1)
##Add Nulls for the passive loadings and sampling rates for the active results
totals$mass_loading_ng <- "NA"
actives$mass_loading_ng <- "NA"
totals$samplingrate <- "NA"
actives$samplingrate <- "NA"
totals$percensored <- "NA"
actives$percensored <- "NA"
totals$Count <- NULL

##Pull all the concentrations together
allconcentrations <- rbind(passives, actives, totals)
##write.csv(allconcentrations, "allconcentrations.csv")   ##check merge if anything is changed

#******Using new censor indicator which is generated individually from passive, active, and total samples 
allconcentrations$Censored <- ifelse(allconcentrations$CenTest=="Yes" | allconcentrations$CenTest=="Yes,Yes",T,F)

##Add a season_study year factor
allconcentrations <- mutate(allconcentrations, Season_StudyYear=paste(Season, Study_Year))

allconcentrations$Sampler_Type <- ifelse(allconcentrations$Sampler_Type %in% "PAHs Air, Passive XAD", "Passive", ifelse(allconcentrations$Sampler_Type %in% "PAHs Air, HiVol XAD", "Gas", ifelse(allconcentrations$Sampler_Type %in% "PAHs Air, HiVol Quartz", "Particle", "Total")))

allconcentrations$location_if_collocated <- as.character(allconcentrations$location_if_collocated)
siteinfo$`location if collocated` <- as.character(siteinfo$`location if collocated`)

siteinfo <- siteinfo[, c("MPCA Codes", "Site Name", "Long Site Name", "City", "location if collocated")]
allconcentrations <- merge(allconcentrations, siteinfo, by.x=c("location_if_collocated", "MPCA_Site_ID"), by.y=c("location if collocated", "MPCA Codes"), all.x=TRUE, all.y=FALSE)
allconcentrations <- unique(allconcentrations)
##write.csv(allconcentrations, "allconcentrations.csv")   ##check merge if anything is changed

##Create a Unified Date Field so that all passive dates are aligned
allconcentrations$Start_Run_Date <- as.POSIXct(strptime(allconcentrations$Start_Run_Date, format="%Y-%m-%d"), tz="GMT")     

allconcentrations <- mutate(allconcentrations, FirstoftheMonth = format(Start_Run_Date, "%Y-%m-01"))
allconcentrations$Start_Run_Date <- as.character(allconcentrations$Start_Run_Date)
allconcentrations <- mutate(allconcentrations, Unified_Date = ifelse(Sampler_Type %in% "Passive", FirstoftheMonth, Start_Run_Date))
allconcentrations$FirstoftheMonth <- NULL

allconcentrations$Unified_Date <- as.POSIXct(strptime(allconcentrations$Unified_Date, format="%Y-%m-%d"), tz="GMT")
allconcentrations$Start_Run_Date <- as.POSIXct(strptime(allconcentrations$Start_Run_Date, format="%Y-%m-%d"), tz="GMT")

pahinfo <- as.data.frame(pahinfo)
pahinfo$ID <- NULL
pahinfo <- pahinfo[, c("CAS", "Chemical Name", "Molecular Weight", "Acronym")]
pahinfo$CAS <- paste0("a", pahinfo$CAS)
names(pahinfo) <- str_replace_all(str_trim(names(pahinfo)), " ", "_")

allconcentrations <- merge(allconcentrations, pahinfo, by.x=c("CAS", "Analyte"), by.y=c("CAS", "Chemical_Name"), all.x=TRUE)
##write_csv(allconcentrations, "allconcentrations.csv") ##check if anything is changed

#allconcentrations$CAS <- gsub("a","", allconcentrations$CAS)
write_csv(allconcentrations, "allPAHconcentrations.csv")

#allconcentrations$CAS <- gsub("a","", allconcentrations$CAS)
#write_csv(allconcentrations, "M:/GregsAnalysisLocation/allPAHconcentrations.csv")

##rmarkdown::render("X:\\Programs\\Air_Quality_Programs\\Air Focus Areas\\EPA Community Scale Monitoring of PAHs Project\\Data Analysis\\Data Processing\\PAH Data Processing.R", output_format = "Rmd_Document")
