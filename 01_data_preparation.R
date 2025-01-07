#Get Ready -----
setwd("C:/Users/User/Desktop/naiades") #home
setwd("//net1.cec.eu.int/jrc-services/IPR-Users/DOUWEHA/Desktop/naiades") #work

installed.packages() #check and install any missing packages
install.packages("callr")
install.packages("remotes")
remotes::install_github("hannaplantain/ipchem-dataretrieval-r-main", force = TRUE)

1library(dplyr)
library(tidyverse)
library(readxl)
library(lubridate)
library(callr)
library(ipchem.dataretrieval)
library(data.table) 
# Chemical data basic prep ----

#Chemical selection
load_summary() # get chem summary from ipchem

#get tox data and priority list - both from Posthuma
Posthuma_tox <- read_excel("Posthuma_List_Filtered.xlsx")
chem_list <- Posthuma_tox$CAS

Posthuma_prio <- read_excel("Posthuma_2020_supptable.xlsx")
prio_chem_list <- Posthuma_prio$CAS

prio_posthuma_tox <- Posthuma_tox %>% #tox data for priority chemicals only
  filter(CAS %in% prio_chem_list)

prio_chem_list <- prio_posthuma_tox$CAS #prio chem list with tox data available

rm(Posthuma_prio)
rm(Posthuma_tox)
rm(chem_list)

#extract relevant chemical data Naiades from IPCHEM database
retrieve_meta("NAIADES")

NAIADES_POST <- chem_summary %>% 
  filter(cas_number %in% prio_chem_list, #change to chem_list for real anals
         dataset == "NAIADES", 
         period %in% 2013:2022) #change period or remove for real analysis

NAIADES_DATA <- NAIADES_POST %>% 
  #retrieve_media(filter = (stringr::str_detect(media, "water"))) %>%
  #retrieve_period(filter=(period= (years))) %>% 
  #retrieve_cas_number(cas_vector_formatted) %>% #select chemicals here
  create_package_url() %>% 
  #head(10) %>% #remove for real analysis, replaced with selected chems
  download_files(exp_jsn = FALSE, json_col = "SOURCE DATA")

# Bind rows
NAIADES_DATA <- lapply(NAIADES_DATA, function(df) {
  df$`Sampling Date` <- as.POSIXct(df$`Sampling Date`, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  return(df)
})

CHEMDF <- dplyr::bind_rows(NAIADES_DATA)

unique(CHEMDF$`Chemical Name`) #have a firstlook
colnames(CHEMDF)
summary(CHEMDF) 

#not same medium/concentration
CHEMDF_sediment <- CHEMDF %>%
  filter(Media == "Water (Sediment in Water)") 
CHEMDF <- CHEMDF %>%
  filter(Media == "Water (Surface Water)") #choosing this one 10x more values

rm(NAIADES_POST)
rm(NAIADES_DATA)

#first cleanup chemical data

#get rid of spaces in colnames
colnames(CHEMDF) <- gsub(" ", "_", colnames(CHEMDF), fixed = TRUE) 

#fix date format
str(CHEMDF$Sampling_Date) 
CHEMDF$Sampling_Date <- as.character(CHEMDF$Sampling_Date)
CHEMDF <- CHEMDF %>% filter(!is.na(Sampling_Date))
CHEMDF$Sampling_Date <- ymd_hms(CHEMDF$Sampling_Date)

CHEMDF <- CHEMDF %>%
  mutate(
    Year = year(Sampling_Date),         
    Month = month(Sampling_Date),     
    Day = day(Sampling_Date),           
    Time = format(Sampling_Date, "%H:%M:%S")  #Extract t in HH:MM:SS 
  )

CHEMDF <- CHEMDF %>%
  mutate(Sampling_Date = as.Date(paste(Year, Month, Day, sep = "-")))

#still to-do: find efficient way to check date formats to include dates that failed to parse

#export CHEMDF to csv, if R has to be terminated, get CHEMDF from here) -----
file_CHEMDF <- file("CHEMDF.csv","r")
chunk_size <- 100000 # choose the best size for you
line_offset <- 0
chunks <- list()
progress_interval <- 500000
total_rows <- 
  
repeat {
  chunk <- fread("CHEMDF.csv", skip = line_offset, nrows = chunk_size)
  if (nrow(chunk) == 0) break
  chunks[[length(chunks) + 1]] <- chunk
  line_offset <- line_offset + nrow(chunk)
  if (line_offset %% progress_interval == 0) {
    cat("Processed", line_offset, "rows\n")
  }
}

CHEMDF <- rbindlist(chunks)

#Import biologicaal data ----
BIO13 <- read.csv("I2M2_13/ResultatsBiologiques.csv", sep=";")
BIO14 <- read.csv("I2M2_14/ResultatsBiologiques.csv", sep=";")
BIO15 <- read.csv("I2M2_15/ResultatsBiologiques.csv", sep=";")
BIO16 <- read.csv("I2M2_16/ResultatsBiologiques.csv", sep=";")
BIO17 <- read.csv("I2M2_17/ResultatsBiologiques.csv", sep=";")
BIO18 <- read.csv("I2M2_18/ResultatsBiologiques.csv", sep=";")
BIO19 <- read.csv("I2M2_19/ResultatsBiologiques.csv", sep=";")
BIO20 <- read.csv("I2M2_20/ResultatsBiologiques.csv", sep=";")
BIO21 <- read.csv("I2M2_21/ResultatsBiologiques.csv", sep=";")
BIO22 <- read.csv("I2M2_22/ResultatsBiologiques.csv", sep=";")

BIO <- dplyr::bind_rows(BIO13, BIO14, BIO15, BIO16, BIO17, BIO18, BIO19, BIO20, BIO21, BIO22)
rm(BIO13, BIO14, BIO15, BIO16, BIO17, BIO18, BIO19, BIO20, BIO21, BIO22)

colnames(BIO)[which(names(BIO) == "CdStationMesureEauxSurface")] <- "Sample_Source_Code"
colnames(BIO)[which(names(BIO) == "LbStationMesureEauxSurface")] <- "Sample_Source_Name"
colnames(BIO)[which(names(BIO) == "ResIndiceResultatBiologique")] <- "I2M2"
colnames(BIO)[which(names(BIO) == "DateDebutOperationPrelBio")] <- "Bio_Sampling_Date"

str(BIO$Bio_Sampling_Date)
BIO$Bio_Sampling_Date <- as.Date(BIO$Bio_Sampling_Date)

BIO <- BIO %>%
  mutate(
    Year = year(Bio_Sampling_Date))  

#fwrite(BIO, "BIO.csv") #export BIO to csv, if R has to be terminated, get BIO data here)
#BIO <- read.csv("BIO.csv")

#matching - before continuing, keep data which contain same Sample_Source_Code ----
#remove leading 0 IPCHEM codes to make Sample_Source_Code match
CHEMDF <- CHEMDF %>%
  mutate(Sample_Source_Code = sub("^.", "", as.character(Sample_Source_Code))) #make sure to run just once

#filter
common_station_codes <- intersect(BIO$Sample_Source_Code, CHEMDF$Sample_Source_Code)
BIO_filtered <- BIO[BIO$Sample_Source_Code %in% common_station_codes, ] 
CHEMDF_filtered <- CHEMDF[CHEMDF$Sample_Source_Code %in% common_station_codes, ]
rm(CHEMDF)

#randomly remove duplicate bio measurements -----
BIO_duplicates <- BIO_filtered %>% #shows how many and which measurements are duplicates (same station, same year)
  group_by(Sample_Source_Code, Year) %>%
  filter(n() > 1)  

BIO_filtered <- BIO_filtered %>% #average measurements in same location on same day
  group_by(Sample_Source_Code, Bio_Sampling_Date) %>%
  mutate(I2M2_is_average = ifelse(n() > 1, TRUE, FALSE)) %>%
  summarise(
    I2M2 = mean(I2M2, na.rm = TRUE),
    
    Sample_Source_Name = first(Sample_Source_Name),
    CdPointEauxSurf = first(CdPointEauxSurf),
    CdSupport = first(CdSupport),
    LbSupport = first(LbSupport),
    CdUniteMesure = first(CdUniteMesure),
    SymUniteMesure = first(SymUniteMesure),
    CdRqIndiceResultatBiologique = first(CdRqIndiceResultatBiologique),
    MnemoRqAna = first(MnemoRqAna),
    CdMethEval = first(CdMethEval),
    RefOperationPrelBio = first(RefOperationPrelBio),
    CdProducteur = first(CdProducteur),
    NomProducteur = first(NomProducteur),
    CdAccredRsIndiceResultatBiologique = first(CdAccredRsIndiceResultatBiologique),
    MnAccredRsIndiceResultatBiologique = first(MnAccredRsIndiceResultatBiologique),
    Date = first(Bio_Sampling_Date),
    Year = first(Year),
    I2M2_is_average = first(I2M2_is_average)  
  ) %>%
  ungroup()

BIO_duplicates <- BIO_filtered %>% 
  group_by(Sample_Source_Code, Year) %>%
  filter(n() > 1)  

# function to remove duplicates randomly---- 
# https://www.r-bloggers.com/2013/01/randomly-deleting-duplicate-rows-from-a-dataframe/
duplicated.random = function(x, y, incomparables = FALSE, ...) { 
  if ( is.vector(x) & is.vector(y) ) { 
    combined = data.frame(station = x, year = y)
    permutation = sample(nrow(combined)) 
    combined.perm = combined[permutation, ] 
    result.perm = duplicated(combined.perm, incomparables, ...) 
    result = result.perm[order(permutation)] 
    
    return(result) 
  } 
  else { 
    stop(paste("duplicated.random() requires both station and year as vectors.")) 
  } 
} 

BIO_filtered$duplicate <- duplicated.random(BIO_filtered$Sample_Source_Code, BIO_filtered$Year)

BIO_unique <- BIO_filtered %>%
  filter(duplicate == FALSE)

# Check for rows in BIO_unique that are duplicate many-to-many matches in CHEMDF_filtered
BIO_duplicates <- BIO_unique %>% #by code - should be 0 now
  group_by(Sample_Source_Code, Year) %>%
  tally() %>%
  filter(n > 1) 

# Select chem data based on bio data----
#add start and end date to BIO_filtered
BIO_unique <- BIO_unique %>%
  mutate(Start_Date = Bio_Sampling_Date - months(3), 
         End_Date = Bio_Sampling_Date + months(0)) 

BIO_unique <- BIO_unique %>%
  mutate(Sample_Source_Code = as.character(Sample_Source_Code))

#keep chem measurements within time window bio measurement
CHEMDF_filtered <- CHEMDF_filtered %>%
  inner_join(BIO_unique,
             by = c("Sample_Source_Code", "Year")) %>%
  mutate(
    Start_Date = as.Date(Start_Date),
    End_Date = as.Date(End_Date),
    Sampling_Date = as.Date(Sampling_Date)  
  ) %>%
  filter(Sampling_Date >= Start_Date & Sampling_Date <= End_Date)

# DON'T FOR NOW - kill measurements below LOQ ----
#colnames(CHEMDF_filtered)
#CHEMDF_LOQ <- CHEMDF_filtered %>%
  #filter(Concentration_Value > LOQ)

# reshape chemical data based on station, select cols for clarity----
CHEMDF_clean <- CHEMDF_filtered %>%
  select(Chemical_Name, CAS_Number, Sample_Source_Code:Media, Latitude, Longitude, Year:Day)

#Get one concentration/chemical/site
## option 1: select highest (when tied, randomly select one ) (lose many data points)
CHEMDF_w_max <- CHEMDF_clean %>%
  group_by(Sample_Source_Code, CAS_Number) %>%  
  arrange(Sample_Source_Code, CAS_Number, desc(Concentration_Value)) %>%  
  slice(1) %>%
  ungroup()

# prep for chemical combination selection & TU conversion -----

#pivot chemical data to wide format 
colnames(CHEMDF_w_max)[which(names(CHEMDF_w_max) == "Sample_Source_Name.x")] <- "Sample_Source_Name"

CHEMDF_max_clean <- CHEMDF_w_max %>%
  select(CAS_Number, Sample_Source_Name, Sample_Source_Code, Year, Concentration_Value)

CHEMDF_duplicates <- CHEMDF_max_clean %>% #check code - should be 0 now
  group_by(CAS_Number, Sample_Source_Code, Year) %>%
  tally() %>%
  filter(n > 1) 

head(CHEMDF_max_clean) #check

CHEMDF_wide <- CHEMDF_max_clean %>%
  pivot_wider(
    names_from = c(CAS_Number),  
    values_from = Concentration_Value,  
    values_fill = list(Concentration_Value = NA)  
  )

head(CHEMDF_wide) #check

#merge dataframes
DF <- CHEMDF_wide %>%
  left_join(BIO_unique, by = c("Sample_Source_Code", "Year"))

#END OF FILE go to 02_mixture_TU
