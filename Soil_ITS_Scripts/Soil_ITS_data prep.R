library(tidyverse)

#import Blast ID df and clean data
Blast_ID<-read_excel('Raw_data/Updated_Data/ABS.MER.fielddata.Feb.2023_Site.Info_AF.xlsx')
Blast_ID$Site= sub(c('ABS00|ABS0'),'',Blast_ID$Site)
Blast_ID$sample_ID<-as.character(Blast_ID$sample_ID)
Blast_ID$transect<-as.character(Blast_ID$transect)
Blast_ID$sample_ID<-sub(c("-"),".",Blast_ID$sample_ID)
Blast_ID$sample_ID<-sub(c("-"),".",Blast_ID$sample_ID)
names(Blast_ID)[9]<-'Veg_Class_Abv'

#Nute and Veg data
Nutrients_Transects<-read_excel('Processed_data/Nutrients_Transect_level.xlsx')
VEG_COVER_Transects <- read_excel("Raw_data/Site_Data/ABS.MER.fielddata.Feb.2023_R.PROCESSED.VEG.COVER_ALL.xlsx", 
                                  sheet = "Transect.Level_Data")
VEG_COVER_Transects$Site= sub(c('ABS00|ABS0'),'',VEG_COVER_Transects$Site)
VEG_COVER_Transects$Transect= sub(c('T'),'',VEG_COVER_Transects$Transect)
Myco_plant_spp<-read.csv('Processed_data/Myco_host_abundance.csv')%>%
  mutate(Site=as.factor(Site),
         Transect=as.factor(Transect))


Blast_ID<-Blast_ID%>%
  dplyr::mutate(Regime = paste(Interval, Severity, sep = "_"))%>%
  rename(Transect=transect)%>%
  select(-`Fire 3`,-`Interval (yrs)...16`,-`FESM severity category`)%>%
  left_join(Myco_plant_spp)%>%
  left_join(Nutrients_Transects)%>%
  left_join(VEG_COVER_Transects)%>%
  #this is selecting for the 30 decomp sites which I have CN data for all of
  filter(!is.na(Carbon))

summary(Blast_ID)#no na;s in df



#loading funguild data
# Load required libraries
library(httr)          # For working with HTTP requests
library(jsonlite)      # For parsing JSON data
library(lubridate)     # For handling dates
library(XML)           # For parsing HTML/XML data

# Function to parse data from FunGuild database
parse_funguild <- function(url = 'http://www.stbates.org/funguild_db.php', tax_name = TRUE) {
  # Parse the HTML content of the webpage
  tmp <- XML::htmlParse(url)
  tmp <- XML::xpathSApply(doc = tmp, path = "//body", fun = XML::xmlValue)
  
  # Convert JSON data from the webpage to a data frame
  db <- jsonlite::fromJSON(txt=tmp)
  
  # Remove unnecessary ID column
  db$`_id` <- NULL
  
  if (tax_name == TRUE) {
    # Define taxonomic level names corresponding to numerical codes
    taxons <- c(
      "keyword",                                                       # 0
      "Phylum", "Subphylum", "Class", "Subclass", "Order", "Suborder", # 3:8
      "Family", "Subfamily", "Tribe", "Subtribe", "Genus",             # 9:13
      "Subgenus", "Section", "Subsection", "Series", "Subseries",      # 15:19
      "Species", "Subspecies", "Variety", "Subvariety", "Form",        # 20:24
      "Subform", "Form Species")
    
    # Create a lookup table for taxonomic levels
    taxmatch <- data.frame(
      TaxID = c(0, 3:13, 15:26),
      Taxon = factor(taxons, levels = taxons)
    )
    
    # Replace numerical taxonomic levels with their corresponding names
    db$taxonomicLevel <- taxmatch[match(db$taxonomicLevel, taxmatch$TaxID), "Taxon"]
  }
  
  # Assign the download date as an attribute to the dataset
  attr(db, "DownloadDate") <- date()
  
  return(db)
}

# Define FunGuild database URL
url  <- "http://www.stbates.org/funguild_db.php"

# Fetch and parse FunGuild data
#FunGuildData_org <- parse_funguild(url)
#write.csv(FunGuildData_org,'Funguild/All_Funguild_10.3.25.csv', row.names = FALSE)
FunGuildData_org<-read.csv('Funguild/All_Funguild_10.3.25.csv')

#taxonomy table
otu<-read.csv('Raw_data/Updated_Data/demultiplexed.cleaned.combined.cf.fasta.blast.i97.a95.csv')
#remove first three rows that do not have taxonomy ids and 2 rows with totals and sample counts 
otu<-otu[-c(1:3),-c(2:3) ]

library(readxl)
Fun_Traits<-read_excel("Processed_data/Polme_etal_2021_FungalTraits.xlsx")
tax <- otu %>%
  separate(taxonomy, into = c("Percentage", "Species_Name", "Accession", "SH_ID", "Reference", "Taxonomy_Details"), sep = "\\|", extra = "drop") %>%
  separate(Taxonomy_Details, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";") %>%
  select(SH_ID:Species)%>%
  mutate(across(Kingdom:Species, ~ gsub('^[a-z]__', '', .))) %>%  # Remove prefixes
  mutate(across(Kingdom:Species, ~ case_when(
    str_detect(.x, "unclassified") ~ NA_character_,  # Change 'unclassified' to NA
    str_detect(.x, "Incertae_sedis") ~ NA_character_,  # Change unknown levels with I_s to NA
    TRUE ~ .x   ))) %>% # Keep original values otherwise
  mutate(Species = ifelse(str_detect(Species, '_sp$'), NA, Species))%>%  # Remove '_sp'
  merge(FunGuildData_org,
        by.x = "Genus",   # 'Genus' column from fungal traits dataset
        by.y = "taxon",   # 'taxon' column from FunGuild dataset
        all.x = TRUE)%>%     # Keep all rows from ft, even if no match in FunGuild
    filter(confidenceRanking == "Highly Probable" | SH_ID %in% c("SH1243176.09FU")) %>%  # Keep "Highly Probable" & specified OTUs
  filter(!str_detect(guild, "-") | SH_ID %in% c("SH1243176.09FU")) %>%  # Exclude multiple guilds except for specified OTUs
  mutate(guild2 = case_when(
    str_detect(guild, "Arbuscular Mycorrhizal") ~ "Arbuscular Mycorrhizal",
    str_detect(guild, "Orchid Mycorrhizal") ~ "Orchid Mycorrhizal",
    str_detect(guild, "Ericoid Mycorrhizal") ~ "Ericoid Mycorrhizal",
    str_detect(guild, "Ectomycorrhizal") ~ "Ectomycorrhizal",
    TRUE ~ "other" ),
    note = if_else(SH_ID %in% c("SH1243176.09FU"), 
                 "Included due to manual BLAST confirmation", 
                 NA_character_))%>%  # Add note for manually included OTUs
    left_join(Fun_Traits%>%rename(Genus=GENUS)%>%
  select(Genus,Ectomycorrhiza_exploration_type_template,Ectomycorrhiza_lineage_template))%>%
  rename(exploration_type=Ectomycorrhiza_exploration_type_template,
         Ecm_lineage=Ectomycorrhiza_lineage_template)%>%
  mutate(exploration_type = na_if(exploration_type, "unknown"))  # Convert "Unknown" to NA
  


temp<-tax%>%filter(is.na(guild))


tax_myco<-tax%>%
  filter(str_detect(guild, "mycorrhizal") | str_detect(guild, "Mycorrhizal"))

write.csv(tax_myco, 'Processed_data/Soil_ITS_myco_tax.csv',row.names = FALSE)

rm(FunGuildData_org,Nutrients_Transects, VEG_COVER_Transects)
#change first col name to something easy
names(otu)[1]<-'SH_ID'

# transpose community table, to sample-taxon, for vegan functions and remove taxonomy col
mat<-otu%>%
  select(-last_col())%>%
  remove_rownames()%>%
  column_to_rownames("SH_ID")%>%
  t()

# assess variation in sampling effort, plotting sample effort curves

# assess variation in sampling effort, plotting sample effort curves
filtered_mat <- mat[rowSums(mat) >= 5000, ]  # Keep only samples with at least 5000 reads

filtered_temp <- rarecurve(filtered_mat, step=1000, tidy=TRUE)

temp <- rarecurve(mat, step=1000, tidy=TRUE)

filter<-Blast_ID %>%
  left_join(filtered_temp %>% rename(sample_ID = Site))%>%
  #this removes the samples that were not properly sequenced (Sites 49 Tran 2- Site 63)
  filter(!str_detect(sample_ID, "^CG124\\.ITS\\.CG(9[6-9]|1[0-1][0-9]|12[0-4])$"))

unfilter<-Blast_ID %>%
  left_join(temp %>% rename(sample_ID = Site))

excluded<-anti_join(unfilter,filter)


depth_ID_filter<-filter%>% #samples from site 37.1 and 31.2 exlcuded because of insufficent depth
  mutate(Site_Color = ifelse(as.numeric(as.character(Site)) %in% 49:63, "red", "black"),
         Site_Color = case_when(
           Site == 49 & Transect == 1 ~ "black",  # Change Site 49, Transect 1 to 'black'
           TRUE ~ Site_Color  # Keep all other values unchanged
         ))%>%
  ggplot(aes(x = Sample, y = Species, colour = Site_Color, group = sample_ID)) + 
  scale_color_identity() +  # Uses the colors directly without a legend
  geom_line() #+
#theme(legend.position = "none")
depth_ID_filter

#plotly::ggplotly(depth_ID_filter)
depth_ID_excluded<-excluded%>%
  mutate(Site_Color = ifelse(as.numeric(as.character(Site)) %in% 49:63, "red", "black"),
         Site_Color = case_when(
           Site == 49 & Transect == 1 ~ "black",  # Change Site 49, Transect 1 to 'black'
           TRUE ~ Site_Color  # Keep all other values unchanged
         ))%>%
  ggplot(aes(x = Sample, y = Species, colour = 'black', group = sample_ID,text=Site)) + 
  scale_color_identity() +  # Uses the colors directly without a legend
  geom_line()+
  theme(legend.position = "none")
depth_ID_excluded
plotly::ggplotly(depth_ID_excluded)

library(patchwork)
depth_ID_excluded + depth_ID_filter

rm(depth_ID_excluded,depth_ID_filter,Myco_plant_spp,Fun_Traits,filtered_temp)



# 1 - focus analysis on mycorrhizal fungal taxa (create vector of OTU ids)
myco_otus <- filter(tax_myco, guild%in% c("mycorrhizal"))$SH_ID
  

mat_long<-filtered_mat %>% 
  as.data.frame() %>% 
  dplyr::select(all_of(myco_otus)) %>% # just myco OTUs
  rownames_to_column('sample_ID')%>%
  pivot_longer(cols = starts_with("SH"), 
               names_to = "SH_ID", 
               values_to = "readcount")%>%
  filter(readcount > 0)  # Remove rows with 0 readcount


mat_long_all<-filtered_mat %>% 
  as.data.frame() %>% 
  rownames_to_column('sample_ID')%>%
  pivot_longer(cols = starts_with("SH"), 
               names_to = "SH_ID", 
               values_to = "readcount")%>%
  filter(readcount > 0)  # Remove rows with 0 readcount


# 2 - join the community table with our metadata table (need to also modify the community table so that it can be joined)
dat_myco <- left_join(mat_long, Blast_ID)%>%
  right_join(tax_myco) %>% 
  #this removes samples that are bad (> CG124.ITS.CG96)
  filter(!str_detect(sample_ID, "^CG124\\.ITS\\.CG(9[6-9]|1[0-1][0-9]|12[0-4])$"))%>%
  #this is selecting for the 30 decomp sites which I have CN data for all of
  filter(!is.na(Carbon))
summary(dat_myco)
#select sites that I will use based on sample ID
Site_Sample_IDs<- dat_myco%>%select(sample_ID)%>%distinct()%>%pull()


dat_all<-left_join(mat_long_all, Blast_ID)%>%
  right_join(tax) %>% 
  #this removes the samples I dont like
  filter(!str_detect(sample_ID, "^CG124\\.ITS\\.CG(9[6-9]|1[0-1][0-9]|12[0-4])$"))%>%
  #this is selecting for the 30 decomp sites which I have CN data for all of
  filter(!is.na(Carbon))
summary(dat_all)#no na;s in df

dat_myco%>%select(Site,Transect,Interval,Severity)%>%distinct()%>%group_by(Interval,Severity)%>%summarise(n())


dat_all%>%select(Site,Transect,Interval,Severity)%>%distinct()%>%group_by(Interval,Severity)%>%summarise(n())

#I know I should have used mutate instead of all of these joins....

#myco
dat_myco_RA_soil<-dat_myco%>%
  group_by(sample_ID)%>%
  summarise(reads_samp=sum(readcount))%>%
  left_join(dat_myco)%>%
  left_join(dat_myco%>%
              group_by(Severity)%>%
              summarise(severity_reads=sum(readcount)))%>%
  left_join(dat_myco%>%
              group_by(Interval)%>%
              summarise(interval_reads=sum(readcount)))%>%
  mutate( myco_reads= sum(readcount),
        RA_samp= readcount/reads_samp,
         RA_total_reads= readcount/myco_reads,
         RA_total_interval= readcount/interval_reads,
         RA_total_severity= readcount/severity_reads)


dat_explo_soil<-dat_myco_RA_soil%>%
  filter(!is.na(exploration_type))%>%
  group_by(Site,Transect,exploration_type)%>%
  summarise(explo_count=sum(readcount))%>%
  left_join(dat_myco_RA_soil%>%select(Site,Transect,Severity,Interval,reads_samp)%>%distinct())%>%
  left_join(dat_myco_RA_soil%>%filter(!is.na(exploration_type))%>%
              group_by(Interval)%>%
              summarise(interval_reads_explo=sum(readcount)))%>%
  left_join(dat_myco_RA_soil%>%filter(!is.na(exploration_type))%>%
              group_by(Severity)%>%
              summarise(Severity_reads_explo=sum(readcount)))%>%
  mutate(RA_explo_Interval=explo_count/interval_reads_explo,
         RA_explo_Severity=explo_count/Severity_reads_explo)
#all fungal types
dat_all_RA_soil<-dat_all%>%
  group_by(sample_ID)%>%
  summarise(reads_samp=sum(readcount))%>%
  left_join(dat_all)%>%
  left_join(dat_all%>%
              group_by(Severity)%>%
              summarise(severity_reads=sum(readcount)))%>%
  left_join(dat_all%>%
              group_by(Interval)%>%
              summarise(interval_reads=sum(readcount)))%>%
  mutate( total_reads= sum(readcount),
        RA_samp= readcount/reads_samp,
         RA_total_reads= readcount/total_reads,
         RA_total_interval= readcount/interval_reads,
         RA_total_severity= readcount/severity_reads)

rm(mat,otu,dat_all, mat_long,mat_long_all)

