library(tidyverse)
library(stringr)

#import Blast ID df and clean data
Blast_ID<-read_excel('Raw_data/Updated_Data/ABS.MER.fielddata.Feb.2023_Site.Info_AF.xlsx')
Blast_ID$Site= sub(c('ABS00|ABS0'),'',Blast_ID$Site)
Blast_ID$sample_ID<-as.character(Blast_ID$sample_ID)
Blast_ID$transect<-as.character(Blast_ID$transect)
Blast_ID$sample_ID<-sub(c("-"),".",Blast_ID$sample_ID)
Blast_ID$sample_ID<-sub(c("-"),".",Blast_ID$sample_ID)
names(Blast_ID)[9]<-'Veg_Class_Abv'


Blast_ID<-Blast_ID%>%
  dplyr::mutate(Regime = paste(Interval, Severity, sep = "_"))%>%
  rename(Transect=transect)%>%
  select(-`Fire 3`,-`Interval (yrs)...16`,-`FESM severity category`)

#funguild output
Myco<-read.csv("Funguild/Myco_guilds.csv")#IDing all OTUs that are mycorrhizal 


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
FunGuildData_org <- parse_funguild(url)

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
  filter(confidenceRanking %in% c('Probable',"Highly Probable"))%>%
  mutate(guild= str_remove(guild,'\\|'),guild= str_remove(guild,'\\|'), #not great code, but easy solution
         trophicMode= str_remove(trophicMode,' '),
         trophicMode= str_replace(trophicMode, 'Pathotroph-Pathotroph-Saprotroph', 'Pathotroph-Saprotroph' ),
         trophicMode= str_replace(trophicMode, 'Saprotroph-Saprotroph-Symbiotroph ', 'Saprotroph-Symbiotroph'))%>%
  mutate(guild2 = case_when(trophicMode == 'Saprotroph' ~ 'Saprotroph',
                            str_detect(guild, "Arbuscular Mycorrhizal") ~ 'Arbuscular Mycorrhizal',
                            str_detect(guild, "Ectomycorrhizal") ~ 'Ectomycorrhizal',
                            str_detect(guild, "Orchid Mycorrhizal") ~ 'Orchid Mycorrhizal',
                            str_detect(guild, "Ericoid Mycorrhizal") ~ 'Ericoid Mycorrhizal',
                            str_detect(guild, "Endophyte") ~ 'Endophyte',
                            str_detect(guild, "Epiphyte") ~ 'Epiphyte',
                            str_detect(Family, "Mortierellaceae") ~ 'Saprotroph',
                            str_detect(guild, "Parasite") ~ 'Parasite',
                            str_detect(guild, "Pathogen") ~ 'Pathogen',
                            TRUE ~ 'unassigned'))%>%
  left_join(Fun_Traits%>%rename(Genus=GENUS)%>%
  select(Genus,Ectomycorrhiza_exploration_type_template,Ectomycorrhiza_lineage_template))%>%
  rename(exploration_type=Ectomycorrhiza_exploration_type_template,
         Ecm_lineage=Ectomycorrhiza_lineage_template)%>%
  mutate(across(c(exploration_type,Ecm_lineage), ~replace_na(.x, "unknown")))# Replace NA with "Unknown" 

tax_myco<-tax%>%
  mutate (guild = case_when(Genus %in% Myco$Genus ~ 'mycorrhizal',
                          TRUE ~ guild ))%>%
  filter(guild=='mycorrhizal')

write.csv(tax_myco, 'Processed_data/Soil_ITS_myco_tax.csv',row.names = FALSE)

rm(FunGuildData_org,Myco)
#change first col name to something easy
names(otu)[1]<-'SH_ID'

# transpose community table, to sample-taxon, for vegan functions and remove taxonomy col
mat<-otu%>%
  select(-last_col())%>%
  remove_rownames()%>%
  column_to_rownames("SH_ID")%>%
  t()

# assess variation in sampling effort, plotting sample effort curves
mat <- mat[rowSums(mat) >= 5000, ]  # Keep only samples with at least 1000 reads


# 1 - focus analysis on mycorrhizal fungal taxa (create vector of OTU ids)
myco_otus <- filter(tax_myco, guild%in% c("mycorrhizal"))$SH_ID

mat_long<-mat %>% 
  as.data.frame() %>% 
  dplyr::select(all_of(myco_otus)) %>% # just myco OTUs
  rownames_to_column('sample_ID')%>%
  pivot_longer(cols = starts_with("SH"), 
               names_to = "SH_ID", 
               values_to = "readcount")

mat_long_all<-mat %>% 
  as.data.frame() %>% 
  rownames_to_column('sample_ID')%>%
  pivot_longer(cols = starts_with("SH"), 
               names_to = "SH_ID", 
               values_to = "readcount")

# 2 - join the community table with our metadata table (need to also modify the community table so that it can be joined)
dat_myco <- left_join(mat_long, Blast_ID)%>%
  left_join(tax_myco) %>% filter(!is.na(guild))%>%
  #this removes samples that are bad
  filter(!str_detect(sample_ID, "^CG124\\.ITS\\.CG(9[6-9]|1[0-1][0-9]|12[0-2])$"))

dat_all<-left_join(mat_long_all, Blast_ID)%>%
  left_join(tax) %>% filter(!is.na(guild))%>%
  #this removes the samples I dont like
  filter(!str_detect(sample_ID, "^CG124\\.ITS\\.CG(9[6-9]|1[0-1][0-9]|12[0-2])$"))


#myco
myco_reads<-dat_myco%>%
  summarise(sum(readcount))%>%
  pull()


dat_myco_RA<-dat_myco%>%
  group_by(sample_ID)%>%
  summarise(reads_samp=sum(readcount))%>%
  left_join(dat_myco)%>%
  left_join(dat_myco%>%
              group_by(Severity)%>%
              summarise(severity_reads=sum(readcount)))%>%
  left_join(dat_myco%>%
              group_by(Interval)%>%
              summarise(interval_reads=sum(readcount)))%>%
  mutate(RA_samp= readcount/reads_samp,
         RA_total_reads= readcount/myco_reads,
         RA_total_interval= readcount/interval_reads,
         RA_total_severity= readcount/severity_reads)

#all fungal types
total_reads<-dat_all%>%
  summarise(sum(readcount))%>%
  pull()

dat_all_RA<-dat_all%>%
  group_by(sample_ID)%>%
  summarise(reads_samp=sum(readcount))%>%
  left_join(dat_all)%>%
  left_join(dat_all%>%
              group_by(Severity)%>%
              summarise(severity_reads=sum(readcount)))%>%
  left_join(dat_all%>%
              group_by(Interval)%>%
              summarise(interval_reads=sum(readcount)))%>%
  mutate(RA_samp= readcount/reads_samp,
         RA_total_reads= readcount/total_reads,
         RA_total_interval= readcount/interval_reads,
         RA_total_severity= readcount/severity_reads)

rm(mat,otu,dat_all,dat_myco, mat_long,mat_long_all)

