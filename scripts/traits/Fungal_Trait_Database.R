# Set working directory (commented out)
# setwd("\\\\ad.uws.edu.au/dfshare/HomesCMB$/90957135/My Documents/R")

# Install necessary packages (commented out)
# install.packages("devtools")
# devtools::install_github("ropenscilabs/datastorr")  # Install 'datastorr' package from GitHub
# devtools::install_github("traitecoevo/fungaltraits")  # Install 'fungaltraits' package from GitHub

# Load required libraries
library(fungaltraits)  # Load fungal trait data
library(ggplot2)       # For visualization
library(httr)          # For working with HTTP requests
library(jsonlite)      # For parsing JSON data
library(lubridate)     # For handling dates
library(XML)           # For parsing HTML/XML data

# Load fungal trait dataset
ft <- fungal_traits()  

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

# Merge fungal trait data with FunGuild data by matching genus names
guil.trait <- merge(ft, FunGuildData_org,
                    by.x = "Genus",   # 'Genus' column from fungal traits dataset
                    by.y = "taxon",   # 'taxon' column from FunGuild dataset
                    all.x = TRUE)     # Keep all rows from ft, even if no match in FunGuild

# Check column names of fungal traits dataset
colnames(ft)

# Find overlapping column names between the merged dataset and FunGuild data
intersect(names(guil.trait), names(FunGuildData_org))

# Load string manipulation library
library(stringr)

# Define a list of guild-related keywords
guild_keywords <- c("Ectomycorrhizal", "Arbuscular Mycorrhizal", "Orchid Mycorrhizal", "Ericoid Mycorrhizal",
                    "Wood Saprotroph", "Soil Saprotroph", "Litter Saprotroph", "Undefiend saprotroph",
                    "Plant Pathogen", "Endophyte", "Fungal Parasite", "Animal Parasite")

# Filter dataset to include only rows where 'guild' contains any of the specified keywords
Guilds <- subset(guil.trait, str_detect(guild, str_c(guild_keywords, collapse = "|")))%>%
  filter(confidenceRanking %in% c("Probable", "Highly Probable"))# Keep only specified values


write.csv(Guilds, "Funguild/All_guilds.csv", row.names = FALSE)


# Define a more specific list of mycorrhizal guilds
myco_keywords <- c("Ectomycorrhizal", "Arbuscular Mycorrhizal", "Orchid Mycorrhizal", "Ericoid Mycorrhizal")


myco_Guilds<-guil.trait%>%  filter(confidenceRanking %in% c('Probable',"Highly Probable"))%>%
  select(Genus,guild,higher_clade,studyName,taxonomicLevel:citationSource)%>%
  filter(str_detect(guild, str_c(myco_keywords, collapse = "|")))       
          

myco_Guilds2<-guil.trait%>%  filter(confidenceRanking %in% c('Probable',"Highly Probable"))%>%
  select(Genus,guild,higher_clade,studyName,taxonomicLevel:citationSource)%>%
filter(str_detect(guild, "mycorrhizal") | str_detect(guild, "Mycorrhizal")) %>%
  filter(!Genus %in% c('Leohumicola', 'Rhizoctonia')) # Exclude specific genera    

# Extract unique guild names from the dataset
xaxis <- Guilds$guild  
unique(guil.trait$guild)

# Save filtered mycorrhizal guilds dataset as a CSV file
write.csv(myco_Guilds, "Funguild/Myco_guilds.csv", row.names = FALSE)


# # function for number of observations 
# give.n <- function(x){
#   return(c(y = median(x)*1.05, label = length(x))) 
#   # experiment with the multiplier to find the perfect position
# }
# 

# stotich<-c("tissue_cn", "tissue_cp","tissue_n", "tissue_np","tissue_c") 
# 
# colnames(Guilds)
# 
# Guilds %>% filter(guild == 'Arbuscular Mycorrhizal') %>% select('tissue_cp','citationSource')
# 
# # ggplot(Guilds, aes(x=xaxis, y=tissue_cn, fill=tissue_cn))+
#   stat_summary(geom='bar', fun=mean,na.rm = T, position=position_dodge()) +
#   stat_summary(geom= 'errorbar', fun.data = mean_cl_normal, na.rm = T,
#                position=position_dodge(width=0.9),width=0.1) +
#   stat_summary(fun.data = give.n, geom = "text", fun.y = median) 
# 
# 
# 
# ggplot(Guilds, aes(x=xaxis, y=tissue_p, fill=tissue_cp))+
#   stat_summary(geom='bar', fun=mean,na.rm = T, position=position_dodge()) +
#   stat_summary(geom= 'errorbar', fun.data = mean_cl_normal, na.rm = T,
#                position=position_dodge(width=0.9),width=0.1)+
# stat_summary(fun.data = give.n, geom = "text", fun.y = median) 
# 
# 
# ggplot(Guilds, aes(x=xaxis, y=tissue_n, fill=tissue_n))+
#   stat_summary(geom='bar', fun=mean,na.rm = T, position=position_dodge()) +
#   stat_summary(geom= 'errorbar', fun.data = mean_cl_normal, na.rm = T,
#                position=position_dodge(width=0.9),width=0.1)
# 
# ggplot(Guilds, aes(x=xaxis, y=tissue_np, fill=tissue_np))+
#   stat_summary(geom='bar', fun=mean,na.rm = T, position=position_dodge()) +
#   stat_summary(geom= 'errorbar', fun.data = mean_cl_normal, na.rm = T,
#                position=position_dodge(width=0.9),width=0.1)
# 
# ggplot(Guilds, aes(x=xaxis, y=tissue_c, fill=tissue_cp))+
#   stat_summary(geom='bar', fun=mean,na.rm = T, position=position_dodge()) +
#   stat_summary(geom= 'errorbar', fun.data = mean_cl_normal, na.rm = T,
#                position=position_dodge(width=0.9),width=0.1)
# 
# ggplot(Guilds, aes(x=xaxis, y=tissue_p, fill=tissue_cp))+
#   stat_summary(geom='bar', fun=mean,na.rm = T, position=position_dodge()) +
#   stat_summary(geom= 'errorbar', fun.data = mean_cl_normal, na.rm = T,
#                position=position_dodge(width=0.9),width=0.1)
