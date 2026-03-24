# Data Description for: Decoupled responses of mycorrhizal fungal communities and function to recurrent wildfire

## Overview
This repository contains data files for analyzing mycorrhizal fungal responses to historical fire severity and frequency in Australian forests.

## File List
1. site_data_Rnd1.csv: Site metadata and environmental variables for the 1st mycelial collection.
2. site_data_Rnd2.csv: Site metadata and environmental variables for the 2nd mycelial collection.
3. taxa_data_Rnd1.csv: OTU/Taxa abundance table for the 1st mycelial collection.
4. taxa_data_Rnd2.csv: OTU/Taxa abundance table for the 2nd mycelial collection.
5. taxonomy_Rnd1.csv: Taxonomic ID's for each OTU for the 1st mycelial collection.
6. taxonomy_Rnd2.csv: Taxonomic ID's for each OTU for the 2nd mycelial collection.


## Variable Descriptions (site_data_Rnd1 & site_data_Rnd2)

### Identifiers
* Site: Unique identifier for the study site.
* Transect: Transect number within the site.
* Location: Sampling location along transect within the site.
* Tube_ID: Unique sample identifier linking site data to taxa data.

### Experimental Treatments
* Fire.Interval: Frequency of fire (e.g., "Long", "Short").
* Fire.Severity: Severity of the fire event (e.g., "High", "Low").
* Fire: Combined factor of Severity x Interval.

### Environmental Variables
* Latitude: GPS Latitude of the sampling point.
* Longitude: GPS Longitude of the sampling point.
* Ortho_P_mg_kg: Soil Orthophosphate concentration (mg/kg).
* Nitrate_mg_kg: Soil Nitrate concentration (mg/kg).
* Ammonia_mg_kg: Soil Ammonia concentration (mg/kg).
* pH: Soil pH.

### Vegetation & Host Metrics
* perc_myco_host_freq: Frequency of mycorrhizal host plants at the site.
* AM_host: Frequency of Arbuscular Mycorrhizal host plants.
* EcM_host: Frequency of Ectomycorrhizal host plants.

### Fungal Metrics
* Biomass_day: Raw fungal biomass divided by days installed (mg/day).
* ratio_myco_nonmyco: Ratio of mycorrhizal to non-mycorrhizal sequencing reads (QC metric).
* myco_reads: Total count of mycorrhizal fungal reads.