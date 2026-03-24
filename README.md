# ABS_FIRE_MYCO

**Decoupled responses of mycorrhizal fungal communities and function to recurrent wildfire**

Analysis of mycorrhizal fungal biomass, community composition, hyphal stoichiometry, and species associations across fire severity and frequency gradients in Australian sclerophyll forests.

## Repository Structure

```
ABS_FIRE_MYCO/
├── Project_script/           # Core reproducible analyses
│   ├── ABS_analysis_summary.R    # Complete analysis pipeline (sections 1–6)
│   ├── ABS_code_SMM.R            # Original working script
│   ├── data/                     # Clean data files for reproduction
│   └── DATA_README.md            # Variable descriptions
│
├── HMSC_ABS/                 # HMSC joint species distribution modelling
│   ├── HMSC_S1–S7 scripts        # Model definition, fitting, evaluation
│   ├── data/                     # HMSC input data
│   ├── functions_hmsc.R          # Helper functions
│   └── models/ & results/        # Fitted models & outputs (git-ignored)
│
├── scripts/                  # Supporting analysis scripts
│   ├── data_prep/                # Data import and joining
│   ├── biomass/                  # Biomass & stoichiometry analyses
│   ├── community/                # Community composition & co-occurrence
│   ├── nutrients/                # Soil nutrient scripts
│   ├── modeling/                 # Statistical modelling
│   ├── traits/                   # Fungal trait databases
│   ├── mapping/                  # Site mapping & vegetation
│   └── misc/                     # Miscellaneous scripts
│
├── Bag_data_seq/             # ITS sequencing data processing
├── Soil_ITS_Scripts/         # Soil ITS community analyses
├── Nutrient_scripts/         # Nutrient analyses & graphing
├── Site_Meta_Scripts/         # Climate & site metadata extraction
├── Stoich_CNP/               # Hyphal C:N:P stoichiometry
├── HMSC_output/              # Post-HMSC association analyses
├── Tables/                   # Statistical output tables
└── ABS_Second_Rnd/           # 2nd mycelial collection pipeline
```

## Reproducing the Analyses

1. Open `ABS_FIRE_MYCO.Rproj` in RStudio
2. Run `Project_script/ABS_analysis_summary.R`
   - Sections 1–5 use data in `Project_script/data/` and `Processed_data/`
   - Section 6 (HMSC) requires fitted models in `HMSC_ABS/models/` — these are large and git-ignored. Re-fitting requires running `HMSC_ABS/HMSC_S1_define_and_fit_models.R`

## Required R Packages

```r
install.packages(c(
  "tidyverse", "lme4", "car", "vegan", "emmeans", "fossil",
  "broom.mixed", "Hmsc", "ape", "coda", "colorspace", "corrplot"
))
```

## Data

See [`Project_script/DATA_README.md`](Project_script/DATA_README.md) for variable descriptions.

Clean data files in `Project_script/data/`:
- `site_data_Rnd1.csv` / `site_data_Rnd2.csv` — Site metadata
- `taxa_data_Rnd1.csv` / `taxa_data_Rnd2.csv` — OTU abundance tables
- `taxonomy_Rnd1.csv` / `taxonomy_Rnd2.csv` — Taxonomic assignments

## Author

Solomon Maerowitz-McMahan  
