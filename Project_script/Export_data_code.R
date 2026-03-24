########################
# 8. Export Data for Publication
########################

# Function to save files
export_clean_data <- function(df, round_name) {
  
  # 1. Create Taxa Data (ID + Taxa cols)
  # Identify the ID column based on the round
  id_col <- if(round_name == "Rnd1") "Tube_ID" else "ID"
  
  taxa_data <- df %>%
    select(all_of(id_col), starts_with("ITSall"))
  
  write.csv(taxa_data, paste0("Project_script/taxa_data_", round_name, ".csv"), row.names = FALSE)
  
  # 2. Create Site Data (All cols EXCEPT Taxa)
  site_data <- df %>%
    select(-starts_with("ITSall"))
  
  write.csv(site_data, paste0("Project_script/site_data_", round_name, ".csv"), row.names = FALSE)
  
  message(paste(round_name, "files created successfully."))
}

# Execute export
export_clean_data(Data_rnd1, "Rnd1")
export_clean_data(Data_rnd2, "Rnd2")
