###############################################################################
# Decoupled responses of mycorrhizal fungal communities and function
# to recurrent wildfire
#
# Author: Solomon Maerowitz-McMahan
# Date:   2025-11-12
#
# This script reproduces the core analyses for the manuscript:
#   1. Fire regime effects on soil nutrients
#   2. Mycorrhizal community structure (alpha/beta diversity, CAP ordination)
#   3. Mycorrhizal biomass production models
#   4. Hyphal stoichiometry (C:N, C:P, N:P) and elemental content
#   5. HMSC joint species distribution modelling
#
# Data files are in Project_script/data/
# See DATA_README.md for variable descriptions.
###############################################################################

library(tidyverse)
library(lme4)
library(car)
library(vegan)
library(emmeans)
library(fossil)
library(broom.mixed)

###############################################################################
# 1. DATA LOADING
###############################################################################

site_rnd1   <- read.csv("Project_script/data/site_data_Rnd1.csv")
taxa_rnd1   <- read.csv("Project_script/data/taxa_data_Rnd1.csv")
tax_id_rnd1 <- read.csv("Project_script/data/taxonomy_Rnd1.csv")

site_rnd2   <- read.csv("Project_script/data/site_data_Rnd2.csv")
taxa_rnd2   <- read.csv("Project_script/data/taxa_data_Rnd2.csv")
tax_id_rnd2 <- read.csv("Project_script/data/taxonomy_Rnd2.csv")

Data_rnd1 <- site_rnd1 %>%
  bind_cols(taxa_rnd1) %>%
  mutate(
    Biomass_day = if_else(ratio_myco_nonmyco < 0.5, NA_real_, Biomass_day),
    across(c(Biomass_day, Ortho_P_mg_kg, Nitrate_mg_kg, Ammonia_mg_kg, pH),
           ~log10(.), .names = "{.col}_log"),
    biomass_kg_ha_day = Biomass_day * (1e+03 / 30),
    Fire = paste(Fire.Severity, "x", Fire.Interval, sep = "\n")
  )

Data_rnd2 <- site_rnd2 %>%
  bind_cols(taxa_rnd2) %>%
  mutate(
    Biomass_day = if_else(ratio_myco_nonmyco < 0.5, NA_real_, Biomass_day),
    across(c(Biomass_day, Ortho_P_mg_kg, pH), ~log10(.), .names = "{.col}_log"),
    biomass_kg_ha_day = Biomass_day * (1e+03 / 30)
  )

###############################################################################
# 2. FIRE REGIME EFFECTS ON SOIL NUTRIENTS
###############################################################################

m1 <- lmer(Ortho_P_mg_kg_log ~ Fire.Severity + Fire.Interval + (1|Site/Transect), data = Data_rnd1)
m2 <- lmer(Nitrate_mg_kg_log ~ Fire.Severity + Fire.Interval + (1|Site/Transect), data = Data_rnd1)
m3 <- lmer(Ammonia_mg_kg_log ~ Fire.Severity + Fire.Interval + (1|Site/Transect), data = Data_rnd1)

Nutrient_1st_Rnd_Anova <- bind_rows(
  Anova(m1, test = 'F') %>% as.data.frame() %>% rownames_to_column("Factor") %>% mutate(Nutrient = 'Orthophosphate', Round = 'First'),
  Anova(m2, test = 'F') %>% as.data.frame() %>% rownames_to_column("Factor") %>% mutate(Nutrient = 'Nitrate', Round = 'First'),
  Anova(m3, test = 'F') %>% as.data.frame() %>% rownames_to_column("Factor") %>% mutate(Nutrient = 'Ammonia', Round = 'First')
)

m_Phos <- lmer(log10(Ortho_P_mg_kg) ~ Fire.Interval + (1|Site/Transect), data = Data_rnd2)
Anova_resin_Phos <- Anova(m_Phos, test = 'F') %>%
  round(2) %>% as.data.frame() %>%
  rownames_to_column("Factor") %>%
  mutate(Nutrient = "Orthophosphate", Round = "Second")

Nutrient_Anova <- bind_rows(Nutrient_1st_Rnd_Anova, Anova_resin_Phos) %>%
  arrange(Factor)

###############################################################################
# 3. MYCORRHIZAL COMMUNITY ANALYSES
###############################################################################

# --- 3a. Alpha diversity (1st collection) ---

alpha_diversity <- Data_rnd1 %>%
  select(Tube_ID, starts_with("ITSall")) %>%
  rowwise() %>%
  mutate(
    Shannon  = diversity(c_across(starts_with("ITSall")), index = "shannon"),
    Simpson  = diversity(c_across(starts_with("ITSall")), index = "simpson"),
    Chao1    = chao1(c_across(starts_with("ITSall"))),
    Observed = sum(c_across(starts_with("ITSall")) > 0),
    Pielou   = if_else(Observed > 1, Shannon / log(Observed), NA_real_)
  ) %>%
  ungroup() %>%
  select(Tube_ID, Shannon, Simpson, Chao1, Observed, Pielou)

alpha_meta <- alpha_diversity %>%
  left_join(Data_rnd1) %>%
  filter(ratio_myco_nonmyco > 0.5)

analyze_alpha <- function(metric) {
  model <- lmer(as.formula(paste0(metric, " ~ Fire.Severity * Fire.Interval + (1|Site/Transect)")), data = alpha_meta)
  Anova(model, test = "F") %>% round(2) %>% as.data.frame() %>% rownames_to_column("Factor") %>% mutate(Metric = metric)
}

anova_results_1st <- bind_rows(
  analyze_alpha("Shannon"),
  analyze_alpha("Pielou"),
  analyze_alpha("Chao1")
) %>% mutate(source = "1st mycelial collection")

# --- 3b. Beta diversity (1st collection) ---

mat_myco <- Data_rnd1 %>% select(starts_with("ITSall"))

permanova_res <- adonis2(mat_myco ~ Fire.Severity * Fire.Interval, data = Data_rnd1,
                         distance = 'robust.aitchison', by = 'margin') %>%
  as.data.frame() %>% rownames_to_column("Factor") %>%
  mutate(Sample_Type = "1st mycelial collection")

# --- 3c. CAP ordination (1st collection) ---

cap_mod <- capscale(mat_myco ~ Fire, data = Data_rnd1, distance = 'robust.aitchison', add = TRUE)
prop_var <- round(cap_mod$CCA$eig / cap_mod$tot.chi * 100, 1)
scrs <- scores(cap_mod, tidy = TRUE)
scrs_site <- filter(scrs, score == "sites")

interval_colors <- c("Long" = "darkred", "Short" = "orange")
severity_linetypes <- c("High" = "solid", "Low" = "dashed")

p3 <- cbind(Data_rnd1, scrs_site) %>%
  ggplot(aes(x = CAP1, y = CAP2)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey70", linetype = 2) +
  geom_point(aes(colour = Fire.Interval, shape = Fire.Severity), size = 4, stroke = 1.5) +
  stat_ellipse(aes(colour = Fire.Interval, linetype = Fire.Severity), linewidth = 1) +
  scale_colour_manual(values = interval_colors) +
  scale_linetype_manual(values = severity_linetypes, guide = "none") +
  scale_shape_manual(values = c("High" = 16, "Low" = 1)) +
  labs(x = paste0("CAP1 (", prop_var[1], "%)"),
       y = paste0("CAP2 (", prop_var[2], "%)"),
       colour = "Fire Frequency", shape = "Fire Severity") +
  theme_classic() +
  theme(legend.position = 'top',
        axis.text = element_text(size = 12, color = 'black'),
        axis.title = element_text(size = 14, color = 'black'),
        legend.text  = element_text(size = 10),
        legend.title = element_text(size = 10)) +
  guides(shape = guide_legend(override.aes = list(size = 8)),
         colour = guide_legend(override.aes = list(size = 8)))

# --- 3d. Alpha diversity (2nd collection) ---

alpha_div_2 <- Data_rnd2 %>%
  rowwise() %>%
  mutate(
    Shannon  = diversity(c_across(starts_with("ITSall")), index = "shannon"),
    Simpson  = diversity(c_across(starts_with("ITSall")), index = "simpson"),
    Chao1    = chao1(c_across(starts_with("ITSall"))),
    Observed = sum(c_across(starts_with("ITSall")) > 0),
    Pielou   = if_else(Observed > 1, Shannon / log(Observed), NA_real_)
  ) %>%
  ungroup()

run_alpha_model <- function(metric, data) {
  formula <- as.formula(paste0(metric, " ~ Fire.Interval + (1|Site/Transect)"))
  model <- lmer(formula, data = data)
  Anova(model, test = "F") %>%
    round(2) %>% as.data.frame() %>%
    rownames_to_column("Factor") %>%
    mutate(Metric = metric)
}

anova_results_2nd <- bind_rows(
  run_alpha_model("Shannon", alpha_div_2),
  run_alpha_model("Pielou", alpha_div_2),
  run_alpha_model("Chao1", alpha_div_2)
) %>% mutate(source = "2nd mycelial collection")

# --- 3e. Beta diversity (2nd collection) ---

Data_rnd2 <- Data_rnd2 %>% filter(!if_all(starts_with("ITSall"), ~ . == 0))
mat_myco_2 <- Data_rnd2 %>% select(starts_with("ITSall"))

adonis_rnd2 <- adonis2(mat_myco_2 ~ Fire.Interval, data = Data_rnd2,
                       distance = "robust.aitchison", by = "margin") %>%
  as.data.frame() %>% rownames_to_column("Factor") %>%
  mutate(Sample_Type = "2nd mycelial collection")

cap_model_2 <- capscale(mat_myco_2 ~ Fire.Interval, data = Data_rnd2,
                        distance = "robust.aitchison", add = TRUE)
proportions_2 <- round(cap_model_2$CCA$eig / cap_model_2$tot.chi * 100, 1)
scrs_2 <- scores(cap_model_2, tidy = TRUE) %>% filter(score == "sites")

p2 <- cbind(Data_rnd2, scrs_2) %>%
  ggplot(aes(x = CAP1, y = MDS1)) +
  geom_vline(xintercept = 0, color = "grey70", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey70", linetype = 2) +
  geom_point(aes(colour = Fire.Interval), size = 8) +
  stat_ellipse(aes(color = Fire.Interval), level = 0.95, linewidth = 2) +
  scale_colour_manual(values = interval_colors) +
  labs(x = paste0("CAP1 (", proportions_2[1], "%)"), y = "MDS1",
       colour = "Fire frequency") +
  xlim(min(scrs_2$CAP1) - 1, max(scrs_2$CAP1) + 0.5) +
  theme_bw(base_size = 20) +
  theme(axis.text = element_text(color = "black", size = 36),
        axis.title = element_text(color = "black", size = 36),
        legend.position = "top",
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30))

###############################################################################
# 4. MYCORRHIZAL BIOMASS PRODUCTION
###############################################################################

# --- 4a. 1st collection ---

m_additive <- lmer(Biomass_day_log ~ Fire.Severity + Fire.Interval +
                     Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg +
                     perc_myco_host_freq + (1 | Site/Transect),
                   data = Data_rnd1)

m_interaction <- lmer(Biomass_day_log ~ Fire.Severity * Fire.Interval +
                        Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg +
                        perc_myco_host_freq + (1 | Site/Transect),
                      data = Data_rnd1)

AIC(m_additive, m_interaction)
m_biomass_day_resins <- m_additive

fixed_effects <- tidy(m_biomass_day_resins, effects = "fixed") %>%
  filter(term != "(Intercept)") %>%
  mutate(term = str_remove(term, 'Low|Short')) %>%
  select(term, estimate, std.error)

anova_df <- Anova(m_biomass_day_resins, test = "F") %>%
  as.data.frame() %>% rownames_to_column("Predictor")

summary_table_1st <- fixed_effects %>%
  rename(Slope = estimate, Std_Error = std.error) %>%
  left_join(anova_df[, c("Predictor", "Df", "Df.res", "F", "Pr(>F)")], by = c("term" = "Predictor")) %>%
  rename(p = `Pr(>F)`) %>%
  mutate(term = term %>% str_replace_all("_mg_kg", "") %>% str_replace_all("_", " ") %>% str_trim()) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  mutate(Round = '1st') %>%
  relocate(Round, term, Slope, Std_Error, F, Df, Df.res, p)

# --- 4b. 2nd collection ---

m_biomass_2 <- lmer(Biomass_day_log ~ Fire.Interval + Ortho_P_mg_kg + perc_myco_host_freq +
                      (1 | Site / Transect), data = Data_rnd2)

Anova_biomass_2 <- Anova(m_biomass_2, test = "F") %>% round(2)

fixed_effects_2 <- tidy(m_biomass_2, effects = "fixed") %>%
  filter(term != "(Intercept)") %>%
  mutate(Predictor = str_remove_all(term, "Low|Short"),
         Predictor = str_replace_all(Predictor, "_mg_kg", ""),
         Predictor = str_replace_all(Predictor, "_", " ")) %>%
  select(Predictor, estimate, std.error)

anova_df_2 <- Anova_biomass_2 %>%
  as.data.frame() %>% rownames_to_column("Predictor") %>%
  mutate(Predictor = str_replace_all(Predictor, "_mg_kg", ""),
         Predictor = str_replace_all(Predictor, "_", " "),
         Predictor = trimws(Predictor))

summary_table_2nd <- fixed_effects_2 %>%
  rename(Slope = estimate, Std_Error = std.error) %>%
  left_join(anova_df_2, by = "Predictor") %>%
  rename(p = `Pr(>F)`) %>%
  mutate(across(where(is.numeric), round, 2), Round = "2nd") %>%
  relocate(Round, Predictor, F, Slope, Std_Error, Df, Df.res, p)

###############################################################################
# 5. HYPHAL STOICHIOMETRY
###############################################################################

CNP_clean <- read.csv("Processed_data/CNP_clean.csv")

Bag_data <- Data_rnd1 %>%
  group_by(Site, Transect) %>%
  summarise(across(c(Ortho_P_mg_kg, Nitrate_mg_kg, Ammonia_mg_kg), ~mean(.x, na.rm = TRUE))) %>%
  left_join(CNP_clean, by = c("Site", "Transect")) %>%
  left_join(Data_rnd1 %>% distinct(Site, Transect, AM_host, EcM_host, perc_myco_host_freq))

stoich_formula <- function(response) {
  as.formula(paste0("log10(", response, ") ~ Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg + Fire.Interval + Fire.Severity + perc_myco_host_freq + (1 | Site)"))
}

C_N_model       <- lmer(stoich_formula("C_N"), data = Bag_data)
C_P_model       <- lmer(stoich_formula("C_P"), data = Bag_data)
N_P_model       <- lmer(stoich_formula("N_P"), data = Bag_data)
Carb_Hyph_model <- lmer(Carbon ~ Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg + Fire.Interval + Fire.Severity + perc_myco_host_freq + (1 | Site), data = Bag_data)
Nitrog_Hyph_model <- lmer(Nitrogen ~ Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg + Fire.Interval + Fire.Severity + perc_myco_host_freq + (1 | Site), data = Bag_data)
Phos_Hyph_model <- lmer(Percent_Phos_ ~ Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg + Fire.Interval + Fire.Severity + perc_myco_host_freq + (1 | Site), data = Bag_data)

stoich_anovas <- list(
  C_N = Anova(C_N_model, test = "F"),
  C_P = Anova(C_P_model, test = "F"),
  N_P = Anova(N_P_model, test = "F"),
  Carbon = Anova(Carb_Hyph_model, test = "F"),
  Nitrogen = Anova(Nitrog_Hyph_model, test = "F"),
  Phosphorus = Anova(Phos_Hyph_model, test = "F")
)

###############################################################################
# 6. HMSC — JOINT SPECIES DISTRIBUTION MODELLING
###############################################################################
#
# This section defines, fits, evaluates, and extracts parameter estimates from
# Hierarchical Models of Species Communities (HMSC; Ovaskainen & Abrego 2020).
# Model fitting is computationally intensive (hours–days); fitted model objects
# are saved to HMSC_ABS/models/ and are not tracked by git.
#
# Required additional packages:
library(Hmsc)
library(ape)
library(coda)
library(colorspace)
library(corrplot)

# --- 6a. Data preparation ---

localDir  <- "HMSC_ABS"
dataDir   <- file.path(localDir, "data")
modelDir  <- file.path(localDir, "models")
resultDir <- file.path(localDir, "results")
if (!dir.exists(resultDir)) dir.create(resultDir)

hmsc_data <- read.csv(file.path(dataDir, "Bag_data.csv"))
PhyData   <- read.csv(file.path(dataDir, "Trait_Phylo_data.csv")) %>%
  mutate(across(everything(), as.factor)) %>%
  select(-exploration_type, -Ecm_lineage, -mean_gs)

TrData <- read.csv(file.path(dataDir, "Trait_Phylo_data.csv")) %>%
  mutate(across(everything(), as.factor)) %>%
  select(OTU:last_col())

otu_order <- TrData$OTU

XData <- hmsc_data %>%
  select(Fire.Interval, Fire.Severity, Ortho_P_mg_kg, Nitrate_mg_kg,
         Ammonia_mg_kg, pH, readcount, log10_biomass_day, perc_myco_host_freq) %>%
  mutate(Fire.Interval = as.factor(Fire.Interval),
         Fire.Severity = as.factor(Fire.Severity),
         Ortho_P_mg_kg = log10(Ortho_P_mg_kg),
         Nitrate_mg_kg = log10(Nitrate_mg_kg),
         Ammonia_mg_kg = log10(Ammonia_mg_kg),
         pH = log10(pH))

YData <- hmsc_data %>%
  select(starts_with('ITSall')) %>%
  select(all_of(otu_order))

Y <- 1 * (as.matrix(YData) > 0)

# --- 6b. Study design & random effects ---

studyDesign <- data.frame(
  Location = hmsc_data$Location,
  Transect = hmsc_data$Transect,
  Site     = hmsc_data$Site
) %>% mutate(across(everything(), as.factor))

rL_location <- HmscRandomLevel(units = studyDesign$Location)
rL_transect <- HmscRandomLevel(units = studyDesign$Transect)
rL_site     <- HmscRandomLevel(units = studyDesign$Site)

# --- 6c. Model specification ---

XFormula <- ~Fire.Severity + Fire.Interval + log(readcount) +
  log10_biomass_day + perc_myco_host_freq +
  Ortho_P_mg_kg + Nitrate_mg_kg + Ammonia_mg_kg + pH

taxonomicTree <- as.phylo(~phylum/class/order/family/genus/OTU, data = PhyData, collapse = FALSE)
taxonomicTree$edge.length <- rep(1, length(taxonomicTree$edge))

TrFormula <- ~exploration_type + Ecm_lineage + mean_gs

Y_cols <- colnames(Y)
TrData <- TrData %>%
  mutate(OTU = factor(OTU, levels = Y_cols)) %>%
  arrange(OTU) %>%
  column_to_rownames('OTU')

m <- Hmsc(Y = Y, XData = XData, XFormula = XFormula,
          phyloTree = taxonomicTree,
          TrData = TrData, TrFormula = TrFormula,
          studyDesign = studyDesign,
          ranLevels = list(Location = rL_location,
                           Transect = rL_transect,
                           Site = rL_site),
          distr = "probit")

models <- list(m)
names(models) <- c("presence-absence model")
save(models, file = file.path(modelDir, "unfitted_models.RData"))

# --- 6d. MCMC fitting ---

load(file = file.path(modelDir, "unfitted_models.RData"))
nm <- length(models)
samples_list <- c(250, 250, 250)
thin_list    <- c(250, 500, 1000)
nChains  <- 2
nParallel <- nChains

for (Lst in seq_along(samples_list)) {
  thin    <- thin_list[Lst]
  samples <- samples_list[Lst]
  filename <- file.path(modelDir, paste0("models_thin_", thin, "_samples_", samples, "_chains_", nChains, ".Rdata"))
  if (file.exists(filename)) {
    message("Model already fitted: ", filename)
    next
  }
  message("Fitting thin=", thin, " samples=", samples, " at ", date())
  for (mi in 1:nm) {
    models[[mi]] <- sampleMcmc(models[[mi]], samples = samples, thin = thin,
                                adaptNf = rep(ceiling(0.4 * samples * thin), models[[mi]]$nr),
                                transient = ceiling(0.5 * samples * thin),
                                nChains = nChains, nParallel = nParallel)
  }
  save(models, file = filename)
}

# --- 6e. Convergence evaluation ---

mpost <- convertToCodaObject(models[[1]], spNamesNumbers = c(TRUE, FALSE),
                             covNamesNumbers = c(TRUE, FALSE))

psrf.beta  <- gelman.diag(mpost$Beta, multivariate = FALSE)$psrf
psrf.omega <- gelman.diag(mpost$Omega[[1]], multivariate = FALSE)$psrf

par(mfrow = c(1, 2))
hist(psrf.beta,  main = "Probit", xlab = "PSRF (Beta)")
hist(psrf.omega, main = "Probit", xlab = "PSRF (Omega)")

# --- 6f. Model fit (cross-validation) ---

nfolds <- 4
for (mi in 1:nm) {
  preds <- computePredictedValues(models[[mi]])
  MF    <- evaluateModelFit(hM = models[[mi]], predY = preds)
  partition <- createPartition(models[[mi]], nfolds = nfolds)
  preds_cv  <- computePredictedValues(models[[mi]], partition = partition, nParallel = nParallel)
  MFCV      <- evaluateModelFit(hM = models[[mi]], predY = preds_cv)
  WAIC      <- computeWAIC(models[[mi]])
}

# --- 6g. Variance partitioning ---

preds <- computePredictedValues(models[[1]])
VP <- computeVariancePartitioning(models[[1]])
plotVariancePartitioning(hM = models[[1]], VP = VP, main = "Variance Partitioning")

# --- 6h. Beta parameter estimates ---

postBeta <- getPostEstimate(models[[1]], parName = "Beta")
plotBeta(models[[1]], post = postBeta, supportLevel = 0.95, param = "Sign",
         plotTree = !is.null(models[[1]]$phyloTree),
         covNamesNumbers = c(TRUE, FALSE),
         spNamesNumbers = c(models[[1]]$ns <= 30, FALSE),
         cex = c(0.6, 0.6, 0.8))

# --- 6i. Species associations (Omega) ---

OmegaCor <- computeAssociations(models[[1]])
for (r in 1:models[[1]]$nr) {
  support_level <- 0.9
  toPlot <- ((OmegaCor[[r]]$support > support_level) +
               (OmegaCor[[r]]$support < (1 - support_level)) > 0) * sign(OmegaCor[[r]]$mean)
  plotOrder <- corrMatOrder(OmegaCor[[r]]$mean, order = "AOE")
  corrplot(toPlot[plotOrder, plotOrder], method = "color",
           col = colorRampPalette(c("blue", "white", "red"))(3),
           mar = c(0, 0, 1, 0),
           main = paste0("Associations: ", names(models[[1]]$ranLevels)[[r]]),
           cex.main = 0.8)
}
