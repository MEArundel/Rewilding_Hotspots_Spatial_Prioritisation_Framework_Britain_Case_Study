############################################################
# Rewilding hotspot viability headroom analysis
############################################################
# Description - This script calculates species-specific 
# viability headroom from hotspot × species population-capacity 
# data, expressing both total viable habitat area and largest 
# viable patch area relative to each species’ Minimum Area 
# Requirement (MAR). It applies a log₂(1 + ratio) transformation 
# to reduce the influence of exceptionally large habitat-to-MAR 
# ratios, summarises the transformed values for each hotspot, 
# and exports species-level and hotspot-level headroom tables.


# Purpose:
#
# 1. Read the hotspot × species population-capacity output
#
# 2. Calculate species-specific viability headroom using:
#
#       MAR ratio = viable habitat area / MAR
#
# 3. Calculate both:
#
#       - total viable habitat headroom
#       - largest viable patch headroom
#
# 4. Apply a log2(1 + ratio) transformation to reduce the
#    influence of exceptionally large habitat-to-MAR ratios
#
# 5. Summarise the transformed headroom values at hotspot
#    level
#
# 6. Save species-level and hotspot-level headroom tables
#
############################################################


# ==========================================================
# 1. Packages
# ==========================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "writexl"
)

packages_to_install <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}

library(dplyr)
library(tidyr)
library(writexl)


# ==========================================================
# 2. Paths
# ==========================================================

input_dir <-
  "./results/Rewilding Hotspots/Barrier Polygon Summary/"

input_file <- file.path(
  input_dir,
  "Hotspot_Species_Population_Capacity_FULL.csv"
)

out_dir <- file.path(
  input_dir,
  "Viability Headroom"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ==========================================================
# 3. Check input
# ==========================================================

if (!file.exists(input_file)) {
  
  stop(
    paste0(
      "Could not find hotspot population-capacity table:\n",
      input_file,
      "\n\nRun the hotspot population-capacity analysis first."
    )
  )
}


# ==========================================================
# 4. Read hotspot × species table
# ==========================================================

dat <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_fields <- c(
  "hotspot_id",
  "country",
  "hotspot_area_km2",
  "species_richness",
  "median_faith_PD",
  "species",
  "n_viable_patches",
  "total_viable_patch_area_km2",
  "largest_patch_area_km2",
  "MAR_km2"
)


missing_fields <- setdiff(
  required_fields,
  names(dat)
)


if (length(missing_fields) > 0) {
  
  stop(
    paste0(
      "Input table is missing required field(s):\n",
      paste(
        missing_fields,
        collapse = "\n"
      )
    )
  )
}


# ==========================================================
# 5. Calculate species-level viability headroom
# ==========================================================

# TOTAL MAR RATIO
#
#          total viable habitat area
# V_ij = --------------------------------
#                    MAR_i
#
#
# LARGEST-PATCH MAR RATIO
#
#                 largest viable patch
# Vmax_ij = --------------------------------
#                         MAR_i
#
#
# Because this input contains only species already supported
# by a hotspot, total MAR ratios should normally be >= 1.
#
# Viability headroom is calculated as:
#
# H_ij = log2(V_ij)
#
# Interpretation:
#
# V = 1  -> H = 0     exactly meets MAR; no headroom
# V = 2  -> H = 1     twice the MAR
# V = 4  -> H = 2     four times the MAR
# V = 8  -> H = 3     eight times the MAR
#
# Positive values therefore indicate increasing habitat
# availability beyond the estimated minimum area requirement.
#
# ==========================================================

headroom_species <- dat %>%
  
  mutate(
    
    # --------------------------------------------------------
    # Total viable habitat relative to MAR
    # --------------------------------------------------------
    
    total_MAR_ratio =
      total_viable_patch_area_km2 /
      MAR_km2,
    
    
    # --------------------------------------------------------
    # Largest individual viable patch relative to MAR
    # --------------------------------------------------------
    
    largest_patch_MAR_ratio =
      largest_patch_area_km2 /
      MAR_km2,
    
    
    # --------------------------------------------------------
    # Primary viability headroom metric
    # --------------------------------------------------------
    
    total_headroom_score =
      log2(
        total_MAR_ratio
      ),
    
    
    # --------------------------------------------------------
    # Sensitivity measure based on largest individual VP
    # --------------------------------------------------------
    
    largest_patch_headroom_score =
      log2(
        largest_patch_MAR_ratio
      )
    
  )


# ==========================================================
# 6. Validation
# ==========================================================

# ----------------------------------------------------------
# Check for non-finite MAR ratios
# ----------------------------------------------------------

if (
  any(
    !is.finite(
      headroom_species$total_MAR_ratio
    )
  )
) {
  
  stop(
    "Non-finite total MAR ratios were generated."
  )
}


# ----------------------------------------------------------
# Check for non-finite primary headroom scores
# ----------------------------------------------------------

if (
  any(
    !is.finite(
      headroom_species$total_headroom_score
    )
  )
) {
  
  stop(
    "Non-finite total headroom scores were generated."
  )
}


# ----------------------------------------------------------
# Check that supported species meet at least one MAR
# ----------------------------------------------------------

# Species in this table have already passed the viability
# threshold. Total MAR ratios substantially below 1 would
# therefore indicate a mismatch in the source data.

below_mar <- headroom_species %>%
  
  filter(
    total_MAR_ratio < 0.999
  )


if (nrow(below_mar) > 0) {
  
  warning(
    nrow(below_mar),
    " hotspot × species records have total MAR ratio < 1."
  )
  
  print(
    below_mar %>%
      select(
        hotspot_id,
        species,
        total_viable_patch_area_km2,
        MAR_km2,
        total_MAR_ratio,
        total_headroom_score
      )
  )
}


# ==========================================================
# 7. Create compact species-level headroom table
# ==========================================================

headroom_species_simple <- headroom_species %>%
  
  select(
    hotspot_id,
    country,
    hotspot_area_km2,
    species_richness,
    median_faith_PD,
    
    species,
    n_viable_patches,
    
    total_viable_patch_area_km2,
    largest_patch_area_km2,
    MAR_km2,
    
    total_MAR_ratio,
    total_headroom_score,
    
    largest_patch_MAR_ratio,
    largest_patch_headroom_score
  ) %>%
  
  arrange(
    hotspot_id,
    species
  )


# ==========================================================
# 8. Summarise viability headroom at hotspot level
# ==========================================================

# The primary hotspot-level metric is the mean species-level
# headroom across all species supported by the hotspot:
#
#                    1
# H_j = ----------------------------- sum(H_ij)
#       number of supported species
#
# Species with H = 0 are retained in this mean because zero is
# a meaningful value representing exactly one MAR.
#
# Species not supported by a hotspot are absent from the
# species-level table and therefore do not enter the mean.
#
# Additional summaries are retained for interpretation:
#
# SUM:
# Total transformed headroom across supported species.
#
# MEDIAN:
# Typical species-level headroom, less sensitive to extremes.
#
# MINIMUM:
# Headroom of the most spatially constrained supported species.
#
# ==========================================================

headroom_hotspot_summary <- headroom_species %>%
  
  group_by(
    hotspot_id
  ) %>%
  
  summarise(
    
    country =
      first(
        country
      ),
    
    hotspot_area_km2 =
      first(
        hotspot_area_km2
      ),
    
    species_richness =
      first(
        species_richness
      ),
    
    median_faith_PD =
      first(
        median_faith_PD
      ),
    
    
    # --------------------------------------------------------
    # Raw total-area MAR ratios
    # --------------------------------------------------------
    
    sum_total_MAR_ratio =
      sum(
        total_MAR_ratio,
        na.rm = TRUE
      ),
    
    mean_total_MAR_ratio =
      mean(
        total_MAR_ratio,
        na.rm = TRUE
      ),
    
    median_total_MAR_ratio =
      median(
        total_MAR_ratio,
        na.rm = TRUE
      ),
    
    minimum_total_MAR_ratio =
      min(
        total_MAR_ratio,
        na.rm = TRUE
      ),
    
    
    # --------------------------------------------------------
    # Primary viability headroom summaries
    # --------------------------------------------------------
    
    headroom_score_sum =
      sum(
        total_headroom_score,
        na.rm = TRUE
      ),
    
    headroom_score_mean =
      mean(
        total_headroom_score,
        na.rm = TRUE
      ),
    
    headroom_score_median =
      median(
        total_headroom_score,
        na.rm = TRUE
      ),
    
    headroom_score_minimum =
      min(
        total_headroom_score,
        na.rm = TRUE
      ),
    
    
    # --------------------------------------------------------
    # Largest-patch sensitivity summaries
    # --------------------------------------------------------
    
    largest_patch_headroom_sum =
      sum(
        largest_patch_headroom_score,
        na.rm = TRUE
      ),
    
    largest_patch_headroom_mean =
      mean(
        largest_patch_headroom_score,
        na.rm = TRUE
      ),
    
    largest_patch_headroom_median =
      median(
        largest_patch_headroom_score,
        na.rm = TRUE
      ),
    
    largest_patch_headroom_minimum =
      min(
        largest_patch_headroom_score,
        na.rm = TRUE
      ),
    
    
    # --------------------------------------------------------
    # Other descriptive information
    # --------------------------------------------------------
    
    total_number_viable_patches =
      sum(
        n_viable_patches,
        na.rm = TRUE
      ),
    
    total_species_viable_area_km2 =
      sum(
        total_viable_patch_area_km2,
        na.rm = TRUE
      ),
    
    .groups = "drop"
    
  ) %>%
  
  arrange(
    hotspot_id
  )


# ==========================================================
# 9. Identify most marginal and strongest species
# ==========================================================

marginal_species <- headroom_species %>%
  
  group_by(
    hotspot_id
  ) %>%
  
  slice_min(
    total_headroom_score,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup() %>%
  
  select(
    hotspot_id,
    weakest_species = species,
    weakest_species_MAR_ratio =
      total_MAR_ratio,
    weakest_species_headroom =
      total_headroom_score
  )


strongest_species <- headroom_species %>%
  
  group_by(
    hotspot_id
  ) %>%
  
  slice_max(
    total_headroom_score,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup() %>%
  
  select(
    hotspot_id,
    strongest_species = species,
    strongest_species_MAR_ratio =
      total_MAR_ratio,
    strongest_species_headroom =
      total_headroom_score
  )


headroom_hotspot_summary <-
  headroom_hotspot_summary %>%
  
  left_join(
    marginal_species,
    by = "hotspot_id"
  ) %>%
  
  left_join(
    strongest_species,
    by = "hotspot_id"
  )


# ==========================================================
# 10. Create species × hotspot MAR-ratio matrix
# ==========================================================

# A value of 0 in this output matrix indicates that the species
# is not supported by that hotspot. These matrix zeros are for
# presentation only and are NOT included in hotspot-level means.

headroom_matrix <- headroom_species %>%
  
  select(
    hotspot_id,
    species,
    total_MAR_ratio
  ) %>%
  
  pivot_wider(
    names_from = species,
    values_from = total_MAR_ratio,
    values_fill = 0
  ) %>%
  
  arrange(
    hotspot_id
  )


# ==========================================================
# 11. Create species × hotspot viability-headroom matrix
# ==========================================================

# A supported species may legitimately have H = 0 when its
# total viable habitat equals exactly one MAR.
#
# Unsupported species are also shown as 0 in this wide matrix,
# so this matrix should not be used to recalculate hotspot
# means. Hotspot means are calculated above directly from the
# species-level table.

headroom_score_matrix <- headroom_species %>%
  
  select(
    hotspot_id,
    species,
    total_headroom_score
  ) %>%
  
  pivot_wider(
    names_from = species,
    values_from = total_headroom_score,
    values_fill = 0
  ) %>%
  
  arrange(
    hotspot_id
  )


# ==========================================================
# 12. Save CSV outputs
# ==========================================================

write.csv(
  headroom_species_simple,
  file.path(
    out_dir,
    "Hotspot_Species_Viability_Headroom.csv"
  ),
  row.names = FALSE
)


write.csv(
  headroom_hotspot_summary,
  file.path(
    out_dir,
    "Hotspot_Viability_Headroom_Summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  headroom_matrix,
  file.path(
    out_dir,
    "Hotspot_Species_MAR_Ratio_Matrix.csv"
  ),
  row.names = FALSE
)


write.csv(
  headroom_score_matrix,
  file.path(
    out_dir,
    "Hotspot_Species_Headroom_Score_Matrix.csv"
  ),
  row.names = FALSE
)


# ==========================================================
# 13. Save Excel workbook
# ==========================================================

writexl::write_xlsx(
  
  list(
    
    Species_Headroom =
      headroom_species_simple,
    
    Hotspot_Summary =
      headroom_hotspot_summary,
    
    MAR_Ratio_Matrix =
      headroom_matrix,
    
    Headroom_Score_Matrix =
      headroom_score_matrix
    
  ),
  
  path = file.path(
    out_dir,
    "Rewilding_Hotspot_Viability_Headroom.xlsx"
  )
)


# ==========================================================
# 14. Print summary
# ==========================================================

message(
  "\n=============================================="
)

message(
  "VIABILITY HEADROOM ANALYSIS"
)

message(
  "==============================================\n"
)


print(
  headroom_hotspot_summary
)


# ==========================================================
# 15. Completion message
# ==========================================================

message(
  
  "\nViability headroom analysis complete.\n\n",
  
  "Primary hotspot-level output:\n",
  
  file.path(
    out_dir,
    "Hotspot_Viability_Headroom_Summary.csv"
  ),
  
  "\n\nSpecies-level output:\n",
  
  file.path(
    out_dir,
    "Hotspot_Species_Viability_Headroom.csv"
  ),
  
  "\n\nOutputs written to:\n",
  
  normalizePath(
    out_dir,
    mustWork = FALSE
  )
)