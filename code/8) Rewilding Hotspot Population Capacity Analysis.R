############################################################
# Identify and characterise rewilding hotspots using
# Faith's Phylogenetic Diversity and potential population
# capacity
############################################################
# Description - This script partitions Great Britain into 
# barrier-defined polygons and identifies the focal species 
# supported by each, assigning viable patches where at least 
# 95% of their area falls within a polygon. It summarises 
# patch numbers, viable habitat area, species-specific MARs 
# and MVPs, MVP equivalents, and potential population sizes. 
# Multi-species assemblages are matched to Faith’s phylogenetic 
# diversity values calculated across 100 dated VertLife mammal 
# phylogenies. Outputs include species-support and hotspot 
# summary tables, patch-level and population-capacity tables, 
# an unranked phylogenetic diversity map, and hotspot polygons 
# in GeoJSON format, linked by the permanent identifier 
# barrier_poly_id.


# Purpose:
#
# 1. Partition Great Britain into barrier-defined polygons
#
# 2. Determine which focal species are supported by each
#    barrier polygon.
#
#    A viable patch is assigned to a barrier polygon where
#    at least 95% of the complete patch area occurs within
#    that polygon.
#
# 3. Retain the individual viable patches assigned to each
#    hotspot and calculate:
#
#       - number of viable patches
#       - total viable patch area
#       - species-specific MAR
#       - species-specific MVP
#       - number of MVP equivalents
#       - potential population size
#
# 4. Match each multi-species assemblage to Faith's
#    Phylogenetic Diversity previously calculated across
#    100 dated VertLife mammal phylogenies.
#
# 5. Produce:
#
#       - complete species-support tables
#       - complete barrier-polygon Faith's PD summary
#       - multi-species hotspot summary
#       - patch-level hotspot-species table
#       - hotspot × species population-capacity table
#       - unranked Faith's PD hotspot map
#       - hotspot polygon GeoJSON
#
# IMPORTANT:
#
# No hotspot ranking is performed in this script.
#
# The permanent spatial identifier is:
#
#       barrier_poly_id
#
############################################################


# Clean session optional
# rm(list = ls())


# ==========================================================
# 1. Packages
# ==========================================================

required_packages <- c(
  "terra",
  "dplyr",
  "tidyr",
  "readxl",
  "writexl"
)

packages_to_install <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}

library(terra)
library(dplyr)
library(tidyr)
library(readxl)
library(writexl)


# ==========================================================
# 2. Helper functions
# ==========================================================


# ----------------------------------------------------------
# Convert SpatVectorCollection to SpatVector
# ----------------------------------------------------------

flatten_to_vect <- function(x) {
  
  if (inherits(x, "SpatVector")) {
    return(x)
  }
  
  if (inherits(x, "SpatVectorCollection")) {
    
    xs <- as.list(x)
    
    xs <- xs[
      !vapply(
        xs,
        is.null,
        logical(1)
      )
    ]
    
    if (length(xs) == 0) {
      return(NULL)
    }
    
    return(
      do.call(
        rbind,
        xs
      )
    )
  }
  
  stop(
    "Unknown object type: ",
    paste(
      class(x),
      collapse = ", "
    )
  )
}


# ----------------------------------------------------------
# Split one polygon using barrier geometries
# ----------------------------------------------------------

split_one_poly <- function(
    poly,
    splitter
) {
  
  out <- tryCatch(
    
    terra::split(
      poly,
      splitter
    ),
    
    error = function(e) {
      poly
    }
  )
  
  flatten_to_vect(out)
}


# ----------------------------------------------------------
# Check that a vector contains a required field
# ----------------------------------------------------------

assert_has_field <- function(
    v,
    field
) {
  
  if (!(field %in% names(v))) {
    
    stop(
      "Missing field '",
      field,
      "' in object. Fields are: ",
      paste(
        names(v),
        collapse = ", "
      )
    )
  }
}


# ----------------------------------------------------------
# Clean species names for matching
# ----------------------------------------------------------

clean_name_key <- function(x) {
  
  x <- gsub(
    "_",
    " ",
    x
  )
  
  x <- gsub(
    "\\s+",
    " ",
    x
  )
  
  tolower(
    trimws(x)
  )
}


# ==========================================================
# 3. Settings and paths
# ==========================================================

set.seed(13)

do_checks <- FALSE

tmp_dir <- "./_terra_tmp"

dir.create(
  tmp_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

terraOptions(
  memfrac = 0.5,
  tempdir = tmp_dir,
  progress = 1
)

path_barriers <-
  "./data/Combined Dispersal Barriers.geojson"

path_boundary <-
  "./data/Britain Boundary Polygon.geojson"

species_dir <-
  "./data/Suitable Habitats with Barriers/"

mar_file <-
  "./data/Species Minimum Area Required.xlsx"

faith_pd_lookup_file <- file.path(
  "./results/Phylogenetic Distances/",
  "Faith_PD_All_Possible_Assemblages_Summary.csv"
)

out_dir <-
  "./results/Rewilding Hotspots/Barrier Polygon Summary/"

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

crs_target <- "EPSG:27700"


# ----------------------------------------------------------
# Patch assignment threshold
# ----------------------------------------------------------
#
# A complete viable patch is assigned to a barrier polygon
# only where at least 95% of the patch lies within it.
#
# Once assigned, the COMPLETE patch area is used in the
# population-capacity calculation.

overlap_threshold <- 0.95

message(
  "Using viable-patch assignment threshold = ",
  overlap_threshold
)


# ----------------------------------------------------------
# Minimum assemblage size
# ----------------------------------------------------------
#
# A rewilding hotspot is retained where at least two focal
# species are supported.
#
# This is an eligibility condition ONLY.
# No richness-based ranking is performed.

minimum_assemblage_size <- 2


# ==========================================================
# 4. Species information
# ==========================================================

species_phylogeny_lookup <- data.frame(
  
  raster_species = c(
    "Beaver",
    "Bison",
    "Brown Bear",
    "Lynx",
    "Elk",
    "Pine Marten",
    "Wild Boar",
    "Wildcat",
    "Wolf",
    "Wolverine"
  ),
  
  common_name = c(
    "Eurasian beaver",
    "European bison",
    "European brown bear",
    "Northern lynx",
    "Eurasian elk",
    "Pine marten",
    "Wild boar",
    "European wildcat",
    "Eurasian wolf",
    "Wolverine"
  ),
  
  tree_species = c(
    "Castor fiber",
    "Bison bonasus",
    "Ursus arctos",
    "Lynx lynx",
    "Alces alces",
    "Martes martes",
    "Sus scrofa",
    "Felis silvestris",
    "Canis lupus",
    "Gulo gulo"
  ),
  
  focal_scientific_name = c(
    "Castor fiber",
    "Bison bonasus",
    "Ursus arctos arctos",
    "Lynx lynx",
    "Alces alces",
    "Martes martes",
    "Sus scrofa",
    "Felis silvestris silvestris",
    "Canis lupus lupus",
    "Gulo gulo"
  ),
  
  stringsAsFactors = FALSE
)


# ----------------------------------------------------------
# Species-specific MVP values
# ----------------------------------------------------------
#
# These are the MVP values corresponding to the MAR values
# used to define the viable patches.

mvp_lookup <- data.frame(
  
  raster_species = c(
    "Beaver",
    "Bison",
    "Brown Bear",
    "Lynx",
    "Elk",
    "Pine Marten",
    "Wild Boar",
    "Wildcat",
    "Wolf",
    "Wolverine"
  ),
  
  MVP = c(
    42,
    263,
    40,
    15,
    50,
    40,
    300,
    40,
    100,
    46
  ),
  
  stringsAsFactors = FALSE
)


# ----------------------------------------------------------
# Flexible aliases for species matching
# ----------------------------------------------------------

species_alias_lookup <- bind_rows(
  
  # Short raster names
  species_phylogeny_lookup %>%
    transmute(
      species_alias = raster_species,
      raster_species,
      common_name,
      tree_species,
      focal_scientific_name
    ),
  
  # Tree species names
  species_phylogeny_lookup %>%
    transmute(
      species_alias = tree_species,
      raster_species,
      common_name,
      tree_species,
      focal_scientific_name
    ),
  
  # Full focal scientific names
  species_phylogeny_lookup %>%
    transmute(
      species_alias = focal_scientific_name,
      raster_species,
      common_name,
      tree_species,
      focal_scientific_name
    ),
  
  # Full common names
  species_phylogeny_lookup %>%
    transmute(
      species_alias = common_name,
      raster_species,
      common_name,
      tree_species,
      focal_scientific_name
    )
  
) %>%
  
  mutate(
    species_key =
      clean_name_key(
        species_alias
      )
  ) %>%
  
  distinct(
    species_key,
    .keep_all = TRUE
  )


# Keep exactly the phylogenetic species ordering used when
# the Faith's PD assemblage lookup was generated.

tree_species_order <-
  species_phylogeny_lookup$tree_species


# ==========================================================
# 5. Load and validate MAR data
# ==========================================================

if (!file.exists(mar_file)) {
  
  stop(
    "MAR file could not be found:\n",
    mar_file
  )
}

mar_data <- readxl::read_excel(
  mar_file
)

required_mar_fields <- c(
  "species",
  "min_area (km2)",
  "min_area (m2)"
)

missing_mar_fields <- setdiff(
  required_mar_fields,
  names(mar_data)
)

if (length(missing_mar_fields) > 0) {
  
  stop(
    paste0(
      "MAR table is missing the following field(s):\n",
      paste(
        missing_mar_fields,
        collapse = "\n"
      )
    )
  )
}


# Normalise MAR species names.

mar_data <- mar_data %>%
  
  mutate(
    
    species_key =
      clean_name_key(
        species
      ),
    
    MAR_km2 =
      as.numeric(
        `min_area (km2)`
      ),
    
    MAR_m2 =
      as.numeric(
        `min_area (m2)`
      )
  )


# Match MAR names to standard raster species names.

mar_alias <- species_alias_lookup %>%
  
  select(
    species_key,
    raster_species
  ) %>%
  
  distinct()


mar_data <- mar_data %>%
  
  left_join(
    mar_alias,
    by = "species_key"
  )


if (any(is.na(mar_data$raster_species))) {
  
  print(
    mar_data %>%
      filter(
        is.na(
          raster_species
        )
      )
  )
  
  stop(
    "One or more species in the MAR table could not be ",
    "matched to the species lookup."
  )
}


# Join MVP values.

species_population_lookup <- mar_data %>%
  
  select(
    raster_species,
    MAR_km2,
    MAR_m2
  ) %>%
  
  left_join(
    mvp_lookup,
    by = "raster_species"
  ) %>%
  
  left_join(
    species_phylogeny_lookup,
    by = "raster_species"
  ) %>%
  
  select(
    raster_species,
    common_name,
    focal_scientific_name,
    tree_species,
    MAR_km2,
    MAR_m2,
    MVP
  )


if (
  any(
    is.na(
      species_population_lookup$MVP
    )
  )
) {
  
  stop(
    "At least one species does not have an MVP value."
  )
}


if (
  any(
    !is.finite(
      species_population_lookup$MAR_km2
    ) |
    species_population_lookup$MAR_km2 <= 0
  )
) {
  
  stop(
    "All MAR values must be finite and greater than zero."
  )
}


message(
  "Species MAR/MVP lookup successfully created."
)

print(
  species_population_lookup
)


# ==========================================================
# 6. Check remaining input files
# ==========================================================

stopifnot(
  file.exists(
    path_barriers
  )
)

stopifnot(
  file.exists(
    path_boundary
  )
)

stopifnot(
  dir.exists(
    species_dir
  )
)

if (
  !file.exists(
    faith_pd_lookup_file
  )
) {
  
  stop(
    paste0(
      "Faith's PD lookup table could not be found:\n",
      faith_pd_lookup_file,
      "\n\nRun the phylogenetic analysis first."
    )
  )
}

message(
  "Using pre-calculated Faith's PD lookup:\n",
  normalizePath(
    faith_pd_lookup_file
  )
)


# ==========================================================
# 7. Load Britain boundary and countries
# ==========================================================

boundary_vect <- terra::vect(
  path_boundary
)

boundary_vect <- terra::project(
  boundary_vect,
  crs_target
)

boundary_vect <- terra::makeValid(
  boundary_vect
)

assert_has_field(
  boundary_vect,
  "shape1"
)

regions <- list(
  
  England =
    boundary_vect[
      boundary_vect$shape1 ==
        "England",
    ],
  
  Scotland =
    boundary_vect[
      boundary_vect$shape1 ==
        "Scotland",
    ],
  
  Wales =
    boundary_vect[
      boundary_vect$shape1 ==
        "Wales",
    ]
)


# Dissolve internal country borders.

boundary_union <- terra::aggregate(
  boundary_vect
)

boundary_union <- terra::makeValid(
  boundary_union
)


# ==========================================================
# 8. Load dispersal barriers
# ==========================================================

barrier_lines_raw <- terra::vect(
  path_barriers
)

barrier_lines_raw <- terra::project(
  barrier_lines_raw,
  crs_target
)

barrier_lines_raw <- terra::makeValid(
  barrier_lines_raw
)

assert_has_field(
  barrier_lines_raw,
  "type"
)


# Add coastline.

coast_line <- terra::as.lines(
  boundary_union
)

coast_line$type <- "coast"

barrier_lines <- rbind(
  barrier_lines_raw,
  coast_line
)

barrier_lines <- terra::makeValid(
  barrier_lines
)


# ==========================================================
# 9. Split Britain into barrier polygons
# ==========================================================

barrier_polys <- boundary_union

non_rail_types <- setdiff(
  unique(
    barrier_lines$type
  ),
  "rail"
)


for (barrier_type in non_rail_types) {
  
  message(
    "Splitting by: ",
    barrier_type
  )
  
  current_barrier <- barrier_lines[
    barrier_lines$type ==
      barrier_type,
  ]
  
  if (
    nrow(
      current_barrier
    ) == 0
  ) {
    next
  }
  
  current_barrier <- terra::makeValid(
    current_barrier
  )
  
  barrier_polys <- flatten_to_vect(
    barrier_polys
  )
  
  new_polys <- vector(
    "list",
    nrow(
      barrier_polys
    )
  )
  
  for (
    polygon_number in
    seq_len(
      nrow(
        barrier_polys
      )
    )
  ) {
    
    new_polys[[polygon_number]] <-
      split_one_poly(
        barrier_polys[
          polygon_number,
        ],
        current_barrier
      )
  }
  
  barrier_polys <- do.call(
    rbind,
    new_polys
  )
}


# Railways represented using a 25 m buffer.

rail_lines <- barrier_lines[
  barrier_lines$type ==
    "rail",
]

if (
  nrow(
    rail_lines
  ) > 0
) {
  
  message(
    "Splitting by: rail (buffered)"
  )
  
  rail_lines <- terra::makeValid(
    rail_lines
  )
  
  rail_buffer <- terra::buffer(
    rail_lines,
    width = 25
  )
  
  barrier_polys <- flatten_to_vect(
    barrier_polys
  )
  
  new_polys <- vector(
    "list",
    nrow(
      barrier_polys
    )
  )
  
  for (
    polygon_number in
    seq_len(
      nrow(
        barrier_polys
      )
    )
  ) {
    
    new_polys[[polygon_number]] <-
      split_one_poly(
        barrier_polys[
          polygon_number,
        ],
        rail_buffer
      )
  }
  
  barrier_polys <- do.call(
    rbind,
    new_polys
  )
}


barrier_polys <- terra::makeValid(
  barrier_polys
)


message(
  "Barrier polygon geometry type: ",
  terra::geomtype(
    barrier_polys
  )
)


if (
  terra::geomtype(
    barrier_polys
  ) != "polygons"
) {
  
  stop(
    "barrier_polys is not polygon geometry. ",
    "Check the barrier input data."
  )
}


# ----------------------------------------------------------
# Assign permanent IDs and calculate polygon areas
# ----------------------------------------------------------

barrier_polys$barrier_poly_id <-
  seq_len(
    nrow(
      barrier_polys
    )
  )


barrier_polys$area_m2 <-
  terra::expanse(
    barrier_polys,
    unit = "m"
  )


barrier_polys$area_km2 <-
  barrier_polys$area_m2 /
  1e6


message(
  "Barrier polygons created: ",
  nrow(
    barrier_polys
  )
)


summary(
  barrier_polys$area_km2
)


if (do_checks) {
  
  plot(
    barrier_polys
  )
  
  plot(
    boundary_union,
    add = TRUE,
    border = "black",
    lwd = 1.5
  )
}


# ==========================================================
# 10. Assign each polygon to one country
# ==========================================================

country_assignments <- lapply(
  
  names(
    regions
  ),
  
  function(country_name) {
    
    intersection <- terra::intersect(
      barrier_polys,
      regions[[country_name]]
    )
    
    if (
      is.null(
        intersection
      ) ||
      nrow(
        intersection
      ) == 0
    ) {
      
      return(NULL)
    }
    
    if (
      !(
        "barrier_poly_id" %in%
        names(
          intersection
        )
      )
    ) {
      
      stop(
        "barrier_poly_id was not retained after ",
        "intersection for ",
        country_name,
        "."
      )
    }
    
    intersection$country <-
      country_name
    
    intersection$overlap_km2 <-
      terra::expanse(
        intersection,
        unit = "km"
      )
    
    intersection[
      ,
      c(
        "barrier_poly_id",
        "country",
        "overlap_km2"
      )
    ]
  }
)


country_assignments <- do.call(
  rbind,
  country_assignments
)


country_df <- as.data.frame(
  country_assignments
) %>%
  
  group_by(
    barrier_poly_id
  ) %>%
  
  slice_max(
    overlap_km2,
    n = 1,
    with_ties = FALSE
  ) %>%
  
  ungroup() %>%
  
  select(
    barrier_poly_id,
    country
  )


# ==========================================================
# 11. Locate final viable-patch rasters
# ==========================================================

species_tifs <- sort(
  
  list.files(
    species_dir,
    pattern =
      "Suitable Habitat Split and Filtered\\.tif$",
    full.names = TRUE
  )
)


if (length(species_tifs) == 0) {
  
  stop(
    "No final species rasters were found in: ",
    species_dir
  )
}


message(
  "Final species rasters found: ",
  length(
    species_tifs
  )
)


# ==========================================================
# 12. Determine species support AND retain assigned patches
# ==========================================================
#
# Two sets of information are retained:
#
# support_tables:
#     one row per species × barrier polygon giving the
#     existing yes/no species-support result.
#
# assigned_patch_tables:
#     one row per qualifying viable patch assigned to a
#     barrier polygon.
#
# A qualifying patch is one for which:
#
#       overlap_prop >= 0.95
#
# The COMPLETE patch area is retained for the subsequent
# population-capacity calculation.

support_tables <- list()

assigned_patch_tables <- list()


for (tif_path in species_tifs) {
  
  
  # --------------------------------------------------------
  # Species name from raster
  # --------------------------------------------------------
  
  raster_species_name <- gsub(
    " Suitable Habitat Split and Filtered\\.tif$",
    "",
    basename(
      tif_path
    )
  )
  
  raster_species_name <- gsub(
    "_",
    " ",
    raster_species_name
  )
  
  raster_species_name <- trimws(
    raster_species_name
  )
  
  
  message(
    "\n--- Species: ",
    raster_species_name,
    " ---"
  )
  
  
  # --------------------------------------------------------
  # Match species identity
  # --------------------------------------------------------
  
  species_key <- clean_name_key(
    raster_species_name
  )
  
  
  lookup_row <- species_alias_lookup[
    species_alias_lookup$species_key ==
      species_key,
    ,
    drop = FALSE
  ]
  
  
  if (nrow(lookup_row) == 0) {
    
    stop(
      "The raster species name '",
      raster_species_name,
      "' could not be matched to the species lookup."
    )
  }
  
  
  if (nrow(lookup_row) > 1) {
    
    stop(
      "The raster species name '",
      raster_species_name,
      "' matched more than one species."
    )
  }
  
  
  raster_species_standard <-
    lookup_row$raster_species[1]
  
  common_name <-
    lookup_row$common_name[1]
  
  tree_species_name <-
    lookup_row$tree_species[1]
  
  focal_scientific_name_current <-
    lookup_row$focal_scientific_name[1]
  
  
  # --------------------------------------------------------
  # Load patch-ID raster
  # --------------------------------------------------------
  
  species_raster <- terra::rast(
    tif_path
  )
  
  
  if (
    is.na(
      terra::crs(
        species_raster
      )
    )
  ) {
    
    stop(
      "Raster has no CRS: ",
      tif_path
    )
  }
  
  
  if (
    !terra::same.crs(
      species_raster,
      crs_target
    )
  ) {
    
    species_raster <- terra::project(
      species_raster,
      crs_target,
      method = "near"
    )
  }
  
  
  # --------------------------------------------------------
  # Check whether any viable patches remain
  # --------------------------------------------------------
  
  raster_non_na <- terra::global(
    !is.na(
      species_raster
    ),
    "sum",
    na.rm = TRUE
  )[1, 1]
  
  
  if (
    is.na(
      raster_non_na
    ) ||
    raster_non_na == 0
  ) {
    
    message(
      "  No viable patches; species skipped."
    )
    
    rm(
      species_raster
    )
    
    gc()
    
    next
  }
  
  
  # --------------------------------------------------------
  # Convert patch IDs to dissolved polygons
  # --------------------------------------------------------
  
  patch_polys <- terra::as.polygons(
    species_raster,
    dissolve = TRUE,
    na.rm = TRUE
  )
  
  
  names(
    patch_polys
  ) <- "patch_id"
  
  
  patch_polys <- terra::makeValid(
    patch_polys
  )
  
  
  # COMPLETE viable-patch area.
  
  patch_polys$patch_area_m2 <-
    terra::expanse(
      patch_polys,
      unit = "m"
    )
  
  
  # --------------------------------------------------------
  # Intersect viable patches with barrier polygons
  # --------------------------------------------------------
  
  pieces <- terra::intersect(
    patch_polys,
    barrier_polys
  )
  
  
  if (
    is.null(
      pieces
    ) ||
    nrow(
      pieces
    ) == 0
  ) {
    
    message(
      "  No intersections with barrier polygons."
    )
    
    rm(
      species_raster,
      patch_polys,
      pieces
    )
    
    gc()
    
    next
  }
  
  
  required_piece_fields <- c(
    "barrier_poly_id",
    "patch_id",
    "patch_area_m2"
  )
  
  
  missing_piece_fields <- setdiff(
    required_piece_fields,
    names(
      pieces
    )
  )
  
  
  if (
    length(
      missing_piece_fields
    ) > 0
  ) {
    
    stop(
      "Missing fields after intersection for ",
      raster_species_standard,
      ": ",
      paste(
        missing_piece_fields,
        collapse = ", "
      )
    )
  }
  
  
  pieces <- terra::makeValid(
    pieces
  )
  
  
  # Area of the portion actually occurring within the
  # barrier polygon.
  
  pieces$piece_area_m2 <-
    terra::expanse(
      pieces,
      unit = "m"
    )
  
  
  pieces_df <- as.data.frame(
    pieces
  )[
    ,
    c(
      "barrier_poly_id",
      "patch_id",
      "patch_area_m2",
      "piece_area_m2"
    )
  ]
  
  
  # --------------------------------------------------------
  # Calculate patch assignment proportions
  # --------------------------------------------------------
  
  patch_in_barrier <- pieces_df %>%
    
    group_by(
      barrier_poly_id,
      patch_id
    ) %>%
    
    summarise(
      
      patch_area_m2 =
        first(
          patch_area_m2
        ),
      
      overlap_area_m2 =
        sum(
          piece_area_m2,
          na.rm = TRUE
        ),
      
      .groups = "drop"
    ) %>%
    
    mutate(
      
      overlap_prop = ifelse(
        patch_area_m2 > 0,
        overlap_area_m2 /
          patch_area_m2,
        NA_real_
      )
    )
  
  
  # --------------------------------------------------------
  # Retain qualifying assigned patches
  # --------------------------------------------------------
  
  assigned_patches <- patch_in_barrier %>%
    
    filter(
      is.finite(
        overlap_prop
      ),
      overlap_prop >=
        overlap_threshold
    ) %>%
    
    mutate(
      
      raster_species =
        raster_species_standard,
      
      common_name =
        common_name,
      
      focal_scientific_name =
        focal_scientific_name_current,
      
      tree_species =
        tree_species_name,
      
      patch_area_km2 =
        patch_area_m2 /
        1e6,
      
      overlap_area_km2 =
        overlap_area_m2 /
        1e6
    ) %>%
    
    select(
      barrier_poly_id,
      raster_species,
      common_name,
      focal_scientific_name,
      tree_species,
      patch_id,
      patch_area_m2,
      patch_area_km2,
      overlap_area_m2,
      overlap_area_km2,
      overlap_prop
    )
  
  
  if (
    nrow(
      assigned_patches
    ) > 0
  ) {
    
    assigned_patch_tables[[
      tree_species_name
    ]] <- assigned_patches
  }
  
  
  # --------------------------------------------------------
  # Existing species-support calculation
  # --------------------------------------------------------
  
  supported_by_polygon <- patch_in_barrier %>%
    
    group_by(
      barrier_poly_id
    ) %>%
    
    summarise(
      
      raster_species =
        raster_species_standard,
      
      common_name =
        common_name,
      
      focal_scientific_name =
        focal_scientific_name_current,
      
      tree_species =
        tree_species_name,
      
      supports_species = any(
        overlap_prop >=
          overlap_threshold,
        na.rm = TRUE
      ),
      
      best_overlap_prop = if (
        all(
          is.na(
            overlap_prop
          )
        )
      ) {
        NA_real_
      } else {
        max(
          overlap_prop,
          na.rm = TRUE
        )
      },
      
      best_patch_id = if (
        all(
          is.na(
            overlap_prop
          )
        )
      ) {
        NA
      } else {
        patch_id[
          which.max(
            overlap_prop
          )
        ][1]
      },
      
      .groups = "drop"
    )
  
  
  support_tables[[
    tree_species_name
  ]] <- supported_by_polygon
  
  
  rm(
    species_raster,
    patch_polys,
    pieces,
    pieces_df,
    patch_in_barrier,
    assigned_patches,
    supported_by_polygon
  )
  
  gc()
}


if (length(support_tables) == 0) {
  
  stop(
    "No species produced support results. ",
    "Check the species rasters and overlap threshold."
  )
}


if (length(assigned_patch_tables) == 0) {
  
  stop(
    "No viable patches satisfied the hotspot assignment ",
    "threshold."
  )
}


# ==========================================================
# 13. Combine species support and assigned patches
# ==========================================================

support_long <- bind_rows(
  support_tables
)


assigned_patches_long <- bind_rows(
  assigned_patch_tables
)


# ----------------------------------------------------------
# Safety check:
# each species patch should be assigned to only one barrier
# polygon under the >=95% rule.
# ----------------------------------------------------------

duplicate_patch_assignments <-
  assigned_patches_long %>%
  
  count(
    tree_species,
    patch_id,
    name = "n_polygon_assignments"
  ) %>%
  
  filter(
    n_polygon_assignments > 1
  )


if (
  nrow(
    duplicate_patch_assignments
  ) > 0
) {
  
  print(
    duplicate_patch_assignments
  )
  
  stop(
    "At least one viable patch was assigned to more than ",
    "one barrier polygon. Inspect the barrier geometries or ",
    "the overlap threshold before continuing."
  )
}


# ----------------------------------------------------------
# Create wide support table
# ----------------------------------------------------------

support_wide <- support_long %>%
  
  select(
    barrier_poly_id,
    tree_species,
    supports_species
  ) %>%
  
  tidyr::pivot_wider(
    
    names_from =
      tree_species,
    
    values_from =
      supports_species,
    
    values_fill =
      FALSE
  )


# Ensure all barrier polygons are represented.

support_wide <- as.data.frame(
  barrier_polys
)[
  ,
  "barrier_poly_id",
  drop = FALSE
] %>%
  
  left_join(
    support_wide,
    by = "barrier_poly_id"
  ) %>%
  
  mutate(
    across(
      where(
        is.logical
      ),
      ~ tidyr::replace_na(
        .x,
        FALSE
      )
    )
  )


species_support_columns <- intersect(
  names(
    support_wide
  ),
  tree_species_order
)


if (
  length(
    species_support_columns
  ) == 0
) {
  
  stop(
    "No phylogenetic species-support columns were created."
  )
}


# ==========================================================
# 14. Build supported assemblage for each polygon
# ==========================================================

supported_species_lists <- lapply(
  
  seq_len(
    nrow(
      support_wide
    )
  ),
  
  function(row_number) {
    
    support_values <- unlist(
      
      support_wide[
        row_number,
        species_support_columns,
        drop = FALSE
      ],
      
      use.names = TRUE
    )
    
    
    support_values <- as.logical(
      support_values
    )
    
    
    species_support_columns[
      support_values
    ]
  }
)


supported_species_text <- vapply(
  
  supported_species_lists,
  
  function(selected_species) {
    
    if (
      length(
        selected_species
      ) == 0
    ) {
      
      return("")
    }
    
    
    paste(
      selected_species,
      collapse = "; "
    )
  },
  
  character(1)
)


# ----------------------------------------------------------
# Create phylogenetic assemblage key
# ----------------------------------------------------------

make_assemblage_key <- function(
    selected_species
) {
  
  selected_species <-
    tree_species_order[
      tree_species_order %in%
        unique(
          selected_species
        )
    ]
  
  
  paste(
    selected_species,
    collapse = " | "
  )
}


assemblage_keys <- vapply(
  supported_species_lists,
  make_assemblage_key,
  character(1)
)


assemblage_richness <- vapply(
  supported_species_lists,
  length,
  integer(1)
)


is_multispecies_assemblage <- vapply(
  supported_species_lists,
  length,
  integer(1)
) >= minimum_assemblage_size


# ==========================================================
# 15. Create barrier-polygon summary
# ==========================================================

hotspot_summary <- support_wide %>%
  
  left_join(
    
    as.data.frame(
      barrier_polys
    )[
      ,
      c(
        "barrier_poly_id",
        "area_km2"
      )
    ],
    
    by = "barrier_poly_id"
  ) %>%
  
  left_join(
    country_df,
    by = "barrier_poly_id"
  ) %>%
  
  mutate(
    
    supported_species =
      supported_species_text,
    
    assemblage_key =
      assemblage_keys,
    
    species_richness =
      assemblage_richness,
    
    is_multispecies_assemblage =
      is_multispecies_assemblage
  ) %>%
  
  relocate(
    barrier_poly_id,
    country,
    area_km2,
    species_richness,
    supported_species,
    assemblage_key,
    is_multispecies_assemblage
  )


# ==========================================================
# 16. Load pre-calculated Faith's PD lookup
# ==========================================================

faith_pd_lookup <- read.csv(
  faith_pd_lookup_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


required_pd_fields <- c(
  "assemblage_key",
  "species_richness",
  "median_faith_PD",
  "lower_95_faith_PD",
  "upper_95_faith_PD",
  "median_relative_PD",
  "lower_95_relative_PD",
  "upper_95_relative_PD",
  "median_PD_percentile",
  "median_SES_PD"
)


missing_pd_fields <- setdiff(
  required_pd_fields,
  names(
    faith_pd_lookup
  )
)


if (
  length(
    missing_pd_fields
  ) > 0
) {
  
  stop(
    paste0(
      "Faith's PD lookup is missing required fields:\n",
      paste(
        missing_pd_fields,
        collapse = "\n"
      )
    )
  )
}


if (
  anyDuplicated(
    faith_pd_lookup$assemblage_key
  )
) {
  
  stop(
    "Faith's PD lookup contains duplicate assemblage keys."
  )
}


message(
  "Faith's PD assemblages available: ",
  nrow(
    faith_pd_lookup
  )
)


if (
  nrow(
    faith_pd_lookup
  ) != 1013
) {
  
  warning(
    "Expected 1013 possible assemblages but found ",
    nrow(
      faith_pd_lookup
    ),
    "."
  )
}


# ==========================================================
# 17. Attach Faith's PD to barrier polygons
# ==========================================================

pd_fields_to_join <- faith_pd_lookup %>%
  
  select(
    assemblage_key,
    median_faith_PD,
    lower_95_faith_PD,
    upper_95_faith_PD,
    median_relative_PD,
    lower_95_relative_PD,
    upper_95_relative_PD,
    median_PD_percentile,
    median_SES_PD
  )


hotspot_summary <- hotspot_summary %>%
  
  left_join(
    pd_fields_to_join,
    by = "assemblage_key"
  )


# ----------------------------------------------------------
# Verify every multi-species assemblage was matched
# ----------------------------------------------------------

unmatched_multispecies <- hotspot_summary %>%
  
  filter(
    is_multispecies_assemblage,
    !is.finite(
      median_faith_PD
    )
  )


if (
  nrow(
    unmatched_multispecies
  ) > 0
) {
  
  print(
    
    unmatched_multispecies %>%
      
      select(
        barrier_poly_id,
        supported_species,
        assemblage_key
      )
  )
  
  
  stop(
    paste0(
      nrow(
        unmatched_multispecies
      ),
      " multi-species barrier polygon(s) could not be ",
      "matched to the pre-calculated Faith's PD lookup."
    )
  )
}


hotspot_summary <- hotspot_summary %>%
  
  relocate(
    median_faith_PD,
    lower_95_faith_PD,
    upper_95_faith_PD,
    median_relative_PD,
    lower_95_relative_PD,
    upper_95_relative_PD,
    median_PD_percentile,
    median_SES_PD,
    .after = assemblage_key
  )


message(
  "All multi-species assemblages successfully matched ",
  "to pre-calculated Faith's PD."
)


# ==========================================================
# 18. Identify multi-species rewilding hotspots
# ==========================================================
#
# There is NO ranking at this stage.
#
# barrier_poly_id remains the permanent hotspot identifier.

rewilding_hotspots <- hotspot_summary %>%
  
  filter(
    is_multispecies_assemblage,
    is.finite(
      median_faith_PD
    )
  ) %>%
  
  arrange(
    barrier_poly_id
  )


if (
  nrow(
    rewilding_hotspots
  ) == 0
) {
  
  stop(
    "No barrier polygons supported a multi-species assemblage."
  )
}


message(
  "\nMulti-species rewilding hotspots identified: ",
  nrow(
    rewilding_hotspots
  )
)


# ==========================================================
# 19. Retain patches belonging to rewilding hotspots
# ==========================================================

hotspot_ids <-
  rewilding_hotspots$barrier_poly_id


hotspot_assigned_patches <- assigned_patches_long %>%
  
  filter(
    barrier_poly_id %in%
      hotspot_ids
  )


# Safety check:
# Every assigned-patch record should correspond to a species
# that is actually recorded as supported by that polygon.

supported_pairs <- support_long %>%
  
  filter(
    supports_species
  ) %>%
  
  distinct(
    barrier_poly_id,
    tree_species
  )


invalid_patch_pairs <- hotspot_assigned_patches %>%
  
  anti_join(
    supported_pairs,
    by = c(
      "barrier_poly_id",
      "tree_species"
    )
  )


if (
  nrow(
    invalid_patch_pairs
  ) > 0
) {
  
  print(
    invalid_patch_pairs
  )
  
  stop(
    "Assigned patch records were found for species not ",
    "recorded as supported in the same barrier polygon."
  )
}


# ==========================================================
# 20. Create hotspot × species viable-patch summary
# ==========================================================

hotspot_species_patch_summary <-
  hotspot_assigned_patches %>%
  
  group_by(
    barrier_poly_id,
    raster_species,
    common_name,
    focal_scientific_name,
    tree_species
  ) %>%
  
  summarise(
    
    n_viable_patches =
      n_distinct(
        patch_id
      ),
    
    total_viable_patch_area_m2 =
      sum(
        patch_area_m2,
        na.rm = TRUE
      ),
    
    total_viable_patch_area_km2 =
      sum(
        patch_area_km2,
        na.rm = TRUE
      ),
    
    total_overlap_area_km2 =
      sum(
        overlap_area_km2,
        na.rm = TRUE
      ),
    
    smallest_patch_area_km2 =
      min(
        patch_area_km2,
        na.rm = TRUE
      ),
    
    largest_patch_area_km2 =
      max(
        patch_area_km2,
        na.rm = TRUE
      ),
    
    mean_patch_area_km2 =
      mean(
        patch_area_km2,
        na.rm = TRUE
      ),
    
    .groups = "drop"
  )


# ==========================================================
# 21. Add MAR and MVP values
# ==========================================================

hotspot_species_population <-
  hotspot_species_patch_summary %>%
  
  left_join(
    
    species_population_lookup %>%
      
      select(
        raster_species,
        MAR_km2,
        MAR_m2,
        MVP
      ),
    
    by = "raster_species"
  )


# Verify all species matched.

missing_population_parameters <-
  hotspot_species_population %>%
  
  filter(
    !is.finite(
      MAR_km2
    ) |
      !is.finite(
        MVP
      )
  )


if (
  nrow(
    missing_population_parameters
  ) > 0
) {
  
  print(
    missing_population_parameters
  )
  
  stop(
    "MAR or MVP values are missing for at least one ",
    "hotspot-species combination."
  )
}


# ==========================================================
# 22. Calculate MVP equivalents and potential population
# ==========================================================
#
# MAR_km2 represents the habitat area required for one MVP.
#
# Therefore:
#
# MVP equivalents =
#
#     total viable habitat area / MAR
#
#
# Potential population =
#
#     MVP equivalents × MVP
#
#
# The unrounded population estimate is retained because it
# should be used in any subsequent statistical/ranking
# analysis.
#
# A rounded whole-individual value is supplied separately
# for presentation.

hotspot_species_population <-
  hotspot_species_population %>%
  
  mutate(
    
    MVP_equivalents =
      total_viable_patch_area_km2 /
      MAR_km2,
    
    potential_population =
      MVP_equivalents *
      MVP,
    
    potential_population_rounded =
      round(
        potential_population
      )
  )


# ==========================================================
# 23. Attach hotspot-level information
# ==========================================================

hotspot_level_fields <-
  rewilding_hotspots %>%
  
  select(
    barrier_poly_id,
    country,
    area_km2,
    species_richness,
    supported_species,
    assemblage_key,
    median_faith_PD,
    lower_95_faith_PD,
    upper_95_faith_PD,
    median_relative_PD,
    median_PD_percentile,
    median_SES_PD
  ) %>%
  
  rename(
    hotspot_id =
      barrier_poly_id,
    
    hotspot_area_km2 =
      area_km2
  )


hotspot_species_population <-
  hotspot_species_population %>%
  
  rename(
    hotspot_id =
      barrier_poly_id,
    
    species =
      common_name,
    
    scientific_name =
      focal_scientific_name
  ) %>%
  
  left_join(
    hotspot_level_fields,
    by = "hotspot_id"
  ) %>%
  
  select(
    
    hotspot_id,
    country,
    hotspot_area_km2,
    species_richness,
    
    median_faith_PD,
    lower_95_faith_PD,
    upper_95_faith_PD,
    median_relative_PD,
    median_PD_percentile,
    median_SES_PD,
    
    species,
    scientific_name,
    raster_species,
    tree_species,
    
    n_viable_patches,
    total_viable_patch_area_km2,
    
    MAR_km2,
    MVP,
    MVP_equivalents,
    potential_population,
    potential_population_rounded,
    
    smallest_patch_area_km2,
    largest_patch_area_km2,
    mean_patch_area_km2,
    
    total_overlap_area_km2,
    
    supported_species,
    assemblage_key
  ) %>%
  
  arrange(
    hotspot_id,
    species
  )


# ==========================================================
# 24. Create simpler publication/inspection table
# ==========================================================
#
# This is the compact table we discussed.
#
# One row =
#
#     one hotspot × one supported species

hotspot_species_population_simple <-
  hotspot_species_population %>%
  
  select(
    hotspot_id,
    hotspot_area_km2,
    species_richness,
    median_faith_PD,
    species,
    n_viable_patches,
    total_viable_patch_area_km2,
    MAR_km2,
    MVP,
    MVP_equivalents,
    potential_population,
    potential_population_rounded
  )


# ==========================================================
# 25. Create patch-level hotspot table
# ==========================================================
#
# This retains the individual patch IDs underlying the
# summary above and is useful for checking/reproducibility.

hotspot_patch_level_table <-
  hotspot_assigned_patches %>%
  
  rename(
    hotspot_id =
      barrier_poly_id,
    
    species =
      common_name,
    
    scientific_name =
      focal_scientific_name
  ) %>%
  
  left_join(
    
    hotspot_level_fields %>%
      select(
        hotspot_id,
        country,
        hotspot_area_km2,
        species_richness,
        median_faith_PD
      ),
    
    by = "hotspot_id"
  ) %>%
  
  left_join(
    
    species_population_lookup %>%
      select(
        raster_species,
        MAR_km2,
        MVP
      ),
    
    by = "raster_species"
  ) %>%
  
  mutate(
    
    patch_MVP_equivalents =
      patch_area_km2 /
      MAR_km2,
    
    patch_potential_population =
      patch_MVP_equivalents *
      MVP
  ) %>%
  
  select(
    hotspot_id,
    country,
    hotspot_area_km2,
    species_richness,
    median_faith_PD,
    species,
    scientific_name,
    raster_species,
    tree_species,
    patch_id,
    patch_area_km2,
    overlap_area_km2,
    overlap_prop,
    MAR_km2,
    MVP,
    patch_MVP_equivalents,
    patch_potential_population
  ) %>%
  
  arrange(
    hotspot_id,
    species,
    patch_id
  )


# ==========================================================
# 26. Validation checks
# ==========================================================


# ----------------------------------------------------------
# Check n_viable_patches against the patch-level table
# ----------------------------------------------------------

patch_count_check <-
  hotspot_patch_level_table %>%
  
  count(
    hotspot_id,
    species,
    name = "calculated_patch_count"
  ) %>%
  
  left_join(
    
    hotspot_species_population_simple %>%
      select(
        hotspot_id,
        species,
        n_viable_patches
      ),
    
    by = c(
      "hotspot_id",
      "species"
    )
  ) %>%
  
  filter(
    calculated_patch_count !=
      n_viable_patches
  )


if (
  nrow(
    patch_count_check
  ) > 0
) {
  
  print(
    patch_count_check
  )
  
  stop(
    "Patch-count validation failed."
  )
}


# ----------------------------------------------------------
# Every included species must have at least one patch
# ----------------------------------------------------------

if (
  any(
    hotspot_species_population_simple$n_viable_patches < 1
  )
) {
  
  stop(
    "At least one hotspot-species record has zero viable patches."
  )
}


# ----------------------------------------------------------
# MVP equivalents should be >= 1 for each assigned patch
# only if each individual patch is at least the MAR.
#
# Because these rasters were already filtered by MAR, this
# should normally be true apart from tiny raster/vector
# numerical differences.
# ----------------------------------------------------------

patch_mar_check <-
  hotspot_patch_level_table %>%
  
  filter(
    patch_area_km2 <
      MAR_km2 * 0.999
  )


if (
  nrow(
    patch_mar_check
  ) > 0
) {
  
  warning(
    nrow(
      patch_mar_check
    ),
    " assigned patch(es) have polygon-derived areas slightly ",
    "below their MAR. Inspect 'patch_mar_check'. This can ",
    "occur because of raster/vector area representation."
  )
}


message(
  "Hotspot-species population table created: ",
  nrow(
    hotspot_species_population_simple
  ),
  " hotspot × species records."
)


# ==========================================================
# 27. Save tables
# ==========================================================


# Complete support results.

write.csv(
  support_long,
  file.path(
    out_dir,
    "BarrierPolygon_SpeciesSupport_LONG.csv"
  ),
  row.names = FALSE
)


# Complete barrier-polygon + PD table.

write.csv(
  hotspot_summary,
  file.path(
    out_dir,
    "BarrierPolygon_FaithPD_Summary_WIDE.csv"
  ),
  row.names = FALSE
)


# Polygons supporting at least one species.

hotspot_summary_nonzero <-
  hotspot_summary[
    hotspot_summary$supported_species != "",
    ,
    drop = FALSE
  ]


write.csv(
  hotspot_summary_nonzero,
  file.path(
    out_dir,
    "BarrierPolygon_FaithPD_Summary_WIDE_nonzero.csv"
  ),
  row.names = FALSE
)


# Multi-species hotspots WITHOUT ranking.

rewilding_hotspots_export <-
  rewilding_hotspots %>%
  
  select(
    -is_multispecies_assemblage
  )


write.csv(
  rewilding_hotspots_export,
  file.path(
    out_dir,
    "Rewilding_Hotspots_FaithPD.csv"
  ),
  row.names = FALSE
)


# Individual assigned viable patches.

write.csv(
  hotspot_patch_level_table,
  file.path(
    out_dir,
    "Hotspot_Species_Viable_Patches_LONG.csv"
  ),
  row.names = FALSE
)


# Full hotspot × species population table.

write.csv(
  hotspot_species_population,
  file.path(
    out_dir,
    "Hotspot_Species_Population_Capacity_FULL.csv"
  ),
  row.names = FALSE
)


# Compact table.

write.csv(
  hotspot_species_population_simple,
  file.path(
    out_dir,
    "Hotspot_Species_Population_Capacity.csv"
  ),
  row.names = FALSE
)


# MAR / MVP lookup used.

write.csv(
  species_population_lookup,
  file.path(
    out_dir,
    "Species_MAR_MVP_Lookup_Used.csv"
  ),
  row.names = FALSE
)


# ----------------------------------------------------------
# Excel workbook containing the main outputs
# ----------------------------------------------------------

writexl::write_xlsx(
  
  list(
    
    Hotspot_Species_Population =
      hotspot_species_population_simple,
    
    Full_Population_Table =
      hotspot_species_population,
    
    Individual_Patches =
      hotspot_patch_level_table,
    
    Hotspot_Summary =
      rewilding_hotspots_export,
    
    MAR_MVP_Lookup =
      species_population_lookup
  ),
  
  path = file.path(
    out_dir,
    "Rewilding_Hotspot_Population_Capacity.xlsx"
  )
)


# ==========================================================
# 28. Prepare hotspot polygons for mapping
# ==========================================================

hotspot_polys <- barrier_polys[
  barrier_polys$barrier_poly_id %in%
    rewilding_hotspots$barrier_poly_id,
]


if (
  is.null(
    hotspot_polys
  ) ||
  nrow(
    hotspot_polys
  ) == 0
) {
  
  stop(
    "No hotspot polygons matched the hotspot IDs."
  )
}


hotspot_order <- match(
  hotspot_polys$barrier_poly_id,
  rewilding_hotspots$barrier_poly_id
)


hotspot_polys$median_faith_PD <-
  rewilding_hotspots$median_faith_PD[
    hotspot_order
  ]


hotspot_polys$species_richness <-
  rewilding_hotspots$species_richness[
    hotspot_order
  ]


hotspot_polys$supported_species <-
  rewilding_hotspots$supported_species[
    hotspot_order
  ]


# ==========================================================
# 29. Assign map colours from Faith's PD
# ==========================================================

pd_min <- min(
  hotspot_polys$median_faith_PD,
  na.rm = TRUE
)


pd_max <- max(
  hotspot_polys$median_faith_PD,
  na.rm = TRUE
)


if (pd_max > pd_min) {
  
  pd_scaled <- (
    hotspot_polys$median_faith_PD -
      pd_min
  ) / (
    pd_max -
      pd_min
  )
  
} else {
  
  pd_scaled <- rep(
    1,
    nrow(
      hotspot_polys
    )
  )
}


pd_palette <- colorRampPalette(
  
  c(
    "#fee8c8",
    "#fdbb84",
    "#e34a33",
    "#7f0000"
  )
  
)(101)


colour_index <- round(
  pd_scaled * 100
) + 1


colour_index <- pmax(
  1,
  pmin(
    101,
    colour_index
  )
)


hotspot_polys$plot_colour <-
  pd_palette[
    colour_index
  ]


# ==========================================================
# 30. Plot UNRANKED Faith's PD hotspot map
# ==========================================================
#
# Numbers are permanent barrier_poly_id values.
# They are NOT ranks.

jpeg(
  
  filename = file.path(
    out_dir,
    "Rewilding_Hotspots_FaithPD_Unranked.jpg"
  ),
  
  width = 3300,
  height = 2500,
  res = 300
)


layout(
  
  matrix(
    c(
      1,
      2
    ),
    nrow = 1
  ),
  
  widths = c(
    8,
    1.3
  )
)


# ----------------------------------------------------------
# Map panel
# ----------------------------------------------------------

par(
  mar = c(
    6,
    1.5,
    3,
    1.5
  )
)


map_extent <- terra::ext(
  boundary_union
)


dx <-
  map_extent$xmax -
  map_extent$xmin

dy <-
  map_extent$ymax -
  map_extent$ymin


xpad <- 0.02 * dx
ypad <- 0.02 * dy


plot(
  
  boundary_union,
  
  col = NA,
  border = "black",
  lwd = 1.5,
  
  xlim = c(
    map_extent$xmin -
      xpad,
    map_extent$xmax +
      xpad
  ),
  
  ylim = c(
    map_extent$ymin -
      ypad,
    map_extent$ymax +
      ypad
  ),
  
  axes = FALSE,
  frame.plot = FALSE
)


plot(
  
  hotspot_polys,
  
  col =
    hotspot_polys$plot_colour,
  
  border = "grey40",
  
  add = TRUE
)


hotspot_centroids <- terra::centroids(
  hotspot_polys
)


title(
  main =
    "Rewilding Hotspots Unranked",
  cex.main = 1.25
)


# ----------------------------------------------------------
# Faith's PD legend
# ----------------------------------------------------------

par(
  mar = c(
    6,
    1,
    3,
    5
  )
)


if (pd_max > pd_min) {
  
  legend_values <- seq(
    pd_min,
    pd_max,
    length.out = length(
      pd_palette
    )
  )
  
} else {
  
  legend_values <- seq(
    pd_min - 0.5,
    pd_max + 0.5,
    length.out = length(
      pd_palette
    )
  )
}


image(
  
  x = 1,
  y = legend_values,
  
  z = matrix(
    legend_values,
    nrow = 1
  ),
  
  col = pd_palette,
  
  zlim = range(
    legend_values
  ),
  
  axes = FALSE,
  xlab = "",
  ylab = "",
  useRaster = TRUE
)


legend_breaks <- pretty(
  
  c(
    pd_min,
    pd_max
  ),
  
  n = 5
)


axis(
  
  side = 4,
  
  at = legend_breaks,
  
  labels = round(
    legend_breaks,
    1
  ),
  
  las = 1,
  cex.axis = 0.8
)


mtext(
  "Median Faith's PD",
  side = 4,
  line = 3,
  cex = 0.85
)


mtext(
  "Lower",
  side = 1,
  line = 2,
  cex = 0.75
)


mtext(
  "Higher",
  side = 3,
  line = 1,
  cex = 0.75
)


dev.off()


# ==========================================================
# 31. Export UNRANKED hotspot polygons
# ==========================================================

# Remove plotting-only colour field before export.

hotspot_polys_export <-
  hotspot_polys


if (
  "plot_colour" %in%
  names(
    hotspot_polys_export
  )
) {
  
  hotspot_polys_export$plot_colour <- NULL
}


terra::writeVector(
  
  hotspot_polys_export,
  
  filename = file.path(
    out_dir,
    "Rewilding_Hotspots_FaithPD.geojson"
  ),
  
  filetype = "GeoJSON",
  overwrite = TRUE
)


# ==========================================================
# 32. Save analysis information
# ==========================================================

analysis_information <- data.frame(
  
  item = c(
    "Patch assignment overlap threshold",
    "Minimum assemblage size",
    "Number of possible phylogenetic assemblages",
    "Number of barrier polygons",
    "Number of multi-species hotspots",
    "Number of hotspot-species records",
    "Number of assigned viable patches",
    "Population calculation",
    "Patch area used for population calculation",
    "Ranking performed",
    "Faith's PD source",
    "MAR source"
  ),
  
  value = c(
    
    overlap_threshold,
    
    minimum_assemblage_size,
    
    nrow(
      faith_pd_lookup
    ),
    
    nrow(
      barrier_polys
    ),
    
    nrow(
      rewilding_hotspots
    ),
    
    nrow(
      hotspot_species_population_simple
    ),
    
    nrow(
      hotspot_patch_level_table
    ),
    
    "Potential population = (total viable patch area / MAR) * MVP",
    
    "Complete area of each patch assigned at >=95% overlap",
    
    "No",
    
    faith_pd_lookup_file,
    
    mar_file
  ),
  
  stringsAsFactors = FALSE
)


write.csv(
  
  analysis_information,
  
  file.path(
    out_dir,
    "Hotspot_Population_Analysis_Information.csv"
  ),
  
  row.names = FALSE
)


# ==========================================================
# 33. Print key output for inspection
# ==========================================================

message(
  "\n=============================================="
)

message(
  "HOTSPOT × SPECIES POPULATION CAPACITY"
)

message(
  "==============================================\n"
)


print(
  hotspot_species_population_simple
)


# ==========================================================
# 34. Completion message
# ==========================================================

message(
  
  "\nHotspot population-capacity analysis complete.\n\n",
  
  "No hotspot ranking has been performed.\n\n",
  
  "Primary population table:\n",
  
  file.path(
    out_dir,
    "Hotspot_Species_Population_Capacity.csv"
  ),
  
  "\n\nExcel workbook:\n",
  
  file.path(
    out_dir,
    "Rewilding_Hotspot_Population_Capacity.xlsx"
  ),
  
  "\n\nOutputs written to:\n",
  
  normalizePath(
    out_dir,
    mustWork = FALSE
  )
)

