############################################################
# Phylogenetic relationships, pairwise distances and
# Faith's Phylogenetic Diversity for ten focal species
############################################################
# Description - R code calculating pairwise phylogenetic distances 
# among ten focal mammals and Faith’s phylogenetic diversity for 
# all 1,013 possible multi-species assemblages across 100 trees. 
# Exports summary tables, comparisons among assemblages of equal 
# species richness, a representative phylogenetic tree, and a distance 
# heatmap.
# ==========================================================
# 1. Packages
# ==========================================================

required_packages <- c(
  "ape"
)

packages_to_install <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}

library(ape)

# ==========================================================
# 2. Input and output directories
# ==========================================================

data_dir <- "./data/"

out_dir <- "./results/Phylogenetic Distances/"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

extract_dir <- file.path(
  data_dir,
  "Tree Pruner Extracted"
)

dir.create(
  extract_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ==========================================================
# 3. Species information
# ==========================================================

species_data <- data.frame(
  
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
  
  # Species-level names used by the VertLife phylogeny
  tree_scientific_name = c(
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
  
  stringsAsFactors = FALSE
)

species_order <- species_data$tree_scientific_name

common_name_lookup <- setNames(
  species_data$common_name,
  species_data$tree_scientific_name
)

focal_name_lookup <- setNames(
  species_data$focal_scientific_name,
  species_data$tree_scientific_name
)

# ==========================================================
# 4. Define all possible assemblages for Faith's PD
# ==========================================================

# Generate every possible assemblage containing at least
# two of the ten focal species.
#
# With ten focal species this gives:
#
# choose(10,2) + choose(10,3) + ... + choose(10,10)
# = 1013 possible multi-species assemblages.
#
# Faith's PD is therefore calculated once for every species
# combination that could subsequently occur within a
# rewilding hotspot.

assemblages <- unlist(
  
  lapply(
    
    2:length(species_order),
    
    function(current_richness) {
      
      combn(
        species_order,
        current_richness,
        simplify = FALSE
      )
    }
  ),
  
  recursive = FALSE
)

names(assemblages) <- sprintf(
  "Assemblage_%04d",
  seq_along(assemblages)
)

message(
  "Number of possible multi-species assemblages: ",
  length(assemblages)
)

stopifnot(
  length(assemblages) == 1013
)

# ----------------------------------------------------------
# Create a deterministic assemblage key
# ----------------------------------------------------------

# Species always occur in species_order. This means the same
# assemblage will always have exactly the same key regardless
# of the order in which species were originally encountered.

make_assemblage_key <- function(selected_species) {
  
  selected_species <- species_order[
    species_order %in%
      unique(selected_species)
  ]
  
  paste(
    selected_species,
    collapse = " | "
  )
}

assemblage_keys <- vapply(
  assemblages,
  make_assemblage_key,
  character(1)
)

if (anyDuplicated(assemblage_keys)) {
  stop("Duplicate assemblage keys were generated.")
}

message(
  "Unique assemblage keys generated: ",
  length(unique(assemblage_keys))
)


# ==========================================================
# 6. Locate the Tree Pruner ZIP file
# ==========================================================

download_stem <-
  "tree-pruner-a7374b0c-0b8a-44f5-bccb-e00714bc3b02"

zip_candidates <- list.files(
  path = data_dir,
  pattern = paste0(
    "^",
    download_stem,
    ".*\\.zip$"
  ),
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(zip_candidates) == 0) {
  
  all_candidates <- list.files(
    path = data_dir,
    pattern = paste0(
      "^",
      download_stem
    ),
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  zip_candidates <- all_candidates[
    vapply(
      all_candidates,
      function(x) {
        !dir.exists(x) &&
          grepl(
            "\\.zip$",
            x,
            ignore.case = TRUE
          )
      },
      logical(1)
    )
  ]
}

if (length(zip_candidates) == 0) {
  
  stop(
    paste0(
      "Could not find the Tree Pruner ZIP file in:\n",
      normalizePath(
        data_dir,
        mustWork = FALSE
      ),
      "\n\nExpected a filename beginning with:\n",
      download_stem
    )
  )
}

if (length(zip_candidates) > 1) {
  
  stop(
    paste0(
      "More than one matching Tree Pruner ZIP file was found:\n",
      paste(
        zip_candidates,
        collapse = "\n"
      )
    )
  )
}

tree_zip_file <- zip_candidates[1]

message(
  "Tree Pruner ZIP file found: ",
  tree_zip_file
)

# ==========================================================
# 7. Extract the Tree Pruner download
# ==========================================================

zip_contents <- unzip(
  tree_zip_file,
  list = TRUE
)

print(zip_contents)

unzip(
  zipfile = tree_zip_file,
  exdir = extract_dir,
  overwrite = TRUE
)

message(
  "Tree files extracted to: ",
  normalizePath(
    extract_dir,
    mustWork = FALSE
  )
)

# ==========================================================
# 8. Locate phylogenetic tree files
# ==========================================================

tree_files <- list.files(
  path = extract_dir,
  pattern = "\\.(tre|tree|trees|nwk|newick|nex|nexus)$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(tree_files) == 0) {
  
  extracted_files <- list.files(
    path = extract_dir,
    full.names = TRUE,
    recursive = TRUE
  )
  
  stop(
    paste0(
      "No recognised phylogenetic tree file was found.\n\n",
      "Extracted files were:\n",
      paste(
        extracted_files,
        collapse = "\n"
      )
    )
  )
}

message(
  "Phylogenetic tree file(s) found:\n",
  paste(
    tree_files,
    collapse = "\n"
  )
)

# ==========================================================
# 9. Read different tree-file formats
# ==========================================================

read_phylogeny_file <- function(file_path) {
  
  extension <- tolower(
    tools::file_ext(file_path)
  )
  
  if (extension %in% c("nex", "nexus")) {
    
    result <- tryCatch(
      ape::read.nexus(file_path),
      error = function(e) {
        stop(
          "Could not read Nexus file:\n",
          file_path,
          "\n\n",
          conditionMessage(e)
        )
      }
    )
    
  } else {
    
    result <- tryCatch(
      ape::read.tree(file_path),
      error = function(e) {
        stop(
          "Could not read Newick tree file:\n",
          file_path,
          "\n\n",
          conditionMessage(e)
        )
      }
    )
  }
  
  result
}

# ==========================================================
# 10. Read all trees
# ==========================================================

tree_objects <- lapply(
  tree_files,
  read_phylogeny_file
)

all_tree_list <- list()

for (object_number in seq_along(tree_objects)) {
  
  current_object <- tree_objects[[object_number]]
  
  if (inherits(current_object, "multiPhylo")) {
    
    all_tree_list <- c(
      all_tree_list,
      unclass(current_object)
    )
    
  } else if (inherits(current_object, "phylo")) {
    
    all_tree_list <- c(
      all_tree_list,
      list(current_object)
    )
    
  } else {
    
    warning(
      "An extracted object was not recognised as phylo or multiPhylo."
    )
  }
}

if (length(all_tree_list) == 0) {
  stop("No readable phylogenetic trees were found.")
}

class(all_tree_list) <- "multiPhylo"

mammal_trees <- all_tree_list

message(
  "Number of phylogenetic trees read: ",
  length(mammal_trees)
)

if (length(mammal_trees) != 100) {
  
  warning(
    "Expected 100 trees, but read ",
    length(mammal_trees),
    ". Inspect 'tree_files' and 'zip_contents'."
  )
}

# ==========================================================
# 11. Clean tip labels
# ==========================================================

clean_tip_labels <- function(x) {
  
  x <- gsub(
    "_",
    " ",
    x
  )
  
  x <- gsub(
    "^['\"]|['\"]$",
    "",
    x
  )
  
  trimws(x)
}

for (tree_number in seq_along(mammal_trees)) {
  
  mammal_trees[[tree_number]]$tip.label <-
    clean_tip_labels(
      mammal_trees[[tree_number]]$tip.label
    )
}

message("Tip labels in the first tree:")

print(
  mammal_trees[[1]]$tip.label
)

# ==========================================================
# 12. Match tree labels to focal species
# ==========================================================

first_tree_labels <- mammal_trees[[1]]$tip.label

missing_species <- setdiff(
  species_order,
  first_tree_labels
)

unexpected_species <- setdiff(
  first_tree_labels,
  species_order
)

if (length(missing_species) > 0) {
  
  stop(
    paste0(
      "The following focal species were not found in the tree:\n",
      paste(
        missing_species,
        collapse = "\n"
      ),
      "\n\nLabels actually present were:\n",
      paste(
        first_tree_labels,
        collapse = "\n"
      )
    )
  )
}

if (length(unexpected_species) > 0) {
  
  message(
    "Additional taxa were present and will be removed:\n",
    paste(
      unexpected_species,
      collapse = "\n"
    )
  )
}

# ==========================================================
# 13. Prune every tree to the ten focal species
# ==========================================================

focal_trees <- lapply(
  mammal_trees,
  function(current_tree) {
    
    ape::keep.tip(
      current_tree,
      species_order
    )
  }
)

class(focal_trees) <- "multiPhylo"

# Rotate tips into a consistent visual order.
# Rotation changes display order but not topology or distances.

for (tree_number in seq_along(focal_trees)) {
  
  focal_trees[[tree_number]] <- ape::rotateConstr(
    focal_trees[[tree_number]],
    species_order
  )
}

# ==========================================================
# 14. Validate branch lengths
# ==========================================================

has_branch_lengths <- vapply(
  focal_trees,
  function(current_tree) {
    
    !is.null(current_tree$edge.length) &&
      all(
        is.finite(
          current_tree$edge.length
        )
      )
  },
  logical(1)
)

if (!all(has_branch_lengths)) {
  
  stop(
    sum(!has_branch_lengths),
    " tree(s) do not contain valid branch lengths."
  )
}

negative_branch_lengths <- vapply(
  focal_trees,
  function(current_tree) {
    
    any(
      current_tree$edge.length < 0,
      na.rm = TRUE
    )
  },
  logical(1)
)

if (any(negative_branch_lengths)) {
  
  stop(
    sum(negative_branch_lengths),
    " tree(s) contain negative branch lengths."
  )
}

ultrametric_checks <- vapply(
  focal_trees,
  ape::is.ultrametric,
  logical(1)
)

message(
  "Ultrametric trees: ",
  sum(ultrametric_checks),
  " of ",
  length(ultrametric_checks)
)

# ==========================================================
# 15. Calculate pairwise distances for every tree
# ==========================================================

n_species <- length(species_order)
n_trees <- length(focal_trees)

distance_array <- array(
  NA_real_,
  dim = c(
    n_species,
    n_species,
    n_trees
  ),
  dimnames = list(
    species_order,
    species_order,
    paste0(
      "Tree_",
      seq_len(n_trees)
    )
  )
)

for (tree_number in seq_along(focal_trees)) {
  
  current_distance_matrix <-
    ape::cophenetic.phylo(
      focal_trees[[tree_number]]
    )
  
  current_distance_matrix <-
    current_distance_matrix[
      species_order,
      species_order
    ]
  
  distance_array[
    ,
    ,
    tree_number
  ] <- current_distance_matrix
}

# ==========================================================
# 16. Summarise pairwise distances across trees
# ==========================================================

median_distance_matrix <- apply(
  distance_array,
  MARGIN = c(1, 2),
  FUN = median,
  na.rm = TRUE
)

lower_95_distance_matrix <- apply(
  distance_array,
  MARGIN = c(1, 2),
  FUN = quantile,
  probs = 0.025,
  na.rm = TRUE
)

upper_95_distance_matrix <- apply(
  distance_array,
  MARGIN = c(1, 2),
  FUN = quantile,
  probs = 0.975,
  na.rm = TRUE
)

message("Median pairwise phylogenetic distances:")

print(
  round(
    median_distance_matrix,
    2
  )
)

# ==========================================================
# 17. Standardise pairwise distances from 0 to 1
# ==========================================================

# These standardised values are retained for visualisation.
#
# 0 = same species
# 1 = most distant pair among the ten focal species

maximum_median_distance <- max(
  median_distance_matrix,
  na.rm = TRUE
)

standardised_distance_matrix <-
  median_distance_matrix /
  maximum_median_distance

print(
  round(
    standardised_distance_matrix,
    3
  )
)

# ==========================================================
# 18. Create one row for each unique species pair
# ==========================================================

pair_indices <- which(
  upper.tri(
    median_distance_matrix
  ),
  arr.ind = TRUE
)

pairwise_distances <- data.frame(
  
  species_1 =
    rownames(median_distance_matrix)[
      pair_indices[, "row"]
    ],
  
  species_2 =
    colnames(median_distance_matrix)[
      pair_indices[, "col"]
    ],
  
  median_phylogenetic_distance =
    median_distance_matrix[pair_indices],
  
  lower_95_phylogenetic_distance =
    lower_95_distance_matrix[pair_indices],
  
  upper_95_phylogenetic_distance =
    upper_95_distance_matrix[pair_indices],
  
  standardised_distance =
    standardised_distance_matrix[pair_indices],
  
  stringsAsFactors = FALSE
)

pairwise_distances$common_name_1 <-
  common_name_lookup[
    pairwise_distances$species_1
  ]

pairwise_distances$common_name_2 <-
  common_name_lookup[
    pairwise_distances$species_2
  ]

pairwise_distances$focal_scientific_name_1 <-
  focal_name_lookup[
    pairwise_distances$species_1
  ]

pairwise_distances$focal_scientific_name_2 <-
  focal_name_lookup[
    pairwise_distances$species_2
  ]

pairwise_distances <- pairwise_distances[
  order(
    pairwise_distances$median_phylogenetic_distance
  ),
]

row.names(pairwise_distances) <- NULL

identifier_columns <- c(
  "common_name_1",
  "focal_scientific_name_1",
  "species_1",
  "common_name_2",
  "focal_scientific_name_2",
  "species_2"
)

pairwise_distances <- pairwise_distances[
  ,
  c(
    identifier_columns,
    setdiff(
      names(pairwise_distances),
      identifier_columns
    )
  )
]

print(pairwise_distances)

# ==========================================================
# 19. Faith's Phylogenetic Diversity function
# ==========================================================

# Faith's PD is the sum of the branch lengths in the subtree
# connecting the selected species.
#
# All assemblages are evaluated against the same focal-tree
# root, making them directly comparable within this analysis.

calculate_faith_pd <- function(tree, selected_species) {
  
  selected_species <- unique(
    selected_species
  )
  
  missing_tips <- setdiff(
    selected_species,
    tree$tip.label
  )
  
  if (length(missing_tips) > 0) {
    
    stop(
      "The following species are absent from the tree: ",
      paste(
        missing_tips,
        collapse = ", "
      )
    )
  }
  
  assemblage_tree <- ape::keep.tip(
    tree,
    selected_species
  )
  
  sum(
    assemblage_tree$edge.length,
    na.rm = TRUE
  )
}

# ==========================================================
# 20. Calculate Faith's PD across all 100 trees
# ==========================================================

assemblage_names <- names(assemblages)

faith_pd_array <- matrix(
  NA_real_,
  nrow = length(assemblages),
  ncol = n_trees,
  dimnames = list(
    assemblage_names,
    paste0(
      "Tree_",
      seq_len(n_trees)
    )
  )
)

for (assemblage_name in assemblage_names) {
  
  for (tree_number in seq_along(focal_trees)) {
    
    faith_pd_array[
      assemblage_name,
      tree_number
    ] <- calculate_faith_pd(
      tree = focal_trees[[tree_number]],
      selected_species =
        assemblages[[assemblage_name]]
    )
  }
}

# ==========================================================
# 21. Standardise PD within species richness
# ==========================================================

# faith_pd_array now already contains every possible
# assemblage of 2-10 species.
#
# Therefore, the distribution of Faith's PD for a given
# species richness can be obtained directly from this matrix
# without recalculating Faith's PD.
#
# Outputs:
#
# relative_PD:
# observed PD / maximum possible PD for that richness
#
# percentile:
# percentage of possible assemblages of the same richness
# whose PD is less than or equal to the observed assemblage
#
# SES_PD:
# standardised effect size relative to all possible
# assemblages of the same richness

assemblage_richness <- vapply(
  assemblages,
  length,
  integer(1)
)

relative_pd_array <- matrix(
  NA_real_,
  nrow = length(assemblages),
  ncol = n_trees,
  dimnames = dimnames(faith_pd_array)
)

percentile_pd_array <- matrix(
  NA_real_,
  nrow = length(assemblages),
  ncol = n_trees,
  dimnames = dimnames(faith_pd_array)
)

ses_pd_array <- matrix(
  NA_real_,
  nrow = length(assemblages),
  ncol = n_trees,
  dimnames = dimnames(faith_pd_array)
)

for (
  current_richness in
  sort(
    unique(
      assemblage_richness
    )
  )
) {
  
  message(
    "Standardising Faith's PD for richness = ",
    current_richness
  )
  
  relevant_rows <- which(
    assemblage_richness ==
      current_richness
  )
  
  richness_pd <- faith_pd_array[
    relevant_rows,
    ,
    drop = FALSE
  ]
  
  for (
    tree_number in
    seq_len(n_trees)
  ) {
    
    current_values <-
      richness_pd[
        ,
        tree_number
      ]
    
    current_mean <- mean(
      current_values,
      na.rm = TRUE
    )
    
    current_sd <- stats::sd(
      current_values,
      na.rm = TRUE
    )
    
    current_max <- max(
      current_values,
      na.rm = TRUE
    )
    
    # --------------------------------------------------------
    # Relative PD
    # --------------------------------------------------------
    
    relative_pd_array[
      relevant_rows,
      tree_number
    ] <-
      current_values /
      current_max
    
    # --------------------------------------------------------
    # PD percentile
    # --------------------------------------------------------
    
    percentile_pd_array[
      relevant_rows,
      tree_number
    ] <-
      vapply(
        
        current_values,
        
        function(current_pd) {
          
          mean(
            current_values <=
              current_pd,
            na.rm = TRUE
          ) * 100
        },
        
        numeric(1)
      )
    
    # --------------------------------------------------------
    # Standardised effect size
    # --------------------------------------------------------
    
    if (
      is.finite(current_sd) &&
      current_sd > 0
    ) {
      
      ses_pd_array[
        relevant_rows,
        tree_number
      ] <-
        (
          current_values -
            current_mean
        ) /
        current_sd
    }
  }
}

# ==========================================================
# 22. Summarise Faith's PD results
# ==========================================================

assemblage_summary <- data.frame(
  
  assemblage =
    names(assemblages),
  
  assemblage_key =
    assemblage_keys,
  
  species_richness =
    assemblage_richness,
  
  supported_species_common = vapply(
    
    assemblages,
    
    function(current_species) {
      
      paste(
        common_name_lookup[
          current_species
        ],
        collapse = " | "
      )
    },
    
    character(1)
  ),
  
  median_faith_PD = apply(
    faith_pd_array,
    1,
    median,
    na.rm = TRUE
  ),
  
  lower_95_faith_PD = apply(
    faith_pd_array,
    1,
    quantile,
    probs = 0.025,
    na.rm = TRUE
  ),
  
  upper_95_faith_PD = apply(
    faith_pd_array,
    1,
    quantile,
    probs = 0.975,
    na.rm = TRUE
  ),
  
  median_relative_PD = apply(
    relative_pd_array,
    1,
    median,
    na.rm = TRUE
  ),
  
  lower_95_relative_PD = apply(
    relative_pd_array,
    1,
    quantile,
    probs = 0.025,
    na.rm = TRUE
  ),
  
  upper_95_relative_PD = apply(
    relative_pd_array,
    1,
    quantile,
    probs = 0.975,
    na.rm = TRUE
  ),
  
  median_PD_percentile = apply(
    percentile_pd_array,
    1,
    median,
    na.rm = TRUE
  ),
  
  median_SES_PD = apply(
    ses_pd_array,
    1,
    median,
    na.rm = TRUE
  ),
  
  stringsAsFactors = FALSE
)

assemblage_summary <- assemblage_summary[
  
  order(
    assemblage_summary$species_richness,
    -assemblage_summary$median_faith_PD
  ),
  
]

row.names(
  assemblage_summary
) <- NULL

print(
  assemblage_summary
)

message(
  "Assemblages summarised: ",
  nrow(
    assemblage_summary
  )
)

stopifnot(
  nrow(
    assemblage_summary
  ) == 1013
)

# ==========================================================
# 23. Create assemblage membership table
# ==========================================================

assemblage_membership <- do.call(
  rbind,
  lapply(
    names(assemblages),
    function(assemblage_name) {
      
      current_species <-
        assemblages[[assemblage_name]]
      
      data.frame(
        assemblage = assemblage_name,
        common_name =
          common_name_lookup[current_species],
        focal_scientific_name =
          focal_name_lookup[current_species],
        tree_scientific_name =
          current_species,
        stringsAsFactors = FALSE
      )
    }
  )
)

row.names(assemblage_membership) <- NULL

# ==========================================================
# 24. Select a representative tree
# ==========================================================

# The representative tree is the sampled tree whose complete
# distance matrix is closest to the median distance matrix.

tree_distance_from_median <- vapply(
  seq_along(focal_trees),
  function(tree_number) {
    
    current_matrix <- distance_array[
      ,
      ,
      tree_number
    ]
    
    sum(
      (
        current_matrix -
          median_distance_matrix
      )^2,
      na.rm = TRUE
    )
  },
  numeric(1)
)

representative_tree_number <- which.min(
  tree_distance_from_median
)

representative_tree <-
  focal_trees[[representative_tree_number]]

representative_tree <- ape::ladderize(
  representative_tree,
  right = FALSE
)

message(
  "Representative tree selected: Tree ",
  representative_tree_number
)

# Clean common-name labels only.
# Individual phylogenetic-isolation scores have been removed.

plot_tip_labels <-
  common_name_lookup[
    representative_tree$tip.label
  ]

# ==========================================================
# 25. Save tables
# ==========================================================

# ----------------------------------------------------------
# Master Faith's PD lookup for all possible assemblages
# ----------------------------------------------------------

write.csv(
  
  assemblage_summary,
  
  file.path(
    out_dir,
    "Faith_PD_All_Possible_Assemblages_Summary.csv"
  ),
  
  row.names = FALSE
)

write.csv(
  species_data,
  file.path(
    out_dir,
    "Focal_Species_Names.csv"
  ),
  row.names = FALSE
)

write.csv(
  round(
    median_distance_matrix,
    6
  ),
  file.path(
    out_dir,
    "Phylogenetic_Distance_Matrix_Median.csv"
  ),
  row.names = TRUE
)

write.csv(
  round(
    lower_95_distance_matrix,
    6
  ),
  file.path(
    out_dir,
    "Phylogenetic_Distance_Matrix_Lower95.csv"
  ),
  row.names = TRUE
)

write.csv(
  round(
    upper_95_distance_matrix,
    6
  ),
  file.path(
    out_dir,
    "Phylogenetic_Distance_Matrix_Upper95.csv"
  ),
  row.names = TRUE
)

write.csv(
  round(
    standardised_distance_matrix,
    6
  ),
  file.path(
    out_dir,
    "Phylogenetic_Distance_Matrix_Standardised_0_to_1.csv"
  ),
  row.names = TRUE
)

write.csv(
  pairwise_distances,
  file.path(
    out_dir,
    "Phylogenetic_Distances_Pairwise.csv"
  ),
  row.names = FALSE
)

write.csv(
  assemblage_summary,
  file.path(
    out_dir,
    "Faith_PD_Assemblage_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  assemblage_membership,
  file.path(
    out_dir,
    "Faith_PD_Assemblage_Membership.csv"
  ),
  row.names = FALSE
)

write.csv(
  faith_pd_array,
  file.path(
    out_dir,
    "Faith_PD_All_100_Trees.csv"
  ),
  row.names = TRUE
)

write.csv(
  relative_pd_array,
  file.path(
    out_dir,
    "Faith_PD_Relative_All_100_Trees.csv"
  ),
  row.names = TRUE
)

write.csv(
  percentile_pd_array,
  file.path(
    out_dir,
    "Faith_PD_Percentiles_All_100_Trees.csv"
  ),
  row.names = TRUE
)

ape::write.tree(
  representative_tree,
  file = file.path(
    out_dir,
    "Representative_Dated_Phylogenetic_Tree.newick"
  )
)

# ==========================================================
# 26. Plot the representative dated tree
# ==========================================================

jpeg(
  filename = file.path(
    out_dir,
    "Focal_Species_Dated_Phylogenetic_Tree.jpg"
  ),
  width = 3600,
  height = 2600,
  res = 300
)

par(
  mar = c(5, 2, 5, 15),
  xpd = TRUE
)

plot(
  representative_tree,
  type = "phylogram",
  direction = "rightwards",
  show.tip.label = FALSE,
  edge.width = 2,
  no.margin = FALSE
)

ape::tiplabels(
  text = plot_tip_labels,
  tip = seq_along(
    representative_tree$tip.label
  ),
  frame = "none",
  adj = c(-0.02, 0.5),
  cex = 0.9
)

ape::axisPhylo(
  backward = TRUE
)

title(
  main = paste(
    "Dated phylogenetic relationships",
    "among focal species"
  ),
  line = 2,
  cex.main = 1.4
)

mtext(
  paste0(
    "Representative tree ",
    representative_tree_number,
    " of ",
    n_trees
  ),
  side = 1,
  line = 3,
  cex = 0.8
)

dev.off()

# ==========================================================
# 27. Plot phylogenetic-distance heatmap with legend
# ==========================================================

scientific_labels <- species_data$focal_scientific_name[
  match(
    rownames(standardised_distance_matrix),
    species_data$tree_scientific_name
  )
]

heatmap_matrix <- standardised_distance_matrix

rownames(heatmap_matrix) <- scientific_labels
colnames(heatmap_matrix) <- scientific_labels

# Colour scale
heatmap_colours <- colorRampPalette(
  c(
    "navy",
    "white",
    "firebrick3"
  )
)(101)

jpeg(
  filename = file.path(
    out_dir,
    "Focal_Species_Phylogenetic_Distance_Heatmap.jpg"
  ),
  width = 3600,
  height = 3000,
  res = 300
)

# Heatmap on the left; colour legend on the right
layout(
  matrix(
    c(1, 2),
    nrow = 1
  ),
  widths = c(
    8,
    1.2
  )
)

# ----------------------------------------------------------
# Heatmap panel
# ----------------------------------------------------------

par(
  mar = c(14, 14, 5, 2)
)

n_heatmap_species <- nrow(
  heatmap_matrix
)

# image() expects the matrix orientation to be transposed.
# Reversing the rows keeps the first species at the top.

heatmap_for_plot <- t(
  heatmap_matrix[
    n_heatmap_species:1,
    ,
    drop = FALSE
  ]
)

image(
  x = seq_len(
    ncol(heatmap_matrix)
  ),
  y = seq_len(
    nrow(heatmap_matrix)
  ),
  z = heatmap_for_plot,
  col = heatmap_colours,
  zlim = c(0, 1),
  axes = FALSE,
  xlab = "",
  ylab = "",
  useRaster = TRUE
)

axis(
  side = 1,
  at = seq_len(
    ncol(heatmap_matrix)
  ),
  labels = colnames(
    heatmap_matrix
  ),
  las = 2,
  tick = FALSE,
  cex.axis = 0.75
)

axis(
  side = 2,
  at = seq_len(
    nrow(heatmap_matrix)
  ),
  labels = rev(
    rownames(
      heatmap_matrix
    )
  ),
  las = 2,
  tick = FALSE,
  cex.axis = 0.75
)

# Add cell borders
abline(
  v = seq(
    1.5,
    n_heatmap_species - 0.5,
    by = 1
  ),
  h = seq(
    1.5,
    n_heatmap_species - 0.5,
    by = 1
  ),
  col = "grey80",
  lwd = 0.5
)

box()

title(
  main = paste(
    "Median standardised phylogenetic distance",
    "across 100 trees"
  ),
  line = 2,
  cex.main = 1.25
)

# Add distance values to cells
for (column_number in seq_len(ncol(heatmap_matrix))) {
  
  for (row_number in seq_len(nrow(heatmap_matrix))) {
    
    original_row_number <-
      nrow(heatmap_matrix) -
      row_number +
      1
    
    current_value <-
      heatmap_matrix[
        original_row_number,
        column_number
      ]
    
    text(
      x = column_number,
      y = row_number,
      labels = sprintf(
        "%.2f",
        current_value
      ),
      cex = 0.55
    )
  }
}

# ----------------------------------------------------------
# Legend panel
# ----------------------------------------------------------

par(
  mar = c(14, 1, 5, 5)
)

legend_values <- seq(
  0,
  1,
  length.out = length(
    heatmap_colours
  )
)

image(
  x = 1,
  y = legend_values,
  z = matrix(
    legend_values,
    nrow = 1
  ),
  col = heatmap_colours,
  zlim = c(0, 1),
  axes = FALSE,
  xlab = "",
  ylab = "",
  useRaster = TRUE
)

axis(
  side = 4,
  at = seq(
    0,
    1,
    by = 0.2
  ),
  labels = sprintf(
    "%.1f",
    seq(
      0,
      1,
      by = 0.2
    )
  ),
  las = 1,
  cex.axis = 0.8
)

mtext(
  "Standardised phylogenetic distance",
  side = 4,
  line = 3,
  cex = 0.85
)

mtext(
  "Closer",
  side = 1,
  line = 2,
  cex = 0.7
)

mtext(
  "More distant",
  side = 3,
  line = 1,
  cex = 0.7
)

dev.off()


# ==========================================================
# 28. Save run information
# ==========================================================

run_information <- data.frame(
  
  item = c(
    "ZIP file",
    "Number of tree files",
    "Number of trees",
    "Number of focal species",
    "Number of assemblages",
    "Number of ultrametric trees",
    "Representative tree number",
    "Maximum median pairwise distance"
  ),
  
  value = c(
    tree_zip_file,
    length(tree_files),
    n_trees,
    n_species,
    length(assemblages),
    sum(ultrametric_checks),
    representative_tree_number,
    maximum_median_distance
  ),
  
  stringsAsFactors = FALSE
)

write.csv(
  run_information,
  file.path(
    out_dir,
    "Phylogenetic_Analysis_Run_Information.csv"
  ),
  row.names = FALSE
)

# ==========================================================
# 29. Completion message
# ==========================================================

message(
  "\nPhylogenetic analysis complete.\n",
  "Outputs saved to:\n",
  normalizePath(
    out_dir,
    mustWork = FALSE
  )
)

