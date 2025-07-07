# ==============================================================================
# Data Setup Script for NMD Gene Analysis
# This script helps set up the data files needed for the analysis pipeline
# ==============================================================================

cat("NMD Gene Analysis - Data Setup\n")
cat("===============================\n\n")

# ==============================================================================
# Configuration
# ==============================================================================

# Define data files and their locations
data_files <- list(
  # Required core data files
  required = list(
    "significant_SNP.RData" = "../significant_SNP.RData",
    "NMD_gene_list.RData" = "../NMD_gene_list.RData",
    "plot_man_anno.csv" = "../plot_man_anno.csv",
    "sig_snp_gene_id_mapped_res.csv" = "../sig_snp_gene_id_mapped_res.csv",
    "NMD_Related_Genes_and_Proteins.csv" = "../NMD_Related_Genes_and_Proteins.csv"
  ),
  
  # Optional data files for full functionality
  optional = list(
    "gwas_catalog_v1.0.2-associations_e112_r2024-09-22.tsv" = "../gwas_catalog_v1.0.2-associations_e112_r2024-09-22.tsv",
    "sig_SNP_anno_res_biomart.csv" = "../sig_SNP_anno_res_biomart.csv",
    "sig_SNP_anno_res.csv" = "../sig_SNP_anno_res.csv"
  ),
  
  # Directory structures
  directories = list(
    "pathway/" = "../pathway/",
    "RBP_map/" = "../RBP_map/",
    "Cancer_KO/" = "../Cancer_KO/", 
    "GWAS_overlap/" = "../GWAS_overlap/",
    "protein_map/" = "../protein_map/"
  )
)

# ==============================================================================
# Helper Functions
# ==============================================================================

# Function to check if file exists
file_exists_with_message <- function(file_path, file_name) {
  if(file.exists(file_path)) {
    cat("✓ Found:", file_name, "\n")
    return(TRUE)
  } else {
    cat("✗ Missing:", file_name, "\n")
    return(FALSE)
  }
}

# Function to copy or link files
setup_file <- function(source_path, dest_path, method = "copy") {
  if(!file.exists(source_path)) {
    return(FALSE)
  }
  
  # Create destination directory if needed
  dest_dir <- dirname(dest_path)
  if(!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
  }
  
  tryCatch({
    if(method == "copy") {
      file.copy(source_path, dest_path, overwrite = TRUE)
    } else if(method == "link") {
      file.symlink(source_path, dest_path)
    }
    return(TRUE)
  }, error = function(e) {
    cat("Error setting up", dest_path, ":", e$message, "\n")
    return(FALSE)
  })
}

# Function to setup directory structure
setup_directory <- function(source_dir, dest_dir, method = "copy") {
  if(!dir.exists(source_dir)) {
    return(FALSE)
  }
  
  if(!dir.exists(dest_dir)) {
    dir.create(dest_dir, recursive = TRUE)
  }
  
  tryCatch({
    if(method == "copy") {
      # Copy entire directory
      file.copy(source_dir, dirname(dest_dir), overwrite = TRUE, recursive = TRUE)
    } else if(method == "link") {
      # Create symbolic link
      file.symlink(source_dir, dest_dir)
    }
    return(TRUE)
  }, error = function(e) {
    cat("Error setting up directory", dest_dir, ":", e$message, "\n")
    return(FALSE)
  })
}

# ==============================================================================
# Data Availability Check
# ==============================================================================

check_data_availability <- function() {
  cat("Checking data availability...\n")
  cat(rep("-", 40), "\n")
  
  # Check required files
  cat("\nRequired Files:\n")
  required_available <- 0
  for(name in names(data_files$required)) {
    path <- data_files$required[[name]]
    if(file_exists_with_message(path, name)) {
      required_available <- required_available + 1
    }
  }
  
  # Check optional files
  cat("\nOptional Files:\n")
  optional_available <- 0
  for(name in names(data_files$optional)) {
    path <- data_files$optional[[name]]
    if(file_exists_with_message(path, name)) {
      optional_available <- optional_available + 1
    }
  }
  
  # Check directories
  cat("\nData Directories:\n")
  dir_available <- 0
  for(name in names(data_files$directories)) {
    path <- data_files$directories[[name]]
    if(file_exists_with_message(path, name)) {
      dir_available <- dir_available + 1
    }
  }
  
  # Summary
  cat(rep("-", 40), "\n")
  cat("Summary:\n")
  cat(sprintf("Required files: %d/%d available\n", required_available, length(data_files$required)))
  cat(sprintf("Optional files: %d/%d available\n", optional_available, length(data_files$optional)))
  cat(sprintf("Directories: %d/%d available\n", dir_available, length(data_files$directories)))
  
  return(list(
    required = required_available == length(data_files$required),
    optional = optional_available,
    directories = dir_available
  ))
}

# ==============================================================================
# Data Setup Functions
# ==============================================================================

setup_data <- function(method = "copy", setup_optional = TRUE) {
  cat(sprintf("\nSetting up data using method: %s\n", method))
  cat(rep("-", 40), "\n")
  
  success_count <- 0
  total_count <- 0
  
  # Setup required files
  cat("\nSetting up required files:\n")
  for(name in names(data_files$required)) {
    source_path <- data_files$required[[name]]
    dest_path <- file.path("data", name)
    
    total_count <- total_count + 1
    if(setup_file(source_path, dest_path, method)) {
      cat("✓ Setup:", name, "\n")
      success_count <- success_count + 1
    } else {
      cat("✗ Failed:", name, "\n")
    }
  }
  
  # Setup optional files
  if(setup_optional) {
    cat("\nSetting up optional files:\n")
    for(name in names(data_files$optional)) {
      source_path <- data_files$optional[[name]]
      dest_path <- file.path("data", name)
      
      total_count <- total_count + 1
      if(setup_file(source_path, dest_path, method)) {
        cat("✓ Setup:", name, "\n")
        success_count <- success_count + 1
      } else {
        cat("⚠ Skipped:", name, "\n")
      }
    }
    
    # Setup directories
    cat("\nSetting up data directories:\n")
    for(name in names(data_files$directories)) {
      source_path <- data_files$directories[[name]]
      dest_path <- file.path("data", name)
      
      if(setup_directory(source_path, dest_path, method)) {
        cat("✓ Setup directory:", name, "\n")
      } else {
        cat("⚠ Skipped directory:", name, "\n")
      }
    }
  }
  
  cat(rep("-", 40), "\n")
  cat(sprintf("Setup completed: %d/%d files processed successfully\n", success_count, total_count))
  
  return(success_count >= length(data_files$required))
}

# ==============================================================================
# Interactive Setup
# ==============================================================================

interactive_setup <- function() {
  cat("\nInteractive Data Setup\n")
  cat("Choose an option:\n")
  cat("1. Check data availability only\n")
  cat("2. Copy files to Git/data/ directory\n")
  cat("3. Create symbolic links to original files\n")
  cat("4. Setup required files only\n")
  cat("5. Exit\n")
  
  choice <- readline(prompt = "Enter your choice (1-5): ")
  
  switch(choice,
    "1" = {
      availability <- check_data_availability()
      if(availability$required) {
        cat("\n✅ All required files are available. You can run the analyses.\n")
      } else {
        cat("\n⚠ Some required files are missing. Please ensure they are available.\n")
      }
    },
    "2" = {
      cat("\nCopying files...\n")
      setup_data(method = "copy", setup_optional = TRUE)
    },
    "3" = {
      cat("\nCreating symbolic links...\n")
      setup_data(method = "link", setup_optional = TRUE)
    },
    "4" = {
      cat("\nSetting up required files only...\n")
      setup_data(method = "copy", setup_optional = FALSE)
    },
    "5" = {
      cat("Exiting setup.\n")
      return()
    },
    {
      cat("Invalid choice. Please run the script again.\n")
    }
  )
}

# ==============================================================================
# Main Execution
# ==============================================================================

# Create data directory if it doesn't exist
if(!dir.exists("data")) {
  dir.create("data")
  cat("Created data directory.\n")
}

# Check if running interactively
if(interactive()) {
  interactive_setup()
} else {
  # Non-interactive mode - just check availability
  cat("Running in non-interactive mode...\n")
  availability <- check_data_availability()
  
  if(!availability$required) {
    cat("\n⚠ WARNING: Some required files are missing!\n")
    cat("Please run this script interactively or ensure all required files are available.\n")
    cat("Required files should be in the parent directory of this Git folder.\n")
  } else {
    cat("\n✅ All required files are available!\n")
    cat("You can now run the analysis pipeline.\n")
  }
}

cat("\nFor more information, see the README.md file.\n") 