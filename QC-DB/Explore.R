library(data.table)
library(dplyr)
library(ggplot2)

results_dir <- "/rna/liha/phylogenomics_practice/QC-DB/Results"
plots_dir <- "/rna/liha/phylogenomics_practice/QC-DB/Plots"
summary_file <- "/rna/liha/phylogenomics_practice/QC-DB/assembly_summary_refseq.txt"

if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

dirs <- list.dirs(results_dir, recursive = FALSE)
dirs <- dirs[!grepl("/\\.", dirs)]

# 1. Read and combine all data from all directories
all_data <- data.frame()

for (d in dirs) {
  dir_name <- basename(d)
  clade <- strsplit(dir_name, "-")[[1]][1]
  
  tsv_path <- file.path(d, "03_stats", "basic_stats.tsv")
  
  if (file.exists(tsv_path)) {
    stats <- tryCatch(
      read.delim(tsv_path, sep = "\t", stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(stats)) {
      stats$clade <- clade

      # BUSCO: pull the C/S/D/F/M percentages out of the short_summary text file
      stats$busco_C <- NA
      stats$busco_S <- NA
      stats$busco_D <- NA
      stats$busco_F <- NA
      stats$busco_M <- NA
      stats$busco_n <- NA

      busco_file <- list.files(file.path(d, "04_busco"),
                                pattern = "^short_summary\\.specific\\..*\\.txt$",
                                recursive = TRUE, full.names = TRUE)
      if (length(busco_file) > 0) {
        busco_lines <- tryCatch(readLines(busco_file[1], warn = FALSE), error = function(e) NULL)
        one_line <- grep("^\\s*C:", busco_lines, value = TRUE)
        if (length(one_line) > 0) {
          m <- regmatches(one_line[1], regexec(
            "C:([0-9.]+)%\\[S:([0-9.]+)%,D:([0-9.]+)%\\],F:([0-9.]+)%,M:([0-9.]+)%,n:([0-9]+)",
            one_line[1]))[[1]]
          if (length(m) == 7) {
            stats$busco_C <- as.numeric(m[2])
            stats$busco_S <- as.numeric(m[3])
            stats$busco_D <- as.numeric(m[4])
            stats$busco_F <- as.numeric(m[5])
            stats$busco_M <- as.numeric(m[6])
            stats$busco_n <- as.numeric(m[7])
          }
        }
      }

      all_data <- rbind(all_data, stats)
    }
  }
}
cat(sprintf("Successfully loaded basic data for %d entries.\n", nrow(all_data)))

# 2. Fetch Genome Quality from assembly_summary_refseq.txt and NCBI FTP
# Initialize the columns beforehand so they exist even if fetching fails
all_data$scaffold_N50 <- NA
all_data$contig_N50 <- NA
all_data$genome_size <- NA
all_data$normalized_scaffold_N50 <- NA

if (file.exists(summary_file) && nrow(all_data) > 0) {
  cat("Fetching genome quality stats from NCBI FTP (this might take a moment)...\n")
  
  # Quickly grep only the accessions we need from the large summary file
  accs <- unique(all_data$accession)
  acc_pattern <- paste(accs, collapse = "|") # Fixed regex pipe
  cmd <- sprintf("grep -E '%s' %s", acc_pattern, summary_file)
  
  sub_summary <- tryCatch(
    read.delim(pipe(cmd), header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = ""),
    error = function(e) NULL
  )
  
  if (!is.null(sub_summary) && nrow(sub_summary) > 0) {
    # V1 is accession, V20 is ftp path, V26 is genome_size
    
    for (i in 1:nrow(all_data)) {
      acc <- all_data$accession[i]
      idx <- match(acc, sub_summary$V1)
      if (!is.na(idx)) {
        # Get genome size
        genome_sz <- as.numeric(sub_summary$V26[idx])
        if (!is.na(genome_sz)) {
          all_data$genome_size[i] <- genome_sz
        }
        
        ftp <- sub_summary$V20[idx]
        if (!is.na(ftp) && ftp != "na" && ftp != "") {
          base_name <- basename(ftp)
          stats_url <- paste0(ftp, "/", base_name, "_assembly_stats.txt")
          
          # Read stats file from FTP
          lines <- tryCatch(readLines(stats_url, warn = FALSE), error = function(e) NULL)
          
          if (!is.null(lines)) {
            # Extract scaffold-N50
            scaf_line <- grep("all\tall\tall\tall\tscaffold-N50", lines, value = TRUE)
            if (length(scaf_line) > 0) {
              all_data$scaffold_N50[i] <- as.numeric(strsplit(scaf_line[1], "\t")[[1]][6])
            }
            # Extract contig-N50
            cont_line <- grep("all\tall\tall\tall\tcontig-N50", lines, value = TRUE)
            if (length(cont_line) > 0) {
              all_data$contig_N50[i] <- as.numeric(strsplit(cont_line[1], "\t")[[1]][6])
            }
          }
        }
      }
    }
    
    # Calculate normalized metric: scaffold_N50 / genome_size
    all_data$normalized_scaffold_N50 <- all_data$scaffold_N50 / all_data$genome_size
    
    cat("Successfully fetched genome quality stats and computed normalized metrics.\n")
  } else {
    cat("Could not find accessions in assembly_summary_refseq.txt.\n")
  }
}

# 3. Define a function to easily plot any column interactively
plot_column <- function(data, col_name, log_y = TRUE) {
  if (!col_name %in% colnames(data)) {
    stop(paste("Column", col_name, "not found in data."))
  }
  
  p <- ggplot(data, aes_string(x = "clade", y = col_name, color = "clade")) +
    geom_jitter(width = 0.2, alpha = 0.7) +
    theme_bw() + # White background theme
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) +
    labs(title = sprintf("Distribution of %s per Clade", col_name),
         x = "Clade",
         y = col_name)
         
  if (log_y) {
    p <- p + scale_y_log10()
  }
  
  output_file <- file.path(plots_dir, sprintf("%s_distribution.png", col_name))
  ggsave(output_file, plot = p, width = 10, height = 6, bg = "white")
  cat(sprintf("Plot saved to %s\n", output_file))
  
  return(p) # Returns the plot object so you can view it interactively
}

# 4. Example Usage (Uncomment or run interactively):
p_n50 <- plot_column(all_data, "protein_length_N50")
p_scaffold_n50 <- plot_column(all_data, "scaffold_N50")







species_info <- fread("/rna/liha/phylogenomics_practice/QC-DB/Scripts/species_info.tsv")
all_data <- left_join(all_data, species_info %>% select(Taxonomy, requested_species, selected_taxon, taxid), by = "taxid")

colnames(all_data)

png("/rna/liha/phylogenomics_practice/QC-DB/Plots/contig_N50_vs_BUSCO_C.png", width = 800, height = 600)
plot(log10(all_data$contig_N50), all_data$busco_C, xlab = "Contig N50 (log10)", ylab = "BUSCO Complete (%)", main = "Contig N50 vs BUSCO Complete")
dev.off()

png("/rna/liha/phylogenomics_practice/QC-DB/Plots/contig_N50_vs_BUSCO_D.png", width = 800, height = 600)
plot(log10(all_data$contig_N50), all_data$busco_D, xlab = "Contig N50 (log10)", ylab = "BUSCO Duplicated (%)", main = "Contig N50 vs BUSCO Duplicated")
dev.off()
 


all_data %>% select(Taxonomy, requested_species, n_proteins, under_50aa, busco_C, busco_D, busco_F)


species <- fread("/rna/liha/phylogenomics_practice/QC-DB/Tmp/assembly_protein_count.tsv")
species %>% filter(in_refseq_protein == "YES") %>% 
filter(!is.na(protein_coding_genes)) %>%
select(Taxonomy, requested_species, protein_coding_genes, assembly_level, contig_n50, genome_size) %>% fwrite("/rna/liha/phylogenomics_practice/QC-DB/Tmp/species_summary.tsv")

