suppressWarnings(suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Biostrings)
})) 

# Define project paths.
base_dir <- "/rna/liha/phylogenomics_practice/SLIT2"
scripts_dir <- file.path(base_dir, "Scripts")

species_file <- file.path(scripts_dir, "species.txt")
fasta_path <- file.path(scripts_dir, "SLIT-homologs.fasta")
output_png <- file.path(scripts_dir, "SLIT_Homolog_Counts.png")








# Load the reference species table used to define plotting order.
species_list <- fread(species_file, header = TRUE)

# Stop early if the filtered homolog FASTA has not been generated yet.
if (!file.exists(fasta_path) || file.info(fasta_path)$size == 0) {
  stop("Missing or empty FASTA file: ", fasta_path)
}

# Read filtered SLIT homolog protein sequences.
seqs <- readAAStringSet(fasta_path, format = "fasta")

# Extract the clade and species label from each FASTA header.
# Example:
# Cephalochordata_Branchiostoma_belcheri_H1_xxx
# becomes:
# Cephalochordata - Branchiostoma belcheri
organism_raw <- sub("_H[0-9]+_.*", "", names(seqs))
extracted_species <- gsub("_", " ", sub("_", " - ", organism_raw))

# Count the number of retained SLIT homologs per species.
df_counts <- as.data.table(table(extracted_species))
setnames(df_counts, c("extracted_species", "N"), c("target_species", "count"))

# Match the plot order to the species reference table.
ordered_levels <- paste(species_list$Target_Clade, species_list$Organism_Name, sep = " - ")
df_counts[, target_species := factor(target_species, levels = ordered_levels)]

# Draw the homolog count bar plot.
plot_homolog_counts <- ggplot(df_counts, aes(x = target_species, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", alpha = 0.8) +
  coord_flip() +
  scale_y_continuous(breaks = function(x) seq(0, max(x, na.rm = TRUE), by = 1)) +
  theme_minimal() +
  labs(
    title = "Counts of SLIT Homologs per Species",
    x = "Clade and Species",
    y = "Number of Homologs"
  ) +
  theme(
    axis.text.y = element_text(face = "italic", size = 10),
    panel.grid.minor.x = element_blank()
  )

# Print and save the plot.
print(plot_homolog_counts)

ggsave(
  filename = output_png,
  plot = plot_homolog_counts,
  width = 10,
  height = 6,
  dpi = 300, 
  bg = "white"
)