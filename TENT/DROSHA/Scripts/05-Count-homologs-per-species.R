suppressWarnings(suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Biostrings)
}))

script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

get_env <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (!is.na(value) && nzchar(value)) value else default
}

binomial_label <- function(x) {
  sub("^([^ ]+ [^ ]+).*", "\\1", x)
}

scripts_dir <- script_dir()
base_dir <- normalizePath(file.path(scripts_dir, ".."), mustWork = TRUE)
output_prefix <- get_env("PIPELINE_OUTPUT_PREFIX", "SLIT")

species_file <- file.path(scripts_dir, "species.txt")
fasta_path <- get_env("PIPELINE_FILTERED_FASTA", file.path(scripts_dir, paste0(output_prefix, "-homologs.fasta")))
output_png <- get_env("PIPELINE_HOMOLOG_COUNTS_PNG", file.path(scripts_dir, paste0(output_prefix, "_Homolog_Counts.png")))

species_list <- fread(species_file, header = TRUE)
species_list[, Organism_Label := binomial_label(Organism_Name)]

if (!file.exists(fasta_path) || file.info(fasta_path)$size == 0) {
  stop("Missing or empty FASTA file: ", fasta_path)
}

seqs <- readAAStringSet(fasta_path, format = "fasta")

organism_raw <- sub("_H[0-9]+_.*", "", names(seqs))
extracted_species <- gsub("_", " ", sub("_", " - ", organism_raw))

df_counts <- as.data.table(table(extracted_species))
setnames(df_counts, c("extracted_species", "N"), c("target_species", "count"))

ordered_levels <- paste(species_list$Target_Clade, species_list$Organism_Label, sep = " - ")
df_counts[, target_species := factor(target_species, levels = ordered_levels)]

plot_homolog_counts <- ggplot(df_counts, aes(x = target_species, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", alpha = 0.8) +
  coord_flip() +
  scale_y_continuous(breaks = function(x) seq(0, max(x, na.rm = TRUE), by = 1)) +
  theme_minimal() +
  labs(
    title = paste("Counts of", output_prefix, "Homologs per Species"),
    x = "Clade and Species",
    y = "Number of Homologs"
  ) +
  theme(
    axis.text.y = element_text(face = "italic", size = 10),
    panel.grid.minor.x = element_blank()
  )

print(plot_homolog_counts)

ggsave(
  filename = output_png,
  plot = plot_homolog_counts,
  width = 10,
  height = 6,
  dpi = 300,
  bg = "white"
)
