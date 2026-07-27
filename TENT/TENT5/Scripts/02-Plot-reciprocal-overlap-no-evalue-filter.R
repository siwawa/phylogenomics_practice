suppressWarnings(suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
}))

base_dir <- "/rna/liha/phylogenomics_practice/TENT/TENT5"
scripts_dir <- file.path(base_dir, "Scripts")
results_dir <- file.path(base_dir, "Blast", "Results")
species_file <- file.path(scripts_dir, "species.txt")
rb_file <- file.path(scripts_dir, "Blast-high-scoring-hits-rb-no-evalue-filter.txt")
output_png <- file.path(scripts_dir, "E-value_Distribution_Overlap_no_filter.png")
query_name <- "TENT5A-B-C-D"

binomial_label <- function(x) {
  sub("^([^ ]+ [^ ]+).*", "\\1", x)
}

species_list <- fread(species_file)
species_list[, organism_label := binomial_label(Organism_Name)]
species_lookup <- unique(species_list[, .(
  organism = paste(
    Target_Clade,
    gsub(" ", "_", organism_label, fixed = TRUE),
    sep = "_"
  ),
  target_species = organism_label
)])

add_species_labels <- function(x) {
  x <- copy(x)
  x[, actual_species := binomial_label(sub(".*\\[(.*?)\\].*", "\\1", stitle))]
  x[species_lookup, target_species := i.target_species, on = .(organism)]
  x
}

blast_cols <- c(
  "qseqid", "saccver", "pident", "length",
  "qcovs", "evalue", "bitscore", "stitle"
)
result_files <- list.files(
  results_dir,
  pattern = glob2rx(paste0("*", query_name, ".txt")),
  recursive = TRUE,
  full.names = TRUE
)
names(result_files) <- result_files

read_result <- function(path) {
  if (file.info(path)$size == 0) return(NULL)
  fread(path, sep = "\t", header = FALSE, col.names = blast_cols, fill = TRUE)
}

df <- rbindlist(lapply(result_files, read_result), idcol = "filepath", fill = TRUE)
df[, source_file := basename(filepath)]
df[, organism := basename(dirname(filepath))]
df <- add_species_labels(df)

unmatched <- unique(df[is.na(target_species), organism])
if (length(unmatched) > 0) {
  stop("No species.txt mapping for result directories: ", paste(unmatched, collapse = ", "))
}
df_clean <- df[actual_species == target_species]

if (!file.exists(rb_file) || file.info(rb_file)$size == 0) {
  stop("Missing no-filter reciprocal table: ", rb_file)
}
df_rb <- add_species_labels(fread(rb_file))
df_rb <- df_rb[!is.na(target_species) & actual_species == target_species]

row_keys <- c(
  "qseqid", "saccver", "pident", "length", "qcovs", "evalue",
  "bitscore", "stitle", "source_file", "organism", "filepath"
)
rb_outside_clean <- df_rb[!unique(df_clean[, ..row_keys]), on = row_keys]
if (nrow(rb_outside_clean) > 0) {
  stop("Reciprocal table contains ", nrow(rb_outside_clean), " rows outside df_clean.")
}

df_clean_plot <- copy(df_clean)
df_rb_plot <- copy(df_rb)
df_clean_plot[evalue == 0, evalue := 1e-300]
df_rb_plot[evalue == 0, evalue := 1e-300]
df_clean_plot[, group := "All species-matched forward hits"]
df_rb_plot[, group := "Reciprocal TENT5 hits"]

plot_data <- rbindlist(list(df_clean_plot, df_rb_plot), use.names = TRUE, fill = TRUE)
plot_data[, log_evalue := -log10(evalue)]

p <- ggplot(plot_data, aes(x = log_evalue, fill = group)) +
  geom_histogram(
    position = "identity",
    alpha = 0.55,
    bins = 50,
    color = "black"
  ) +
  scale_fill_manual(values = c(
    "All species-matched forward hits" = "lightblue",
    "Reciprocal TENT5 hits" = "darkblue"
  )) +
  theme_minimal() +
  labs(
    title = "E-value Distribution Without the 1e-10 Filter",
    subtitle = paste0(
      "All species-matched rows: ", nrow(df_clean_plot),
      " | Reciprocal TENT5 rows: ", nrow(df_rb_plot)
    ),
    x = "-log10(E-value)",
    y = "Count",
    fill = NULL
  )

ggsave(output_png, p, width = 10, height = 6, dpi = 300, bg = "white")

message("All species-matched forward rows: ", nrow(df_clean))
message("Species-matched reciprocal rows: ", nrow(df_rb))
message("Overlap plot: ", output_png)
