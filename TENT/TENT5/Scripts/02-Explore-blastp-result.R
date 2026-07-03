suppressWarnings(suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(dplyr)
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

read_fasta_lengths <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Missing or empty query FASTA file: ", path)
  }

  lines <- readLines(path, warn = FALSE)
  lengths <- integer()
  current_length <- 0L
  seen_header <- FALSE

  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line)) next

    if (startsWith(line, ">")) {
      if (seen_header) lengths <- c(lengths, current_length)
      current_length <- 0L
      seen_header <- TRUE
    } else if (seen_header) {
      current_length <- current_length + nchar(gsub("\\s+", "", line))
    }
  }

  if (seen_header) lengths <- c(lengths, current_length)
  if (length(lengths) == 0 || any(lengths <= 0)) {
    stop("No non-empty protein sequences found in query FASTA file: ", path)
  }

  lengths
}

regex_from_suffix <- function(suffix) {
  glob2rx(paste0("*", suffix))
}

scripts_dir <- script_dir()
base_dir <- normalizePath(file.path(scripts_dir, ".."), mustWork = TRUE)
blast_results_dir <- file.path(base_dir, "Blast", "Results")
query_name <- get_env("PIPELINE_QUERY_NAME", "SLIT1-2-3")
output_prefix <- get_env("PIPELINE_OUTPUT_PREFIX", "SLIT")
query_fasta <- file.path(base_dir, "Blast", "Query", paste0(query_name, ".fasta"))
query_lengths <- read_fasta_lengths(query_fasta)
reference_length <- mean(query_lengths)
reference_label <- sprintf(
  "Mean query length from %s (%.1f aa; n=%d)",
  basename(query_fasta),
  reference_length,
  length(query_lengths)
)
reference_species <- get_env("PIPELINE_REFERENCE_SPECIES", "Homo sapiens")

species_list <- fread(file.path(scripts_dir, "species.txt"), header = TRUE)
species_list[, Organism_Label := binomial_label(Organism_Name)]

file_list <- list.files(
  path = blast_results_dir,
  pattern = regex_from_suffix(paste0(query_name, ".txt")),
  recursive = TRUE,
  full.names = TRUE
)

if (length(file_list) == 0) {
  stop("No BLAST result files found for query name '", query_name, "' under ", blast_results_dir)
}

names(file_list) <- file_list
blast_cols <- c("qseqid", "saccver", "pident", "length", "qcovs", "evalue", "bitscore", "stitle")

safe_fread <- function(f) {
  if (file.info(f)$size == 0) return(NULL)
  fread(f, col.names = blast_cols, fill = TRUE, sep = "\t", header = FALSE)
}

df <- rbindlist(lapply(file_list, safe_fread), idcol = "filepath", fill = TRUE)
if (nrow(df) == 0) {
  stop("All BLAST result files are empty for query name '", query_name, "'.")
}

df[, source_file := basename(filepath)]
df[, organism := basename(dirname(filepath))]

df[, actual_species := sub(".*\\[(.*?)\\].*", "\\1", stitle)]
df[, target_species := gsub("_", " ", sub("^[^_]+_", "", organism))]
df_clean <- df[actual_species == target_species]

if (nrow(df_clean) == 0) {
  stop("No BLAST rows matched the target species names after filtering.")
}

df_bit <- df_clean[bitscore >= 300]
p1 <- ggplot(df_bit, aes(x = bitscore)) +
  geom_histogram(fill = "steelblue", color = "black", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution of Bit Scores (>= 300)",
    subtitle = paste("Total filtered alignments:", nrow(df_bit)),
    x = "Bit Score",
    y = "Count"
  )
print(p1)

df_plot <- copy(df_clean[evalue < 1e-10])
df_plot[evalue == 0, evalue := 1e-300]

p3 <- ggplot(df_plot, aes(x = -log10(evalue))) +
  geom_histogram(fill = "steelblue", color = "black", bins = 30) +
  theme_minimal() +
  labs(
    title = "Distribution of E-Values (< 1e-10)",
    subtitle = paste("Total filtered alignments:", nrow(df_plot), "(E-value=0 adjusted to 1e-300)"),
    x = "-log10(E-Value)",
    y = "Count"
  )
print(p3)
ggsave(filename = file.path(scripts_dir, "E-value_Distribution.png"), plot = p3, width = 8, height = 5, bg = "white")

df_plot4 <- copy(df_clean[evalue < 1e-80])
df_plot4[evalue == 0, evalue := 1e-300]

p4 <- ggplot(df_plot4, aes(x = length, y = -log10(evalue))) +
  geom_point(color = "navy", alpha = 0.4, size = 2) +
  geom_point(data = df_plot4[target_species == reference_species], color = "red", alpha = 0.9, size = 3) +
  theme_minimal() +
  labs(
    title = "E-value vs Aligned Length (E-value < 1e-80)",
    subtitle = paste("Dashed line:", reference_label, "| Red dots:", reference_species),
    x = "Aligned Length (aa)",
    y = "-log10(E-Value)"
  ) +
  theme(plot.title = element_text(size = 14))

p4 <- p4 + geom_vline(xintercept = reference_length, color = "darkred", linetype = "dashed", linewidth = 0.8)

print(p4)
ggsave(filename = file.path(scripts_dir, "E-value_vs_Aligned_Length.png"), plot = p4, width = 8, height = 5, bg = "white")

df_clean[, target_species := factor(target_species, levels = species_list$Organism_Label)]
df_hits <- df_clean[evalue <= 1e-300]

total_hits <- ggplot(df_hits, aes(x = target_species)) +
  geom_bar(fill = "lightblue", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  coord_flip() +
  labs(
    title = paste("Number of High-Scoring Hits (E-value <= 1e-300) per Organism for", output_prefix),
    y = "Count"
  ) +
  scale_x_discrete(drop = FALSE)
print(total_hits)
ggsave(filename = file.path(scripts_dir, paste0(output_prefix, "_Hit_Counts.png")), plot = total_hits, width = 10, height = 6, bg = "white")

df_hits <- df_hits %>% select(qseqid, saccver, pident, length, qcovs, evalue, bitscore, stitle, source_file, organism, filepath)
fwrite(df_hits, file = file.path(scripts_dir, "Blast-high-scoring-hits.txt"), sep = "\t")
