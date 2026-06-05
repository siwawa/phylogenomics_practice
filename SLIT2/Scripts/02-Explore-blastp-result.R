library(data.table)
library(ggplot2)
library(dplyr)

# ==========================================
# 1. Data Loading & Preprocessing
# ==========================================

# Setup directories and load reference lists
base_dir <- "/rna/liha/phylogenomics_practice/SLIT2/Blast/Results"
species_list <- fread("/rna/liha/phylogenomics_practice/SLIT2/Scripts/species.txt", header = TRUE)

# Retrieve BLAST result files and assign names for idcol extraction
file_list <- list.files(path = base_dir, pattern = "*SLIT1-2-3.txt$", recursive = TRUE, full.names = TRUE)
names(file_list) <- file_list

# Define BLAST outfmt 6 columns
blast_cols <- c("qseqid", "saccver", "pident", "length", "qcovs", "evalue", "bitscore", "stitle")

# Safe read function to skip empty files
safe_fread <- function(f) {
  if (file.info(f)$size == 0) return(NULL) 
  fread(f, col.names = blast_cols, fill = TRUE, sep = "\t", header = FALSE)
}

# Bind all files and extract directory metadata
df <- rbindlist(lapply(file_list, safe_fread), idcol = "filepath", fill = TRUE)
df[, source_file := basename(filepath)]
df[, organism := basename(dirname(filepath))]


# Filter noise (remove clustered DB artifacts by matching species names)
df[, actual_species := sub(".*\\[(.*?)\\].*", "\\1", stitle)]
# Urochordata_Thalia_democratica -> Thalia democratica  
df[, target_species := gsub("_", " ", sub("^[^_]+_", "", organism))]
df_clean <- df[actual_species == target_species]



# ==========================================
# 2. Data Visualization
# ==========================================

# Plot 1: Bit Score Distribution
df_bit <- df_clean[bitscore >= 300]
p1 <- ggplot(df_bit, aes(x = bitscore)) +
  geom_histogram(fill = "steelblue", color = "black", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of Bit Scores (>= 300)",
       subtitle = paste("Total filtered alignments:", nrow(df_bit)),
       x = "Bit Score",
       y = "Count")
print(p1)

# Plot 2: E-value Distribution
df_plot <- df_clean[evalue < 1e-10]
df_plot[evalue == 0, evalue := 1e-300] # Adjust absolute zero values for log scale

p3 <- ggplot(df_plot, aes(x = -log10(evalue))) +
  geom_histogram(fill = "steelblue", color = "black", bins = 30) +
  theme_minimal() +
  labs(title = "Distribution of E-Values (< 1e-10)",
       subtitle = paste("Total filtered alignments:", nrow(df_plot), 
                        "(E-value=0 adjusted to 1e-300)"),
       x = "-log10(E-Value)",
       y = "Count")
print(p3) 

# Plot 3: E-value vs Aligned Length
df_plot4 <- df_clean[evalue < 1e-80]
df_plot4[evalue == 0, evalue := 1e-300]

p4 <- ggplot(df_plot4, aes(x = length, y = -log10(evalue))) +
  geom_point(color = "navy", alpha = 0.4, size = 2) +
  geom_point(data = df_plot4[organism == "Homo sapiens"], 
             color = "red", alpha = 0.9, size = 3) +
  geom_vline(xintercept = 1529, color = "darkred", linetype = "dashed", size = 0.8) +
  theme_minimal() +
  labs(title = "E-value vs Aligned Length (E-value < 1e-80)",
       subtitle = "Dashed line: Human SLIT2 length (1529 aa) | Red dots: Homo sapiens",
       x = "Aligned Length (aa)",
       y = "-log10(E-Value)") +
  theme(plot.title = element_text(size = 14))
print(p4)


# ==========================================
# 3. Organism Hit Counts
# ==========================================

# Align target_species order with the phylogenetic reference list
df_clean[, target_species := factor(target_species, levels = species_list$Organism_Name)]

# Isolate hits with absolute maximum significance
df_hits <- df_clean[evalue <= 1e-300]

# Plot 4: Hit counts for all SLIT proteins
total_hits <- ggplot(df_hits, aes(x = target_species)) +
  geom_bar(fill = "lightblue", color = "black") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  coord_flip() +
  labs(title = "Number of High-Scoring Hits (E-value <= 1e-300) per Organism",
       y = "Count") + 
  scale_x_discrete(drop = FALSE)  # Display organisms with zero hits 
print(total_hits)
ggsave(filename = "/rna/liha/phylogenomics_practice/SLIT2/Scripts/SLIT_Hit_Counts.png", plot = total_hits, width = 10, height = 6)



df_hits <- df_hits %>% select(qseqid, saccver, pident, length, qcovs, evalue, bitscore, stitle, source_file, organism, filepath)
fwrite(df_hits, file = "/rna/liha/phylogenomics_practice/SLIT2/Scripts/Blast-high-scoring-hits.txt", sep = "\t")












# PLOT 7: Counts of hits per species
library(Biostrings)

fasta_path <- "/rna/liha/phylogenomics_practice/SLIT2/Scripts/SLIT-homologs.fasta"
a <- readAAStringSet(fasta_path, format = "fasta")

# 2. 헤더 파싱 (Clade를 지우지 않고 살려서 라벨로 만듭니다)
organism_raw <- sub("_H[0-9]+_.*", "", names(a))

# 첫 번째 언더바만 " - " 로 바꾸고, 나머지 언더바는 공백으로 바꿈
# 결과: "Cephalochordata_Branchiostoma_belcheri" -> "Cephalochordata - Branchiostoma belcheri"
extracted_species <- gsub("_", " ", sub("_", " - ", organism_raw))

# 3. 카운트
df_counts <- as.data.table(table(extracted_species))
setnames(df_counts, c("extracted_species", "N"), c("target_species", "count"))

# 4. ★ 수정된 핵심 로직: 계통수 레퍼런스 리스트 순서에 맞게 Factor 지정 ★
# species_list에서도 위와 똑같은 형태("Clade - Genus species")로 글자를 합쳐서 레벨을 만듭니다.
ordered_levels <- paste(species_list$Target_Clade, species_list$Organism_Name, sep = " - ")
df_counts[, target_species := factor(target_species, levels = ordered_levels)]

# 5. 막대 그래프 생성
plot_7 <- ggplot(df_counts, aes(x = target_species, y = count)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", alpha = 0.8) +
  coord_flip() + 
  scale_y_continuous(breaks = function(x) seq(0, max(x, na.rm = TRUE), by = 1)) + 
  theme_minimal() +
  labs(title = "Counts of SLIT Homologs per Species",
       x = "Clade & Species",
       y = "Number of Homologs") +
  theme(axis.text.y = element_text(face = "italic", size = 10), 
        panel.grid.minor.x = element_blank()) 

# 그래프 출력
print(plot_7)



