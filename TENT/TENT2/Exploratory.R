library(data.table)
library(ggplot2)
source("/rna/liha/tools/ggplot_theme/set_ggplot_theme.R")
source("/rna/liha/phylogenomics_practice/TENT/TENT2/Exploratory-helper.R")

base_dir <- "/rna/liha/phylogenomics_practice/TENT/TENT2"
figure_dir <- file.path(base_dir, "Figures")
species_info <- fread(file.path(base_dir, "Scripts", "species_info.tsv"))
blast_hits <- fread(file.path(base_dir, "BLASTp-results", "TENT2_BLASTp_results_isoform_resolved.tsv"))
reciprocal_result <- fread(file.path(base_dir, "BLASTp-results", "TENT2_reciprocal_blast_result.tsv"))

human_colors <- c("Homo sapiens" = "red", "Other" = "grey40")
highlighted_species <- c("Homo sapiens", "Bolinopsis microptera", "Rhizophagus irregularis", "Oryza sativa")
species_colors <- c("Homo sapiens" = "red", "Bolinopsis microptera" = "blue", "Rhizophagus irregularis" = "forestgreen", "Oryza sativa" = "orange", "Other" = "grey70")
paralogs_of_interest <- c("TENT2", "TENT4A", "TENT4B", "Tut4", "Tut7")

# Shows the number of nuclear RefSeq proteins available per sampled species across broad taxonomic groups.
plot_refseq_counts(species_info, file.path(figure_dir, "TENT2_RefSeq_by_Taxonomy_dotplot.png"))

# Compares forward-BLAST percent identity among taxonomic groups ordered from close to distant from humans.
plot_identity_by_taxonomic_distance(blast_hits, file.path(figure_dir, "TENT2_pident_by_taxonomic_distance.png"))

# Compares forward-BLAST bitscores among taxonomic groups ordered from close to distant from humans.
plot_bitscore_by_taxonomic_distance(blast_hits, file.path(figure_dir, "TENT2_bitscore_by_taxonomic_distance.png"))

# Compares forward-BLAST query coverage among taxonomic groups ordered from close to distant from humans.
plot_qcov_by_taxonomic_distance(blast_hits, file.path(figure_dir, "TENT2_qcov_by_taxonomic_distance.png"))

# Shows whether forward hits with stronger significance also cover more of TENT2, with human hits highlighted.
plot_qcov_vs_evalue(blast_hits, "Homo sapiens", human_colors, file.path(figure_dir, "TENT2_BLASTp_query_coverage_vs_evalue_human.png"))

# Shows the relationship between forward-hit percent identity and significance, with human hits highlighted.
plot_pident_vs_evalue(blast_hits, "Homo sapiens", human_colors, file.path(figure_dir, "TENT2_BLASTp_pident_vs_evalue.png"))

# Shows query coverage against forward-hit significance while highlighting four representative species.
plot_qcov_vs_evalue(blast_hits, highlighted_species, species_colors, file.path(figure_dir, "TENT2_BLASTp_query_coverage_vs_evalue.png"))

# Compares forward-hit identity and significance among reciprocal TENT2, other-TENT, and non-TENT assignments.
plot_forward_hits_by_reciprocal_category(blast_hits, reciprocal_result, file.path(figure_dir, "TENT2_forward_pident_vs_evalue_by_category.png"))

# Displays the human reverse-BLAST candidates for one forward accession and highlights its strongest TENT assignments.
plot_reverse_bitscore_profile(reciprocal_result, "XP_642359.1", file.path(figure_dir, "TENT2_reverse_bitscore_XP_642359.1.png"))

# Shows the bitscore advantage of TENT2 over the second-best distinct reverse-hit group when TENT2 ranks first.
tent2_gap <- plot_tent2_bitscore_gap(blast_hits, reciprocal_result, file.path(figure_dir, "TENT2_RBH_bitscore_gap_histogram.png"))

# Compares reciprocal bitscores across human TENT paralogs for non-human accessions ordered by taxonomy.
plot_tent_family_dispersion(blast_hits, reciprocal_result, file.path(figure_dir, "TENT_family_bitscore_dispersion_by_taxonomy.png"))

# Counts Diphodia accessions whose top reciprocal TENT paralog is TENT4A or TENT4B in each taxonomy.
tent4_high <- plot_top_paralog_counts(blast_hits, reciprocal_result, c("TENT4A", "TENT4B"), "forestgreen", "Number of Accessions with TENT4A/TENT4B as Top Reciprocal BLAST Hit, by Taxonomy", file.path(figure_dir, "TENT4AB_high_count_by_taxonomy.png"))

# Counts Diphodia accessions whose top reciprocal TENT paralog is Tut4 or Tut7 in each taxonomy.
tut47_high <- plot_top_paralog_counts(blast_hits, reciprocal_result, c("Tut4", "Tut7"), "orchid", "Number of Accessions with Tut4/Tut7 as Top Reciprocal BLAST Hit, by Taxonomy", file.path(figure_dir, "Tut4_7_high_count_by_taxonomy.png"))
print(tut47_high[, .(qseqid, top_paralog, bitscore, requested_species, sscinames, Taxonomy)])
print(tut47_high[, .N, by = Taxonomy][order(-N)])
print(tut47_high[, .N, by = requested_species][order(-N)])

# Shows the presence or absence of each selected top-hit TENT paralog across sampled Diphodia species.
plot_paralog_presence_heatmap(blast_hits, reciprocal_result, paralogs_of_interest, file.path(figure_dir, "Diphodia_paralog_presence_heatmap.png"))
