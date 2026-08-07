# Shared data-preparation helpers ------------------------------------------------

# Assign each detailed taxonomy path to the broad groups used in the figures.
taxonomy_group_order <- c("Gnathostomes", "Basal chordates", "Echinoderms + hemichordates", "Protostomes", "Deep Metazoa", "Other Holozoa", "Holomycota", "Apusomonadida", "Amoebozoa", "Diphodia", "Asgard Archaea (outgroup)")

assign_taxonomy_group <- function(taxonomy, requested_species = rep(NA_character_, length(taxonomy))) {
  group <- rep(NA_character_, length(taxonomy))

  group[grepl("^Obazoa-Holozoa-Metazoa", taxonomy)] <- "Deep Metazoa"
  group[grepl("^Obazoa-Holozoa-Metazoa-Protostomia", taxonomy)] <- "Protostomes"
  group[grepl("^Obazoa-Holozoa-Metazoa-(Echinodermata|Hemichordata)", taxonomy)] <- "Echinoderms + hemichordates"
  group[grepl("^Obazoa-Holozoa-Metazoa-Chordata-(Cephalochordata|Tunicata|Cyclostomata)", taxonomy)] <- "Basal chordates"
  group[grepl("^Obazoa-Holozoa-Metazoa-Chordata-(Chondrichthyes|Sarcopterygii)", taxonomy)] <- "Gnathostomes"
  group[requested_species == "Homo sapiens"] <- "Gnathostomes"
  group[grepl("^Obazoa-Holozoa-(Choanoflagellatea|Filasterea|Ichthyosporea)", taxonomy)] <- "Other Holozoa"
  group[grepl("^Obazoa-Holomycota", taxonomy)] <- "Holomycota"
  group[grepl("^Obazoa-Apusomonadida", taxonomy)] <- "Apusomonadida"
  group[grepl("^Amoebozoa", taxonomy)] <- "Amoebozoa"
  group[grepl("^(Discoba|SAR|Archaeplastida|Cryptista|Haptophyta)", taxonomy)] <- "Diphodia"
  group[grepl("^Asgard_Archaea", taxonomy)] <- "Asgard Archaea (outgroup)"

  factor(group, levels = taxonomy_group_order)
}

# Collapse human TENT isoforms by gene name while retaining non-TENT accessions separately.
collapse_reciprocal_hits <- function(reciprocal_result) {
  collapsed <- copy(reciprocal_result)
  collapsed[, group_key := ifelse(!is.na(is_human_TENT) & is_human_TENT == TRUE & !is.na(human_TENT) & human_TENT != "", human_TENT, saccver)]
  collapsed <- collapsed[, .(bitscore = max(bitscore)), by = .(qseqid, group_key, human_TENT)]

  collapsed <- collapsed[order(qseqid, -bitscore, group_key)]
  collapsed[, rank := seq_len(.N), by = qseqid]

  collapsed
}

# Build a ranked table of TENT paralogs for forward hits with at least two TENT matches.
prepare_top_paralogs <- function(blast_hits, reciprocal_result) {
  plot_dt <- reciprocal_result[!is.na(is_human_TENT) & is_human_TENT == TRUE & !is.na(human_TENT), .(bitscore = max(bitscore)), by = .(qseqid, human_TENT)]
  multi_qseqids <- plot_dt[, .N, by = qseqid][N >= 2, qseqid]
  plot_dt <- plot_dt[qseqid %in% multi_qseqids]

  # Human forward hits are references rather than comparative observations.
  human_qseqids <- blast_hits[requested_species == "Homo sapiens", unique(saccver)]
  plot_dt <- plot_dt[!qseqid %in% human_qseqids]

  # Attach the species and taxonomy information from the forward search.
  qseqid_taxonomy <- unique(blast_hits[!is.na(saccver), .(qseqid = saccver, Taxonomy, requested_species, sscinames)])
  plot_dt <- merge(plot_dt, qseqid_taxonomy, by = "qseqid", all.x = TRUE)
  plot_dt[, Taxonomy_group := assign_taxonomy_group(Taxonomy, requested_species)]

  plot_dt <- plot_dt[order(qseqid, -bitscore, human_TENT)]
  plot_dt[, rank := seq_len(.N), by = qseqid]
  top_paralog_dt <- plot_dt[rank == 1, .(qseqid, top_paralog = human_TENT)]

  merge(plot_dt, top_paralog_dt, by = "qseqid")
}

# Forward-BLAST plots -----------------------------------------------------------

# Plot the usable RefSeq protein count for every sampled non-human species.
plot_refseq_counts <- function(species_info, output_file) {
  filtered <- copy(species_info)[(NP_ > 10 | XP_ > 10) & requested_species != "Homo sapiens"]
  filtered[, RefSeq_count := NP_ + XP_]
  filtered[, Taxonomy_top := sub("[-_].*", "", Taxonomy)]

  tax_order <- filtered[, .(total = sum(RefSeq_count)), by = Taxonomy_top][order(-total), Taxonomy_top]
  filtered[, Taxonomy_top := factor(Taxonomy_top, levels = tax_order)]

  png(output_file, width = 1200, height = 700)
  print(ggplot(filtered, aes(x = Taxonomy_top, y = RefSeq_count)) +
    geom_jitter(width = 0.15, height = 0, size = 3, alpha = 0.7, color = "black") +
    labs(x = "Taxonomy", y = "RefSeq Protein Accessions per Species (NP_ + XP_)", title = "RefSeq Protein Accessions per Species by Taxonomy") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.text.y = element_text(size = 15), axis.title = element_text(size = 15), plot.title = element_text(size = 20, face = "bold")))
  dev.off()

  invisible(filtered)
}

# Compare forward-hit identity across broad taxonomic distances from humans.
plot_identity_by_taxonomic_distance <- function(blast_hits, output_file) {
  plot_dt <- copy(blast_hits)[!is.na(Taxonomy) & !is.na(pident)]
  plot_dt[, Taxonomy_group := assign_taxonomy_group(Taxonomy, requested_species)]
  plot_dt <- plot_dt[!is.na(Taxonomy_group)]
  plot_dt[, point_group := ifelse(requested_species == "Homo sapiens", "Homo sapiens", "Other")]

  png(output_file, width = 1300, height = 700)
  print(ggplot(plot_dt, aes(x = Taxonomy_group, y = pident)) +
    geom_boxplot(outlier.shape = NA, fill = "grey85", color = "grey40", alpha = 0.6, width = 0.6) +
    geom_jitter(aes(color = point_group), width = 0.15, height = 0, size = 2.2, alpha = 0.8) +
    scale_color_manual(values = c("Homo sapiens" = "red", "Other" = "steelblue"), name = "") +
    labs(x = "Taxonomic group (ordered: close to human → far from human)", y = "Percent Identity (%)", title = "TENT2 BLASTp Hit Percent Identity by Taxonomic Distance from Human") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.text.y = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12)))
  dev.off()

  invisible(plot_dt)
}

# Compare forward-hit bitscores across broad taxonomic distances from humans.
plot_bitscore_by_taxonomic_distance <- function(blast_hits, output_file) {
  plot_dt <- copy(blast_hits)[!is.na(Taxonomy) & !is.na(bitscore)]
  plot_dt[, Taxonomy_group := assign_taxonomy_group(Taxonomy, requested_species)]
  plot_dt <- plot_dt[!is.na(Taxonomy_group)]
  plot_dt[, point_group := ifelse(requested_species == "Homo sapiens", "Homo sapiens", "Other")]

  png(output_file, width = 1300, height = 700)
  print(ggplot(plot_dt, aes(x = Taxonomy_group, y = bitscore)) +
    geom_boxplot(outlier.shape = NA, fill = "grey85", color = "grey40", alpha = 0.6, width = 0.6) +
    geom_jitter(aes(color = point_group), width = 0.15, height = 0, size = 2.2, alpha = 0.8) +
    scale_color_manual(values = c("Homo sapiens" = "red", "Other" = "steelblue"), name = "") +
    labs(x = "Taxonomic group (ordered: close to human → far from human)", y = "BLASTp bitscore", title = "TENT2 BLASTp Hit Bitscore by Taxonomic Distance from Human") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.text.y = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12)))
  dev.off()

  invisible(plot_dt)
}

# Compare forward-hit query coverage across broad taxonomic distances from humans.
plot_qcov_by_taxonomic_distance <- function(blast_hits, output_file) {
  plot_dt <- copy(blast_hits)[!is.na(Taxonomy) & !is.na(qcovs)]
  plot_dt[, Taxonomy_group := assign_taxonomy_group(Taxonomy, requested_species)]
  plot_dt <- plot_dt[!is.na(Taxonomy_group)]
  plot_dt[, point_group := ifelse(requested_species == "Homo sapiens", "Homo sapiens", "Other")]

  png(output_file, width = 1300, height = 700)
  print(ggplot(plot_dt, aes(x = Taxonomy_group, y = qcovs)) +
    geom_boxplot(outlier.shape = NA, fill = "grey85", color = "grey40", alpha = 0.6, width = 0.6) +
    geom_jitter(aes(color = point_group), width = 0.15, height = 0, size = 2.2, alpha = 0.8) +
    scale_color_manual(values = c("Homo sapiens" = "red", "Other" = "steelblue"), name = "") +
    labs(x = "Taxonomic group (ordered: close to human → far from human)", y = "Query coverage (%)", title = "TENT2 BLASTp Query Coverage by Taxonomic Distance from Human") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.text.y = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12)))
  dev.off()

  invisible(plot_dt)
}

# Plot query coverage against forward-BLAST significance with selected species highlighted.
plot_qcov_vs_evalue <- function(blast_hits, highlighted_species, colors, output_file) {
  species_group <- ifelse(blast_hits$requested_species %in% highlighted_species, blast_hits$requested_species, "Other")
  point_col <- colors[species_group]
  point_cex <- ifelse(species_group == "Other", 0.8, 1.2)

  png(output_file, width = 800, height = 600)
  plot(blast_hits$qcovs, -log10(blast_hits$evalue), xlab = "Query Coverage (%)", ylab = expression(-log[10]("E-value")), main = "Query Coverage vs -log10(E-value)", col = point_col, pch = 19, cex = point_cex)
  legend("topleft", legend = names(colors), col = colors, pch = 19, bty = "o", cex = 0.8)
  dev.off()
}

# Plot percent identity against forward-BLAST significance with selected species highlighted.
plot_pident_vs_evalue <- function(blast_hits, highlighted_species, colors, output_file) {
  species_group <- ifelse(blast_hits$requested_species %in% highlighted_species, blast_hits$requested_species, "Other")
  point_col <- colors[species_group]
  point_cex <- ifelse(species_group == "Other", 0.8, 1.2)

  png(output_file, width = 800, height = 600)
  plot(blast_hits$pident, -log10(blast_hits$evalue), xlab = "Percent identity (%)", ylab = expression(-log[10]("E-value")), main = "Percent identity vs -log10(E-value)", col = point_col, pch = 19, cex = point_cex)
  legend("topleft", legend = names(colors), col = colors, pch = 19, bty = "o", cex = 0.8)
  dev.off()
}

# Reciprocal-BLAST plots --------------------------------------------------------

# Compare forward-hit quality among TENT2, other-TENT, and non-TENT reciprocal assignments.
plot_forward_hits_by_reciprocal_category <- function(blast_hits, reciprocal_result, output_file) {
  # Use the strongest reciprocal hit to assign each forward accession to a category.
  reciprocal_best <- reciprocal_result[order(evalue, -bitscore), .SD[1], by = qseqid]
  merged <- merge(blast_hits, reciprocal_best[, .(qseqid, is_human_TENT, human_TENT)], by.x = "saccver", by.y = "qseqid", all.x = TRUE)
  merged[, neg_log10_evalue_fwd := -log10(evalue)]

  merged[, category := "Other"]
  merged[!is.na(is_human_TENT) & is_human_TENT == TRUE, category := "TENT family"]
  merged[!is.na(is_human_TENT) & is_human_TENT == TRUE & human_TENT == "TENT2", category := "TENT2"]
  merged[, category := factor(category, levels = c("Other", "TENT family", "TENT2"))]

  png(output_file, width = 900, height = 700)
  print(ggplot(merged[order(category)], aes(x = pident, y = neg_log10_evalue_fwd, color = category)) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_color_manual(values = c("Other" = "grey70", "TENT family" = "blue", "TENT2" = "red"), name = "") +
    labs(x = "Forward BLAST Percent Identity (%)", y = expression("Forward BLAST"~-log[10]("E-value")), title = "Forward BLAST: Percent Identity vs E-value by Reciprocal Category") +
    theme_minimal() +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12)))
  dev.off()

  invisible(merged)
}

# Show the ranked human reverse hits for one selected forward accession.
plot_reverse_bitscore_profile <- function(reciprocal_result, target_accession, output_file) {
  rev_hits <- copy(reciprocal_result)[qseqid == target_accession]

  # Retain one representative isoform for each identified TENT paralog.
  rev_hits[, group_key := ifelse(!is.na(is_human_TENT) & is_human_TENT == TRUE & !is.na(human_TENT) & human_TENT != "", human_TENT, saccver)]
  rev_hits <- rev_hits[, .SD[which.max(bitscore)], by = group_key][order(-bitscore, saccver)]

  # Label the identified TENT paralogs and distinguish the best reverse hit.
  rev_hits[, hit_label := ifelse(!is.na(is_human_TENT) & is_human_TENT == TRUE & !is.na(human_TENT) & human_TENT != "", paste0(saccver, " (", human_TENT, ")"), saccver)]
  rev_hits[, hit_label := factor(hit_label, levels = unique(hit_label))]
  rev_hits[, highlight := "Other reverse hit"]
  rev_hits[human_TENT == "TENT2", highlight := "TENT2 hit"]
  rev_hits[1, highlight := "Reciprocal best hit"]
  rev_hits[, xpos := seq_len(.N)]
  n_slots <- max(5, nrow(rev_hits))

  png(output_file, width = 1200, height = 600)
  print(ggplot(rev_hits, aes(x = xpos, y = bitscore, fill = highlight)) +
    geom_col(width = 0.6) +
    scale_x_continuous(breaks = rev_hits$xpos, labels = rev_hits$hit_label, limits = c(0.5, n_slots + 0.5)) +
    scale_fill_manual(values = c("Other reverse hit" = "grey60", "Reciprocal best hit" = "red", "TENT2 hit" = "blue"), name = "") +
    labs(x = "Reverse BLAST hit (accession, TENT paralog if identified)", y = "Bitscore", title = paste("Reverse BLAST hits for", target_accession)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9), axis.text.y = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(size = 15, face = "bold"), legend.text = element_text(size = 12)))
  dev.off()

  invisible(rev_hits)
}

# Plot the bitscore gap between a top-ranked TENT2 hit and its nearest competitor.
plot_tent2_bitscore_gap <- function(blast_hits, reciprocal_result, output_file) {
  collapsed <- collapse_reciprocal_hits(reciprocal_result)
  top1 <- collapsed[rank == 1]

  # Keep TENT2-best queries that also have a distinct second reverse-hit group.
  tent2_only_qseqids <- collapsed[, .N, by = qseqid][N == 1, qseqid]
  competitor_qseqids <- top1[human_TENT == "TENT2" & !qseqid %in% tent2_only_qseqids, qseqid]
  sub <- collapsed[qseqid %in% competitor_qseqids & rank <= 2]
  diff_dt <- sub[, .(bitscore_diff = bitscore[rank == 1] - bitscore[rank == 2]), by = qseqid]

  # Exclude human reference accessions from the comparative distribution.
  human_qseqids <- blast_hits[requested_species == "Homo sapiens", unique(saccver)]
  diff_dt <- diff_dt[!qseqid %in% human_qseqids]

  png(output_file, width = 900, height = 600)
  print(ggplot(diff_dt, aes(x = bitscore_diff)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "white", alpha = 0.85) +
    labs(x = "Bitscore difference in reciprocal BLAST (TENT2 best hit − 2nd best hit)", y = "Number of forward BLAST hits", title = "Bitscore difference between TENT2 and Runner-up") +
    theme_minimal() +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(size = 14, face = "bold")))
  dev.off()

  invisible(diff_dt)
}

# Display the reciprocal bitscores of multiple TENT paralogs along a taxonomy-ordered x-axis.
plot_tent_family_dispersion <- function(blast_hits, reciprocal_result, output_file) {
  plot_dt <- prepare_top_paralogs(blast_hits, reciprocal_result)
  plot_dt <- plot_dt[!is.na(Taxonomy_group)]

  # Within each taxonomy, group accessions by their top reciprocal TENT paralog.
  top_paralog_order <- c("TENT2", "mtPAP", "Tut1", "Tut7", "Tut4", "TENT4A", "TENT4B", "TENT5A", "TENT5B", "TENT5C", "TENT5D")
  ordering_dt <- plot_dt[rank == 1, .(qseqid, Taxonomy_group, top_paralog, top_bitscore = bitscore)]
  ordering_dt[, top_paralog := factor(top_paralog, levels = top_paralog_order)]
  ordering_dt <- ordering_dt[order(Taxonomy_group, top_paralog, -top_bitscore, qseqid)]
  qseqid_order <- ordering_dt$qseqid
  plot_dt[, qseqid := factor(qseqid, levels = qseqid_order)]

  # Define one shaded background block for each broad taxonomic group.
  group_blocks <- ordering_dt[, .(xmin = min(match(qseqid, qseqid_order)) - 0.5, xmax = max(match(qseqid, qseqid_order)) + 0.5), by = Taxonomy_group]
  group_colors <- c("Gnathostomes" = "#F7CAC9", "Basal chordates" = "#FADBD8", "Echinoderms + hemichordates" = "#FDE2C5", "Protostomes" = "#FCECC9", "Deep Metazoa" = "#FDE0DC", "Other Holozoa" = "#FDEBC8", "Holomycota" = "#FFF6C8", "Apusomonadida" = "#E3F2CE", "Amoebozoa" = "#CDEEDD", "Diphodia" = "#CDE7F0", "Asgard Archaea (outgroup)" = "#DCD6F0")

  png(output_file, width = 1400, height = 700)
  print(ggplot() +
    geom_rect(data = group_blocks, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Taxonomy_group), alpha = 0.35, show.legend = FALSE) +
    scale_fill_manual(values = group_colors) +
    geom_point(data = plot_dt, aes(x = qseqid, y = bitscore, color = human_TENT), size = 2.5, alpha = 0.9) +
    labs(x = "Forward BLAST accession (ordered by taxonomic distance from human)", y = "Reciprocal BLAST bitscore", color = "TENT family", title = "Bitscore Dispersion Across TENT Family Paralogs, Ordered by Taxonomic Group\n(Non-human accessions with ≥2 distinct TENT family hits)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), axis.text.y = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 11)))
  dev.off()

  invisible(plot_dt)
}

# Count accessions whose strongest TENT-family match belongs to a selected paralog set.
plot_top_paralog_counts <- function(blast_hits, reciprocal_result, paralogs, fill, title, output_file) {
  plot_dt <- prepare_top_paralogs(blast_hits, reciprocal_result)
  selected <- plot_dt[Taxonomy_group == "Diphodia" & rank == 1 & top_paralog %in% paralogs][order(-bitscore)]

  counts <- selected[, .N, by = Taxonomy][order(-N)]
  counts[, Taxonomy := factor(Taxonomy, levels = Taxonomy)]

  png(output_file, width = 1000, height = 700)
  print(ggplot(counts, aes(x = Taxonomy, y = N)) +
    geom_col(fill = fill, alpha = 0.85) +
    labs(x = "Taxonomy", y = "Count", title = title) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11), axis.text.y = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(size = 14, face = "bold")))
  dev.off()

  invisible(selected)
}

# Summarize the presence of selected top-hit TENT paralogs across Diphodia species.
plot_paralog_presence_heatmap <- function(blast_hits, reciprocal_result, paralogs, output_file) {
  plot_dt <- prepare_top_paralogs(blast_hits, reciprocal_result)
  diphodia_top <- plot_dt[Taxonomy_group == "Diphodia" & rank == 1 & top_paralog %in% paralogs]

  # Expand to every species-paralog combination so absences are explicitly plotted.
  presence_dt <- unique(diphodia_top[, .(present = TRUE), by = .(requested_species, top_paralog)])
  species_taxonomy <- unique(blast_hits[, .(requested_species, Taxonomy)])
  full_grid <- CJ(requested_species = sort(unique(diphodia_top$requested_species)), top_paralog = paralogs)
  heatmap_dt <- merge(full_grid, presence_dt, by = c("requested_species", "top_paralog"), all.x = TRUE)
  heatmap_dt[, presence := ifelse(is.na(present), "No", "Yes")]
  heatmap_dt <- merge(heatmap_dt, species_taxonomy, by = "requested_species", all.x = TRUE)
  heatmap_dt[, top_paralog := factor(top_paralog, levels = paralogs)]

  # Arrange species by their detailed taxonomy within the faceted heatmap.
  species_order <- unique(heatmap_dt[, .(requested_species, Taxonomy)])[order(Taxonomy, requested_species)]
  heatmap_dt[, requested_species := factor(requested_species, levels = species_order$requested_species)]

  png(output_file, width = 600, height = 1000)
  print(ggplot(heatmap_dt, aes(x = top_paralog, y = requested_species, fill = presence)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_manual(values = c("Yes" = "black", "No" = "grey90"), name = "") +
    facet_grid(Taxonomy ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = "Top reciprocal BLAST paralog", y = NULL, title = "Presence/Absence of Top-Hit TENT in Diphodia") +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), axis.title = element_text(size = 14), plot.title = element_text(size = 14, face = "bold"), strip.text.y.left = element_text(angle = 0, size = 12, face = "bold"), strip.placement = "outside", panel.spacing = grid::unit(0.3, "lines")))
  dev.off()

  invisible(heatmap_dt)
}
