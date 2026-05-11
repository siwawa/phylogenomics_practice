library(data.table)
library(Biostrings) 

base_dir <- "/rna/liha/phylogenomics_practice/UCE"
fasta_dir <- file.path(base_dir, "harvest_UCE", "05-fasta-output")
lastz_dir <- file.path(base_dir, "harvest_UCE", "04-lastz-output")
species_file <- file.path(base_dir, "harvest_UCE", "00-species", "species.tsv")
out_uce_dir <- file.path(base_dir, "Unaligned_fasta") 

if(!dir.exists(out_uce_dir)) dir.create(out_uce_dir, recursive = TRUE)

species_dt <- fread(species_file, sep = "\t", stringsAsFactors = FALSE)
species_dt[, file_name := gsub(" ", "_", Organism_Name)]



MIN_IDENTITY <- 82.0 

get_uce_identities <- function(f_name) {
  lastz_path <- file.path(lastz_dir, paste0("metazoan_probes-70.fasta_v_", f_name, ".lastz.clean"))
  
  if (file.exists(lastz_path)) {
    dt <- fread(lastz_path, header = FALSE, fill = TRUE, stringsAsFactors = FALSE)
    
    if (!is.null(dt) && nrow(dt) > 0) {
      dt[, V15_num := as.numeric(V15)]
      filtered_dt <- dt[!is.na(V15_num) & V15_num >= MIN_IDENTITY]
      
      if (nrow(filtered_dt) > 0) {
        filtered_dt[, UCE_ID := regmatches(V7, regexpr("uce-[0-9]+", V7))]
        
        # 한 종에서 같은 UCE에 두 번 매칭된 경우(Paralog 등), 가장 높은 점수(max)를 취함
        res <- filtered_dt[, .(Identity = max(V15_num)), by = UCE_ID]
        res[, file_name := f_name]
        return(res)
      }
    }
  }
  # 매칭된 것이 없으면 빈 데이터 테이블 반환
  return(data.table(UCE_ID = character(0), Identity = numeric(0), file_name = character(0))) 
}

list_of_dts <- lapply(species_dt$file_name, get_uce_identities)
long_dt <- rbindlist(list_of_dts) 




# long_dt에 Target_Clade 정보 병합 (나중에 Deep/Shallow 분기 계산용)
long_dt <- merge(long_dt, species_dt[, .(file_name, Target_Clade)], by = "file_name", all.x = TRUE)

# 3. 변경됨: long_dt를 기반으로 각 종별 UCE 개수 카운트 업데이트
sp_counts <- long_dt[, .(UCE_Count = .N), by = file_name]
species_dt <- merge(species_dt, sp_counts, by = "file_name", all.x = TRUE)
species_dt[is.na(UCE_Count), UCE_Count := 0] # NA를 0으로 처리 
 




# 4. Generate Summary Table by Target Clade
summary_table <- species_dt[, .(
  Total_Species = .N,                              
  Species_with_UCEs = sum(UCE_Count > 0),          
  Avg_UCEs_per_Species = round(mean(UCE_Count), 1),
  Total_UCE_Records = sum(UCE_Count)               
), by = Target_Clade]

print(summary_table[order(-Total_UCE_Records)]) 
fwrite(summary_table, file = file.path(base_dir, "06-summary_table.tsv"), sep = "\t", row.names = FALSE)







# 6. Generate UCE-centric Statistics 
uce_summary <- long_dt[, .(
  Total_Species_Count = .N, 
  Deep_Species_Count = sum(Target_Clade != "Euteleostomi"), 
  Shallow_Species_Count = sum(Target_Clade == "Euteleostomi"),
  Median_Identity = round(median(Identity, na.rm = TRUE), 2) 
), by = UCE_ID]

# Sort the dataframe
setorder(uce_summary, -Deep_Species_Count, -Total_Species_Count)
setorder(uce_summary, -Shallow_Species_Count, -Total_Species_Count)

uce_stats_file <- file.path(base_dir, "uce_stats_all.tsv")
fwrite(uce_summary, file = uce_stats_file, sep = "\t", row.names = FALSE)

uce_stats_file_top10 <- file.path(base_dir, "uce_stats_top10.tsv") 
fwrite(head(uce_summary, 10), file = uce_stats_file_top10, sep = "\t", row.names = FALSE) 






# 7. Extract FASTA for Top UCEs 
target_uces <- head(uce_summary, 10)$UCE_ID 

for (target_uce in target_uces) {
  out_file <- file.path(out_uce_dir, paste0(target_uce, ".fasta"))
  if (file.exists(out_file)) {
    file.remove(out_file)
  }
} 

valid_species <- species_dt[UCE_Count > 0]  

for (i in 1:nrow(valid_species)) {
  sp_name <- tolower(valid_species$file_name[i])
  sp_clade <- gsub(" ", "_", valid_species$Target_Clade[i])
  fasta_path <- file.path(fasta_dir, paste0(sp_name, ".fasta"))
  
  if (file.exists(fasta_path)) {
    sp_seqs <- readDNAStringSet(fasta_path)
    sp_headers <- names(sp_seqs)
    extracted_ids <- regmatches(sp_headers, regexpr("uce-[0-9]+", sp_headers))
    
    for (target_uce in target_uces) {
      idx <- match(target_uce, extracted_ids)
      
      if (!is.na(idx)) {
        target_seq <- sp_seqs[idx]
        names(target_seq) <- paste0(sp_clade, "_", sp_name)
        out_file <- file.path(out_uce_dir, paste0(target_uce, ".fasta"))
        writeXStringSet(target_seq, filepath = out_file, append = TRUE)
      }
    }
  }
}

cat(sprintf("Successfully saved %d UCE alignments in %s\n", length(target_uces), out_uce_dir))