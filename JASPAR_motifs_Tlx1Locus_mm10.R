#---------------------------
# Install & load dependencies
#---------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("JASPAR2024", "TFBSTools", "Biostrings", "RSQLite"))

library(JASPAR2024)
library(TFBSTools)
library(Biostrings)
library(RSQLite)

#---------------------------
# 1. Load your FASTA
#---------------------------
fasta_file <- "/Users/veralaub/Documents/postdoc/collaboration/Maurizio/JASPAR_Tlx1/Tlx1_mm10_+-2kb_up+downstream.fasta"
seqs <- readDNAStringSet(fasta_file)
sequence <- seqs[[1]]  # assume one sequence

#---------------------------
# 2. Open the JASPAR2024 SQLite DB
#---------------------------
jaspar_obj <- JASPAR2024()
conn <- dbConnect(SQLite(), db(jaspar_obj))

#---------------------------
# 3. Find all motifs whose NAME starts with "Klf"
#---------------------------
query <- "
SELECT BASE_ID, VERSION, NAME
FROM MATRIX
WHERE NAME LIKE 'Klf%'
"
motif_df <- dbGetQuery(conn, query)

if(nrow(motif_df) == 0) {
  stop("No Klf motifs found in the JASPAR2024 database.")
}

#---------------------------
# 4. Scan sequence for each PWM
#---------------------------
all_hits <- lapply(seq_len(nrow(motif_df)), function(i) {
  row <- motif_df[i, ]
  base_id <- row$BASE_ID
  tf_name <- row$NAME
  
  # getMatrixByID expects only the BASE_ID string
  pfm <- getMatrixByID(conn, ID = base_id)
  pwm <- toPWM(pfm, type = "log2probratio")
  
  # scan the sequence
  matches <- searchSeq(pwm, sequence, min.score = "80%", strand = "*")
  df <- as.data.frame(matches)
  
  # add motif metadata
  df$BASE_ID <- base_id
  df$TF_NAME <- tf_name
  df
})

#---------------------------
# 5. Combine & save
#---------------------------
hits_df <- do.call(rbind, all_hits)
hits_df <- hits_df[, c("BASE_ID","TF_NAME","start","end","strand","score","sequence")]

write.csv(hits_df, "/Users/veralaub/Documents/postdoc/collaboration/Maurizio/JASPAR_Tlx1/Klf_family_motif_hits_mm10.csv", row.names = FALSE)
print(head(hits_df))


#---------------------------
# 6. Disconnect
#---------------------------
dbDisconnect(conn)


################################################################################
########################### Visualization of results ###########################
#---------------------------
# Load packages
#---------------------------
library(dplyr)
library(ggplot2)

#---------------------------
# 0. Filter high-confidence hits
#---------------------------
filtered_hits <- hits_df %>%
  filter(relScore >= 0.85)

#---------------------------
# 1. Define Klf order: mouse then human for each family member
#---------------------------
klf_order <- c(
  "Klf1","KLF1",
  "Klf2","KLF2",
  "Klf3","KLF3",
  "Klf4","KLF4",
  "Klf5","KLF5",
  "Klf6","KLF6",
  "Klf7","KLF7",
  "Klf8","KLF8",
  "Klf9","KLF9",
  "Klf10","KLF10",
  "Klf11","KLF11",
  "Klf12","KLF12",
  "Klf13","KLF13",
  "Klf14","KLF14"
)

# Apply ordering to TF_NAME factor
filtered_hits$TF_NAME <- factor(filtered_hits$TF_NAME, levels = klf_order)

#---------------------------
# 2. Summarize hits per TF
#---------------------------
summary_df <- filtered_hits %>%
  group_by(TF_NAME) %>%
  summarize(
    n_hits = n(),
    median_score = median(relScore),
    mean_score = mean(relScore)
  )

print(summary_df)

#---------------------------
# 3. Prepare hit counts for labels
#---------------------------
hit_counts <- filtered_hits %>%
  group_by(TF_NAME) %>%
  summarize(n_hits = n())

#---------------------------
# 4. Create boxplot
#---------------------------
p <- ggplot(filtered_hits, aes(x = TF_NAME, y = relScore, fill = TF_NAME)) +
  geom_boxplot(outlier.shape = 16, alpha = 0.7) +
  geom_text(
    data = hit_counts,
    aes(x = TF_NAME, y = max(filtered_hits$relScore) + 0.02, label = n_hits),
    inherit.aes = FALSE
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Distribution of high-confidence Klf motif hits +/-2kb of Tlx1 locus (mouse & human)",
    x = "Klf TF",
    y = "Relative PWM score"
  ) +
  guides(fill = FALSE)

# Display plot
print(p)

#---------------------------
# 5. Save plot
#---------------------------
output_folder <- dirname("/Users/veralaub/Documents/postdoc/collaboration/Maurizio/JASPAR_Tlx1/Klf_family_motif_hits_mm10.csv")
ggsave(filename = file.path(output_folder, "Klf_motif_Tlx1Locus_boxplot_mouse_human.png"),
       plot = p, width = 12, height = 6, dpi = 300)
