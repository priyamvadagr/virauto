# Get all significant SNPs at thresholds of 0.001, 0.005, 0.01, and 0.05 for MS, T1D, and SLE
library(data.table)

gwas_dir <- '/ix/djishnu/Priyamvada/people/Varun/GWAS_data'
out_dir <- '/ix/djishnu/Priyamvada/auto_immune/1000G_LD/gwas_snps'
disease_names <- c('MS', 'SLE', 'T1D')
file_names <- c('MS_summary_stats.tsv.bgz', 'SLE_summary_stats.tsv.bgz', 'E4_DM1_UKBB_gwas.tsv.bgz')
thresh <- c(0.001, 0.005, 0.01, 0.05)

for (f in seq_along(file_names)) {
  file_path <- file.path(gwas_dir, file_names[f])
  disease <- disease_names[f]
  
  # Read entire GWAS file
  gwas_df <- fread(file_path)
  
  # Create output directory if it doesn't exist
  dir.create(file.path(out_dir, disease), recursive = TRUE, showWarnings = FALSE)
  
  for (t in thresh) {
    # Filter by -log10(p-value)
    gwas_df_flt <- gwas_df[neglog10_pval_EUR >= -log10(t)]
    
    # Construct CHR:POS:REF:ALT ID
    gwas_df_flt[, ID := paste0(chr, ":", pos, ":", ref, ":", alt)]
    
    # Write to file
    out_file <- file.path(out_dir, disease, paste0(disease, "_", t, ".txt"))
    fwrite(gwas_df_flt[, .(ID)], file = out_file, col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
}

