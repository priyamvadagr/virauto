#!/bin/bash
# ===========================================================================
# Aggregate NetMHCpan Predictions into Long-Format Parquet (Partitioned by Allele)
# ===========================================================================
# File: run_annotate_allele_binding.sh
#SBATCH -J annotate_allele_binding
##SBATCH --cluster smp
#SBATCH --constraint=genoa         # Use AMD EPYC 9374F nodes (fastest CPU)
#SBATCH -N 1
#SBATCH -c 64                      # 64 CPU cores on AMD node
#SBATCH -t 0-12:00
#SBATCH --mem=100G                 # Use high memory for large batch processing
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

# ===========================================================================
# Configuration
# ===========================================================================
CLASS="Class_I"
TYPE="A"                            # MHC class filter (A/B/C/all)
K_MER=9                            # Peptide length
BATCH_SIZE=100                     # Files per write batch
FORCE_REPROCESS=true              # Set to true to ignore manifest

DATA_DIR="/ix/djishnu/Priyamvada/virauto/results/netmhcpan/virscan/"
RESULTS_DIR="/ix/djishnu/Priyamvada/virauto/results/mimicry_analysis/virscan/"

PYTHON_SCRIPT="/ix/djishnu/Priyamvada/virauto/analyses/mimicry_analysis/python/annotate_delta_binding.py"
LOGDIR="/ix/djishnu/Priyamvada/virauto/analyses/mimicry_analysis/logs"
mkdir -p "$LOGDIR"

WORKERS=60                         # Use slightly fewer than total cores

# ===========================================================================
# Logging setup
# ===========================================================================
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTFILE="${LOGDIR}/delta_binding_${CLASS}_k${K_MER}_${TIMESTAMP}.out"
exec > "$OUTFILE" 2>&1

echo "=================================================="
echo "       NetMHCpan Aggregation by Allele Started"
echo "=================================================="
echo "Job ID:        $SLURM_JOB_ID"
echo "Node:          $SLURMD_NODENAME (AMD EPYC 9374F)"
echo "CPU Cores:     $SLURM_CPUS_PER_TASK"
echo "Workers:       $WORKERS"
echo "Batch Size:    $BATCH_SIZE"
echo "K-mer length:  $K_MER"
echo "Class Filter:  $CLASS"
echo "Started at:    $(date)"
echo "=================================================="

# ===========================================================================
# Command
# ===========================================================================
python "$PYTHON_SCRIPT" \
    --class "$CLASS" \
    --kmer "$K_MER" \
    --hla_type "$TYPE" \
    --netmhc_dir "$DATA_DIR" \
    --out_dir "$RESULTS_DIR" \
    --batch_size "$BATCH_SIZE" \
    --workers "$WORKERS" \
    $( $FORCE_REPROCESS && echo "--force_reprocess" )

# ===========================================================================
# Completion
# ===========================================================================
echo "=================================================="
echo "Completed at: $(date)"
echo "=================================================="
