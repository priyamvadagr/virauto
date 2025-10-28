#!/bin/bash
# ===========================================================================
# OPTION 1: Single High-Performance Job on AMD EPYC node (FASTEST for single chunk)
# ===========================================================================
# File: run_single_chunk_amd.sh
#SBATCH -J aggregate_netmhc
#SBATCH --cluster htc
#SBATCH --constraint=genoa         # Use AMD EPYC 9374F nodes (fastest)
#SBATCH -N 1
#SBATCH -c 64                      # Use all cores on AMD node
#SBATCH -t 0-12:00
#SBATCH --mem=512G                 # Use most of the 768GB available
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

# Configuration
TYPE="Type1_NR"
CLASS="all"
CHUNK='chunk_1'
K_MER=9
WORKERS=60                         # Slightly less than cores for overhead
BATCH_SIZE=100                     # Large batch size with ample memory

# Directory for all logs
LOGDIR="/ix/djishnu/Priyamvada/virauto/analyses/virscan/logs"
mkdir -p "$LOGDIR"

# Generate timestamp for unique log naming
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
OUTFILE="${LOGDIR}/${TYPE}_${CLASS}_${K_MER}_aggregate_${TIMESTAMP}.out"

exec > "$OUTFILE" 2>&1

echo "=================================================="
echo "      NetMHCpan Aggregation Job Started"
echo "=================================================="
echo "Job ID:        $SLURM_JOB_ID"
echo "Node:          $SLURMD_NODENAME (AMD EPYC 9374F)"
echo "CPU Cores:     $SLURM_CPUS_PER_TASK"
echo "Memory:        512GB"
echo "Workers:       $WORKERS"
echo "Batch Size:    $BATCH_SIZE"
echo "Chunk:         $CHUNK"
echo "Started at:    $(date)"
echo "=================================================="

python /ix/djishnu/Priyamvada/virauto/analyses/virscan/python/aggregate_allele_specific_results.py \
    --type $TYPE \
    --class $CLASS \
    --chunk $CHUNK \
    --kmer $K_MER \
    --workers $WORKERS \
    --batch_size $BATCH_SIZE

echo "=================================================="
echo "Completed at: $(date)"
echo "=================================================="

