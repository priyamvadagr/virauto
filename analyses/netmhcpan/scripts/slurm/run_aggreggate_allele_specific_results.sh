#!/bin/bash
#SBATCH -J aggregate_netmhc
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -t 1-00:00
#SBATCH --mem=100G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu
#SBATCH --cluster=htc

# ===================================================
# Setup and Timestamped Logging
# ===================================================
TYPE="Type1_NR"
CLASS="all"
CHUNK='chunk_3'
K_MER=9

# Directory for all logs
LOGDIR="/ix/djishnu/Priyamvada/virauto/analyses/virscan/logs"
mkdir -p "$LOGDIR"

# Generate timestamp for unique log naming
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

# Construct log file paths
OUTFILE="${LOGDIR}/${TYPE}_${CLASS}_${K_MER}_aggregate_${TIMESTAMP}.out"
ERRFILE="${LOGDIR}/${TYPE}_${CLASS}_${K_MER}_aggregate_${TIMESTAMP}.err"


# Redirect all output (stdout + stderr) to the combined log file
exec > "$OUTFILE" 2>&1

# ===================================================
# Job Metadata and Start Banner
# ===================================================
echo "=================================================="
echo "      NetMHCpan Aggregation Job Started"
echo "=================================================="
echo "Job ID:        $SLURM_JOB_ID"
echo "Timestamp:     $TIMESTAMP"
echo "Node:          $SLURMD_NODENAME"
echo "MHC Type:      $TYPE"
echo "MHC Class:     $CLASS"
echo "Peptide k-mer: $K_MER"
echo "Output Log:    $OUTFILE"
echo "Error Log:     $ERRFILE"
echo "Started at:    $(date)"
echo "=================================================="
echo

# ===================================================
# Run Python Aggregation Script
# ===================================================
python /ix/djishnu/Priyamvada/virauto/analyses/virscan/python/aggregate_allele_specific_results.py \
    --type $TYPE \
    --class $CLASS \
    --chunk $CHUNK \
    --kmer $K_MER

# ===================================================
# Job Completion Banner
# ===================================================
echo
echo "=================================================="
echo "      NetMHCpan Aggregation Job Completed"
echo "=================================================="
echo "Finished at:   $(date)"
echo "Job ID:        $SLURM_JOB_ID"
echo "Output saved:  $OUTFILE"
echo "=================================================="
