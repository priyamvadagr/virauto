# virauto

The repository is organized into **modular analyses** (e.g., LD expansion, NetMHCpan,
BLAST), each self-contained with configuration, scripts, notebooks, and logs.  
Outputs are stored centrally under `results/`.

-------

## Repository Layout

analyses/ # per-analysis modules (e.g. ld_expansion, netmhcpan, blast)  
results/ # centralized outputs (per module subfolders)  
data/refs/ # shared reference files (genomes, HLA lists, databases)  
src/virauto/ # reusable helper code  
tests/ # smoke/unit tests  
docs/ # notes, documentation  
analyses/ # per-analysis modules (e.g. ld_expansion, netmhcpan, blast)  
results/ # centralized outputs (per module subfolders)  
data/refs/ # shared reference files (genomes, HLA lists, databases)  
src/virauto/ # reusable helper code  
tests/ # smoke/unit tests  
docs/ # notes, documentation  

Each analysis module follows a standard structure: 

analyses/<module>/    
config/ # YAML/TSV configs, manifests    
scripts/ # Python/R/Bash/SLURM entrypoints    
notebooks/ # Jupyter / RMarkdown exploration     
logs/ # runtime logs    
slurm/ # HPC job scripts    
