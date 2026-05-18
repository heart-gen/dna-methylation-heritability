#!/bin/bash
#SBATCH --account=b1042
#SBATCH --partition=genomics
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=12G
#SBATCH --mail-type=FAIL
#SBATCH --job-name=regctx_plot
#SBATCH --mail-user=kynon.benjamin@northwestern.edu
#SBATCH --output=logs/regctx_plot.%j.log

log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

log_message "**** Job starts ****"

echo "**** QUEST info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"

module purge
module list

mkdir -p ../_m/logs
source /projects/p32505/opt/miniforge3/etc/profile.d/conda.sh
conda activate rnaseq

log_message "Plotting regulatory context"

RUNS=(
  "SHARED shared"
  "AA AA_only"
  "EA EA_only"
)

for RUN in "${RUNS[@]}"; do
  read -r POP VMR_SET <<< "$RUN"
  Rscript "../_h/05.plot_regulatory_context.R" all_individuals "$POP" "$VMR_SET"
done

conda deactivate
log_message "**** Job ends ****"
