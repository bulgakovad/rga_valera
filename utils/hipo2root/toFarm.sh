#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name=H2R_RGA_Valera
#SBATCH --mail-type=ALL
#SBATCH --mail-user=aleksandr.bulgakov1999@gmail.com
#SBATCH --output=/farm_out/%u/%x-%A_%a-%j-%N.out
#SBATCH --error=/farm_out/%u/%x-%A_%a-%j-%N.err
#SBATCH --partition=production
#SBATCH --account=clas12
#SBATCH --mem-per-cpu=5000
#SBATCH --time=24:00:00



# Navigate to the working directory (if needed)
cd /w/hallb-scshelf2102/clas12/bulgakov/projects/rga_valera/utils/hipo2root

module use /scigroup/cvmfs/hallb/clas12/sw/modulefiles
module load clas12root


# Run the executable with the job array index (if needed)
srun clas12root -q -b ana12GeVShortFCQA.C  --in=pass1_rga_inclusive_runs_valera.dat >> ../../logs/pass1_rga_inclusive_runs_valera.log