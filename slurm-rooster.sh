#!/bin/bash
#SBATCH --job-name=nacie-up    # Job Name
#SBATCH --partition=daily    # Using 'daily' will grant higher priority than 'general'
#SBATCH --time=23:00:00    # Time needed for running the job. Must match with 'partition' limits.
#SBATCH --nodes=1            # Number of nodes
#SBATCH --ntasks=1          # Number of tasks
#SBATCH --cpus-per-task=1    # Double if hyperthreading enabled
#SBATCH --ntasks-per-core=1  # Run one task per core
#SBATCH --hint=nomultithread # Disable Hyperthreading
#SBATCH --error=slurm-%j.err # Define a file for standard error messages
##SBATCH --exclusive         # Uncomment if you want exclusive usage of the nodes

# module use unstable
# module load ANSYS/2020R1-1
module load psi-python39/2021.11

# [Optional:BEGIN] Specify your license server if this is not 'lic-ansys.psi.ch'
# LICENSE_SERVER=<your_license_server>
# export ANSYSLMD_LICENSE_FILE=1055@$LICENSE_SERVER
# export ANSYSLI_SERVERS=2325@$LICENSE_SERVER
# [Optional:END]


python3 A_rooster.py > log


