#!/bin/bash

# Get the original directory
original_dir=$(pwd)

for file_name in $(ls data/); do
    sbatch <<EOF
#!/bin/bash
#SBATCH --time=05:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=1024M   # memory per CPU core

# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
python submision_script.py $file_name #> ${file_name%.pickle}.out

# Copy the output file back to the original directory
#cp ${file_name%.pickle}.out $original_dir
EOF
done
