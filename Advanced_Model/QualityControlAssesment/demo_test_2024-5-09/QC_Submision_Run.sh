#!/bin/bash

# Define the parent directory
parent_dir="QC"

# Navigate to the parent directory
#cd $parent_dir

## Find the most recent directory
#latest_dir=$(ls -td -- */ | head -n 1)

# Navigate to the most recent directory
#cd $latest_dir

#cp ../../"NW_functions.py" .

# Get the original directory (now the most recent directory)
#original_dir=$(pwd)
echo "Test 1"
for file_name in $(ls); do
    sbatch <<EOF
#!/bin/bash
#SBATCH --time=10:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=1024M   # memory per CPU core
echo "TEST2"
# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE
echo "TEST3"
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
python QC_Submision_Run.py $file_name #> ${file_name%.pickle}.out
echo "TEST4"
# Copy the output file back to the original directory
#cp ${file_name%.pickle}.out $original_dir
EOF
done
