#!/bin/bash

# Specify the source directory
src_dir="./test_setup"

# Create a new directory with the date and time
new_dir="./test_$(date +%Y-%m-%d_%H:%M:%S)"
mkdir -p "$new_dir"

# Copy the contents of the source directory to the new directory
cp -r -v "$src_dir/"* "$new_dir/"

echo "Contents of $src_dir have been copied to $new_dir"

cd "$new_dir/"

dos2unix *

chmod 777 *

echo "All files have been converted"

echo "Please remember to copy in the new functions module as NW_functions.py before executing the test, otherwise the default NW_funtions.py fount in test_setup will be used "

python "QC_Test_Construction.py" 

echo "Test construction has begun, to verify progress use the command squeue -u enteryourusername"
echo "Test setup and execution complete"``

