import os
import shutil
import subprocess
from datetime import datetime

# Specify the source directory
src_dir = "./batch_setup"

# Create a new directory with the date and time
new_dir = "./batch_" + datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
os.makedirs(new_dir, exist_ok=True)

# Copy the contents of the source directory to the new directory
for item in os.listdir(src_dir):
    s = os.path.join(src_dir, item)
    d = os.path.join(new_dir, item)
    if os.path.isdir(s):
        shutil.copytree(s, d, False, None)
    else:
        shutil.copy2(s, d)

print(f"Contents of {src_dir} have been copied to {new_dir}")

os.chdir(new_dir)

# Convert all files to Unix format and change permissions
#for file in os.listdir("."):
#    subprocess.run(["dos2unix", file])
#    os.chmod(file, 0o777)

#print("All files have been converted")

print("Please remember to fill out the create_simulation.py, following the format demostrateted in example.py before executing the batch")

# Run the Python script


print("Batch setup  complete")
