from functions import *
import subprocess

#### INSERT CREATION SCRIPT HERE ####
  ## EXAMPLE CAN BE FOUND IN example.py ##


def run_python_script():
    try:
        # The command you want to execute
        command = "bash submision_script.sh"
        
        # Use subprocess to execute the command
        process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        
        # Get the output
        output, error = process.communicate()
        
        # Print the output
        print('Output: ', output.decode())
        
        # Print the error if there is any
        if error:
            print('Error: ', error.decode())
            
    except Exception as e:
        print("An error occurred: ", str(e))

# Call the function
run_python_script()
