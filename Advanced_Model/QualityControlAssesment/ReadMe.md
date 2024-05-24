# Quality Control Assessment Tool

The quality control assessment is a developer's tool used to compare iterations of the code. This is done by creating a test folder whose data is generated using a new iteration of the functions.py. Once the data has been generated it is compared to a base set of data to check its performance. In this description, we will explore the structure of this directory, explain the use of each file, and describe the process for conducting an assessment.

## Directory Structure 
* QCbase
  - QCbase is a directory for the collection of data generated with an older, more established version of the code. Before it can be used the data must be generated. this can be done by first calling QC_Control.py by typing *python QC_Control.py* in the command line. This generates a batch of files for each of the following:

      * high concentration (10 mol% UCl3)
>  low concentration (0.1 wt% UCl3)
>  fast scan rate (1 V/s)
>  low scan rate (10 mV/s)
 

-- Diffusion only (no migration, ideal activites, no ohmic losses)  ###0000
#- Test solid activity only                                         ###0001
#- Test ion activity only                                           ###0010
#- Test migration only                                              ###0100
#- Test migration with complexation                                 ###0100b
#- Test ohmic losses                                                ###1000
#- Test ohmic losses with 85% IR compensation                       ###1000b
#- Test solid activity, ion activity, migration with complexation, and ohmic losses with IR compensation. ###1111
 
* test_setup
  - explanation  
* demo_test_2024-5-09
  - explanation 

## Instructions for Conducting an Assessment 
