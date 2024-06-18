# Quality Control Assessment Tool

The quality control assessment is a developer's tool used to compare iterations of the code. This is done by creating a test folder whose data is generated using a new iteration of the functions.py. Once the data has been generated it is compared to a base set of data to check its performance. In this description, we will explore the structure of this directory, explain the use of each file, and describe the process for conducting an assessment.

## Directory Structure 
* ### QCbase
  - QCbase is a directory for the collection of data generated with an older, more established version of the code. Before it can be used the data must be generated. this can be done by first calling QC_Control.py by typing *python QC_Control.py* in the command line. This generates a batch of files for every combination of the following:

      * high concentration (10 mol% UCl3)
      * low concentration (0.1 wt% UCl3)
      * fast scan rate (1 V/s)
      * low scan rate (10 mV/s)
        *  Diffusion only (no migration, ideal activites, no ohmic losses) ###0000
        *  Test solid activity only                                        ###0001
        *  Test ion activity only                                          ###0010
        *  Test migration only                                             ###0100
        *  Test migration with complexation                                ###0100b
        *  Test ohmic losses                                               ###1000
        *  Test ohmic losses with 85% IR compensation                      ###1000b
        *  Test solid activity, ion activity, migration with complexation, and ohmic losses with IR compensation.                                         ###1111

   Once those are created, the runs may be started by calling *bash QC_base_check.sh*. This will submit and run each pickle data file to the supercomputer. Once complete, the run will:
          1) update the pickle file with a full data set
          2) print a *.out file for each data set with information on the speed at which data points were collected
          3) create a png file for the simulation.

    After the QC base files have all been run successfully we can move forward with testing.
* ### test_setup
  - This directory contains ll the raw files needed to create a test file these include:
     * NW_functions.py : A file containing the new functions.py that is to be compared during the quality assessment. Before each test make sure that this is the most recent updated file.
     * OutProcessing.py : A python script that will compare the completed test runs against those generated previously in the QCbase directory.
     * QC_Submision_Run.py : The python submission script for a run.
     * QC_Submision_Run.sh : The bash submission script for a run.
     * QC_Test_Construction.py : The primary executable which will create new pickle data files for the diffrent tests and then submit said files to the supercomputer using the other scripts in the directory.
  - 
 
-     
* ### demo_test_2024-5-09
  - This directory contains a selection of files that can be found in a test directory that was succesfuly run. It should be noted that only a few files are included for demonstration purposes and that no data files are included as they are too large for demonstration.
* ### test_creation.sh
  - A Bash script that when called will create a new directory time-stamped and ready to run.

## Instructions for Conducting an Assessment 
Note: These instructions have been prepared as if these tests were being conducted on the BYU super computer.

1) Ensure that a complete QCbase folder has been created using the correct (validated) version of functions.py.
2) Copy the new version of the functions.py script into the test_setup directory using the name NW_functions.py.
3) Call *bash test_creation.sh* , this will generate a new test directory with the files we need. It will also initiate the test by calling the job submission script.
4) Check on the progress of each test by calling *squeue -u enteryourusername*, using your own username.
5) Once all tests are complete, navigate to the new directory and call *python OutProcessing.py*, which will create various data files, including comparative performance plots for each test and a QA_Output.txt file which summarizes the Quality control Assessment Results. 
