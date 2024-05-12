# Advanced Simultation Tool for Cyclic Voltammetry of Electrodeposition Reactions in High-Temperature Molten Salts

This Code was First Uploaded on 5/10/2024 by Braden Clayton. It was designed and developed by Michael Stoddard as part of his graduate studies while at Brigham Young University.

## Acknowledgement and Attribution
This code is freely available, however it is **required** to cite the following when using this code:
Stoddard, M. (2024). *The Advancement of Experimental and Computation Tools for the Study of Molten Salt Chemistry to Facilitate the Extraction  of Strategic Elements in Nuclear Applications* [Ph.D. Dissertation]. Brigham Young University.

## Purpose 
Intended to be an advanced simulation tool applicable to a variety of electrochemical systems, this code is still under development. The code is being developed specifically for voltammetry of electrodeposition in high-temperature molten salts (e.g., LiCl-KCl, FLiBe). It's intended to relax the assumptions involved in the Berzins-Delahay model, which is commonly used in analyzing electrodeposition peaks in high-temperature molten salts. However, the Berzins-Delahay model assumes:
  -  unit activity of the deposits (i.e., pure),
  - ideal behavior of analytes ions (i.e., activity coefficient of 1), and
  - mass transport by diffusion only
In it current form, this advanced model, has include the following missing in the Berzins-Delahay model:
  - non-unit activity of initial deposits on a forgien substrate (e.g., U deposition onto W)
  - activity coeffient as a function of composition based on colbalt(II) chloride data
  - migration of charged particles due to an electric field
  - uncompensated resistance

Users are encourage to use this model and provide feedback to accelerate its development and validation

## Structure
The code currently (5/10/2024) consists of two tools. The first is called batch creation and the second quality control assessment. 
Batch creation can be considered a tool for general applications of the code.
  Files and directories pertaining to this tool include:
    > batch_creation.sh
        A bash script that creates a new directory with the appropriate files.
    > batch_creation.py 
        A Python script that creates a new directory with the appropriate files.
    > batch_setup/ 
        A directory containing the basic files required for each batch. See the directory's readme for further detail.
    > demo_batch/
        A directory containing a selection of files from a successful batch. Raw data could not be included due to file size limitations.
    > batch_2024-05-10_14-17-51
        A Sample directory was created to demonstrate the code which was never run.
  
    The batch function is initiated by calling either of the creation scripts. This can be don by executing the file in a python editor such as pycharm,
    or by using the comandline while in this directory and calling "python batch_creation.py" or "bash batch_creation.sh".
    Once this is done, a new directory will be made with the appropriate base files.

    The next step will be to create the data files that will be executed as part of the simulation. Details on this process are included in the readme file 
    included in the batch_setup directory which will also be copied to the new directory. 
    
Quality control assessment is a developer's tool used to compare iterations of the code.
  A more detailned explanation of this code can be found in the directory's readme file.

