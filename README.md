# HEMStoEC
This repository contains code and data relevant to HEMStoEC: Home Energy Management Systems to Energy Communities DataSet - https://zenodo.org/record/8096648

All code for the generation of the dataset was written in Matlab R2022. 

Daily information is received by the data acquisition system in a zipped file, which should be placed in the same directory (denoted as root directory) of the function files. 

A sample can be found in 2023_06_11_00_00_00.zip. The README and the VARS files provide information about the format of the files enclosed in the zip file. 

Matlab data is extracted from the unzipped file using the Matlab function extract_quadro_10.m. The command extract_quadro_10('2023_06_11_00_00_00') creates a Matlab data file 2023_06_11_00_00_00.mat inside the 2023_06_11_00_00_00 directory.  

Gaps are identified and data is interpolated using the function Validate_Quadro_4.m. A data file 2023_06_11_00_00_00_cor.mat is created, again inside the 2023_06_11_00_00_00 directory, upon the command Validate_Quadro_4('2023_06_11_00_00_00','2023_06_11_23_23_59').

Data with a common time basis is achieved using the Matlab function convert_quadro_10_cor.m. Using the command convert_quadro_10_cor('2023_06_11_00_00_00','2023_06_11_23_23_59','', minutes(15),hours(1)), the data file 2023_06_11_00_00_00 to 2023_06_11_23_23_59 excl  pst 15 min est 1 hr_cor.mat is created, this time in the root directory. 

A matlab file, Factor.mat, needs to be placed in the root directory.