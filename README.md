# eORCA1-BGC-PDAF
Scripts to set up NEMO runs on the eORCA1 domain for data assimilation runs using PDAF. Configured to use either MEDUSA or ERSEM biogeochemical models

# Code compiliation
To compile code, cd into code directory. Update the work directory in 0_set_environment.sh, then tun scripts 1-5 to clone repositories and compile code for xios, fabm, pdaf and nemo. 

# Link Larger run files
cd to INPUTS directory and execute link_inputs.sh. This will create directories with links to all the files needed to run. 

# To run a simulation
cd to cycle and update parameters in set_environment.sh. Then submit a job using the command:
sbatch --export=year=XXXX cycle_year.slurm
where XXXX is the year you want to start, e.g. 2000. The run will submit 12 months on a single job and will cycle each year up until the end year specified in set_environment.sh. 
