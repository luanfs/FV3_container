# FV3_container

----------------------------------------------------------------------------------------------------------------------
# Container
Steps to run the container:
 1) Pull the docker image: $docker pull gfdlfv3/shield-dev
 2) cd to sh and run ./init_container.sh
 3) It starts a command line like: root@dbdb48866381:
 
# Compilation
Steps to compile the code in the container:
 1) Change the directory to build: root@dbdb48866381:/# cd SHiELD_build/Build/
 2) Compile the code  (it should take a while): root@dbdb48866381:/SHiELD_build/Build# ./COMPILE solo sw debug gnu 64bit cleanall
 3) Exit the container and save its image: $docker commit dbdb48866381 sw

# Running a shallow-water or advection simulation
 1) Run the container saved image using the script ./run_sw_container.sh
 2) cd to SHiELD_SRC/test/
 3) Run ./run_adv.sh and adjust this file for you preferences 
 4) Plots scripts are available in /plot (outputs are saved in /graphs).
 
You should recompile the code using compile.sh.

----------------------------------------------------------------------------------------------------------------------
References:
- Mouallem, J., Harris, L., & Chen, X. (2023). Implementation of the novel Duo-Grid in GFDL's FV3 dynamical core. Journal of Advances in Modeling Earth Systems, 15, e2023MS003712. https://doi.org/10.1029/2023MS003712 
- Mouallem, J. (2023). Implementation of the novel duo-grid in GFDL’s FV3 dynamical core - code and simulations files [Software] [Dataset]. Zenodo. https://doi.org/10.5281/zenodo.8327578
