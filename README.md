# New_Code
A new code for N-body/Hydrodynamics simulations in astrophysics

This new code use a novel design for the adaptative mesh refinement method. 
The code is focused to small scale devices as desk computers or laptos.  
Right now the code is functional with a novel design in the AMR part, but it only can work in a secuencial way and with a sharing time-step.

The future updating will include the paralelization in both, CPU and GPU. Moreover, we will include the variable time-step for the different level of refinements.

###########################################################################################################################################################

The code is written in c, and compiled and linked with gcc. 
The make file is inside of the bin folder. 

Console Input: Does not required
Static Input: The code use constant parameters which are defined in the "global_variables.c" source file and in "commen.h" header file. 
Runtime Input: Currently there are 2 options. First a random imput given in the "input.c" source file in the function called                                              "input_particle_initialization". The another input type is using a .txt file given the position, x,y,z and the velocities vx,vy and vz per                  each particle. This model is in the "input.c" source file in the function called "input_plummer_model".  
               
Output: The code relases binary files with the information about the mass, position, velocity, energies, and time of the simulation of the whole bunch of           particles. the source file which contain that is called "output_snapshots.c". Moreover, the code generates a .dat file called "Parameters.dat"             which contains the main parameter values used in the simulation. Finally, there is another information by console as an approximation of the memory         used and the time-clock of each part of the simulation code.               
