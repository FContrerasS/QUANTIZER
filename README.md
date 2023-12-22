# QUANTIZER
QUANTIZER is a collisionless N-body code for simulations in astrophysics

#############################################################################################

Author: Felipe Contreras (contreras.sep.felipe@gmail.com).

Date: 2022-05-26.

Code version: beta-v1.0.0

#############################################################################################

This new code uses a novel design for the adaptative mesh refinement method. 

The code is focused to small scale devices as desk computers or laptos. 

Right now the code is functional with a novel design in the AMR part, but it only can work in a secuencial way and with a sharing time-step.

The future updating will include the paralelization in both, CPU and GPU. Moreover, we will include the variable time-step for the different level of refinements.

#############################################################################################

QUANTIZER is written in c, and compiled and linked with gcc. 

The make file is inside of the bin folder. 

Console Input: Does not required

Static Input: The code use constant parameters which are defined in the "global_variables.c" source file and in "commen.h" header file. 

Runtime Input: In the input file put your .txt file containing the position in x,y,z, and velocity x,y,z and mass of every paticle per row. In the source folder you will find the file  "input.c". There you can add the name of your file to run it. Moreover, it is necessary to fill some variables as the number of particles, the boxsize, between others. They appears in the file "global_variables.c"
               
Output: The code relases binary files with the information about the mass, position, velocity, energies, and time of the simulation of the whole bunch of particles. the source file which contain that is called "output_snapshots.c". Moreover, the code generates a .dat file called "Parameters.dat" which contains the main parameter values used in the simulation. The output of QUANTIZER, is localized in the "output" folder.  Finally, there is another information by console as an approximation of the memory used and the time-clock of each part of the simulation code.

An example appears in the folder "examples". To visualize it use the python file "example_plot.py"

#############################################################################################

To run the code use the makefile localized in the "bin" folder. 

Compilation:

> make

run the code

> make run








               
