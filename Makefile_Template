# Build a executable called vaca


# the compiler: gcc for C program, define as g++ for C++
CC = gcc

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall

#Optimization Flags
OPTIMIZATION  =  -O3

GPROF = -pg

# Headers directory
HEADERS = ../headers

# Initial Conditions 
SOURCE = ../source/global_variables.c ../source/input.c ../source/main.c ../source/output_snapshots.c ../source/tree_construction.c ../source/initialize_node.c ../source/potential_head_node.c	../source/initialization.c ../source/poisson_error.c	../source/grid_density.c ../source/potential.c ../source/grid_acceleration.c ../source/particle_acceleration.c ../source/timestep_1.c ../source/timestep_2.c ../source/particle_updating_A.c ../source/particle_updating_B.c ../source/terminal_print.c ../source/observables.c ../source/output_main_parameters.c ../source/reset.c ../source/initial_potential_and_acceleration_head.c ../source/space_check.c

OBJ = global_variables.o input.o main.o output_snapshots.o tree_construction.o initialize_node.o potential_head_node.o initialization.o poisson_error.o grid_density.o potential.o grid_acceleration.o particle_acceleration.o timestep_1.o timestep_2.o particle_updating_A.o particle_updating_B.o terminal_print.o observables.o output_main_parameters.o reset.o initial_potential_and_acceleration_head.o space_check.o
 
 
.PHONY: run

# the build target executable:

# Main.cpp Initialization.cpp Jacobi.cpp GS.cpp Fine2coarse.cpp Coarse2fine.cpp Vgrid.cpp Fgrid.cpp Wgrid.cpp Error.cpp

# the build target executable:
EXECUTABLE = nbody.exe
	
all:
	make compile -s
	make link -s

compile:
	$(CC) $(CFLAGS) -I$(HEADERS) $(OPTIMIZATION) -c $(SOURCE) 

link:
	$(CC)  $(OBJ) -o $(EXECUTABLE) -lm
	
gprof: *.cpp
	$(CC) $(CFLAGS) -I $(HEADERS) $(INICOND) $(OPTIMIZATION) $(GPROF) -o $(EXECUTABLE).gprof $(SOURCE) -lm
	
clean:
	$(RM) -f $(EXECUTABLE) *.o
	
run:
	./$(EXECUTABLE)	
	
test:
	make clean -s
	make respaldo -s
	make run -s
	make clean -s
	
respaldo: 
	$(CC) $(CFLAGS) -I $(HEADERS) $(OPTIMIZATION) -o $(EXECUTABLE) $(SOURCE) -lm	
	
