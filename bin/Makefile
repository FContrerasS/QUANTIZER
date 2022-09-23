# Build a executable called vaca


# the compiler: gcc for C program, define as g++ for C++
CC := gcc

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS := -g -Wall

#Optimization Flags
OPTIMIZATION := -O3

GPROF := -pg

# Headers directory
HEADERS = ../headers

#Source Directories
FMDir := ../source/force_methods
GDir := ../source

#Sources
FMSources := conjugate_gradient.c conjugate_gradient_multigrid.c jacobi.c gauss_saidel.c multigrid.c
FMSources += multigrid2.c
##
GSources := global_variables.c input.c main.c output_snapshots.c tree_construction.c initialize_node.c 
GSources += potential_head_node.c initialization.c poisson_error.c	grid_density.c potential.c 
GSources += grid_acceleration.c particle_acceleration.c timestep.c particle_updating_A.c 
GSources += particle_updating_B.c terminal_print.c observables.c output_main_parameters.c reset.c 
GSources += initial_potential_and_acceleration_head.c space_check.c tree_adaptation.c new_node.c 
GSources += add_node_to_stack.c initialize_cell_struct.c ptcl_idx_to_box_idx.c garbage_collector.c 
GSources += potential_branches.c 

FORCE_METHODS = $(patsubst %,$(FMDir)/%,$(FMSources))
GENERAL = $(patsubst %,$(GDir)/%,$(GSources))

#All sources
SOURCES = $(FORCE_METHODS) $(GENERAL)


#Objects
#ODIR = ob3

OBJ := global_variables.o input.o main.o output_snapshots.o tree_construction.o initialize_node.o 
OBJ += potential_head_node.o initialization.o poisson_error.o grid_density.o potential.o grid_acceleration.o 
OBJ += particle_acceleration.o timestep.o particle_updating_A.o particle_updating_B.o 
OBJ += terminal_print.o observables.o output_main_parameters.o reset.o initial_potential_and_acceleration_head.o 
OBJ += space_check.o tree_adaptation.o new_node.o add_node_to_stack.o initialize_cell_struct.o 
OBJ += ptcl_idx_to_box_idx.o garbage_collector.o potential_branches.o conjugate_gradient.o 
OBJ += conjugate_gradient_multigrid.o jacobi.o gauss_saidel.o multigrid.o multigrid2.o
OBJ += 

#OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

.PHONY: run

# the build target executable:

# Main.cpp Initialization.cpp Jacobi.cpp GS.cpp Fine2coarse.cpp Coarse2fine.cpp Vgrid.cpp Fgrid.cpp Wgrid.cpp Error.cpp

# the build target executable:
EXECUTABLE := nbody.exe
	
all:
	make compile -s
	make link -s

compile:
	$(CC) $(CFLAGS) -I$(HEADERS) $(OPTIMIZATION) -c $(SOURCES)

link:
	$(CC)  $(OBJ) -o $(EXECUTABLE) -lm
	
gprof: *.cpp
	$(CC) $(CFLAGS) -I $(HEADERS) $(INICOND) $(OPTIMIZATION) $(GPROF) -o $(EXECUTABLE).gprof $(SOURCES) -lm
	
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
	$(CC) $(CFLAGS) -I $(HEADERS) $(OPTIMIZATION) -o $(EXECUTABLE) $(SOURCES) -lm	
	
