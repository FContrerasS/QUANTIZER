New N-body code overview (architecture, input/output, general principles)
=========================================================================

**VERSION INFORMATION**: Felipe Contreras, 2022-10-01, version 1.0.


# Overall architecture of the N-body code #

## Files and directories ##

After downloading the new N-body code, one can see the following files in the
root directory contains:

- `bin/` contains the `Makefile` with which you can compile the code by typing
one of this two options with the bash in the `bin` folder: 
  -# `make clean; make; make run` 
  -# `make test`

- `headers/` contains all the headers files with a .h suffix.

- `source/` contains the C files for each modulce of the new N-body code module,
and 2 folders, the `force_methods`, and the `initial_conditions`.

- `source/force_methods/` contains the numerical methods to solve the Poisson
equation.

- `source/initial_conditions/` contains several possible particle initial
condition as input of the code.

- OTHERS files, this part is still under construction.

# The code flow #

## Key Concepts ##

In this secction we are goind to define some useful concepts definitions wich
are used through the whole manual.

- \anchor Key_Concepts_Capacity [a] **Capacity**:  Maximum amount of elements
that can be held in an array.

- \anchor Key_Concepts_Size [b] **Size**: Current number of elements in an
array.

- \anchor Key_Concepts_Smallest_Box [c] **Smallest Box**:  Logical array with
box geometry containing all the existing cells of the node, that is, cells with
status \f$ > -4 \f$, using the smallest possible space.

- \anchor Key_Concepts_Code_Space [d] **Code Space**: In a typical simulation,
the user chooses the coarsest level of refinement \c lmin (\f$ l_{min}\f$), the
maximum level \c lmax (\f$ l_{max}\f$), and the simulation box with any length
unit, and any coordinate system. But, the code always transforms the user box in
the coarsest level of refinement \c lmin to fit in a cube of side equal to \f$
2^{l_{min}} \f$ localized at the position *(0,0,0)* in a coordinate system, so
the cube can be described as the set of points which belong to \f$
[0,2^{l_{min}})\times [0,2^{l_{min}})\times [0,2^{l_{min}})\f$. The
3-Dimensional space which goes from \f$ (-\infty,\infty)\times (-\infty,
\infty)\times (-\infty,\infty)\f$ contains this cube of side \f$ 2^{l_{min}} \f$
localized at *(0,0,0)* is called the \f${\color{red} \mathbf{ Code\ Space\ of\
refinement\ l_{min} }}\f$. For any level of refinement *l*, the box is localized
in the same coordinate *(0,0,0)* but using a space equal to \f$ [0,2^l)^3\f$.

- \anchor Key_Concepts_Logical_Space [e] **Logical Space**: Because of the
enormous memory space required to store a full level of refinement *l* of size
\f$ (2^{l})^3 \f$, the concept of boxes is implemented in every node. Boxes are
logical entities that represent only a small piece of space of the *Code Space
of refinement level l*, where *l* is the level of refinement of the node. This
representation of the space is called the \f${\color{red} \mathbf{ Logical\
Space}}\f$. For every node at every level of refinement, this *Logical Space*
has dimensions of \f$ [0,A)\times[0,B)\times[0,C) \f$, where *A, B, C* are the
dimensions of the box of the node.

## Overview ##

### Vision ###

The New N-body code is a code designed to perform simulations of N-body without collisions in the area of astrophysics, particularly in cosmology. The focus of the New N-body code is to be run effectively in small-scale devices such as desktop computers and laptops, with a code design that is efficient, modular, accurate and user friendly. 

To run fast the New N-body on low-power machines, we decided to use C as the main language, and to build an optimum code design to include vector instructions (AVX), parallelization on CPUs, and mainly the use of GPUs. 

To achieve modularity, the code is implemented in serval modules which are almost independent between them. We also decided to use the scheme of \f$``{\color{red} \mathbf{ grid\ based}}"\f$ method to compute the acceleration of the particles through the Poisson equation. This choice is because we intend to implement in the future other types of gravity solvers for MOND, f(R), GR extension, etc, in an interchangeable way.
  
The accuracy of the code is reached by performing a finite difference method of 7-points for the 3D-Laplacian, 3 or 5 stencil points (user option) for the computation of the acceleration, and with a \e "Cloud In Cell" (CIC) scheme for the mass and acceleration transfer between grid points and particles. Moreover, we are using the Leapfrog method as the integrator of the time steps using the "kick-drift-kick" form which is separated in two steps. First is the "predictor step":

\f[ v^{n+1/2} = v^{n} + a^{n} \cdot {dt \over 2} ,\f]
\f[ x^{n+1} = x^{n} + v^{n+1/2} \cdot dt .\f]

And, after that, there is a "corrector step":

\f[ v^{n+1} = v^{n+1/2} + a^{n+1} \cdot {dt \over 2} .\f]

With this scheme is possible to move forward the time of the simulation for long periods with a low loss on the accuracy of the solution. 

It is really hard to accommodate a design of a code that has these three characteristics: efficiency, modularity, and accuracy at the same time. The efficiency and the accuracy go hand in hand, and in several parts, we decide to sacrifice a bit of performance in favor of getting better modularity. Therefore, incorporating the "User Friendly" characteristic into the code seems to be an almost impossible task. We know, that most astrophysics use Python as the main or one of the main languages in your day-to-day. So initially, we thought to implement the code in Python, but the efficiency diminished markedly. Moreover, the use of GPUs is very well designed to be implemented in C/C++ but not in Python. For this reason, we decided to do not to write the code in Python and sacrifice "User Friendly" in favor of efficiency. However, if you are astrophysics dedicated to working in programming mainly, and you require to modify the code, you probably will know basic programming languages such as C/C++ or Fortran, which will be enough for you to perform the modifications. But if you are only a user of the code, it will not be necessary for you to know about the C language because you only will need to modify a .dat file with the initial parameters. Moreover, we have planned to implement a GUI user interface in the future. 

To carry out our vision of the code, we created a new design based on a \e "tree structure" with \e "nodes". 

### Tree Structure, Nodes and AMR ###

We implemented a grid-based scheme with levels of refinement organized through a tree, where every level of refinement corresponds to a level of the tree. The figure \ref FIG--Tree_of_Refinements "FIG: Tree of Refinements" shows a representation of this. 

\anchor FIG--Tree_of_Refinements \image latex FIG--Tree_of_Refinements.pdf "Tree of Refinements" 

In every level of the tree, there are several \f${\color{red} \mathbf{Nodes}}\f$ which are the minimum tree structure. The Parent node is connected with its child nodes, and the child nodes are connected with their parent node, but there is no connection between sibling nodes. The Nodes represent refinement zones separated physically from the other refinement zones at the same level. Every node contains cells, particles, grid properties (acceleration, density, and potential), boundary flags (boolean parameters), type of grid point (interior, boundary, simulation boundary), box logical representation, and many other parameters (for a complete review of the content of the nodes see \ref node "struct node").

 This tree is adapted in every time step in the well-known scheme "Adaptive Mesh Refinement" (AMR). In every adaptation, some new nodes are created, others removed and other adapted to accommodate to the new particle configuration. 

The moment when the tree is adapted and the whole flow of the simulation structure can be seen in the following chapter. 

### Code Flow ###

In this section, we are going to show an overview of the operation of the New    N-body code. 

This overview can be observed in the flow chart \ref FIG--Global_Flow_Chart "FIG: Global Flow Chart". Here, the code flow, can be separated into two big parts separated by the central question \e "Is dt and N° steps lower than the Threshold?". <b>(a)</b>  The initial part is always executed and this one is in charge of initializing the input, global and user parameters, the initial tree of refinement levels, the computation of the cero time step, i.e, when \f$ dt=0 \f$, and the exportation of the first snapshot with the data of the initial particle configuration and observables. <b>(b)</b> The second part is repeated in a loop that ends when the time reached by the simulation, or the number of iterations exceeds the threshold previously defined by the user. In this loop, the process to compute the acceleration of the particles is repeated, but using the new particle configuration and a new tree configuration of the refinement levels. This last part is the heart of the code, which allows us to perform the well known \e "Adaptive Mesh Refinement" scheme (AMR), and the details of its operation can be seen in the module \ref tree_adaptation.c. 

\anchor FIG--Global_Flow_Chart \image latex FIG--Global_Flow_Chart.pdf "Global Flow Chart" 

### How the Nodes work? ###


