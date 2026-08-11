For any of the directories:

#cd $dir;

#make -f Makefile.tsp veryclean deps tsp

or 

#make -f Makefile.ising veryclean deps ising

then, the executable is in 

tsp/tsp_sa

or

tsp/tsp_sa.mpi

or

ising/ising_sa

or 

ising/ising_sa.mpi

With the mpi suffix if it is a parallel run, no mpi suffix for serial runs.

The hierarchical tsp_sa has the capacity to do runs where the initial states 
of all processors are specified by the user. 

The serial_ising prints state files and can be restarted from them,
however at the moment it appears that it does not properly recover the
random number generation and does not properly deal with old log files.

I am working on fixing that, but for the time being I would be cautious
about using such state files.

As of now, it is required that both a .instance and .params file
exist.

In the examples, the parallel files should work for 20 processors with
either hcdr or cdr.

For [h]cdr runs, you must run with

mpirun -np 20 [dir]/cdr_ising/ising/ising_sa.mpi parallel_tsp_example

where [dir] is where you keep this repository.
The program will automatically locate the instance and param files, and produce
parallel_tsp_example.log, parallel_tsp_example.output, and if you fix state
files it will also produce a set of parallel_tsp_example_XX.state
for each of the processors used.
