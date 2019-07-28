# Conway's Game of Life using MPI

## Description
- **serial.cpp** is a serial implementation without use of MPI
- **mpi.cpp** is a MPI based implementation without use of ghost cells.
- **mpi_ghost.cpp** is MPI based implementation using ghost cells.

Sample input/output files are provided. *#* - live cell. . *-* dead cell. In case of MPI version use **MPIC++** compiler. You can use **--oversubscribe** if testing on local machine for big processor numbers.

The main grid was partionened into smaller 2D subgrids.
