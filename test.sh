#!/bin/bash

# Name of the file containing the definition of the game grid for the Life game
life_grid_file=$1;

# Number of iterations to be performed
num_iterations=$2;

# Set the number of processors to be used
# Number of processors to be used = number of rows in the grid
num_processors=$(sed -n '$=' "$life_grid_file")

# Compile the source code
mpic++ -o life life.cpp

# Run the program with the given parameters
mpirun --oversubscribe -np $num_processors life $life_grid_file $num_iterations

# Cleanup
rm -f life