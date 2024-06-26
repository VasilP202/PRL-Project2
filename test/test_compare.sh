#!/bin/bash

# Name of the file containing the definition of the game grid for the Life game
life_grid_file=$1;

# Number of iterations to be performed
num_iterations=$2;

# Number of processors to be used = number of rows in the game grid
num_processors=$(sed -n '$=' "$life_grid_file")

# Compile the source code
mpic++ -o life life.cpp

# Run the program
mpirun --oversubscribe -np $num_processors life $life_grid_file $num_iterations > c.out

python3 script.py $life_grid_file $num_iterations > py.out

# Compare the outputs
diff c.out py.out

# Cleanup
rm -f life c.out py.out