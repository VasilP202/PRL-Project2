#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

using namespace std;

const int GRID_SIZE = 8;

void initializeGridFromFile(vector < vector < int > > & grid, const string & filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    string line;
    int row = 0;
    while (getline(file, line)) {
        if (row >= GRID_SIZE) {
            break;
        }
        for (int col = 0; col < GRID_SIZE; ++col) {
            if (col >= line.size()) {
                break;
            }
            if (line[col] == '1') {
                grid[row][col] = 1;
            } else {
                grid[row][col] = 0;
            }
        }
        ++row;
    }
}

void updateGrid(vector < vector < int > > & grid) {
}

void printGrid(const vector < vector < int > > & grid) {
    for (int row = 0; row < grid.size(); ++row) {
        for (int col = 0; col < grid[row].size(); ++col) {
            cout << grid[row][col] << " ";
        }
        cout << endl;
    }
}

int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    string filename = argv[1];
    int num_steps = atoi(argv[2]);
    const int grid_size = num_procs;
    
    vector < vector < int > > grid(num_procs, vector < int > (num_procs, 0));

    if (argc != 3) {
        if (rank == 0) {
            cerr << "Usage: " << argv[0] << " <input_file> <number_of_steps>" << endl;
        }
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
        initializeGridFromFile(grid, filename);
        printGrid(grid);
    }

    int row[grid_size];
    int neighbors[8]; 

    MPI_Scatter(&grid[0][0], grid_size, MPI_INT, row, grid_size, MPI_INT, 0, MPI_COMM_WORLD);

    int left_rank = (rank - 1 + num_procs) % num_procs;
    int right_rank = (rank + 1) % num_procs;

    MPI_Sendrecv(&row, grid_size, MPI_INT, left_rank, 0,
                 &neighbors[0], grid_size, MPI_INT, right_rank, 0,
                 MPI_COMM_WORLD, &status);

    MPI_Sendrecv(&row, grid_size, MPI_INT, right_rank, 0,
                 &neighbors[4], grid_size, MPI_INT, left_rank, 0,
                 MPI_COMM_WORLD, &status);


    // Print the neighbors
    cout << "Rank " << rank << " has neighbors: ";
    for (int i = 0; i < 8; ++i) {
        cout << neighbors[i] << " ";
    }

    /* for (int step = 0; step < num_steps; ++step) {
        updateGrid(grid);
        MPI_Barrier(MPI_COMM_WORLD);
    } */

    MPI_Finalize();
    return 0;
}