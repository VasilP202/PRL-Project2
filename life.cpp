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
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (argc != 3) {
        if (rank == 0) {
            cerr << "Usage: " << argv[0] << " <input_file> <number_of_steps>" <<
                endl;
        }
        MPI_Finalize();
        return 1;
    } else {
        if (rank == 0) {
            cout << "Number of processes: " << num_procs << endl;
            cout << "Input file: " << argv[1] << endl;
            cout << "Number of steps: " << argv[2] << endl;
        }
    }
    string filename = argv[1];
    int num_steps = atoi(argv[2]);

    vector < vector < int > > grid(num_procs, vector < int > (num_procs, 0));
    initializeGridFromFile(grid, filename);

    if (rank == 0) {
        printGrid(grid);
    }

    /* for (int step = 0; step < num_steps; ++step) {
        updateGrid(grid);

        MPI_Barrier(MPI_COMM_WORLD);
    } */

    MPI_Finalize();
    return 0;
}