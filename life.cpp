#include <cstdlib>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <vector>

using namespace std;

const int GRID_SIZE = 8; // TODO remove
const int NUM_NEIGHBORS = 8;
const int CELL_DEAD = 0;
const int CELL_ALIVE = 1;


void initializeGridFromFile(vector<vector<int>> & grid, const string & filename, int grid_size) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    string line;
    int row = 0;
    while (getline(file, line)) {
        if (row >= grid_size) {
            break;
        }
        for (int col = 0; col < grid_size; ++col) {
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

void printGrid(const vector < vector < int > > & grid) {
    for (int row = 0; row < grid.size(); ++row) {
        for (int col = 0; col < grid[row].size(); ++col) {
            cout << grid[row][col] << " ";
        }
        cout << endl;
    }
}

int updateCellValue(int cell_value, vector<int> & neighbors) {
    int live_neighbors = 0;
    int dead_neighbors = 0;

    for (int i = 0; i < neighbors.size(); ++i) {
        if (neighbors[i] == CELL_ALIVE)
            ++live_neighbors;
        else // CELL_DEAD
            ++dead_neighbors;

        // Cell dies due to overpopulation if > 3 neighbors are alive
        if (cell_value == CELL_ALIVE && live_neighbors > 3)
            return CELL_DEAD;         
    }

    // Cell dies due to underpopulation if < 2 neighbors are alive
    if (cell_value == CELL_ALIVE && live_neighbors < 2)
        return CELL_DEAD;
    
    // Cell becomes alive if exactly 3 neighbors are alive
    if (cell_value == CELL_DEAD && live_neighbors == 3)
        return CELL_ALIVE;

    // If no conditions are met, cell value remains the same    
    return cell_value;
}

int main(int argc, char ** argv) {
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    const int grid_size = num_procs;
    int num_steps = atoi(argv[2]);

    vector<vector<int>> grid(num_procs, vector<int> (num_procs, 0));
    vector<int> flat_grid(grid_size * grid_size);

    if (argc != 3) {
        if (rank == 0) {
            cerr << "Usage: " << argv[0] << " <input_file> <number_of_steps>" << endl;
        }
        MPI_Finalize();
        return 1;
    }

    if (rank == 0) {
        string filename = argv[1];
        initializeGridFromFile(grid, filename, grid_size);
        for (int i = 0; i < grid_size; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                flat_grid[i * grid_size + j] = grid[i][j];
            }
        }
    }

    // Vector to store the row of the grid
    vector<int> row(grid_size);

    // Distribute the grid to all processes
    // Each processor gets a row of the grid based on its rank
    // Thus a single processor is going to handle a single row of a grid
    MPI_Scatter(flat_grid.data(), grid_size, MPI_INT, row.data(), grid_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Get the rank of the neighbor processors - Wrap-around implementation
    int neighbor_up_rank = (rank - 1 + num_procs) % num_procs;
    int neighbor_down_rank = (rank + 1) % num_procs;


    for (int step = 0; step < num_steps; ++step) {
        // Vectors to store the rows from the neighbors
        vector<int> row_up(grid_size);
        vector<int> row_down(grid_size);
        // Simultaenously send my row and receive rows from neighbors
        MPI_Sendrecv(
            row.data(), grid_size, MPI_INT,  neighbor_up_rank, 0, 
            row_down.data(), grid_size, MPI_INT, neighbor_down_rank, 0,
            MPI_COMM_WORLD, &status
        );
        MPI_Sendrecv(
            row.data(), grid_size, MPI_INT, neighbor_down_rank, 0, 
            row_up.data(), grid_size, MPI_INT, neighbor_up_rank, 0,
            MPI_COMM_WORLD, &status
        );
        // Vector to store the updated row
        vector<int> updated_row(grid_size);
        for (int i = 0; i < grid_size; ++i) {
            vector<int> neighbors(NUM_NEIGHBORS);
            // Get the neighbors of the current cell - Wrap-around implementation
            neighbors[0] = row_up[i];                                   // N
            neighbors[1] = row_up[(i + 1) % grid_size];                 // NE
            neighbors[2] = row_up[(i - 1 + grid_size) % grid_size];     // NW
            neighbors[3] = row[(i - 1 + grid_size) % grid_size];        // W
            neighbors[4] = row[(i + 1) % grid_size];                    // E
            neighbors[5] = row_down[i];                                 // S
            neighbors[6] = row_down[(i + 1) % grid_size];               // SE
            neighbors[7] = row_down[(i - 1 + grid_size) % grid_size];   // SW

            // Update the cell value
            updated_row[i] = updateCellValue(row[i], neighbors);
        }
        // Update the row with the new values
        row = updated_row;

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Vector to store all rows in the root process
    vector<int> all_rows;
    if (rank == 0) {
        all_rows.resize(grid_size * num_procs);
    }

    // Gather all rows to the root process
    MPI_Gather(row.data(), grid_size, MPI_INT, all_rows.data(), grid_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Print all rows in the root process
    int total_population = 0;
    if (rank == 0) {
        for (int i = 0; i < num_procs; ++i) {
            for (int j = 0; j < grid_size; ++j) {
                total_population += all_rows[i * grid_size + j] ? 1 : 0;
                cout << all_rows[i * grid_size + j];
            }
            cout << endl;
        }
        //cout << "Total population: " << total_population << endl;
    }

    MPI_Finalize();
    return 0;
}