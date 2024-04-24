import numpy as np
import sys 

def random_grid(grid_size):
    grid = np.random.randint(2, size=(grid_size, grid_size))
    # Write the grid to a file
    output_filename = f"random{grid_size}"
    with open(output_filename, 'w') as f:
        for i in range(grid_size):
            for j in range(grid_size):
                f.write(str(grid[i][j]))
            f.write('\n')

random_grid(int(sys.argv[1]))
