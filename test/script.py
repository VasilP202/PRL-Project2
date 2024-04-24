import numpy as np
import sys


def game_of_life(grid, num_steps):
    grid_size = len(grid)
    for step in range(num_steps):
        new_grid = np.zeros((grid_size, grid_size), dtype=int)
        for i in range(grid_size):
            for j in range(grid_size):
                # Count the number of live neighbors
                live_neighbors = 0
                for x in range(i-1, i+2):
                    for y in range(j-1, j+2):
                        # Wrap around the grid
                        x = (x + grid_size) % grid_size
                        y = (y + grid_size) % grid_size
                        if (x != i or y != j) and grid[x][y] == 1:
                            live_neighbors += 1
                # Apply the rules of the game
                if grid[i][j] == 1:
                    if live_neighbors < 2 or live_neighbors > 3:
                        new_grid[i][j] = 0
                    else:
                        new_grid[i][j] = 1
                else:
                    if live_neighbors == 3:
                        new_grid[i][j] = 1
        grid = new_grid
    
    return grid        


# Read the filename from the command line
filename = sys.argv[1]
num_steps = int(sys.argv[2])

# Load the data
f = open(filename)

line = f.readline().strip()

data = []
while line:
    data.append([int(c) for c in line])
    line = f.readline().strip()

grid_size = len(data)
grid = game_of_life(data, num_steps)

for i in range(grid_size):
    print(f"{i}: ", end='')
    for j in range(grid_size):
        print(grid[i][j], end='')
    print()