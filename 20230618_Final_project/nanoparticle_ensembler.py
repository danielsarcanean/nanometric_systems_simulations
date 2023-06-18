##################################################
# SPHERICAL CERIUM STACKING MONTECARLO SIMULATOR #
# MADE BY: DANIEL SARCANEAN                      #
##################################################


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import imageio

# Constants
k_B = 1.38e-23  # Boltzmann constant in J/K
T = 1000  # Temperature in K

# Lennard-Jones potential parameters for cerium (KIM Models)
epsilon = 4.11772195e-19    # Depth of the potential well in J
sigma = 3.63487e-10  # Distance at which the potential is zero in m

# Initialize the 3D grid and select the middle point as a catalyzer/initiator
grid_size = 100
grid = np.zeros((grid_size, grid_size, grid_size))
grid[grid_size//2, grid_size//2, grid_size//2] = 1
center = np.array([grid_size // 2, grid_size // 2, grid_size // 2])
radius = 4.5479 # Real radius is 15 - normalization factor (r/3.65)

# Monte Carlo steps (Choose accordingly for better coverage)
n_steps = 2000
# List to store file names
filenames = []

def lennard_jones_potential(r):
    return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)

# Return fictice temperature factor 
def calculate_temperature_factor(x, y, z, center, max_radius):
    distance_to_center = np.linalg.norm(np.array([x, y, z]) - center)
    # Normalize the distance to the range [0, 1]
    normalized_distance = distance_to_center / max_radius
    # Calculate the temperature factor as a function of the normalized distance
    # For example, you could use a linear function, but other functions are possible
    temperature_factor = normalized_distance
    return temperature_factor

     
def get_neighbors(x, y, z):
    # Returns a list of the coordinates of the neighbors of the point (x, y, z)
    neighbors = []
    for dx in [-1, 0, 1]:
        for dy in [-1, 0, 1]:
            for dz in [-1, 0, 1]:
                if dx != 0 or dy != 0 or dz != 0:  # Exclude the point itself
                    neighbors.append(((x + dx) % grid_size, (y + dy) % grid_size, (z + dz) % grid_size))
    return neighbors

''' 
# Only use to see different values of epsilon and sigma in LJ potential affinity
def calculate_crystallinity(x, y, z):
    # Calculate the "crystallinity" of the point (x, y, z) based on the number of nearest neighbors
    return sum(grid[i, j, k] for i, j, k in get_neighbors(x, y, z))
'''

for step in range(n_steps):
    # Create a list of all the grid points that are adjacent to existing atoms
    adjacent_points = []
    for x in range(grid_size):
        for y in range(grid_size):
            for z in range(grid_size):
                if grid[x, y, z] == 1:
                    adjacent_points.extend(get_neighbors(x, y, z))
                    
    # Remove duplicates and points that already contain an atom
    adjacent_points = list(set(adjacent_points) - set(zip(*np.where(grid == 1))))

    # Randomly select a grid point from the list of adjacent points, or from the entire grid if there are no adjacent points
    if adjacent_points:
        x, y, z = adjacent_points[np.random.randint(0, len(adjacent_points))]
    else:
        x, y, z = np.random.randint(0, grid_size, size=3)
        
    # Check if the addition of the new atom would keep the particle spherical
    distance_to_center = np.linalg.norm(np.array([x, y, z]) - center)
    if distance_to_center > radius:
        continue  # Skip this iteration and do not add the atom
    
    # Calculate the energy change due to the addition of another atom
    delta_E = 0
    for i, j, k in get_neighbors(x, y, z):
        if grid[i, j, k] == 1:
            r = np.sqrt((x - i)**2 + (y - j)**2 + (z - k)**2)
            delta_E += lennard_jones_potential(r)
    #delta_E -= calculate_crystallinity(x, y, z)  # Decrease the energy when the number of nearest neighbors increases

    # Metropolis-Hastings algorithm
    if delta_E < 0:
        # If the energy decreases, accept the change
        grid[x, y, z] = 1
    else:
        # If the energy increases, accept the change with a probability given by the Boltzmann factor
        if np.random.uniform(0, 1) < np.exp(-delta_E / (k_B * T)):
            grid[x, y, z] = 1

    # Visualize the grid
    if step % 50 == 0:  # Update the plot every 1000 steps
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        x, y, z = np.where(grid == 1)
        ax.scatter(x, y, z)
        filename = f'step_{step}.png'
        plt.savefig(filename)
        filenames.append(filename)
        plt.close()

# Create a gif from the images
with imageio.get_writer('simulation.gif', mode='I') as writer:
    for filename in filenames:
        image = imageio.v2.imread(filename)
        writer.append_data(image)

grid_spacing = 3.65 # Angstroms (an extended factor is used to mimick 
# experimentally accurate data and prevent automatic bonding by VMD)

# Save protein data bank file
with open('structure1.pdb', 'w') as f:
    f.write('MODEL\n')
    x, y, z = np.where(grid == 1)
    num_atoms = np.sum(grid)
    for i in range(len(x)):
        x_coord = "{:8.3f}".format(x[i] * grid_spacing)
        y_coord = "{:8.3f}".format(y[i] * grid_spacing)
        z_coord = "{:8.3f}".format(z[i] * grid_spacing)
        e = "{:4d}".format(i+1)
        h = "{:5d}".format(i+1)
        temp_factor = 10*calculate_temperature_factor(x[i], y[i], z[i], center, radius)
        temp_factor_str = "{:6.2f}".format(temp_factor)
        f.write(f'ATOM  {h}  CE  ICEM {e}    {x_coord}{y_coord}{z_coord}  1.00{temp_factor_str}      CEM\n')
    
    f.write(f'TER  {num_atoms+1:4.0f}      ICEM  {num_atoms + 1:0f}\n')
    f.write('ENDMDL\n')


# Save grid data in .xyz format
with open('structure1.xyz', 'w') as f:
    num_atoms = np.sum(grid)  # Count the number of atoms in the grid
    f.write(f'{num_atoms}\n')
    f.write('Atoms\n')
    x, y, z = np.where(grid == 1)
    for i in range(len(x)):
        x_coord = x[i] * grid_spacing
        y_coord = y[i] * grid_spacing
        z_coord = z[i] * grid_spacing
        f.write(f'Ce {x_coord} {y_coord} {z_coord}\n')
