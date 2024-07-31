import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Vector3D:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

class ReadParseData:
    def __init__(self):
        self.grid_data = {}
        self.mesh_elements = {}

    def parse_grid_line(self, line):
        if line.startswith("GRID"):
            node_id = int(line[8:16].strip())
            x = float(line[24:32].strip())
            y = float(line[32:40].strip())
            z = float(line[40:48].strip())
            self.grid_data[node_id] = Vector3D(x, y, z)

    def parse_mesh_line(self, line):
        if line.startswith("CTRIA3"):
            parts = line.split()
            element_id = int(parts[1])
            node1 = int(parts[3])
            node2 = int(parts[4])
            node3 = int(parts[5])
            self.mesh_elements[element_id] = [node1, node2, node3]

    def read_from_file(self, filename):
        with open(filename, 'r') as file:
            for line in file:
                self.parse_grid_line(line)
                self.parse_mesh_line(line)

    def plot_geometry(self, ax):
        # Plot nodes
        x = [node.x for node in self.grid_data.values()]
        y = [node.y for node in self.grid_data.values()]
        z = [node.z for node in self.grid_data.values()]
        ax.scatter(x, y, z, c='b', marker='o', s=10, label='FEM Nodes')

        # Plot mesh elements
        for element in self.mesh_elements.values():
            for i in range(3):
                node1 = self.grid_data[element[i]]
                node2 = self.grid_data[element[(i+1)%3]]
                ax.plot([node1.x, node2.x], [node1.y, node2.y], [node1.z, node2.z], 'r-', linewidth=0.5)

    def get_coordinate_ranges(self):
        x = [node.x for node in self.grid_data.values()]
        y = [node.y for node in self.grid_data.values()]
        z = [node.z for node in self.grid_data.values()]
        return (min(x), max(x)), (min(y), max(y)), (min(z), max(z))

# Define the equation parameters
h = 4.725783
k = -0.011001
a = 0.677848
b = 1.065240
z0 = 0.796295

# Set up the plot
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')

# Read and plot FEM geometry
analyzer = ReadParseData()
analyzer.read_from_file("data/cylinder.fem")
analyzer.plot_geometry(ax)

# Get coordinate ranges from FEM data
(x_min, x_max), (y_min, y_max), (z_min, z_max) = analyzer.get_coordinate_ranges()

# Define the coordinate ranges for the elliptic paraboloid
x_range = np.linspace(x_min, x_max, 100)
y_range = np.linspace(y_min, y_max, 100)
x, y = np.meshgrid(x_range, y_range)

# Calculate z based on the elliptic paraboloid equation
z = ((x - h)**2 / a**2) + ((y - k)**2 / b**2) + z0

# Plot the elliptic paraboloid surface
surf = ax.plot_surface(x, y, z, cmap='viridis', alpha=0.7, label='Fitted Surface')

# Set the coordinate ranges
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
ax.set_zlim(z_min, z_max)

# Label the axes
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Add a title
ax.set_title('FEM Geometry and Fitted Elliptic Paraboloid')

# Add a color bar
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)

# Add a legend
ax.legend()

# Show the plot
plt.show()