import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Vector3D:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __sub__(self, other):
        return Vector3D(self.x - other.x, self.y - other.y, self.z - other.z)

class ReadParseData:
    def __init__(self):
        self.gridData = {}
        self.meshElements = {}

    def parse_grid_line(self, line):
        if line.startswith("GRID"):
            nodeId = int(line[8:16].strip())
            x = float(line[24:32].strip())
            y = float(line[32:40].strip())
            z = float(line[40:48].strip())
            self.gridData[nodeId] = Vector3D(x, y, z)

    def parse_mesh_line(self, line):
        if line.startswith("CTRIA3"):
            parts = line.split()
            elementId = int(parts[1])
            node1 = int(parts[3])
            node2 = int(parts[4])
            node3 = int(parts[5])
            self.meshElements[elementId] = [node1, node2, node3]

    def read_from_file(self, filename):
        with open(filename, 'r') as file:
            for line in file:
                self.parse_grid_line(line)
                self.parse_mesh_line(line)

    def calculate_surface_area(self):
        total_area = 0.0
        for element, nodes in self.meshElements.items():
            if len(nodes) == 3:
                p1 = self.gridData[nodes[0]]
                p2 = self.gridData[nodes[1]]
                p3 = self.gridData[nodes[2]]
                v1 = p2 - p1
                v2 = p3 - p1
                cross_product = np.cross([v1.x, v1.y, v1.z], [v2.x, v2.y, v2.z])
                area = 0.5 * np.linalg.norm(cross_product)
                total_area += area
        return total_area

    def plot_mesh(self, ax):
        for element, nodes in self.meshElements.items():
            x = [self.gridData[node].x for node in nodes]
            y = [self.gridData[node].y for node in nodes]
            z = [self.gridData[node].z for node in nodes]
            x.append(x[0])  # Close the loop
            y.append(y[0])
            z.append(z[0])
            ax.plot(x, y, z, color='b')

        # Plot nodal points
        x = [point.x for point in self.gridData.values()]
        y = [point.y for point in self.gridData.values()]
        z = [point.z for point in self.gridData.values()]
        ax.scatter(x, y, z, color='r', s=10)  # Nodal points in red

# Function to plot doubly curved shell
def plot_doubly_curved_shell(ax, coeffs, x_range, y_range):
    x = np.linspace(x_range[0], x_range[1], 100)
    y = np.linspace(y_range[0], y_range[1], 100)
    X, Y = np.meshgrid(x, y)
    Z = (coeffs[0] + coeffs[1] * X + coeffs[2] * Y +
         coeffs[3] * X**2 + coeffs[4] * X * Y + coeffs[5] * Y**2)
    ax.plot_surface(X, Y, Z, color='r', alpha=0.6, edgecolor='none')

# Function to plot singly curved shell
def plot_singly_curved_shell(ax, coeffs, x_range, y_range):
    x = np.linspace(x_range[0], x_range[1], 100)
    y = np.linspace(y_range[0], y_range[1], 100)
    X, Y = np.meshgrid(x, y)
    Z = coeffs[0] * (X - coeffs[1])**2 + coeffs[2] * Y + coeffs[3]
    ax.plot_surface(X, Y, Z, color='g', alpha=0.6, edgecolor='none')

# Function to plot flat panel surface
def plot_flat_panel(ax, coeffs, x_range, y_range):
    x = np.linspace(x_range[0], x_range[1], 100)
    y = np.linspace(y_range[0], y_range[1], 100)
    X, Y = np.meshgrid(x, y)
    Z = coeffs[0] + coeffs[1] * X + coeffs[2] * Y
    ax.plot_surface(X, Y, Z, color='y', alpha=0.6, edgecolor='none')

# Usage
if __name__ == "__main__":
    read_parse_data = ReadParseData()
    read_parse_data.read_from_file('data/testGeo.fem')
    surface_area = read_parse_data.calculate_surface_area()
    print(f"Total surface area: {surface_area:.6f}")

    # Coefficients and ranges for the surfaces (example values, replace with actual values)
    #doubly_curved_coeffs = [60.474667 , -25.030457, -0.000741, 2.622955, 0.006062, 1.198502]
    doubly_curved_coeffs = [4.254311 , -1.281356, -0.374673, 0.121153, 0.082139, 0.087465]

    singly_curved_coeffs = [0.118743, -1.252326, 0.014910, 4.176812]
    flat_panel_coeffs = [1.524741, -0.129304, 0.015091]
    x_range = [4.470440, 4.996350]
    y_range = [-0.372765, 0.377163]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    read_parse_data.plot_mesh(ax)
    #plot_doubly_curved_shell(ax, doubly_curved_coeffs, x_range, y_range)
    plot_singly_curved_shell(ax, singly_curved_coeffs, x_range, y_range)
    #plot_flat_panel(ax, flat_panel_coeffs, x_range, y_range)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()