import re
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

def fix_scientific_notation(value_str):
    """
    Fixes scientific notation in cases like '123-03' to '123e-03'.
    """
    sci_notation_pattern = re.compile(r'([-+]?\d*\.?\d+)([-+]\d+)')
    match = sci_notation_pattern.search(value_str)
    
    if match:
        # If a match is found, reconstruct the correct format
        return match.group(1) + 'e' + match.group(2)
    
    # If no match, return the original string
    return value_str

class ReadParseData:
    def __init__(self):
        self.grid_data = {}
        self.mesh_elements = {}



    def parse_grid_line(self, line):
        if line.startswith("GRID"):
            node_id = int(line[8:16].strip())
            x_str = line[24:32].strip()
            y_str = line[32:40].strip()
            z_str = line[40:48].strip()

            # Apply the fix_scientific_notation to x, y, and z strings
            x_str = fix_scientific_notation(x_str)
            y_str = fix_scientific_notation(y_str)
            z_str = fix_scientific_notation(z_str)

            # Convert to floats
            x = float(x_str)
            y = float(y_str)
            z = float(z_str)

            # Store the parsed coordinates
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
        try:
            with open(filename, 'r') as file:
                for line in file:
                    self.parse_grid_line(line)
                    self.parse_mesh_line(line)
        except FileNotFoundError:
            print(f"Error: File '{filename}' not found.")
        except Exception as e:
            print(f"Error: {e}")

    def calculate_surface_area(self):
        total_area = 0.0
        for element, nodes in self.mesh_elements.items():
            if len(nodes) == 3:
                p1 = self.grid_data[nodes[0]]
                p2 = self.grid_data[nodes[1]]
                p3 = self.grid_data[nodes[2]]
                v1 = p2 - p1
                v2 = p3 - p1
                cross_product = np.cross([v1.x, v1.y, v1.z], [v2.x, v2.y, v2.z])
                area = 0.5 * np.linalg.norm(cross_product)
                total_area += area
        return total_area

    def plot_mesh(self, ax):
        for element, nodes in self.mesh_elements.items():
            x = [self.grid_data[node].x for node in nodes]
            y = [self.grid_data[node].y for node in nodes]
            z = [self.grid_data[node].z for node in nodes]
            x.append(x[0])  # Close the loop
            y.append(y[0])
            z.append(z[0])
            ax.plot(x, y, z, color='b')

        # Plot nodal points
        x = [point.x for point in self.grid_data.values()]
        y = [point.y for point in self.grid_data.values()]
        z = [point.z for point in self.grid_data.values()]
        ax.scatter(x, y, z, color='r', s=10)  # Nodal points in red

def plot_doubly_curved_shell(ax, coeffs, x_range, y_range):
    x = np.linspace(x_range[0], x_range[1], 100)
    y = np.linspace(y_range[0], y_range[1], 100)
    X, Y = np.meshgrid(x, y)
    Z = (coeffs[0] + coeffs[1] * X + coeffs[2] * Y +
         coeffs[3] * X**2 + coeffs[4] * X * Y + coeffs[5] * Y**2)
    ax.plot_surface(X, Y, Z, color='r', alpha=0.6, edgecolor='none')

def plot_singly_curved_shell(ax, coeffs, x_range, y_range):
    x = np.linspace(x_range[0], x_range[1], 100)
    y = np.linspace(y_range[0], y_range[1], 100)
    X, Y = np.meshgrid(x, y)
    Z = coeffs[0] * (X - coeffs[1])**2 + coeffs[2] * Y + coeffs[3]
    ax.plot_surface(X, Y, Z, color='g', alpha=0.6, edgecolor='none')

def plot_singly_curved_shell2(ax, coeffs, x_range, y_range):
    x = np.linspace(x_range[0], x_range[1], 100)
    y = np.linspace(y_range[0], y_range[1], 100)
    X, Y = np.meshgrid(x, y)
    Z =coeffs[0] * Y**2 + coeffs[1] * X 
    ax.plot_surface(X, Y, Z, color='g', alpha=0.6, edgecolor='none')

def plot_doubly_curved_shell2(ax, coeffs, x_range, y_range):
    x = np.linspace(x_range[0], x_range[1], 100)
    y = np.linspace(y_range[0], y_range[1], 100)
    X, Y = np.meshgrid(x, y)
    Z = coeffs[0] * X**2 + coeffs[1] * X + coeffs[2] * Y + coeffs[3] * Y**2
    ax.plot_surface(X, Y, Z, color='g', alpha=0.6, edgecolor='none')

def plot_flat_panel(ax, coeffs, x_range, y_range):
    x = np.linspace(x_range[0], x_range[1], 100)
    y = np.linspace(y_range[0], y_range[1], 100)
    X, Y = np.meshgrid(x, y)
    Z = coeffs[0] + coeffs[1] * X + coeffs[2] * Y
    ax.plot_surface(X, Y, Z, color='y', alpha=0.6, edgecolor='none')

if __name__ == "__main__":
    filename = 'data/doublyCurvedShell_2.fem'  # Manually specify the filename here

    read_parse_data = ReadParseData()
    read_parse_data.read_from_file(filename)
    surface_area = read_parse_data.calculate_surface_area()
    print(f"Total surface area: {surface_area:.6f}")

    # Coefficients and ranges for the surfaces (example values, replace with actual values)
    doubly_curved_coeffs = [60.474667, -25.030457, -0.000741, 2.622955, 0.006062, 1.198502]
    singly_curved_coeffs = [0.427566, -0.024748, -0.163008, 1.655713]
    flat_panel_coeffs = [1.524741, -0.129304, 0.015091]

    #doubly_curved_coeffs2 = [-0.072603, 0.531280, 0.020841, 0.401458]
    doubly_curved_coeffs2 = [-0.117619, 0.410272, -0.000002, -7.154899]
    singly_curved_coeffs2 = [-4.317719, 0.127801]

    x_range = [2.093541, 2.527030]
    y_range = [-0.136805, 0.136805]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    read_parse_data.plot_mesh(ax)
    # Uncomment the surface you want to plot
    #plot_doubly_curved_shell(ax, doubly_curved_coeffs, x_range, y_range)
    #plot_doubly_curved_shell2(ax, doubly_curved_coeffs2, x_range, y_range)
    #plot_flat_panel(ax, flat_panel_coeffs, x_range, y_range)
    plot_singly_curved_shell2(ax, singly_curved_coeffs2, x_range, y_range)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
