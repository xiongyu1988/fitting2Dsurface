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

class ReadParseData:
    def __init__(self):
        self.grid_data = {}
        self.mesh_elements = {}

    def parse_grid_line(self, line):
        if line.startswith("GRID"):
            node_id = int(line[8:16].strip())
            x = float(line[24:32].strip())
            y = float(line[32:40].strip())
            z_str = line[40:48].strip()
            if 'e' not in z_str and ('+' in z_str or '-' in z_str[1:]):
                z_str = z_str[:-3] + 'e' + z_str[-3:]
            z = float(z_str)
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
    Z = coeffs[0] * (Y - coeffs[1])**2 + coeffs[2] * X + coeffs[3]
    ax.plot_surface(X, Y, Z, color='g', alpha=0.6, edgecolor='none')

def plot_flat_panel(ax, coeffs, x_range, y_range):
    x = np.linspace(x_range[0], x_range[1], 100)
    y = np.linspace(y_range[0], y_range[1], 100)
    X, Y = np.meshgrid(x, y)
    Z = coeffs[0] + coeffs[1] * X + coeffs[2] * Y
    ax.plot_surface(X, Y, Z, color='y', alpha=0.6, edgecolor='none')

if __name__ == "__main__":
    filename = 'data/testGeo.fem'  # Manually specify the filename here

    read_parse_data = ReadParseData()
    read_parse_data.read_from_file(filename)
    surface_area = read_parse_data.calculate_surface_area()
    print(f"Total surface area: {surface_area:.6f}")

    # Coefficients and ranges for the surfaces (example values, replace with actual values)
    doubly_curved_coeffs = [-11.209817, -0.000904, -2.865284, 0.000055, 0.000078, -0.093105]
    singly_curved_coeffs = [-0.093109, -15.383293, -0.001024, 10.832004]
    flat_panel_coeffs = [5.703799, -0.003804, -0.001578]

    x_range = [0.0, 20.0]
    y_range = [-25.3859, -5.38594]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    read_parse_data.plot_mesh(ax)
    # Uncomment the surface you want to plot
    # plot_doubly_curved_shell(ax, doubly_curved_coeffs, x_range, y_range)
    # plot_singly_curved_shell(ax, singly_curved_coeffs, x_range, y_range)
    # plot_flat_panel(ax, flat_panel_coeffs, x_range, y_range)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()
