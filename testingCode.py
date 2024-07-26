import numpy as np
from scipy.linalg import lstsq

def parse_coordinate(coord_str):
    if '.' in coord_str[1:]:
        parts = coord_str.split('.')
        if len(parts) == 3:
            mantissa = float(parts[0] + '.' + parts[1])
            exponent = int(parts[2])
            return mantissa * (10 ** -exponent)
    return float(coord_str)

def read_fem_file(filename):
    grid_data = {}
    mesh_elements = {}
    
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('GRID'):
                parts = line.split()
                node_id = int(parts[1])
                x = parse_coordinate(line[24:32].strip())
                y = parse_coordinate(line[32:40].strip())
                z = parse_coordinate(line[40:48].strip())
                grid_data[node_id] = np.array([x, y, z])
            elif line.startswith('CTRIA3'):
                parts = line.split()
                element_id = int(parts[1])
                node1, node2, node3 = map(int, parts[3:6])
                mesh_elements[element_id] = [node1, node2, node3]
    
    return grid_data, mesh_elements

def calculate_surface_area(grid_data, mesh_elements):
    total_area = 0
    for element_id, nodes in mesh_elements.items():
        p1, p2, p3 = [grid_data[node] for node in nodes]
        v1 = p2 - p1
        v2 = p3 - p1
        cross = np.cross(v1, v2)
        area = 0.5 * np.linalg.norm(cross)
        total_area += area
        if element_id <= 5:
            print(f"Element {element_id} area: {area}")
    return total_area

def fit_doubly_curved_shell(grid_data):
    points = np.array(list(grid_data.values()))
    x, y, z = points[:, 0], points[:, 1], points[:, 2]
    A = np.column_stack([np.ones_like(x), x, y, x**2, x*y, y**2])
    coeffs, _, _, _ = lstsq(A, z)
    return coeffs

def main():
    filename = "data/geo1.fem"
    grid_data, mesh_elements = read_fem_file(filename)
    
    print("Sample Grid Data:")
    for i, (node_id, coords) in enumerate(grid_data.items()):
        print(f"Node {node_id}: ({coords[0]:.6f}, {coords[1]:.6f}, {coords[2]:.6f})")
        if i >= 4:
            break
    
    print("\nSample Mesh Elements:")
    for i, (element_id, nodes) in enumerate(mesh_elements.items()):
        print(f"Element {element_id}: {' '.join(map(str, nodes))}")
        if i >= 4:
            break
    
    surface_area = calculate_surface_area(grid_data, mesh_elements)
    print(f"\nMesh surface area: {surface_area:.6f}")
    
    shell_coeffs = fit_doubly_curved_shell(grid_data)
    print("\nShell coefficients:")
    print(shell_coeffs)

if __name__ == "__main__":
    main()