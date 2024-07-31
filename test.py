import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define the mesh grid for u and v
u = np.linspace(-10, 10, 1000)
v = np.linspace(-20, 20, 1000)
u, v = np.meshgrid(u, v)

R1 = 20  # Known radius of curvature
R2 = 40  #

# Define the parameters a and b
a = np.sqrt(2*R1)  # a = sqrt(2*R1) where R1 is unknown
b = np.sqrt(2*R2)  # b = sqrt(2*R2) where R2 is unknown

# Define the elliptic paraboloid equation
z = (u**2 / a**2) + (v**2 / b**2)

# Calculate partial derivatives
zu = 2 * u / a**2
zv = 2 * v / b**2

# First fundamental form coefficients
E = 1 + zu**2
G = 1 + zv**2

# Second derivatives
zuu = 2 / a**2
zvv = 2 / b**2

# Normal vector components
nx = -zu
ny = -zv
nz = np.ones_like(z)
n = np.sqrt(nx**2 + ny**2 + nz**2)
nx /= n
ny /= n
nz /= n

# Second fundamental form coefficients
L = zuu / n
N = zvv / n

# Principal curvatures
K1 = L / E
K2 = N / G

# Radii of curvature
R1_values = 1 / K1
R2_values = 1 / K2

# Angles relative to the normal vector
angles = np.degrees(np.arccos(nz))

# Plotting the surface
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(u, v, z, cmap='viridis', edgecolor='none')

# Set plot limits
ax.set_xlim([-10, 10])
ax.set_ylim([-20, 20])
ax.set_zlim([0, max(20, 40)])

# Labels and title
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
ax.set_title('Elliptic Paraboloid with Calculated Radii of Curvature')

# Show the plot
plt.show()

# Display results
print(f'Calculated Mean Radius of Curvature R1: {np.mean(R1_values)}')
print(f'Calculated Mean Radius of Curvature R2: {np.mean(R2_values)}')
print(f'Mean Angle Relative to Normal Vector: {np.mean(angles)} degrees')
