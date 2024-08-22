import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Function to plot the elliptic paraboloid with principal curvatures highlighted
def plot_elliptic_paraboloid_with_curvatures(Rx=5, Ry=3):
    # Calculate a and b based on the desired radii Rx and Ry
    a = 1 / (2 * Rx)
    b = 1 / (2 * Ry)
    
    # Generate mesh grid
    x = np.linspace(-5, 5, 100)
    y = np.linspace(-5, 5, 100)
    x, y = np.meshgrid(x, y)
    z = a * x**2 + b * y**2

    # Plotting the surface
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x, y, z, cmap='viridis', edgecolor='none', alpha=0.8)

    # Highlighting the principal curvatures along the x and y axes
    x_curve = np.linspace(-5, 5, 100)
    y_curve = np.zeros_like(x_curve)
    z_x_curve = a * x_curve**2
    ax.plot(x_curve, y_curve, z_x_curve, color='red', linewidth=3, label=f'Curvature in X (Rx={Rx})')

    y_curve = np.linspace(-5, 5, 100)
    x_curve = np.zeros_like(y_curve)
    z_y_curve = b * y_curve**2
    ax.plot(x_curve, y_curve, z_y_curve, color='blue', linewidth=3, label=f'Curvature in Y (Ry={Ry})')

    # Adding labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Elliptic Paraboloid with Rx={Rx}, Ry={Ry}')

    # Adding legend
    ax.legend()

    # Annotating the radii of curvature
    ax.text(0, 0, a * (5**2), f'Rx = {Rx}', color='red', fontsize=12)
    ax.text(0, 0, b * (5**2), f'Ry = {Ry}', color='blue', fontsize=12)

    plt.show()

# Plot with specific Rx and Ry values
plot_elliptic_paraboloid_with_curvatures(Rx=5, Ry=3)
