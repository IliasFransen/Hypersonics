import numpy as np
import matplotlib.pyplot as plt

def generate_heatshield_points(D, H, num_points=20):
    # Parametric equations for an ellipse
    dth = 0.0000001
    R=D/2
    th = np.linspace(np.pi+dth, 2*np.pi-dth, num_points)
    x = R*np.cos(th)
    y = H*np.sin(th)
    # Apply curvature to the points
    return x, y

def plot_heatshield(x, y):
    plt.scatter(x, y)
    plt.title("Apollo Command Module Heatshield")
    plt.xlabel("X-axis")
    plt.ylabel("Y-axis")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True)
    plt.show()

# Example usage
# D = 2  # Replace with your actual value
# H = 0.5  # Replace with your actual value
#
# x_points, y_points = generate_heatshield_points(D, H, 11)
# plot_heatshield(x_points,y_points)