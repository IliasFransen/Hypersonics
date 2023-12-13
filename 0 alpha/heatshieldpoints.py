import numpy as np
import matplotlib.pyplot as plt

def generate_heatshield_points(D, H, num_points=101):
    # Parametric equations for an ellipse
    R = D / 2
    np.linspace(-R, R, num_points)
    th = np.linspace(np.pi, 2 * np.pi, num_points)
    x = R * np.cos(th)
    y = H * np.sin(th)
    # Apply curvature to the points
    return x, y

"""
def plot_apollo_heatshield(x, y):
     plt.scatter(x, y)
     plt.title("Apollo Command Module Heatshield")
     plt.xlabel("X-axis")
     plt.ylabel("Y-axis")
     plt.gca().set_aspect('equal', adjustable='box')
     plt.grid(True)
     plt.show()
"""
# # Example usage
# radius_of_curvature = 6  # Replace with your actual value
# diameter = 4  # Replace with your actual value
# height = 0.5  # Replace with your actual value
#
# x_points, y_points = generate_heatshield_points(radius_of_curvature, diameter, height)
# plot_apollo_heatshield(x_points,y_points)