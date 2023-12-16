import numpy as np
import matplotlib.pyplot as plt

<<<<<<< Updated upstream
def generate_heatshield_points(radius_of_curvature, diameter, height, num_points=100):
=======
def generate_heatshield_points(D, H, num_points=51):
>>>>>>> Stashed changes
    # Parametric equations for an ellipse
    theta = np.linspace(0, 2*np.pi, num_points)
    x = diameter/2 * np.cos(theta)
    y = height/2 * np.sin(theta)

    # Apply curvature to the points
<<<<<<< Updated upstream
    x_curved = x + radius_of_curvature - diameter/2
    x,y = x_curved[y <= 0], y[y <= 0]
    return x, y

# def plot_apollo_heatshield(x, y):
#     plt.scatter(x, y)
#     plt.title("Apollo Command Module Heatshield")
#     plt.xlabel("X-axis")
#     plt.ylabel("Y-axis")
#     plt.gca().set_aspect('equal', adjustable='box')
#     plt.grid(True)
#     plt.show()

=======
    #plot_apollo_heatshield(x,y)
    return x, y


def plot_apollo_heatshield(x, y):
     plt.scatter(x, y)
     plt.title("Apollo Command Module Heatshield")
     plt.xlabel("X-axis")
     plt.ylabel("Y-axis")
     plt.gca().set_aspect('equal', adjustable='box')
     plt.grid(True)
     plt.show()

>>>>>>> Stashed changes
# # Example usage
# radius_of_curvature = 6  # Replace with your actual value
# diameter = 4  # Replace with your actual value
# height = 0.5  # Replace with your actual value
#
# x_points, y_points = generate_heatshield_points(radius_of_curvature, diameter, height)
# plot_apollo_heatshield(x_points,y_points)