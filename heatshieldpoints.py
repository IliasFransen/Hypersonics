# import numpy as np
# import matplotlib.pyplot as plt
#
# def generate_heatshield_points(radius_of_curvature, height, diameter):
#     # Calculate the semi-axes of the ellipse
#     a = diameter / 2
#     b = height
#
#     # Generate points along the ellipse
#     theta = np.linspace(0, np.pi, 1000)
#     x = a * np.cos(theta)
#     y = b * np.sin(theta) - b + radius_of_curvature
#
#     return x, y
#
# # Example usage:
# radius_of_curvature = 6.03504
# height = 12
# diameter = 5.03
#
# x, y = generate_heatshield_points(radius_of_curvature, height, diameter)
#
# # Plot the points
# plt.figure(figsize=(6,6))
# plt.plot(x, y)
# plt.xlabel('X')
# plt.ylabel('Y')
# plt.title('Heatshield Profile')
# plt.grid(True)
# plt.show()

import numpy as np
import matplotlib.pyplot as plt

def generate_heatshield_points(radius, height, diameter):
    # Calculate the x points for a circle
    theta = np.linspace(0, 2*np.pi, 100)
    x_circle = radius * np.cos(theta)

    # Calculate the y points for a circle
    y_circle = radius * np.sin(theta)

    # Scale the y points based on the height and diameter
    y_scaled = y_circle * (height / diameter)

    # Plot the heatshield
    plt.plot(x_circle, y_scaled, label='Heatshield')
    plt.title('Space Crew Capsule Heatshield')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.legend()
    plt.axis('equal')
    plt.grid(True)
    plt.show()

# Example usage
radius_of_curvature = 10  # Adjust as needed
height_of_heatshield = 5  # Adjust as needed
diameter_of_capsule = 8   # Adjust as needed

generate_heatshield_points(radius_of_curvature, height_of_heatshield, diameter_of_capsule)