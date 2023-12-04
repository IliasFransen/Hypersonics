import matplotlib.pyplot as plt
import numpy as np

def create_curved_plate(radius, width, height):
    # Generate theta values for the curve
    theta = np.linspace(0, 2*np.pi, 200)

    # Parametric equations for a circle
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)

    # Translate the circle to create the plate
    y = y + radius-height
    x = x[y <= 0]
    y = y[y <= 0]
    return x, y

def plot_curved_plate(x, y):
    plt.plot(x, y, label='Curved Plate')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.legend()
    plt.axis('equal')
    plt.grid(True)
    plt.show()

# Example usage:
radius = 4.69392
width = 3.9116
height = 0.635

x_plate, y_plate = create_curved_plate(radius, width, height)
plot_curved_plate(x_plate, y_plate)
