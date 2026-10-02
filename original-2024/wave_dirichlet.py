import numpy as np
import matplotlib.pyplot as plt

# Defining Parameters
max_t = 10  # maximum time
max_x = max_y = 10.0  # spatial domain
dt = 0.01  # time step
dx = dy = 0.1  # spatial step
num_t = int(max_t / dt) + 1  # number of time steps
num_x = int(max_x / dx) + 1  # number of x steps
num_y = int(max_y / dy) + 1  # number of y steps

# Initialize the grid
w = np.zeros((num_t, num_x, num_y))

# Defining a fixed Initial condition with Gaussian pulse
for j in range(num_x):
    for k in range(num_y):
        x = j * dx - dx / 2
        y = k * dy
        w[0, j, k] = np.exp(-((x - 5)**2 + (y - 5)**2))

# Copying the initial condition to the first time step
w[1, :, :] = w[0, :, :]

# Define Updating Wave Function with Finite difference method
for n in range(1, num_t - 1):
    for j in range(1, num_x - 1):
        for k in range(1, num_y - 1):
            w[n+1, j, k] = (2 * w[n, j, k] - w[n-1, j, k] +
                            dt**2 * ((w[n, j+1, k] - 2*w[n, j, k] + w[n, j-1, k]) / dx**2 +
                                     (w[n, j, k+1] - 2*w[n, j, k] + w[n, j, k-1]) / dy**2))

    # Apply Dirichlet boundary condition at x = 0
    w[n+1, 0, :] = 0

# Plotting the result at the final time step
plt.imshow(w[-1, :, :], extent=[0, max_x, 0, max_y])
plt.colorbar()
plt.title("Wave at t = {:.2f}".format(max_t))
plt.xlabel("x")
plt.ylabel("y")

# Save the figure
plt.savefig("Wave_simulation t=10.png")

plt.show()
