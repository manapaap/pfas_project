# -*- coding: utf-8 -*-
"""
Plotting a 3D surface for the paper
"""

import numpy as np
import matplotlib.pyplot as plt


def fun(x, y):
    return x**2 + y**2


def bettercircle(radius=0.5):
    t = np.linspace(0, 2*np.pi, 100)
    xc = -2
    yc = 2
    r = radius

    x = r*np.cos(t) + xc
    y = r*np.sin(t) + yc
    z = fun(x, y)
    return x, y, z


def line(start=1.75, end=1.3):
    x = np.arange(-start, -end, 0.05)
    y = -x
    z = fun(x, y)
    return x, y, z


def arrowhead(start=1.6, end=1.25):
    x = np.arange(-start, -end, 0.05)
    y = np.zeros_like(x) + end

    y2 = np.arange(end, start, 0.05)
    x2 = np.zeros_like(y2) - end

    x_vals = np.concatenate((x, x2), axis=0)
    y_vals = np.concatenate((y, y2), axis=0)
    z_vals = fun(x_vals, y_vals)

    return x_vals, y_vals, z_vals


def main():
    # Plot the surface
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    # x = y = np.arange(-3.0, 3.0, 0.05)
    x = np.arange(-3, 2, 0.05)
    y = np.arange(-2, 3, 0.05)

    X, Y = np.meshgrid(x, y)
    zs = np.array(fun(np.ravel(X), np.ravel(Y)))
    Z = zs.reshape(X.shape)

    ax.plot_surface(X, Y, Z, cmap='viridis', alpha=0.72)

    # Add the circle
    x_cir, y_cir, z_cir = bettercircle(0.2)
    ax.plot3D(x_cir, y_cir, z_cir, 'red', antialiased=False, linewidth=2)

    # Add the arrow line
    x_lin, y_lin, z_lin = line(end=1.3)
    ax.plot3D(x_lin, y_lin, z_lin, 'red',
              antialiased=False, linewidth=1)

    # Add arrowhead
    x_hd, y_hd, z_hd = arrowhead(end=1.28, start=1.6)
    ax.plot3D(x_hd, y_hd, z_hd, 'red',
              antialiased=False, linewidth=1)

    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_zticklabels([])

    plt.show()


if __name__ == '__main__':
    main()
