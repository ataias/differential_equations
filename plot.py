#!/usr/bin/env python3

import numpy as np
import struct
import os
import argparse
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

def read_matrix(filename):
    f = open(filename, 'rb')
    n = struct.unpack('i', f.read(4))[0]
    print("Started reading file")
    mesh = np.zeros((n,n), dtype=np.double)
    for i in range(0, n):
        for j in range(0, n):
            mesh[i,j] = struct.unpack('d',f.read(8))[0]
    f.close()
    print("Finished reading file")
    return n, mesh

def plot(n, mesh):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    X = np.linspace(0, 1, n)
    Y = np.linspace(0, 1, n)
    X, Y = np.meshgrid(X, Y)
    R = np.sqrt(X**2 + Y**2)
    Z = mesh
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,
                           linewidth=0, antialiased=False)

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()

def main():
    parser = argparse.ArgumentParser(description='Read 2D array of doubles and plot it considering x in (0,1) and y in (0,1).')
    parser.add_argument('--filename', dest='filename', type=str,
                    help='the binary file that is going to be used', required=True)
    args = parser.parse_args()
    n, mesh = read_matrix(args.filename)

    plot(n, mesh)

if __name__ == '__main__':
    main()
