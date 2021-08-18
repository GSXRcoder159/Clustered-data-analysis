# -*-coding: utf-8-*-

# THIS PROGRAM IS USED TO PLOT THE POINTS IN A CONVEX HULL OF A CLUSTER ONTO A SCATTER PLOT


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


file = open('ABSOLUTE PATH TO THE SOURCE FILE', 'r')

number = ""

x = []
y = []
a = []
b = []

for line in file:
    number = ""
    for character in line:
        if character == ' ':
            a.append(int(float(number)))
            number = ""
        elif character == '#':
            x.append(np.array(a))
            y.append(np.array(b))
            a.clear()
            b.clear()
        elif character == '\n' and number != '':
            b.append(int(float(number)))
            number = ""
        else:
            number += character

#GET THE X-COORDINATES AND Y-COORDINATES OF THE POINTS ON THE CONVEX HULL OF A CLUSTER
#n - the index of the cluster in l1_clusters
def get_cluster_coordinates(n):
    file = open(
        'C:\\Users\\adamk\\OneDrive\\Plocha\\Clusers_dimensions.txt', 'r')

    number = ""
    lines = file.readlines()

    x = []
    y = []

    number = ""
    for character in lines[n]:
        if character == ' ':
            x.append(int(number))
            number = ""
        elif character == '\t':
            y.append(int(number))
            number = ""
        else:
            number += character

    a = np.array(x)
    b = np.array(y)

    return a, b

#PLOTS THE ACQUIRED DATA
def plot(n):
    x0 = x[n]
    y0 = y[n]

    x1, y1 = get_cluster_coordinates(n)

    plt.subplot(1, 2, 1)

    plt.scatter(x0, y0)

    plt.title('Border plot: ' + str(n))
    plt.xlabel('X - value')
    plt.ylabel('Y - value')

    plt.subplot(1, 2, 2)
    plt.scatter(x1, y1)

    plt.title('Cluster plot: ' + str(n))
    plt.xlabel('X - value')
    plt.ylabel('Y - value')

    plt.show()

#PLOT THE 20TH CLUSTER (INDEX 19)
plot(8)
