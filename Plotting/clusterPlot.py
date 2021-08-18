
# THIS PROGRAM IS USED TO PLOT PIXELS FROM A GIVEN CLUSTER ONTO A SCATTER PLOT

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


file = open('ABSOLUTE PATH TO THE SOURCE FILE', 'r')

number = ""
lines = file.readlines()

# LOAD A FILE AND FROM IT THE X-COORDINATES AND Y-COORDINATES OF PIXELS FROM A GIVEN CLUSTER
# n - the index of the cluster in l1_clusters
def get_cluster_coordinates(n):
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
            number = number + character

    a = np.array(x)
    b = np.array(y)

    return a, b

#PLOTS THE ACQUIRED DATA
def plot_cluster(n):
    x = np.array([])
    y = np.array([])

    x, y = get_cluster_coordinates(n)

    plt.scatter(x, y)
    plt.title('Cluster plot: ' + str(n))
    plt.xlabel('X - value')
    plt.ylabel('Y - value')

    plt.show()

#PLOT THE 20TH CLUSTER (INDEX 19)
plot_cluster(19)
