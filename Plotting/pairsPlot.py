import numpy as np
import matplotlib.pyplot as plt
from numpy.core.defchararray import index
import pandas as pd
import seaborn as sns
import decimal as decimal
import matplotlib.cm as cm
import math as math
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator


# THIS PROGRAM IS USED TO PLOT CLUSTERS AND PAIRS OF CLUSTERS WITH THEIR ENERGIES ONTO AN "IMSHOW"


#LOAD THE X-COORDINATE, Y-COORDINATE AND TIME OF ALL PIXELS OF ALL CLUSTERS
def get_cluster_coordinates():
    x = []
    y = []
    t = []
    a1 = []
    b1 = []
    c1 = []
    a2 = []
    b2 = []
    c2 = []
    i = 0

    file = open('THE ABSOLUTE PATH TO THE SOURCE FILE', 'r')
    number = ""

    for line in file:
        number = ""
        for character in line:
            if character == ' ':
                if(i == 0):
                    x.append(int(number))
                    number = ""
                    i = 1
                elif(i == 1):
                    y.append(int(number))
                    number = ""
                    i = 2
            elif character == '\n':
                if number != '':
                    t.append(float(number))
                    number = ""
                i = 0
            elif character == '$':
                a = np.array(x)
                a1.append(a)
                b = np.array(y)
                b1.append(b)
                c = np.array(t)
                c1.append(c)
                x.clear()
                y.clear()
                t.clear()
            elif character == '#':
                a = np.array(x)
                a2.append(a)
                b = np.array(y)
                b2.append(b)
                c = np.array(t)
                c2.append(c)
                x.clear()
                y.clear()
                t.clear()
            else:
                number += character

    return a1, b1, c1, a2, b2, c2

#PLOT A SINGLE CLUSTER FROM LAYER 1
# n - the index of the cluster in l1_clusters
def orig_plot(n):
    a1 = []
    b1 = []
    c1 = []
    a2 = []
    b2 = []
    c2 = []

    a1, b1, c1, a2, b2, c2 = get_cluster_coordinates()

    a_1 = []
    a_1 = a1[n].tolist()
    b_1 = []
    b_1 = b1[n].tolist()
    c_1 = []
    c_1 = c1[n].tolist()

    a_2 = []
    a_2 = a2[n].tolist()
    b_2 = []
    b_2 = b2[n].tolist()
    c_2 = []
    c_2 = c2[n].tolist()

    cols = [[0]*255 for _ in range(255)]
    cols1 = [[0]*255 for _ in range(255)]
    cols3 = [[0]*255 for _ in range(255)]

    length = len(a_1)
    temp = np.array(c_1)
    m = temp.max()

    length1 = len(a_2)
    temp1 = np.array(c_2)
    m1 = temp1.max()

    for index in range(length):
        c_1[index] = math.pow((c_1[index] / m), 10000000000)
    for i in range(length):
        cols[b_1[i]][a_1[i]] = c_1[i]

    for index in range(length1):
        c_2[index] = math.pow((c_2[index] / m1), 10000000000)
    for i in range(length1):
        cols1[b_2[i]][a_2[i]] = c_2[i]

    for i in range(len(cols)):
        for j in range(len(cols[i])):
            cols3[i][j] = ((cols[i][j]/2) + (cols[i][j]/2))

    t1 = np.array(cols)
    t2 = np.array(cols1)
    t = np.array(cols3)

    np.set_printoptions(threshold=np.inf)

    # START
    plt.subplot(211)
    nrows, ncols = 255, 255
    grid = t1.reshape((ncols, nrows))
    plt.imshow(grid,
               interpolation='nearest', cmap=cm.hot, origin='lower')
    plt.xlim([0, 255])
    plt.ylim([0, 255])
    plt.title("Layer 1")

    plt.subplot(212)
    nrows, ncols = 255, 255
    grid = t2.reshape((ncols, nrows))
    plt.imshow(grid,
               interpolation='nearest', cmap=cm.hot, origin='lower')
    plt.xlim([0, 255])
    plt.ylim([0, 255])
    plt.title("Layer 2")

    plt.subplots_adjust(bottom=0.05, right=0.8, top=0.95)
    cax = plt.axes([0.85, 0.1, 0.075, 0.8])
    plt.colorbar(cax=cax)


#MAKE TWO PLOTS ONE OF A CLUSTER FORM THE LAYER 1 AND SECOND OF THE CORESPONDING PAIR FORM LAYER 2
# n - the index of the cluster in l1_clusters
def double_plot(n):
    a1 = []
    b1 = []
    c1 = []
    a2 = []
    b2 = []
    c2 = []

    a1, b1, c1, a2, b2, c2 = get_cluster_coordinates()

    a_1 = []
    a_1 = a1[n].tolist()
    b_1 = []
    b_1 = b1[n].tolist()
    c_1 = []
    c_1 = c1[n].tolist()

    a_2 = []
    a_2 = a2[n].tolist()
    b_2 = []
    b_2 = b2[n].tolist()
    c_2 = []
    c_2 = c2[n].tolist()

    cols = [[0]*255 for _ in range(255)]
    cols1 = [[0]*255 for _ in range(255)]
    cols3 = [[0]*255 for _ in range(255)]

    length = len(a_1)
    temp = np.array(c_1)
    m = temp.max()

    length1 = len(a_2)
    temp1 = np.array(c_2)
    m1 = temp1.max()

    for index in range(length):
        c_1[index] = math.pow((c_1[index] / m), 10000000000)
    for i in range(length):
        cols[b_1[i]][a_1[i]] = c_1[i]

    for index in range(length1):
        c_2[index] = math.pow((c_2[index] / m1), 10000000000)
    for i in range(length1):
        cols1[b_2[i]][a_2[i]] = c_2[i]

    for i in range(len(cols)):
        for j in range(len(cols[i])):
            cols3[i][j] = ((cols[i][j]/2) + (cols[i][j]/2))

    t1 = np.array(cols)
    t2 = np.array(cols1)
    t = np.array(cols3)

    np.set_printoptions(threshold=np.inf)

    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.suptitle('Cluster pairs ploted', fontsize=20)

    # Plot number 1
    ax1.set_title('Layer 1')
    nrows, ncols = 255, 255
    grid = t1.reshape((ncols, nrows))
    im1 = ax1.imshow(grid,
                     interpolation='nearest', cmap=cm.hot, origin='lower', aspect='auto')
    plt.xlim([0, 255])
    plt.ylim([0, 255])

    # Plot number 2
    ax2.set_title('Layer 2')
    nrows, ncols = 255, 255
    grid = t2.reshape((ncols, nrows))
    im2 = ax2.imshow(grid,
                     interpolation='nearest', cmap=cm.hot, origin='lower', aspect='auto')
    plt.xlim([0, 255])
    plt.ylim([0, 255])
    divider2 = make_axes_locatable(ax2)
    cax2 = divider2.append_axes("right", size="20%", pad=0.05)
    cbar2 = plt.colorbar(
        im2, cax=cax2)

    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.show()

#PLOT A CLUSTER FROM LAYER 1 AND ITS CORESPONDING PAIR FROM LAYER 2
# n - the index of the cluster in l1_clusters
def both_plot(n):
    a1 = []
    b1 = []
    c1 = []
    a2 = []
    b2 = []
    c2 = []

    a1, b1, c1, a2, b2, c2 = get_cluster_coordinates()

    a_1 = []
    a_1 = a1[n].tolist()
    b_1 = []
    b_1 = b1[n].tolist()
    c_1 = []
    c_1 = c1[n].tolist()

    a_2 = []
    a_2 = a2[n].tolist()
    b_2 = []
    b_2 = b2[n].tolist()
    c_2 = []
    c_2 = c2[n].tolist()

    cols = [[0]*255 for _ in range(255)]
    cols1 = [[0]*255 for _ in range(255)]
    cols3 = [[0]*255 for _ in range(255)]

    length = len(a_1)
    temp = np.array(c_1)
    m = temp.max()

    length1 = len(a_2)
    temp1 = np.array(c_2)
    m1 = temp1.max()

    for index in range(length):
        c_1[index] = math.pow((c_1[index] / m), 10000000000)
    for i in range(length):
        cols[b_1[i]][a_1[i]] = c_1[i]

    for index in range(length1):
        c_2[index] = math.pow((c_2[index] / m1), 10000000000)
    for i in range(length1):
        cols1[b_2[i]][a_2[i]] = c_2[i]

    for i in range(len(cols)):
        for j in range(len(cols[i])):
            cols3[i][j] = ((cols[i][j]/2) + (cols1[i][j]/2))

    t1 = np.array(cols)
    t2 = np.array(cols1)
    t = np.array(cols3)

    np.set_printoptions(threshold=np.inf)

    plt.title('Both layers')
    nrows, ncols = 255, 255
    grid = t.reshape((ncols, nrows))
    plt.imshow(grid,
               interpolation='nearest', cmap=cm.hot, origin='lower', aspect='auto')
    plt.xlim([0, 255])
    plt.ylim([0, 255])
    plt.colorbar()

    plt.tight_layout()
    plt.subplots_adjust(top=0.85)
    plt.show()


possible_pairs = {0, 1, 2, 3, 4, 5}

for items in possible_pairs:
    both_plot(items)

# orig_plot(n)
# double_plot(n)
# both_plot(n)
