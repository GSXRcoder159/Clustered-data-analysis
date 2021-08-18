# -*-coding: utf-8-*-

# THIS PROGRAM IS USED TO PLOT DIFFERENCES IN SIZES BETWEEN CLUSTER PAIRS INTO A HISTOGRAM

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


file = open('ABSOLUTE PATH TO THE SOURCE FILE', 'r')

number = ""
n = 0

l1 = [] #LIST OF THE INDICES OF CLUSTERS FROM LAYER 1
l2 = [] #LIST OF THE INDICES OF CLUSTERS FROM LAYER 2
dsize = []  #THE DIFFERENCE IN SIZE OF THE TWO CLUSTERS

for line in file:
    number = ""
    for character in line:
        if character == ' ' and n == 0:
            l1.append(int(number))
            number = ""
            n = 1
        elif character == ' ' and n == 1:
            l2.append(int(number))
            number = ""
            n = 0
        elif character == '\n':
            dsize.append(int(number))
            number = ""
        else:
            number += character

l1 = np.array(l1)
l2 = np.array(l2)
dsize = np.array(dsize)


# Density Plot and Histogram of all arrival delays
sns.distplot(dsize, hist=False, kde=True,
             kde_kws={'linewidth': 4})

plt.title('Differences in cluster sizes')
plt.xlabel('Difference in size (px)')
plt.ylabel('Density')

plt.show()
