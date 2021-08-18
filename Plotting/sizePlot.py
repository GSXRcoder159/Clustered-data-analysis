# -*-coding: utf-8-*-

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# THIS PROGRAM IS USED TO PLOT THE SIZES OF CLUSTERS INTO A HISTOGRAM

file = open('ABSOLUTE PATH TO THE SOURCE FILE', 'r')

number = ""

size = []
time = []

for line in file:
    number = ""
    for character in line:
        if character == ' ':
            size.append(int(number))
            number = ""
        elif character == '\n':
            time.append(float(number))
            number = ""
        else:
            number = number + character

size = np.array(size)
time = np.array(time)

# size

n, bins, patches = plt.hist(
    x=size, bins=int(size.size / 12), color='#8c73ff', alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel("Size [in pixels]")
plt.ylabel("Count")
plt.title("Size")
maxfreq = n.max()
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)

plt.show()
