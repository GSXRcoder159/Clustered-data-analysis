import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# THIS PROGRAM IS USED TO PLOT ENERGIES OF CLUSTERS INTO A HISTOGRAM

file = open('ABSOLUTE PATH TO THE SOURCE FILE', 'r')

number = ""

energy = []

for line in file:
    number = ""
    for character in line:
        if character == '\n':
            energy.append(float(number))
            number = ""
        else:
            number = number + character

size = np.array(energy)

# size

n, bins, patches = plt.hist(
    x=size, bins=size.size, color='#8c73ff', alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel("Total energy [in keV]")
plt.ylabel("Count")
plt.title("Total energy")
maxfreq = n.max()
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)


plt.show()
