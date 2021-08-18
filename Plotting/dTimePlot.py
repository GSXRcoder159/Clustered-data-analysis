import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# THIS PROGRAM IS USED TIME DIFFERENCES OF A CLUSTERS INTO A HISTOGRAM

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

first_edge, last_edge = size.min(), size.max()

n_equal_bins = 10
bin_edges = np.linspace(start=first_edge, stop=last_edge,
                        num=n_equal_bins + 1, endpoint=True)

bcounts = np.bincount(size)
hist, _ = np.histogram(size, range=(0, size.max()), bins=size.max() + 1)

np.array_equal(hist, bcounts)

dict(zip(np.unique(size), bcounts[bcounts.nonzero()]))

n, bins, patches = plt.hist(
    x=time, bins=200, color='#8c73bb', alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel("delta time [in ns]")
plt.ylabel("Count")
plt.title("Delta time")
maxfreq = n.max()
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)


plt.show()
