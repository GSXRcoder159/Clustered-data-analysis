import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# ???

file = open('ABSOLUTE PATH TO THE SOURCE FILE', 'r')

number = ""

x_data = []
y_data = []

for line in file:
    number = ""
    for character in line:
        if character == ' ':
            y_data.append(int(number))
            number = ""
        elif character == '\n':
            x_data.append(int(number))
            number = ""
        else:
            number = number + character

x_data = np.array(x_data)
y_data = np.array(y_data)
# print(x_data)

first_edge, last_edge = x_data.min(), x_data.max()

n_equal_bins = 10
bin_edges = np.linspace(start=first_edge, stop=last_edge,
                        num=n_equal_bins + 1, endpoint=True)

bcounts = np.bincount(x_data)
hist, _ = np.histogram(x_data, range=(0, x_data.max()), bins=x_data.max() + 1)

np.array_equal(hist, bcounts)

dict(zip(np.unique(x_data), bcounts[bcounts.nonzero()]))


n, bins, patches = plt.hist(
    x=x_data, bins=100, color='#8c73bb', alpha=0.7, rwidth=0.85)
plt.grid(axis='y', alpha=0.75)
plt.xlabel("Size [in pixels]")
plt.ylabel("Count")
plt.title("")
maxfreq = n.max()
plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)

# # plt.show()

# sns.set_style('darkgrid')
# sns.displot(size)
plt.show()
