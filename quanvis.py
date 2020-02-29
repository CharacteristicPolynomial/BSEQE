import numpy as np
import matplotlib.pyplot as plt

q1 = np.loadtxt("quant_2.list")

num_bins = 100
print(np.mean(q1), np.std(q1))
plt.figure()
plt.hist(q1, num_bins, facecolor='red', alpha=0.3)
plt.show()