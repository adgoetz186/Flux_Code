import numpy as np

sample_points = np.loadtxt("524288_sec",delimiter=",")
distance = 0
for i in range(len(sample_points)-1):
	distance += np.sqrt(np.sum((sample_points[i]-sample_points[i+1])**2))
print(distance/(len(sample_points)-1))
print(np.shape(sample_points))
