import numpy as np

# np.zeros(10),10**(np.arange(10)/3)
sample = np.random.random((100000, 2))
z = (sample[:,0]+sample[:,1])/2
z = np.reshape(z,(-1,1))
print(z)
sample = np.hstack((sample,z))
distance = 0
for i in range(len(sample) - 1):
	distance += np.sqrt(np.sum((sample[i] - sample[i + 1]) ** 2))
print(distance / (len(sample) - 1))