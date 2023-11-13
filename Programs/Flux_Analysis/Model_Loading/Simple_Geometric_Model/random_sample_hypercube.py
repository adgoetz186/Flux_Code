import numpy as np

#np.zeros(10),10**(np.arange(10)/3)
sample = np.random.random((1000000,10))
for i in range(np.shape(sample)[1]):
	print(10**(i/3))
	sample[:,i] *= 10**(i/3)
	
distance = 0
for i in range(len(sample)-1):
	distance += np.sqrt(np.sum((sample[i]-sample[i+1])**2))
print(distance/(len(sample)-1))