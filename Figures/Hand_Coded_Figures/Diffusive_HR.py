import matplotlib.pyplot as plt
import numpy as np

time_between_points = [1,2,4,8,16,32,64,128,256,512,1024]
average_distance = np.array([79,9.4,152,205,243,288.9,360,346,387,389.9,411.9])
true_random = 411

plt.hlines(1,1,1024,color = "black", linestyles="dashed",label = "Truly random points")
plt.scatter(time_between_points,average_distance/true_random,label = "Unconstrained 10-D box, hit and run (CD)")



time_between_points = [1,2,4,8,16,32,64,128,256,512,1024]
average_distance = np.array([0.3428438201278746,0.4597912091795292,0.5281095773121651,0.5492830022118915,0.5949885185718863,0.5781098578563868,0.5777820847347932,0.558967504502721,0.5815214118245637, 0.5728365408477571,0.5748064554450056])
true_random = 0.5823775339799527

#plt.hlines(1,1,512,color = "orange",label = "Truly random points")
plt.scatter(time_between_points,average_distance/true_random,label = "2-D plane in 3-D box, hit and run (CD)")
plt.xlabel("Iterations between saved points")
plt.ylabel("Distance between saved points")
plt.legend()
plt.show()

time_between_points = [1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,118202]
average_distance = [741.5618085,961.5147288,1526.756,2106.708618,3069.560047,4370.757162,5914.104576,8647.69611,12261.88455,16882.27467,23570.56496,32037.75822,43324.62799,57525.29482,74509.79455,93954,109221.2782,118202.8384]


plt.scatter(time_between_points,average_distance,label = "Recon2.2 model, hit and run (CD)")
plt.xlabel("Iterations between saved points")
plt.ylabel("Distance between saved points")
plt.legend()
plt.show()