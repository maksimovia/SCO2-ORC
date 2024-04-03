import pandas as pd
import prop
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt


df = pd.read_csv('comp790.txt', sep='\t',header=None)
#print(df)
N = 200
M = 20
pi = np.linspace(1.1,7,N)
pk = np.linspace(7.5,8.5,M)
A = np.zeros((N,M),dtype='float32')
for i in range(0,N):
    for j in range(0, M):
        A[i,j] = df.iloc[j+i*M,-1]
#print(df.iloc[0:20,0:4])
#print(A[0,:])
# plt.contourf(pk,pi,A,cmap="Blues",levels=25)
# plt.show()
B = np.zeros(N)
for i in range(0,N):
    B[i]=max(A[i,:])
plt.plot(pi,B)
plt.show()