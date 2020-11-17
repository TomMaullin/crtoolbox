import numpy as np
import time

a = np.random.randn(100,500,2)

b1 = a[:,:,0]
b2 = a[:,:,1]

t1 = time.time()
a_std = np.std(a,axis=0,ddof=1)
t2 = time.time()
print(t2-t1)


t1 = time.time()
b1_std = np.std(b1,axis=0,ddof=1)
b2_std = np.std(b2,axis=0,ddof=1)
t2 = time.time()
print(t2-t1)

print(a_std.shape)
print(b1_std.shape)
print(b2_std.shape)
