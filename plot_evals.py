
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math

from loadW import loadW

N=3000
p=50/N

W45=loadW(N,p,7)
D45 = np.linalg.eigvals(W45)

W200=loadW(N,p,8)
D200 = np.linalg.eigvals(W200)

W1000=loadW(N,p,9)
D1000 = np.linalg.eigvals(W1000)

matplotlib.rcParams.update({'font.size': 22})
plt.rc('xtick', labelsize=32)
plt.rc('ytick', labelsize=32)


plt.figure()
plt.scatter(np.real(D45), np.imag(D45), color="blue")
plt.hold(True)
plt.scatter(np.real(D200), np.imag(D200), color="green")

plt.scatter(np.real(D1000), np.imag(D1000), color="red")

plt.legend(["L=45", "L=200", "L=1000"])
plt.xlabel('$Re(\lambda)$', fontsize=32)
plt.ylabel('$Im(\lambda)$', fontsize=32)
plt.tight_layout()



