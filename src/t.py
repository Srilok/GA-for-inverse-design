import numpy as np
from slope import obj_slope
ref= np.loadtxt('ideal.dat')
from scipy.interpolate import UnivariateSpline
M = UnivariateSpline(ref[:,0]*1e-6,ref[:,1],k=2,s=0.0)

a= np.loadtxt('/home/share/cnm50256/Disk_opt/1D/Lumerical_runs/phase/phase_31.dat')
b= np.loadtxt('/home/share/cnm50256/Disk_opt/1D/Lumerical_runs/phase/phase_24.dat')
ix = a[:,0]>=0.0
x = a[ix,0]
pred1 = a[ix,1]
pred1 -= pred1[0]
pred2 = b[ix,1]
pred2 -= pred2[0]
target = M(x)

print(10.*np.mean(np.abs(target-pred1)))
print(10.*np.mean(np.abs(target-pred2)))

print(obj_slope(target,pred1))
M1 = UnivariateSpline(x,pred1,k=2,s=1.5000)
print(obj_slope(target,pred2))
M2 = UnivariateSpline(x,pred2,k=2,s=1.5000)

#b= np.loadtxt('/home/share/cnm50256/Disk_opt/1D/phase/phase_2.dat')
#print(obj_slope(b[:,1],a[:,1]))
#print(ref.shape)
#print(a.shape)
#print(b.shape)
#print(obj_slope(ref[:,1],a[:,1]))
##print(obj_slope(ref[:,1],b[:,1]))

#from scipy.interpolate import UnivariateSpline
##data = np.loadtxt('/home/share/cnm50256/Disk_opt/1D/Lumerical_runs/phase/phase_12.dat')
#data = np.loadtxt('/home/share/cnm50256/Disk_opt/1D/ideal.dat')
#x = data[:,0]
#y = data[:,1]
#M = UnivariateSpline(x,y,k=2,s=1.5)
#
import matplotlib.pyplot as plt
plt.scatter(ref[:,0]*1e-6,ref[:,1],color='grey')
plt.scatter(x,target,color='orange')
plt.scatter(x,pred1,color='blue')
plt.scatter(x,M1(x),color='cyan')
plt.scatter(x,pred2,color='green')
plt.scatter(x,M2(x),color='yellow')
plt.xlim([0,6e-6])  
plt.ylim([-40,40])
plt.show()
