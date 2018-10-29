from sys import argv 
import numpy as np 
import matplotlib.pyplot as plt 
from slope import obj_slope 
from CalcObjective import refCylin

def compare_target(argv): 
        nfiles = len(argv) -1
        for i in range(nfiles): 
              test = np.loadtxt(argv[i+1])
              target = np.loadtxt("ideal_unwrap.dat")
              
              ix = test[:,0] >= 0.0 
              
              x_data = test[ix,0]*1e+6 
              y_data = test[ix,1] 
              y_data -= y_data[0]  
              
              # base_val = test[0,0]
              # test[:,0] = test[:,0] - base_val #0.5e-6
              # test[:,0] = test[:,0]*1e+6

              # base_val = test[0,1] 
              # test[:,1] = test[:,1] - base_val 

              target_obj = [refCylin(x) for x in test[ix,0]]
              target_obj = np.array(target_obj) 

              obj = 0.0
              obj += 10*np.mean(np.abs(target_obj-y_data))
              
              #obj += obj_slope(target_obj, y_data)

              print("Objective: %.2f"%obj)
              
              fig = plt.figure()
              ax = fig.add_subplot(211)

              
              plt.plot(x_data,y_data,'o') 
              plt.plot(target[:,0],target[:,1]) 

              plt.xlim([0,6])
              plt.ylim([-20,4])

              ax1 = fig.add_subplot(212)

              #print(y_data)
              #print(len(y_data))
              ix = y_data < -3.14
              while True in ix:
                ix = y_data < -3.14
                y_data[ix] += 2*3.14
        #        print len(ix)

              ix = target[:,1] < -3.14
              while True in ix:
                ix = target[:,1] < -3.14
                target[ix,1] += 2*3.14

              ix = y_data > 3.14
              # while True in ix:
              #         ix = y_data > 3.14
              #         y_data[ix] -= -2*3.14
        #     #          print len(ix)

              # ix = target[:,1] > 3.14
              # while True in ix:
              #         ix = target[:,1] > 3.14
              #         target[ix,1] -= 2*3.14


         #    for i in range(len(y_data)):
         #    if y_data[i] < -3.14:
         #       shift = -3.14 - y_data[i]
         #       y_data[i] += 2*3.14 + shift 
         #    if y_data[i] > 3.14: 
         #       shift = y_data -3.14
         #       y_data[i] -= 2*3.14 

         #    for i in range(len(target)):
         #    if target[i,1] < -3.14:
         #       shift = -3.14 - target[i,1] 
         #       target[i,1] += 2*3.14 + shift  
         #    if target[i,1] > 3.14: 
         #       shift = target[i,1] -3.14
         #       target[i,1] -= 2*3.14 +shift  

        #     plt.plot(x_data,y_data,'o') 
              plt.plot(target[:,0],target[:,1]) 

              plt.xlim([0,6])
              plt.ylim([-4,4])


              plt.plot(x_data,y_data)
        plt.show() 

if __name__ == "__main__": 
        compare_target(argv)
