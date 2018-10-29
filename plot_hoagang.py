from sys import argv 
import numpy as np 
import matplotlib.pyplot as plt 
from slope import obj_slope 
from CalcObjective import refCylin

def compare_target(filename): 
        test = np.loadtxt(filename)
        target = np.loadtxt("ideal.dat")
        
        ix = test[:,0] >= 0.0 
        
        x_data = test[ix,0]*1e+6 
        y_data = test[ix,1] 

        target_obj = [refCylin(x) for x in test[ix,0]]
        target_obj = np.array(target_obj) 

        obj = 0.0
        obj += 10*np.mean(np.abs(target_obj-y_data))
        
#        obj += obj_slope(target_obj, y_data)

        print("Objective: %.2f"%obj)
        
        fig = plt.figure()
        ax = fig.add_subplot(211)

        
        plt.plot(x_data,y_data,'o') 
        plt.plot(target[:,0],target[:,1]+y_data[0]) 

        plt.xlim([0,6])
        plt.ylim([-20,20])

        ax1 = fig.add_subplot(212)

        target[:,1] +=3.14
        y_data -=y_data[0] + 3.14
        #print(y_data)
        #print(len(y_data))
        ix = y_data < -3.14
        while True in ix:
                ix = y_data < -3.14
                y_data[ix] += 2*3.14
#                print len(ix)

        ix = target[:,1] < -3.14
        while True in ix:
                ix = target[:,1] < -3.14
                target[ix,1] += 2*3.14

        ix = y_data > 3.14
        while True in ix:
                ix = y_data > 3.14
                y_data[ix] -= -2*3.14
                print len(ix)

        ix = target[:,1] > 3.14
        while True in ix:
                ix = target[:,1] > 3.14
                target[ix,1] -= 2*3.14



#        plt.plot(x_data,y_data,'o') 
        plt.plot(target[:,0],target[:,1]) 

        plt.xlim([0,6])
        plt.ylim([-4,4])


        plt.plot(x_data,y_data)
        plt.show() 

if __name__ == "__main__": 
        filename = argv[1]
        compare_target(filename)
