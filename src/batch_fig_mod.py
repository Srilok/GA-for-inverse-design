from sys import argv 
import numpy as np 
import os 
from os.path import isfile,isdir, join 
import matplotlib.pyplot as plt
import errno


def batch_plot_fig(foldername): 

        target = np.loadtxt("ideal.dat")
        
        only_folders = [ f for f in os.listdir(foldername) if isfile(join(foldername,f))]
    
        fig_folder = "./FitFig_1D_Try3" 
        try:
                os.makedirs(fig_folder)
        except OSError as exception:
                if exception.errno != errno.EEXIST:
                        raise
        for folder in only_folders:
                try:
                    if isfile(foldername +"/"+ folder):
                    #if isfile(join(foldername + folder,"phase.dat")):
                        print folder 
                        test  = np.loadtxt(foldername + "/" +folder)
                        #test = np.loadtxt(foldername + folder + "/phase.dat")
                        ix = test[:,0] >=0 
                        plot_test_x = test[ix,0]*1e+6
                        plot_test_phase = test[ix,1] 
                        plot_test_phase -= plot_test_phase[0]
                        #base_val = test[0,0]
                        #test[:,0] = test[:,0] - base_val 
                        #test[:,0] = test[:,0]*1e+6
                        fig = plt.figure() 
                        ax = fig.add_subplot(211)

                        plt.plot(plot_test_x,plot_test_phase,'o',label = folder) 
                        plt.plot(target[:,0],target[:,1],linewidth = 2) 
                        plt.xlim([0,6])
                        plt.ylim([-20,2])
                        plt.xlabel("x (um)") 
                        plt.ylabel("Phase")
                        ax1 = fig.add_subplot(212)


                        ix = plot_test_phase < -3.14
                        while True in ix:
                            ix = plot_test_phase < -3.14
                            plot_test_phase[ix] += 2*3.14
#                           print len(ix)

                        ix = target[:,1] < -3.14
                        while True in ix:
                            ix = target[:,1] < -3.14
                            target[ix,1] += 2*3.14
                        
                        plt.plot(plot_test_x,plot_test_phase,'o') 
                        plt.plot(target[:,0],target[:,1]) 

                        plt.xlim([0,6])
                        plt.ylim([-4,4])


                        plt.plot(plot_test_x,plot_test_phase)

                        plt.savefig(fig_folder + "/" + folder+ ".png")
                        plt.close()
                except: 
                        raise
                        continue 



if __name__ == "__main__": 
        foldername = argv[1]
        batch_plot_fig(foldername)
