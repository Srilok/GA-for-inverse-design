from sys import argv 
import numpy as np 
import os 
from os.path import isfile,isdir, join 
import matplotlib.pyplot as plt
import errno


def batch_plot_fig(foldername): 

        target = np.loadtxt("ideal.dat")
        
        only_folders = [ f for f in os.listdir(foldername) if isfile(join(foldername,f))]
    
        fig_folder = "./FitFig_1D" 
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
                        base_val = test[0,0]
                        test[:,0] = test[:,0] - base_val 
                        test[:,0] = test[:,0]*1e+6
                        plt.plot(test[:,0],test[:,1],'o',label = folder) 
                        plt.plot(target[:,0],target[:,1],linewidth = 2) 
                        plt.xlim([0,6])
                        plt.xlabel("x (um)") 
                        plt.ylabel("Phase")
                        plt.savefig(fig_folder + "/" + folder+ ".png")
                        plt.close()
                except: 
                        raise
                        continue 



if __name__ == "__main__": 
        foldername = argv[1]
        batch_plot_fig(foldername)
