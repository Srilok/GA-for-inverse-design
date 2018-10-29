from scipy.interpolate import UnivariateSpline
from scipy.misc import derivative 
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt 
import numpy as np 
from math import pi, sqrt
from slope import obj_slope

from find_peaks import find_peaks
from sim_setup import run_lumerical

lam = 532e-9
foc = 10e-6
tpi = pi*2.0
shift = 0.0

#=========================================  
def refCylin(x):
    val = (tpi/lam) * (foc - sqrt(x*x + foc*foc))
#    while val > tpi:
#        val -= tpi
    #while val < -pi:
    #    val += tpi
#        shift += tpi
    return val
#==========================================
def refCylin_wrapped(x):
    val = (tpi/lam) * (foc - sqrt(x*x + foc*foc))
#    while val > tpi:
#        val -= tpi
    while val < -pi:
        val += tpi
#        shift += tpi
    return val
#==========================================  
def refSphere(x,y):
    val = (tpi/lam) * (foc - sqrt(x*x + y*y + foc*foc)) - pi
#    while val > tpi:
#        val -= tpi
    while val < -pi:
        val += tpi
    return val
#=========================================  
def testfunc():
    dx = 200.0/10000.0
    for x in range(10000):
        print x*dx, refCylin(x*dx)
#=========================================  
def set_target():
    from random import random
    outfile1 = open("ideal_wrap.dat","w")
    dx = 200.0/10000.0
    dx = dx*1e-6
    for x in range(10000):
        #print(' '.join(str(x) for x in [x*dx*1e6, refCylin(x*dx), "\n"]))
        outfile1.write(' '.join(str(x) for x in [x*dx*1e6, refCylin_wrapped(x*dx), "\n"]))
    outfile1.close()

    outfile2 = open("ideal_unwrap.dat","w") 
    for x in range(10000):
        #print(' '.join(str(x) for x in [x*dx*1e6, refCylin(x*dx), "\n"]))
        outfile2.write(' '.join(str(x) for x in [x*dx*1e6, refCylin(x*dx), "\n"]))
    outfile2.close()

    peaks = find_peaks("ideal_wrap.dat")
    return peaks 
#=========================================  


#---------------------  
def get_slope(x,phase, end_point,k = 2, s=0, plot_spline=False): 
        from scipy.interpolate import UnivariateSpline
        index = x<end_point
        
        fit_x = x[index]
        fit_y = phase[index]
        
        spline = UnivariateSpline(fit_x,fit_y,k=k,s=s)
        slope = spline.derivative()
        slope_x = slope(fit_x)

        
        #plot_spline = True
        
        if plot_spline: 
                
                plt.plot(fit_x,fit_y,'o')
                plt.plot(fit_x,spline(fit_x))
                plt.xlim([0,fit_x[-1]])
                plt.ylim([-4,4])
                plt.show() 
                plt.savefig('test_spline')
                plt.close()

        return slope_x
#---------------------

def analyze(individual,l1,l2,peaks):

    from math import fabs
    import numpy as np
    
    # Sometimes Lumerical crashes with a licensing error 
    # Try for say 5 times and if it still fails then assign 
    # very large value for objetive and 0 for transmission
    try: 
        run_dir,index = run_lumerical(l1,l2,individual)
        #Load Output from Lumerical
        infile = open(run_dir + 
                "/phase/phase_"+str(index)+".dat","r")
    except:

        err = 1000000000 
        trans = 0.0
        return (err,trans) 
 #       while True:
 #           print "WARNING:RERUNNING LUMERICAL.."
 #           run_dir,index = run_lumerical(l1,l2,individual) 

 #           
 #           try: 
 #               infile = open(run_dir + 
 #                   "/phase/phase_"+str(index)+".dat","r")
 #               break
 #           except: 
 #               continue
    #numfile = open("ndisk.dat", "r")
    #peakfile = open(run_dir+ "peak.dat","r")
    #ndisk = numfile.readline()
    #ndisk = float(ndisk)
    data = []
    for line in infile:
        col = line.split()
        if col[0] == "#": 
            trans = float(col[-1])
        #col = [float(x) for x in col]
        elif float(col[0]) >= -1.0e-10: 
            data.append([float(col[0]),float(col[1])])
        #data.append(col)

    xmin = data[0][0]
    ymin = data[0][1]
    for i, item in enumerate(data):
        data[i][0] -= xmin
        data[i][1] -= ymin
    data = np.array(data)
    
    peaks = np.array(peaks)
    #Get the slope of phase till the first peak
    #point_slope = get_slope(data[:,0],data[:,1],end_point=peaks[0]) 
    
    target = [] 
    predicted = []

    target_seg1 = []
    predicted_seg1 = []
    cnt = 0
    err1,err2,err3,err4,err5, err = 0.0,0.0,0.0,0.0,0.0,0.0
    #with open("test.data","w") as f: 
    
    #Fitting absolute values  
    for i, item in enumerate(data):
        refVal = refCylin(item[0])
        target.append(refVal)
        predicted.append(item[1])
        peak_flag = peaks - item[0]
        
        #Peak Values
        #if any(abs(peak_flag) <= 0.2e-6):
        #    err += 100*fabs(refVal-item[1])
            #    err += 100000*fabs(refVal - item[1])
        #if item[0] < 2.32e-6: 
        #    target.append(refVal)
        #    predicted.append(item[1])
            #if 0 < peak_flag[0] <= 0.2e-6:  #First peak fitting
            #    err += 100000*fabs(refVal-item[1])
            #else:
        if item[0] < 3.0e-6: 
            err1 += 10*fabs(refVal-item[1])

            target_seg1.append(refVal) 
            predicted_seg1.append(item[1])
        else:
            err1 += fabs(refVal-item[1])


    err1= 1000* err1/len(data)

    target = np.array(target) 
    predicted = np.array(predicted)

    #Overall Slope 
    #err2 += obj_slope(target,predicted,dlevel=1)
    err2 += obj_slope(target,predicted,k=2,s=0,dlevel=1)
    
    #Overall Curvature 
    #err3 += obj_slope(target,predicted,dlevel=2) 
    err3 += obj_slope(target,predicted,k=2,s=0,dlevel=2)
    
    #Initial Slope 
    err4 += obj_slope(target_seg1,predicted_seg1,k=2,s=0,dlevel=1)

    # Initial Curvature 
    err5 += obj_slope(target_seg1,predicted_seg1,dlevel=2)
    
    
    err = err1 + err2 + err3 + err4 +err5
    with open("obj_split.dat","a") as f: 
            f.write("Index: %d \t AbsDiff: %.5f \t Slope: %.5f \t Curvature: %.5f \t InitialSlope: %.5f \t InitialCurvature: %.5f \t Total: %.5f \n "%(individual.index, 
                    err1,err2, err3, err4, err5, err))
    return (err,trans)
#-----------------------------------------------------

#Calculates the gaussian function
def Gauss(x, a, x0, sigma):
        return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def analyzeFullSim(individual,l1,l2,peaks):

    from math import fabs
    import numpy as np
    
    [xMin,yMin,rMin,pMin] = l1
    [xMax,yMax,rMax,pMax] = l2

    # Sometimes Lumerical crashes with a licensing error 
    # Try for say 5 times and if it still fails then assign 
    # very large value for objetive and 0 for transmission
    try: 
        run_dir,index = run_lumerical(l1,l2,individual)
        #Load Output from Lumerical
        infile = open(run_dir + 
                "/phase/foc_"+str(index)+".dat","r")
    except:

        err = 10000000
        trans = 0.0
        raise
        #return (err,trans) 

    x = []
    focal_line = []
    read_intensity = False
    for line in infile:
        col = line.split()
        if col[0] == "#": 
            trans = float(col[2])
            err = -1*float(col[-1])
            read_intensity = True 
        else: 
            x.append(float(col[0]))
            focal_line.append(float(col[1]))
    
    x = np.array(x) 
    focal_line = np.array(focal_line)
    
    # Select only the centre portion of the curve to fit 
    # Minimizing the effect of secondary peaks in the focal line 
    I = np.abs(x) < ((xMax-xMin)/6.0)

    mean = 0
    sigma = np.sqrt(sum(focal_line[I] * (x[I] - mean)**2) / sum(focal_line[I]))
     
    #Curve Fitting
    popt,pcov = curve_fit(Gauss, x[I], focal_line[I], p0=[max(focal_line), mean, sigma])
    
    # For-debugging
    
    #print "Paramters: ", popt
    #print "Efficiency:",np.sqrt(2)*popt[0]*popt[2]*np.sqrt(np.pi)
    #
    a = abs(popt[0]) 
    c = min(abs(popt[2]),1.0e-6)
    #eff =  np.sqrt(2)*abs(popt[0])*abs(popt[2])*np.sqrt(np.pi)/(2*(xMax-xMin))
    eff =  np.sqrt(2)*a*c*np.sqrt(np.pi)/(2*(xMax-xMin))

    err = -1*eff
    
    if not read_intensity: 
        err = 10000000
        trans = 0.0
        return (err,trans) 
    else: 
        return (err,trans)

#---------------------------------------------------
#=========================================  
if __name__ == "__main__":
    analyze()


