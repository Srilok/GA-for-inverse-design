from sys import argv


def find_peaks(filename): 
        x =[] 
        phase = []
        peaks = []
        rfile = open(filename,'r')
        for line in rfile: 
                col = line.split() 
                x.append(float(col[0])*1.0e-6)
                phase.append(float(col[1]))

        for i in range(len(x)-1):

                if abs(phase[i+1] - phase[i]) >= 3.12:
                        peaks.append(x[i])

        rfile.close()
        outfile = open("peak.dat","w") 

        for peak in peaks: 
                outfile.write(str(peak) + "\n")

        outfile.close()
        return peaks 



if __name__ =="__main__":

        filename = argv[1]
        find_peaks(filename)
