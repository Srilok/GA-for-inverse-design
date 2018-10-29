from dump import dump
import numpy as np 

def load_dump(filename,l1,l2): 

    #Read the initial populaiton file
    d = dump(filename)
    d.tselect.all() 
    ids = d.time()
    print ids
    [xMin,yMin,rMin,pMin] = l1 
    [xMax,yMax,rMax,pMax] = l2 

    start_pool = []
    note = []
    for i in range(len(ids)):
        ind = []
        x,y,r = d.vecs(ids[i],"x","y","radius")
        
        x = np.array(x) 
        y = np.array(y) 
        r = np.array(r)

        x = (x - xMin)/(xMax-xMin) 
        y = (y - yMin)/(yMax-yMin) 
        r = (r - rMin)/(rMax-rMin) 

#        r = (r-rMin)/(rMax-rMin)

        ndisks = len(x)
        for j in range(ndisks): 
            ind.append([x[j],y[j],r[j]])
        start_pool.append(ind)
        if ids[i] == 1: 
                note.append("Haogang Struct")
        elif ids[i] == 1344: 
                note.append("Good fit bad focus")
        elif ids[i] == 5962: 
                note.append("5962")
        elif ids[i] == 7242: 
                note.append("7242")
        else: 
                note.append("seeded_strcut")
    return start_pool, note

if __name__ == "__main__": 
    load_dump("init_population.dump", l1 = [0.0,0.0,110e-9], l2=[6.0e-6,310e-9,150e-9])
