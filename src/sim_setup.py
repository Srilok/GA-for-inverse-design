from ScriptTemplate import InputTemplate
import os
import errno
import numpy as np
from itertools import combinations
from SimulationClassDef import Simulation
import subprocess32 as subproc
from math import fabs,copysign

#------------------------------------------------------------
def get_box_input(input_file): 

        
    with open(input_file) as f:
        
        line = f.readline(); 
        col = line.split() 
        nMin = int(col[0])
        nMax = int(col[1]) 
        
        line = f.readline(); 
        col = line.split() 
        pMin = float(col[0])
        pMax = float(col[1]) 

        line = f.readline(); 
        col = line.split() 
        xMin = float(col[0])
        xMax = float(col[1]) 

    
        line = f.readline(); 
        col = line.split() 
        yMin = float(col[0])
        yMax = float(col[1]) 
        
        line = f.readline(); 
        col = line.split() 
        rMin = float(col[0])
        rMax = float(col[1]) 

    return [nMin,nMax,pMin,pMax ,xMin,xMax,yMin,yMax,rMin,rMax]
#----------------------------------------------------------
def make_lumerical_input(index,individual,parent_dir,run_dir,phase_dir,l1,l2,mesh_override=False): 
    
    mesh_override = True 

    [xMin,yMin,rMin,pMin] = l1 
    [xMax,yMax,rMax,pMax] = l2 
    
    disktemp = InputTemplate()
    disktemp.readTemplate(parent_dir + "/disktemp.txt")
    diskStr = ""

    meshtemp = InputTemplate()
    meshtemp.readTemplate(parent_dir +"/meshtemp.txt") 
    meshStr= ""

    #posttemp = InputTemplate() 
    #posttemp.readTemplate(parent_dir + "/post.lsf") 
    #posttemp.clearReplacement()
    #posttemp.addReplacement('write("' + phase_dir + '/phase_' + 
    #                str(index) + '.dat",str);')
    #posttemp.addReplacement('write("' + phase_dir + '/phase_' + 
    #                str(index) + '.dat",str);')
    #posttemp.generateNewFile(run_dir+"/post.lsf")
   

    posttemp = InputTemplate() 
    posttemp.readTemplate(parent_dir + "/get_focaline.lsf") 
    posttemp.clearReplacement()
    posttemp.addReplacement('write("' + phase_dir + '/foc_' + 
                    str(index) + '.dat",str);')
    posttemp.addReplacement('write("' + phase_dir + '/foc_' + 
                    str(index) + '.dat",str);')
    posttemp.generateNewFile(run_dir+"/get_focaline.lsf")
    
    coords =[] 
    # Translate the reduced cell coordinates into Lumerical box coordinates.
    for i,disk in enumerate(individual):
        x,y,r = tuple(disk)
        
        #Make symmetric copies of disks centered around x=0 y=0
        xbox1 = (xMax - xMin)*x # + xMin
        ybox1 = (yMax - yMin)*(y-0.5) #+ yMin
        rbox1 = (rMax - rMin)*r + rMin 
        coords.append([xbox1,ybox1,rbox1])
        
        if x!=0:
            xbox2 = -1*xbox1 #Making a symmetric copy of each disk
            ybox2 = ybox1
            rbox2 = rbox1 
            coords.append([xbox2,ybox2,rbox2])
    coords.sort()
    coords = np.array(coords)

    # Load the Lumerical commands into a string variable for later output. 
    for i in range(len(coords)):
        disktemp.clearReplacement()
        x,y,r = tuple(coords[i])
        disktemp.addReplacement("set(\"x\"," + str(x) + "); ") 
        disktemp.addReplacement("set(\"y\"," + str(y) + "); ") 
        disktemp.addReplacement("set(\"radius\"," + str(r) + "); ") 
        diskStr = diskStr + disktemp.getstring()
        diskStr = diskStr + "\n"
        meshStr = ""
        
        if mesh_override == True and i < len(coords)-1 and x>=0: 

            meshtemp.clearReplacement()
            x1,y1,r1 = tuple(coords[i])
            x2,y2,r2 = tuple(coords[i+1])
            gap = (x2 -x1 -r1 -r2)
            if gap <= 50e-9:
                mesh_size = max(gap/5,10e-9)
                xc = x1 + r1 + gap/2.0 
                meshtemp.addReplacement('set("name","mesh'+str(i)+'");')
                meshtemp.addReplacement('set("x",'+str(xc)+');')
                meshtemp.addReplacement('set("x span",'+str(gap)+');')
                meshtemp.addReplacement('set("dx",'+str(mesh_size)+');')
                meshtemp.addReplacement('set("dy",'+str(mesh_size)+');')
                meshtemp.addReplacement('set("dz",'+str(mesh_size)+');')
                meshStr = meshStr+meshtemp.getstring()
                diskStr = diskStr + meshStr+"\n"
            
                
    # Write the new lumerical script
    boxLx = "sidex=" + str(2*((xMax-xMin)+1.0e-6)) + ";\n" 
    boxLy = "sidey=" + str(yMax-yMin) + ";\n"
    filetemp = InputTemplate()
    filetemp.readTemplate(parent_dir + "/temp50Disk.lsf")
    filetemp.addReplacement(boxLx) 
    filetemp.addReplacement(boxLy) 
    filetemp.addReplacement(diskStr)
    #filetemp.addReplacement('save("rad100gap50_'+str(individual.index)+'");')
    filetemp.generateNewFile(run_dir+"/OutGA.lsf")
    return None 
    
#----------------------------------------------------------
def write_dump_file(index,individual,l1,l2,filename="disk.dump"): 
    
    [xMin,yMin,rMin,prMin] = l1 
    [xMax,yMax,rMax,prMax] = l2  

    with open(filename,"a") as f: 
        
        header = ["ITEM: TIMESTEP"]
        #Header of dump file 
        header.append(str(index))
        header.append("ITEM: NUMBER OF ATOMS")
        header.append(str(len(individual)))
        header.append("ITEM: BOX BOUNDS pp pp pp")
        header.append(str(xMin) +" "+ str(xMax)) 
        header.append(str(yMin) +" "+ str(yMax)) 
        header.append("0.00 1.0e-6")
        header.append("ITEM: ATOMS id type x y radius")
        header_str = "\n".join(header) 
        f.write(header_str)
        f.write("\n") 

        #Body of dump file 
        body = []
        for i,disk in enumerate(individual): 
            
            x,y,r = tuple(disk)
            xbox = (xMax - xMin)*x + xMin 
            ybox = (yMax - yMin)*y + yMin 
            if r <0: 
                print "flag",r, individual.index
            rbox = (rMax - rMin)*r + rMin 

            line = []
            line.append(str(i+1)) 
            line.append("1") 
            line.append(str(xbox)) 
            line.append(str(ybox)) 
            line.append(str(rbox))

            body.append(" ".join(line)) 

        body_str = "\n".join(body)
        f.write(body_str)
        f.write("\n")

def run_lumerical(l1,l2,individual):
    
     
    parent_dir = os.getcwd() 
    redirect_file = "SIM_OUTPUT"
    #index = -1
    index =  individual.index
    #Create temp folder for lumerical run if not already present
    run_dir = parent_dir + "/Lumerical_runs"
    phase_dir = run_dir +"/phase"

    try:
            os.makedirs(run_dir)
    except OSError as exception:
            if exception.errno != errno.EEXIST:
                    raise
    
    try:
            os.makedirs(phase_dir)
    except OSError as exception:
            if exception.errno != errno.EEXIST:
                    raise
    os.system("cp " + parent_dir+ "/blank.fsp " + run_dir)
    make_lumerical_input(index, individual,parent_dir,run_dir,phase_dir,l1,l2)
    write_dump_file(index,individual,l1,l2)
    LS = Simulation()
    LS.setWorkDir(run_dir)
    command = ["time","-p","fdtd-solutions","blank.fsp", "-nw","-run","OutGA.lsf"]
    LS.setProgArg(command)
    LS.setWorkDir(run_dir)

    # Try till the file is made 
    # Lumerical license issue: Can only run one instance 
    # Needed when running parallel lumerical GA
    make_file=False
    while make_file==False:
        try:
            LS.runProgram(redirect_file = redirect_file) 
            make_file=True
        except: 
            print "Making File..Lumerical busy"
            continue
    LS.clearProgArg()
    #command = ["time","-p","mpirun","fdtd-engine-ompi-lcl","rad100gap50_"+str(individual.index)+".fsp"]
    command = ["time","-p","mpirun","fdtd-engine-ompi-lcl","rad100gap50.fsp"]
    LS.setProgArg(command)

    run_sim=False
    while run_sim==False:
        try:
            LS.runProgram(redirect_file = redirect_file) 
            run_sim=True
        except:
            print "waiting to run simulation..Lumerical busy"
            continue

    LS.clearProgArg() 
    #command = ["time","-p","fdtd-solutions","rad100gap50_"+str(individual.index)+".fsp", "-nw","-run","post.lsf"]
    #command = ["time","-p","fdtd-solutions","rad100gap50.fsp", "-nw","-run","post.lsf"]
    command = ["time","-p","fdtd-solutions","rad100gap50.fsp", "-nw","-run","get_focaline.lsf"]
    LS.setProgArg(command)
    
    make_output=False
    while make_output==False:
        try: 
            LS.runProgram(redirect_file = redirect_file) 
            make_output=True
        except: 
            print "making output..Lumerical busy"
            continue

    return run_dir, index
#==============================================================

if __name__ == "__main__":
    main()
