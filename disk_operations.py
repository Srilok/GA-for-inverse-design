import numpy as np 
import random 
from sim_setup import write_dump_file 

from itertools import combinations 
from deap import base,creator,tools
from itertools import repeat
from collections import Sequence 


#--------------------------------------------
# -- define individual and population --
def guess_random_disks(nMin,nMax,l1,l2):
    
    [xMin,yMin,rMin] = l1 
    [xMax,yMax,rMax] = l2  
    
    #print("Creating Initial random disks")
    individual = []
    ndisk = int(round((nMax - nMin)*np.random.random())) +nMin 
    for disk in range(ndisk): 
        r = np.random.random() 
        rbox = (rMax - rMin)*r + rMin

        #Buffer Values along the x and y dimensions to ensure disks
        #don't fall outside the box 
        ybuffer = rbox/(yMax- yMin)
        xbuffer = rbox/(xMax -xMin)
    
        #y = (1-2*ybuffer)*np.random.sample() +ybuffer
        y = 0.5
        x = (1-2*xbuffer)*np.random.sample() +xbuffer
    
        individual.append([x,y,r])

    correct_if_overlap(individual,l1,l2,mode ="delete")
    return individual
#-----------------------------------------------------
def guess_chain_disks(l1,l2):
    
    [xMin,yMin,rMin,pMin] = l1 
    [xMax,yMax,rMax,pMax] = l2  
    
    #print("Creating Initial random disks")
    individual = []
    xchain = 0
    cnt = 0 
    rchain = 0
    
    #Reduced values of pitch (normalized with box length along x)
    prMax = pMax/(xMax-xMin)
    prMin = pMin/(xMax-xMin)

    while xchain+rchain < 1:

        #Inserting the 1st disk
        if cnt ==0:
            r = np.random.random()
            rbox = (rMax - rMin)*r + rMin

            #Buffer Values along the x and y dimensions to ensure disks
            #don't fall outside the box 
            ybuffer = rbox/(yMax- yMin)
            xbuffer = rbox/(xMax -xMin)
                
            #y = (1-2*ybuffer)*np.random.sample() +ybuffer
            y = 0.5
            x = 0.0
            xchain = x
            rchain = xbuffer 

        else: 
            r = np.random.random()
            rbox = (rMax - rMin)*r + rMin

            #Buffer Values along the x and y dimensions to ensure disks
            #don't fall outside the box 
            ybuffer = rbox/(yMax- yMin)
            xbuffer = rbox/(xMax -xMin)
                
            
            y = 0.5
            gap = (prMax-prMin)*np.random.random() + prMin
            x = xchain + rchain + xbuffer + gap
                
            xchain = x  
            rchain = xbuffer 

        cnt += 1 
        if xchain + rchain < 1:
            individual.append([x,y,r])

    sort_x(individual)
    return individual


# -----------------------------------------
# Constrain the disks with the user specified bounds 

def constrain_disk(individual,nMin,nMax,l1,l2):

#    [xMin,yMin,rMin] = l1 
#    [xMax,yMax,rMax] = l2  
   
    [xMin,yMin,rMin,pMin] = l1 
    [xMax,yMax,rMax,pMax] = l2  

    ndisk = len(individual)
    overlap = False
#    constrain_met = None 
#    if ndisk < nMin:
#        nAdd = nMin - ndisk
#        individual = add_disk(individual,nAdd,l1,l2)
#        constrain_met = False 
#    if ndisk > nMax:
#        nDel = ndisk - nMax
#        individual = remove_disk(individual,nDel) 
#        constrain_met = False 
    
    for disk in individual:
     
        r = disk[2] 
        x = disk[0] 
        y = disk[1]
            
        if r <=0 or r>=1:
            r = np.random.random()
        rbox = (rMax - rMin)*r + rMin

        #Buffer Values along the x and y dimensions to ensure disks
        #don't fall outside the box 
        ybuffer = rbox/(yMax- yMin)
        xbuffer = rbox/(xMax -xMin)
     
        disk[2] = r
        
        if x < xbuffer or x < 0: 
            disk[0] = 0
        else:
            disk[0] = min(x,1-xbuffer) 
        ## 1D CASE
        if y !=  0.5: 
            disk[1] = 0.5
        
        
        ## 2D CASE
        #if y < ybuffer or y <0: 
        #    disk[1] = 0
        #else:
        #    disk[1] = min(y,1-ybuffer)
        
        if disk[0] < 0 or disk[1] < 0:
            print "WARNING: Negative co-ordiantes after constraining"
            print xbuffer 
            print ybuffer 

    overlap = correct_if_overlap(individual,l1,l2,mode = "delete")
    add_disk_uniform(individual,1,l1,l2)
    
    #

    if overlap: 
        constrain_met = False
    else: 
        constrain_met = True

    return constrain_met 
     
#-----------------------------------------
def disk_mutate(individual,indpb,mu,sigma):

    #print("Mutating")
    size = len(individual)
    num_coords = np.shape(np.array(individual))[-1]

    if not isinstance(mu, Sequence):
        mu_seq = [mu for _ in range(num_coords)]
        mu = mu_seq
    elif len(mu) < num_coords:
        raise IndexError("mu must be at least the size of number of coordinates: %d < %d" 
                            % (len(mu), num_coords))
    if not isinstance(sigma, Sequence):
        sigma_seq = [sigma for _ in range(num_coords)]
        sigma = sigma_seq
    elif len(sigma) < num_coords:
        raise IndexError("sigma must be at least the size of number of coordinates: %d < %d" 
                            % (len(sigma), num_coords))
    
    for i in range(size):
        rand = random.random() 
    #    print rand, indpb, rand < indpb 
        if rand < indpb:
            for j,m,s in zip(range(num_coords),mu,sigma):
                ran = random.gauss(m,s)
     #           print ran,s
                #individual[i][j] += random.gauss(m, s)
                individual[i][j] += ran
            individual[i][1] = 0.5
    sort_x(individual)
    return individual
#-------------------------------------------------
def disk_mutate_x(individual,indpb,mu,sigma):

    #print("Mutating")
    size = len(individual)
    num_coords = np.shape(np.array(individual))[-1]

    #if not isinstance(mu, Sequence):
    #    mu_seq = [mu for _ in range(num_coords)]
    #    mu = mu_seq
    #elif len(mu) < num_coords:
    #    raise IndexError("mu must be at least the size of number of coordinates: %d < %d" 
    #                        % (len(mu), num_coords))
    #if not isinstance(sigma, Sequence):
    #    sigma_seq = [sigma for _ in range(num_coords)]
    #    sigma = sigma_seq
    #elif len(sigma) < num_coords:
    #    raise IndexError("sigma must be at least the size of number of coordinates: %d < %d" 
    #                        % (len(sigma), num_coords))
    

    #xfix = 5.0e-6/12.0e-6 
    xfix = 0.0 
    for i in xrange(size):
        if individual[i][0] > xfix:
            if random.random() < indpb:
                individual[i][0] += random.gauss(mu, sigma)
                #Check to see if the mutation takes the disk into the first half 
                while individual[i][0] < xfix: 
                        individual[i][0] += random.gauss(mu,sigma)
    sort_x(individual)
    return individual
#-------------------------------------------------------------
def disk_mutate_r(individual,indpb,mu,sigma):

    #print("Mutating")
    size = len(individual)
    num_coords = np.shape(np.array(individual))[-1]

    #if not isinstance(mu, Sequence):
    #    mu_seq = [mu for _ in range(num_coords)]
    #    mu = mu_seq
    #elif len(mu) < num_coords:
    #    raise IndexError("mu must be at least the size of number of coordinates: %d < %d" 
    #                        % (len(mu), num_coords))
    #if not isinstance(sigma, Sequence):
    #    sigma_seq = [sigma for _ in range(num_coords)]
    #    sigma = sigma_seq
    #elif len(sigma) < num_coords:
    #    raise IndexError("sigma must be at least the size of number of coordinates: %d < %d" 
    #                        % (len(sigma), num_coords))
    
    #xfix = 5.0e-6/12.0e-6 
    xfix = 0.0  
    for i in xrange(size):
        if individual[i][0] > xfix:
            if random.random() < indpb:
                individual[i][2] += random.gauss(mu, sigma)
                #Check to see if the mutation takes the disk into the first half 
    sort_x(individual)
    return individual


#-------------------------------------------
def tail_mutate(individual,indpb,mu,sigma):

    #print("Mutating")
    size = len(individual)
    num_coords = np.shape(np.array(individual))[-1]
    
    if not isinstance(mu, Sequence):
        mu_seq = [mu for _ in range(num_coords)]
        mu = mu_seq
    elif len(mu) < num_coords:
        raise IndexError("mu must be at least the size of number of coordinates: %d < %d" 
                            % (len(mu), num_coords))
    if not isinstance(sigma, Sequence):
        sigma_seq = [sigma for _ in range(num_coords)]
        sigma = sigma_seq
    elif len(sigma) < num_coords:
        raise IndexError("sigma must be at least the size of number of coordinates: %d < %d" 
                            % (len(sigma), num_coords))
 
    #xfix = 5.0e-6/12.0e-6 
    xfix = 0.0  
    for i in xrange(size):
        if individual[i][0] > xfix: 
            if random.random() < indpb:
                for j,m,s in zip(range(num_coords),[0.0,0.0,0.0],sigma):
                    individual[i][j] += random.gauss(m,s)

                    #Check to see if the mutation takes the disk into the first half 
                    while individual[i][0] < xfix: 
                        individual[i][0] += random.gauss(m,s)
        individual[i][1] = 0.5
    sort_x(individual)
    return individual


#------------------------------------------
def remove_disk(individual,nDel):

    #print("Removing %d disks"%nDel)
    ndisk = len(individual)
    for _ in range(nDel):
        individual.pop()
    return individual
#------------------------------------------
# Add a disk without causing disk overlap
def add_disk(individual,nAdd,l1,l2):

    #print("Adding %d disks"%nAdd)

    [xMin,yMin,rMin] = l1 
    [xMax,yMax,rMax] = l2  
    cnt = 0
    for i in range(nAdd): 
        while cnt < 100:
            cnt += 1
            r = np.random.random() 
            rbox = (rMax -rMin)*r + rMin 

            #Buffer Values along the x and y dimensions to ensure disks
            #don't fall outside the box 
            ybuffer = rbox/(yMax- yMin)
            xbuffer = rbox/(xMax -xMin)

            #y = (1-2*ybuffer)*np.random.random() +ybuffer
            y = 0.5
            x = (1-2*xbuffer)*np.random.random() +xbuffer
            
            check_species = individual[:] 
            check_species.append([x,y,r])
            
            # Check if the added disk overlaps 
            # If overlap is true then continue 
            # If overlap is false break and return individual
            if not check_overlap(check_species,l1,l2): 
                individual.append([x,y,r])
                #print "Succesfully added disks"
                break
    
    return individual

#------------------------------------------
# Add a disk without causing disk overlap
def add_disk2(individual,nAdd,l1,l2):

    #print("Adding %d disks"%nAdd)

    [xMin,yMin,rMin,pMin] = l1 
    [xMax,yMax,rMax,pMax] = l2

    prMin = pMin/(xMax-xMin) 
    prMax = pMax/(xMax-xMin) 

    rmin_red = rMin/(xMax-xMin)
    rmax_red = rMax/(xMax-xMin)
    
    cnt = 0

    while True: 

        sort_x(individual)
        flag = False 
        
        ndisk = len(individual)
        for i in range(ndisk-1): 
            x1 = individual[i][0] 
            y1 = individual[i][1] 
            r1 = individual[i][2]
        
            x2 = individual[i+1][0] 
            y2 = individual[i+1][1] 
            r2 = individual[i+1][2]
        
            rbox1 = (rMax-rMin)*r1 +rMin
            rbox2 = (rMax-rMin)*r2 +rMin
        
            xbuffer1 = rbox1/(xMax -xMin)
            xbuffer2 = rbox2/(xMax -xMin)
        
            #Check to see if a disk can be added at the begining 
            if i ==0: 
                if x1 > rmin_red + prMin + xbuffer1: 
                    flag = True
                    gap_avail = x1 - xbuffer1 - prMin 
        
                    r_guess = np.random.random() 
                    xbuffer = (min(rmax_red,gap_avail)-rmin_red)*r_guess + rmin_red 
                    rbox = xbuffer*(xMax-xMin)
                    if rbox < rMin: 
                            print "flag",gap_avail,x1,xbuffer1,rmin_red, prMin,xbuffer  
                    r = (rbox-rMin)/(rMax-rMin)
                    x = 0
                    y = 0.5
                    individual.append([x,y,r])
                    break
                    #print "Adding disk at the begining of %d"%individual.index
            
            #Check to see if a disk can be added at the end 
            if i+1 ==ndisk-1: 
                if x2 < 1- (2*rmin_red + prMin + xbuffer2): 
                    flag = True 
                    gap_avail = 1 - (x2 + xbuffer2 + prMin)
        
                    r_guess = np.random.random() 
                    xbuffer = (min(rmax_red,gap_avail/2)-rmin_red)*r_guess + rmin_red 
                    rbox = xbuffer*(xMax-xMin)

                    rbuffer = min(gap_avail/2-xbuffer,prMax)
                
                    r = (rbox-rMin)/(rMax-rMin)
                    x = x2 + xbuffer2 + xbuffer +rbuffer*np.random.random() + prMin 
                    y = 0.5
                    individual.append([x,y,r])
                    break
                    #print "Adding disk at the end of %d"%individual.index

            #Check to see if there is enough spacing to insert a disk 
            if (x2-x1) > xbuffer1 + xbuffer2 + 2*prMin + 2*rmin_red:  
                flag = True 
                gap_avail = x2 - x1 - xbuffer1 - xbuffer2 -2*prMin
        
                #random number used to assign radius 
                r_guess = np.random.random()
                xbuffer = (min(rmax_red,gap_avail/2)-rmin_red)*r_guess + rmin_red
                rbox = xbuffer*(xMax-xMin)

                rbuffer = min(gap_avail/2-xbuffer,prMax)

                r = (rbox -rMin)/(rMax-rMin)
                y = 0.5
                x = x1 + xbuffer1+ xbuffer + rbuffer*np.random.random() +prMin

                #Inserting disk.
                individual.append([x,y,r])   
                #print "Adding disk at the middle of %d"%individual.index
                break

            # If the distance between disks are more than the maximum allowed
            # shift the right disk to the left. 
            elif (x2-x1) > xbuffer1 + xbuffer2 + prMax*1.05:
                    flag = True 
                    individual[i+1][0] = x1 + xbuffer1 + prMax +xbuffer2
                    break 

        if flag == False: 
            #print "Could not add any disks, disks already close packed!"
            break

    sort_x(individual)
    return individual

#-------------------------------------------
# Add a disk without causing disk overlap
def add_disk_uniform(individual,nAdd,l1,l2):

    #print("Adding %d disks"%nAdd)

    [xMin,yMin,rMin,pMin] = l1 
    [xMax,yMax,rMax,pMax] = l2

    prMin = pMin/(xMax-xMin) 
    prMax = pMax/(xMax-xMin) 

    rmin_red = rMin/(xMax-xMin)
    rmax_red = rMax/(xMax-xMin)
    
    cnt = 0

    while True: 

        sort_x(individual)
        flag = False 
        
        ndisk = len(individual)
        for i in range(ndisk-1): 
            x1 = individual[i][0] 
            y1 = individual[i][1] 
            r1 = individual[i][2]
        
            x2 = individual[i+1][0] 
            y2 = individual[i+1][1] 
            r2 = individual[i+1][2]
        
            rbox1 = (rMax-rMin)*r1 +rMin
            rbox2 = (rMax-rMin)*r2 +rMin
        
            xbuffer1 = rbox1/(xMax -xMin)
            xbuffer2 = rbox2/(xMax -xMin)
        
            # Check to see if 1st disk satisfies the gap constrain 
            # If not push it to within the gap range 
            if i ==0: 
                if (x1-xbuffer1) > prMax/2.0: 
                    flag = True
        
                    x = 0
                    r = np.random.random() 
                    y = 0.5
                    individual[i][0] =x 
                    individual[i][1] =y 
                    individual[i][2] =r
                    break
                    #print "Adding disk at the begining of %d"%individual.index
            
            #Check to see if a disk can be added at the end 
            if i+1 ==ndisk-1: 
                if x2 < 1- (2*rmax_red + prMax + xbuffer2): 
                    flag = True 

                    r = np.random.random() 
                    rbox = (rMax-rMin)*r + rMin 
                    xbuffer = rbox/(xMax-xMin)
                    
                    x = x2 + xbuffer2 + xbuffer + np.random.random()*(prMax-prMin)+prMin
                    y = 0.5
                    individual.append([x,y,r])
                    break
                    #print "Adding disk at the end of %d"%individual.index

            #Check to see if there is enough spacing to insert a disk 
            if (x2-x1) > xbuffer1 + xbuffer2 + 2*prMax + 2*rmax_red:  
                flag = True 
        
                #random number used to assign radius 
                r = np.random.random() 
                rbox = (rMax-rMin)*r + rMin 
                xbuffer = rbox/(xMax-xMin)
                    
                y = 0.5
                
                #Choose reference disk to guess to calculate the gap
                choice = random.choice([1,2])
                if choice == 1: 
                    x = x1 + (xbuffer1 + xbuffer + np.random.random()*(prMax-prMin) + prMin)
                if choice == 2:
                    x = x2 - (xbuffer2 + xbuffer + np.random.random()*(prMax-prMin) + prMin)

                #Inserting disk.
                individual.append([x,y,r])   
                #print "Adding disk at the middle of %d"%individual.index
                break

            # If the distance between disks are more than the maximum allowed
            # shift the right disk to the left. 
            if (x2-x1) > xbuffer1 + xbuffer2 + prMax:
                flag = True 
                individual[i+1][0] = x1 + xbuffer1 + xbuffer2 + np.random.random()*(prMax - prMin) + prMin
                break 

        if flag == False: 
            #print "Could not add any disks, disks already close packed!"
            break

    sort_x(individual)
    return individual


#-------------------------------------------
# To check if disk overlaps 
# Return ture if overlaps 
# Return False if no overlap

def check_overlap(individual,l1,l2):

#    [xMin,yMin,rMin] = l1 
#    [xMax,yMax,rMax] = l2  

    [xMin,yMin,rMin,pMin] = l1 
    [xMax,yMax,rMax,pMax] = l2  

    overlap = False
    ndisk = len(individual)
    
    for i in range(ndisk):

        x1,y1,r1 = individual[i] 

        xbox1 = (xMax - xMin)*x1 + xMin 
        ybox1 = (yMax - yMin)*y1 + yMin 
        rbox1 = (rMax - rMin)*r1 + rMin 
        
        for j in range(i+1,ndisk):

            x2,y2,r2 = individual[j] 
        
            xbox2 = (xMax - xMin)*x2 + xMin 
            ybox2 = (yMax - yMin)*y2 + yMin 
            rbox2 = (rMax - rMin)*r2 + rMin

            dx = xbox1 - xbox2
            dy = ybox1 - ybox2
            d = (dx*dx + dy*dy)**0.5
            
            R = rbox1 + rbox2
            if d < R + pMin:
                overlap = True
                return overlap 

    return overlap
#------------------------------------------------------
# Corrects the radius of the disk if there is a overlap
def correct_if_overlap(individual,l1,l2, mode ="delete"): 

    
#    [xMin,yMin,rMin] = l1 
#    [xMax,yMax,rMax] = l2  

    [xMin,yMin,rMin,pMin] = l1 
    [xMax,yMax,rMax,pMax] = l2  

    constrain_conflict = False
    overlap = False
    ndisk = len(individual)
    remove_disk_index = []
#    ########################################################
    if mode == "adjust":
        print("Mode adjust")
        for i in range(ndisk):

            x1,y1,r1 = individual[i] 

            xbox1 = (xMax - xMin)*x1 + xMin 
            ybox1 = (yMax - yMin)*y1 + yMin 
            rbox1 = (rMax - rMin)*r1 + rMin 
        
            for j in range(i+1,ndisk):

                x2,y2,r2 = individual[j] 
        
                xbox2 = (xMax - xMin)*x2 + xMin 
                ybox2 = (yMax - yMin)*y2 + yMin 
                rbox2 = (rMax - rMin)*r2 + rMin

                dx = xbox1 - xbox2
                dy = ybox1 - ybox2
                d = (dx*dx + dy*dy)**0.5

                if d <= (rbox1 + rbox2+pMin):
                    overlap = True
                    print "Correcting overlap..."

                    # Checks for the special case where disks are too close
                    # Radius adjustment is not possible without violating constrains
                    if d < 2*rMin:
                        constrain_conflict = True
                        choice = random.choice([i,j])
                        if choice not in remove_disk_index:  
                            remove_disk_index.append(choice)
                    else: 
                        if rbox1 > rbox2: 
                            if rMin > rbox1 - (rbox1 + rbox2 - d):
                                new_rbox1 = rMin
                                new_rbox2 = d - rMin 
                            else: 
                                new_rbox1 = rbox1 - (rbox1+ rbox2 -d) 
                                new_rbox2 = rbox2 
                        else:
                            if rMin > rbox2 - (rbox1 + rbox2 -d):
                                new_rbox2 = rMin 
                                new_rbox1 = d - rMin 
                            else:
                                new_rbox2 = rbox2 - (rbox1+rbox2 -d) 
                                new_rbox1 = rbox1

                        new_r1 = (new_rbox1 - rMin)/(rMax - rMin) - 0.01 
                        new_r2 = (new_rbox2 - rMin)/(rMax - rMin) - 0.01
                        
                        individual[i][2] = new_r1
                        individual[j][2] = new_r2
        
        # If constrain conflict is true remove the disks causing the conflict and add 
        # new disks without overlap
        if constrain_conflict: 
            rmdisk =[individual[i] for i in remove_disk_index]
            [individual.remove(disk) for disk in rmdisk]
            
            if not check_overlap(individual,l1,l2):
                #print "Disk Removed"
                #for k in range(len(remove_disk_index)):
                #add_disk(individual,len(remove_disk_index),l1,l2)
                add_disk_uniform(individual,len(remove_disk_index),l1,l2)
            else: 
                raise StandardError 
            
        sort_x(individual)
#        ####################################################

    if mode == "delete": 
        for i in range(ndisk):

            x1,y1,r1 = tuple(individual[i]) 

            xbox1 = (xMax - xMin)*x1 + xMin 
            ybox1 = (yMax - yMin)*y1 + yMin 
            rbox1 = (rMax - rMin)*r1 + rMin 
        
            for j in range(i+1,ndisk):

                x2,y2,r2 = tuple(individual[j]) 
        
                xbox2 = (xMax - xMin)*x2 + xMin 
                ybox2 = (yMax - yMin)*y2 + yMin 
                rbox2 = (rMax - rMin)*r2 + rMin

                dx = xbox1 - xbox2
                dy = ybox1 - ybox2
                d = (dx*dx + dy*dy)**0.5
                

                if d < (rbox1 + rbox2+pMin):
                    overlap = True
                    #print "Correcting overlap...delete"
                    choice = random.choice([i,j])
                    if choice not in remove_disk_index:  
                        remove_disk_index.append(choice)
       
        # If constrain conflict is true remove the disks causing the conflict and add 
        # new disks without overlap
        if overlap:
            rmdisk =[individual[i] for i in remove_disk_index]
            [individual.remove(disk) for disk in rmdisk]

            if not check_overlap(individual,l1,l2):
                #print "Disk Removed"
                #for k in range(len(remove_disk_index)):
                #add_disk(individual,len(remove_disk_index),l1,l2)
                add_disk_uniform(individual,len(remove_disk_index),l1,l2)
            else: 
                raise StandardError 
            
            #individual = disk_removed_ind[:]
        sort_x(individual)

    return overlap
#--------------------------------------------
          
def crossover(individual1,individual2,l1,l2,mode='single'): 
    
    if individual1 == individual2: 
        print "WARNING: Crossing over identical individuals...."            
        return
    
    if mode == 'single': 
            individual1, individual2 = crossover_1pt(individual1,individual2,l1,l2)
            return individual1, individual2
    if mode == 'double': 
            individual1, individual2 = crossover_2pt(individual1,individual2,l1,l2)
            return individual1, individual2
    if mode == 'triple': 
            individual1, individual2 = crossover_3pt(individual1,individual2,l1,l2)
            return individual1, individual2


#    write_dump_file(individual1.index,individual1,l1,l2,"crossed_over.dump")
#    write_dump_file(individual2.index,individual2,l1,l2,"crossed_over.dump")

#    
#
##    with open("crossover.dat","a") as f:
##        f.write("before crossover \n")
##        f.write("child1 \n")
##        for disk in individual1: 
##            f.write(str(disk))
##            f.write("\n") 
##        f.write("child2 \n")
##        for disk in individual2:
##            f.write(str(disk)) 
##            f.write("\n")
##        f.write("Equation of line:: x1:" + str(x1) + 
##                             "y1:" +str(y1) + "x2:" + str(x2) + "y2:" + str(y2) +"\n")
##        f.write("After crossover \n")
##        f.write("child1 \n")
##        for disk in individual1: 
##            f.write(str(disk))
##            f.write("\n") 
##        f.write("child2 \n")
##        for disk in individual2:
##            f.write(str(disk)) 
##            f.write("\n")

    return individual1, individual2
#--------------------------------------------------------
def crossover_1pt(individual1,individual2,l1,l2): 
    
    ind1 = np.array(individual1)
    xmin1, xmax1 = np.min(ind1[:,0]), np.max(ind1[:,0])
    ind2 = np.array(individual2)
    xmin2, xmax2 = np.min(ind2[:,0]), np.max(ind2[:,0])
    
    xmin = min(xmin1,xmin2)
    xmin = 5.0e-6/12.e-6
    xmax = max(xmax1,xmax2)

    
    #Generate points for random line    
    x1 = np.random.random()*(xmax-xmin) + xmin

    #with open("cross_over_debug.dat","a") as f: 
    #        f.write("x1: %.5f,  ind1:%d, ind2:%d"%(x1,
    #                individual1.index,individual2.index))
    #        f.write("\n")

    
    N_ind1 = len(individual1)
    N_ind2 = len(individual2)

    #Make a temp copy of the ind disk list
    temp_copy1 = individual1[:]
    temp_copy2 = individual2[:]

    #Empty individuals
    for _ in range(N_ind1):
        individual1.pop() 
    for _ in range(N_ind2):
        individual2.pop() 


    if individual1 !=[] or individual2!=[]: 
            raise 

    for disk in temp_copy1: 
        if disk[0] <= x1:
            individual1.append(disk) 
        else: 
            individual2.append(disk)

    for disk in temp_copy2: 
        if disk[0] <=x1: 
            individual2.append(disk) 
        else: 
            individual1.append(disk)

    
    sort_x(individual1) 
    sort_x(individual2)

    return individual1, individual2

def crossover_2pt(individual1,individual2,l1,l2): 
    
    ind1_seg1,ind1_seg2,ind1_seg3 = 0,0,0
    ind2_seg1,ind2_seg2,ind2_seg3 = 0,0,0

    ind1 = np.array(individual1)
    xmin1, xmax1 = np.min(ind1[:,0]), np.max(ind1[:,0])
    ind2 = np.array(individual2)
    xmin2, xmax2 = np.min(ind2[:,0]), np.max(ind2[:,0])
    
    xmin = min(xmin1,xmin2)
    xmin = 5.0e-6/12.e-6
    xmax = max(xmax1,xmax2)
    
    # While loop ensures atleast one disk between random cuts
    while True:
        #Generate points for random line    
        xr1 = np.random.random()*(xmax-xmin) + xmin
        xr2 = np.random.random()*(xmax-xmin) + xmin

        x1 = min(xr1,xr2) 
        x2 = max(xr1,xr2)
        
        for disk in individual1: 
            if disk[0] <=x1: 
                ind1_seg1+=1
            elif x1 < disk[0] <=x2: 
                ind1_seg2+=1 
            elif x2 < disk[0]: 
                ind1_seg3+=1 

        for disk in individual2: 
            if disk[0] <=x1: 
                ind2_seg1+=1
            elif x1 < disk[0] <=x2: 
                ind2_seg2+=1 
            elif x2 < disk[0]: 
                ind2_seg3+=1 

        seg_count = np.array([ind1_seg1,ind1_seg2,ind1_seg3,
                            ind2_seg1,ind2_seg2,ind2_seg3])


        if all(seg_count!= 0): 
                break

    #Make a temp copy of the ind disk list
    temp_copy1 = individual1[:]
    temp_copy2 = individual2[:]

    N_ind1 = len(individual1)
    N_ind2 = len(individual2)

    #Empty individuals
    for _ in range(N_ind1):
        individual1.pop() 
    for _ in range(N_ind2):
        individual2.pop() 

    if individual1 !=[] or individual2!=[]: 
        raise 

    for disk in temp_copy1: 
        if disk[0] <= x1 or disk[0] >=x2:
            individual1.append(disk) 
        else:
            individual2.append(disk)

    for disk in temp_copy2: 
        if disk[0] <=x1 or disk[0] >= x2: 
            individual2.append(disk) 
        else: 
            individual1.append(disk)

    sort_x(individual1) 
    sort_x(individual2)
    #print individual1,individual2
    return individual1, individual2

def crossover_3pt(individual1,individual2,l1,l2): 

    
    ind1_seg1,ind1_seg2,ind1_seg3,ind1_seg4 = 0,0,0,0
    ind2_seg1,ind2_seg2,ind2_seg3,ind2_seg4 = 0,0,0,0


    ind1 = np.array(individual1)
    xmin1, xmax1 = np.min(ind1[:,0]), np.max(ind1[:,0])
    ind2 = np.array(individual2)
    xmin2, xmax2 = np.min(ind2[:,0]), np.max(ind2[:,0])
    
    xmin = min(xmin1,xmin2)
    xmin = 5.0e-6/12.e-6
    xmax = max(xmax1,xmax2)

    # While loop ensures atleast one disk between random cuts
    while True:
        #Generate points for random line    
        xr1 = np.random.random()*(xmax-xmin) + xmin
        xr2 = np.random.random()*(xmax-xmin) + xmin
        xr3 = np.random.random()*(xmax-xmin) + xmin
        
        guess = [xr1,xr2,xr3]
        guess.sort()
        [x1,x2,x3] = guess

        for disk in individual1: 
                if disk[0] <=x1: 
                        ind1_seg1+=1
                elif x1< disk[0] <=x2: 
                        ind1_seg2+=1
                elif x2 < disk[0] <=x3: 
                        ind1_seg3+=1
                else: 
                        ind1_seg4+=1

        for disk in individual2:
                if disk[0] <=x1: 
                        ind2_seg1+=1
                elif x1< disk[0] <=x2: 
                        ind2_seg2+=1
                elif x2 < disk[0] <=x3: 
                        ind2_seg3+=1
                else: 
                        ind2_seg4+=1

        seg_count = np.array([ind1_seg1,ind1_seg2,ind1_seg3,ind1_seg4,
                            ind2_seg1,ind2_seg2,ind2_seg3,ind2_seg4])

        if all(seg_count!= 0): 
                break

    #Make a temp copy of the ind disk list
    temp_copy1 = individual1[:]
    temp_copy2 = individual2[:]

    N_ind1 = len(individual1)
    N_ind2 = len(individual2)

    #Empty individuals
    for _ in range(N_ind1):
        individual1.pop() 
    for _ in range(N_ind2):
        individual2.pop() 


    if individual1 !=[] or individual2!=[]: 
        raise 

    for disk in temp_copy1: 
        if disk[0] <= x1 or x2 <= disk[0] <=x3:
            individual1.append(disk) 
        else:
            individual2.append(disk)

    for disk in temp_copy2: 
        if disk[0] <=x1 or x2<= disk[0] <= x3: 
            individual2.append(disk) 
        else: 
            individual1.append(disk)

    
    #print "%d %d: %.5f %5f %.5f"%(individual1.index, individual2.index,x1,x2,x3) 
    sort_x(individual1) 
    sort_x(individual2)

    return individual1, individual2
#----------------------------------------------------------------------------------------
def invert(individual): 
    ndisk = len(individual)

    for i in range(ndisk): 
            individual[i][0] *=-1 
            individual[i][0] += 1

    return individual

def replace_random(individual,l1,l2): 
        
    [xMin,yMin,rMin] = l1 
    [xMax,yMax,rMax] = l2  

    ndisk = len(individual)

    for _ in range(ndisk): 
        individual.pop() 
        #print "Individual emptied" 
        #print individual

    for _ in range(ndisk): 
        cnt = 0 
        while cnt < 1000: 
            cnt +=1 
            r = np.random.random() 
            rbox = (rMax - rMin)*r + rMin 

            #Buffer Values along the x and y dimensions to ensure disks
            #don't fall outside the box 
            ybuffer = rbox/(yMax- yMin)
            xbuffer = rbox/(xMax -xMin)

            #y = (1-2*ybuffer)*np.random.random() +ybuffer
            y = 0.5
            x = (1-2*xbuffer)*np.random.random() +xbuffer

            individual.append([x,y,r])
    
            # Check if the added disk overlaps 
            # If overlap is true then continue 
            # If overlap is false break and return individual
            if not check_overlap(individual,l1,l2): 
                break
            else: 
                individual.remove([x,y,r])

    return individual 
#--------------------------------------------------------
def chain_replace_random(individual,l1,l2): 
        
    [xMin,yMin,rMin,pMin] = l1 
    [xMax,yMax,rMax,pMax] = l2  

    ndisk = len(individual)

    for _ in range(ndisk): 
        individual.pop() 
        #print "Individual emptied" 
        #print individual

    guessed_disks = guess_chain_disks(l1,l2)

    for disk in guessed_disks: 
        individual.append(disk)

    sort_x(individual)
    return individual 
#---------------------------------------------------------------
def tail_replace_random(individual,l1,l2): 
        
    sort_x(individual)
    [xMin,yMin,rMin,pMin] = l1 
    [xMax,yMax,rMax,pMax] = l2  

    ndel = 0 

    #xfix = 5.0e-6/12.0e-6 
    xfix = 0.0 
    for i,disk in enumerate(individual): 
        if disk[0] > xfix:  
           ndel+=1

    for _ in range(ndel): 
        individual.pop() 

    guessed_disks = guess_chain_disks(l1,l2)

    for disk in guessed_disks:
        if disk[0] > xfix: 
            individual.append(disk)

    sort_x(individual)
    return individual 
#---------------------------------------------------------------
#Sort disk according to the x co-ordinate
def sort_x(individual): 
    individual.sort() 
    return individual 


if __name__ == "__main__": 
        crossover() 
