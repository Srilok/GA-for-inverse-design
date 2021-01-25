#!/home/share/cnm50256/bin/miniconda2/bin/python
"""
Srilok: GA based on DEAP to optimize focusing efficiency of flat lenses (1D probelm) 
Built on Henry's code for force field fitting 
"""
import sys
from datetime import datetime
import random, pickle, shutil
import numpy as np
import os 
from itertools import combinations
from deap import base,creator,tools

# ----- Default for optional setup -----
##
##
init_population_file = None 
expand_from_file = None
read_restart_file = None
write_restart_file = None 

# ----- User define ------

#from cal_obj import objective
from disk_operations import crossover, correct_if_overlap,invert
from disk_operations import check_overlap, guess_chain_disks 
from disk_operations import  remove_disk
#from disk_operations import add_disk
from disk_operations import add_disk_uniform as add_disk
from disk_operations import disk_mutate,constrain_disk
from disk_operations import sort_x 
#from disk_operations import replace_random  
from disk_operations import chain_replace_random as replace_random 
from sim_setup import  get_box_input, write_dump_file 
from read_dump import load_dump
from CalcObjective import set_target
from CalcObjective import analyzeFullSim as objective 

#-----------------------------------------------------------------------
# For debugging 
#def objective(individual,l1,l2,peaks):
#    
#    write_dump_file(individual.index,individual,l1,l2)
#    return np.random.random(),individual.index 
#    err = 0.0
#    for disk in individual: 
#        err += disk[0]*disk[2]

#    return err, individual.index
#-----------------------------------------------------------------------

index = 0 
parent_dir = os.getcwd()
BOX_settings = "in.box" 
max_generation = 5000
population_size = 80
number_of_ind_promoted =  24
probability_crossover = 0.3
probability_large_mutation = 0.4
probability_invert = 0.1
probability_small_mutation = 0.4
probability_add_disk = 0.0
probability_remove_disk = 0.0
probability_replace_random = 0.1
number_of_contestants = population_size - number_of_ind_promoted
#----------------------------------------------------------------------
# Check to see if the probabilities add upto 1 
p_tot = probability_large_mutation + probability_small_mutation \
                + probability_add_disk + probability_remove_disk \
                + probability_replace_random

#try: 
#    if not (0.99 < p_tot < 1.01):
#        raise
#except: 
#    print("Sum of the probabilities of operations does not add upto 1") 
#    print("Exiting..")
#    sys.exit()
#----------------------------------------------------------------------

init_population_file = 'restart.dump'
#init_population_file = None 
offspring_info_file = 'offspring.dat'
redundancy_monitor = 'redundancy.dat'
best_objective_file = 'best_objective.dat'
pool_evolution_file = 'pool.dat'
expand_from_file = 'foreigner.txt'
index_file = 'index_count'
#read_restart_file = 'checkpt.ga.restart'
write_restart_file = 'checkpt.ga'
write_restart_every = 1 # unit: generations
hof = tools.HallOfFame(12)


large_sigma_x = 200.0 /6000.0 
large_sigma_r = 20.0/75.0

small_sigma_x = 20.0/6000.0 
small_sigma_r = 10.0/75.0

large_sigma = [large_sigma_x,0.0,large_sigma_r] 
small_sigma = [small_sigma_x,0.0,small_sigma_r] 

# ------------------------------------------
def setup_ga():

    
    #Set the target phase profile and get the peaks 
    peaks = set_target()

    # -1 for minimizing, 1 for maximizing, multi values for multi-obj
    creator.create("optimizer", base.Fitness, weights=(-1,))
    creator.create("ind_type", list, index= int, history= "Initial",
                    trans = 0.0,fitness=creator.optimizer)
    
    #Read the box constrains from input  
    [nMin, nMax,pMin,pMax, xMin, xMax, yMin, yMax, rMin, rMax] = get_box_input(BOX_settings)
    global l1, l2
    l1 = [xMin,yMin,rMin,pMin] 
    l2 = [xMax,yMax,rMax,pMax]
    l12 = [ max_val - min_val for min_val,max_val in zip(l1,l2)]

    # -- define ga operations --
    wf.register("crossover", crossover,l1=l1,l2=l2)
    wf.register("large_mutate", disk_mutate, indpb=0.5,
                mu=0.0, sigma=large_sigma)
    
    wf.register("small_mutate",disk_mutate, indpb = 1.0, 
                mu = 0.0,sigma=small_sigma)

    wf.register("AddDisk", add_disk, nAdd =1, l1=l1,l2=l2)

    wf.register("RemoveDisk",remove_disk, nDel = 1)
    
    wf.register("invert",invert)

    def selUnique(individuals, k):
        return random.sample(individuals,k)
    
    #wf.register("select", tools.selTournament, tournsize=3)
    
    #wf.register("select", tools.selRoulette)
    
    def selRoulette_fitness(individuals, k, fit_attr="fitness"):
        obj_list =[]
        for ind in individuals:
                obj = ind.fitness.values[0]
                obj_list.append(obj)
        #obj_list.append(float(ind.fitness.values[0]) for ind in individuals)
        best_obj = min(obj_list) 
        worst_obj = max(obj_list)
        RelFitness =[]
        for ind in individuals: 
            RelFitness.append((ind.fitness.values[0] - worst_obj)/(best_obj - worst_obj))
            
        sum_fits = sum(relfit for relfit in RelFitness)
        chosen = []
        for i in xrange(k):
            u = random.random() * sum_fits
            sum_ = 0 
            for j,ind in enumerate(individuals):
                sum_ += RelFitness[j]
                if sum_ > u:
                    chosen.append(ind)
                    break
#        with open("chosen","a") as g:
#                for ind in chosen:
#                    g.write("%d"%ind.index) 
#                    g.write("\n") 
#
        return chosen 
   
    #wf.register("select", selRoulette_fitness)
    wf.register("select", selUnique)

    wf.register("evaluate", objective,l1=l1,l2=l2,peaks=peaks)

    wf.register("CorrectOverlap",correct_if_overlap,l1=l1,l2=l2,mode = "delete")

    wf.register("sort_best", tools.selBest)

    wf.register("random_sel",tools.selRandom)                               
    
    wf.register("limit_range",constrain_disk,
                    nMin=nMin, nMax=nMax,l1=l1,l2=l2)

    #wf.register("guess_random_disks", guess_random_disks, 
    #                nMin=nMin, nMax=nMax, l1=l1, l2=l2)
    
    wf.register("guess_random_disks", guess_chain_disks, 
                    l1=l1, l2=l2)
                
                
    wf.register("replace_by_random", replace_random,l1 =l1,l2 =l2)
    
    wf.register("sortX",sort_x)

    wf.register("individual", tools.initIterate, creator.ind_type,
                    wf.guess_random_disks)

    wf.register("pop_type", tools.initRepeat, list, wf.individual)

    #---------------Assign index------------------
    def assign_index(individual): 
        global index
        index += 1
        individual.index = index
    wf.register("assign_index",assign_index)


    # -- define procedure to replace an individual --
    def replace_individual(individual,replacement,note=None):
        ndisk = len(individual)
        # Make individual empty 
        for _ in range(ndisk): 
            individual.pop() 
        # Add replacement 
        for disk in replacement: 
            individual.append(disk)
        if not note is None: 
            individual.history = note 
        return individual
    wf.register("replace_individual",replace_individual)

    # -- define procedure to expand a population --
    def expand_population(population,new_size):
        new_population = wf.pop_type(n=new_size-len(population))
        return population + new_population
    wf.register("expand_population",expand_population)


    # -- define procedure to check redundancy --
    # Returns True if redundant
    def check_redundancy(ind1,ind2,tol): 
        n1 = len(ind1) 
        n2 = len(ind2)
        test1 = np.array(ind1) 
        test2 = np.array(ind2)
        tol_x = tol[0] 
        tol_r = tol[1]
        if n1 == n2: 
            if np.max(np.abs(test1[:,0]-test2[:,0])) > tol_x \
                    or np.max(np.abs(test1[:,2]-test2[:,2])) > tol_r: 
                return False 
            else: 
                return True
        else: 
            return False
    wf.register("is_same",check_redundancy)

# ------------------------------------------
def run_ga():

    global hof, read_restart_file, expand_from_file

    print('-'*60)
    #---------RESET RUN INFO------------------    
    #Reset objective info 
    with open("obj_split.dat","w") as f: 
        f.write("# Objective split up of individuals \n")

    
    #Reset index count 
    #if read_restart_file is None and init_population_file is None: 
    with open(index_file,"w") as f: 
        f.write("0")

    #Reset offspring file 
    with open(offspring_info_file,"w") as f: 
        f.write("#Offspring fitness and history")
        f.write("\n")

    #Reset objective file 
    with open(best_objective_file,"w") as f: 
        f.write("#Evolution of the best objective with generation")
        f.write("\n")
        f.write("#Generation \t Objective")
        f.write("\n")

    #Reset redundancy file  
    with open(redundancy_monitor,"w") as f: 
        f.write("#Redundancy Monitor")
        f.write("\n")

    #Reset Pool summary 
    with open(pool_evolution_file,"w") as f: 
        f.write("# Pool Summary ")
        f.write("\n")
    #------------------------------------------
    
    # restart an evolution
    if not read_restart_file is None:
        try:
            gen, population, hof =\
            restart(read_restart_file, expand_from_file)
            if not init_population_file is None:
                replace_from_file(population,init_population_file,
                                l1=l1,l2=l2,addmsg = ',start new')
            
        except:
            print("cannot restart from %s, start new"%
                  read_restart_file)
            read_restart_file = None

    # create a new evolution
    if read_restart_file is None:
        gen = 0
        population = wf.pop_type(n=population_size)
        print("create a new population of size %d"%population_size)
        replace_from_file(population,init_population_file,l1=l1,l2=l2,
                          addmsg=', start new')
        print("run %d generations"%max_generation)

    print('-'*60)
    sys.stdout.flush()

    if gen >= max_generation:
        sys.exit()

    wf.map(wf.assign_index,population)
    wf.map(wf.sortX,population)

    fitnesses = list(wf.map(wf.evaluate, population))

    for ind, obj in zip(population, fitnesses):
        ind.fitness.values = obj[0],
        ind.trans = obj[1]

    # Extract all the fitnesses
    objs = [ind.fitness.values[0] for ind in population]

    # Begin the evolution
    while gen < max_generation:
        
        population = wf.sort_best(population, population_size)
        write_pool(population,gen)
        print("-- Generation %i -- %s"%(gen,
            datetime.now().strftime("%d %B, %Y at %H:%M:%S.")))
        sys.stdout.flush()

        # Select the next generation individuals 
        offspring = wf.select(population, len(population)-number_of_ind_promoted)
        
        #Clone the selected individuals
        offspring = list(wf.map(wf.clone, offspring))
        
        # Apply crossover and mutation on the offspring
        for child1, child2 in zip(offspring[::2], offspring[1::2]):
            if random.random() < probability_crossover:
                # If the the parents are same, mutate one of them before cross-over
                if wf.is_same(child1,child2,[1e-3,1e-3]): 
                    if np.random.random() < 0.5:
                        wf.large_mutate(child1)
                        wf.limit_range(child1)
                    else: 
                        wf.small_mutate(child1) 
                        wf.limit_range(child1)
                rand = random.random()
                if rand <= 0.4: 
                    wf.crossover(child1, child2, mode ='single')
                    history = "1Pt crossover(%d & %d)"%(child1.index,child2.index)
                    child1.history = history 
                    child2.history = history
                elif  0.4 <= rand < 0.9: 
                    wf.crossover(child1, child2, mode ='double')
                    history = "2Pt crossover(%d & %d)"%(child1.index,child2.index)
                    child1.history = history 
                    child2.history = history
                else:
                    wf.crossover(child1, child2, mode ='triple')
                    history = "3Pt crossover(%d & %d)"%(child1.index,child2.index)
                    child1.history = history 
                    child2.history = history

                del child1.fitness.values
                del child2.fitness.values
                #write_dump_file(child1.index,child1,l1,l2,"crossed_over.dump")
                #write_dump_file(child2.index,child2,l1,l2,"crossed_over.dump")

        # Mutation, disk removal and addition
        for mutant in offspring:
            red_test = wf.clone(mutant)
            rand = random.random() 
            prob = probability_large_mutation 
            if mutant.fitness.valid:
                if rand < prob:
                    wf.large_mutate(mutant)
                    if not wf.is_same(mutant,red_test,[1e-3,1e-3]):
                        del mutant.fitness.values
                        history = "Large mutation( %d )"%mutant.index
                        mutant.history = history 
                        wf.sortX(mutant)

            prob += probability_small_mutation
            if mutant.fitness.valid: 
                if rand < prob:
                    wf.small_mutate(mutant)
                    if not wf.is_same(mutant,red_test,[1e-3,1e-3]):
                        del mutant.fitness.values
                        history = "Small mutation (%d)"%(mutant.index)
                        mutant.history = history
                        wf.sortX(mutant)

           # prob += probability_add_disk
           # if mutant.fitness.valid:
           #     if rand < prob: 
           #         print mutant.index,"Add Disk"
           #         wf.AddDisk(mutant)
           #         if np.max(red_test - np.array(mutant)) > 1e-3:
           #         #if red_test!= mutant:
           #             del mutant.fitness.values
           #             history = "Added random disk(%d)"%(mutant.index)
           #             mutant.history = history 
           #             wf.sortX(mutant)
           #             write_dump_file(mutant.index,ind,l1,l2,filename="after_operation.dump")

           # prob += probability_remove_disk
           # if mutant.fitness.valid:
           #     if rand < prob:
           #         print mutant.index,"Remove disk"
           #         wf.RemoveDisk(mutant)
           #         if np.max(red_test - np.array(mutant)) > 1e-3:
           #         #if red_test!= mutant:
           #             del mutant.fitness.values
           #             history = "Removed random disk(%d)"%(mutant.index)
           #             mutant.history = history 
           #             wf.sortX(mutant)
           #             write_dump_file(mutant.index,ind,l1,l2,filename="after_operation.dump")

            prob += probability_replace_random
            if mutant.fitness.valid:
                if rand < prob:
                    wf.replace_by_random(mutant)
                    del mutant.fitness.values
                    history = "Replaced by random individual(%d)"%(mutant.index)
                    mutant.history = history 
                    wf.sortX(mutant)
            
            prob += probability_invert
            if mutant.fitness.valid:
                if rand < prob:
                    wf.invert(mutant)
                    del mutant.fitness.values
                    history = "Inverted(%d)"%(mutant.index)
                    mutant.history = history 
                    wf.sortX(mutant)
        for ind in offspring: 
            if ind.fitness.valid: 
                wf.replace_by_random(ind)
                del ind.fitness.values
                history = "Replaced by random individual(%d) AFTER-check"%(ind.index)
                ind.history = history 

        # constrain disks within bounds
        for ind in offspring:
            if not ind.fitness.valid:
                constrain_met = wf.limit_range(ind)
                if not constrain_met:
                    history = "\t {Adjusted for constrains}"
                    ind.history += history
                    wf.sortX(mutant)
        
        #Redundancy gaurd 
        test_list = list(wf.map(wf.clone,population))
        for new_ind in offspring: 
            for test_ind in test_list: 
                    if wf.is_same(new_ind,test_ind,[1e-3, 1e-2]): 
                        wf.replace_by_random(new_ind)
                        del new_ind.fitness.values
                        history = "Replaced by random individual(%d) RED "%(new_ind.index)
                        new_ind.history += history 

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]

        #Reassign index
        wf.map(wf.assign_index,invalid_ind)
        fitnesses = list(wf.map(wf.evaluate, invalid_ind))
        
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit[0], 
            ind.trans = fit[1]
        ##########################################################
        # UPDATE POOL
        write_offspring(offspring,gen) 
        
        # Candidates for next generation
        candidates = list(wf.map(wf.clone,offspring + population[number_of_ind_promoted:]))
        # Sort candidates according to fitness
        candidates = wf.sort_best(candidates, len(candidates))
        
        #Update population
        population[number_of_ind_promoted:] = list(wf.map(wf.clone, candidates[:number_of_contestants]))
        #population[number_of_ind_promoted:] = list(wf.map(wf.clone, offspring))
        hof.update(population)
        print_gen(population,gen,hof)
        gen += 1

        #Pool dump files 
        for ind in population: 
            write_dump_file(ind.index,ind,l1,l2,filename="pool.dump")

        # write restart file
        if (not write_restart_file is None) and (gen % write_restart_every == 0):
            cp = dict(population=population, generation=gen, 
                      halloffame=hof, rndstate=random.getstate())
            with open(write_restart_file+'.tmp', "wb") as f:
                pickle.dump(cp,f)
            shutil.move(write_restart_file+'.tmp',write_restart_file)

# ------------------------------------------
def read_next_line(fileobj,comments='#'):

    # read next line in file, skip comments
    for line in fileobj:
        ix = line.find(comments)
        if ix >= 0:
            line = line[:ix]
        line = line.strip()
        if not line:
            continue
        return line.split()
    return None

# ------------------------------------------
def print_gen(population,gen,hof=None):

    # Gather all the fitnesses in one list and print the stats
    objs = [ind.fitness.values[0] for ind in population]
    pop_size = len(population)
    mean = sum(objs) / pop_size
    sum2 = sum(x*x for x in objs)
    std = abs(sum2 / pop_size - mean**2)**0.5
    with open(best_objective_file,"a") as f: 
        f.write("%d \t %.5f"%(gen,min(objs)))
        f.write("\n")
    print("  Min %s  Max %s  Avg %s  Std %s" % 
          (min(objs),max(objs),mean,std))
    if not hof is None:
        print("Hall Of Hame")
        print("============")
        for i,hero in enumerate(hof):
                print("Index:%d \t Fitness:%.2f \t Trans: %.2f \t History:%s"%(hero.index, 
                    hero.fitness.values[0],hero.trans,hero.history)) 
    return

# ------------------------------------------
def replace_from_file(population,load_file,l1,l2,i=0,addmsg=''):

    # replace individuals by individuals from a file
    if load_file is None:
        return
    try:
        start_pool,note = load_dump(load_file,l1,l2) 
        len_pop = len(population)
        if len(start_pool) ==1: 
            i =0
            for ind in population:
                wf.replace_individual(ind,start_pool[0],note=note[0])
                i += 1
                added = i
                print('added %d using %s'%(added,load_file))
        else:
            i = 0 
            for i,ind in enumerate(start_pool): 
                wf.replace_individual(population[len_pop-i-1],ind,note=note[i])
                i += 1
                added = i
                print('added %d using %s'%(added,load_file))
    except:
        print('cannot load %s'%load_file+addmsg)

# ------------------------------------------
def restart(restart_file,expand_from_file):

    # load restart file
    with open(restart_file, "r") as f:
        cp = pickle.load(f)
    gen = cp["generation"]
    population = cp["population"]
    random.setstate(cp["rndstate"])
    hof = cp["halloffame"]
    old_size = len(population)
    run_gen = max(max_generation-gen,0)
    print("restart from %s , gen %d (population size %d)"%
          (restart_file,gen,old_size))
    print('-'*60)
    print("-- Generation %i --"%gen)
    print_gen(population,gen,hof)
    print('-'*60)
    #if run_gen == 0:
    #    print("run %d new generations"%run_gen)
    #else:
    #    if old_size < population_size:
    #        print("expand population size %d to %d"%
    #              (old_size,population_size))
    #        population = wf.expand_population(population,population_size)
    #        replace_from_file(population,expand_from_file,old_size,
    #                          addmsg=',  random only')
    #        print("run %d new generations"%run_gen)
    #    elif old_size > population_size:
    #        print("shrink population size %d to %d "%
    #              (old_size,population_size)+
    #              ", run %d new generations"%run_gen)
    #        population = wf.select(population,population_size)
    #    else:
    #        print("kept population size %d , run %d new generations"%
    #              (population_size,run_gen))
    #wf.map(wf.limit_range,population)
    return (gen,population,hof)

#-------------------------------------------
def check_redundancy(offspring,gen): 
    with open(redundancy_monitor,"a") as f: 
        f.write("-"*60)
        f.write("Generation:%d"%gen)
        f.write("\n")
        for i in range(len(offspring)): 
            for j in range(i+1,len(offspring)): 
                if offspring[i] == offspring[j]:
                    f.write("Offspring %d and offspring %d are identical" \
                                    %(offspring[i].index, offspring[j].index))
                    f.write("\n")
    return 
# ------------------------------------------
def write_offspring(offspring,gen): 
    check_redundancy(offspring,gen)
    with open(offspring_info_file,"a") as f: 
        f.write("-"*60+ "\n")
        f.write("Generation: %d"%gen)
        f.write("\n")
        for ind in offspring:
            f.write("Index:%d \t Fitness:%.5f \t Trans:%.5f \t History:%s" \
                            %(ind.index, ind.fitness.values[0],ind.trans,ind.history))
            f.write("\n")
#-------------------------------------------
def write_pool(population,gen): 
    #check_redundancy(offspring,gen)
    with open(pool_evolution_file,"a") as f: 
        f.write("-"*60+ "\n")
        f.write("Generation: %d"%gen)
        f.write("\n")
        for ind in population:
            f.write("Index:%d  \t Fitness:%.5f \t Trans:%.5f \t History:%s" \
                            %(ind.index, ind.fitness.values[0],ind.trans, ind.history))
            f.write("\n")
#---------------------------------------------

wf = base.Toolbox()
setup_ga()

# ------------------------------------------
if __name__ == "__main__":
    try:
        backend = int(sys.argv[1])
        if not backend in (0,1,2,-1,-2):
            raise
    except:
        print("for serial,\n  python ga.py 0\n")
        print("for one node (use multiprocessing),\n"+
        "  python ga.py 1 Nworker\n")
        print("for many nodes (use scoop),\n"+
        "  python -m scoop --hostfile hosts -n Nworker ga.py 2\n")
        sys.exit()


    if backend == -1: 
        #Set the target phase profile and get the peaks 
        peaks = set_target()
        [nMin, nMax, xMin, xMax, yMin, yMax, rMin, rMax] = get_box_input(BOX_settings)
        global l1, l2 
        l1 = [xMin,yMin,rMin] 
        l2 = [xMax,yMax,rMax]
        start_pool = load_dump(init_population_file,l1=l1,l2=l2)
        print(wf.evaluate(start_pool[0]))
   
    if backend == -2:

        [nMin, nMax, xMin, xMax, yMin, yMax, rMin, rMax] = get_box_input(BOX_settings)
        l1 = [xMin,yMin,rMin] 
        l2 = [xMax,yMax,rMax]
        num = 20 
        population = wf.pop_type(num)
        print np.shape(np.array(population[0]))[-1]
        for i in range(num):  
            start_pool = load_dump(init_population_file,l1=l1,l2=l2)
            wf.replace_individual(population[i],start_pool[0])
            rand = random.random() 
            if rand < 0.4: 
                wf.small_mutate(population[i]) 
            if 0.4 < rand < 0.6: 
                wf.RemoveDisk(population[i])
            if 0.6 < rand < 0.8: 
                wf.AddDisk(population[i])
            if rand > 0.8: 
                rand_ind = wf.guess_random_disks()
                print rand_ind
                wf.crossover(population[i],rand_ind)

            wf.limit_range(population[i])
            write_dump_file(population[i].index,population[i],l1=l1,l2=l2)         
    
    if backend == 0: # serial
        wf.register("map",map)
        run_ga()

    if backend == 1: # multiprocessing
        import multiprocessing as mp
        Nproc = min(int(sys.argv[2]),population_size)
        pool = mp.Pool(processes=Nproc) #,initializer=pool_init,initargs=
        wf.register("map",pool.map)
        run_ga()

    elif backend == 2: # scoop
        # usage: python -m scoop --hostfile hosts -n Nproc ga.py
        from scoop import futures
        wf.register("map",futures.map)
        run_ga()

    print("-- Done -- %s"%(
        datetime.now().strftime("%d %B, %Y at %H:%M:%S.")))
    sys.stdout.flush()
