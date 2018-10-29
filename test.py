import ga 
from disk_operations import guess_chain_disks ,add_disk2
from disk_operations import disk_mutate
from sim_setup import write_dump_file

print ga.l1 
print ga.l2 

for i in range(1000): 
    individual = guess_chain_disks(ga.l1,ga.l2)
    individual = disk_mutate(individual,0.5,0,0.05)
    write_dump_file(i+1,individual,ga.l1,ga.l2,"test.dump")

print len(individual)
disk = individual[2]
print disk
individual.remove(disk)
print len(individual)

add_disk2(individual,1,ga.l1,ga.l2)
print len(individual) 

write_dump_file(1,individual,ga.l1,ga.l2,"add_test.dump")
