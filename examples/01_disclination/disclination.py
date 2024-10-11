#!/usr/bin/env python
# coding: utf-8


#import the code
import pymembrane as mb
import numpy as np


#create a system 
system = mb.System()
#read the mesh
N = 10 #pentagon size
vertex_file = 'vertices_N' + str(N) + '.inp'
face_file = 'faces_N' + str(N) + '.inp'
system.read_mesh_from_files(files={'vertices':vertex_file, 'faces':face_file})


#save the mesh to display
#create dumper
dump = system.dump() 
dump.vtk("initial mesh")


#add the evolver class where the potentials and integrators are added
evolver = mb.Evolver(system)


#add the potentials in this case we will add bending and stretching 
# stretching 
evolver.add_force("Mesh>Harmonic", {"k":{"0":str(50.0)}, 
                                    "l0":{"0":str(1.0)}})

# limit potential
evolver.add_force("Mesh>Limit", {"lmin":{"0":str(0.7)}, 
                                 "lmax":{"0":str(1.3)}})
# bending potential
kappa = 1.0
evolver.add_force("Mesh>Bending", {"kappa":{"1":str(kappa)}
                                       })

print(evolver.get_force_info())


#add the monte carlo integrator
# vertex move:
# first we need to know the edge length to move it appropiate:
compute = system.compute_mesh()
edge_lengths = compute.compute_edge_lengths()
avg_edge_length= np.mean(edge_lengths)
print("avg_edge_length = ", avg_edge_length)

evolver.add_integrator("Mesh>MonteCarlo>vertex>move", {"dr":"0.005"})


## Now we want to have 100 snapshots every 1000 steps so then
snapshots = 100
run_steps = 1000

## then we want to run the simulation for a temperature 1e-6
evolver.set_global_temperature(str(1e-4))

dump.vtk("pentagon_t0")
for snapshot in range(snapshots):
    evolver.evolveMC(MCS=run_steps)
    dump.vtk("pentagon_t" + str(snapshot*run_steps))

edge_lengths = compute.compute_edge_lengths()
avg_edge_length= np.mean(edge_lengths)
print("avg_edge_length = ", avg_edge_length)



