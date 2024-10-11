import sys,os
# sys.path.append('../')
import pymembrane as mb 
import numpy as np
import pandas as pd
from termcolor import *
import json


system = mb.System()
fac=[0.4,0.5]
system.read_mesh_from_files(files={'vertices':'vertices_T192.inp','faces':'faces_T192.inp'})

# vertices = system.getVertices()
# Nv=len(vertices)
# id0s=np.random.choice(range(Nv),round(fac[0] * Nv),replace=False)
# id1s=np.random.choice(list(set(range(Nv))-set(id0s)),round(fac[1] * Nv),replace=False)
# id2s=list(set(range(Nv))-set(id0s)-set(id1s))

# for vi in id0s: vertices[vi].type=0
# for vi in id1s: vertices[vi].type=1
# for vi in id2s: vertices[vi].type=2

# system.setVertices(vertices)
# print("id 0:",id0s)
# print("id 1:",id1s)
# print("id 2:",id2s)


# simulation = mb.MembraneFactory(trimesh)
# simulation.dump("vtk","mesh")
dump = system.dump() 
dump.vtk("initial mesh")


evolver = mb.Evolver(system)
# harmonic bond potential
evolver.add_force("Mesh>Harmonic", 
{"k":{"0": str(10.)}, 
"l0":{"0": str(1.0)}})

# Bending Potential
evolver.add_force("Mesh>BendingGK", 
{"H0":{"0": "0."}, 
 "kappaH":{"0": "4."},
 "kappaG":{"0": "-0.66667*4."}})


# limit edge potential
# evolver.add_force("Mesh>Limit", 
# {"lmin":{"0": "1.3"}, 
#  "lmax":{"0": "0.7"}})

print(evolver.get_force_info())
evolver.add_integrator("Mesh>Brownian>Positions", {"T":"0.0"})

vertices = system.getVertices()
faces = system.getFaces()
for fi in range(len(faces)):
	print("original metric:",faces[fi].get_metric([vertices[faces[fi].v1].r,vertices[faces[fi].v2].r,vertices[faces[fi].v3].r]))
	faces[fi].set_refmetric(1.2*np.array(faces[fi].get_metric([vertices[faces[fi].v1].r,vertices[faces[fi].v2].r,vertices[faces[fi].v3].r])))
	print("after:",faces[fi].get_refmetric())
system.setFaces(faces)

compute = system.compute_mesh()
energy = compute.compute_mesh_energy(evolver)


frame = 0
run_steps = 20000
dump_period = 1000

for frame in range(run_steps):
    evolver.evolveMD()
    # acc=simulation.simulate_monte_carlo(dump_period, temp)
    # print("accepted ratio:", acc[0][1] / (trimesh.getNumVertices + trimesh.getNumVertices * (1. - sum([x**2 for x in [fac[0],fac[1],1.-fac[0]-fac[1]]]))) / dump_period)
    if frame % dump_period ==0:
        print(evolver.get_force_info())
        print("mesh energy", sum(energy['vertices']))
        dump.vtk("growth_" + str(frame))
# simulation.dump('json','shell')