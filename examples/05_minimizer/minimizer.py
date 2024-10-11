#First we need to import the module for the simulations
import pymembrane as mb
#numpy
import numpy as np
import sys

#create a system 
system = mb.System()
#read the mesh
vertex_file = 'vertices_T2500.inp'
face_file = 'faces_T2500.inp'
system.read_mesh_from_files(files = {'vertices': vertex_file, 'faces': face_file})

dump = system.dump() 
dump.vtk("initial mesh")

kb = 0.01
print("kb:", kb)
#add the evolver class where the potentials and integrators are added
evolver = mb.Evolver(system)

evolver.add_force("Mesh>Harmonic", {"k": {"0": str(1.)},
                                    "l0": {"0": str(1.0)}})

evolver.add_force("Mesh>Bending", {"kappa": {"0": str(kb)}})


evolver.add_minimizer("Mesh>Fire", {"dt": "0.1", "max_iter": "100000", "ftol": "0.0001", "etol": "0.0000001"})

print(evolver.get_force_info())
print(evolver.get_minimizer_info())

computer = system.compute_mesh()
energies = computer.compute_mesh_energy(evolver)
print("initial total energy:", sum(energies["vertices"]) + sum(energies["edges"]) + sum(energies["faces"]))

step, max_step = 0, 1000
while True:
    evolver.minimize()
    info = evolver.get_minimizer_info()
    print(step, info)
    dump.vtk(f"fire_kb_{kb}_step_{step}")
    energies = computer.compute_mesh_energy(evolver)
    print("vertex energy:", sum(energies["vertices"]), "edge energy:", sum(energies["edges"]), "face energy:", sum(energies["faces"]))
    print("total energy:", sum(energies["vertices"]) + sum(energies["edges"]) + sum(energies["faces"]))

    if (info["Mesh>Fire"]["converge"] == 'true') or step >= max_step: 
        break
    step += 1
dump.vtk("fire_kb_" + str(kb))    
