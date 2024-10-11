#First we need to import the module for the simulations
import pymembrane as mb
import numpy as np

def run(fac0, fac1):
    print(fac0, fac1)
    #create a system 
    system = mb.System()
    #read the mesh
    vertex_file = 'vertices_T192.inp'
    face_file = 'faces_T192.inp'
    system.read_mesh_from_files(files = {'vertices': vertex_file, 'faces': face_file})

    vertices = system.getVertices()
    Nv=len(vertices)
    id0s = np.random.choice(range(Nv),round(fac0 * Nv),replace = False)
    id1s = np.random.choice(list(set(range(Nv))-set(id0s)), round(fac1 * Nv), replace = False)
    id2s = list(set(range(Nv)) - set(id0s) - set(id1s))

    for vi in id0s: vertices[vi].type = 0
    for vi in id1s: vertices[vi].type = 1
    for vi in id2s: vertices[vi].type = 2

    system.setVertices(vertices)

    #create dumper
    dump = system.dump() 
    dump.vtk("initial mesh")

    #add the evolver class where the potentials and integrators are added
    evolver = mb.Evolver(system)

    evolver.add_force("Mesh>Harmonic", {"k": {"0": str(10.), "1": str(10.), "2": str(10.)},
                                        "l0": {"0": str(1.0), "1": str(1.0), "2": str(1.0)}})

    evolver.add_force("Mesh>BendingGK", {"H0": {"0": "0.", "1": "0.", "2": "0."},
                                        "kappaH": {"0": "4.", "1": "1.", "2": "0.8"},
                                        "kappaG": {"0": str(-0.66667 * 4.), "1": str(-0.66667 * 1.), "2": str(-0.66667 * 0.8)}})

    evolver.add_force("Mesh>Line Tension", {"gamma": {"0": "0.0004", "1": "0.0004", "2": "0.0004"}}) 

    evolver.add_force("Mesh>Limit", {"lmax": {"0": "1.3", "1": "1.3"},
                                    "lmin": {"0": "0.7", "1": "0.7"}})

    evolver.add_integrator("Mesh>MonteCarlo>vertex>move", {"dr": "0.05"})
    evolver.add_integrator("Mesh>MonteCarlo>vertex>swap", {"every step": "1"})

    print(evolver.get_force_info())

    computer = system.compute_mesh()
    energies = computer.compute_mesh_energy(evolver)
    print("initial total energy:", sum(energies["vertices"]) + sum(energies["edges"]) + sum(energies["faces"]))

    ## Now we want to run 20000 steps and take snapshot every 10000 steps so then
    frame = 0
    run_steps = 20
    templist1 = [0.01 * i for i in range(10, 0, -1)] + [10**-7]
    templist2 = [0.001 * i for i in range(10, 0, -1)] + [10**-7]
    templist = templist1 * 5 + templist2

    for tempindex, temp in enumerate(templist):
        ## we step the temperature first
        evolver.set_global_temperature(str(temp))
        frame += 1
        acc = evolver.evolveMC(run_steps)
        print("vertex move accepted ratio:", acc['Mesh>MonteCarlo>vertex>move'] / system.Numvertices / run_steps)
        print("vertex swap accepted ratio:", acc['Mesh>MonteCarlo>vertex>swap'] / system.Numvertices / run_steps)
        energies = computer.compute_mesh_energy(evolver)
        Ev, Ee, Ef = sum(energies["vertices"]), sum(energies["edges"]), sum(energies["faces"])
        print(f"frame:{frame}, vertex energy: {Ev}, edge energy: {Ee}, face energy:, {Ef}")
        dump.vtk("vertex_swap_" + str(frame))

    dump.vtk(f"vertex_swap_fraction_({fac0}, {fac1}, {1-fac0-fac1:.1f})")

for f0, f1 in [(0.1, 0.8), (0.2, 0.7), (0.4, 0.5), (0.5, 0.4), (0.6, 0.3), (0.8, 0.1)]:
    run(f0, f1)