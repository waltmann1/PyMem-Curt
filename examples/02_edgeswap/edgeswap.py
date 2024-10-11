#First we need to import the module for the simulations
import pymembrane as mb
#numpy
import numpy as np


def run(fac):
    #create a system 
    system = mb.System()
    #read the mesh
    vertex_file = 'vertices_T189.inp'
    face_file = 'faces_T189.inp'
    system.read_mesh_from_files(files = {'vertices': vertex_file, 'faces': face_file})
    # fac = 0.7
    edges = system.getEdges()
    Ne = len(edges)
    id0s = np.random.choice(range(Ne),round(fac * Ne), replace = False)
    id1s = list(set(range(Ne)) - set(id0s))
    for li in id1s: 
        edges[li].type = 1
    for li in id0s: 
        edges[li].type = 0
    system.setEdges(edges)

    #create dumper
    dump = system.dump() 
    dump.vtk("initial mesh")

    #add the evolver class where the potentials and integrators are added
    evolver = mb.Evolver(system)

    evolver.add_force("Mesh>Harmonic", {"k":{"0": "5.77", "1": "5.77"}, 
                                        "l0":{"0": "1.0", "1": "1.0"}})

    evolver.add_force("Mesh>Bending", {"kappa":{"0": "0.06", "1": "28.9"}})

    evolver.add_force("Mesh>Limit", {"lmax":{"0": "1.3", "1": "1.3"},
                                    "lmin":{"0": "0.7", "1": "0.7"},})

    evolver.add_integrator("Mesh>MonteCarlo>vertex>move", {"dr":"0.05"})
    evolver.add_integrator("Mesh>MonteCarlo>edge>swap", {"every step": "1"})

    print(evolver.get_force_info())

    computer = system.compute_mesh()
    energies = computer.compute_mesh_energy(evolver)
    print("initial total energy:", sum(energies["vertices"]) + sum(energies["edges"]) + sum(energies["faces"]))

    ## Now we want to run 100000 steps and take snapshot every 10000 steps so then
    frame = 0
    run_steps = 20000
    anneal_steps = 5
    templist = np.linspace(0.1, 10**-7, 10)

    for anneal in range(anneal_steps):
        for tempindex, temp in enumerate(templist):
            evolver.set_global_temperature(str(temp))
            frame += 1
            acc = evolver.evolveMC(run_steps)
            print("edge swap accepted ratio:", acc['Mesh>MonteCarlo>edge>swap'] / system.Numedges / run_steps)
            print("vertexmove accepted ratio:", acc['Mesh>MonteCarlo>vertex>move'] / system.Numvertices / run_steps)
            energies = computer.compute_mesh_energy(evolver)
            print("frame:", frame, "vertex energy:", sum(energies["vertices"]), "edge energy:", sum(energies["edges"]), "face energy:", sum(energies["faces"]))
            # dump.edge_vtk("edge_swap_frac_" + str(frame))

    dump.edge_vtk("edge_swap_fraction_" + str(fac) + "_" + str(frame))

# for f in [0.1, 0.2, 0.5, 0.8]:
#     run(f)

for f in [0., 0.015, 0.023, 0.025, 0.035, 0.045, 0.01, 0.02, 0.03, 0.05, 0.15, 0.3, 1.]:
    run(f)