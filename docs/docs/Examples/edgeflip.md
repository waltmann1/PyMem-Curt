---
layout: default
title: Fluid Membranes
parent: Examples
nav_order: 3
---

# Edgeflip elastic shell

{: .no_toc }

In this example we simulate the structure of a fluid membrane by using edge-flip.


<dl>
  <dt>By</dt>
  <dd>Siyu Li</dd>
  <dt>Type</dt>
  <dd>Monte Carlo</dd>
</dl>

---

```python
#First we need to import the module for the simulations
import pymembrane as mb
#numpy
import numpy as np
```

#### Loading Mesh

We first load the icosahedral mesh, the triangulation number used in the paper is T=189 with (h, k)=(12, 3). The input file can be found in the end of this section.

```python
#create a system
system = mb.System()
#read the mesh
Nv = 72
vertex_file = 'vertices_Nv_' + str(Nv) + '.inp'
face_file = 'faces_Nv_' + str(Nv) + '.inp'
system.read_mesh_from_files(files = {'vertices': vertex_file, 'faces': face_file})
```

```python
#save the mesh to display
#create dumper
dump = system.dump()
dump.vtk("initial_mesh")
```

```python
#add the evolver class where the potentials and integrators are added
evolver = mb.Evolver(system)
```

```python
evolver.add_force("Mesh>Harmonic", {"k": {"0": str(10.)},
                                    "l0": {"0": str(1.0)}})

evolver.add_force("Mesh>Limit", {"lmax": {"0": "1.5"},
                                "lmin": {"0": "0.5"}})
```

```python
evolver.add_integrator("Mesh>MonteCarlo>vertex>move", {"dr": "0.05", "spherical_move": "true"})
evolver.add_integrator("Mesh>MonteCarlo>edge>flip", {"every step": "1"})
```

```python
print(evolver.get_force_info())
computer = system.compute_mesh()
energies = computer.compute_mesh_energy(evolver)
print("initial vertex energy:", sum(energies["vertices"]), "edge energy:", sum(energies["edges"]), "face energy:", sum(energies["faces"]))
print("initial total energy:", sum(energies["vertices"]) + sum(energies["edges"]) + sum(energies["faces"]))
```

```python
run_steps = 50000
templist = np.linspace(0.1, 10**-7, 50)
frame = 0

for tempindex, temp in enumerate(templist):
    evolver.set_global_temperature(str(temp))
    acc = evolver.evolveMC(run_steps)
    print(f"frame:{frame}")
    print("vertexmove accepted ratio:", acc['Mesh>MonteCarlo>vertex>move'] / system.Numvertices / run_steps)
    print("edge flip accepted ratio:", acc['Mesh>MonteCarlo>edge>flip'] / system.Numedges / run_steps)
    frame += 1
    dump.vtk(f"edge_flip_Nv_{system.Numvertices}_{frame}")
    energies = computer.compute_mesh_energy(evolver)
    print("vertex energy:", sum(energies["vertices"]), "edge energy:", sum(energies["edges"]), "face energy:", sum(energies["faces"]))
    print("total energy:", sum(energies["vertices"]) + sum(energies["edges"]) + sum(energies["faces"]))

dump.vtk(f"edge_flip_Nv_{system.Numvertices}")
```

Download the [**script**](../../attached/edgeflip.py), [**initial vertice file**](../../attached/vertices_Nv_32.inp), [**initial face file**](../../attached/faces_Nv_32.inp).

#### Visualization in Paraview

![](Nv32.pdf)
