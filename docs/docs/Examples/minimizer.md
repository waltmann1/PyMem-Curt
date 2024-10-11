---
layout: default
title: Energy minimization
parent: Examples
nav_order: 3
---

# Energy minimization of a closed shell

{: .no_toc }

In this example we reproduce the energy minimization of an icosahedral shell [Lidmar, et al](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.68.051910)

## Table of contents

<dl>
  <dt>By</dt>
  <dd>Siyu Li</dd>
  <dt>Type</dt>
  <dd>Minimizer</dd>
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
vertex_file = 'vertices_T2500.inp'
face_file = 'faces_T2500.inp'
system.read_mesh_from_files(files = {'vertices': vertex_file, 'faces': face_file})
```

```python
dump = system.dump()
dump.vtk("initial mesh")
```

```python
kb = 0.01
print("kb:", kb)
```

    kb: 0.01

```python
evolver = mb.Evolver(system)

evolver.add_force("Mesh>Harmonic", {"k": {"0": str(1.)},
                                    "l0": {"0": str(1.0)}})

evolver.add_force("Mesh>Bending", {"kappa": {"0": str(kb)}})
```

After assign all the potentials to the system we now need to add the **_Monte Carlo_** integrator. In the following example, each Monte Carlo run will contain $N_v$ attempts move of vertice and $N_v$ attempts of edge swaps.

```python
evolver.add_minimizer("Mesh>Fire", {"dt": "0.1", "max_iter": "10000", "ftol": "0.01", "etol": "0.0000001"})
```

```python
print(evolver.get_force_info())
print(evolver.get_minimizer_info())
computer = system.compute_mesh()
energies = computer.compute_mesh_energy(evolver)
print("initial total energy:", sum(energies["vertices"]) + sum(energies["edges"]) + sum(energies["faces"]))
```

    [{'kappa': '0.01 1 1 1 1 1 1 1 1 1 ', 'name': 'bending', 'type': 'edge'}, {'k': '1 0 0 0 0 0 0 0 0 0 ', 'l0': '1 1 1 1 1 1 1 1 1 1 ', 'name': 'harmonic', 'type': 'edge'}]
    {'Mesh>Fire': {'alpha': '0.1', 'converge': 'false', 'dE': '0', 'dt': '0.01', 'energy': '0', 'etol': '1e-07', 'fnorm': '0', 'ftol': '0.01', 'max_iter': '10', 'n_since_start': '0', 'name': 'fire'}}
    initial total energy: 4.410563588330755

Finally we use annealing method to integrate the monte carlo simulation. Five cooling cycles are used below where each cycle containes a linear decreasing of temperature.

```python
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
```

Download the [**script**](../../attached/minimizer.py), [**initial vertice file**](../../attached/vertices_T2500.inp), [**initial face file**](../../attached/faces_T2500.inp).

#### Visualization in Paraview

![twocomp_fac](../../images/twocomp_fac.pdf)