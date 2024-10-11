---
layout: default
title: Disclination
parent: Examples
nav_order: 1
---

# Disclination

{: .no_toc }

This example reproduce the disclination in a pentagon of 'Radius = N'. For more information see: [Seung, H. S., and David R. Nelson. "Defects in flexible membranes with crystalline order." Physical Review A 38.2 (1988): 1005.](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.38.1005). The color code in the images represent the mean curvature of the mesh.


<dl>
  <dt>By</dt>
  <dd>Daniel Matoz Fernandez</dd>
  <dt>Type</dt>
  <dd>Monte Carlo</dd>
</dl>

---

```python
#import the code
import pymembrane as mb
import numpy as np
```

```python
#create a system
system = mb.System()
#read the mesh
N = 10 #pentagon size
vertex_file = 'vertices_N' + str(N) + '.inp'
face_file = 'faces_N' + str(N) + '.inp'
system.read_mesh_from_files(files={'vertices':vertex_file, 'faces':face_file})
```

```python
#save the mesh to display
#create dumper
dump = system.dump()
dump.vtk("initial mesh")
```

![disclination_init](../../images/01_disclination_init.png)

```python
#add the evolver class where the potentials and integrators are added
evolver = mb.Evolver(system)
```

```python
# first we need to know the edge length to move it appropiate:
compute = system.compute_mesh()
edge_lengths = compute.compute_edge_lengths()
avg_edge_length= np.mean(edge_lengths)
print("avg_edge_length = ", avg_edge_length)

#add the potentials in this case we will add bending and stretching
# stretching potential
evolver.add_force("Mesh>Harmonic", {"k":{"0":str(50.0)},
                                    "l0":{"0":str(1.0)}})

# limit potential
evolver.add_force("Mesh>Limit", {"lmin":{"0":str(0.7)},
                                 "lmax":{"0":str(1.3)}})
# bending potential
evolver.add_force("Mesh>Bending", {"kappa":{"1":str(1.0)}})


#edges = system.getEdges()
#for edge in edges:
#    print(edge.type)
```

```python
#add the monte carlo integrator
# vertex move:


evolver.add_integrator("Mesh>MonteCarlo>vertex>move", {"dr":str(0.005),
                                                       "seed":str(123949)})


```

```python
## Now we want to have 100 snapshots every 1000 steps so then
snapshots = 100
run_steps = 2000

## then we want to run the simulation for a temperature 1e-6
evolver.set_global_temperature(str(1e-3))

dump.vtk("pentagon_t0")
for snapshot in range(snapshots):
    acceptance = evolver.evolveMC(MCS=run_steps)
    print(acceptance)
    dump.vtk("pentagon_t" + str(snapshot*run_steps))

```

![01_disclination](../../images/01_disclination.png)

Download the [**script**](../../attached/disclination.py), [**initial vertice file**](../../attached/vertices_N10.inp), [**initial face file**](../../attached/faces_N10.inp).
