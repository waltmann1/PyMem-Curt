---
layout: default
title: Bending
parent: Potentials
nav_order: 2
---

# Bending 
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Definition
---
<dl>
  <dt>Name</dt>
  <dd>Mesh>Bending</dd>
  <dt>Type</dt>
  <dd>Bending</dd>
  <dt>Defined in</dt>
  <dd>Edges</dd>
</dl>
---
Use normals that belong to the two traingles that share the edge and calculate bengind energy as

<div align="center">
$E = \frac{1}{2} \kappa \vert \hat{n}_k- \hat{n}_l \vert^2$, 
</div>

where $ \kappa $ is the bending rigidity and $ \hat{n}_i, \hat{n}_j $ are the normal unit vectors of the triangles meeting at the edge.


## Python calling

```python
evolver.add_force("Mesh>Bending", {"kappa":{"1":str(kappa_value)}
                                       })
```                                
