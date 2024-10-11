---
layout: default
title: Line Tension
parent: Potentials
nav_order: 4
---

# Line Tension
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
  <dd>Mesh>Line Tension</dd>
  <dt>Type</dt>
  <dd>Limit</dd>
  <dt>Defined in</dt>
  <dd>Vertices</dd>
</dl>
---
Iterate over all neighbors of a given vertex. If a current neighbor is of a different type add gamma to the line tension energy. Since in principle we can associate different $\gamma$ to each vertex  we use a simple algebraic men to determine the line tension between vertices with different local gamma. This will be useful only if we have a three or more components system.

## Python calling

```python
evolver.add_force("Mesh>Limit", {"gamma":{"1":str(gamma_value)}})
```                                
