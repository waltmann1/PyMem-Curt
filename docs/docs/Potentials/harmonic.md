---
layout: default
title: Harmonic
parent: Potentials
nav_order: 3
---

# Harmonic
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
  <dd>Mesh>Harmonic</dd>
  <dt>Type</dt>
  <dd>Harmonic</dd>
  <dt>Defined in</dt>
  <dd>Edges</dd>
</dl>
---
Calculate harmonic energy of the edge using
<div align="center">
$ E = \frac{1}{2} k \left( l - l_0 \right)^2 $
</div>
where $ k $ is the spring constant, $ l_0 $ is length of unstretched bond (edge) and $ l $ is current edge length.

## Python calling

```python
evolver.add_force("Mesh>Harmonic", {"k":{"1":str(spring_constant_value)},
                                    "l0":{"1":str( unstretched_bond_value)}
                                       })
```                                
