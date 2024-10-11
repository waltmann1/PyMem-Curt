---
layout: default
title: Limit
parent: Potentials
nav_order: 4
---

# Limit
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
  <dd>Limit</dd>
  <dt>Defined in</dt>
  <dd>Edges</dd>
</dl>
---
If edge is longer or shorter than a given values, energy is infinity, otherwise it is zero,
<div align="center">
\begin{equation}
E =
\begin{cases}
\infty\qquad \text{if} \qquad r < l_{min}\\
\infty\qquad \text{if} \qquad r > l_{max}\\
0 \qquad \text{otherwise}
\end{cases}.
\end{equation}
</div>

## Python calling

```python
evolver.add_force("Mesh>Limit", {"lmin":{"1":str(min_constant_value)},
                                    "lmax":{"1":str(max_bond_value)}
                                       })
```                                
