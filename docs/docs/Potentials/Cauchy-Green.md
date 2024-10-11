---
layout: default
title: Cauchy-Green
parent: Potentials
nav_order: 1
---

# Cauchy-Green
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
  <dd>Mesh>Cauchy-Green</dd>
  <dt>Type</dt>
  <dd>Stretching</dd>
  <dt>Defined in</dt>
  <dd>Faces</dd>
</dl>
---

Stretching energy density is calculated as
<div align="center">
$E=\frac{Y}{2}\frac{1}{1+\nu}\left[Tr\left(C^2\right) + \frac{\nu}{1-\nu}\left(Tr C\right)^2\right]$,
</div>
where $ k $ is the elastic constant, $ \nu $ if the Poisson ratio, and $ C = \frac{1}{2}\left(FQ-I\right) $, with
$ F $ being Gram matrix for strained triangle, $ Q $ being the inverse of the Gram matrix for the unstretched 
triangle and $ I $ is the identity matrix. 




## Python calling

```python
evolver.add_force("Mesh>Cauchy-Green", {"E":{"1":str(E_value)},  #Elastic constant
                                        "h":{"1":str(h_value)},  #Thickness 
                                        "nu":{"1":str(nu_value)} #Poisson ratio
                                       })
```                                
