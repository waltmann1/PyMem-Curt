---
layout: default
title: Constant Face Area
parent: Potentials
nav_order: 4
---

# Constant Face Area 
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
  <dd>Mesh>Constant Area</dd>
  <dt>Type</dt>
  <dd>Stretching</dd>
  <dt>Defined in</dt>
  <dd>Faces</dd>
</dl>
---

For a given triangle t compute contribution from the surface tension caused by change in the surface area. 
<div align="center">
    $ E = \frac{1}{2} \sigma \left(A - A_0\right)^2 $
</div>
where $ A $ is current area of the triangle, $ A_0 $  is the native (given) triangle area, and $ \sigma $ if the surface tension.




## Python calling

```python
evolver.add_force("Mesh>Constant Area", {"sigma":{"1":str(E_value)}
                                       })
```                                
