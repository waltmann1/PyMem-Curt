---
layout: default
title: Home
nav_order: 1
description: "Continuum simulation of elastic membranes."
permalink: /
---

![](/images/logo-white.png)


# Continuum simulation of elastic membranes
{: .fs-9 }

PyMemebrane is **python** package for the simulation of elasticity in thin membranes. The code is written in **C++** and extensively uses [Standard Template Library](https://en.wikipedia.org/wiki/Standard_Template_Library) and [PyBind11](https://pybind11.readthedocs.io/en/stable/index.html) which provides a modular object-oriented design allowing simple extensions and maintenance.

{: .fs-6 .fw-300 }

[Get started now](#getting-started){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 } [View it on GitHub](https://github.com/fdmatoz/pymembrane){: .btn .fs-5 .mb-4 .mb-md-0 }

---

## Getting started

### Dependencies

We recommend using Anaconda and creating a PYMEMB environment. To create a new conda environment type:

```bash
conda create -n PYMEMB python=3.7 numpy scipy matplotlib jupyter jupyterlab vtk pip
```

This should install all packages needed to run the code. Otherwise, in the existing conda environment (PYMEMB) please install these additional packages:

* cmake
* ipympl
* nodejs
* VTK

```bash
conda install -c anaconda cmake
```
```bash
conda install -c conda-forge ipympl
```
```bash
conda install -c conda-forge nodejs
```
Please also install the following extensions for **JupyterLab**:

```bash
jupyter labextension install @jupyter-widgets/jupyterlab-manager jupyter-matplotlib
```

In addition we assume that you have a working C/C++ compiler that it supports the **C++14** or newer standard and [vtk](https://vtk.org/download/) installed*.

<small>*Note: You can installed from source or follow the recommendations of your distribution.</small>

We also recommend that you install [Paraview](https://www.paraview.org/) which will allow you to visualize the mesh files and their attributes.


### Local installation using python

In order to make the Python modules visible, please install the pymembrane module. From the **``pymembrane``**  directory type:

``python setup.py install``

## Quick start: Desclination problem

Please refer to the [basic example](./docs/Examples/desclination.md).

---

## About the project

PyMemebrane is &copy; 2019-{{ "now" | date: "%Y" }} by [Daniel Matoz Fernandez](http://www.danielmatoz.com).

### License

PyMemebrane is distributed by an [MIT license](LICENSE).

