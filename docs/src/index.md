# PsiTK.jl

PsiTK is a Julia package that serves as a bridge between the plane wave-based [Density-Functional Toolkit (DFTK.jl)](https://github.com/JuliaMolSim/DFTK.jl) and electron correlation solvers, such as the massively parallel coupled cluster code [Cc4s](https://gitlab.cc4s.org/cc4s/cc4s).

<p align="left"> <img src="/assets/workflow.png" alt="workflow" height="80"> </p>

Building on DFTK's infrastructure for self-consistent field calculations of 3D models under periodic boundary conditions, PsiTK provides the necessary capabilities to prepare electron correlation methods and interface data to the solver.
This encompasses a range of functions, including the definition and compression of virtual orbital manifolds and density fitting for the calculation of electron-repulsion integrals (Coulomb vertex).

Furthermore, PsiTK natively includes algorithms for MP2 and RPA calculations.


## API Reference

```@autodocs
Modules = [PsiTK]
```
