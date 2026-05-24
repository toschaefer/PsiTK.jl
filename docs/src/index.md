# PsiTK.jl

PsiTK is a Julia package that serves as a bridge between the plane wave-based [Density-Functional Toolkit (DFTK.jl)](https://docs.dftk.org) and electron correlation solvers, such as the massively parallel coupled cluster code [Cc4s](https://gitlab.cc4s.org/cc4s/cc4s).

![workflow](assets/workflow.png)

Building on DFTK's infrastructure for self-consistent field calculations of 3D models under periodic boundary conditions, PsiTK provides the necessary capabilities to prepare electron correlation methods and interface data to the solver.
This encompasses a range of functions, including the definition and compression of virtual orbital manifolds and density fitting for the calculation of electron-repulsion integrals (Coulomb vertex).

## Installation

PsiTK is currently unregistered. You can install it directly from GitHub using Julia's package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/toschaefer/PsiTK.jl")
```

## Quick Start / Tutorials

To see PsiTK in action—from running a Hartree-Fock calculation to generating Density Specific Virtuals and computing the Coulomb vertex—we recommend looking at the `examples/` directory in the repository. A complete, runnable example for a solid Helium system is available at `examples/helium_cc4s/helium.jl`.

## Testing

To run the full test suite, enter the Julia package manager (`]`) and run:
```julia
pkg> test
```
To filter and run specific test items, pass the corresponding tags as test arguments:
```bash
julia --project -e 'using Pkg; Pkg.test(test_args=["coulomb_vertex"])'
```
