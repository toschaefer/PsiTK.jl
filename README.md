<p align="left"> <img src="docs/src/logo/PsiTK+name.png" alt="PostDFTK logo" height="100"> </p>

# Wavefunction Toolkit

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://toschaefer.github.io/PsiTK.jl/)
[![Build Status](https://github.com/toschaefer/PsiTK.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/toschaefer/PsiTK.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Status](https://img.shields.io/badge/Status-Experimental-red)](https://github.com/toschaefer/PsiTK.jl)



**PsiTK** is a lightweight framework for developing correlated wavefunction workflows in the plane wave basis.

Built on top of [DFTK.jl](https://github.com/JuliaMolSim/DFTK.jl), it acts as the wavefunction-based counterpart for post-Hartree–Fock and beyond-DFT methods.

Features:

* virtual orbital space compression (e.g. natural orbitals)
* MP2 ground-state energy
* interface to the high-performance coupled cluster (CC) solver [Cc4s](https://gitlab.cc4s.org/cc4s/cc4s), enabling CC calculations via the workflow DFTK &rarr; PsiTK &rarr; Cc4s.

## Motivation

While established electronic codes (VASP, Quantum Espresso, ABINIT, CP2K, ...) are highly optimized for production calculations, their sheer scale can make implementing novel theoretical ideas challenging. 
PsiTK fills this gap by minimizing the distance between scientific intuition and executable code, offering a flexible testbed without sacrificing readability or efficiency.

## Usage (TOD)

```julia
using DFTK
using PseudoPotentialData
using PsiTK

pd_pbe_family = PseudoFamily("dojo.nc.sr.pbe.v0_5.stringent.upf") 
Si = ElementPsp(:Si, pd_pbe_family)
a = 10.26  # Silicon lattice constant in Bohr
lattice = a / 2 * [[0 1 1.];
                   [1 0 1.];
                   [1 1 0.]]
atoms     = [Si, Si]
positions = [ones(3)/8, -ones(3)/8]

model  = model_HF(lattice, atoms, positions)
basis  = PlaneWaveBasis(model; basis.Ecut, basis.kgrid)
scfres = self_consistent_field(basis)
```
