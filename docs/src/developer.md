```@meta
CurrentModule = PsiTK
```

# Developer Guide

This page is a map of the code: it names the parts, shows how they fit together and states
the conventions shared by all of them. The details (signatures, formulas, costs) live in
the docstrings, collected in the [Code Reference](@ref), and the code itself remains the
authoritative documentation. When in doubt, follow the nearest existing example.

## Nouns and verbs

PsiTK builds on [DFTK](https://docs.dftk.org)'s plane-wave infrastructure (e.g.
`PlaneWaveBasis`, FFTs, Hamiltonian, exchange operators) and borrows its pattern: *nouns*
are plain structs that carry either data or configuration, rather than both; *verbs* are
functions that take the data they need as arguments and the configuration as a dispatch
argument.

| noun | role | lives in |
|---|---|---|
| [`OrbitalSpace`](@ref) | data: orbitals, energies, occupation on a basis; what orbital-space verbs and interfaces consume and produce | `src/orbital_spaces.jl` |
| [`DensityFitting`](@ref) | data: the Coulomb vertex `Γ` together with its auxiliary basis (G vectors, kernel, singular vectors) | `src/coulomb_vertex.jl` |
| target, e.g. [`DensitySpecificVirtuals`](@ref) | physics: *which* eigenvalue problem defines the virtual orbitals | `src/virtual_orbitals.jl` |
| solver, e.g. [`LOBPCG`](@ref) | numerics: *how* an eigenvalue problem is solved | `src/eigensolvers/` |
| strategy, e.g. [`CoulombGramian`](@ref) | numerics: *how* the Coulomb vertex is compressed | `src/coulomb_vertex.jl` |
| operator, e.g. [`LevelShiftedOperator`](@ref) | numerics: matrix-free wrapper of a DFTK operator, applied with `mul!` | `src/operators.jl` |
| callback, e.g. [`ShowProgress`](@ref) | output: reports the progress of long computations | `src/callbacks.jl` |

The verbs transform the two data nouns into each other:

| verb | transforms | lives in |
|---|---|---|
| [`OrbitalSpace`](@ref)`(scfres)` | DFTK SCF result → `OrbitalSpace` | `src/orbital_spaces.jl` |
| [`split_occupied_virtual`](@ref), [`select_orbitals`](@ref), [`merge_spaces`](@ref), [`canonicalize_orbitals`](@ref) | `OrbitalSpace` → `OrbitalSpace` | `src/orbital_spaces.jl` |
| [`generate_orbitals`](@ref) | `OrbitalSpace` (occupied) → `OrbitalSpace` (virtual) | `src/virtual_orbitals.jl` |
| [`compute_overlap_densities`](@ref), [`compute_coulomb_vertex`](@ref) | `OrbitalSpace` → `DensityFitting` | `src/coulomb_vertex.jl` |
| [`compress_coulomb_vertex`](@ref) | `DensityFitting` → `DensityFitting` | `src/coulomb_vertex.jl` |
| [`compute_delta_integrals`](@ref) | orbitals → grid integrals for the interfaces | `src/delta_integrals.jl` |
| [`dump_cc4s_files`](@ref) | `OrbitalSpace` + `DensityFitting` → files | `src/interfaces/` |

A typical workflow composes the verbs like this:

```julia
occ_space, _ = split_occupied_virtual(OrbitalSpace(scfres))              # DFTK enters here
dsv_space    = generate_orbitals(DensitySpecificVirtuals(n_orbitals=50), occ_space, ham)
active_space = canonicalize_orbitals(merge_spaces(occ_space, dsv_space), ham)
fitting      = compress_coulomb_vertex(compute_coulomb_vertex(active_space), CoulombGramian())
dump_cc4s_files(active_space, fitting)                                    # src/interfaces/
```

`src/PsiTK.jl` includes the files in dependency order and declares the exports next to the
`include` of the file defining them. Every export needs a docstring (`docs/make.jl` builds
with `checkdocs = :exports`).

## Conventions

- **k-points:** per-k-point data are `Vector`s over k-points, e.g. `ψ[ik]` is an
  ``N_\text{pw} \times N`` matrix of orbital coefficients and `eigenvalues[ik]` a vector.
  Verbs loop over k-points; the Coulomb vertex and the interfaces are Γ-point only so far.
- **Vertex layout:** the Coulomb vertex `Γ` has the shape `(Nk, Nbra, Nk, Nket, NG)`,
  or `(Nk, Nbra, Nk, Nket, NF)` once compressed.
- **Orthonormality:** an `OrbitalSpace` records whether its orbitals are orthonormal
  (`is_orthonormal`). Orbitals of a generalized eigenvalue problem, e.g.
  [`DensitySpecificVirtuals`](@ref), are not; [`canonicalize_orbitals`](@ref) makes them
  orthonormal, and the interfaces require it.
- **Callbacks:** long computations take a `callback` keyword, called with a `NamedTuple`
  after each step; `callback=identity` silences it. Step-wise loops pass
  `(; step, total_steps)`, which [`ShowProgress`](@ref) turns into a progress bar.

## Extension points

Extending the code usually means adding one noun and its methods of an extension point.
The docstring of each extension point states its contract.

### Virtual-orbital targets and eigensolvers

[`generate_orbitals`](@ref) splits physics from numerics, so that targets and solvers
combine freely. Per k-point it

1. takes the eigenvalue problem `(; A, B, ε_offset)` of the target from
   [`_eigenproblems`](@ref). The operators are matrix-free wrappers of DFTK operators
   ([`LevelShiftedOperator`](@ref), [`ProjectedShiftedOperator`](@ref)) that push the
   occupied space out of the low end of the spectrum;
2. builds a stochastic guess `X0` of [`_n_guess`](@ref)`(solver, n)` vectors orthogonal to
   the occupied space, where `n` is the number of requested orbitals;
3. solves for the lowest `n` eigenpairs with [`_solve`](@ref)`(A, B, X0, n, P, solver)`;
4. returns an `OrbitalSpace`, orthonormal or not as given by [`_is_orthonormal`](@ref).

- **A new target** is a configuration struct with an `n_orbitals` field, added to the
  [`VirtualOrbitalTarget`](@ref) union, plus methods of [`_eigenproblems`](@ref) and
  [`_is_orthonormal`](@ref) in `src/virtual_orbitals.jl`.
- **A new eigensolver** is a configuration struct plus a method of [`_solve`](@ref)
  (and optionally of [`_n_guess`](@ref)) in `src/eigensolvers/`.

### Compression strategies

[`compress_coulomb_vertex`](@ref) delegates to [`_compress_coulomb_vertex`](@ref) and
accumulates the transformations of repeated compressions in the
[`DensityFitting`](@ref). **A new strategy** is a configuration struct plus a method of
[`_compress_coulomb_vertex`](@ref) in `src/coulomb_vertex.jl`.

### Correlation-solver interfaces

**A new interface** is a file in `src/interfaces/` consuming an [`OrbitalSpace`](@ref) and a
[`DensityFitting`](@ref), so that only that file knows the solver's file format. It checks
its preconditions (e.g. orthonormal canonical orbitals) itself, see
[`dump_cc4s_files`](@ref).

## Computational cost

A few functions dominate the runtime and memory of a calculation once there are a few
hundred virtual orbitals. To keep their cost visible, every public function and every
extension point (target, solver, compression strategy) whose time or memory grows faster
than linearly with the system size, i.e. with the number of orbitals or plane waves (log
factors from FFTs do not count), documents its cost in a `# Cost` section at the end of
its docstring:

- **time:** the leading terms in O-notation and the step they come from;
- **memory:** the peak memory in complex numbers (16 bytes each).

A change of the scaling updates this section. Following the split between physics and
numerics, a target documents the cost of one application of its operators, and a solver
its cost in terms of operator applications. For example, for `compute_coulomb_vertex`:

```
# Cost
Those of [`compute_overlap_densities`](@ref), which the vertex overwrites in place:
- time: ``O(N_\text{bra} N_\text{ket} N_r \log N_r)``
- memory: ``N_\text{bra} N_\text{ket} N_G`` complex numbers for `Γ`
```

Counts use the following symbols. In code, the name is the symbol without the underscore
of its subscript (``N_G`` ↔ `NG`); counts without a symbol are snake_case (`n_orbitals`).
Each `# Cost` section states which of its arguments the symbols refer to.

| docstring | code | meaning |
|---|---|---|
| ``N_\text{occ}``, ``N_\text{virt}``, ``N`` | `Nocc`, `Nvirt`, `N` | occupied, virtual and all orbitals |
| ``N_\text{bra}``, ``N_\text{ket}`` | `Nbra`, `Nket` | orbitals on the bra and ket side of a pair |
| ``N_k`` | `Nk` | k-points |
| ``N_\text{pw}`` | `Npw` | plane waves of the orbital basis |
| ``N_G`` | `NG` | plane waves of the overlap densities, set by `Ecut_ratio` |
| ``N_r`` | `Nr` | points of the FFT grid |
| ``N_F`` | `NF` | compressed auxiliary index of the Coulomb vertex |
| ``N_{pp}`` | `Npp` | orbital pairs of the Coulomb vertex, ``N_k^2 N_\text{bra} N_\text{ket}`` |
