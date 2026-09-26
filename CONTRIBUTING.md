# Contributing

Contributions via pull request are welcome. Please:

- Add tests for new functionality as `@testitem`s with a tag (see `test/runtests.jl`).
- Add a docstring to every exported name; the documentation build fails otherwise.
- Follow the existing code style ([Blue Style](https://github.com/invenia/BlueStyle),
  lines of about 92 characters).
- Keep pull requests focused on a single change or feature.

Questions? Open an issue.

## Code design

PsiTK builds on DFTK's plane-wave infrastructure (`PlaneWaveBasis`, FFTs, Hamiltonian and
exchange operators) and follows its pattern: *nouns* are plain structs that carry either data
or configuration, never both; *verbs* are functions that take the data as arguments and the
configuration as a dispatch argument.

| noun | role | lives in |
|---|---|---|
| `OrbitalSpace` | data: orbitals, energies, occupations on a basis; the object every verb consumes and produces | `src/orbital_spaces.jl` |
| `DensityFitting` | data: the Coulomb vertex `Γ` together with its auxiliary basis (G vectors, kernel, singular vectors) | `src/coulomb_vertex.jl` |
| target, e.g. `DensitySpecificVirtuals` | physics: *which* eigenvalue problem defines the virtual orbitals | `src/virtual_orbitals.jl` |
| solver, e.g. `LOBPCG` | numerics: *how* an eigenvalue problem is solved | `src/eigensolvers/` |
| strategy, e.g. `CoulombGramian` | numerics: *how* the Coulomb vertex is compressed | `src/coulomb_vertex.jl` |

The verbs compose into the workflow

```julia
occ_space, _ = split_occupied_virtual(OrbitalSpace(scfres))              # DFTK enters here only
dsv_space    = generate_orbitals(DensitySpecificVirtuals(n_orbitals=50), occ_space, ham)
active_space = canonicalize_orbitals(merge_spaces(occ_space, dsv_space), ham)
fitting      = compress_coulomb_vertex(compute_coulomb_vertex(active_space), CoulombGramian())
dump_cc4s_files(active_space, fitting)                                    # src/interfaces/
```

`generate_orbitals` separates physics from numerics: `_eigenproblems(target, occ_space, ham)`
returns the operators `(; A, B, ε_offset)` per k-point, `_solve(A, B, X0, P, solver)` returns
the lowest eigenpairs, so every target works with every solver.

To add

- **a virtual-orbital target:** a configuration struct with an `n_orbitals` field, added to
  the `VirtualOrbitalTarget` union, plus `_eigenproblems` and `_is_orthonormal` methods in
  `src/virtual_orbitals.jl`;
- **an eigensolver:** a configuration struct plus a `_solve` method in `src/eigensolvers/`;
- **a compression strategy:** a struct plus `_compress_coulomb_vertex(Γ, strategy) -> (Γ_F, U)`
  in `src/coulomb_vertex.jl`;
- **a correlation-solver interface:** a file in `src/interfaces/` consuming an `OrbitalSpace`
  and a `DensityFitting`; nothing else knows a solver's file format.
