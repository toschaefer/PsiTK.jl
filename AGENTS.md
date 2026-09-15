# AGENTS.md — AI Assistant Guide

You are an expert Julia developer. Never write "Pythonic" Julia. 
Read the `README.md` and linked documentation before contributing.

## Core Philosophy
* **Performance First:** High-performance is the absolute top priority—even slightly edging out readability. Fast, zero-allocation, and type-stable code comes before all else.
* **Readable & Simple:** Simplicity over feature count. While performance is the ultimate priority, prioritize transparent and hackable implementations over clever brevity.
* **Code as Documentation:** The code itself is the primary documentation. 
  * Keep functions short and trackable.
  * Write docstrings ONLY for the public API or highly non-obvious user-facing functions. **When adding a new public function, you must also append it to `docs/src/code_reference.md`.**
  * Use inline comments ONLY to explain physics/math reasoning, subtle normalizations, or workarounds. Never restate what the code mechanically does.

## Device-Agnostic Code (GPU/CPU)
The exact same code path must run on both CPU and GPU.
* **Allocation:** Never hardcode `zeros(...)` or `CuArray(...)`. Allocate using `similar(X)` or `fill!(similar(X), 0)` so the result inherits the device of `X`.
* **Data Transfer:** Never use vendor-specific constructors like `CuArray(A)`. For instance (since we depend on `DFTK.jl`), leverage its existing infrastructure to move data (e.g., `DFTK.to_device(arch, A)` and `DFTK.to_cpu(A)`)
* **Scalar Indexing:** No scalar indexing into GPU arrays in hot loops. If you must touch a single element, explicitly wrap it in `GPUArraysCore.@allowscalar`.
* **Kernels:** Functions running on the device must take `isbits` arguments. Do not close over non-`isbits` values.

## Style & Conventions
Follow [Julia Blue Style](https://github.com/invenia/BlueStyle). Line length: ~92 characters.

* **Signatures:** Break long function signatures vertically (one argument per line) with a trailing comma.
* **Variables:** Use explicit NamedTuples `(; var=val)`, not `(var=val)`. Prefer readable, descriptive names over abbreviated ones.
* **Loops:** Use `=` for ranges (`for i = 1:10`) and `in` for collections (`for item in array`).
* **Types:** Always use explicit braces for where clauses: `where {T <: AbstractFloat}`.
* **Arguments:** Keyword arguments must be explicit. No implicit positional-to-keyword promotion.
* **Internals & Placeholders:** Prefix internal helpers with `_` (e.g., `_compute_density`). Use `identity` as a placeholder for empty callbacks.

## Units
Use atomic units throughout. Lengths are in Bohr, energies in Hartree.
```julia
using Unitful, UnitfulAtomic
austrip(10u"eV")    # Convert 10 eV → Hartree
auconvert(u"Å", 1.2) # Convert 1.2 Bohr → Ångström

```

## Testing & Quality Assurance

* **Strict Quality Checks:** We use `Aqua.jl` to enforce strict code quality. All new code must be free of type ambiguities, unbound arguments, and pirated methods. Verify against `test/aqua.jl`.
* **Test Structure:** Place tests in the appropriate files within the `test/` directory. If using `@testitem`, ensure the test block is fully self-contained.

## Git & CI Workflow

* **Main branch:** `main` (all PRs target this branch).
* **Workflow:** Fork → Branch → PR to `main`.
* **Releases & Dependencies:** Releases are tagged via `TagBot`. Dependency `[compat]` bounds in `Project.toml` are managed automatically by `Dependabot`.
