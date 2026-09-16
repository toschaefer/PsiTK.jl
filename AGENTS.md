# AGENTS.md — AI Assistant Guide

You are an expert Julia developer. Never write "Pythonic" Julia. 
Read the `README.md` and linked documentation before contributing.

## Core Philosophy
* **Performance First:** High-performance is the absolute top priority—even slightly edging out readability. Fast, zero-allocation, and type-stable code comes before all else.
* **Readable & Simple:** Simplicity over feature count. While performance is the ultimate priority, prioritize transparent and hackable implementations over clever brevity.
* **Code as Documentation:** The code itself is the primary documentation. 
  * Keep functions short and trackable.
  * Write docstrings ONLY for the public API or highly non-obvious user-facing functions. Every exported name needs a docstring (`docs/make.jl` builds with `checkdocs = :exports` and fails otherwise); `docs/src/code_reference.md` picks them up automatically via `@autodocs`, so never list functions there by hand. When one docstring covers several methods (e.g. a convenience form), its argument list must explain every signature.
  * Use inline comments ONLY to explain physics/math reasoning, subtle normalizations, or workarounds. Never restate what the code mechanically does.

## Device-Agnostic Code (GPU/CPU)
**Planned, not yet enforced.** The code is CPU-only; the rules below apply to new kernels and
to explicit GPU-porting tasks. Never convert existing code or add GPU dependencies unasked.
The exact same code path must run on both CPU and GPU.
* **Allocation:** Never hardcode `zeros(...)` or `CuArray(...)`. Allocate using `similar(X)` or `fill!(similar(X), 0)` so the result inherits the device of `X`.
* **Data Transfer:** Never use vendor-specific constructors like `CuArray(A)`. For instance (since we depend on `DFTK.jl`), leverage its existing infrastructure to move data (e.g., `DFTK.to_device(arch, A)` and `DFTK.to_cpu(A)`)
* **Scalar Indexing:** No scalar indexing into GPU arrays in hot loops. If you must touch a single element, explicitly wrap it in `GPUArraysCore.@allowscalar`.
* **Kernels:** Functions running on the device must take `isbits` arguments. Do not close over non-`isbits` values.

## Style & Conventions
Follow [Julia Blue Style](https://github.com/invenia/BlueStyle). Line length: ~92 characters.


## Testing & Quality Assurance

* **Strict Quality Checks:** We use `Aqua.jl` to enforce strict code quality. All new code must be free of type ambiguities, unbound arguments, and pirated methods. Verify against `test/aqua.jl`.
* **Test Structure:** Place tests in the appropriate files within the `test/` directory. A `@testitem` must not rely on state from other items; shared fixtures go into `@testmodule`s (e.g. `TestSystems` in `test/systems.jl`, via `setup=[TestSystems]`). Tag items (`tags=[:name]`) so they can be selected with `Pkg.test(test_args=["name"])`.

## Git & CI Workflow

* **Main branch:** `main` (all PRs target this branch).
* **Workflow:** Fork → Branch → PR to `main`.
* **Releases & Dependencies:** Releases are tagged via `TagBot`. Dependency `[compat]` bounds in `Project.toml` are managed automatically by `Dependabot`.
