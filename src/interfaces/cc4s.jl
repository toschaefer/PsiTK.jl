# This file initially based on code of the experimental "cc4s" branch in DFTK
# written by Michael Herbst

"""
    write_cc4s_tensor(folder, name, tensor_data; kwargs...)

Generic function to write a Tensor object in Cc4s format, including its YAML
metadata file and its binary/text elements file.
"""
function write_cc4s_tensor(
    folder::AbstractString,
    name::AbstractString,
    tensor_data;
    dimensions::Vector{<:Dict},
    scalarType::String="Real64",
    elementType::String="IeeeBinaryFile",
    metaData::Dict=Dict{String,Any}(),
    unit::Float64=1.0,
)
    @assert elementType in ("TextFile", "IeeeBinaryFile")
    yamlfile, elementsfile = _cc4s_files(folder, name)

    data = Dict(
        "version" => 100,
        "type" => "Tensor",
        "scalarType" => scalarType,
        "dimensions" => dimensions,
        "elements" => Dict("type" => elementType),
        "unit" => unit,
        "metaData" => metaData,
    )
    open(fp -> YAML.write(fp, data), yamlfile, "w")

    open(elementsfile, "w") do fp
        if elementType == "TextFile"
            for val in tensor_data
                println(fp, val)
            end
        else
            for val in tensor_data
                write(fp, val)
            end
        end
    end

    return [yamlfile, elementsfile]
end

# Write EigenEnergies.yaml and EigenEnergies.elements
function write_eigenenergies(
    folder::AbstractString,
    eigenvalues::AbstractVector,
    εF::Number,
)
    @assert length(eigenvalues) == 1
    εk = eigenvalues[1]

    # Note: Eigenenergies need to be ordered in *non-decreasing* order for cc4s !
    @assert maximum(abs, sort(εk) - εk) < 1e-10

    dimensions = [Dict("length" => length(εk), "type" => "State")]
    metaData = Dict("fermiEnergy" => εF, "energies" => εk)

    return write_cc4s_tensor(
        folder,
        "EigenEnergies",
        εk;
        dimensions,
        scalarType="Real64",
        elementType="TextFile",
        metaData,
    )
end

@doc raw"""
    dump_cc4s_files(
        active_space::OrbitalSpace,
        fitting::DensityFitting;
        folder::AbstractString=joinpath(pwd(), "cc4s"),
        force=false
    )

Write Cc4s input files (*.yaml and *.elements):
- EigenEnergies
- CoulombVertex
- DeltaIntegralsHH
- DeltaIntegralsPPHH
- GridVectors
- CoulombPotential
- CoulombVertexSingularVectors (only if `fitting` is compressed)

Requires a single k-point at Γ and orthonormal canonical orbitals (see
[`canonicalize_orbitals`](@ref)) with eigenvalues in non-decreasing order and integer
occupations. Throws an `ArgumentError` otherwise, or if a file to write exists and `force`
is false, and a `DimensionMismatch` if `fitting` does not match `active_space`; all before
anything is computed or written.

# Arguments
- `active_space`: the `OrbitalSpace` containing the bands (occupation and eigenvalues)
- `fitting`: the [`DensityFitting`](@ref) of the active space, uncompressed from
  [`compute_coulomb_vertex`](@ref) or compressed by [`compress_coulomb_vertex`](@ref)
- `folder`: the target folder (created if missing)
- `force`: if true existing files will be overwritten

# Returns
The list of written file paths.

# Cost
Dominated by the Delta integrals, see [`compute_delta_integrals`](@ref):
- time: ``O(N_r N_\text{virt}^2 N_\text{occ}^2)``
- memory: ``N_r (N_\text{virt}^2 + N_\text{occ}^2)`` complex numbers

The files hold ``N^2 N_F`` complex numbers for the Coulomb vertex.
"""
function dump_cc4s_files(
    active_space::OrbitalSpace,
    fitting::DensityFitting;
    folder::AbstractString=joinpath(pwd(), "cc4s"),
    force=false,
)
    _check_cc4s_input(active_space, fitting, folder, force)
    mkpath(folder)

    # --- dump Eigenvalues
    # For cc4s we just pass the eigenvalues from the active space
    basis = active_space.basis
    eigenvalues = active_space.eigenvalues
    εF = active_space.εF  # Fermi level

    files_ene = write_eigenenergies(folder, eigenvalues, εF)

    # --- dump Coulomb Vertex
    files_coul = write_coulomb_vertex(folder, fitting.Γ)

    # --- Split space by Fermi Energy for Delta Integrals
    idx_holes = findall(ε -> ε <= εF, eigenvalues[1])
    idx_parts = findall(ε -> ε > εF, eigenvalues[1])

    hole_space = select_orbitals(active_space, idx_holes)
    particle_space = select_orbitals(active_space, idx_parts)

    # --- dump DeltaIntegralsHH
    DeltaIntegralsHH = compute_delta_integrals(basis, hole_space.ψ, Val(:HH))
    Nocc = size(DeltaIntegralsHH, 1)

    # Row-major ordering for C++: i, j
    tensor_data_hh = (
        convert(Complex{Cdouble}, DeltaIntegralsHH[i, j])
        for i in 1:Nocc for j in 1:Nocc
    )
    dim_hh = [
        Dict("length" => Nocc, "type" => "Hole"),
        Dict("length" => Nocc, "type" => "Hole"),
    ]
    files_hh = write_cc4s_tensor(
        folder,
        "DeltaIntegralsHH",
        tensor_data_hh;
        dimensions=dim_hh,
        scalarType="Complex64",
        elementType="IeeeBinaryFile",
    )

    # --- dump DeltaIntegralsPPHH
    DeltaIntegralsPPHH = compute_delta_integrals(
        basis,
        particle_space.ψ,
        hole_space.ψ,
        Val(:PPHH),
    )
    Nvirt = size(DeltaIntegralsPPHH, 1)

    # Row-major ordering for C++: a, b, i, j
    tensor_data_pphh = (
        convert(Complex{Cdouble}, DeltaIntegralsPPHH[a, b, i, j])
        for a in 1:Nvirt for b in 1:Nvirt for i in 1:Nocc for j in 1:Nocc
    )
    dim_pphh = [
        Dict("length" => Nvirt, "type" => "Particle"),
        Dict("length" => Nvirt, "type" => "Particle"),
        Dict("length" => Nocc, "type" => "Hole"),
        Dict("length" => Nocc, "type" => "Hole"),
    ]
    files_pphh = write_cc4s_tensor(
        folder,
        "DeltaIntegralsPPHH",
        tensor_data_pphh;
        dimensions=dim_pphh,
        scalarType="Complex64",
        elementType="IeeeBinaryFile",
    )

    # --- dump Grid Vectors
    files_grid = write_grid_vectors(folder, basis, fitting.G_vectors)

    # --- dump Coulomb Potential
    files_pot = write_coulomb_potential(folder, fitting.kernel_fourier)

    # --- dump Coulomb Vertex Singular Vectors
    U = fitting.singular_vectors
    files_u = isnothing(U) ? String[] : write_singular_vectors(folder, U)

    return vcat(files_ene, files_coul, files_hh, files_pphh, files_grid, files_pot, files_u)
end

# All requirements of the Cc4s interface, which the writers below only assert
function _check_cc4s_input(active_space, fitting, folder, force)
    kpoints = active_space.basis.kpoints
    length(kpoints) == 1 && iszero(only(kpoints).coordinate) ||
        throw(ArgumentError("Cc4s interface supports only a single k-point at Γ. Run " *
                            "the SCF with `kgrid=[1, 1, 1]` and without spin " *
                            "polarization."))
    # Cc4s expects canonical HF orbitals: orthonormal, with Fock eigenvalues
    active_space.is_orthonormal ||
        throw(ArgumentError("Cc4s interface requires orthonormal canonical orbitals. " *
                            "Apply `canonicalize_orbitals` to the active space first."))
    ε = only(active_space.eigenvalues)
    maximum(abs, sort(ε) - ε; init=0.0) < 1e-10 ||
        throw(ArgumentError("Cc4s interface requires eigenvalues in non-decreasing " *
                            "order. Apply `canonicalize_orbitals` to the active space " *
                            "first."))
    # Cc4s only supports gapped systems with integer occupancies
    for occ in only(active_space.occupation)
        occ ≈ 0.0 || occ ≈ 1.0 || occ ≈ 2.0 ||
            throw(ArgumentError("Cc4s interface requires integer occupations, found " *
                                "$occ. Use a gapped reference without smearing."))
    end
    N = size(only(active_space.ψ), 2)
    size(fitting.Γ)[1:4] == (1, N, 1, N) ||
        throw(DimensionMismatch("The Coulomb vertex of size $(size(fitting.Γ)[1:4]) " *
                                "does not match the active space of $N orbitals. " *
                                "Compute `fitting` from `active_space` itself."))
    force && return nothing
    for name in _cc4s_tensor_names(fitting), file in _cc4s_files(folder, name)
        isfile(file) &&
            throw(ArgumentError("$file exists. Pass `force=true` to overwrite it or " *
                                "choose another `folder`."))
    end
    return nothing
end

function _cc4s_tensor_names(fitting)
    names = ["EigenEnergies", "CoulombVertex", "DeltaIntegralsHH", "DeltaIntegralsPPHH",
             "GridVectors", "CoulombPotential"]
    isnothing(fitting.singular_vectors) || push!(names, "CoulombVertexSingularVectors")
    return names
end

function _cc4s_files(folder, name)
    return joinpath(folder, "$name.yaml"), joinpath(folder, "$name.elements")
end

# Write CoulombVertex.yaml and CoulombVertex.elements
function write_coulomb_vertex(
    folder::AbstractString,
    ΓnmF::AbstractArray{T,5},
) where {T}
    Nk = size(ΓnmF, 1)
    N = size(ΓnmF, 2)
    NF = size(ΓnmF, 5)
    @assert Nk == size(ΓnmF, 3)
    @assert N == size(ΓnmF, 4)
    @assert Nk == 1  # 1 kpt is hard-coded for now (see write_eigenenergies)

    dimensions = [
        Dict("length" => NF, "type" => "AuxiliaryField"),
        Dict("length" => Nk * N, "type" => "State"),
        Dict("length" => Nk * N, "type" => "State"),
    ]
    metaData = Dict("halfGrid" => 0)  # Complex integrals

    # Cc4s is written in C++
    # C++ is row-major, julia is column-major. Therefore we write
    # ΓnmF in a stream using chunks of all field Fs for given (n,m)
    tensor_data = (
        convert(Vector{Complex{Cdouble}}, vec(ΓnmF[1, n, 1, m, :]))
        for n in 1:N for m in 1:N
    )

    return write_cc4s_tensor(
        folder,
        "CoulombVertex",
        tensor_data;
        dimensions,
        scalarType="Complex64",
        elementType="IeeeBinaryFile",
        metaData,
    )
end

# Write GridVectors.yaml and GridVectors.elements
function write_grid_vectors(
    folder::AbstractString,
    basis::PlaneWaveBasis,
    G_vectors::AbstractVector,
)
    # The GridVectors object contains the grid vectors of the employed plane-wave basis set
    model = basis.model
    # Convert integer G-vectors to Cartesian coordinates
    G_cartesian = [model.recip_lattice * G for G in G_vectors]

    dimensions = [
        Dict("length" => 3, "type" => "Vector"),
        Dict("length" => length(G_cartesian), "type" => "Momentum"),
    ]
    # Unit is 1.0 (Bohr^-1)

    # Gi, Gj, Gk for reference only, they are not needed for Cartesian coordinates
    metaData = Dict(
        "Gi" => model.recip_lattice[:, 1],
        "Gj" => model.recip_lattice[:, 2],
        "Gk" => model.recip_lattice[:, 3],
    )

    tensor_data = (G[i] for G in G_cartesian for i in 1:3)

    return write_cc4s_tensor(
        folder,
        "GridVectors",
        tensor_data;
        dimensions,
        scalarType="Real64",
        elementType="TextFile",
        metaData,
        unit=1.0,
    )
end

# Write CoulombPotential.yaml and CoulombPotential.elements
function write_coulomb_potential(
    folder::AbstractString,
    kernel_fourier::AbstractVector,
)
    dimensions = [Dict("length" => length(kernel_fourier), "type" => "Momentum")]

    return write_cc4s_tensor(
        folder,
        "CoulombPotential",
        kernel_fourier;
        dimensions,
        scalarType="Real64",
        elementType="TextFile",
        unit=1.0,
    )
end

# Write CoulombVertexSingularVectors.yaml and CoulombVertexSingularVectors.elements
function write_singular_vectors(
    folder::AbstractString,
    coulomb_vertex_singular_vectors::AbstractMatrix{T},
) where {T}
    # coulomb_vertex_singular_vectors has dimensions (NG, NF)
    NG, NF = size(coulomb_vertex_singular_vectors)
    dimensions = [
        Dict("length" => NF, "type" => "AuxiliaryField"),
        Dict("length" => NG, "type" => "Momentum"),
    ]

    # row-major write: loop over F then G
    tensor_data = (
        convert(Complex{Cdouble}, coulomb_vertex_singular_vectors[iG, iF])
        for iF in 1:NF for iG in 1:NG
    )

    return write_cc4s_tensor(
        folder,
        "CoulombVertexSingularVectors",
        tensor_data;
        dimensions,
        scalarType="Complex64",
        elementType="IeeeBinaryFile",
        unit=1.0,
    )
end
