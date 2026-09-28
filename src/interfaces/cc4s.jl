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
    force=false,
)
    yamlfile = joinpath(folder, "$name.yaml")
    elementsfile = joinpath(folder, "$name.elements")
    if !force && (isfile(yamlfile) || isfile(elementsfile))
        error("Generated files $yamlfile and/or $elementsfile exists.")
    end

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
        elseif elementType == "IeeeBinaryFile"
            for val in tensor_data
                write(fp, val)
            end
        else
            error("Unsupported elementType $elementType")
        end
    end

    return [yamlfile, elementsfile]
end

# Write EigenEnergies.yaml and EigenEnergies.elements
function write_eigenenergies(
    folder::AbstractString,
    eigenvalues::AbstractVector,
    εF::Number;
    force=false,
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
        force,
    )
end

"""
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

Requires a Gamma-only calculation with integer occupations.

# Arguments
- `active_space`: the `OrbitalSpace` containing the bands (occupation and eigenvalues)
- `fitting`: the [`DensityFitting`](@ref) of the active space, uncompressed from
  [`compute_coulomb_vertex`](@ref) or compressed by [`compress_coulomb_vertex`](@ref)
- `folder`: the target folder (created if missing)
- `force`: if true existing files will be overwritten

# Returns
The list of written file paths.
"""
function dump_cc4s_files(
    active_space::OrbitalSpace,
    fitting::DensityFitting;
    folder::AbstractString=joinpath(pwd(), "cc4s"),
    force=false,
)
    # Cc4s expects canonical HF orbitals: orthonormal, with Fock eigenvalues
    if !active_space.is_orthonormal
        error("Cc4s interface requires orthonormal canonical orbitals. " *
              "Apply `canonicalize_orbitals` to the active space first.")
    end
    size(fitting.Γ, 2) == size(fitting.Γ, 4) == size(active_space.ψ[1], 2) ||
        error("The Coulomb vertex does not match the active space. Compute `fitting` " *
              "from `active_space` itself.")
    mkpath(folder)

    # --- dump Eigenvalues
    # For cc4s we just pass the eigenvalues from the active space
    # (Assuming 1 kpoint for now)
    basis = active_space.basis
    eigenvalues = active_space.eigenvalues
    εF = active_space.εF  # Fermi level

    # Cc4s only supports gapped systems with integer occupancies
    for occ in active_space.occupation[1]
        if !(occ ≈ 0.0 || occ ≈ 1.0 || occ ≈ 2.0)
            error("Cc4s interface requires integer occupancies. " *
                  "Fractional occupations detected.")
        end
    end

    files_ene = write_eigenenergies(folder, eigenvalues, εF; force)

    # --- dump Coulomb Vertex
    files_coul = write_coulomb_vertex(folder, fitting.Γ; force)

    # --- Split space by Fermi Energy for Delta Integrals
    idx_holes = findall(ε -> ε <= εF, eigenvalues[1])
    idx_parts = findall(ε -> ε > εF, eigenvalues[1])

    hole_space = select_orbitals(active_space, idx_holes)
    particle_space = select_orbitals(active_space, idx_parts)

    # --- dump DeltaIntegralsHH
    DeltaIntegralsHH = compute_delta_integrals(basis, hole_space.ψ, Val(:HH))
    N_occ = size(DeltaIntegralsHH, 1)

    # Row-major ordering for C++: i, j
    tensor_data_hh = (
        convert(Complex{Cdouble}, DeltaIntegralsHH[i, j])
        for i in 1:N_occ for j in 1:N_occ
    )
    dim_hh = [
        Dict("length" => N_occ, "type" => "Hole"),
        Dict("length" => N_occ, "type" => "Hole"),
    ]
    files_hh = write_cc4s_tensor(
        folder,
        "DeltaIntegralsHH",
        tensor_data_hh;
        dimensions=dim_hh,
        scalarType="Complex64",
        elementType="IeeeBinaryFile",
        force,
    )

    # --- dump DeltaIntegralsPPHH
    DeltaIntegralsPPHH = compute_delta_integrals(
        basis,
        particle_space.ψ,
        hole_space.ψ,
        Val(:PPHH),
    )
    N_virt = size(DeltaIntegralsPPHH, 1)

    # Row-major ordering for C++: a, b, i, j
    tensor_data_pphh = (
        convert(Complex{Cdouble}, DeltaIntegralsPPHH[a, b, i, j])
        for a in 1:N_virt for b in 1:N_virt for i in 1:N_occ for j in 1:N_occ
    )
    dim_pphh = [
        Dict("length" => N_virt, "type" => "Particle"),
        Dict("length" => N_virt, "type" => "Particle"),
        Dict("length" => N_occ, "type" => "Hole"),
        Dict("length" => N_occ, "type" => "Hole"),
    ]
    files_pphh = write_cc4s_tensor(
        folder,
        "DeltaIntegralsPPHH",
        tensor_data_pphh;
        dimensions=dim_pphh,
        scalarType="Complex64",
        elementType="IeeeBinaryFile",
        force,
    )

    # --- dump Grid Vectors
    files_grid = write_grid_vectors(folder, basis, fitting.G_vectors; force)

    # --- dump Coulomb Potential
    files_pot = write_coulomb_potential(folder, fitting.kernel_fourier; force)

    # --- dump Coulomb Vertex Singular Vectors
    U = fitting.singular_vectors
    files_u = isnothing(U) ? String[] : write_singular_vectors(folder, U; force)

    return vcat(files_ene, files_coul, files_hh, files_pphh, files_grid, files_pot, files_u)
end

# Write CoulombVertex.yaml and CoulombVertex.elements
function write_coulomb_vertex(
    folder::AbstractString,
    ΓnmF::AbstractArray{T,5};
    force=true,
) where {T}
    n_kpt = size(ΓnmF, 1)
    n_bands = size(ΓnmF, 2)
    n_aux_field = size(ΓnmF, 5)
    @assert n_kpt == size(ΓnmF, 3)
    @assert n_bands == size(ΓnmF, 4)
    @assert n_kpt == 1  # 1 kpt is hard-coded for now (see write_eigenenergies)

    dimensions = [
        Dict("length" => n_aux_field, "type" => "AuxiliaryField"),
        Dict("length" => n_kpt * n_bands, "type" => "State"),
        Dict("length" => n_kpt * n_bands, "type" => "State"),
    ]
    metaData = Dict("halfGrid" => 0)  # Complex integrals

    # Cc4s is written in C++
    # C++ is row-major, julia is column-major. Therefore we write
    # ΓnmF in a stream using chunks of all field Fs for given (n,m)
    tensor_data = (
        convert(Vector{Complex{Cdouble}}, vec(ΓnmF[1, n, 1, m, :]))
        for n in 1:n_bands for m in 1:n_bands
    )

    return write_cc4s_tensor(
        folder,
        "CoulombVertex",
        tensor_data;
        dimensions,
        scalarType="Complex64",
        elementType="IeeeBinaryFile",
        metaData,
        force,
    )
end

# Write GridVectors.yaml and GridVectors.elements
function write_grid_vectors(
    folder::AbstractString,
    basis::PlaneWaveBasis,
    G_vectors::AbstractVector;
    force=false,
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
        force,
    )
end

# Write CoulombPotential.yaml and CoulombPotential.elements
function write_coulomb_potential(
    folder::AbstractString,
    kernel_fourier::AbstractVector;
    force=false,
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
        force,
    )
end

# Write CoulombVertexSingularVectors.yaml and CoulombVertexSingularVectors.elements
function write_singular_vectors(
    folder::AbstractString,
    coulomb_vertex_singular_vectors::AbstractMatrix{T};
    force=false,
) where {T}
    # coulomb_vertex_singular_vectors has dimensions (N_G, N_F)
    N_G, N_F = size(coulomb_vertex_singular_vectors)
    dimensions = [
        Dict("length" => N_F, "type" => "AuxiliaryField"),
        Dict("length" => N_G, "type" => "Momentum"),
    ]

    # row-major write: loop over F then G
    tensor_data = (
        convert(Complex{Cdouble}, coulomb_vertex_singular_vectors[iG, iF])
        for iF in 1:N_F for iG in 1:N_G
    )

    return write_cc4s_tensor(
        folder,
        "CoulombVertexSingularVectors",
        tensor_data;
        dimensions,
        scalarType="Complex64",
        elementType="IeeeBinaryFile",
        unit=1.0,
        force,
    )
end
