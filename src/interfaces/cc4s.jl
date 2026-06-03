# This file initially based on code of the experimental "cc4s" branch in DFTK written by Michael Herbst

using YAML
using DFTK

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
    scalarType::String = "Real64",
    elementType::String = "IeeeBinaryFile",
    metaData::Dict = Dict{String,Any}(),
    unit::Float64 = 1.0,
    force = false
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
    force = false,
)
    @assert length(eigenvalues) == 1
    εk = eigenvalues[1]

    # Note: Eigenenergies need to be ordered in *non-decreasing* order for cc4s !
    @assert maximum(abs, sort(εk) - εk) < 1e-10

    dimensions = [Dict("length" => length(εk), "type" => "State")]
    metaData = Dict("fermiEnergy" => εF, "energies" => εk)

    return write_cc4s_tensor(
        folder, "EigenEnergies", εk;
        dimensions = dimensions,
        scalarType = "Real64",
        elementType = "TextFile",
        metaData = metaData,
        force = force
    )
end

"""
    dump_cc4s_files(
        active_space::OrbitalSpace,
        ΓmnG::AbstractArray,
        G_vectors::AbstractVector = SVector{3, Int}[],
        kernel_fourier::AbstractVector = Float64[];
        coulomb_vertex_singular_vectors::Union{AbstractMatrix, Nothing} = nothing,
        folder::AbstractString=joinpath(pwd(), "cc4s"),
        force=false
    )

Write Cc4s input files (*.yaml and *.elements):
- EigenEnergies
- CoulombVertex
- DeltaIntegralsHH
- DeltaIntegralsPPHH
- GridVectors (if G_vectors is provided)
- CoulombPotential (if kernel_fourier is provided)
- CoulombVertexSingularVectors (if coulomb_vertex_singular_vectors is provided)

# Arguments
- `active_space`: the `OrbitalSpace` containing the bands (occupations and eigenvalues)
- `ΓmnG`: the (compressed) Coulomb vertex
- `G_vectors`: corresponding plane-wave vectors
- `kernel_fourier`: the evaluated Coulomb potential at the given G_vectors
- `coulomb_vertex_singular_vectors`: transformation matrix from Coulomb vertex compression
- `folder`: the target folder
- `force`: if true existing files will be overwritten
"""
function dump_cc4s_files(
    active_space::OrbitalSpace,
    ΓmnG::AbstractArray,
    G_vectors::AbstractVector,
    kernel_fourier::AbstractVector;
    coulomb_vertex_singular_vectors::Union{AbstractMatrix, Nothing} = nothing,
    folder::AbstractString = joinpath(pwd(), "cc4s"),
    force = false,
)
    mkpath(folder)

    # --- dump Eigenvalues
    # For cc4s we just pass the eigenvalues from the active space
    # (Assuming 1 kpoint for now)
    eigenvalues = active_space.eigenvalues
    εF = active_space.εF  # Fermi level

    # Cc4s only supports gapped systems with integer occupancies
    for occ in active_space.occupations[1]
        if !(occ ≈ 0.0 || occ ≈ 1.0 || occ ≈ 2.0)
            error("Cc4s interface requires integer occupancies. Fractional occupations detected.")
        end
    end

    files_ene = write_eigenenergies(folder, eigenvalues, εF; force)

    # --- dump Coulomb Vertex
    files_coul = write_coulomb_vertex(folder, ΓmnG; force)

    # --- Split space by Fermi Energy for Delta Integrals
    idx_holes = findall(ε -> ε <= εF, eigenvalues[1])
    idx_parts = findall(ε -> ε > εF, eigenvalues[1])

    hole_space = select_orbitals(active_space, idx_holes)
    particle_space = select_orbitals(active_space, idx_parts)

    # --- dump DeltaIntegralsHH
    DeltaIntegralsHH = compute_delta_integrals(hole_space, Val(:HH))
    N_occ = size(DeltaIntegralsHH, 1)
    
    # Row-major ordering for C++: i, j
    tensor_data_hh = (convert(Complex{Cdouble}, DeltaIntegralsHH[i, j]) for i=1:N_occ for j=1:N_occ)
    dim_hh = [
        Dict("length" => N_occ, "type" => "Hole"),
        Dict("length" => N_occ, "type" => "Hole")
    ]
    files_hh = write_cc4s_tensor(
        folder, "DeltaIntegralsHH", tensor_data_hh;
        dimensions = dim_hh,
        scalarType = "Complex64",
        elementType = "IeeeBinaryFile",
        force = force
    )

    # --- dump DeltaIntegralsPPHH
    DeltaIntegralsPPHH = compute_delta_integrals(particle_space, hole_space, Val(:PPHH))
    N_virt = size(DeltaIntegralsPPHH, 1)
    
    # Row-major ordering for C++: a, b, i, j
    tensor_data_pphh = (
        convert(Complex{Cdouble}, DeltaIntegralsPPHH[a, b, i, j]) 
        for a=1:N_virt for b=1:N_virt for i=1:N_occ for j=1:N_occ
    )
    dim_pphh = [
        Dict("length" => N_virt, "type" => "Particle"),
        Dict("length" => N_virt, "type" => "Particle"),
        Dict("length" => N_occ, "type" => "Hole"),
        Dict("length" => N_occ, "type" => "Hole")
    ]
    files_pphh = write_cc4s_tensor(
        folder, "DeltaIntegralsPPHH", tensor_data_pphh;
        dimensions = dim_pphh,
        scalarType = "Complex64",
        elementType = "IeeeBinaryFile",
        force = force
    )

    # --- dump Grid Vectors
    files_grid = write_grid_vectors(folder, active_space.basis, G_vectors; force)

    # --- dump Coulomb Potential
    files_pot = write_coulomb_potential(folder, kernel_fourier; force)

    # --- dump Coulomb Vertex Singular Vectors
    files_u = isnothing(coulomb_vertex_singular_vectors) ? String[] : write_singular_vectors(folder, coulomb_vertex_singular_vectors; force)

    return vcat(files_ene, files_coul, files_hh, files_pphh, files_grid, files_pot, files_u)
end


# Write CoulombVertex.yaml and CoulombVertex.elements
function write_coulomb_vertex(
    folder::AbstractString,
    ΓnmF::AbstractArray{T,5};
    force = true,
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

    # Cc4s is writte in C++
    # C++ is row-major, julia is column-major. Therefore we write
    # ΓnmF in a stream using chunks of all field Fs for given (n,m)
    tensor_data = (
        convert(Vector{Complex{Cdouble}}, vec(ΓnmF[1, n, 1, m, :])) 
        for n = 1:n_bands for m = 1:n_bands
    )

    return write_cc4s_tensor(
        folder, "CoulombVertex", tensor_data;
        dimensions = dimensions,
        scalarType = "Complex64",
        elementType = "IeeeBinaryFile",
        metaData = metaData,
        force = force
    )
end

# Write GridVectors.yaml and GridVectors.elements
function write_grid_vectors(
    folder::AbstractString,
    basis::PlaneWaveBasis,
    G_vectors::AbstractVector;
    force = false,
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
    
    # We still provide the Gi, Gj, Gk for reference, though not strictly needed for Cartesian
    metaData = Dict(
        "Gi" => model.recip_lattice[:, 1],
        "Gj" => model.recip_lattice[:, 2],
        "Gk" => model.recip_lattice[:, 3],
    )

    tensor_data = (
        G[i] for G in G_cartesian for i = 1:3
    )

    return write_cc4s_tensor(
        folder, "GridVectors", tensor_data;
        dimensions = dimensions,
        scalarType = "Real64",
        elementType = "TextFile",
        metaData = metaData,
        unit = 1.0,
        force = force
    )
end

# Write CoulombPotential.yaml and CoulombPotential.elements
function write_coulomb_potential(
    folder::AbstractString,
    kernel_fourier::AbstractVector;
    force = false,
)
    dimensions = [
        Dict("length" => length(kernel_fourier), "type" => "Momentum")
    ]

    return write_cc4s_tensor(
        folder, "CoulombPotential", kernel_fourier;
        dimensions = dimensions,
        scalarType = "Real64",
        elementType = "TextFile",
        unit = 1.0,
        force = force
    )
end

# Write CoulombVertexSingularVectors.yaml and CoulombVertexSingularVectors.elements
function write_singular_vectors(
    folder::AbstractString,
    coulomb_vertex_singular_vectors::AbstractMatrix{T};
    force = false,
) where {T}
    # coulomb_vertex_singular_vectors has dimensions (N_G, N_F)
    N_G, N_F = size(coulomb_vertex_singular_vectors)
    dimensions = [
        Dict("length" => N_F, "type" => "AuxiliaryField"),
        Dict("length" => N_G, "type" => "Momentum")
    ]
    
    # row-major write: loop over F then G
    tensor_data = (
        convert(Complex{Cdouble}, coulomb_vertex_singular_vectors[iG, iF]) for iF = 1:N_F for iG = 1:N_G
    )

    return write_cc4s_tensor(
        folder, "CoulombVertexSingularVectors", tensor_data;
        dimensions = dimensions,
        scalarType = "Complex64",
        elementType = "IeeeBinaryFile",
        unit = 1.0,
        force = force
    )
end
