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
        "unit" => 1.0,  # DFTK using Hartree as well
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

"""
    dump_cc4s_files(
        active_space::OrbitalSpace,
        ΓmnG::AbstractArray,
        folder::AbstractString=joinpath(pwd(), "cc4s");
        force=false
    )

Write Cc4s input files (*.yaml and *.elements):
- EigenEnergies
- CoulombVertex
- DeltaIntegralsHH
- DeltaIntegralsPPHH

# Arguments
- `active_space`: the `OrbitalSpace` containing the bands (occupations and eigenvalues)
- `ΓmnG`: the (compressed) Coulomb vertex
- `folder`: the target folder
- `force`: if true existing files will be overwritten
"""
function dump_cc4s_files(
    active_space::OrbitalSpace,
    ΓmnG::AbstractArray,
    folder::AbstractString = joinpath(pwd(), "cc4s");
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

    hole_space = select_orbitals(active_space, IndexSelection(idx_holes))
    particle_space = select_orbitals(active_space, IndexSelection(idx_parts))

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

    return vcat(files_ene, files_coul, files_hh, files_pphh)
end
