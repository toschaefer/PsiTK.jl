# This file initially based on code of the experimental "cc4s" branch in DFTK written by Michael Herbst

using YAML

# Write EigenEnergies.yaml and EigenEnergies.elements
function write_eigenenergies(
    folder::AbstractString,
    eigenvalues::AbstractVector,
    εF::Number;
    force=false
)
    @assert length(eigenvalues) == 1
    εk = eigenvalues[1]

    # Note: Eigenenergies need to be ordered in *non-decreasing* order for cc4s !
    @assert maximum(abs, sort(εk) - εk) < 1e-10

    yamlfile     = joinpath(folder, "EigenEnergies.yaml")
    elementsfile = joinpath(folder, "EigenEnergies.elements")
    if !force && (isfile(yamlfile) || isfile(elementsfile))
        error("Generated files $yamlfile and/or $elementsfile exists.")
    end

    metadata = Dict(
        "fermiEnergy" => εF,
        "energies"    => εk,
    )
    data = Dict(
        "version"    => 100,
        "type"       => "Tensor",
        "scalarType" => "Real64",
        "dimensions" => [Dict("length" => length(εk),
                              "type"   => "State")],
        "elements"   => Dict("type" => "TextFile"),
        "unit"       => 1.0,  # DFTK using Hartree as well
        "metaData"   => metadata,
    )
    open(fp -> YAML.write(fp, data), yamlfile, "w")

    open(elementsfile, "w") do fp
        for ε in εk
            println(fp, ε)
        end
    end

    [yamlfile, elementsfile]
end


# Write CoulombVertex.yaml and CoulombVertex.elements
function write_coulomb_vertex(
    folder::AbstractString, 
    ΓnmF::AbstractArray{T, 5};
    force=true
) where {T}
    n_kpt   = size(ΓnmF, 1)
    n_bands = size(ΓnmF, 2)
    n_aux_field = size(ΓnmF, 5)
    @assert n_kpt   == size(ΓnmF, 3)
    @assert n_bands == size(ΓnmF, 4)
    @assert n_kpt  == 1  # 1 kpt is hard-coded for now (see write_eigenenergies)

    yamlfile     = joinpath(folder, "CoulombVertex.yaml")
    elementsfile = joinpath(folder, "CoulombVertex.elements")
    if !force && (isfile(yamlfile) || isfile(elementsfile))
        error("Generated files $yamlfile and/or $elementsfile exists.")
    end

    dimensions = [
        Dict("length" => n_aux_field,     "type" => "AuxiliaryField"),
        Dict("length" => n_kpt * n_bands, "type" => "State"),
        Dict("length" => n_kpt * n_bands, "type" => "State"),
    ]
    data = Dict(
        "version"    => 100,
        "type"       => "Tensor",
        "scalarType" => "Complex64",
        "dimensions" => dimensions,
        "elements"   => Dict("type" => "IeeeBinaryFile"),
        "unit"       => 1.0,  # DFTK using Hartree as well
        "metaData"   => Dict("halfGrid" => 0),  # Complex integrals
    )
    open(fp -> YAML.write(fp, data), yamlfile, "w")

    # Cc4s is writte in C++
    # C++ is row-major, julia is column-major. Therefore we write
    # ΓnmF in a stream using chunks of all field Fs for given (n,m)
    open(elementsfile, "w") do fp
        for n = 1:n_bands, m = 1:n_bands
           binary = convert(Vector{Complex{Cdouble}}, vec(ΓnmF[1,n,1,m,:]))
           write(fp, binary)
        end # n,m
    end

    [yamlfile, elementsfile]
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

# Arguments
- `active_space`: the `OrbitalSpace` containing the bands (occupations and eigenvalues)
- `ΓmnG`: the (compressed) Coulomb vertex
- `folder`: the target folder
- `force`: if true existing files will be overwritten
"""
function dump_cc4s_files(
    active_space::OrbitalSpace,
    ΓmnG::AbstractArray,
    folder::AbstractString=joinpath(pwd(), "cc4s");
    force=false
)
    mkpath(folder)

    # --- dump Eigenvalues
    # For cc4s we just pass the eigenvalues from the active space
    # (Assuming 1 kpoint for now)
    eigenvalues = active_space.eigenvalues
    εF = active_space.εF  # Fermi level
    
    files_ene = write_eigenenergies(folder, eigenvalues, εF; force)

    # --- dump Coulomb Vertex
    files_coul = write_coulomb_vertex(folder, ΓmnG; force)

    append!(files_ene, files_coul)
end
