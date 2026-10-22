# This file initially based on code of the experimental "cc4s" branch in DFTK
# written by Michael Herbst

@doc raw"""
    write_cc4s_tensor(folder, name, tensor_data; kwargs...)

Generic function to write a Tensor object in Cc4s format, including its YAML
metadata file and its binary/text elements file. `tensor_data` iterates the elements with
the first index of `dimensions` fastest, as Julia stores arrays: Cc4s fills its tensors by
the global index ``g = i_1 + n_1 (i_2 + n_2 (i_3 + \dots))``.
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
                _print_cc4s_element(fp, val)
            end
        else
            for val in tensor_data
                write(fp, val)
            end
        end
    end

    return [yamlfile, elementsfile]
end

# Cc4s text files hold one element per line, complex numbers as (re,im)
_print_cc4s_element(io, x::Real) = println(io, x)
_print_cc4s_element(io, z::Complex) = println(io, '(', real(z), ',', imag(z), ')')

@doc raw"""
    read_cc4s_tensor(folder, name)

Read the tensor `name` in Cc4s format from the files `name.yaml` and `name.elements` in
`folder`, as written by Cc4s or by [`dump_cc4s_files`](@ref). For a set of tensors, e.g.
the blocks of a density matrix, use [`read_cc4s_tensor_set`](@ref).

The indices of the returned array follow the order of the `dimensions` in the YAML file.
Cc4s stores the first index fastest, as Julia does, unless the elements are marked as
`rowMajor`, in which case the array is permuted accordingly.

Throws an `ArgumentError` if a file is missing, describes a tensor set or uses a format
that is not supported (a version other than 100, block-sparse tensors, the `halfGrid`
format of real Coulomb vertices), and a `DimensionMismatch` if the number of elements in
the file does not match the dimensions.

# Returns
A `NamedTuple` `(; tensor, dimensions, unit, metaData)` with
- `tensor`: an `Array{Float64}` for the scalar type `Real64` or an `Array{ComplexF64}`
  for `Complex64`, not scaled by `unit`;
- `dimensions`: a `Vector` of `Dict`s, one per index, with its `"length"` and optionally
  its `"type"`, e.g. `"Hole"`, `"Particle"` or `"State"`;
- `unit`: the unit of the elements;
- `metaData`: the `Dict` of meta data, empty if the file has none.

# Cost
Time and memory are linear in the number of elements, e.g. ``N^2 N_F`` complex numbers
for the Coulomb vertex.
"""
function read_cc4s_tensor(folder::AbstractString, name::AbstractString)
    node = _load_cc4s_yaml(folder, name)
    _is_cc4s_tensor_set(node) &&
        throw(ArgumentError("$name.yaml describes a Cc4s tensor set. Read it with " *
                            "`read_cc4s_tensor_set`."))
    return _read_cc4s_tensor_node(node, joinpath(folder, name), "$name.yaml")
end

@doc raw"""
    read_cc4s_tensor_set(folder, name)

Read the set of tensors `name` in Cc4s format from the file `name.yaml` in `folder`, which
lists the tensors under `components`, and the files `name.components.<key>.elements`
holding the elements of the component `<key>`. Cc4s writes e.g. the one-particle reduced
density matrix as the set `rdm` of its blocks `"hh"`, `"hp"`, `"ph"` and `"pp"`, which
[`merge_hole_particle_blocks`](@ref) assembles into the full matrix.

Throws the exceptions of [`read_cc4s_tensor`](@ref), and an `ArgumentError` if
`name.yaml` does not describe a tensor set.

# Returns
A `Dict` from the key of each component to its `NamedTuple`
`(; tensor, dimensions, unit, metaData)`, see [`read_cc4s_tensor`](@ref).

# Cost
Time and memory are linear in the number of elements of all components.
"""
function read_cc4s_tensor_set(folder::AbstractString, name::AbstractString)
    node = _load_cc4s_yaml(folder, name)
    _is_cc4s_tensor_set(node) ||
        throw(ArgumentError("$name.yaml does not describe a Cc4s tensor set. Read it " *
                            "with `read_cc4s_tensor`."))
    path = joinpath(folder, "$name.components")
    return Dict{String,Any}(
        string(key) => _read_cc4s_tensor_node(
            component, "$path.$key", "$name.yaml (component $key)",
        )
        for (key, component) in node["components"]
    )
end

@doc raw"""
    merge_hole_particle_blocks(components)

Assemble a tensor over all states from its blocks over holes and particles, e.g. the
one-particle reduced density matrix ``D_{pq}`` from its blocks `"hh"`, `"hp"`, `"ph"` and
`"pp"`, as returned by [`read_cc4s_tensor_set`](@ref). The key of a block has one letter
per index, `h` for the ``N_\text{occ}`` holes and `p` for the ``N_\text{virt}``
particles. In each index of the result, the holes come first, followed by the particles,
which is the order of the states written by [`dump_cc4s_files`](@ref).

Throws an `ArgumentError` unless the keys are exactly all ``2^r`` combinations of `h` and
`p` for a tensor of rank ``r`` (e.g. the Coulomb integrals of Cc4s lack blocks), and a
`DimensionMismatch` if the sizes of the blocks are inconsistent.

# Returns
The `Array` of size ``(N_\text{occ} + N_\text{virt})^r``.

# Cost
Time and memory are linear in the number of elements, ``(N_\text{occ} + N_\text{virt})^r``.
"""
function merge_hole_particle_blocks(components::AbstractDict)
    rank = length(first(keys(components)))
    blocks = Set(join(letters) for letters in Iterators.product(ntuple(_ -> "hp", rank)...))
    Set(keys(components)) == blocks ||
        throw(ArgumentError("The blocks $(sort!(collect(keys(components)))) are not all " *
                            "combinations of h and p for a tensor of rank $rank."))
    Nocc = size(components["h"^rank].tensor, 1)
    Nvirt = size(components["p"^rank].tensor, 1)
    T = mapreduce(c -> eltype(c.tensor), promote_type, values(components))
    tensor = zeros(T, ntuple(_ -> Nocc + Nvirt, rank))
    for (key, component) in components
        ranges = Tuple(letter == 'h' ? (1:Nocc) : (Nocc + 1:Nocc + Nvirt) for letter in key)
        size(component.tensor) == length.(ranges) ||
            throw(DimensionMismatch("The block $key has the size " *
                                    "$(size(component.tensor)), expected " *
                                    "$(length.(ranges)) for $Nocc holes and $Nvirt " *
                                    "particles."))
        tensor[ranges...] = component.tensor
    end
    return tensor
end

function _load_cc4s_yaml(folder, name)
    yamlfile = joinpath(folder, "$name.yaml")
    isfile(yamlfile) || throw(ArgumentError("Cc4s file $yamlfile not found."))
    return YAML.load_file(yamlfile)
end

# Cc4s writes a tensor set with `type: TensorSet`. Without a type, Cc4s reads the map
# entry by entry, which reaches the same element files.
function _is_cc4s_tensor_set(node)
    return haskey(node, "components") && get(node, "type", "TensorSet") == "TensorSet"
end

# Read the tensor described by `node` from the file `path.elements`; `source` names the
# YAML node in error messages
function _read_cc4s_tensor_node(node, path, source)
    elementsfile = "$path.elements"
    isfile(elementsfile) || throw(ArgumentError("Cc4s file $elementsfile not found."))
    _check_cc4s_tensor_node(node, source)

    T = node["scalarType"] == "Real64" ? Float64 : ComplexF64
    dimensions = _cc4s_list(node["dimensions"])
    lens = Tuple(Int(dim["length"]) for dim in dimensions)
    elements = if node["elements"]["type"] == "TextFile"
        _read_cc4s_elements_text(elementsfile, T, prod(lens))
    else
        _read_cc4s_elements_binary(elementsfile, T, prod(lens))
    end
    tensor = if get(node["elements"], "rowMajor", false)
        permutedims(reshape(elements, reverse(lens)), length(lens):-1:1)
    else
        reshape(elements, lens)
    end
    metaData = something(get(node, "metaData", nothing), Dict{String,Any}())
    return (; tensor, dimensions, unit=Float64(node["unit"]), metaData)
end

# All format requirements of a tensor, which the element readers only assert
function _check_cc4s_tensor_node(node, source)
    get(node, "type", nothing) == "Tensor" ||
        throw(ArgumentError("$source does not describe a Cc4s tensor."))
    get(node, "version", nothing) == 100 ||
        throw(ArgumentError("$source has the Cc4s format version " *
                            "$(get(node, "version", nothing)), only 100 is supported."))
    node["scalarType"] in ("Real64", "Complex64") ||
        throw(ArgumentError("$source has the scalar type $(node["scalarType"]), only " *
                            "Real64 and Complex64 are supported."))
    node["elements"]["type"] in ("TextFile", "IeeeBinaryFile") ||
        throw(ArgumentError("$source has elements of type " *
                            "$(node["elements"]["type"]), only TextFile and " *
                            "IeeeBinaryFile are supported."))
    haskey(node, "nonZeroCondition") &&
        throw(ArgumentError("$source describes a block-sparse tensor, which is not " *
                            "supported. Write it densely from Cc4s."))
    metaData = something(get(node, "metaData", nothing), Dict())
    get(metaData, "halfGrid", 0) == 1 &&
        throw(ArgumentError("$source uses the halfGrid format of real Coulomb " *
                            "vertices, which is not supported."))
    return nothing
end

# Cc4s writes YAML lists as maps with the keys 0, 1, ...
_cc4s_list(node::AbstractVector) = node
_cc4s_list(node::AbstractDict) = [node[i] for i in 0:(length(node) - 1)]

function _read_cc4s_elements_binary(file, ::Type{T}, n) where {T}
    filesize(file) == n * sizeof(T) ||
        throw(DimensionMismatch("$file holds $(filesize(file)) bytes, but its " *
                                "dimensions require $n elements of $(sizeof(T)) bytes."))
    return read!(file, Vector{T}(undef, n))
end

function _read_cc4s_elements_text(file, ::Type{T}, n) where {T}
    elements = Vector{T}(undef, n)
    count = 0
    for line in eachline(file)
        isempty(strip(line)) && continue
        count += 1
        count <= n && (elements[count] = _parse_cc4s_element(T, line))
    end
    count == n ||
        throw(DimensionMismatch("$file holds $count elements, but its dimensions " *
                                "require $n."))
    return elements
end

_parse_cc4s_element(::Type{Float64}, s) = parse(Float64, s)
function _parse_cc4s_element(::Type{ComplexF64}, s)
    re, im = split(strip(s, (' ', '\t', '(', ')')), ',')
    return complex(parse(Float64, re), parse(Float64, im))
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

    tensor_data_hh = (convert(Complex{Cdouble}, Δ) for Δ in DeltaIntegralsHH)
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

    tensor_data_pphh = (convert(Complex{Cdouble}, Δ) for Δ in DeltaIntegralsPPHH)
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

    # Γ[F, n, m] in chunks of all F, which avoids a permuted copy of the whole vertex.
    # Cc4s conjugates the first orbital index n, as PsiTK does.
    tensor_data = (
        convert(Vector{Complex{Cdouble}}, ΓnmF[1, n, 1, m, :])
        for m in 1:N for n in 1:N
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

    # U[F, G], the transpose of the stored (NG, NF) matrix
    tensor_data = (
        convert(Complex{Cdouble}, coulomb_vertex_singular_vectors[iG, iF])
        for iG in 1:NG for iF in 1:NF
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
