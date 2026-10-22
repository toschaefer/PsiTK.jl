@testitem "Cc4s Tensor IO" tags=[:cc4s] begin
    using PsiTK: write_cc4s_tensor

    folder = mktempdir()
    dims(lens...) = [Dict("length" => l, "type" => "Index$d") for (d, l) in enumerate(lens)]

    # Round trips through write_cc4s_tensor, which streams the first index fastest
    A = rand(3, 5)
    write_cc4s_tensor(folder, "A", vec(A); dimensions=dims(3, 5), elementType="TextFile",
                      metaData=Dict("answer" => 42), unit=0.5)
    (; tensor, dimensions, unit, metaData) = read_cc4s_tensor(folder, "A")
    @test tensor == A && tensor isa Matrix{Float64}
    @test dimensions == dims(3, 5)
    @test unit == 0.5 && metaData == Dict("answer" => 42)

    B = rand(ComplexF64, 2, 3, 4)
    write_cc4s_tensor(folder, "B", vec(B); dimensions=dims(2, 3, 4), scalarType="Complex64")
    @test read_cc4s_tensor(folder, "B").tensor == B
    write_cc4s_tensor(folder, "Bt", vec(B); dimensions=dims(2, 3, 4),
                      scalarType="Complex64", elementType="TextFile")
    @test read_cc4s_tensor(folder, "Bt").tensor == B

    # Tensors as written by Cc4s: dimensions as maps with keys 0, 1, ..., no meta data,
    # complex text elements as (re,im) and optionally the last index fastest (rowMajor)
    C = [complex(i, -j) for i in 1:2, j in 1:3]
    function write_cc4s_style(name, rowMajor, elements)
        open(joinpath(folder, "$name.yaml"), "w") do io
            println(io, """
                dimensions:
                  0:
                    length: 2
                    type: Hole
                  1:
                    length: 3
                    type: Particle
                elements:
                  type: TextFile
                  rowMajor: $rowMajor
                scalarType: Complex64
                type: Tensor
                unit: 1
                version: 100""")
        end
        open(io -> foreach(z -> println(io, "($(real(z)),$(imag(z)))"), elements),
             joinpath(folder, "$name.elements"), "w")
    end
    write_cc4s_style("rdm.hp", false, vec(C))
    rdm = read_cc4s_tensor(folder, "rdm.hp")
    @test rdm.tensor == C && rdm.unit == 1.0 && isempty(rdm.metaData)
    @test [d["type"] for d in rdm.dimensions] == ["Hole", "Particle"]
    write_cc4s_style("rdm.hp.rowMajor", true, vec(permutedims(C)))
    @test read_cc4s_tensor(folder, "rdm.hp.rowMajor").tensor == C

    # Missing files, mismatching element counts and unsupported formats are rejected
    @test_throws ArgumentError read_cc4s_tensor(folder, "missing")
    write_cc4s_tensor(folder, "short", 1:5; dimensions=dims(2, 3), elementType="TextFile")
    @test_throws DimensionMismatch read_cc4s_tensor(folder, "short")
    write_cc4s_tensor(folder, "long", rand(7); dimensions=dims(2, 3))
    @test_throws DimensionMismatch read_cc4s_tensor(folder, "long")
    write_cc4s_tensor(folder, "half", rand(ComplexF64, 4); dimensions=dims(4),
                      scalarType="Complex64", metaData=Dict("halfGrid" => 1))
    @test_throws ArgumentError read_cc4s_tensor(folder, "half")

    # Tensor sets as written by Cc4s, e.g. the density matrix as its hole/particle blocks
    # in rdm.components.<key>.elements, with or without the type TensorSet in rdm.yaml
    Nocc, Nvirt = 2, 3
    D = rand(ComplexF64, Nocc + Nvirt, Nocc + Nvirt)
    states = Dict('h' => 1:Nocc, 'p' => (Nocc + 1):(Nocc + Nvirt))
    function write_cc4s_set(name, type)
        components = Dict{String,Any}()
        for key in ("hh", "hp", "ph", "pp")
            block = D[states[key[1]], states[key[2]]]
            write(joinpath(folder, "$name.components.$key.elements"), vec(block))
            components[key] = Dict(
                "dimensions" => Dict(0 => Dict("length" => size(block, 1)),
                                     1 => Dict("length" => size(block, 2))),
                "elements" => Dict("type" => "IeeeBinaryFile"),
                "scalarType" => "Complex64",
                "type" => "Tensor",
                "unit" => 3.388234215746401e-316,
                "version" => 100,
            )
        end
        node = Dict{String,Any}("components" => components)
        isnothing(type) || (node["type"] = type)
        open(io -> PsiTK.YAML.write(io, node), joinpath(folder, "$name.yaml"), "w")
    end
    for (name, type) in (("rdm", nothing), ("rdm_typed", "TensorSet"))
        write_cc4s_set(name, type)
        set = read_cc4s_tensor_set(folder, name)
        @test sort!(collect(keys(set))) == ["hh", "hp", "ph", "pp"]
        @test set["hp"].tensor == D[1:Nocc, (Nocc + 1):end]
        @test merge_hole_particle_blocks(set) == D
    end

    # Tensors and tensor sets are not mixed up, and blocks must be complete and consistent
    @test_throws ArgumentError read_cc4s_tensor(folder, "rdm")
    @test_throws ArgumentError read_cc4s_tensor_set(folder, "A")
    blocks = read_cc4s_tensor_set(folder, "rdm")
    @test_throws DimensionMismatch merge_hole_particle_blocks(
        merge(blocks, Dict("hp" => (; tensor=zeros(Nocc, Nocc)))),
    )
    delete!(blocks, "ph")
    @test_throws ArgumentError merge_hole_particle_blocks(blocks)
end
