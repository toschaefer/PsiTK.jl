module PsiTKJLD2Ext

using DFTK
using JLD2
using PsiTK

function PsiTK.save_orbitals(::Val{:jld2}, filename::AbstractString, space::OrbitalSpace)
    JLD2.jldopen(filename * ".new", "w") do jld
        jld["ψ"] = space.ψ
        jld["eigenvalues"] = space.eigenvalues
        jld["occupation"] = space.occupation
        jld["εF"] = space.εF
        jld["is_orthonormal"] = space.is_orthonormal
        for (key, value) in pairs(basis_fingerprint(space.basis))
            jld[string(key)] = value
        end
    end
    mv(filename * ".new", filename; force=true)
    return nothing
end

function PsiTK.load_orbitals(::Val{:jld2}, filename::AbstractString, basis::PlaneWaveBasis)
    return JLD2.jldopen(filename, "r") do jld
        for (key, value) in pairs(basis_fingerprint(basis))
            saved = jld[string(key)]
            saved == value ||
                throw(ArgumentError("The basis does not match the one the orbitals in " *
                                    "$filename were saved with: $key is $saved in the " *
                                    "file, but $value in the basis."))
        end
        OrbitalSpace(basis, jld["ψ"], jld["eigenvalues"], jld["occupation"], jld["εF"],
                     jld["is_orthonormal"])
    end
end

# Properties of the basis the orbital coefficients refer to, stored instead of the basis
function basis_fingerprint(basis::PlaneWaveBasis)
    return (; basis.Ecut, basis.fft_size, basis.kcoords_global,
            Npw=[length(G_vectors(basis, kpt)) for kpt in basis.kpoints])
end

end
