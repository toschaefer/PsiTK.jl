using ProgressMeter

raw"""
Compute the Coulomb vertex
```math
Γ_{km,\tilde{k}n,G} = ∫_Ω \sqrt{\frac{4π}{|G|^2} e^{-ir\cdot G} ψ_{km}^∗ ψ_{\tilde{k}n} dr
```
where `n_bands` is the number of bands to be considered.
"""
function compute_coulomb_vertex(
    basis,
    ψ::AbstractVector{<:AbstractArray{T}};
    n_bands=size(ψ[1], 2)
) where {T}
    mpi_nprocs(basis.comm_kpts) > 1 && error("Cannot use mpi")
    if length(basis.kpoints) > 1 && basis.use_symmetries_for_kpoint_reduction
        error("Cannot use symmetries right now.")
        # This requires appropriate insertion of kweights
    end

    kpt   = basis.kpoints[1]
    n_G   = length(G_vectors(basis, kpt)) # works only for 1-kpoint
    n_kpt = length(basis.kpoints)

    # create index to map G to -G
    Gs  = G_vectors(basis, kpt)
    G_to_idx    = Dict(Gs[i] => i for i in 1:n_G)
    idx_minus_G = [G_to_idx[-Gs[i]] for i in 1:n_G]

    # allocate coulomb vertex
    ΓmnG  = zeros(complex(T), n_kpt, n_bands, n_kpt, n_bands, n_G)

    # TODO: we need a callback structure for output
    # show progress via ProgressMeter
    progress = Progress(
         n_bands*(n_bands-1)÷2*size(basis.kpoints,1)^2; 
         desc="Compute Coulomb vertices", 
         dt=0.5, 
         barlen=20, 
         color=:black
    )

    @views for (ikn, kptn) in enumerate(basis.kpoints), n = 1:n_bands
        ψnk_real = ifft(basis, kptn, ψ[ikn][:, n])

        for (ikm, kptm) in enumerate(basis.kpoints)
            q = kptn.coordinate - kptm.coordinate
            kernel_sqrt = sqrt.(DFTK.compute_coulomb_kernel(basis; q))

            for m in 1:n_bands
                if m > n continue end # only for Gamma-only

                ψmk_real = ifft(basis, kptm, ψ[ikm][:, m])

                # kptn has to be some qptn (but works for Gamma-only)
                overlap_density = fft(basis, kptn, conj(ψmk_real) .* ψnk_real)
                ΓmnG[ikm, m, ikn, n, :] .= kernel_sqrt .* overlap_density

                # TODO (Gamma only, no spin polarization)
                # fill lower triangle via ΓmnG = ΓnmG*
                if m != n
                    mn_values = ΓmnG[ikm, m, ikn, n, :]
                    ΓmnG[ikn, n, ikm, m, :] .= conj.(mn_values[idx_minus_G[:]])
                end
                
                next!(progress) # update ProgressMeter
            end  
        end 
    end  
    ΓmnG
end
function compute_coulomb_vertex(scfres::NamedTuple)
    compute_coulomb_vertex(scfres.basis, scfres.ψ; n_bands=scfres.n_bands_converge)
end


# CoulombGramian E(G,G') defined as E = -Γ^† * Γ, where Γ(ab,G) = <a|G|b>
# see Eq. (10) in Hummel et al., JCTC (doi.org/10.1063/1.4977994).
# The operator CoulombGramian enables efficient application to a vector
# E*v = -Γ^† * (Γ*v) without full construction of E(G,G')
# in order to diagonalize E through iterative methods.
struct CoulombGramian{T}
    Γmat::T
end
function LinearAlgebra.mul!(Y, op::CoulombGramian, X)
    T = eltype(op)
    Ywork = zeros(T, size(op.Γmat,1), size(X,2))
    mul!(Ywork, op.Γmat, X, -1.0, 0.0)
    mul!(Y, op.Γmat', Ywork)
    return Y
end
function Base.:*(op::CoulombGramian, X::AbstractMatrix)
    T_out = promote_type(eltype(op), eltype(X))
    Y = similar(X, T_out)
    mul!(Y, op, X) 
    return Y
end
Base.size(op::CoulombGramian) = (size(op.Γmat, 2), size(op.Γmat, 2))
Base.eltype(op::CoulombGramian) = eltype(op.Γmat)
LinearAlgebra.ishermitian(op::CoulombGramian) = true

function svdcompress_coulomb_vertex(
    ΓmnG::AbstractArray{T,5}; 
    thresh=1e-6 # in Hartree
) where {T}
    Γmat = reshape(ΓmnG, prod(size(ΓmnG)[1:4]), size(ΓmnG, 5))
    NFguess = round(Int, 10*size(Γmat,1)^0.5)
    NG = size(Γmat,2)
    ϕk = randn(ComplexF64, NG, NFguess)
    for a in 1:NFguess
        ϕk[:,a] ./= norm(ϕk[:,a]) # normalize
    end
    
    E_GG = CoulombGramian(Γmat)

    lobpcg_thresh = thresh/10 # thresh for LOBPCG should be smaller than thresh
    println("Compress Coulomb vertices.")
    @time FFF = DFTK.LOBPCG(E_GG, ϕk, I, I, lobpcg_thresh, 500, display_progress=true)
    nkeep = findlast(s -> abs(s) > thresh, FFF.λ)
    println("Compressed Coulomb vertices from NG=$(NG) to NF=$(nkeep).")
    @views if isnothing(nkeep)
        ΓmnG
    else
        cΓmat = Γmat * FFF.X[:, 1:nkeep] 
        reshape(cΓmat, size(ΓmnG)[1:4]..., nkeep)
    end
     
    #@time U, S, V = tsvd(Γmat, NFguess; tolconv=tolconv, maxiter=maxiter)
    ##@time res = eigen(Γmat' * Γmat)
    #@time F = svd(Γmat)

    #Serror = abs.(S - F.S[1:NFguess])
    ##println(Serror)
    ##println(" ")
    #println("max error at: ", findmax(Serror))

    #tol = sqrt(thresh) # singular values are sqrt of energies
    #nkeep = findlast(s -> abs(s) > tol, F.S)
    #@views if isnothing(nkeep)
    #    ΓmnG
    #else
    #    cΓmat = F.U[:, 1:nkeep] * Diagonal(F.S[1:nkeep])
    #    reshape(cΓmat, size(ΓmnG)[1:4]..., nkeep)
    #end
end
