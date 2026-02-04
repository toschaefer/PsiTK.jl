raw"""
Compute the Coulomb vertex
```math
Γ_{\bm km, \bm kn \bm G} = \int_Ω \sqrt{\frac{4π}{\bm G^2}} e^{-i\bm r \bm G} \psi_{\bm km}^∗ \psi_{\bm k' n} dr
```
where `n_bands` is the number of bands to be considered.
"""
function compute_coulomb_vertex(scfres::NamedTuple, callback=nothing)
    basis = scfres.basis
    n_bands = scfres.n_bands_converge
    nk = length(basis.kpoints) 

    @withprogress name="Compute Coulomb Vertex" begin 
        # === set up callback ===
        total_steps = (n_bands*(n_bands+1)÷2)*nk^2 # only upper triangle of ΓmnG
        internal_callback = make_coulomb_vertex_callback(total_steps, @progressid)
        full_callback = chain_callbacks(internal_callback, callback)

        # === compute Coulomb Vertex ===
        _compute_coulomb_vertex(basis, scfres.ψ; n_bands, callback=full_callback)
    end
end

function _compute_coulomb_vertex(
    basis,
    ψ::AbstractVector{<:AbstractArray{T}};
    n_bands=size(ψ[1], 2),
    callback=nothing
) where {T}
    kpt   = basis.kpoints[1]
    n_G   = length(G_vectors(basis, kpt))
    n_kpt = length(basis.kpoints)

    # === Create index to map G to -G ===
    Gs  = G_vectors(basis, kpt)
    G_to_idx    = Dict(Gs[i] => i for i in 1:n_G)
    idx_minus_G = [G_to_idx[-Gs[i]] for i in 1:n_G]

    # allocate coulomb vertex
    ΓmnG  = zeros(complex(T), n_kpt, n_bands, n_kpt, n_bands, n_G)

    # Callback
    step_counter = 0
    if !isnothing(callback)
        callback(CoulombVertexInfo(step_counter))
    end

    # === Calculate Coulomb Vertex ΓmnG ===
    @views for (ikn, kptn) in enumerate(basis.kpoints), n = 1:n_bands
        # Prepare ψnk(r)
        ψnk_real = ifft(basis, kptn, ψ[ikn][:, n])

        for (ikm, kptm) in enumerate(basis.kpoints)
            # Compute momentum transfer q and Coulomb kernel
            q = kptn.coordinate - kptm.coordinate
            kernel_sqrt = sqrt.(DFTK.compute_coulomb_kernel(basis; q))

            for m in 1:n_bands
                # Compute upper triangle only (m <= n)
                # The lower triangle is filled via Hermitian conjugation below.
                if m > n continue end 

                # Prepare ψmk(r)
                ψmk_real = ifft(basis, kptm, ψ[ikm][:, m])

                # Calcualte overlap density ρ_nm(r) = ψm*(r)ψn(r) and FFT to reciprocal space
                overlap_density = fft(basis, kptn, conj(ψmk_real) .* ψnk_real)

                # store entry of Coulomb Vertex 
                ΓmnG[ikm, m, ikn, n, :] .= kernel_sqrt .* overlap_density

                # Fill lower triangle via Γmn(-G) = conjg(ΓnmG)
                if m != n
                    value = ΓmnG[ikm, m, ikn, n, :]
                    ΓmnG[ikn, n, ikm, m, :] .= conj.(value[idx_minus_G[:]])
                end
                
                # Callback
                step_counter += 1
                if !isnothing(callback)
                    callback(CoulombVertexInfo(step_counter))
                end
            end  
        end 
    end  
    ΓmnG
end


# Info struct and callback generator for compute_coulomb_vertex
struct CoulombVertexInfo <: AbstractAlgoInfo 
    step::Int
end
function make_coulomb_vertex_callback(total_steps, progress_id)
    return function(info::CoulombVertexInfo)
        fraction = info.step / total_steps
        @logprogress fraction message="Computing Coulomb Vertex..." _id=progress_id
    end
end


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



# Info struct and callback generator for compress_coulomb_vertex
struct CompressCoulombVertexInfo <: AbstractAlgoInfo 
    iter::Int
end
function make_compress_vertex_callback(progress_id)
    fraction = NaN # indeterminate state (required iterations unknown)
    @logprogress fraction message="start compression" _id=progress_id
    return function(info::CompressCoulombVertexInfo)
        #message = "Iter: $(info.iter) | Time: $(round(info.time_taken, digits=4)) s"
        message = "Iter: $(info.iter)"
        @logprogress fraction message _id=progress_id
    end
end
