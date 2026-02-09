@doc raw"""
Compute the Coulomb vertex
```math
Γ_{\bm km, \bm kn \bm G} = \int_Ω \sqrt{\frac{4π}{\bm G^2}} e^{-i\bm r \bm G} \psi_{\bm km}^∗ \psi_{\bm k' n} dr
```
where `n_bands` is the number of bands to be considered.
"""
function compute_coulomb_vertex(scfres::NamedTuple)
    basis = scfres.basis
    n_bands = scfres.n_bands_converge
    nk = length(basis.kpoints) 

    # === set up callback ===
    total_steps = (n_bands*(n_bands+1)÷2)*nk^2 # only upper triangle of ΓmnG
    callback = make_coulomb_vertex_callback(total_steps)

    # === compute Coulomb Vertex ===
    _compute_coulomb_vertex(basis, scfres.ψ; n_bands, callback)
end

# This function initially based on code of the experimental "cc4s" branch in DFTK written by Michael Herbst
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

    # TODO:
    # Idea is to make some outer loop over the m-slices
    # for m_slice in m_slices
    #     precalculate ψmk_real for all m in this slice
    # end

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
                # TODO: pre-calculate some of them (not all because virtual space can be large)
                ψmk_real = ifft(basis, kptm, ψ[ikm][:, m])

                # Calcualte overlap density ρ_nm(r) = ψm*(r)ψn(r) and FFT to reciprocal space
                overlap_density = fft(basis, kptn, conj.(ψmk_real) .* ψnk_real)

                # store entry of Coulomb Vertex 
                ΓmnG[ikm, m, ikn, n, :] .= kernel_sqrt .* overlap_density

                # Fill lower triangle via Γmn(-G) = conjg(ΓnmG)
                if m != n
                    #value = ΓmnG[ikm, m, ikn, n, :]
                    #ΓmnG[ikn, n, ikm, m, :] .= conj.(@view value[idx_minus_G])
                    ΓmnG[ikn, n, ikm, m, :] .= conj.(view(ΓmnG, ikm, m, ikn, n, idx_minus_G))
                end
                
                # Callback
                if !isnothing(callback)
                    callback()
                end
            end  
        end 
    end  
    ΓmnG
end

function make_coulomb_vertex_callback(total_steps)
    p = Progress(
        total_steps; 
        desc="Computing Coulomb Vertex",
        dt=0.5,
        barlen=20,
        barglyphs=BarGlyphs(' ', '━', '╸', '─', ' '),
        color=:normal
    )
    return function()
        next!(p)
    end
end



@doc raw"""
    ΓCompressionStrategy

Abstract type for different strategies of compressing the 
Coulomb vertex $\Gamma_{nm}^G$.

Available models:
- [`CoulombGramian`](@ref) 
- [`AdaptiveRandomizedSVD`](@ref) (default)
"""
abstract type ΓCompressionStrategy end


@doc raw"""
    compress_coulomb_vertex(ΓmnG; thresh=1e-6, 

TODO

"""
function compress_coulomb_vertex(
    ΓmnG::AbstractArray{T,5}; 
    thresh=1e-6, # in Hartree
    compression_strategy::ΓCompressionStrategy=AdaptiveRandomizedSVD()
) where {T}
    _compress_coulomb_vertex(ΓmnG, thresh, compression_strategy)
end


@doc raw"""
    CoulombGramian <: ΓCompressionStrategy

AdaptiveRandomizedSVD
"""
struct CoulombGramian <: ΓCompressionStrategy end
function _compress_coulomb_vertex(
    ΓmnG::AbstractArray{T,5},
    thresh, # in Hartree
    ::CoulombGramian
) where {T}
    Γmat = reshape(ΓmnG, prod(size(ΓmnG)[1:4]), size(ΓmnG, 5))
    Npp, NG = size(Γmat)

    E_GG = -Hermitian(Γmat' * Γmat)         # Gramian in full PW basis
    λ, U = eigen(E_GG)                      # diagonalize
    NF = findlast(s -> abs(s) > thresh, λ)  # truncate based on thresh
    if isnothing(NF)
        ΓmnG
    else
        ΓmnF = Γmat * U[:, 1:NF]           # rotate
        reshape(ΓmnF, size(ΓmnG)[1:4]..., NF)
    end
end


@doc raw"""
    AdaptiveRandomizedSVD <: ΓCompressionStrategy

Compressing the Coulomb vertex $\Gamma_{mn}^{G}$ through an adaptive randomized SVD.

This algorithm approximates the range of the row space of $\Gamma$ (the orbital indices
are considered as superindex) through a thin basis Q, such that
```math
\Gamma \approx \Gamma Q Q^\dagger
```
where $\Gamma$ is a $N_{pp} \times N_G$ and $Q$ a $N_G \times N_F$ matrix.
This is done through a stochastic Q and a diagonalization of
```math
E = -\tilde \Gamma^\dagger \tilde \Gamma = U D U^\dagger
```    
where $\tilde \Gamma = \Gamma Q$. 
The compressed $\Gamma$ is then obtained via $\Gamma_\text{compressed} = \tilde \Gamma U$.

The dimension $N_F$ is found by a preceding adaptive range finder. 
This finder iteratively increases the columns of Q (i.e. $N_F$) in steps of $\sqrt{N_{pp}}$ 
and stops when the error for a stochastic test vector $\omega$
```math
\varepsilon = \frac{ \Vert (1 - QQ^\dagger)\Gamma^\dagger \omega \Vert}{\Vert \omega \Vert}
```
is smaller than thresh/10.
"""
struct AdaptiveRandomizedSVD <: ΓCompressionStrategy end
function _compress_coulomb_vertex(
    ΓmnG::AbstractArray{T,5},
    thresh, # in Hartree
    ::AdaptiveRandomizedSVD
) where {T}
    Γmat = reshape(ΓmnG, prod(size(ΓmnG)[1:4]), size(ΓmnG, 5))
    Npp, NG = size(Γmat)
    
    # === Adaptive Range Finder for NF ===
    
    # Matrix for the stochastic guess basis (intially empty)
    Q = Matrix{T}(undef, NG, 0) 

    # Step size for increasing the basis = (√Npp)/10
    column_block_size = round(Int, Npp^0.5)  
    
    # Stochastic test vector for error estimation
    ω = randn(T, Npp) 
    ω_norm = norm(ω)
    
    # target error a little smaller than √thresh
    target_error = sqrt(thresh)/10

    # set current error initially larger than stop criterion 
    current_error = 2 * target_error
    
    # Iterate until convergence
    while current_error > target_error && size(Q, 2) < NG
        Ω = randn(T, Npp, column_block_size) # D@doc raw a new random block
        Y_block = Γmat' * Ω                  # Project Γ onto Ω
        
        # Orthogonalize Y_block against existing Q (Gram-Schmidt)
        if size(Q, 2) > 0
            coeffs = Q' * Y_block
            Y_block .-= Q * coeffs
        end
        
        Q_block = Matrix(qr(Y_block).Q) # Orthonormalize block itself (QR)
        Q = hcat(Q, Q_block)            # Update stochastic basis Q
        
        # current_error = || Γ' * ω - Q * (Q' * Γ' * ω) || / ||ω||
        proj_ω = Γmat' * ω
        coeffs_ω = Q' * proj_ω
        rem_ω = proj_ω - Q * coeffs_ω
        current_error = norm(rem_ω) / ω_norm
    end

    # === Compression Step ===
    Γ_proj = Γmat * Q                       # Project Γ onto Q 
    E_GG = -Hermitian(Γ_proj' * Γ_proj)     # Gramian in Q basis
    λ, U = eigen(E_GG)                      # diagonalize
    NF = findlast(s -> abs(s) > thresh, λ)  # truncate based on thresh
    if isnothing(NF)
        ΓmnG
    else
        ΓmnF = Γ_proj * U[:, 1:NF]           # rotate
        reshape(ΓmnF, size(ΓmnG)[1:4]..., NF)
    end
end


