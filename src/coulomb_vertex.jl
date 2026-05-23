@doc raw"""
    compute_coulomb_vertex(
        active_space::OrbitalSpace;
        n_bands=size(active_space.ψ[1], 2)
    )

Compute the Coulomb vertex
```math
Γ_{mn \bm G} = \int_Ω  \; \sqrt{v(\bm G)} \; \psi_{m}(\bm r)^∗ \psi_{n}(\bm r)  \; e^{-i\bm r \bm G}  \; d^3 r
```
where $v(\bm G)$ is the Coulomb potential, e.g.
```math
v(\bm G) = \sqrt{\frac{4π}{\bm G^2}}
```

# Arguments
- `active_space`: the active orbital space 
- `n_bands`: number of bands to be considered

"""
function compute_coulomb_vertex(
    active_space::OrbitalSpace;
    interaction_kernel = DFTK.Coulomb(DFTK.ProbeCharge()),
    n_bands = size(active_space.ψ[1], 2),
    compression = nothing,
    Ecut_ratio = 2/3,
)
    basis = active_space.basis
    nk = length(basis.kpoints)

    # === set up callback ===
    total_steps = (n_bands*(n_bands+1)÷2)*nk^2 # only upper triangle of ΓmnG
    callback = make_coulomb_vertex_callback(total_steps)

    # === compute Coulomb Vertex ===
    ΓmnG = _compute_coulomb_vertex(
        basis,
        interaction_kernel,
        active_space.ψ;
        n_bands,
        callback,
    )

    # === Filter G vectors (Ecut_ratio) ===
    if Ecut_ratio !== nothing
        # reduce the plane wave cutoff
        # this only works for Gamma-only now
        Ecut_reduced = basis.Ecut * Ecut_ratio
        kpt = basis.kpoints[1]
        model = basis.model
        G_mask =
            [sum(abs2, model.recip_lattice * G) / 2 <= Ecut_reduced for G in kpt.G_vectors]
        G_indices = findall(G_mask)
        nG_reduced = length(G_indices)
        nk1, nb1, nk2, nb2, nG = size(ΓmnG)
        ΓmnG_reduced = zeros(eltype(ΓmnG), nk1, nb1, nk2, nb2, nG_reduced)
        ΓmnG_reduced[:, :, :, :, :] = ΓmnG[:, :, :, :, G_indices]
        ΓmnG = ΓmnG_reduced
    end

    if !isnothing(compression)
        return compress_coulomb_vertex(ΓmnG, compression)
    end

    ΓmnG
end

# This function initially based on code of the experimental "cc4s" branch in DFTK written by Michael Herbst
function _compute_coulomb_vertex(
    basis,
    interaction_kernel,
    ψ::AbstractVector{<:AbstractArray{T}};
    n_bands = size(ψ[1], 2),
    callback = nothing,
) where {T}
    kpt = basis.kpoints[1]
    n_G = length(G_vectors(basis, kpt))
    n_kpt = length(basis.kpoints)

    # === Create index to map G to -G ===
    Gs = G_vectors(basis, kpt)
    G_to_idx = Dict(Gs[i] => i for i = 1:n_G)
    idx_minus_G = [G_to_idx[-Gs[i]] for i = 1:n_G]

    # allocate coulomb vertex
    ΓmnG = zeros(complex(T), n_kpt, n_bands, n_kpt, n_bands, n_G)

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
            kernel_sqrt = sqrt.(DFTK.compute_kernel_fourier(interaction_kernel, basis; q))

            for m = 1:n_bands
                # Compute upper triangle only (m <= n)
                # The lower triangle is filled via Hermitian conjugation below.
                if m > n
                    continue
                end

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
                    ΓmnG[ikn, n, ikm, m, :] .=
                        conj.(view(ΓmnG, ikm, m, ikn, n, idx_minus_G))
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





@doc raw"""
    CoulombGramian

This strategy compresses the Coulomb vertex $\Gamma_{mn}^{G}$ through the largest
eigenvalues of the Coulomb Gramian
```math
H = - \Gamma^\dagger \Gamma = U \Lambda U^\dagger
```    
The compressed $\Gamma$ is then obtained via $\Gamma_\text{compressed} = \Gamma U$,
where the columns of $U$ are restricted such that $\lambda >$ `thresh`.
"""
Base.@kwdef struct CoulombGramian
    thresh::Float64 = 1e-6
end
function compress_coulomb_vertex(
    ΓmnG::AbstractArray{T,5},
    strategy::CoulombGramian,
) where {T}
    thresh = strategy.thresh
    Γmat = reshape(ΓmnG, prod(size(ΓmnG)[1:4]), size(ΓmnG, 5))
    Npp, NG = size(Γmat)

    H = -Hermitian(Γmat' * Γmat)         # Gramian in full PW basis
    λ, U = eigen(H)                      # diagonalize
    NF = findlast(s -> abs(s) > thresh, λ)  # truncate based on thresh
    if isnothing(NF)
        ΓmnG
    else
        ΓmnF = Γmat * U[:, 1:NF]           # rotate
        reshape(ΓmnF, size(ΓmnG)[1:4]..., NF)
    end
end


@doc raw"""
    AdaptiveRandomizedSVD

This strategy compresses the Coulomb vertex $\Gamma_{mn}^{G}$ via an adaptive randomized SVD.

The algorithm approximates the range of the row space of $\Gamma$ (the orbital indices
are considered as superindex) through a thin basis Q, such that
```math
\Gamma \approx \Gamma Q Q^\dagger
```
where $\Gamma$ is a $N_{pp} \times N_G$ and $Q$ a $N_G \times N_F$ matrix.
This is done through a stochastic Q and a diagonalization of
```math
H = -\tilde \Gamma^\dagger \tilde \Gamma = U \Lambda U^\dagger
```    
where $\tilde \Gamma = \Gamma Q$. 
The compressed $\Gamma$ is then obtained via $\Gamma_\text{compressed} = \tilde \Gamma U$.

The dimension $N_F$ is found by a preceding adaptive range finder. 
This finder iteratively increases the columns of Q (i.e. $N_F$) in steps of $2\sqrt{N_{pp}}$ 
and stops when the error for a stochastic test vector $\omega$
```math
\varepsilon =  \Vert (1 - QQ^\dagger)\Gamma^\dagger \omega \Vert
```
is smaller than thresh/2.
"""
Base.@kwdef struct AdaptiveRandomizedSVD
    thresh::Float64 = 1e-6
end
function compress_coulomb_vertex(
    ΓmnG::AbstractArray{T,5},
    strategy::AdaptiveRandomizedSVD,
) where {T}
    thresh = strategy.thresh
    Γmat = reshape(ΓmnG, prod(size(ΓmnG)[1:4]), size(ΓmnG, 5))
    Npp, NG = size(Γmat)

    # === Adaptive Range Finder for NF ===

    # Store blocks of the stochastic guess basis
    Q_blocks = Matrix{T}[]

    # Step size for increasing the basis = 2*(√Npp)
    column_block_size = round(Int, 2*Npp^0.5)

    # Stochastic test vector for error estimation
    ω = randn(T, Npp)

    # target error a little smaller than √thresh
    target_error = sqrt(thresh)/2

    # set current error initially larger than stop criterion 
    current_error = 2 * target_error

    # Precompute projection of the test vector onto Γ
    proj_ω = Γmat' * ω
    rem_ω = copy(proj_ω)

    current_cols = 0

    # Iterate until convergence
    while current_error > target_error && current_cols < NG
        current_block_size = min(column_block_size, NG - current_cols)
        Ω = randn(T, Npp, current_block_size) # Draw a new random block
        Y_block = Γmat' * Ω                   # Project Γ onto Ω

        # Orthogonalize Y_block against existing Q: we do iterated Gram-Schmidt 
        # to preserve orthogonality (assuming "twice is enough" rule)
        if !isempty(Q_blocks)
            # first pass
            for Qb in Q_blocks
                coeffs1 = Qb' * Y_block
                Y_block .-= Qb * coeffs1
            end
            # second pass
            for Qb in Q_blocks
                coeffs2 = Qb' * Y_block
                Y_block .-= Qb * coeffs2
            end
        end

        Q_block = Matrix(qr(Y_block).Q) # Orthonormalize block itself (QR)
        push!(Q_blocks, Q_block)        # Update stochastic basis blocks
        current_cols += current_block_size

        # Update current_error incrementally
        coeffs_ω = Q_block' * rem_ω
        rem_ω .-= Q_block * coeffs_ω
        current_error = norm(rem_ω)
    end

    # Combine blocks to form the full basis
    Q = reduce(hcat, Q_blocks)

    # === Compression Step ===
    Γ_proj = Γmat * Q                       # Project Γ onto Q
    H = -Hermitian(Γ_proj' * Γ_proj)        # Gramian in Q basis
    λ, U = eigen(H)                         # diagonalize
    NF = findlast(s -> abs(s) > thresh, λ)  # truncate based on thresh
    if isnothing(NF)
        ΓmnG
    else
        ΓmnF = Γ_proj * U[:, 1:NF]           # rotate
        reshape(ΓmnF, size(ΓmnG)[1:4]..., NF)
    end
end
