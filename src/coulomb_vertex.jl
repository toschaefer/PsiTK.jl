@doc raw"""
    compute_overlap_densities(
        bra_space::OrbitalSpace,
        ket_space::OrbitalSpace;
        n_bands_bra=size(bra_space.ψ[1], 2),
        n_bands_ket=size(ket_space.ψ[1], 2),
        Ecut_ratio=1.0,
        callback=identity
    )
    compute_overlap_densities(space::OrbitalSpace; n_bands=size(space.ψ[1], 2), kwargs...)

Compute the overlap densities in reciprocal space
```math
ρ_{mn \bm G} = \int_Ω \; \psi_{m}(\bm r)^∗ \psi_{n}(\bm r)  \; e^{-i\bm r \bm G}  \; d^3 r
```

# Arguments
- `bra_space`: the bra orbital space (e.g. occupied space)
- `ket_space`: the ket orbital space (e.g. virtual space)
- `n_bands_bra`: number of bands to be considered from bra_space
- `n_bands_ket`: number of bands to be considered from ket_space
- `Ecut_ratio`: ratio to reduce the plane-wave cutoff for the densities (default: 1.0),
  `nothing` keeps the full plane-wave grid
- `callback`: called after each orbital pair with `(; step, total_steps)`,
  e.g. `callback=ShowProgress()` for a progress bar (default: no output)

# Returns
A tuple `(ρmnG, G_vectors)`:
- `ρmnG`: the overlap densities as a tensor of shape `(nk, n_bands_bra, nk, n_bands_ket, nG)`.
- `G_vectors`: the corresponding plane-wave vectors.
"""
function compute_overlap_densities(
    bra_space::OrbitalSpace,
    ket_space::OrbitalSpace;
    n_bands_bra = size(bra_space.ψ[1], 2),
    n_bands_ket = size(ket_space.ψ[1], 2),
    Ecut_ratio = 1.0,
    callback = identity,
)
    basis = bra_space.basis
    G_indices = _G_indices_within_cutoff(basis, Ecut_ratio)
    ρmnG = _compute_overlap_densities(
        basis,
        bra_space.ψ,
        ket_space.ψ;
        n_bands_bra,
        n_bands_ket,
        G_indices,
        callback,
    )

    return ρmnG, G_vectors(basis, basis.kpoints[1])[G_indices]
end

function compute_overlap_densities(space::OrbitalSpace; n_bands = size(space.ψ[1], 2), kwargs...)
    return compute_overlap_densities(space, space; n_bands_bra=n_bands, n_bands_ket=n_bands, kwargs...)
end

# Indices of the G vectors of the first k-point within the reduced cutoff Ecut * Ecut_ratio
# (Gamma-only for now)
function _G_indices_within_cutoff(basis, Ecut_ratio)
    Gs = G_vectors(basis, basis.kpoints[1])
    isnothing(Ecut_ratio) && return eachindex(Gs)
    recip_lattice = basis.model.recip_lattice
    Ecut_reduced = basis.Ecut * Ecut_ratio
    return findall(G -> sum(abs2, recip_lattice * G) / 2 <= Ecut_reduced, Gs)
end

@doc raw"""
    compute_coulomb_vertex(
        bra_space::OrbitalSpace,
        ket_space::OrbitalSpace;
        interaction_kernel=DFTK.Coulomb(DFTK.ProbeCharge()),
        n_bands_bra=size(bra_space.ψ[1], 2),
        n_bands_ket=size(ket_space.ψ[1], 2),
        Ecut_ratio=2/3,
        callback=identity
    )
    compute_coulomb_vertex(space::OrbitalSpace; n_bands=size(space.ψ[1], 2), kwargs...)

Compute the Coulomb vertex
```math
Γ_{mn \bm G} = \sqrt{v(\bm G)} \; ρ_{mn \bm G}
```
where $ρ_{mn \bm G}$ are the overlap densities (see [`compute_overlap_densities`](@ref))
and $v(\bm G)$ is the interaction kernel, e.g. the Coulomb potential
```math
v(\bm G) = \frac{4π}{\bm G^2}
```

# Arguments
- `bra_space`: the bra orbital space (e.g. occupied space)
- `ket_space`: the ket orbital space (e.g. virtual space)
- `interaction_kernel`: the DFTK interaction kernel to use (default: Coulomb)
- `n_bands_bra`: number of bands to be considered from bra_space
- `n_bands_ket`: number of bands to be considered from ket_space
- `Ecut_ratio`: ratio to reduce the plane-wave cutoff for the vertex (default: 2/3)
- `callback`: called after each orbital pair with `(; step, total_steps)`,
  e.g. `callback=ShowProgress()` for a progress bar (default: no output)

# Returns
A tuple `(ΓmnG, G_vectors, kernel_fourier)`:
- `ΓmnG`: the Coulomb vertex tensor in the uncompressed plane-wave basis.
- `G_vectors`: the corresponding plane-wave vectors.
- `kernel_fourier`: the evaluated interaction kernel at the returned G vectors.
"""
function compute_coulomb_vertex(
    bra_space::OrbitalSpace,
    ket_space::OrbitalSpace;
    interaction_kernel = DFTK.Coulomb(DFTK.ProbeCharge()),
    n_bands_bra = size(bra_space.ψ[1], 2),
    n_bands_ket = size(ket_space.ψ[1], 2),
    Ecut_ratio = 2/3,
    callback = identity,
)
    ρmnG, G_vectors = compute_overlap_densities(
        bra_space,
        ket_space;
        n_bands_bra,
        n_bands_ket,
        Ecut_ratio,
        callback,
    )

    basis = bra_space.basis
    G_indices = _G_indices_within_cutoff(basis, Ecut_ratio)
    kernel_fourier = DFTK.compute_kernel_fourier(interaction_kernel, basis)[G_indices]

    # Γ = √v ⊙ ρ along the G axis; ρmnG is not needed anymore, so scale in place
    ΓmnG = ρmnG
    ΓmnG .*= reshape(sqrt.(kernel_fourier), 1, 1, 1, 1, :)

    return ΓmnG, G_vectors, kernel_fourier
end

function compute_coulomb_vertex(space::OrbitalSpace; n_bands = size(space.ψ[1], 2), kwargs...)
    return compute_coulomb_vertex(space, space; n_bands_bra=n_bands, n_bands_ket=n_bands, kwargs...)
end

# This function initially based on code of the experimental "cc4s" branch in DFTK written by Michael Herbst
function _compute_overlap_densities(
    basis,
    ψ_bra::AbstractVector{<:AbstractArray{T}},
    ψ_ket::AbstractVector{<:AbstractArray{T}};
    n_bands_bra = size(ψ_bra[1], 2),
    n_bands_ket = size(ψ_ket[1], 2),
    G_indices = eachindex(G_vectors(basis, basis.kpoints[1])),
    callback = identity,
) where {T}
    kpt = basis.kpoints[1]
    n_kpt = length(basis.kpoints)

    # === Create index to map each stored G to -G on the full grid ===
    Gs = G_vectors(basis, kpt)
    G_to_idx = Dict(Gs[i] => i for i in eachindex(Gs))
    idx_minus_G = [G_to_idx[-Gs[i]] for i in G_indices]

    # allocate overlap densities
    ρmnG = zeros(complex(T), n_kpt, n_bands_bra, n_kpt, n_bands_ket, length(G_indices))

    is_symmetric = (ψ_bra === ψ_ket)
    if is_symmetric
        total_steps = (n_bands_bra*(n_bands_bra+1)÷2)*n_kpt^2 # only upper triangle of ρmnG
    else
        total_steps = n_bands_bra * n_bands_ket * n_kpt^2
    end
    step = 0

    # TODO:
    # Idea is to make some outer loop over the m-slices
    # for m_slice in m_slices
    #     precalculate ψmk_real for all m in this slice
    # end

    # === Calculate overlap densities ρmnG ===
    @views for (ikn, kptn) in enumerate(basis.kpoints), n = 1:n_bands_ket
        # Prepare ψnk(r)
        ψnk_real = ifft(basis, kptn, ψ_ket[ikn][:, n])

        for (ikm, kptm) in enumerate(basis.kpoints)
            for m = 1:n_bands_bra
                # Compute upper triangle only (m <= n) if spaces are symmetric
                # The lower triangle is filled via Hermitian conjugation below.
                if is_symmetric && m > n
                    continue
                end

                # Prepare ψmk(r)
                # TODO: pre-calculate some of them (not all because virtual space can be large)
                ψmk_real = ifft(basis, kptm, ψ_bra[ikm][:, m])

                # Calcualte overlap density ρ_nm(r) = ψm*(r)ψn(r) and FFT to reciprocal space
                overlap_density = fft(basis, kptn, conj.(ψmk_real) .* ψnk_real)

                # store entry of the overlap densities
                ρmnG[ikm, m, ikn, n, :] .= overlap_density[G_indices]

                # Fill lower triangle via ρmn(-G) = conjg(ρnmG)
                if is_symmetric && m != n
                    ρmnG[ikn, n, ikm, m, :] .= conj.(overlap_density[idx_minus_G])
                end

                step += 1
                callback((; step, total_steps))
            end
        end
    end
    ρmnG
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

Returns a tuple `(ΓmnF, coulomb_vertex_singular_vectors)`, where `coulomb_vertex_singular_vectors` is the applied transformation matrix.
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
        return ΓmnG, LinearAlgebra.I(size(Γmat, 2))
    else
        ΓmnF = Γmat * U[:, 1:NF]           # rotate
        return reshape(ΓmnF, size(ΓmnG)[1:4]..., NF), U[:, 1:NF]
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

Returns a tuple `(ΓmnF, coulomb_vertex_singular_vectors)`, where `coulomb_vertex_singular_vectors` is the effective transformation matrix.

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
        return ΓmnG, LinearAlgebra.I(size(Γmat, 2))
    else
        coulomb_vertex_singular_vectors = Q * U[:, 1:NF]
        ΓmnF = Γ_proj * U[:, 1:NF]           # rotate
        return reshape(ΓmnF, size(ΓmnG)[1:4]..., NF), coulomb_vertex_singular_vectors
    end
end
