
function compute_mp2_dummy(
    occ_space::OrbitalSpace,
    virt_space::OrbitalSpace,
    ΓhpF::AbstractArray,
    kernel_fourier::AbstractVector
)
    # We currently focus only on one k-point (Gamma-point only)
    ik = 1
    n_bands_occ = size(occ_space.ψ[ik], 2)
    n_bands_virt = size(virt_space.ψ[ik], 2)
    mp2_pair_energies = zeros(Float64, n_bands_occ, n_bands_occ)

    εh = occ_space.eigenvalues[ik]
    εp = virt_space.eigenvalues[ik]

    for i = 1:n_bands_occ, j = 1:n_bands_occ
        Δ_ijpp = [
            (εh[i] + εh[j] - εp[a] - εp[b]) 
            for a=1:n_bands_virt, b=1:n_bands_virt
        ]
        
        v_ijpp =  real(ΓhpF[ik,i,ik,:,:]) * real(conj.(ΓhpF[ik,j,ik,:,:]))
        v_ijpp += imag(ΓhpF[ik,i,ik,:,:]) * imag(conj.(ΓhpF[ik,j,ik,:,:]))

        v_jipp =  real(ΓhpF[ik,j,ik,:,:]) * real(conj.(ΓhpF[ik,i,ik,:,:]))
        v_jipp += imag(ΓhpF[ik,j,ik,:,:]) * imag(conj.(ΓhpF[ik,i,ik,:,:]))
        
        t_ijpp = v_ijpp ./ Δ_ijpp

        #v_ppij = conj.(v_ijpp)
        #v_ppji = conj.(v_jipp)

        #mp2_direct   =  2 * dot(t_ijpp, v_ijpp)
        #mp2_exchange = -1 * dot(t_ijpp, v_jipp)

        mp2_direct   =  2 * (v_ijpp ⋅ t_ijpp)
        mp2_exchange = -1 * (v_jipp ⋅ t_ijpp)

        mp2_pair_energies[i,j] = mp2_direct + mp2_exchange
    end
    mp2_pair_energies
end
