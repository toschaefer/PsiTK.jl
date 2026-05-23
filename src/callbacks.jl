
"""
Default callback function to dump status of the DFTK.LOBPCG function

# Arguments
- `thresh`: the target threshold for the LOBPCG solver
- `description`: a prefix string for the output
"""
function make_lobpcg_callback(thresh; description=nothing)
    start_time = time()

    prefix = ""
    if !isnothing(description)
        prefix = "$description | "
    end

    return function(info)
        niter = info.niter
        nlocked = info.nlocked
        n_conv_check = info.n_conv_check
        resid_norm = norm(info.resid_history[1:n_conv_check, niter+1])
        time_str = TimerOutputs.prettytime((time()-start_time)*1e9)
        @printf("\r\e[2K") # Carriage return
        @printf(
            "%sIteration: %d | Converged = %d / %d | Residual = %.2e (target: %.1e) | Elapsed time = %s", 
            prefix, 
            niter, 
            nlocked, 
            n_conv_check, 
            resid_norm, 
            thresh, 
            time_str)

        # line break if LOBPCG finished
        if nlocked >= n_conv_check 
            @printf("\n") 
        end
    end
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
