"""
Combines multiple callbacks into one.
"""
function chain_callbacks(callbacks...)
    # Filter out 'nothing' to avoid errors
    valid_callbacks = filter(!isnothing, callbacks)
    
    if isempty(valid_callbacks)
        return nothing
    elseif length(valid_callbacks) == 1
        return first(valid_callbacks)
    else
        # Return a functor that iterates through them
        return function(info)
            for cb in valid_callbacks
                cb(info)
            end
        end
    end
end


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
