"""
Abstract parent type for all algorithm states for callbacks.
Any new algorithm can define a struct that inherits from this.
"""
abstract type AbstractAlgoInfo end


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
        return function(info::AbstractAlgoInfo)
            for cb in valid_callbacks
                cb(info)
            end
        end
    end
end



