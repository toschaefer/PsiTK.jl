"""
    ShowProgress(; desc="Progress")

Callback printing a progress bar for step-wise computations (e.g.
[`compute_overlap_densities`](@ref)). The callback is called with a NamedTuple
`info` providing at least `info.step` and `info.total_steps`; the bar is created lazily
on the first call.
"""
struct ShowProgress
    desc::String
    progress::Ref{Union{Nothing, Progress}}
end
function ShowProgress(; desc = "Progress")
    return ShowProgress(desc, Ref{Union{Nothing, Progress}}(nothing))
end

function (callback::ShowProgress)(info)
    if isnothing(callback.progress[])
        callback.progress[] = Progress(
            info.total_steps;
            desc = callback.desc,
            dt = 0.5,
            barlen = 20,
            barglyphs = BarGlyphs(' ', '━', '╸', '─', ' '),
            color = :normal,
        )
    end
    next!(callback.progress[])
end
