using TestItemRunner

# Default tags to exclude
excluded = [:slow]
included = Symbol[]

# Parse command line ARGS or environment variables
args = Symbol.(ARGS)
if "PSITK_TEST_ARGS" in keys(ENV) && isempty(ARGS)
    args = Symbol.(split(ENV["PSITK_TEST_ARGS"], "-"))
end

for arg in args
    if arg == :all
        empty!(excluded)
    elseif startswith(string(arg), "no")
        push!(excluded, Symbol(string(arg)[3:end])) # e.g. "noslow" -> excludes :slow
    else
        push!(included, arg)
    end
end

# If we are specifically including tags, don't arbitrarily exclude :slow unless it's explicitly in `excluded`
if length(included) > 0 && !(:slow in excluded)
    empty!(excluded)
end

function psitk_testfilter(ti)
    if any(in(ti.tags), excluded)
        return false
    elseif isempty(included)
        return true
    elseif any(in(ti.tags), included)
        return true
    else
        return false
    end
end

println("Running PsiTK tests")
if !isempty(excluded)
    println("    Excluded tags: $(join(excluded, ", "))")
end
if !isempty(included)
    println("    Included tags: $(join(included, ", "))")
end

@run_package_tests filter=psitk_testfilter verbose=true
