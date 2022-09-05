using Plots, Latexify
import Plots.plot, Plots.plot!; export plot, plot!
Plots.default(label=nothing, linewidth=2)
Plots.default(fontfamily="Computer Modern", titlefont="Computer Modern" , tickfont="Computer Modern")

function plot(res::Result; kwargs...)
    p = length(size(res.solutions)) == 1 ? plot1D(res; kwargs...) : error("only 1D plots supported now")
end

plot!(res::Result; kwargs...) = plot(res; add=true, kwargs...)


"""
$(TYPEDSIGNATURES)

Goes over a solution and an equally-sized array (a "mask") of booleans. 
true  -> solution unchanged
false -> changed to NaN (omitted from plotting)
"""
function _apply_mask(solns::Array{Vector{ComplexF64}},  booleans)
    factors = replace.(booleans, 0 => NaN)
    map(.*, solns, factors)
end


function plot1D(res::Result; x::String, y::String, class="physical", not_class=[], add=false, kwargs...)

    length(size(res.solutions)) != 1 && error("1D plots of not-1D datasets are usually a bad idea.")
    X, Y = transform_solutions(res, x), transform_solutions(res, y) # first transform, then filter
    
    Y = class !== "all" ? _apply_mask(Y, _get_mask(res, class, not_class)) : Y
    branches = _realify(Y)
          
    # start a new plot if needed
    p = add ? Plots.plot!() : Plots.plot()

    # colouring is matched to branch index - matched across plots    
    for k in findall(x -> !all(isnan.(x)), branches[1:end]) # skip NaN branches but keep indices
        l = _is_labeled(p, k) ? nothing : k
        Plots.plot!(_realify.(getindex.(X, k)), branches[k];  color=k, label=l, xlabel=latexify(x), ylabel=latexify(y), kwargs...)
    end
    
    return p
end

function _realify(x)
    !is_real(x) && !isnan(x) ? (@warn "Values with non-negligible complex parts have been projected on the real axis!") : nothing
    real(x)
end

_realify(v::Vector) = [_realify.(getindex.(v, k)) for k in 1:length(v[1])]

_vectorise(s::Vector) = s
_vectorise(s) = [s]


# returns an array of bools to mark solution which fall into classes but not not_classes
function _get_mask(res, classes, not_classes)
    bools = vcat([res.classes[c] for c in _vectorise(classes)], [map(.!, res.classes[c]) for c in _vectorise(not_classes)])
    map(.*, bools...)
end

# return true if p already has a label for branch index idx
_is_labeled(p::Plots.Plot, idx::Int64) = in(string(idx), [sub[:label] for sub in p.series_list])







