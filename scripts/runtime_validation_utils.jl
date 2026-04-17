#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PhysiologicallyBasedDemographicModels
using Test

include(joinpath(@__DIR__, "extract_code.jl"))

function vignette_source_path(qmd_name::AbstractString)
    base = replace(basename(qmd_name), ".qmd" => "")
    nested = joinpath(@__DIR__, "..", "vignettes", base, base * ".qmd")
    legacy = joinpath(@__DIR__, "..", "vignettes", qmd_name)

    if isfile(nested)
        return nested
    elseif isfile(legacy)
        return legacy
    end

    throw(ArgumentError("Vignette not found: $qmd_name"))
end

"""
    run_vignette_runtime(qmd_name; min_figures=0)

Extract a Quarto vignette into a temporary Julia script, execute it, and
return a small summary. The vignette globals remain available in `Main` for
additional assertions in the caller.
"""
function run_vignette_runtime(qmd_name::AbstractString; min_figures::Int=0)
    qmd_path = vignette_source_path(qmd_name)

    return mktempdir() do tmpdir
        run_path = joinpath(tmpdir, replace(basename(qmd_name), ".qmd" => "_run.jl"))
        extract_julia_code(qmd_path, run_path)
        Base.include(Main, run_path)

        fig_count = if isdefined(Main, :figdir)
            vignette_figdir = getfield(Main, :figdir)
            isdir(vignette_figdir) ? length(readdir(vignette_figdir)) : 0
        else
            0
        end

        fig_count >= min_figures || error("Expected at least $min_figures figures from $qmd_name, found $fig_count")
        return (; qmd=qmd_name, fig_count=fig_count)
    end
end
