#!/usr/bin/env julia

using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const VIGNETTE_ROOT = joinpath(ROOT, "vignettes")

function vignette_sources()
    dirs = filter(name -> isdir(joinpath(VIGNETTE_ROOT, name)) && occursin(r"^\d{2}_", name),
                  readdir(VIGNETTE_ROOT))
    sort!(dirs)
    return [joinpath(VIGNETTE_ROOT, d, d * ".qmd") for d in dirs if isfile(joinpath(VIGNETTE_ROOT, d, d * ".qmd"))]
end

function render_vignette(qmd::String)
    quarto = Sys.which("quarto")
    quarto === nothing && error("quarto is required to render vignettes")
    run(`$quarto render $qmd`)
end

function main()
    failures = String[]
    for qmd in vignette_sources()
        @printf("=== Rendering %s ===\n", basename(dirname(qmd)))
        try
            render_vignette(qmd)
        catch err
            push!(failures, string(basename(dirname(qmd)), ": ", err))
        end
    end

    if isempty(failures)
        println("All vignette renders completed successfully.")
        return
    end

    println("\nRender failures:")
    for failure in failures
        println(" - ", failure)
    end
    error("One or more vignette renders failed")
end

main()
