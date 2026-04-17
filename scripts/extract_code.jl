#!/usr/bin/env julia
"""Extract Julia code blocks from a .qmd file and write a runnable .jl script."""

function extract_julia_code(qmd_path::String, output_path::String)
    lines = readlines(qmd_path)
    out = IOBuffer()
    in_block = false
    
    vname = replace(basename(qmd_path), ".qmd" => "")
    println(out, "# Auto-extracted from $(basename(qmd_path))")
    println(out, "figdir = joinpath(@__DIR__, \"figures\", \"$(vname)\")")
    println(out, "mkpath(figdir)")
    println(out, "fig_counter = Ref(0)")
    println(out, "")
    
    for line in lines
        if startswith(strip(line), "```{julia}")
            in_block = true
            continue
        elseif in_block && startswith(strip(line), "```")
            in_block = false
            println(out, "")
            continue
        end
        if in_block
            println(out, line)
        end
    end
    
    write(output_path, String(take!(out)))
    println("Extracted code from $(basename(qmd_path)) -> $(basename(output_path))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) >= 2
        extract_julia_code(ARGS[1], ARGS[2])
    elseif length(ARGS) == 1
        outpath = replace(ARGS[1], ".qmd" => "_run.jl")
        extract_julia_code(ARGS[1], outpath)
    else
        println("Usage: julia extract_code.jl <input.qmd> [output.jl]")
    end
end
