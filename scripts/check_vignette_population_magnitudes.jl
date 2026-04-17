#!/usr/bin/env julia

using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const VIGNETTE_ROOT = joinpath(ROOT, "vignettes")
const DEFAULT_THRESHOLD = 1e10
const MAX_HITS_PER_FILE = 5

struct NumericHit
    line::Int
    raw::String
    value::Float64
    context::String
end

function vignette_markdown_outputs()
    dirs = filter(name -> isdir(joinpath(VIGNETTE_ROOT, name)) && occursin(r"^\d{2}_", name),
                  readdir(VIGNETTE_ROOT))
    sort!(dirs)
    return [joinpath(VIGNETTE_ROOT, d, d * ".md") for d in dirs if isfile(joinpath(VIGNETTE_ROOT, d, d * ".md"))]
end

function magnitude_threshold()
    isempty(ARGS) && return DEFAULT_THRESHOLD
    try
        return parse(Float64, ARGS[1])
    catch
        error("Invalid magnitude threshold '$(ARGS[1])'. Pass a numeric value, e.g. 1e10.")
    end
end

function extract_hits(line::String, threshold::Float64)
    hits = NumericHit[]

    for m in eachmatch(r"(?<![A-Za-z0-9_.])([-+]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?(?:[eE][+-]?\d+)?)(?!%)", line)
        raw = m.captures[1]
        value = try
            parse(Float64, replace(raw, "," => ""))
        catch
            NaN
        end
        if isfinite(value) && abs(value) >= threshold
            push!(hits, NumericHit(0, raw, abs(value), strip(line)))
        end
    end

    for m in eachmatch(r"(?<![A-Za-z0-9_.])([-+]?\d+(?:\.\d+)?)\s*[×x]\s*10\^\{?([+-]?\d+)\}?(?!%)", line)
        coeff = try
            parse(Float64, m.captures[1])
        catch
            NaN
        end
        exponent = try
            parse(Int, m.captures[2])
        catch
            0
        end
        value = coeff * 10.0^exponent
        if isfinite(value) && abs(value) >= threshold
            push!(hits, NumericHit(0, m.match, abs(value), strip(line)))
        end
    end

    hits
end

function scan_markdown(path::String, threshold::Float64)
    hits = NumericHit[]
    in_fence = false

    for (line_no, line) in enumerate(eachline(path))
        stripped = lstrip(line)
        if startswith(stripped, "```")
            in_fence = !in_fence
            continue
        end
        in_fence && continue

        for hit in extract_hits(line, threshold)
            push!(hits, NumericHit(line_no, hit.raw, hit.value, hit.context))
        end
    end

    sort!(hits, by = hit -> hit.value, rev = true)
    hits
end

function main()
    threshold = magnitude_threshold()
    offenders = Tuple{String, Vector{NumericHit}}[]

    for md in vignette_markdown_outputs()
        hits = scan_markdown(md, threshold)
        isempty(hits) || push!(offenders, (md, hits))
    end

    if isempty(offenders)
        println("All rendered vignette outputs are below the magnitude threshold of ",
                @sprintf("%.6g", threshold), ".")
        return
    end

    sort!(offenders, by = item -> first(item[2]).value, rev = true)

    println("Vignette output magnitude failures (threshold = ", @sprintf("%.6g", threshold), "):")
    for (md, hits) in offenders
        vignette = basename(dirname(md))
        println(" - ", vignette, ": ", length(hits), " hit(s), max = ", @sprintf("%.6g", first(hits).value))
        for hit in Iterators.take(hits, MAX_HITS_PER_FILE)
            println("   L", hit.line, ": ", hit.raw, " :: ", first(hit.context, 220))
        end
    end

    error("Rendered vignette outputs exceed the population magnitude threshold")
end

main()
