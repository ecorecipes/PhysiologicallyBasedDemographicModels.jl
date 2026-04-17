#!/usr/bin/env julia

using Printf

const ROOT = normpath(joinpath(@__DIR__, ".."))
const VIGNETTE_ROOT = joinpath(ROOT, "vignettes")
const REFERENCE_BIB = joinpath(VIGNETTE_ROOT, "references.bib")

function vignette_sources()
    dirs = filter(name -> isdir(joinpath(VIGNETTE_ROOT, name)) && occursin(r"^\d{2}_", name),
                  readdir(VIGNETTE_ROOT))
    sort!(dirs)
    return [joinpath(VIGNETTE_ROOT, d, d * ".qmd") for d in dirs if isfile(joinpath(VIGNETTE_ROOT, d, d * ".qmd"))]
end

function has_reference_path(text::String)
    occursin(r"bibliography:\s*(\.\./references\.bib|\n\s*-\s*\.\./references\.bib)", text)
end

function strip_code(text::String)
    text = replace(text, r"(?ms)^```.*?^```" => "")
    replace(text, r"`[^`\n]+`" => "")
end

function bibliography_keys()
    bib = read(REFERENCE_BIB, String)
    Set(m.captures[1] for m in eachmatch(r"@\w+\{([^,]+),", bib))
end

function is_bibliography_key(key::AbstractString)
    !any(startswith(key, prefix) for prefix in
         ("fig-", "tbl-", "eq-", "sec-", "thm-", "lem-", "cor-", "prp-", "alg-", "lst-", "apx-"))
end

function citation_keys(text::String)
    clean = strip_code(text)
    keys = Set{String}()
    for m in eachmatch(r"(?<![A-Za-z0-9_])@([A-Za-z][A-Za-z0-9:_-]+)", clean)
        key = rstrip(m.captures[1], ':')
        is_bibliography_key(key) && push!(keys, key)
    end
    keys
end

function main()
    failures = String[]
    bib_keys = bibliography_keys()
    for qmd in vignette_sources()
        text = read(qmd, String)
        cited_keys = citation_keys(text)
        has_reference_path(text) || push!(failures, "$(basename(dirname(qmd))): missing ../references.bib bibliography")
        isempty(cited_keys) && push!(failures, "$(basename(dirname(qmd))): no citation found in body")
        missing_keys = sort!(collect(setdiff(cited_keys, bib_keys)))
        isempty(missing_keys) || push!(failures,
                                       "$(basename(dirname(qmd))): missing bibliography keys: " * join(missing_keys, ", "))
    end

    if isempty(failures)
        println("All vignette metadata checks passed.")
        return
    end

    println("Metadata failures:")
    for failure in failures
        println(" - ", failure)
    end
    error("Vignette metadata checks failed")
end

main()
