using Test

@testset "vignette outputs — rendered .md sanity checks" begin
    vdir = joinpath(@__DIR__, "..", "vignettes")
    error_markers = [
        "ERROR:",
        "LoadError:",
        "MethodError",
        "UndefVarError",
        "DimensionMismatch",
        "Stacktrace:",
    ]

    rendered = String[]
    for d in readdir(vdir; join=true)
        isdir(d) || continue
        basename(d) in ("_freeze", ".quarto", "figures") && continue
        name = basename(d)
        md = joinpath(d, "$(name).md")
        isfile(md) && push!(rendered, md)
    end

    @test length(rendered) >= 30  # at least most vignettes have rendered md

    for md in rendered
        content = read(md, String)
        @testset "$(basename(md))" begin
            # Basic sanity: rendered file should not be trivially empty.
            @test length(content) > 500

            # Rendered output should not carry Julia error tracebacks.
            lines = split(content, '\n')
            in_output = false
            found_error = false
            for line in lines
                if startswith(strip(line), "```")
                    in_output = !in_output
                    continue
                end
                if in_output
                    for marker in error_markers
                        if occursin(marker, line)
                            found_error = true
                            break
                        end
                    end
                end
                found_error && break
            end
            @test !found_error
        end
    end
end
