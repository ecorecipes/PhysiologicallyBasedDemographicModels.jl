using Documenter
using PhysiologicallyBasedDemographicModels
using StructuredPopulationCore

const TUTORIAL_BASENAMES = [
    "01_getting_started",
    "02_cotton_plant",
    "03_coffee_berry_borer",
    "04_grapevine",
    "05_lobesia_overwintering",
    "06_bt_cotton_resistance",
    "07_pesticide_resistance",
    "08_screwworm_sit",
    "09_bt_cotton_india",
    "10_tsetse_ecosocial",
    "11_olive_climate",
]

function sync_tutorials_from_vignettes()
    quarto = Sys.which("quarto")
    quarto === nothing && error("Quarto is required to build the tutorial docs from vignettes/*/*.qmd")

    repo_root = normpath(joinpath(@__DIR__, ".."))
    vignette_dir = joinpath(repo_root, "vignettes")
    tutorial_dir = joinpath(@__DIR__, "src", "tutorials")
    mkpath(tutorial_dir)

    for base in TUTORIAL_BASENAMES
        qmd = joinpath(vignette_dir, base, base * ".qmd")
        md = base * ".md"
        run(`$quarto render $qmd --to gfm --output $md --output-dir $tutorial_dir`)

        rendered_in_src = joinpath(@__DIR__, "src", md)
        rendered_in_tutorials = joinpath(tutorial_dir, md)
        if isfile(rendered_in_src)
            mv(rendered_in_src, rendered_in_tutorials; force=true)
        end
    end
end

sync_tutorials_from_vignettes()

# --- HTML build ---
makedocs(;
    modules = [PhysiologicallyBasedDemographicModels, StructuredPopulationCore],
    warnonly = true,
    authors = "Simon Frost",
    sitename = "PhysiologicallyBasedDemographicModels.jl",
    remotes = nothing,
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://ecorecipes.github.io/PhysiologicallyBasedDemographicModels.jl"),
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Getting Started" => "tutorials/01_getting_started.md",
            "Cotton Plant Model" => "tutorials/02_cotton_plant.md",
            "Coffee Berry Borer" => "tutorials/03_coffee_berry_borer.md",
            "Grapevine C/N Model" => "tutorials/04_grapevine.md",
            "Lobesia Overwintering" => "tutorials/05_lobesia_overwintering.md",
            "Bt Cotton Resistance" => "tutorials/06_bt_cotton_resistance.md",
            "Pesticide Resistance Optimization" => "tutorials/07_pesticide_resistance.md",
            "Screwworm SIT Eradication" => "tutorials/08_screwworm_sit.md",
            "Indian Bt Cotton Bioeconomics" => "tutorials/09_bt_cotton_india.md",
            "Tsetse Ecosocial Model" => "tutorials/10_tsetse_ecosocial.md",
            "Olive Climate Economics" => "tutorials/11_olive_climate.md",
        ],
        "API Reference" => [
            "Types & Traits" => "api/types.md",
            "Weather & Forcing" => "api/weather.md",
            "Dynamics" => "api/dynamics.md",
            "Problem & Solution" => "api/problems.md",
            "Analysis" => "api/analysis.md",
            "Multi-Species Interactions" => "api/interactions.md",
            "Economics" => "api/economics.md",
            "Genetics" => "api/genetics.md",
            "Epidemiology" => "api/epidemiology.md",
            "Utilities" => "api/utilities.md",
        ],
    ])

# --- PDF build (LaTeX) ---
makedocs(;
    modules = [PhysiologicallyBasedDemographicModels, StructuredPopulationCore],
    warnonly = true,
    authors = "Simon Frost",
    sitename = "PhysiologicallyBasedDemographicModels.jl",
    remotes = nothing,
    format = Documenter.LaTeX(platform = "none"),
    build = "build_pdf",
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Getting Started" => "tutorials/01_getting_started.md",
            "Cotton Plant Model" => "tutorials/02_cotton_plant.md",
            "Coffee Berry Borer" => "tutorials/03_coffee_berry_borer.md",
            "Grapevine C/N Model" => "tutorials/04_grapevine.md",
            "Lobesia Overwintering" => "tutorials/05_lobesia_overwintering.md",
            "Bt Cotton Resistance" => "tutorials/06_bt_cotton_resistance.md",
            "Pesticide Resistance Optimization" => "tutorials/07_pesticide_resistance.md",
            "Screwworm SIT Eradication" => "tutorials/08_screwworm_sit.md",
            "Indian Bt Cotton Bioeconomics" => "tutorials/09_bt_cotton_india.md",
            "Tsetse Ecosocial Model" => "tutorials/10_tsetse_ecosocial.md",
            "Olive Climate Economics" => "tutorials/11_olive_climate.md",
        ],
        "API Reference" => [
            "Types & Traits" => "api/types.md",
            "Weather & Forcing" => "api/weather.md",
            "Dynamics" => "api/dynamics.md",
            "Problem & Solution" => "api/problems.md",
            "Analysis" => "api/analysis.md",
            "Multi-Species Interactions" => "api/interactions.md",
            "Economics" => "api/economics.md",
            "Genetics" => "api/genetics.md",
            "Epidemiology" => "api/epidemiology.md",
            "Utilities" => "api/utilities.md",
        ],
    ])

deploydocs(;
    repo = "github.com/ecorecipes/PhysiologicallyBasedDemographicModels.jl.git",
)
