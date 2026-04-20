using Documenter
using PhysiologicallyBasedDemographicModels
using StructuredPopulationCore

# (basename, page-title) for each vignette. Generated from vignette YAML titles.
const TUTORIALS = [
    ("01_getting_started", "Getting Started with PBDMs"),
    ("02_cotton_plant", "Cotton Plant Model"),
    ("03_coffee_berry_borer", "Coffee Berry Borer Lifecycle"),
    ("04_grapevine", "Grapevine Carbon-Nitrogen Model"),
    ("05_lobesia_overwintering", "Lobesia Overwintering and Diapause"),
    ("06_bt_cotton_resistance", "Bt Cotton Resistance Evolution"),
    ("07_pesticide_resistance", "Pesticide Resistance Optimization"),
    ("08_screwworm_sit", "Screwworm SIT Eradication"),
    ("09_bt_cotton_india", "Indian Bt Cotton Bioeconomics"),
    ("10_tsetse_ecosocial", "Tsetse-Cattle-Human Ecosocial Model"),
    ("11_olive_climate", "Mediterranean Olive Climate Economics"),
    ("12_cassava_mealybug", "Cassava Mealybug Biocontrol in Africa"),
    ("13_medfly_invasion", "Mediterranean Fruit Fly Invasive Potential"),
    ("14_east_coast_fever", "East Coast Fever in African Livestock"),
    ("15_rice_weed_competition", "Rice–Weed Competition"),
    ("16_aedes_albopictus", "Aedes albopictus Invasion Risk in Europe"),
    ("17_pink_bollworm", "Pink Bollworm Climate Limits"),
    ("18_vine_mealybug", "Vine Mealybug Biocontrol"),
    ("19_bmsb_biocontrol", "Brown Marmorated Stink Bug Biocontrol"),
    ("20_cabbage_maggot", "Cabbage Root Fly Diapause Dynamics"),
    ("21_tomato_ipm", "Processing Tomato Crop-Pest Management"),
    ("22_tuta_absoluta", "Tuta absoluta Invasion Risk Assessment"),
    ("23_yellow_starthistle", "Yellow Starthistle Biological Control"),
    ("24_asian_citrus_psyllid", "Asian Citrus Psyllid and Citrus Greening Disease"),
    ("25_bemisia_tabaci", "Bemisia tabaci Invasion Risk in Europe"),
    ("26_light_brown_apple_moth", "Light Brown Apple Moth Invasion Potential"),
    ("27_spotted_alfalfa_aphid", "Spotted Alfalfa Aphid Biological Control"),
    ("28_olive_bactrocera", "Olive–Olive Fly System Under Climate Warming"),
    ("29_cowpea_thrips", "Cowpea–Thrips Agroecosystem in West Africa"),
    ("30_bean_growth", "Common Bean Growth Types I–III: Yield Prediction"),
    ("31_cotton_boll_weevil", "Cotton–Boll Weevil Interaction in Brazil"),
    ("32_oleander_scale", "Oleander Scale Regulation by Competing Parasitoids"),
    ("33_rice_fish_agroecosystem", "Rice–Fish Integrated Agroecosystem"),
    ("34_plant_aphid_parasite", "Plant–Aphid–Parasitoid Tritrophic Dynamics"),
    ("35_apple_tree", "Golden Delicious Apple Tree Growth Model"),
    ("36_fusarium_nematode", "Fusarium Wilt and Root-Knot Nematode in Cotton"),
    ("37_whitefly_autoparasitoid", "Cotton-Whitefly-Autoparasitoid Dynamics"),
    ("38_tropical_fruit_flies", "Tropical Fruit Fly Invasive Potential Under Climate Change"),
    ("39_china_bt_cotton", "Economics of Bt Cotton in China"),
    ("40_spodoptera_frugiperda", "Fall Armyworm Establishment Risk in Europe"),
    ("41_olive_fly_ode", "Olive Fruit Fly Population Dynamics"),
    ("42_spodoptera_biocontrol", "Biological Control of Spodoptera exigua"),
    ("43_philaenus_phenology", "Philaenus spumarius Phenology Model"),
    ("44_risk_index", "Physiological Risk Index for Pest Establishment"),
    ("45_dsuzukii_pde", "Drosophila suzukii Adult Male Dynamics"),
    ("46_medfly_kolmogorov", "Mediterranean Fruit Fly Under Climate Change"),
    ("47_consumer_resource", "Stage-Structured Consumer-Resource Dynamics"),
    ("48_xylella_ecoepi", "Xylella fastidiosa Transmission in Olive Groves"),
    ("49_bombus_hsp", "Genetically Modified Thermal Tolerance in Bombus terrestris"),
    ("50_type_hierarchy", "Type Hierarchy Reference"),
    ("51_state_variables", "State Variables, Bulk Populations, and Phase Callbacks"),
    ("52_extended_rules_events", "Extended Rules and Scheduled Events"),
    ("53_theory_helpers", "Theoretical Helpers — Compensation, Isoclines, and Assembly"),
    ("54_continuous_pspm", "Continuous-Time and PSPM Formulations"),
    ("55_management_economics", "Management Optimisation and Economics"),
    ("56_ensembles_misc", "Ensembles, Filters, and Weather Helpers"),
    ("57_verticillium_dp", "Verticillium Wilt — Multi-Season Dynamic-Programming Management"),
    ("58_cbb_bioeconomics", "Coffee Berry Borer — Bio-Economic Analysis of Control Tactics"),
    ("59_olive_climate", "Olive / Olive-Fly Bio-Economics under Climate Warming"),
    ("60_tuta_absoluta_invasion", "Invasion-Risk Assessment for *Tuta absoluta* — a Mechanistic PBDM"),
    ("61_lobesia_voltinism", "Voltinism Shifts of *Lobesia botrana* under Climate Warming"),
    ("62_screwworm_sit", "SIT Overflooding Ratio for Screwworm Eradication"),
    ("63_bayesian_mortality", "Bayesian Inference of Stage-Specific Mortality"),
    ("64_bmsb_tritrophic", "Tritrophic Biocontrol Design for Brown Marmorated Stink Bug"),
]

"""
Sync vignette outputs into `docs/src/tutorials/`.

For each vignette, copy the pre-rendered GFM markdown (`vignettes/<n>/<n>.md`)
and its companion figure directory (`<n>_files/`) into the documentation
source tree. If the markdown does not yet exist, fall back to running
`quarto render` to produce it.
"""
function sync_tutorials_from_vignettes(; force_rerender::Bool = false)
    repo_root = normpath(joinpath(@__DIR__, ".."))
    vignette_dir = joinpath(repo_root, "vignettes")
    tutorial_dir = joinpath(@__DIR__, "src", "tutorials")
    mkpath(tutorial_dir)

    for (base, _title) in TUTORIALS
        srcdir   = joinpath(vignette_dir, base)
        src_md   = joinpath(srcdir, base * ".md")
        src_figs = joinpath(srcdir, base * "_files")
        dst_md   = joinpath(tutorial_dir, base * ".md")
        dst_figs = joinpath(tutorial_dir, base * "_files")

        need_render = force_rerender || !isfile(src_md)
        if need_render
            quarto = Sys.which("quarto")
            quarto === nothing && error("Quarto is required to (re-)render vignette '$base'")
            qmd = joinpath(srcdir, base * ".qmd")
            run(`$quarto render $qmd --to gfm`)
        end

        cp(src_md, dst_md; force = true)
        if isdir(src_figs)
            isdir(dst_figs) && rm(dst_figs; recursive = true)
            cp(src_figs, dst_figs)
        end

        # Quarto's gfm output emits multi-line `<img …>` blocks for figures.
        # Documenter's CommonMark parser escapes these as text, so rewrite
        # them to standard `![](path)` markdown image syntax.
        normalize_img_tags!(dst_md)
    end
end

function normalize_img_tags!(path::AbstractString)
    text = read(path, String)
    rx = r"<img\s+([^>]*?)/?>"s
    fixed = replace(text, rx => function (m)
        attrs = match(rx, m).captures[1]
        srcm = match(r"src=\"([^\"]+)\"", attrs)
        srcm === nothing && return m
        "![]($(srcm.captures[1]))"
    end)
    fixed === text || write(path, fixed)
    return nothing
end

sync_tutorials_from_vignettes(; force_rerender = get(ENV, "PBDM_RERENDER", "false") == "true")

const TUTORIAL_PAGES = [title => "tutorials/$(base).md" for (base, title) in TUTORIALS]

const API_PAGES = [
    "Types & Traits"            => "api/types.md",
    "Weather & Forcing"         => "api/weather.md",
    "Dynamics"                  => "api/dynamics.md",
    "Problem & Solution"        => "api/problems.md",
    "Analysis"                  => "api/analysis.md",
    "Multi-Species Interactions" => "api/interactions.md",
    "Economics"                 => "api/economics.md",
    "Genetics"                  => "api/genetics.md",
    "Epidemiology"              => "api/epidemiology.md",
    "Utilities"                 => "api/utilities.md",
]

makedocs(;
    modules   = [PhysiologicallyBasedDemographicModels, StructuredPopulationCore],
    warnonly  = true,
    authors   = "Simon Frost",
    sitename  = "PhysiologicallyBasedDemographicModels.jl",
    remotes   = nothing,
    format    = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical  = "https://ecorecipes.github.io/PhysiologicallyBasedDemographicModels.jl",
        repolink   = "https://github.com/ecorecipes/PhysiologicallyBasedDemographicModels.jl",
        size_threshold = 1_000_000,
        size_threshold_warn = 500_000,
    ),
    pages = [
        "Home"          => "index.md",
        "Tutorials"     => TUTORIAL_PAGES,
        "API Reference" => API_PAGES,
    ],
)

if get(ENV, "PBDM_BUILD_PDF", "false") == "true"
    makedocs(;
        modules   = [PhysiologicallyBasedDemographicModels, StructuredPopulationCore],
        warnonly  = true,
        authors   = "Simon Frost",
        sitename  = "PhysiologicallyBasedDemographicModels.jl",
        remotes   = nothing,
        format    = Documenter.LaTeX(platform = "none"),
        build     = "build_pdf",
        pages = [
            "Home"          => "index.md",
            "Tutorials"     => TUTORIAL_PAGES,
            "API Reference" => API_PAGES,
        ],
    )
end

deploydocs(;
    repo = "github.com/ecorecipes/PhysiologicallyBasedDemographicModels.jl.git",
)
