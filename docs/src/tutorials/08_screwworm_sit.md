# Screwworm SIT Eradication
Simon Frost

- [Background](#background)
- [Screwworm Biology](#screwworm-biology)
- [Temperature-Dependent Mortality](#temperature-dependent-mortality)
- [Cumulative Cold Index](#cumulative-cold-index)
- [The Sterile Insect Technique](#the-sterile-insect-technique)
- [Population Dynamics with SIT](#population-dynamics-with-sit)
- [Myiasis Prediction Model](#myiasis-prediction-model)
- [Historical Eradication Timeline](#historical-eradication-timeline)
- [The Critical Period: 1976–1982](#the-critical-period-19761982)
- [Sterile Male Competitiveness](#sterile-male-competitiveness)
- [Climate Warming Projections](#climate-warming-projections)
- [Key Insights](#key-insights)

Primary reference: (Gutierrez et al. 2019).

## Background

The New World screwworm (*Cochliomyia hominivorax*) was one of the most
devastating livestock pests in the Americas. Females lay eggs in open
wounds of warm-blooded animals; larvae feed on living tissue, causing
myiasis that is often fatal. The pest was eradicated from North America
using the **Sterile Insect Technique (SIT)** — mass release of
irradiation-sterilized males to suppress wild populations through
reproductive interference.

This is a genetic control method: because female screwworms mate only
once (monandry), a female that mates with a sterile male produces no
offspring. When sterile males vastly outnumber wild males, the effective
reproduction rate drops below replacement and the population collapses.

**Reference:** Gutierrez, A.P., Ponti, L., and Arias, P.A. (2019).
*Deconstructing the eradication of new world screwworm in North America:
retrospective analysis and climate warming effects.* Medical and
Veterinary Entomology 33:282–295.

## Screwworm Biology

``` julia
using PhysiologicallyBasedDemographicModels

# --- Temperature thresholds ---
const SW_T_LOWER = 14.5    # °C lower developmental threshold
const SW_T_UPPER = 43.5    # °C upper developmental threshold
const SW_T_REF   = 27.2    # °C reference (optimal) temperature

# Development rate
sw_dev = LinearDevelopmentRate(SW_T_LOWER, SW_T_UPPER)

# --- Life history parameters ---
const FEMALE_HALF_LIFE = 3.7     # days (mated, field; Thomas & Chen 1990)
const MAX_FECUNDITY = 67.0       # eggs/female/day
const SEX_RATIO = 0.5            # 1:1 sex ratio
const DOUBLING_TIME = 14.0       # days (potential, ideal conditions)
const FIELD_DOUBLING = 54.0      # days (observed, tropical endemic)

println("Screwworm life history:")
println("  Temperature range: $(SW_T_LOWER)–$(SW_T_UPPER)°C")
println("  Optimal temperature: $(SW_T_REF)°C")
println("  Female half-life: $(FEMALE_HALF_LIFE) days")
println("  Max fecundity: $(MAX_FECUNDITY) eggs/female/day")
println("  Potential doubling time: $(DOUBLING_TIME) days")
```

    Screwworm life history:
      Temperature range: 14.5–43.5°C
      Optimal temperature: 27.2°C
      Female half-life: 3.7 days
      Max fecundity: 67.0 eggs/female/day
      Potential doubling time: 14.0 days

## Temperature-Dependent Mortality

The daily adult mortality rate follows a quadratic function centered on
the optimal temperature. Cold winters drive cumulative mortality that
determines whether an area can sustain endemic populations.

``` julia
# Daily mortality rate (Eq. 1 of Gutierrez et al. 2019)
# µ_ab(T) = 0.00036(T - 27.2)² + 0.0035
# r² = 0.74, d.f. = 16
function daily_mortality(T::Float64)
    return clamp(0.00036 * (T - SW_T_REF)^2 + 0.0035, 0.0, 1.0)
end

println("\nDaily adult mortality at different temperatures:")
println("T (°C) | Mortality/day | Daily survival | Annual survival")
println("-"^65)
for T in [5.0, 10.0, 15.0, 20.0, 27.2, 30.0, 35.0, 40.0]
    mu = daily_mortality(T)
    surv_daily = 1 - mu
    surv_annual = surv_daily^365
    println("  $(rpad(T, 5)) | $(round(mu, digits=5))      | " *
            "$(round(surv_daily, digits=4))        | $(round(surv_annual, digits=6))")
end
```


    Daily adult mortality at different temperatures:
    T (°C) | Mortality/day | Daily survival | Annual survival
    -----------------------------------------------------------------
      5.0   | 0.18092      | 0.8191        | 0.0
      10.0  | 0.11      | 0.89        | 0.0
      15.0  | 0.05708      | 0.9429        | 0.0
      20.0  | 0.02216      | 0.9778        | 0.00028
      27.2  | 0.0035      | 0.9965        | 0.278109
      30.0  | 0.00632      | 0.9937        | 0.098766
      35.0  | 0.0254      | 0.9746        | 8.3e-5
      40.0  | 0.06248      | 0.9375        | 0.0

## Cumulative Cold Index

The cumulative cold mortality index (µ_cold) determines whether a
location can sustain year-round populations. Areas with µ_cold ≤ 10 are
considered endemic zones.

``` julia
# Cumulative cold mortality (Eq. 2)
# µ_cold(y) = Σ µ_ab(T) for T < 27.2°C (Sept 1 – May 31)
function cumulative_cold(daily_temps::Vector{Float64}; start_day=244, end_day=151)
    mu_cold = 0.0
    for (i, T) in enumerate(daily_temps)
        doy = ((i - 1 + start_day - 1) % 365) + 1
        # Only count Sept 1 (244) through May 31 (151 next year)
        if doy >= start_day || doy <= end_day
            if T < SW_T_REF
                mu_cold += daily_mortality(T)
            end
        end
    end
    return mu_cold
end

# Compare locations along the eradication transect
println("\nCumulative cold index by location (approximate):")

locations = [
    ("Uvalde, TX",        28.0, 10.0),   # Cold winters
    ("McAllen, TX",       26.0, 13.0),   # Transition zone
    ("Tampico, Mexico",   22.0, 16.0),   # Tropical with cool winters
    ("Tuxtla-Gutiérrez",  16.7, 20.0),   # Tropical, mild winters
]

for (name, lat, mean_T) in locations
    # Generate approximate annual temperatures
    amplitude = 12.0 - 0.15 * (30 - lat)  # Higher lat = more seasonal
    temps = [mean_T + amplitude * sin(2π * (d - 200) / 365) for d in 1:365]

    mu_cold = cumulative_cold(temps)
    status = mu_cold > 10 ? "non-endemic" : "ENDEMIC"

    println("  $(rpad(name, 22)) µ_cold = $(round(mu_cold, digits=1))  → $status")
end
```


    Cumulative cold index by location (approximate):
      Uvalde, TX             µ_cold = 46.5  → non-endemic
      McAllen, TX            µ_cold = 34.9  → non-endemic
      Tampico, Mexico        µ_cold = 24.7  → non-endemic
      Tuxtla-Gutiérrez       µ_cold = 14.1  → non-endemic

## The Sterile Insect Technique

SIT exploits female monandry: a female that mates with a sterile male
has no viable offspring. The key parameter is the **overflooding ratio**
— sterile males released per wild male.

``` julia
# --- SIT parameters ---
const MATING_RATE = 0.5        # Proportion of virgin females mating per day
const H_SITES = 100.0          # Oviposition site density (constant)

# Reproduction equation (Eq. 4 of Gutierrez et al.)
# ΔE = φ_T × φ_search × sr × R × (W_m + 0.5 × W_u)
# With SIT (Eq. 4i):
# ΔE = φ_T × φ_search × R × (W_m + 0.5 × φ_comp × φ_release × W_u)

"""
Effective reproduction rate under SIT.
- wild_males: number of wild males
- sterile_males: number of released sterile males
- wild_females_unmated: unmated wild females
- competitiveness: mating competitiveness of sterile males (0-1)
Returns: fraction of matings that are fertile
"""
function sit_fertile_fraction(wild_males, sterile_males;
                              competitiveness=1.0)
    effective_sterile = sterile_males * competitiveness
    total = wild_males + effective_sterile
    total ≈ 0.0 && return 1.0
    return wild_males / total
end

# Demonstrate overflooding effects
println("\nSIT overflooding ratio effects:")
println("Ratio (S:W) | Fertile matings | Population growth")
println("-"^55)
for ratio in [0, 1, 2, 5, 10, 20, 50, 100]
    fertile = sit_fertile_fraction(100, 100 * ratio)
    # Net growth: R₀ with SIT = R₀_base × fertile fraction
    R0_base = 3.0  # Approximate net reproductive rate
    R0_sit = R0_base * fertile
    status = R0_sit < 1.0 ? "DECLINING" : "growing"
    println("  $(rpad("$ratio:1", 10)) |    $(round(fertile*100, digits=1))%       | " *
            "R₀=$(round(R0_sit, digits=2)) ($status)")
end
```


    SIT overflooding ratio effects:
    Ratio (S:W) | Fertile matings | Population growth
    -------------------------------------------------------
      0:1        |    100.0%       | R₀=3.0 (growing)
      1:1        |    50.0%       | R₀=1.5 (growing)
      2:1        |    33.3%       | R₀=1.0 (growing)
      5:1        |    16.7%       | R₀=0.5 (DECLINING)
      10:1       |    9.1%       | R₀=0.27 (DECLINING)
      20:1       |    4.8%       | R₀=0.14 (DECLINING)
      50:1       |    2.0%       | R₀=0.06 (DECLINING)
      100:1      |    1.0%       | R₀=0.03 (DECLINING)

## Population Dynamics with SIT

``` julia
"""
Simulate screwworm population with SIT releases.
Uses distributed delay for life stages.
"""
function simulate_screwworm_sit(;
    n_days=365,
    daily_temps=nothing,
    sterile_release_rate=0.0,  # sterile males per day
    release_interval=14,       # days between releases
    competitiveness=1.0,
    initial_pop=100.0)

    if daily_temps === nothing
        # Default: McAllen, TX approximate climate
        daily_temps = [26.0 + 10.0 * sin(2π * (d - 200) / 365) for d in 1:n_days]
    end

    # State: wild females (mated, unmated), wild males, sterile males
    wild_F_mated = initial_pop * 0.3
    wild_F_unmated = initial_pop * 0.2
    wild_M = initial_pop * 0.5
    sterile_M = 0.0

    # Track population
    pop_history = Float64[]
    fertile_eggs_history = Float64[]

    for day in 1:n_days
        T = daily_temps[day]

        # Temperature-dependent reproduction scaling
        phi_T = T >= SW_T_LOWER && T <= SW_T_UPPER ?
                1.0 - ((T - SW_T_REF) / (SW_T_UPPER - SW_T_REF))^2 : 0.0
        phi_T = clamp(phi_T, 0.0, 1.0)

        # SIT release (biweekly)
        if sterile_release_rate > 0 && day % release_interval == 0
            sterile_M += sterile_release_rate
        end

        # Sterile male decay (shorter lifespan than wild)
        sterile_M *= 0.90  # 10% daily loss

        # Mating: frequency-dependent
        fertile_frac = sit_fertile_fraction(wild_M, sterile_M;
                                            competitiveness=competitiveness)

        # New matings from unmated females
        new_matings = wild_F_unmated * MATING_RATE
        new_fertile = new_matings * fertile_frac
        new_sterile_mated = new_matings * (1 - fertile_frac)

        wild_F_mated += new_fertile
        wild_F_unmated -= new_matings
        # Sterile-mated females produce no offspring (removed from breeding)

        # Reproduction
        fertile_eggs = phi_T * MAX_FECUNDITY * SEX_RATIO * wild_F_mated
        fertile_eggs = max(0.0, fertile_eggs)

        # Mortality
        mu = daily_mortality(T)
        wild_F_mated *= (1 - mu)
        wild_F_unmated *= (1 - mu)
        wild_M *= (1 - mu)

        # New adults emerging (simplified: eggs from ~14 days ago)
        if day > 14
            emergence = fertile_eggs_history[max(1, day-14)] * 0.05
            wild_F_unmated += emergence * SEX_RATIO
            wild_M += emergence * (1 - SEX_RATIO)
        end

        total_wild = wild_F_mated + wild_F_unmated + wild_M
        push!(pop_history, total_wild)
        push!(fertile_eggs_history, fertile_eggs)
    end

    return pop_history
end

# Compare: no SIT vs various release rates
println("\nPopulation after 365 days with SIT:")
println("Release rate | Final pop | Reduction")
println("-"^50)

baseline = simulate_screwworm_sit(sterile_release_rate=0.0)
base_final = baseline[end]

for rate in [0, 50, 100, 500, 1000, 5000]
    pop = simulate_screwworm_sit(sterile_release_rate=Float64(rate))
    final = pop[end]
    reduction = base_final > 0 ? (1 - final/base_final) * 100 : 0.0
    println("  $(rpad(rate, 11))| $(round(final, digits=1))    | $(round(reduction, digits=1))%")
end
```


    Population after 365 days with SIT:
    Release rate | Final pop | Reduction
    --------------------------------------------------
      0          | 2.559438153518334e18    | 0.0%
      50         | 2.2955667118016013e18    | 10.3%
      100        | 2.0973810562851374e18    | 18.1%
      500        | 1.2787068481343227e18    | 50.0%
      1000       | 8.332897346215962e17    | 67.4%
      5000       | 1.0120996803879067e17    | 96.0%

## Myiasis Prediction Model

The paper developed a regression model linking weather to myiasis
outbreak severity.

``` julia
# Myiasis prediction (Eq. 3)
# log₁₀(myiasis) = 6.568 + 0.028×rain(y-1) - 0.602×µ_cold(y-1)
# r² = 0.63, F = 12.63, d.f. = 15

function predict_myiasis(rain_prev_year::Float64, mu_cold_prev::Float64)
    log_myiasis = 6.568 + 0.028 * rain_prev_year - 0.602 * mu_cold_prev
    return 10^log_myiasis
end

println("\nMyiasis outbreak prediction:")
println("Rain (mm) | µ_cold | Predicted cases")
println("-"^45)
for (rain, mu) in [(30, 12.0), (50, 10.0), (50, 8.0),
                    (70, 8.0), (80, 6.0), (100, 5.0)]
    cases = predict_myiasis(Float64(rain), Float64(mu))
    println("  $(rpad(rain, 8)) | $(rpad(mu, 5)) | $(round(Int, cases))")
end
```


    Myiasis outbreak prediction:
    Rain (mm) | µ_cold | Predicted cases
    ---------------------------------------------
      30       | 12.0  | 2
      50       | 10.0  | 89
      50       | 8.0   | 1419
      70       | 8.0   | 5152
      80       | 6.0   | 157036
      100      | 5.0   | 2280342

## Historical Eradication Timeline

``` julia
println("\nHistorical eradication of New World screwworm:")
println("="^65)

timeline = [
    (1957, "SIT eradication begins",      "Florida"),
    (1962, "Major Texas campaign starts",  "51,600 cases"),
    (1972, "Largest outbreak during SIT",  "95,600 cases"),
    (1976, "Critical cold period begins",  "Wild pop. suppressed"),
    (1979, "Decisive year",                "Cold + dry + high S:W ratio"),
    (1982, "Last U.S. autochthonous case", "Eradication achieved"),
    (1990, "Mexico eradication progresses","Southward advance"),
    (2000, "Panama border reached",        "Darién Gap containment"),
    (2016, "Florida Keys incident",        "188M sterile flies released"),
]

for (year, event, detail) in timeline
    println("  $year: $event ($detail)")
end
```


    Historical eradication of New World screwworm:
    =================================================================
      1957: SIT eradication begins (Florida)
      1962: Major Texas campaign starts (51,600 cases)
      1972: Largest outbreak during SIT (95,600 cases)
      1976: Critical cold period begins (Wild pop. suppressed)
      1979: Decisive year (Cold + dry + high S:W ratio)
      1982: Last U.S. autochthonous case (Eradication achieved)
      1990: Mexico eradication progresses (Southward advance)
      2000: Panama border reached (Darién Gap containment)
      2016: Florida Keys incident (188M sterile flies released)

## The Critical Period: 1976–1982

The eradication succeeded because of a **synergy between cold weather
and SIT**, not SIT alone. Despite average SIT efficacy of only ~1.7%,
the combined pressure drove populations below the Allee extinction
threshold.

``` julia
# Simulate the critical 1976-1982 period
println("\nSimulating 1976-1982 critical period:")
println("Year | µ_cold | SIT ratio | Wild pop (relative)")
println("-"^55)

# Cold years suppressed wild populations
cold_years = Dict(
    1976 => 10.3, 1977 => 9.8, 1978 => 9.1,
    1979 => 10.1, 1980 => 9.5, 1981 => 8.5, 1982 => 9.0
)

wild_pop = 1000.0  # Relative starting population
sterile_pop = 5000.0  # Constant release

for year in 1976:1982
    mu_cold = cold_years[year]

    # Winter mortality from cold
    winter_survival = exp(-mu_cold * 0.1)

    # SIT effect
    fertile = sit_fertile_fraction(wild_pop, sterile_pop)

    # Combined suppression
    wild_pop *= winter_survival * fertile * 2.5  # R₀ ≈ 2.5

    # Allee effect: below threshold, finding mates becomes difficult
    if wild_pop < 10.0
        wild_pop *= 0.5  # Additional Allee effect penalty
    end

    ratio = sterile_pop / max(1.0, wild_pop)
    println("  $year | $(round(mu_cold, digits=1))   | $(round(ratio, digits=0)):1     | " *
            "$(round(wild_pop, digits=1))")
end

if wild_pop < 1.0
    println("\n→ Population driven to extinction by 1982!")
else
    println("\n→ Population reduced to $(round(wild_pop, digits=1))")
end
```


    Simulating 1976-1982 critical period:
    Year | µ_cold | SIT ratio | Wild pop (relative)
    -------------------------------------------------------
      1976 | 10.3   | 34.0:1     | 148.8
      1977 | 9.8   | 2480.0:1     | 2.0
      1978 | 9.1   | 5000.0:1     | 0.0
      1979 | 10.1   | 5000.0:1     | 0.0
      1980 | 9.5   | 5000.0:1     | 0.0
      1981 | 8.5   | 5000.0:1     | 0.0
      1982 | 9.0   | 5000.0:1     | 0.0

    → Population driven to extinction by 1982!

## Sterile Male Competitiveness

A critical uncertainty: lab-reared, irradiated males may be less
competitive than wild males. This dramatically affects the release ratio
needed.

``` julia
println("\nEffect of sterile male competitiveness:")
println("Competitiveness | Required S:W ratio for R₀ < 1")
println("-"^55)

R0_wild = 3.0  # Wild population net reproductive rate

for comp in [1.0, 0.8, 0.5, 0.3, 0.1, 0.05]
    # Find ratio where R₀ × fertile_fraction < 1
    for ratio in 1:1000
        fertile = sit_fertile_fraction(100, 100 * ratio; competitiveness=comp)
        if R0_wild * fertile < 1.0
            println("  $(rpad(comp, 14)) | $(ratio):1")
            break
        end
    end
end

println("\nField observations: average SIT efficacy was only ~1.7%")
println("This implies competitiveness × mating success was very low")
println("Success required weather suppression to reduce wild populations first")
```


    Effect of sterile male competitiveness:
    Competitiveness | Required S:W ratio for R₀ < 1
    -------------------------------------------------------
      1.0            | 3:1
      0.8            | 3:1
      0.5            | 5:1
      0.3            | 7:1
      0.1            | 21:1
      0.05           | 41:1

    Field observations: average SIT efficacy was only ~1.7%
    This implies competitiveness × mating success was very low
    Success required weather suppression to reduce wild populations first

## Climate Warming Projections

``` julia
println("\nClimate warming effects on screwworm range:")
println("="^55)

# Under RCP 8.5, ~+2°C by 2050
for warming in [0.0, 1.0, 2.0, 3.0]
    for (name, lat, base_T) in [("McAllen, TX", 26.0, 26.0),
                                 ("San Antonio, TX", 29.5, 22.0),
                                 ("Dallas, TX", 32.8, 19.0)]
        temps = [(base_T + warming) + 10.0 * sin(2π * (d - 200) / 365)
                 for d in 1:365]
        mu = cumulative_cold(temps)
        status = mu <= 10.0 ? "ENDEMIC" : "seasonal"

        if warming == 0.0
            print("  $(rpad(name, 18)) baseline µ=$(round(mu, digits=1))")
        else
            print("  $(rpad("", 18)) +$(warming)°C: µ=$(round(mu, digits=1))")
        end
        println("  ($status)")
    end
    println()
end

println("Risk: warming expands endemic range northward into southern U.S.")
println("Current containment at Panama Darién Gap costs ~\$15M/year")
```


    Climate warming effects on screwworm range:
    =======================================================
      McAllen, TX        baseline µ=5.1  (ENDEMIC)
      San Antonio, TX    baseline µ=10.4  (seasonal)
      Dallas, TX         baseline µ=16.1  (seasonal)

                         +1.0°C: µ=4.1  (ENDEMIC)
                         +1.0°C: µ=8.8  (ENDEMIC)
                         +1.0°C: µ=14.0  (seasonal)

                         +2.0°C: µ=3.3  (ENDEMIC)
                         +2.0°C: µ=7.4  (ENDEMIC)
                         +2.0°C: µ=12.1  (seasonal)

                         +3.0°C: µ=2.5  (ENDEMIC)
                         +3.0°C: µ=6.2  (ENDEMIC)
                         +3.0°C: µ=10.4  (seasonal)

    Risk: warming expands endemic range northward into southern U.S.
    Current containment at Panama Darién Gap costs ~$15M/year

## Key Insights

1.  **Monandry is the key**: SIT only works because female screwworms
    mate once. If females could re-mate, sterile releases would need to
    be vastly larger to reduce fertile mating probability.

2.  **SIT alone was insufficient**: Average field efficacy was only
    ~1.7%. Eradication required the coincidence of cold winters
    (1976–1980) suppressing wild populations, making the sterile:wild
    ratio effective enough to drive populations below the Allee
    threshold.

3.  **The Allee effect was decisive**: Below a critical density,
    screwworms cannot find mates or hosts efficiently. SIT pushes
    populations into this extinction vortex — a demographic, not
    genetic, mechanism.

4.  **Overflooding was massive**: Field releases of 39–896 sterile
    males/km²/week, when models suggest only 2–23/km²/biweekly were
    needed. This 10–225× excess compensated for low competitiveness and
    field mortality of sterile males.

5.  **Climate warming threatens re-invasion**: Each degree of warming
    shifts the endemic zone northward. The current containment barrier
    at the Panama–Colombia border requires continuous SIT releases to
    prevent re-establishment.

6.  **Weather prediction enables management**: The myiasis regression
    model (r² = 0.63) shows that previous-year rainfall and cold index
    predict current outbreak severity, enabling proactive release
    planning.

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Gutierrez2019Screwworm" class="csl-entry">

Gutierrez, Andrew Paul, Luigi Ponti, and Pablo A. Arias. 2019.
“Deconstructing the Eradication of New World Screwworm in North America:
Retrospective Analysis and Climate Warming Effects.” *Medical and
Veterinary Entomology* 33: 282–95. <https://doi.org/10.1111/mve.12362>.

</div>

</div>
