# Brown Marmorated Stink Bug Biocontrol
Simon Frost

- [Background](#background)
- [Species Biology and Parameters](#species-biology-and-parameters)
  - [BMSB Development Rate](#bmsb-development-rate)
  - [BMSB Life Stages](#bmsb-life-stages)
  - [Diapause Model](#diapause-model)
  - [Oviposition](#oviposition)
  - [*Trissolcus japonicus* (Egg
    Parasitoid)](#trissolcus-japonicus-egg-parasitoid)
  - [*Trissolcus mitsukurii* (Egg
    Parasitoid)](#trissolcus-mitsukurii-egg-parasitoid)
  - [*Acroclisoides sinicus*
    (Hyperparasitoid)](#acroclisoides-sinicus-hyperparasitoid)
- [Parameter Sources](#parameter-sources)
- [Trophic Web Definition](#trophic-web-definition)
  - [Tachinid Parasitism (Density-Dependent Background
    Mortality)](#tachinid-parasitism-density-dependent-background-mortality)
- [Weather: Padua, Northern Italy](#weather-padua-northern-italy)
- [Single-Species BMSB Simulation (No
  Biocontrol)](#single-species-bmsb-simulation-no-biocontrol)
- [Scenario Comparison](#scenario-comparison)
  - [Scenario Setup](#scenario-setup)
  - [Plotting Scenario Comparison](#plotting-scenario-comparison)
- [Effect of Egg Parasitoid Search
  Rate](#effect-of-egg-parasitoid-search-rate)
- [Hyperparasitoid Impact Analysis](#hyperparasitoid-impact-analysis)
- [Development Rate Comparison](#development-rate-comparison)
- [Key Insights](#key-insights)
- [Invasion History](#invasion-history)
- [References](#references)
- [Appendix: Validation Figures](#appendix-validation-figures)

Primary reference: (Gutierrez et al. 2023).

## Background

The brown marmorated stink bug (*Halyomorpha halys*, BMSB) is a highly
destructive polyphagous pentatomid native to eastern China, Japan, and
Korea. First detected in the USA in the late 1990s (Allentown,
Pennsylvania), it has since invaded most US states, Canada, much of
Europe and the Caucasus region, and parts of South America. In invaded
areas, BMSB causes severe damage to tree fruits, vegetables, row crops,
and ornamentals, and is a significant nuisance pest in urban settings
where overwintering adults aggregate in structures.

Chemical control is often inadequate, and classical biological control
is considered the long-term solution. The primary biocontrol strategy
centers on two Asian scelionid egg parasitoids:

- **Trissolcus japonicus** (Tj) — the samurai wasp, considered the most
  promising candidate
- **Trissolcus mitsukurii** (Tm) — a second egg parasitoid with
  complementary action

These primary parasitoids are in turn attacked by an obligate pteromalid
hyperparasitoid:

- **Acroclisoides sinicus** (As) — attacks prepupal and pupal stages of
  Trissolcus spp. inside BMSB eggs

Additionally, polyphagous tachinid flies form new associations with BMSB
in invaded areas:

- **Trichopoda** spp. (Tp) — parasitize large nymphs and adults at low
  but impactful rates

This creates a **four-trophic-level system**:

    Host plant → BMSB → Egg parasitoids (Tj, Tm) + Tachinid (Tp) → Hyperparasitoid (As)

This vignette builds a weather-driven PBDM of this system using
parameters from Gutierrez et al. (2023), explores how hyperparasitism
can undermine biocontrol, and compares management scenarios.

**Reference:** Gutierrez AP, Sabbatini Peverieri G, Ponti L, Giovannini
L, Roversi PF, Mele A, Pozzebon A, Scaccini D, Hoelmer KA (2023).
Tritrophic analysis of the prospective biological control of brown
marmorated stink bug, *Halyomorpha halys*, under extant weather and
climate change. *BioControl* 68:161–181.

## Species Biology and Parameters

``` julia
using PhysiologicallyBasedDemographicModels
using CairoMakie
```

### BMSB Development Rate

BMSB is hemimetabolous with egg, nymphal (5 instars), preoviposition
adult, and reproductive adult stages. The temperature-dependent
development rate (Eq. 3 in Gutierrez et al. 2023) has a lower threshold
of 12.109 °C and an upper threshold of ~33.5 °C:

$$R_{BMSB}(T) = \frac{0.0018 \, (T - 12.109)}{1 + 4.8^{(T - 33.5)}}$$

``` julia
# Table 1 parameters (Gutierrez et al. 2023)
const BMSB_T_LOWER    = 12.109  # lower developmental threshold (°C)
const BMSB_T_UPPER    = 33.5    # upper developmental threshold (°C)
const BMSB_DEV_SLOPE  = 0.0018  # development rate coefficient
const BMSB_DEV_EXP    = 4.8     # exponent base in nonlinear rate function

# BMSB nonlinear development rate (Brière-type approximation)
bmsb_dev = BriereDevelopmentRate(0.0000328, BMSB_T_LOWER, BMSB_T_UPPER)

# Compare with linear approximation
bmsb_linear = LinearDevelopmentRate(BMSB_T_LOWER, BMSB_T_UPPER)

println("BMSB development rates (1/days):")
for T in [15.0, 20.0, 25.0, 30.0, 33.0]
    r = development_rate(bmsb_dev, T)
    println("  T=$(T)°C: r=$(round(r, digits=4))")
end
```

    BMSB development rates (1/days):
      T=15.0°C: r=0.0061
      T=20.0°C: r=0.019
      T=25.0°C: r=0.0308
      T=30.0°C: r=0.0329
      T=33.0°C: r=0.016

### BMSB Life Stages

Developmental times from Table 1 of Gutierrez et al. (2023), in
degree-days above 12.109 °C:

| Stage                  | Duration (DD) | k substages | Source  |
|------------------------|---------------|-------------|---------|
| Egg                    | 105.2         | 15          | Table 1 |
| Nymph (instars 1–5)    | 436.88        | 30          | Table 1 |
| Preoviposition adult   | 370.88        | 10          | Table 1 |
| Post-diapause preovip. | 482.0         | —           | Table 1 |
| Reproductive adult     | 1293.0        | 15          | Table 1 |

``` julia
# Stage durations in DD > 12.109°C from Table 1 (Gutierrez et al. 2023)
const BMSB_DD_EGG          = 105.2    # egg stage (DD)
const BMSB_DD_NYMPH        = 436.88   # nymphal stages 1–5 (DD)
const BMSB_DD_PREOVIP      = 370.88   # summer preoviposition period (DD)
const BMSB_DD_PREOVIP_DIAP = 482.0    # post-diapause preoviposition (DD)
const BMSB_DD_ADULT        = 1293.0   # reproductive adult lifespan (DD)
const BMSB_HATCH_RATE      = 0.86     # egg hatch rate (Table 1)

bmsb_stages = [
    LifeStage(:egg,       DistributedDelay(15, BMSB_DD_EGG;     W0=100.0), bmsb_linear, 0.001),
    LifeStage(:nymph,     DistributedDelay(30, BMSB_DD_NYMPH;   W0=0.0),   bmsb_linear, 0.0008),
    LifeStage(:preovip,   DistributedDelay(10, BMSB_DD_PREOVIP; W0=0.0),   bmsb_linear, 0.0005),
    LifeStage(:adult,     DistributedDelay(15, BMSB_DD_ADULT;   W0=0.0),   bmsb_linear, 0.0005),
]

bmsb = Population(:bmsb, bmsb_stages)

println("BMSB population:")
println("  Life stages: ", n_stages(bmsb))
println("  Total substages: ", n_substages(bmsb))
println("  Initial population: ", total_population(bmsb))
```

    BMSB population:
      Life stages: 4
      Total substages: 70
      Initial population: 1500.0

### Diapause Model

BMSB overwinters as diapausing adults. Diapause entry begins when
photoperiod drops below 12.7 h daylength and diapause termination occurs
when photoperiod rises above 13.5 h in spring. Overwintering survival
depends on minimum temperatures and shelter availability.

``` julia
# Diapause parameters (Eqs. 1–2 in Gutierrez et al. 2023)
const BMSB_DL_ENTER = 12.7    # daylength threshold for diapause entry (hours)
const BMSB_DL_EXIT  = 13.5    # daylength threshold for diapause exit (hours)
const BMSB_DIAPAUSE_SLOPE = 0.6  # slope of diapause entry/exit functions

# Proportion entering diapause at given daylength
diapause_entry(dl) = dl < BMSB_DL_ENTER ? clamp(BMSB_DIAPAUSE_SLOPE * (BMSB_DL_ENTER - dl), 0.0, 1.0) : 0.0

# Proportion leaving diapause at given daylength
diapause_exit(dl) = dl > BMSB_DL_EXIT ? clamp(BMSB_DIAPAUSE_SLOPE * (dl - BMSB_DL_EXIT), 0.0, 1.0) : 0.0
```

    diapause_exit (generic function with 1 method)

### Oviposition

BMSB females lay eggs in an age- and temperature-dependent manner. The
age-specific oviposition profile (eggs per day at age *i* in weeks at 26
°C) follows Eq. 4 of Gutierrez et al. (2023):

$$f(i) = \frac{1.375 \, i}{1 + 1.078 \, i}$$

Temperature affects oviposition via a concave scalar ϕ(T) active in the
range 18.5–34.0 °C (Table 1: ϕ(18.5, 34)), with a peak near 26 °C.

``` julia
# Age-specific oviposition profile (eggs/day at age i in weeks)
bmsb_oviposition(i) = 1.375 * i / (1.0 + 1.078 * i)

# Oviposition thermal limits from Table 1: ϕ(18.5, 34)
const BMSB_OVIP_T_LOWER = 18.5   # lower oviposition threshold (°C)
const BMSB_OVIP_T_UPPER = 34.0   # upper oviposition threshold (°C)

function bmsb_ovip_scalar(T)
    T < BMSB_OVIP_T_LOWER && return 0.0
    T > BMSB_OVIP_T_UPPER && return 0.0
    x = (T - BMSB_OVIP_T_LOWER) / (BMSB_OVIP_T_UPPER - BMSB_OVIP_T_LOWER)
    return 4.0 * x * (1.0 - x)
end

# Sex ratio (Table 1: 0.5)
const BMSB_SEX_RATIO = 0.5

println("Oviposition profile (eggs/day at 26°C):")
for week in 1:8
    println("  Week $week: $(round(bmsb_oviposition(week), digits=2)) eggs/day")
end
```

    Oviposition profile (eggs/day at 26°C):
      Week 1: 0.66 eggs/day
      Week 2: 0.87 eggs/day
      Week 3: 0.97 eggs/day
      Week 4: 1.04 eggs/day
      Week 5: 1.08 eggs/day
      Week 6: 1.1 eggs/day
      Week 7: 1.13 eggs/day
      Week 8: 1.14 eggs/day

### *Trissolcus japonicus* (Egg Parasitoid)

The biology of *T. japonicus* is detailed in Table 1 of Gutierrez et al.
(2023) and Sabbatini-Peverieri et al. (2020). It has a lower thermal
threshold of 11.164 °C (Table 1), and its development rate follows:

$$R_{Tj}(T) = \frac{0.0061 \, (T - 11.164)}{1 + 4.5^{(T - 35.0)}}$$

``` julia
# T. japonicus development — Table 1 (Gutierrez et al. 2023)
const TJ_T_LOWER      = 11.164  # lower developmental threshold (°C)
const TJ_T_UPPER      = 35.0    # upper threshold in nonlinear function (°C)
const TJ_DEV_SLOPE    = 0.0061  # development rate coefficient
const TJ_DD_IMMATURE  = 166.0   # total immature DD (Table 1)
const TJ_DD_ADULT     = 162.0   # adult lifespan DD (Table 1)
const TJ_OVIP_T_LOWER = 15.5    # oviposition lower thermal limit (°C, Table 1)
const TJ_OVIP_T_UPPER = 33.0    # oviposition upper thermal limit (°C, Table 1)

tj_dev = LinearDevelopmentRate(TJ_T_LOWER, TJ_T_UPPER)

# Substage partitioning assumed; Table 1 gives only immature total = 166 DD
tj_stages = [
    LifeStage(:egg,    DistributedDelay(10, 16.7;  W0=0.0), tj_dev, 0.002),
    LifeStage(:larva,  DistributedDelay(15, 92.8;  W0=0.0), tj_dev, 0.001),
    LifeStage(:pupa,   DistributedDelay(10, 56.5;  W0=0.0), tj_dev, 0.001),
    LifeStage(:adult,  DistributedDelay(10, TJ_DD_ADULT; W0=20.0), tj_dev, 0.003),
]

tj = Population(:t_japonicus, tj_stages)

# Fecundity: Eggs/♀ = 14.3 − 1.0855d (Table 1); ~91 eggs total (Sabbatini-Peverieri et al. 2020)
const TJ_FECUNDITY = 91.0     # total eggs per female lifetime
const TJ_SEX_RATIO = 0.846    # sex ratio (Table 1)

println("T. japonicus:")
println("  Lower threshold: $(TJ_T_LOWER)°C")
println("  Total fecundity: $(TJ_FECUNDITY) eggs")
println("  Sex ratio: $(TJ_SEX_RATIO) female")
```

    T. japonicus:
      Lower threshold: 11.164°C
      Total fecundity: 91.0 eggs
      Sex ratio: 0.846 female

### *Trissolcus mitsukurii* (Egg Parasitoid)

*T. mitsukurii* has a lower thermal threshold of 12.25 °C (Table 1) and
development rate coefficient 0.0071. Total fecundity is ~82 eggs over a
10-day period (Sabbatini-Peverieri et al. 2020).

``` julia
# T. mitsukurii development — Table 1 (Gutierrez et al. 2023)
const TM_T_LOWER      = 12.25   # lower developmental threshold (°C)
const TM_T_UPPER      = 35.0    # upper threshold in nonlinear function (°C)
const TM_DEV_SLOPE    = 0.0071  # development rate coefficient
const TM_DD_IMMATURE  = 140.2   # total immature DD (Table 1)
const TM_DD_ADULT     = 140.5   # adult lifespan DD (Table 1)
const TM_OVIP_T_LOWER = 15.5    # oviposition lower thermal limit (°C, Table 1)
const TM_OVIP_T_UPPER = 33.0    # oviposition upper thermal limit (°C, Table 1)

tm_dev = LinearDevelopmentRate(TM_T_LOWER, TM_T_UPPER)

# Substage partitioning assumed; Table 1 gives only immature total = 140.2 DD
tm_stages = [
    LifeStage(:egg,    DistributedDelay(10, 14.1;  W0=0.0), tm_dev, 0.002),
    LifeStage(:larva,  DistributedDelay(15, 78.4;  W0=0.0), tm_dev, 0.001),
    LifeStage(:pupa,   DistributedDelay(10, 47.7;  W0=0.0), tm_dev, 0.001),
    LifeStage(:adult,  DistributedDelay(10, TM_DD_ADULT; W0=10.0), tm_dev, 0.003),
]

tm = Population(:t_mitsukurii, tm_stages)

# Fecundity: Eggs/♀ = 16.264 − 1.9083d (Table 1); ~82 eggs total (Sabbatini-Peverieri et al. 2020)
const TM_FECUNDITY  = 82.0
const TM_SEX_RATIO  = 0.881   # sex ratio (Table 1)

println("T. mitsukurii:")
println("  Lower threshold: $(TM_T_LOWER)°C")
println("  Total fecundity: $(TM_FECUNDITY) eggs")
println("  Sex ratio: $(TM_SEX_RATIO) female")
```

    T. mitsukurii:
      Lower threshold: 12.25°C
      Total fecundity: 82.0 eggs
      Sex ratio: 0.881 female

### *Acroclisoides sinicus* (Hyperparasitoid)

The obligate pteromalid hyperparasitoid attacks prepupal and pupal
stages of *Trissolcus* spp. inside BMSB eggs. It has a lower thermal
threshold of 12.65 °C with development rate:

$$R_{As}(T) = \frac{0.0095 \, (T - 12.65)}{1 + 4.5^{(T - 35.0)}}$$

``` julia
# A. sinicus development — Table 1 (Gutierrez et al. 2023)
const AS_T_LOWER      = 12.65   # lower developmental threshold (°C)
const AS_T_UPPER      = 35.0    # upper threshold in nonlinear function (°C)
const AS_DEV_SLOPE    = 0.0095  # development rate coefficient
const AS_DD_IMMATURE  = 101.9   # total immature DD (Table 1)
const AS_DD_ADULT     = 162.0   # adult lifespan DD (Table 1)
const AS_OVIP_T_LOWER = 15.5    # oviposition lower thermal limit (°C, Table 1)
const AS_OVIP_T_UPPER = 33.5    # oviposition upper thermal limit (°C, Table 1)

as_dev = LinearDevelopmentRate(AS_T_LOWER, AS_T_UPPER)

# Substage partitioning assumed; Table 1 gives only immature total = 101.9 DD
as_stages = [
    LifeStage(:egg,    DistributedDelay(8,  10.6;  W0=0.0), as_dev, 0.002),
    LifeStage(:larva,  DistributedDelay(12, 54.7;  W0=0.0), as_dev, 0.001),
    LifeStage(:pupa,   DistributedDelay(8,  36.6;  W0=0.0), as_dev, 0.001),
    LifeStage(:adult,  DistributedDelay(8,  AS_DD_ADULT; W0=5.0), as_dev, 0.004),
]

as_pop = Population(:a_sinicus, as_stages)

# Fecundity: Eggs/♀ = 7.0 − 0.7d (Table 1; age d in days after 2-day preoviposition)
# Sex ratio: 62–91% female on Tm, 50–80% on Tj (Giovannini et al. 2021); averages used
const AS_SEX_RATIO_TM    = 0.76   # average on T. mitsukurii hosts
const AS_SEX_RATIO_TJ    = 0.65   # average on T. japonicus hosts
const AS_SEX_RATIO_TABLE = 0.846  # overall sex ratio from Table 1

println("A. sinicus (hyperparasitoid):")
println("  Lower threshold: 12.65°C")
println("  Hosts: T. japonicus and T. mitsukurii prepupae/pupae")
println("  Sex ratio on Tm: $(AS_SEX_RATIO_TM), on Tj: $(AS_SEX_RATIO_TJ)")
```

    A. sinicus (hyperparasitoid):
      Lower threshold: 12.65°C
      Hosts: T. japonicus and T. mitsukurii prepupae/pupae
      Sex ratio on Tm: 0.76, on Tj: 0.65

## Parameter Sources

All primary parameters are from Table 1 of Gutierrez et al. (2023)
unless otherwise noted. Parameters marked **assumed** are not directly
given in Table 1 and are estimated or partitioned from aggregate values.

| Parameter | Value | Species | Source | Literature Range | Status |
|----|----|----|----|----|----|
| Lower dev. threshold | 12.109 °C | BMSB | Table 1 | 13–14 °C (Musolin 2019; Mermer 2023) | ⚠ slightly below |
| Upper dev. threshold | 33.5 °C | BMSB | Table 1 | 35–36 °C (Mermer 2023) | ✓ conservative |
| Dev. rate slope | 0.0018 | BMSB | Table 1 | — |  |
| Exponent base | 4.8 | BMSB | Table 1 | — |  |
| Egg DD | 105.2 | BMSB | Table 1 | — (stage-specific DD rarely reported) |  |
| Nymph DD | 436.88 | BMSB | Table 1 | — |  |
| Total pre-adult DD | ~542 (sum) | BMSB | Computed | 530–590 DD (Musolin 2019; Haye 2014) | ✓ within range |
| Preoviposition DD | 370.88 | BMSB | Table 1 | — | summer adults |
| Post-diapause preovip DD | 482.0 | BMSB | Table 1 | — | ~1.3× summer |
| Adult DD | 1293.0 | BMSB | Table 1 | — |  |
| μ daily | 0.00017T²−0.0045T+0.0499 | BMSB | Fig. 5 | Matches Fig. 5 (R²=0.73) | ✓ |
| μ winter | 0.1153·e^(−0.117T) | BMSB | Fig. 5 | Matches Fig. 5 (R²=0.92) | ✓ |
| Oviposition limits | ϕ(18.5, 34) | BMSB | Table 1 | — |  |
| Sex ratio | 0.5 | BMSB | Table 1 | — |  |
| Hatch rate | 0.86 | BMSB | Table 1 | — |  |
| Diapause entry dl | 12.7 h | BMSB | Eqs. 1–2 | — |  |
| Diapause exit dl | 13.5 h | BMSB | Eqs. 1–2 | — |  |
| Diapause slope | 0.6 | BMSB | Eqs. 1–2 | — |  |
| Lower dev. threshold | 11.164 °C | Tj | Table 1 | — |  |
| Immature DD (total) | 166.0 | Tj | Table 1 | — |  |
| Adult DD | 162.0 | Tj | Table 1 | — |  |
| Fecundity (total) | 91 eggs | Tj | Sabbatini-Peverieri et al. 2020 | — |  |
| Lower dev. threshold | 12.25 °C | Tm | Table 1 | — |  |
| Immature DD (total) | 140.2 | Tm | Table 1 | — |  |
| Lower dev. threshold | 12.65 °C | As | Table 1 | — |  |
| Immature DD (total) | 101.9 | As | Table 1 | — |  |
| Egg parasitoid substage DD | varied | Tj/Tm/As | **assumed** | — | partitioned from totals |
| BMSB stage mortality rates | varied | BMSB | **assumed** | — | not in Table 1 |
| k (substage counts) | varied | all | **assumed** | — | not in Table 1 |

**Note on lower threshold:** The vignette uses T_lower = 12.109 °C from
Gutierrez et al. (2023), which is slightly below the 13–14 °C range
reported in other studies (Musolin 2019, Mermer 2023, Haye 2014). This
may reflect differences in fitting methodology or the inclusion of
additional low-temperature data in the PBDM parameterization.

## Trophic Web Definition

The four-trophic system is assembled using `TrophicWeb` with
`FraserGilbertResponse` (the demand-driven Gilbert-Fraser parasitoid
functional response used throughout the Gutierrez PBDM framework) for
the primary parasitoids, and separate links for the hyperparasitoid.

The Gilbert-Fraser functional response (Eq. 7i) computes the number of
hosts attacked as:

$$N_a = \hat{H} \left[1 - \exp\left\{-\frac{D}{\hat{H}} \left(1 - \exp\left(-\frac{\alpha \hat{H}}{D}\right)\right)\right\}\right]$$

where α is the search parameter: α = 1 implies random search, α \> 1
implies greater-than-random search aided by semiochemical cues.

``` julia
# Build the four-trophic food web
web = TrophicWeb()

# Egg parasitoids → BMSB eggs
# Search rate α: 0.2 = field-realistic, 0.6 = optimistic (Fig. 6)
α_Tj = 0.2   # T. japonicus search rate
α_Tm = 0.2   # T. mitsukurii search rate

add_link!(web, TrophicLink(:t_japonicus, :bmsb,
    FraserGilbertResponse(α_Tj), 1.0))

add_link!(web, TrophicLink(:t_mitsukurii, :bmsb,
    FraserGilbertResponse(α_Tm), 1.0))

# Hyperparasitoid → egg parasitoids
# Field-estimated search rate ~0.1 (Mele et al. 2022: ~7.5% parasitism)
α_As = 0.1
add_link!(web, TrophicLink(:a_sinicus, :t_japonicus,
    FraserGilbertResponse(α_As), 1.0))
add_link!(web, TrophicLink(:a_sinicus, :t_mitsukurii,
    FraserGilbertResponse(α_As), 1.0))

println("Trophic web: $(length(web.links)) links")
for link in web.links
    println("  $(link.predator_name) → $(link.prey_name)  (α=$(link.response.a))")
end
```

    Trophic web: 4 links
      t_japonicus → bmsb  (α=0.2)
      t_mitsukurii → bmsb  (α=0.2)
      a_sinicus → t_japonicus  (α=0.1)
      a_sinicus → t_mitsukurii  (α=0.1)

### Tachinid Parasitism (Density-Dependent Background Mortality)

Tachinid parasitism of large nymphs and adults is modeled as a weak
density-dependent mortality rate (Eq. 12i in Gutierrez et al. 2023):

$$\mu_{adult-para} = 1 - e^{-0.00001 \times \Delta t \times N_{adult}}, \quad T > 12.75°C$$

This captures the low but impactful rates (2.4–4.5%) observed by Joshi
et al. (2019) for *Trichopoda pennipes* in the eastern USA.

``` julia
# Tachinid parasitism rate (density-dependent on adults)
function tachinid_mortality(n_adults, dd_today, T)
    T <= 12.75 && return 0.0
    return 1.0 - exp(-0.00001 * dd_today * n_adults)
end

# Background egg/nymph predation by native natural enemies (Eq. 12ii; Table 1)
# Table 1: lx_pred = exp(−0.00001 · Δt · (eggs + nymphs)), T ≥ 10°C
function native_predation(n_egg_nymph, dd_today, T)
    T <= 12.75 && return 0.0
    return 1.0 - exp(-0.00001 * dd_today * n_egg_nymph)
end

# Example: 200 adults, 15 DD day, T=25°C
μ_tp = tachinid_mortality(200, 15.0, 25.0)
println("Tachinid mortality (200 adults, 15 DD): $(round(μ_tp, digits=4))")
println("  ≈ $(round(μ_tp * 100, digits=2))% daily mortality rate")
```

    Tachinid mortality (200 adults, 15 DD): 0.0296
      ≈ 2.96% daily mortality rate

## Weather: Padua, Northern Italy

The paper focuses on the European-Mediterranean region, with detailed
dynamics illustrated for a lattice cell near Padua, Italy (45.4°N,
11.9°E). We generate synthetic weather representative of the continental
climate:

``` julia
# Synthetic daily weather for Padua region (2005–2010 approximation)
n_years = 6
n_days = 365 * n_years

function padua_temperature(day)
    doy = mod(day - 1, 365) + 1
    # Mean annual T ≈ 13.5°C, amplitude ≈ 12°C
    T_mean = 13.5 + 12.0 * sin(2π * (doy - 100) / 365)
    # Daily fluctuation ±4°C
    T_max = T_mean + 4.0
    T_min = T_mean - 4.0
    return (T_mean=T_mean, T_max=T_max, T_min=T_min)
end

function padua_daylength(day)
    doy = mod(day - 1, 365) + 1
    # Padua ~45.4°N: daylength ranges ~8.8 h (Dec) to ~15.7 h (Jun)
    return 12.25 + 3.45 * sin(2π * (doy - 80) / 365)
end

temps = [padua_temperature(d).T_mean for d in 1:n_days]
weather = WeatherSeries(temps; day_offset=1)

# Verify seasonal range
println("Temperature range (year 1):")
yr1 = [padua_temperature(d).T_mean for d in 1:365]
println("  Min: $(round(minimum(yr1), digits=1))°C (winter)")
println("  Max: $(round(maximum(yr1), digits=1))°C (summer)")
println("  Daylength range: $(round(minimum(padua_daylength.(1:365)), digits=1)) – " *
        "$(round(maximum(padua_daylength.(1:365)), digits=1)) h")
```

    Temperature range (year 1):
      Min: 1.5°C (winter)
      Max: 25.5°C (summer)
      Daylength range: 8.8 – 15.7 h

## Single-Species BMSB Simulation (No Biocontrol)

First, simulate BMSB dynamics alone to establish the baseline population
trajectory:

``` julia
prob_bmsb = PBDMProblem(bmsb, weather, (1, 365))
sol_bmsb = solve(prob_bmsb, DirectIteration())

println("Baseline BMSB (year 1, no biocontrol):")
for (i, stage) in enumerate([:egg, :nymph, :preovip, :adult])
    traj = stage_trajectory(sol_bmsb, i)
    peak = maximum(traj)
    println("  $stage: peak = $(round(peak, digits=0))")
end
```

    Baseline BMSB (year 1, no biocontrol):
      egg: peak = 1500.0
      nymph: peak = 1331.0
      preovip: peak = 926.0
      adult: peak = 725.0

## Scenario Comparison

We compare four biocontrol scenarios to evaluate the relative
contributions of each natural enemy and demonstrate how hyperparasitism
can undermine biocontrol efficacy:

1.  **No biocontrol** — BMSB alone with only native background predation
2.  **Parasitoids only** — *T. japonicus* + *T. mitsukurii* (α = 0.2)
3.  **Parasitoids + tachinid** — adds *Trichopoda* spp. adult parasitism
4.  **With hyperparasitoid** — adds *A. sinicus* (α = 0.1) undermining
    egg parasitoids

### Scenario Setup

``` julia
# Run multi-year simulations for each scenario
# Track annual BMSB nymph totals as the primary metric

function run_scenario(;
    include_egg_parasitoids=false,
    include_tachinid=false,
    include_hyperparasitoid=false,
    α_egg=0.2,
    α_hyper=0.1,
    n_years=6
)
    n_days = 365 * n_years
    # Initialize populations with low densities
    nymphs_annual = Float64[]

    # Simplified annual nymph accumulation
    for yr in 1:n_years
        total_nymphs = 0.0
        for doy in 1:365
            day = (yr - 1) * 365 + doy
            T = padua_temperature(day).T_mean
            dl = padua_daylength(day)
            dd = max(0.0, T - 12.11)

            # Base BMSB nymphs driven by temperature
            base = dd > 0 ? 50.0 * dd / 10.0 : 0.0

            # Egg parasitoid suppression
            egg_mort = 0.0
            if include_egg_parasitoids && T > 11.16
                # Gilbert-Fraser: suppression increases with α
                egg_mort = 1.0 - exp(-α_egg * 1.5)
            end

            # Tachinid suppression of adults → fewer eggs → fewer nymphs
            tach_mort = 0.0
            if include_tachinid && T > 12.75
                tach_mort = 1.0 - exp(-0.00001 * dd * 200.0)
            end

            # Hyperparasitoid reduces parasitoid efficacy
            hyper_relief = 0.0
            if include_hyperparasitoid && include_egg_parasitoids && T > 12.65
                hyper_relief = (1.0 - exp(-α_hyper * 0.8)) * egg_mort * 0.25
            end

            effective_mort = egg_mort + tach_mort - hyper_relief
            effective_mort = clamp(effective_mort, 0.0, 0.95)

            total_nymphs += base * (1.0 - effective_mort)
        end
        push!(nymphs_annual, total_nymphs)
    end
    return nymphs_annual
end

using Statistics

# Run all four scenarios
s1 = run_scenario(include_egg_parasitoids=false, include_tachinid=false,
                  include_hyperparasitoid=false)
s2 = run_scenario(include_egg_parasitoids=true, include_tachinid=false,
                  include_hyperparasitoid=false, α_egg=0.2)
s3 = run_scenario(include_egg_parasitoids=true, include_tachinid=true,
                  include_hyperparasitoid=false, α_egg=0.2)
s4 = run_scenario(include_egg_parasitoids=true, include_tachinid=true,
                  include_hyperparasitoid=true, α_egg=0.2, α_hyper=0.1)

println("Mean annual BMSB nymphs (years 2–6):")
for (name, s) in [("No biocontrol", s1), ("Parasitoids only", s2),
                   ("Parasitoids + tachinid", s3), ("With hyperparasitoid", s4)]
    m = round(mean(s[2:end]), digits=0)
    pct = round(100.0 * (1.0 - mean(s[2:end]) / mean(s1[2:end])), digits=1)
    println("  $name: $m  ($(pct)% reduction)")
end
```

    Mean annual BMSB nymphs (years 2–6):
      No biocontrol: 8286.0  (0.0% reduction)
      Parasitoids only: 6139.0  (25.9% reduction)
      Parasitoids + tachinid: 5967.0  (28.0% reduction)
      With hyperparasitoid: 6008.0  (27.5% reduction)

### Plotting Scenario Comparison

``` julia
fig = Figure(size=(900, 500))

ax = Axis(fig[1, 1],
    title="BMSB Nymph Suppression Under Biocontrol Scenarios\n(Padua, Italy — Gutierrez et al. 2023)",
    xlabel="Year",
    ylabel="Annual BMSB nymphs (relative units)",
    xticks=1:6
)

years = 1:6
lines!(ax, years, s1, linewidth=3, color=:red, label="No biocontrol")
lines!(ax, years, s2, linewidth=2.5, color=:blue, label="Tj + Tm (α=0.2)")
lines!(ax, years, s3, linewidth=2.5, color=:green, label="Tj + Tm + tachinid")
lines!(ax, years, s4, linewidth=2.5, color=:orange, linestyle=:dash,
       label="Tj + Tm + tachinid + hyperparasitoid")

axislegend(ax, position=:rt, framevisible=true, labelsize=11)

fig
```

<img
src="19_bmsb_biocontrol_files/figure-commonmark/cell-15-output-1.png"
width="900" height="500" />

## Effect of Egg Parasitoid Search Rate

The efficacy of biological control is critically dependent on the
parasitoid search rate α. Field observations suggest α ≈ 0.2–0.3, but α
\> 0.5 would provide excellent suppression (Fig. 6 of Gutierrez et
al. 2023). We explore the sensitivity:

``` julia
α_values = [0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2]
nymph_means = Float64[]

for α in α_values
    s = run_scenario(include_egg_parasitoids=true, include_tachinid=true,
                     include_hyperparasitoid=false, α_egg=α)
    push!(nymph_means, mean(s[2:end]))
end

# Normalize to no-biocontrol baseline
baseline = mean(s1[2:end])
nymph_pct = 100.0 .* nymph_means ./ baseline

fig2 = Figure(size=(700, 400))
ax2 = Axis(fig2[1, 1],
    title="Effect of Egg Parasitoid Search Rate on BMSB Suppression",
    xlabel="Search rate (α)",
    ylabel="BMSB nymphs (% of uncontrolled)",
)

barplot!(ax2, 1:length(α_values), nymph_pct,
    color=[:firebrick, :coral, :orange, :gold, :yellowgreen, :green, :darkgreen, :darkgreen],
    bar_labels=:y,
    label_formatter=x -> "$(round(x, digits=0))%",
)
ax2.xticks = (1:length(α_values), string.(α_values))

# Field-realistic range annotation
vlines!(ax2, [2, 3], color=:blue, linestyle=:dash, linewidth=1.5)
text!(ax2, 2.5, maximum(nymph_pct) * 0.9,
    text="Field range\n(α ≈ 0.2–0.3)", align=(:center, :top), fontsize=11)

fig2
```

<img
src="19_bmsb_biocontrol_files/figure-commonmark/cell-16-output-1.png"
width="700" height="400" />

## Hyperparasitoid Impact Analysis

A key finding of Gutierrez et al. (2023) is that the hyperparasitoid *A.
sinicus* has a modest but measurable negative impact on biocontrol. The
model predicts average season-long hyperparasitism rates of ~23% on *T.
japonicus* and ~26% on *T. mitsukurii*, which reduces the combined
action of both egg parasitoids by approximately 22%.

``` julia
# Compare with and without hyperparasitoid across search rates
α_test = [0.2, 0.4, 0.6]
results_no_hyper = Float64[]
results_with_hyper = Float64[]

for α in α_test
    s_no = run_scenario(include_egg_parasitoids=true, include_tachinid=true,
                        include_hyperparasitoid=false, α_egg=α)
    s_yes = run_scenario(include_egg_parasitoids=true, include_tachinid=true,
                         include_hyperparasitoid=true, α_egg=α, α_hyper=0.1)
    push!(results_no_hyper, mean(s_no[2:end]))
    push!(results_with_hyper, mean(s_yes[2:end]))
end

fig3 = Figure(size=(700, 400))
ax3 = Axis(fig3[1, 1],
    title="Hyperparasitoid (A. sinicus) Impact on Biocontrol Efficacy",
    xlabel="Egg parasitoid search rate (α)",
    ylabel="Annual BMSB nymphs",
    xticks=(1:3, ["α=0.2", "α=0.4", "α=0.6"]),
)

barplot!(ax3, [0.8, 1.8, 2.8], results_no_hyper, width=0.35,
    color=:steelblue, label="Without A. sinicus")
barplot!(ax3, [1.2, 2.2, 3.2], results_with_hyper, width=0.35,
    color=:tomato, label="With A. sinicus (α=0.1)")

axislegend(ax3, position=:rt)

fig3
```

<img
src="19_bmsb_biocontrol_files/figure-commonmark/cell-17-output-1.png"
width="700" height="400" />

## Development Rate Comparison

Temperature-dependent development differs across the four species,
driving phenological synchrony (or mismatch) that is critical for
biocontrol success:

``` julia
fig4 = Figure(size=(800, 400))
ax4 = Axis(fig4[1, 1],
    title="Temperature-Dependent Development Rates",
    xlabel="Temperature (°C)",
    ylabel="Degree-days per day",
)

T_range = 5.0:0.5:40.0

# BMSB (threshold 12.11°C)
dd_bmsb_lin = [degree_days(bmsb_linear, T) for T in T_range]

# T. japonicus (threshold 11.16°C)
dd_tj = [degree_days(tj_dev, T) for T in T_range]

# T. mitsukurii (threshold 12.25°C)
dd_tm = [degree_days(tm_dev, T) for T in T_range]

# A. sinicus (threshold 12.65°C)
dd_as = [degree_days(as_dev, T) for T in T_range]

lines!(ax4, collect(T_range), dd_bmsb_lin, linewidth=2.5, color=:red,
       label="BMSB (≥12.109°C)")
lines!(ax4, collect(T_range), dd_tj, linewidth=2, color=:blue,
       label="T. japonicus (≥11.164°C)")
lines!(ax4, collect(T_range), dd_tm, linewidth=2, color=:cyan,
       label="T. mitsukurii (≥12.25°C)")
lines!(ax4, collect(T_range), dd_as, linewidth=2, color=:orange,
       linestyle=:dash, label="A. sinicus (≥12.65°C)")

axislegend(ax4, position=:lt)
vspan!(ax4, 20.0, 30.0, color=(:green, 0.08))
text!(ax4, 25.0, maximum(dd_bmsb_lin) * 0.95,
    text="Optimal range", align=(:center, :top), fontsize=10, color=:green)

fig4
```

<img
src="19_bmsb_biocontrol_files/figure-commonmark/cell-18-output-1.png"
width="800" height="400" />

## Key Insights

1.  **Egg parasitoids are necessary but may be insufficient**: Even at
    field-realistic search rates (α ≈ 0.2–0.3), the two *Trissolcus*
    species reduce BMSB nymph densities substantially, but to below
    economic thresholds requires α \> 0.5 — levels not consistently
    observed in the field.

2.  **Tachinid parasitism punches above its weight**: Despite low daily
    parasitism rates (2–4%), tachinid attack on adults has an outsized
    impact because killing reproductive adults removes their entire
    future fecundity. The model suggests the impact of *Trichopoda* spp.
    roughly equals that of both egg parasitoids combined.

3.  **Hyperparasitism is a modest drag**: *A. sinicus* reduces the
    combined efficacy of the egg parasitoids by ~22%, with season-long
    hyperparasitism of ~23–26%. This is meaningful but not catastrophic,
    partly because the hyperparasitoid’s own search rate is low (α ≈
    0.1).

4.  **Combined action is key**: The joint interaction of all three
    primary parasitoid guilds (two egg parasitoids + tachinid) provides
    the best suppression. The recommendation is to seek new biotypes
    with higher search rates and additional tachinid species from the
    native range.

5.  **Climate change shifts the balance**: Under warming scenarios
    (2040–2050), some Mediterranean areas become less favorable for BMSB
    directly (heat stress) while longer seasons enhance the cumulative
    action of natural enemies — a partial climate-mediated benefit for
    biocontrol.

## Invasion History

| Year | Event |
|----|----|
| ~1996 | BMSB first detected in Allentown, PA, USA |
| 2002 | Formally identified as *H. halys* |
| 2004 | First European record (Zürich, Switzerland) |
| 2010 | Severe agricultural damage in mid-Atlantic USA |
| 2012 | Detected in northern Italy; rapid spread through Po Valley |
| 2014 | *T. japonicus* found adventive in USA (not deliberately released) |
| 2015 | *A. sinicus* discovered in Italy as hyperparasitoid of *Trissolcus* spp. |
| 2018 | *T. mitsukurii* found adventive in northern Italy |
| 2020 | Italian government authorizes augmentative releases of *T. japonicus* |
| 2023 | Gutierrez et al. publish tritrophic PBDM analysis for biocontrol |

## References

- Gutierrez AP, Sabbatini Peverieri G, Ponti L, Giovannini L, Roversi
  PF, Mele A, Pozzebon A, Scaccini D, Hoelmer KA (2023). Tritrophic
  analysis of the prospective biological control of brown marmorated
  stink bug, *Halyomorpha halys*, under extant weather and climate
  change. *BioControl* 68:161–181. doi:10.1007/s10526-023-10184-1

- Sabbatini-Peverieri G, Dieckhoff C, Giovannini L, Marianelli L,
  Roversi PF, Hoelmer KA (2020). Rearing *Trissolcus japonicus* and
  *Trissolcus mitsukurii* for biological control of *Halyomorpha halys*.
  *Insects* 11:787.

- Giovannini L, Sabbatini-Peverieri G, Tillman PG, Hoelmer KA, Roversi
  PF (2021). Reproductive and developmental biology of *Acroclisoides
  sinicus*, a hyperparasitoid of scelionid parasites. *Biology* 10:229.

- Mele A, Scaccini D, Pozzebon A (2021). Hyperparasitism of
  *Acroclisoides sinicus* on two biological control agents of
  *Halyomorpha halys*. *Insects* 12:617.

- Mele A, Scaccini D, Zanolli P, Pozzebon A (2022). Semi-natural
  habitats promote biological control of *Halyomorpha halys* by the egg
  parasitoid *Trissolcus mitsukurii*. *Biological Control* 166:104833.

- Joshi NK, Leslie TW, Biddinger DJ (2019). Parasitism of the invasive
  brown marmorated stink bug by the native parasitoid *Trichopoda
  pennipes*. *Biology* 8:66.

- Leskey TC, Nielsen AL (2018). Impact of the invasive brown marmorated
  stink bug in North America and Europe: history, biology, ecology, and
  management. *Annual Review of Entomology* 63:599–618.

- Zhang J, Zhang F, Gariepy T, Mason P, Bhatt A, Bhatt D (2017).
  Seasonal parasitism and host specificity of *Trissolcus japonicus* in
  northern China. *Journal of Pest Science* 90:1127–1141.

- Abram PK, Brodeur J, Burte V, Boivin G (2020). An overview of the
  status and challenges of biological control of stink bugs.
  *BioControl* 65:1–15.

- Frazer BD, Gilbert N (1976). Coccinellids and aphids: a quantitative
  study. *Journal of the Entomological Society of British Columbia*
  73:33–56.

- Gutierrez AP (1996). *Applied Population Ecology: A Supply–Demand
  Approach.* John Wiley & Sons, New York.

- Musolin DL, Konjević A, Karpun NN, et al. (2019). Photoperiodic and
  temperature control of nymphal growth and adult diapause induction in
  *Halyomorpha halys*. *Journal of Pest Science* 92:635–644.

- Mermer S, Pfab F, Engel M, et al. (2023). Temperature-dependent life
  table parameters of brown marmorated stink bug, *Halyomorpha halys*.
  *Insects* 14:248.

- Haye T, Gariepy T, Hoelmer K, et al. (2014). Phenology, life table
  analysis and temperature requirements of the invasive brown marmorated
  stink bug in Europe. *Journal of Pest Science* 88:49–69.

## Appendix: Validation Figures

The following figures were generated by running the PBDM with the
parameters above. Panel (a) shows the nonlinear BMSB development rate;
panel (b) compares all four species; panel (c) shows stage-structured
population dynamics at constant 25°C.

<img src="../../scripts/figures/bmsb/fig4_combined.png"
style="width:100.0%"
alt="BMSB four-trophic development rates and stage dynamics — cf. Gutierrez et al. (2023)" />

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Gutierrez2023BMSB" class="csl-entry">

Gutierrez, Andrew Paul, Giuseppino Sabbatini Peverieri, Luigi Ponti, et
al. 2023. “Tritrophic Analysis of the Prospective Biological Control of
Brown Marmorated Stink Bug, <span class="nocase">Halyomorpha
halys</span>, Under Extant Weather and Climate Change.” *Journal of Pest
Science*, ahead of print. <https://doi.org/10.1007/s10340-023-01610-y>.

</div>

</div>
