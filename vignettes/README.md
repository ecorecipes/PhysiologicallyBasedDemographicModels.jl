# Vignettes — `PhysiologicallyBasedDemographicModels.jl`

This directory contains 64 worked vignettes (Quarto `.qmd`).
Most reproduce a published PBDM study end-to-end; vignettes 50–56 are
API/feature documentation, and #49 is an original speculative scenario.

Each vignette is rendered in three formats. The committed GFM markdown is
linked below; HTML and PDF are gitignored but reproducible via
`bash scripts/render_all_vignettes.sh {gfm|html|pdf}`.

## Fidelity grades

Grades come from the comprehensive audit (see `CHANGELOG.md` and the
`vign_audit_*.md` reports in the workspace root). **Bold** grades indicate a
post-audit fix; the parenthetical shows the pre-fix grade.

- **A**: faithful reproduction — parameters and functional forms verified.
- **B**: mostly faithful, minor pedagogical simplifications documented.
- **C**: significant simplifications or transcription issues (now mostly fixed).
- **N**: API documentation or original scenario, not a paper reproduction.

## Vignette index

| # | Title | Source | Grade |
|---|---|---|---|
| 01 | [Getting Started with PBDMs](01_getting_started/01_getting_started.md) | `@Manetsch1976` | A |
| 02 | [Cotton Plant Model](02_cotton_plant/02_cotton_plant.md) | `@Gutierrez1984Cotton` | A |
| 03 | [Coffee Berry Borer Lifecycle](03_coffee_berry_borer/03_coffee_berry_borer.md) | `@Cure2020Coffee` | A |
| 04 | [Grapevine Carbon-Nitrogen Model](04_grapevine/04_grapevine.md) | `@Wermelinger1991Grapevine` | A |
| 05 | [Lobesia Overwintering and Diapause](05_lobesia_overwintering/05_lobesia_overwintering.md) | `@Baumgartner2012Lobesia` | B |
| 06 | [Bt Cotton Resistance Evolution](06_bt_cotton_resistance/06_bt_cotton_resistance.md) | `@Gutierrez2006BtCottonChina` | B |
| 07 | [Pesticide Resistance Optimization](07_pesticide_resistance/07_pesticide_resistance.md) | `@Gutierrez1979PesticideResistance` | A |
| 08 | [Screwworm SIT Eradication](08_screwworm_sit/08_screwworm_sit.md) | `@Gutierrez2019Screwworm` | A |
| 09 | [Indian Bt Cotton Bioeconomics](09_bt_cotton_india/09_bt_cotton_india.md) | `@Gutierrez2020IndianBtCotton` | B |
| 10 | [Tsetse-Cattle-Human Ecosocial Model](10_tsetse_ecosocial/10_tsetse_ecosocial.md) | `@Baumgartner2008EcoSocial` | B |
| 11 | [Mediterranean Olive Climate Economics](11_olive_climate/11_olive_climate.md) | `@Ponti2014OliveClimate` | B+ |
| 12 | [Cassava Mealybug Biocontrol in Africa](12_cassava_mealybug/12_cassava_mealybug.md) | `@Gutierrez1988CassavaMealybug` | A− |
| 13 | [Mediterranean Fruit Fly Invasive Potential](13_medfly_invasion/13_medfly_invasion.md) | `@Gutierrez2011Medfly` | A |
| 14 | [East Coast Fever in African Livestock](14_east_coast_fever/14_east_coast_fever.md) | `@Gilioli2009EastCoastFever` | **A** (was C+) |
| 15 | [Rice–Weed Competition](15_rice_weed_competition/15_rice_weed_competition.md) | `@Graf1990RiceWeeds` | A |
| 16 | [Aedes albopictus Invasion Risk in Europe](16_aedes_albopictus/16_aedes_albopictus.md) | `@Pasquali2020Aedes` | A |
| 17 | [Pink Bollworm Climate Limits](17_pink_bollworm/17_pink_bollworm.md) | `@Gutierrez2006PinkBollworm` | A |
| 18 | [Vine Mealybug Biocontrol](18_vine_mealybug/18_vine_mealybug.md) | `@Gutierrez2008VineMealybug` | A |
| 19 | [Brown Marmorated Stink Bug Biocontrol](19_bmsb_biocontrol/19_bmsb_biocontrol.md) | `@Gutierrez2023BMSB` | A |
| 20 | [Cabbage Root Fly Diapause Dynamics](20_cabbage_maggot/20_cabbage_maggot.md) | `@Johnsen1997Cabbage` | A |
| 21 | [Processing Tomato Crop-Pest Management](21_tomato_ipm/21_tomato_ipm.md) | `@Wilson1986TOMDAT` | B |
| 22 | [Tuta absoluta Invasion Risk Assessment](22_tuta_absoluta/22_tuta_absoluta.md) | `@Ponti2021Tuta` | A |
| 23 | [Yellow Starthistle Biological Control](23_yellow_starthistle/23_yellow_starthistle.md) | `@Gutierrez2005YellowStarthistle` | A |
| 24 | [Asian Citrus Psyllid and Citrus Greening Disease](24_asian_citrus_psyllid/24_asian_citrus_psyllid.md) | `@Gutierrez2013ACP` | B+ |
| 25 | [Bemisia tabaci Invasion Risk in Europe](25_bemisia_tabaci/25_bemisia_tabaci.md) | `@Gilioli2014Bemisia` | A |
| 26 | [Light Brown Apple Moth Invasion Potential](26_light_brown_apple_moth/26_light_brown_apple_moth.md) | — | B |
| 27 | [Spotted Alfalfa Aphid Biological Control](27_spotted_alfalfa_aphid/27_spotted_alfalfa_aphid.md) | — | A |
| 28 | [Olive–Olive Fly System Under Climate Warming](28_olive_bactrocera/28_olive_bactrocera.md) | `@Ponti2009OliveFly` | A |
| 29 | [Cowpea–Thrips Agroecosystem in West Africa](29_cowpea_thrips/29_cowpea_thrips.md) | `@Tamo1993Cowpea` | **B** (was C) |
| 30 | [Common Bean Growth Types I–III: Yield Prediction](30_bean_growth/30_bean_growth.md) | — | **B** (was C) |
| 31 | [Cotton–Boll Weevil Interaction in Brazil](31_cotton_boll_weevil/31_cotton_boll_weevil.md) | — | B |
| 32 | [Oleander Scale Regulation by Competing Parasitoids](32_oleander_scale/32_oleander_scale.md) | `@Gutierrez2007OleanderScale` | B |
| 33 | [Rice–Fish Integrated Agroecosystem](33_rice_fish_agroecosystem/33_rice_fish_agroecosystem.md) | `@DOultremont2002RiceFish` | B |
| 34 | [Plant–Aphid–Parasitoid Tritrophic Dynamics](34_plant_aphid_parasite/34_plant_aphid_parasite.md) | `@Gilbert1973PlantAphidParasite` | **A** (was B) |
| 35 | [Golden Delicious Apple Tree Growth Model](35_apple_tree/35_apple_tree.md) | `@Baumgartner1986Apple` | B |
| 36 | [Fusarium Wilt and Root-Knot Nematode in Cotton](36_fusarium_nematode/36_fusarium_nematode.md) | — | B |
| 37 | [Cotton-Whitefly-Autoparasitoid Dynamics](37_whitefly_autoparasitoid/37_whitefly_autoparasitoid.md) | — | A |
| 38 | [Tropical Fruit Fly Invasive Potential Under Climate Change](38_tropical_fruit_flies/38_tropical_fruit_flies.md) | — | B |
| 39 | [Economics of Bt Cotton in China](39_china_bt_cotton/39_china_bt_cotton.md) | — | B |
| 40 | [Fall Armyworm Establishment Risk in Europe](40_spodoptera_frugiperda/40_spodoptera_frugiperda.md) | `@Gilioli2023Spodoptera` | A |
| 41 | [Olive Fruit Fly Population Dynamics](41_olive_fly_ode/41_olive_fly_ode.md) | `@Rossini2022BactroceraODE` | A |
| 42 | [Biological Control of Spodoptera exigua](42_spodoptera_biocontrol/42_spodoptera_biocontrol.md) | `@GarayNarvaez2015Spodoptera` | B |
| 43 | [Philaenus spumarius Phenology Model](43_philaenus_phenology/43_philaenus_phenology.md) | `@Sweidan2025PhilaenusPhenology` | A |
| 44 | [Physiological Risk Index for Pest Establishment](44_risk_index/44_risk_index.md) | `@Ndjomatchoua2024RiskIndex` | A |
| 45 | [Drosophila suzukii Adult Male Dynamics](45_dsuzukii_pde/45_dsuzukii_pde.md) | `@Rossini2020DSuzukiiMales` | A |
| 46 | [Mediterranean Fruit Fly Under Climate Change](46_medfly_kolmogorov/46_medfly_kolmogorov.md) | `@Gilioli2022CeratitisNonlinear` | B |
| 47 | [Stage-Structured Consumer-Resource Dynamics](47_consumer_resource/47_consumer_resource.md) | `@DeRoos2008PSPM` | A |
| 48 | [Xylella fastidiosa Transmission in Olive Groves](48_xylella_ecoepi/48_xylella_ecoepi.md) | `@Weber2026Xylella` | B |
| 49 | [Genetically Modified Thermal Tolerance in Bombus terrestris](49_bombus_hsp/49_bombus_hsp.md) | — | N |
| 50 | [Type Hierarchy Reference](50_type_hierarchy/50_type_hierarchy.md) | — | N |
| 51 | [State Variables, Bulk Populations, and Phase Callbacks](51_state_variables/51_state_variables.md) | — | N |
| 52 | [Extended Rules and Scheduled Events](52_extended_rules_events/52_extended_rules_events.md) | — | N |
| 53 | [Theoretical Helpers — Compensation, Isoclines, and Assembly](53_theory_helpers/53_theory_helpers.md) | — | N |
| 54 | [Continuous-Time and PSPM Formulations](54_continuous_pspm/54_continuous_pspm.md) | — | N |
| 55 | [Management Optimisation and Economics](55_management_economics/55_management_economics.md) | — | N |
| 56 | [Ensembles, Filters, and Weather Helpers](56_ensembles_misc/56_ensembles_misc.md) | — | N |
| 57 | [Verticillium Wilt — Multi-Season Dynamic-Programming Management](57_verticillium_dp/57_verticillium_dp.md) | — | **A** (was C) |
| 58 | [Coffee Berry Borer — Bio-Economic Analysis of Control Tactics](58_cbb_bioeconomics/58_cbb_bioeconomics.md) | — | **A** (was B) |
| 59 | [Olive / Olive-Fly Bio-Economics under Climate Warming](59_olive_climate/59_olive_climate.md) | — | B |
| 60 | [Invasion-Risk Assessment for *Tuta absoluta* — a Mechanistic PBDM](60_tuta_absoluta_invasion/60_tuta_absoluta_invasion.md) | — | B |
| 61 | [Voltinism Shifts of *Lobesia botrana* under Climate Warming](61_lobesia_voltinism/61_lobesia_voltinism.md) | — | **A** (was B) |
| 62 | [SIT Overflooding Ratio for Screwworm Eradication](62_screwworm_sit/62_screwworm_sit.md) | — | B |
| 63 | [Bayesian Inference of Stage-Specific Mortality](63_bayesian_mortality/63_bayesian_mortality.md) | — | **B** (was C) |
| 64 | [Tritrophic Biocontrol Design for Brown Marmorated Stink Bug](64_bmsb_tritrophic/64_bmsb_tritrophic.md) | — | **B** (was C) |

## Building

Render one vignette:
```bash
cd vignettes/01_getting_started
quarto render 01_getting_started.qmd --to gfm   # or html, pdf
```

Render all 64 in batch:
```bash
bash scripts/render_all_vignettes.sh gfm    # or html, pdf
```
