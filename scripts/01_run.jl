# Auto-extracted from 01_getting_started.qmd
figdir = joinpath(@__DIR__, "figures", "01_getting_started")
mkpath(figdir)
fig_counter = Ref(0)

using Pkg
Pkg.develop(path="path/to/PhysiologicallyBasedDemographicModels.jl")

using PhysiologicallyBasedDemographicModels

# 1. Define a development rate model
# Eggs develop above 10°C, with upper limit at 35°C
dev_rate = LinearDevelopmentRate(10.0, 35.0)

# 2. Create a life stage with distributed delay
# k=20 substages, mean developmental time τ=100 DD, initial population=1000
egg_delay = DistributedDelay(20, 100.0; W0=1000.0)
egg_stage = LifeStage(:egg, egg_delay, dev_rate, 0.005)  # μ=0.005 per DD

# 3. Create a population with just eggs
pop = Population(:insect, [egg_stage])

println("Substages: ", n_substages(pop))      # 20
println("Initial population: ", total_population(pop))  # 20000 (20 × 1000)
println("Delay variance: ", delay_variance(egg_delay))  # 500.0

# At 25°C: 25 - 10 = 15 degree-days per day
println(degree_days(dev_rate, 25.0))  # 15.0

# At 5°C: below threshold, no development
println(degree_days(dev_rate, 5.0))   # 0.0

# At 40°C: capped at T_upper - T_lower = 25
println(degree_days(dev_rate, 40.0))  # 25.0

# Constant 25°C for 30 days
weather = WeatherSeries(fill(25.0, 30); day_offset=1)

# Create and solve the problem
prob = PBDMProblem(pop, weather, (1, 30))
sol = solve(prob, DirectIteration())

println(sol)  # PBDMSolution(30 days, 1 stages, retcode=Success)

# Cumulative degree-day accumulation
cdd = cumulative_degree_days(sol)
println("Total DD after 30 days: ", round(cdd[end], digits=1))  # 435.0

# Population trajectory over time
tp = total_population(sol)
println("Initial: ", round(tp[1], digits=1))
println("Final:   ", round(tp[end], digits=1))

# Per-step growth rates
println("Mean λ: ", round(mean(sol.lambdas), digits=4))

dev = LinearDevelopmentRate(10.0, 35.0)

stages = [
    LifeStage(:egg,    DistributedDelay(15, 80.0;  W0=500.0), dev, 0.003),
    LifeStage(:larva,  DistributedDelay(25, 200.0; W0=0.0),   dev, 0.008),
    LifeStage(:pupa,   DistributedDelay(10, 70.0;  W0=0.0),   dev, 0.002),
    LifeStage(:adult,  DistributedDelay(5,  150.0; W0=0.0),   dev, 0.015),
]
pop = Population(:moth, stages)

# Or use the convenience constructor
pop2 = make_population(:moth,
    [(:egg, 15, 80.0), (:larva, 25, 200.0), (:pupa, 10, 70.0), (:adult, 5, 150.0)],
    dev; mortality=0.005)

# Simulate over a growing season with synthetic weather
weather = SinusoidalWeather(20.0, 10.0; phase=200.0)
prob = PBDMProblem(pop, weather, (1, 365))
sol = solve(prob, DirectIteration())

# When does 50% of maturation occur? (phenological midpoint)
mid = phenology(sol; threshold=0.5)
println("Phenological midpoint: day ", mid)

# Stage-specific trajectories
for i in 1:4
    traj = stage_trajectory(sol, i)
    peak_day = sol.t[argmax(traj)]
    println("$(pop.stages[i].name): peak at day $peak_day")
end

fr = FraserGilbertResponse(0.7)

# When supply greatly exceeds demand → acquire ≈ demand
println(acquire(fr, 1000.0, 10.0))   # ≈ 10.0

# When supply is scarce → acquire << demand
println(acquire(fr, 5.0, 100.0))     # ≈ 3.4

# The supply/demand ratio φ scales vital rates
φ = supply_demand_ratio(fr, 50.0, 100.0)
println("φ = ", round(φ, digits=3))  # ≈ 0.295

pool = MetabolicPool(
    80.0,                                         # supply
    [25.0, 30.0, 40.0, 20.0],                    # demands
    [:respiration, :growth, :reproduction, :reserves]
)

alloc = allocate(pool)
println("Respiration: ", alloc[1])   # 25.0 (full)
println("Growth:      ", alloc[2])   # 30.0 (full)
println("Reproduction:", alloc[3])   # 25.0 (partial!)
println("Reserves:    ", alloc[4])   # 0.0  (nothing left)

φ = supply_demand_index(pool)
println("φ = ", round(φ, digits=3))  # 0.696

# Leaf respiration: 3% of dry mass/day at 25°C, Q₁₀ = 2.3
leaf_resp = Q10Respiration(0.03, 2.3, 25.0)

println("At 15°C: ", round(respiration_rate(leaf_resp, 15.0), digits=4))
println("At 25°C: ", round(respiration_rate(leaf_resp, 25.0), digits=4))
println("At 35°C: ", round(respiration_rate(leaf_resp, 35.0), digits=4))

