using PhysiologicallyBasedDemographicModels

seedling_dev = LinearDevelopmentRate(8.0, 30.0)
adult_dev = LinearDevelopmentRate(12.0, 35.0)

pop = Population(:heterogeneous_dd, [
    LifeStage(:seedling, DistributedDelay(15, 120.0; W0=2.0), seedling_dev, 0.02),
    LifeStage(:adult, DistributedDelay(15, 250.0; W0=0.5), adult_dev, 0.01),
])

weather = WeatherSeries(fill(20.0, 14); day_offset=1)
prob = PBDMProblem(pop, weather, (1, 14))
sol = solve(prob, DirectIteration())

println("Stage-specific degree-day example")
println("  cumulative DD (stage 1): ", round(cumulative_degree_days(sol; stage_idx=1)[end], digits=1))
println("  cumulative DD (stage 2): ", round(cumulative_degree_days(sol; stage_idx=2)[end], digits=1))
println("  day-1 DD by stage: ", round.([stage_degree_days(sol, 1)[1], stage_degree_days(sol, 2)[1]], digits=1))
