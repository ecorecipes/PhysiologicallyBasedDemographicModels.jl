using PhysiologicallyBasedDemographicModels

dev = LinearDevelopmentRate(10.0, 35.0)
pop = Population(:demo, [
    LifeStage(:juvenile, DistributedDelay(20, 180.0; W0=5.0), dev, 0.001),
    LifeStage(:adult, DistributedDelay(20, 320.0; W0=1.0), dev, 0.0005),
])

weather = WeatherSeries(fill(24.0, 20); day_offset=1)
prob = PBDMProblem(pop, weather, (1, 20))
sol = solve(prob, DirectIteration())

println("Basic delay example")
println("  final population: ", round(total_population(sol)[end], digits=3))
println("  cumulative degree-days: ", round(cumulative_degree_days(sol)[end], digits=1))
println("  stage totals on day 20: ", round.(sol.u[end], digits=3))
