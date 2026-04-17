using PhysiologicallyBasedDemographicModels

dev = LinearDevelopmentRate(12.0, 35.0)
response = FraserGilbertResponse(0.7)
resp = Q10Respiration(0.016, 2.3, 25.0)

bdf = BiodemographicFunctions(dev, response, resp; label=:cotton_bdf)
mp = MetabolicPool(1.0, [0.8, 0.5, 1.2], [:leaf, :root, :fruit])
hybrid = CoupledPBDMModel(bdf, mp; label=:cotton_hybrid)

cotton = Population(:cotton, [
    LifeStage(:leaf, DistributedDelay(25, 700.0; W0=0.5), dev, 0.0002),
    LifeStage(:root, DistributedDelay(25, 180.0; W0=0.2), dev, 0.0005),
    LifeStage(:fruit, DistributedDelay(25, 800.0; W0=0.0), dev, 0.0002),
])

weather_days = [DailyWeather(24.0 + 4sin(2π * d / 90), 20.0, 30.0; radiation=18.0) for d in 1:90]
weather = WeatherSeries(weather_days; day_offset=1)
prob = PBDMProblem(hybrid, cotton, weather, (1, 90))
sol = solve(prob, DirectIteration())

println("Hybrid cotton example")
println("  approach family: ", approach_family(hybrid))
println("  mean daily lambda: ", round(net_growth_rate(sol), digits=4))
println("  final stage totals: ", round.(sol.u[end], digits=4))
