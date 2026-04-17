@testset "Coupled Population Systems" begin
    dev = LinearDevelopmentRate(10.0, 35.0)

    function make_test_pop(name, w0_juv, w0_adult)
        Population(name, [
            LifeStage(:juvenile, DistributedDelay(10, 100.0; W0=w0_juv), dev, 0.01),
            LifeStage(:adult,    DistributedDelay(10, 200.0; W0=w0_adult), dev, 0.005),
        ])
    end

    @testset "PopulationSystem construction" begin
        pop1 = make_test_pop(:a, 50.0, 25.0)
        pop2 = make_test_pop(:b, 30.0, 15.0)

        sys = PopulationSystem(
            :comp_a => pop1,
            :comp_b => pop2,
        )
        @test length(sys) == 2
        @test collect(keys(sys)) == [:comp_a, :comp_b]
        @test haskey(sys, :comp_a)
        @test !haskey(sys, :comp_c)
        @test total_population(sys) ≈ 750.0 + 450.0
        @test component_total(sys, :comp_a) ≈ 750.0
        @test component_total(sys, :comp_b) ≈ 450.0

        totals = component_totals(sys)
        @test totals[:comp_a] ≈ 750.0
        @test totals[:comp_b] ≈ 450.0

        # Duplicate name should error
        @test_throws ArgumentError PopulationSystem(:x => pop1, :x => pop2)
        # Empty should error
        @test_throws ArgumentError PopulationSystem()
    end

    @testset "PopulationComponent metadata" begin
        pop = make_test_pop(:test, 10.0, 5.0)
        comp = PopulationComponent(pop; species=:insect, type=:wild, patch=:field1)
        @test comp.species == :insect
        @test comp.type == :wild
        @test comp.patch == :field1

        sys = PopulationSystem(
            :wild => PopulationComponent(pop; species=:sw, type=:wild),
            :sterile => PopulationComponent(make_test_pop(:s, 5.0, 2.0); species=:sw, type=:sterile),
        )
        wild_comps = by_type(sys, :wild)
        @test length(wild_comps) == 1
        @test wild_comps[1][1] == :wild

        sw_comps = by_species(sys, :sw)
        @test length(sw_comps) == 2
    end

    @testset "inject! and remove_fraction!" begin
        pop = make_test_pop(:test, 100.0, 50.0)
        sys = PopulationSystem(:p => pop)

        before = component_total(sys, :p)
        inject!(sys, :p, 20.0)
        @test component_total(sys, :p) ≈ before + 20.0

        inject!(sys, :p, 2, 10.0)  # into adult stage
        @test component_total(sys, :p) ≈ before + 30.0

        total_before = component_total(sys, :p)
        remove_fraction!(sys, :p, 0.1)
        @test component_total(sys, :p) ≈ total_before * 0.9 atol=0.1

        @test_throws BoundsError inject!(sys, :p, 5, 10.0)
        @test_throws ArgumentError remove_fraction!(sys, :p, 1.5)
    end

    @testset "PulseRelease event" begin
        pop = make_test_pop(:t, 10.0, 5.0)
        sys = PopulationSystem(:t => pop)

        ev = PulseRelease(:t, 100.0, 7; start_day=1)

        # Day 1 fires (offset 0 from start_day)
        @test apply_event!(ev, sys, nothing, 1, nothing) == true
        # Day 2 does not fire
        @test apply_event!(ev, sys, nothing, 2, nothing) == false
        # Day 8 fires (offset 7 from start_day)
        @test apply_event!(ev, sys, nothing, 8, nothing) == true
    end

    @testset "SingleDayRelease event" begin
        pop = make_test_pop(:t, 10.0, 5.0)
        sys = PopulationSystem(:t => pop)

        ev = SingleDayRelease(:t, 50.0, 10)
        before = component_total(sys, :t)
        @test apply_event!(ev, sys, nothing, 9, nothing) == false
        @test apply_event!(ev, sys, nothing, 10, nothing) == true
        @test component_total(sys, :t) ≈ before + 50.0
        @test apply_event!(ev, sys, nothing, 11, nothing) == false
    end

    @testset "SprayEvent" begin
        pop1 = make_test_pop(:a, 100.0, 50.0)
        pop2 = make_test_pop(:b, 100.0, 50.0)
        sys = PopulationSystem(:a => pop1, :b => pop2)

        ev = SprayEvent([:a, :b], [0.5, 0.3], [5])
        before_a = component_total(sys, :a)
        before_b = component_total(sys, :b)

        @test apply_event!(ev, sys, nothing, 4, nothing) == false
        @test apply_event!(ev, sys, nothing, 5, nothing) == true
        @test component_total(sys, :a) ≈ before_a * 0.5 atol=0.1
        @test component_total(sys, :b) ≈ before_b * 0.7 atol=0.1
    end

    @testset "ReproductionRule" begin
        pop = make_test_pop(:t, 10.0, 5.0)
        sys = PopulationSystem(:t => pop)
        w = DailyWeather(25.0, 20.0, 30.0)

        rule = ReproductionRule(:t, (sys, w, day, p) -> 10.0)
        before = component_total(sys, :t)
        result = apply_rule!(rule, sys, w, 1, nothing)
        @test result.offspring ≈ 10.0
        @test component_total(sys, :t) ≈ before + 10.0
    end

    @testset "MortalityRule" begin
        pop = make_test_pop(:t, 100.0, 50.0)
        sys = PopulationSystem(:t => pop)
        w = DailyWeather(25.0, 20.0, 30.0)

        rule = MortalityRule(:t, (sys, w, day, p) -> 0.1)
        before = component_total(sys, :t)
        result = apply_rule!(rule, sys, w, 1, nothing)
        @test result.mortality ≈ 0.1
        @test component_total(sys, :t) ≈ before * 0.9 atol=1.0
    end

    @testset "TransferRule" begin
        pop1 = make_test_pop(:src, 100.0, 50.0)
        pop2 = make_test_pop(:dst, 10.0, 5.0)
        sys = PopulationSystem(:src => pop1, :dst => pop2)
        w = DailyWeather(25.0, 20.0, 30.0)

        rule = TransferRule(:src, :dst, (sys, w, day, p) -> 0.5)
        src_before = delay_total(sys[:src].population.stages[1].delay)
        dst_before = delay_total(sys[:dst].population.stages[1].delay)

        result = apply_rule!(rule, sys, w, 1, nothing)
        @test result.transferred ≈ src_before * 0.5 atol=1.0
    end

    @testset "CustomRule" begin
        pop = make_test_pop(:t, 10.0, 5.0)
        sys = PopulationSystem(:t => pop)
        w = DailyWeather(25.0, 20.0, 30.0)

        rule = CustomRule(:my_rule, (sys, w, day, p) -> (metric=day * 2.0,))
        result = apply_rule!(rule, sys, w, 5, nothing)
        @test result.metric ≈ 10.0
    end

    @testset "Coupled solver - basic" begin
        pop1 = make_test_pop(:prey, 100.0, 50.0)
        pop2 = make_test_pop(:pred, 10.0, 5.0)

        sys = PopulationSystem(
            :prey => PopulationComponent(pop1; species=:prey),
            :pred => PopulationComponent(pop2; species=:predator),
        )

        days = [DailyWeather(25.0, 20.0, 30.0; radiation=20.0, photoperiod=14.0) for _ in 1:30]
        weather = WeatherSeries(days)

        repro = ReproductionRule(:prey,
            (sys, w, day, p) -> delay_total(sys[:prey].population.stages[2].delay) * 0.05)

        obs = Observable(:total, (sys, w, day, p) -> total_population(sys))

        prob = PBDMProblem(MultiSpeciesPBDMNew(), sys, weather, (1, 30);
            rules = AbstractInteractionRule[repro],
            events = AbstractScheduledEvent[],
            observables = [obs]
        )

        sol = solve(prob, DirectIteration())
        @test sol.retcode == :Success
        @test length(sol.t) == 30
        @test length(sol[:prey]) == 30
        @test length(sol[:pred]) == 30
        @test length(sol.observables[:total]) == 30
        @test all(sol[:prey] .>= 0)
        @test all(sol[:pred] .>= 0)
    end

    @testset "Coupled solver - with events" begin
        pop = make_test_pop(:target, 100.0, 50.0)
        sys = PopulationSystem(:target => pop)

        days = [DailyWeather(25.0, 20.0, 30.0; radiation=20.0) for _ in 1:20]
        weather = WeatherSeries(days)

        release = PulseRelease(:target, 50.0, 5; start_day=1)

        prob = PBDMProblem(MultiTypePBDM(), sys, weather, (1, 20);
            events = AbstractScheduledEvent[release]
        )

        sol = solve(prob, DirectIteration())
        @test sol.retcode == :Success
        # PulseRelease fires on days 1, 6, 11, 16 => 4 events
        @test length(sol.event_log) == 4
    end

    @testset "Coupled solver - predation" begin
        pop_prey = Population(:prey, [
            LifeStage(:imm, DistributedDelay(10, 100.0; W0=200.0), dev, 0.005),
            LifeStage(:adult, DistributedDelay(10, 200.0; W0=100.0), dev, 0.003),
        ])
        pop_pred = Population(:pred, [
            LifeStage(:imm, DistributedDelay(8, 80.0; W0=20.0), dev, 0.005),
            LifeStage(:adult, DistributedDelay(8, 150.0; W0=10.0), dev, 0.003),
        ])

        sys = PopulationSystem(:prey => pop_prey, :pred => pop_pred)
        days = [DailyWeather(25.0, 20.0, 30.0; radiation=20.0) for _ in 1:50]
        weather = WeatherSeries(days)

        pred_rule = PredationRule(:pred, :prey,
            HollingTypeII(0.3, 0.01); conversion=0.2)

        prob = PBDMProblem(MultiSpeciesPBDMNew(), sys, weather, (1, 50);
            rules = AbstractInteractionRule[pred_rule]
        )

        sol = solve(prob, DirectIteration())
        @test sol.retcode == :Success
        # Prey should decline faster with predation
        @test sol[:prey][end] < sol[:prey][1]
    end

    @testset "Coupled solver - convenience constructor" begin
        pop = make_test_pop(:x, 10.0, 5.0)
        sys = PopulationSystem(:x => pop)
        days = [DailyWeather(25.0, 20.0, 30.0) for _ in 1:10]
        weather = WeatherSeries(days)

        # Default to MultiTypePBDM
        prob = PBDMProblem(sys, weather, (1, 10))
        sol = solve(prob, DirectIteration())
        @test sol.retcode == :Success
        @test length(sol.t) == 10
    end

    @testset "StressRule construction" begin
        sr = StressRule(:test_stress, (sys, w, day, p) -> Dict{Symbol, Vector{Float64}}())
        @test sr.name == :test_stress
        @test sr isa AbstractStressRule
    end

    @testset "StressRule - solver integration" begin
        pop1 = make_test_pop(:a, 100.0, 50.0)
        pop2 = make_test_pop(:b, 100.0, 50.0)

        # Run WITHOUT stress for baseline
        sys_no = PopulationSystem(:a => pop1, :b => pop2)
        days = [DailyWeather(25.0, 20.0, 30.0) for _ in 1:20]
        weather = WeatherSeries(days)

        prob_no = PBDMProblem(MultiSpeciesPBDMNew(), sys_no, weather, (1, 20))
        sol_no = solve(prob_no, DirectIteration())

        # Run WITH stress on component :a
        pop1s = make_test_pop(:a, 100.0, 50.0)
        pop2s = make_test_pop(:b, 100.0, 50.0)
        sys_s = PopulationSystem(:a => pop1s, :b => pop2s)

        stress_rule = StressRule(:heavy_stress, (sys, w, day, p) ->
            Dict(:a => [0.5, 0.5])  # heavy stress on component :a
        )

        prob_s = PBDMProblem(MultiSpeciesPBDMNew(), sys_s, weather, (1, 20);
                             stress_rules=AbstractStressRule[stress_rule])
        sol_s = solve(prob_s, DirectIteration())

        @test sol_s.retcode == :Success
        # Stressed population should be smaller
        @test sol_s[:a][end] < sol_no[:a][end]
        # Unstressed population should be similar
        @test sol_s[:b][end] ≈ sol_no[:b][end] atol=1e-6
    end

    @testset "StressRule - multiple rules merge" begin
        pop = make_test_pop(:x, 100.0, 50.0)
        sys = PopulationSystem(:x => pop)
        days = [DailyWeather(25.0, 20.0, 30.0) for _ in 1:10]
        weather = WeatherSeries(days)

        sr1 = StressRule(:light, (sys, w, day, p) -> Dict(:x => [0.1, 0.2]))
        sr2 = StressRule(:nitrogen, (sys, w, day, p) -> Dict(:x => [0.3, 0.1]))

        prob = PBDMProblem(MultiSpeciesPBDMNew(), sys, weather, (1, 10);
                           stress_rules=AbstractStressRule[sr1, sr2])
        sol = solve(prob, DirectIteration())
        @test sol.retcode == :Success
        # With merged stress (max), population should decline faster
        @test sol[:x][end] < 150.0  # initial total = 150
    end

    # ========================================================================
    # StateVariable tests
    # ========================================================================

    @testset "ScalarState construction and access" begin
        ss = ScalarState(:temperature, 25.0)
        @test get_state(ss) == 25.0
        @test ss.name == :temperature
        @test !has_auto_update(ss)

        set_state!(ss, 30.0)
        @test get_state(ss) == 30.0

        @test snapshot(ss) == 30.0
    end

    @testset "ScalarState auto-update" begin
        ss = ScalarState(:cum_dd, 0.0; update=(val, sys, w, day, p) -> val + 10.0)
        @test has_auto_update(ss)

        pop = make_test_pop(:a, 50.0, 25.0)
        sys = PopulationSystem(:a => pop)
        w = DailyWeather(25.0)
        update_state!(ss, sys, w, 1, nothing)
        @test get_state(ss) ≈ 10.0
        update_state!(ss, sys, w, 2, nothing)
        @test get_state(ss) ≈ 20.0
    end

    @testset "ArrayState construction and access" begin
        as = ArrayState(:carbon, [1.0, 2.0, 3.0, 4.0])
        @test get_state(as) == [1.0, 2.0, 3.0, 4.0]
        @test as.name == :carbon
        @test !has_auto_update(as)
        @test length(snapshot(as)) == 4
        # snapshot should be a copy
        s = snapshot(as)
        s[1] = 999.0
        @test get_state(as)[1] == 1.0
    end

    @testset "ArrayState auto-update" begin
        as = ArrayState(:pools, [10.0, 20.0]; update=(v, sys, w, day, p) -> begin
            v .+= 1.0
        end)
        @test has_auto_update(as)

        pop = make_test_pop(:a, 50.0, 25.0)
        sys = PopulationSystem(:a => pop)
        w = DailyWeather(25.0)
        update_state!(as, sys, w, 1, nothing)
        @test get_state(as) ≈ [11.0, 21.0]
    end

    @testset "DictState construction and access" begin
        ds = DictState(:patch_water, Dict(:north => 100.0, :south => 50.0))
        @test get_state(ds)[:north] == 100.0
        @test ds.name == :patch_water
        @test !has_auto_update(ds)

        s = snapshot(ds)
        s[:north] = 0.0
        @test get_state(ds)[:north] == 100.0
    end

    @testset "DictState auto-update" begin
        ds = DictState(:nutrient, Dict(:N => 50.0); update=(d, sys, w, day, p) -> begin
            d[:N] -= 1.0
        end)
        @test has_auto_update(ds)

        pop = make_test_pop(:a, 50.0, 25.0)
        sys = PopulationSystem(:a => pop)
        w = DailyWeather(25.0)
        update_state!(ds, sys, w, 1, nothing)
        @test get_state(ds)[:N] ≈ 49.0
    end

    # ========================================================================
    # PopulationSystem with state variables
    # ========================================================================

    @testset "PopulationSystem with state" begin
        pop = make_test_pop(:a, 50.0, 25.0)
        ss = ScalarState(:allele_freq, 0.01)
        as = ArrayState(:carbon, [1.0, 2.0])

        sys = PopulationSystem(:a => pop; state=[ss, as])
        @test has_state(sys, :allele_freq)
        @test has_state(sys, :carbon)
        @test !has_state(sys, :nonexistent)
        @test get_state(sys, :allele_freq) ≈ 0.01
        set_state!(sys, :allele_freq, 0.05)
        @test get_state(sys, :allele_freq) ≈ 0.05

        # show includes state info
        s = sprint(show, sys)
        @test occursin("state vars", s)
    end

    # ========================================================================
    # BulkPopulation tests
    # ========================================================================

    @testset "BulkPopulation construction" begin
        bp = BulkPopulation(:inoculum, 100.0)
        @test total_population(bp) ≈ 100.0
        @test n_stages(bp) == 1
        @test bp.name == :inoculum
    end

    @testset "BulkPopulation inject and remove" begin
        bp = BulkPopulation(:pest, 50.0)
        inject!(bp, 20.0)
        @test total_population(bp) ≈ 70.0

        inject!(bp, 1, 10.0)  # stage arg ignored
        @test total_population(bp) ≈ 80.0

        remove_fraction!(bp, 0.5)
        @test total_population(bp) ≈ 40.0

        remove_fraction!(bp, 1, 0.25)  # stage arg ignored
        @test total_population(bp) ≈ 30.0
    end

    @testset "BulkPopulation carrying capacity" begin
        bp = BulkPopulation(:scale, 90.0; K=100.0)
        inject!(bp, 50.0)
        @test total_population(bp) ≈ 100.0  # capped at K
    end

    @testset "BulkPopulation growth_fn stepping" begin
        # Logistic growth: N_new = N * (1 + r*(1 - N/K))
        bp = BulkPopulation(:logistic, 10.0;
            growth_fn=(N, w, d, p) -> N * (1.0 + 0.1*(1.0 - N/1000.0)),
            K=1000.0)
        w = DailyWeather(25.0)
        step_bulk!(bp, w, 1, nothing)
        @test total_population(bp) > 10.0

        # No growth_fn → no-op
        bp2 = BulkPopulation(:static, 50.0)
        step_bulk!(bp2, w, 1, nothing)
        @test total_population(bp2) ≈ 50.0
    end

    @testset "BulkPopulation in PopulationSystem" begin
        pop = make_test_pop(:arthropod, 50.0, 25.0)
        bp = BulkPopulation(:inoculum, 200.0)

        sys = PopulationSystem(:arthropod => pop, :inoculum => bp)
        @test length(sys) == 2
        @test component_total(sys, :inoculum) ≈ 200.0
        @test total_population(sys) ≈ 750.0 + 200.0

        inject!(sys, :inoculum, 50.0)
        @test component_total(sys, :inoculum) ≈ 250.0

        remove_fraction!(sys, :inoculum, 0.4)
        @test component_total(sys, :inoculum) ≈ 150.0
    end

    # ========================================================================
    # WeatherConditionalEvent tests
    # ========================================================================

    @testset "WeatherConditionalEvent fires on predicate" begin
        pop = make_test_pop(:pest, 100.0, 0.0)
        sys = PopulationSystem(:pest => pop)

        # Only fire when photoperiod < 11
        fired_count = Ref(0)
        wce = WeatherConditionalEvent(:short_day,
            (w, day, sys, p) -> w.photoperiod < 11.0,
            (sys, w, day, p) -> begin
                fired_count[] += 1
                true
            end
        )

        w_long = DailyWeather(25.0, 20.0, 30.0; photoperiod=14.0)
        @test apply_event!(wce, sys, w_long, 1, nothing) == false
        @test fired_count[] == 0

        w_short = DailyWeather(25.0, 20.0, 30.0; photoperiod=10.0)
        @test apply_event!(wce, sys, w_short, 2, nothing) == true
        @test fired_count[] == 1
    end

    @testset "WeatherConditionalEvent with rainfall" begin
        pop = make_test_pop(:pathogen, 10.0, 0.0)
        sys = PopulationSystem(:pathogen => pop)

        wce = WeatherConditionalEvent(:rain_trigger,
            (w, day, sys, p) -> w.rainfall > 5.0,
            (sys, w, day, p) -> begin
                inject!(sys, :pathogen, 50.0)
                true
            end
        )

        w_dry = DailyWeather(25.0, 20.0, 30.0; rainfall=1.0)
        @test apply_event!(wce, sys, w_dry, 1, nothing) == false

        w_wet = DailyWeather(25.0, 20.0, 30.0; rainfall=10.0)
        @test apply_event!(wce, sys, w_wet, 2, nothing) == true
        @test component_total(sys, :pathogen) > 50.0  # injected 50 on top of existing
    end

    # ========================================================================
    # PhaseCallback tests
    # ========================================================================

    @testset "PhaseCallback construction" begin
        cb = PhaseCallback(:my_hook, PRE_STEP, (sys, w, d, p) -> nothing)
        @test cb.name == :my_hook
        @test cb.phase == PRE_STEP
    end

    @testset "PhaseCallback ordering in solver" begin
        pop = make_test_pop(:x, 50.0, 25.0)
        weather = SinusoidalWeather(25.0, 5.0)
        execution_order = Symbol[]

        sys = PopulationSystem(:x => pop)

        cbs = PhaseCallback[
            PhaseCallback(:pre_ev, PRE_EVENT, (s, w, d, p) -> push!(execution_order, :pre_event)),
            PhaseCallback(:post_ev, POST_EVENT, (s, w, d, p) -> push!(execution_order, :post_event)),
            PhaseCallback(:pre_st, PRE_STEP, (s, w, d, p) -> push!(execution_order, :pre_step)),
            PhaseCallback(:post_st, POST_STEP, (s, w, d, p) -> push!(execution_order, :post_step)),
            PhaseCallback(:eod, END_OF_DAY, (s, w, d, p) -> push!(execution_order, :end_of_day)),
        ]

        prob = PBDMProblem(MultiTypePBDM(), sys, weather, (1, 1);
                           callbacks=cbs)
        sol = solve(prob, DirectIteration())
        @test sol.retcode == :Success

        # Verify phase ordering for a single day
        @test execution_order == [:pre_event, :post_event, :pre_step, :post_step, :end_of_day]
    end

    # ========================================================================
    # Solver integration tests for new features
    # ========================================================================

    @testset "Solver with state variables records history" begin
        pop = make_test_pop(:x, 50.0, 25.0)
        weather = SinusoidalWeather(25.0, 5.0)

        cum_dd_state = ScalarState(:cum_dd, 0.0;
            update=(val, sys, w, day, p) -> val + max(0.0, w.T_mean - 10.0))

        sys = PopulationSystem(:x => pop; state=[cum_dd_state])

        prob = PBDMProblem(MultiTypePBDM(), sys, weather, (1, 30))
        sol = solve(prob, DirectIteration())
        @test sol.retcode == :Success
        @test haskey(sol.state_history, :cum_dd)
        @test length(sol.state_history[:cum_dd]) == 30
        # cumulative DD should be monotonically increasing
        h = sol.state_history[:cum_dd]
        @test all(h[i] <= h[i+1] for i in 1:29)
        @test h[end] > 0
    end

    @testset "Solver with BulkPopulation component" begin
        pop = make_test_pop(:arthropod, 50.0, 25.0)
        bp = BulkPopulation(:inoculum, 100.0;
            growth_fn=(N, w, d, p) -> N * 1.02,  # 2% daily growth
            K=10000.0)
        weather = SinusoidalWeather(25.0, 5.0)

        sys = PopulationSystem(:arthropod => pop, :inoculum => bp)

        prob = PBDMProblem(MultiSpeciesPBDMNew(), sys, weather, (1, 50))
        sol = solve(prob, DirectIteration())
        @test sol.retcode == :Success
        @test sol[:inoculum][end] > 100.0  # grew
        @test sol[:inoculum][end] < 10000.0  # didn't hit cap yet in 50 days
        @test length(sol.component_names) == 2
        # Stage totals for BulkPopulation: 1 row
        @test size(sol.component_stage_totals[:inoculum]) == (1, 50)
    end

    @testset "Solver with WeatherConditionalEvent" begin
        pop = make_test_pop(:moth, 100.0, 50.0)
        weather = WeatherSeries([
            DailyWeather(25.0, 20.0, 30.0; photoperiod=14.0),
            DailyWeather(25.0, 20.0, 30.0; photoperiod=14.0),
            DailyWeather(15.0, 10.0, 20.0; photoperiod=10.0),  # short day
            DailyWeather(15.0, 10.0, 20.0; photoperiod=10.0),
            DailyWeather(15.0, 10.0, 20.0; photoperiod=10.0),
        ])

        sys = PopulationSystem(:moth => pop)
        wce = WeatherConditionalEvent(:diapause_kill,
            (w, day, sys, p) -> w.photoperiod < 11.0,
            (sys, w, day, p) -> begin
                remove_fraction!(sys, :moth, 0.3)
                true
            end
        )

        prob = PBDMProblem(MultiTypePBDM(), sys, weather, (1, 5);
                           events=AbstractScheduledEvent[wce])
        sol = solve(prob, DirectIteration())
        @test sol.retcode == :Success
        # Should have logged WeatherConditionalEvent firings on days 3, 4, 5
        wce_events = filter(e -> e[2] == :diapause_kill, sol.event_log)
        @test length(wce_events) == 3
        @test all(e[1] >= 3 for e in wce_events)
    end

    @testset "Solver with all new features combined" begin
        pop = make_test_pop(:pest, 50.0, 25.0)
        bp = BulkPopulation(:soil_N, 100.0)
        weather = SinusoidalWeather(25.0, 5.0)

        n_state = ScalarState(:cumulative_N_uptake, 0.0;
            update=(val, sys, w, day, p) -> val + 0.5)

        callback_log = Int[]
        cb = PhaseCallback(:log_day, POST_STEP, (sys, w, d, p) -> push!(callback_log, d))

        sys = PopulationSystem(
            :pest => pop,
            :soil_N => bp;
            state=[n_state]
        )

        # Rule that transfers N from soil to uptake
        rule = CustomRule(:uptake, (sys, w, day, p) -> begin
            if component_total(sys, :soil_N) > 1.0
                remove_fraction!(sys, :soil_N, 0.01)
            end
            (uptake=0.01,)
        end)

        prob = PBDMProblem(MultiSpeciesPBDMNew(), sys, weather, (1, 20);
                           rules=AbstractInteractionRule[rule],
                           callbacks=PhaseCallback[cb])
        sol = solve(prob, DirectIteration())
        @test sol.retcode == :Success
        @test length(callback_log) == 20
        @test sol.state_history[:cumulative_N_uptake][end] ≈ 10.0
        @test sol[:soil_N][end] < 100.0  # some N was removed
        @test haskey(sol.rule_log, :uptake)
    end

    # ================================================================
    # New extension types
    # ================================================================

    @testset "GenomeState" begin
        locus = DialleleicLocus(0.05, 0.5)
        gs = GenomeState(:resistance, locus)
        @test get_state(gs) ≈ 0.05
        freqs = get_genotypes(gs)
        @test freqs.RR ≈ 0.05^2
        @test freqs.SR ≈ 2 * 0.05 * 0.95
        @test freqs.SS ≈ 0.95^2
        @test freqs.RR + freqs.SR + freqs.SS ≈ 1.0
        set_state!(gs, 0.3)
        @test get_state(gs) ≈ 0.3
        @test snapshot(gs) ≈ 0.3
        @test get_locus(gs) === locus
        @test !has_auto_update(gs)

        # With auto-update
        gs2 = GenomeState(:r2, DialleleicLocus(0.1, 0.0);
            update=(locus, sys, w, day, p) -> (locus.R = min(1.0, locus.R + 0.01)))
        @test has_auto_update(gs2)

        # update_state! clamps to [0, 1] even if user update_fn overshoots
        gs_overshoot = GenomeState(:r3, DialleleicLocus(0.5, 0.0);
            update=(locus, sys, w, day, p) -> (locus.R = 2.0))
        bp = BulkPopulation(:pest, 100.0)
        sys_clamp = PopulationSystem(:pest => bp; state=[gs_overshoot])
        update_state!(sys_clamp.state[:r3], sys_clamp, DailyWeather(20.0, 10.0, 15.0), 1, nothing)
        @test get_state(sys_clamp.state[:r3]) == 1.0

        # SelectionRule: Hardy-Weinberg selection wired into coupled rule phase
        gs_sel = GenomeState(:resistance_sel, DialleleicLocus(0.1, 0.5))
        bp_sel = BulkPopulation(:pest, 100.0)
        sys_sel = PopulationSystem(:pest => bp_sel; state=[gs_sel])
        # Strongly favor RR (R-homozygous)
        fitness_fn = (sys, w, day, p) -> GenotypeFitness(0.1, 0.5, 1.0)
        rule = SelectionRule(:resistance_sel, fitness_fn; name=:sel)
        result = apply_rule!(rule, sys_sel, DailyWeather(20.0, 10.0, 15.0), 1, nothing)
        @test haskey(result, :R_before)
        @test haskey(result, :R_after)
        @test haskey(result, :w_bar)
        @test result.R_before ≈ 0.1
        @test result.R_after > result.R_before
        @test 0.0 <= result.R_after <= 1.0

        # SelectionRule errors when state missing
        sys_empty = PopulationSystem(:pest => BulkPopulation(:pest, 1.0))
        rule_missing = SelectionRule(:missing_locus, fitness_fn)
        @test_throws ArgumentError apply_rule!(rule_missing, sys_empty,
            DailyWeather(20.0, 10.0, 15.0), 1, nothing)

        # SelectionRule errors when fitness_fn returns wrong type
        bad_fn = (sys, w, day, p) -> 0.5
        rule_bad = SelectionRule(:resistance_sel, bad_fn)
        @test_throws ArgumentError apply_rule!(rule_bad, sys_sel,
            DailyWeather(20.0, 10.0, 15.0), 1, nothing)
    end

    @testset "DiapauseState + DiapauseRule" begin
        ds = DiapauseState(:winter_pool, 100.0; cold_dd=50.0,
            induction_fn=(dl, T) -> dl < 12.0 ? 0.1 : 0.0,
            emergence_fn=(dd) -> dd > 200.0 ? 0.05 : 0.0,
            cold_survival_fn=(cdd) -> max(0.0, 1.0 - 0.001 * cdd))
        @test get_state(ds) ≈ 100.0
        @test snapshot(ds).pool ≈ 100.0
        @test snapshot(ds).cold_dd ≈ 50.0
        @test !has_auto_update(ds)

        # Build a system with BulkPopulation + DiapauseState
        bp = BulkPopulation(:moths, 500.0)
        sys = PopulationSystem(:moths => bp; state=[ds])
        @test has_state(sys, :winter_pool)

        rule = DiapauseRule(:moths, :winter_pool; T_cold_base=10.0)
        w_cold = DailyWeather(5.0, 2.0, 8.0)
        w_cold_short = DailyWeather(5.0, 2.0, 8.0)  # photoperiod defaults to 12.0

        result = apply_rule!(rule, sys, w_cold_short, 1, nothing)
        @test haskey(result, :pool)
        @test haskey(result, :emerged)
        @test haskey(result, :survival)
        @test haskey(result, :entering)
        @test result.entering >= 0.0
    end

    @testset "PhenologyState" begin
        ps = PhenologyState(:cotton,
            [(:squaring, 450.0), (:flowering, 700.0), (:boll_fill, 1000.0)];
            base_temp=12.0)
        @test get_state(ps) ≈ 0.0
        @test get_phase(ps) == :pre
        @test !past_milestone(ps, :squaring)
        @test resource_availability(ps, :flowering) ≈ 0.0

        # Simulate auto-update
        bp = BulkPopulation(:pest, 10.0)
        sys = PopulationSystem(:pest => bp; state=[ps])
        w = DailyWeather(30.0)  # 18 DD above base_temp=12
        for d in 1:30
            update_state!(ps, sys, w, d, nothing)
        end
        @test get_state(ps) ≈ 30 * 18.0
        @test get_phase(ps) == :squaring
        @test past_milestone(ps, :squaring)
        @test !past_milestone(ps, :flowering)
        @test resource_availability(ps, :squaring; ramp_dd=100.0) ≈ 0.9

        @test has_auto_update(ps)
        @test snapshot(ps).phase == :squaring

        # PhenologyState with custom dd_fn (e.g. a nonlinear dev rate)
        # Here dd_fn returns a constant 5 DD/day regardless of weather.
        ps_custom = PhenologyState(:custom, [(:m1, 50.0), (:m2, 100.0)];
            dd_fn=(w, day, p) -> 5.0)
        sys_custom = PopulationSystem(:pest => BulkPopulation(:pest, 1.0);
            state=[ps_custom])
        for d in 1:15
            update_state!(ps_custom, sys_custom, DailyWeather(0.0), d, nothing)
        end
        @test get_state(ps_custom) ≈ 75.0
        @test get_phase(ps_custom) == :m1
        @test !past_milestone(ps_custom, :m2)
    end

    @testset "SoilState" begin
        ss = SoilState(:pathogen; inoculum=10.0, virulence=0.2)
        @test get_inoculum(ss) ≈ 10.0
        @test get_virulence(ss) ≈ 0.2
        @test get_state(ss) == (inoculum=10.0, virulence=0.2)
        @test !has_auto_update(ss)

        set_state!(ss, (inoculum=5.0, virulence=0.3))
        @test get_inoculum(ss) ≈ 5.0
        @test get_virulence(ss) ≈ 0.3
        @test snapshot(ss) == (inoculum=5.0, virulence=0.3)

        # With auto-update
        ss2 = SoilState(:decay; inoculum=100.0, virulence=0.1,
            update=(inoc, vir, sys, w, day, p) -> (inoc * 0.99, vir))
        @test has_auto_update(ss2)
    end

    @testset "SpatialGrid" begin
        bp1 = BulkPopulation(:pest, 100.0)
        bp2 = BulkPopulation(:pest, 50.0)
        bp3 = BulkPopulation(:pest, 200.0)

        sys1 = PopulationSystem(:pest => bp1)
        sys2 = PopulationSystem(:pest => bp2)
        sys3 = PopulationSystem(:pest => bp3)

        grid = SpatialGrid([:p1 => sys1, :p2 => sys2, :p3 => sys3])
        @test length(grid) == 3
        @test grid[1] === sys1
        @test grid[:p2] === sys2

        # DispersalRule
        dr = DispersalRule(:pest;
            emigration_fn=(N, phi, w, day, p) -> N > 80.0 ? 0.1 : 0.0,
            patch_finding=0.8)
        w = DailyWeather(25.0)
        result = apply_dispersal!(dr, grid, w, 1, nothing)
        @test result.total_emigrants > 0
        @test result.survivors > 0
    end

    @testset "EnsemblePBDMProblem" begin
        bp = BulkPopulation(:pop, 10.0; growth_fn=(N, w, d, p) -> N * (1 + p.r))
        sys = PopulationSystem(:pop => bp)
        wx = WeatherSeries([DailyWeather(25.0) for _ in 1:30])
        prob = PBDMProblem(sys, wx, (1, 30); p=(r=0.05,))

        ens = EnsemblePBDMProblem(prob;
            prob_func = (prob, i, _) -> begin
                r_val = 0.01 * i
                bp_new = BulkPopulation(:pop, 10.0; growth_fn=(N, w, d, p) -> N * (1 + p.r))
                sys_new = PopulationSystem(:pop => bp_new)
                PBDMProblem(sys_new, prob.weather, prob.tspan; p=(r=r_val,))
            end,
            output_func = (sol, i) -> (sol[:pop][end], false))

        results = solve(ens, DirectIteration(); trajectories=5)
        @test length(results) == 5
        @test results[1] < results[5]  # higher r → higher final pop
        @test !results.converged
    end
end
