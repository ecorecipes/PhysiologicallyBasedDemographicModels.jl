@testset "Theory — Analytical Tools" begin

    @testset "Functional response traits" begin
        fg = FraserGilbertResponse(0.8)
        h2 = HollingTypeII(0.5, 0.1)
        h3 = HollingTypeIII(0.5, 10.0, 0.1)

        @test is_ratio_dependent(fg) == true
        @test is_ratio_dependent(h2) == false
        @test is_ratio_dependent(h3) == false
        @test apparency(fg) == 0.8
    end

    @testset "Compensation point" begin
        resp = Q10Respiration(0.02, 2.0, 25.0)

        # At reference temperature, R(25) = 0.02
        φ = compensation_point(resp, 25.0, 0.5)
        @test φ ≈ 0.02 / 0.5  # = 0.04

        # With conversion efficiency
        φ_ε = compensation_point(resp, 25.0, 0.5; conversion_efficiency=0.5)
        @test φ_ε ≈ 0.02 / (0.5 * 0.5)  # = 0.08

        # Higher temperature → higher respiration → higher MCP
        φ_hot = compensation_point(resp, 35.0, 0.5)
        @test φ_hot > φ

        # BDF convenience method
        bdf = BiodemographicFunctions(
            LinearDevelopmentRate(10.0, 35.0),
            FraserGilbertResponse(0.8),
            resp)
        φ_bdf = compensation_point(bdf, 25.0; demand_rate=0.5)
        @test φ_bdf ≈ φ

        # Edge cases
        @test_throws ArgumentError compensation_point(resp, 25.0, 0.0)
        @test_throws ArgumentError compensation_point(resp, 25.0, 0.5;
                                                       conversion_efficiency=0.0)
    end

    @testset "Life history strategy" begin
        # r-selected: growth before reproduction
        pool_r = MetabolicPool(100.0, [20.0, 40.0, 40.0],
                               [:respiration, :growth, :reproduction])
        @test life_history_strategy(pool_r) == :r_selected

        # K-selected: reproduction before growth
        pool_K = MetabolicPool(100.0, [20.0, 40.0, 40.0],
                               [:respiration, :reproduction, :growth])
        @test life_history_strategy(pool_K) == :K_selected

        # Unknown: missing labels
        pool_unk = MetabolicPool(100.0, [50.0, 50.0],
                                 [:respiration, :storage])
        @test life_history_strategy(pool_unk) == :unknown

        # CoupledPBDMModel convenience
        bdf = BiodemographicFunctions(
            LinearDevelopmentRate(10.0, 35.0),
            FraserGilbertResponse(0.8),
            Q10Respiration(0.02, 2.0, 25.0))
        model = CoupledPBDMModel(bdf, pool_r)
        @test life_history_strategy(model) == :r_selected
    end

    @testset "Consumer isocline" begin
        fr = FraserGilbertResponse(0.8)
        resp = Q10Respiration(0.02, 2.0, 25.0)
        M1_range = 0.0:10.0:500.0

        iso = consumer_isocline(fr, resp, 25.0;
                                demand_rate=0.5, M1_range=M1_range)
        @test iso.isocline_type == :consumer
        @test iso.slope !== nothing
        @test iso.slope > 0

        # Line through origin: M₂(0) = 0
        @test iso.consumer_biomass[1] ≈ 0.0
        # Linear: M₂ = slope * M₁
        @test iso.consumer_biomass[end] ≈ iso.slope * M1_range[end]
        # Lengths match
        @test length(iso.resource_biomass) == length(M1_range)

        # Population that cannot persist (MCP ≥ 1)
        resp_hot = Q10Respiration(5.0, 2.0, 25.0)  # Very high respiration
        iso_nan = consumer_isocline(fr, resp_hot, 25.0;
                                    demand_rate=0.5, M1_range=M1_range)
        @test all(isnan, iso_nan.consumer_biomass)
        @test iso_nan.slope === nothing
    end

    @testset "Resource isocline" begin
        fr = FraserGilbertResponse(0.8)
        M1_range = 1.0:5.0:995.0

        iso = resource_isocline(fr;
                                intrinsic_rate=0.1,
                                carrying_capacity=1000.0,
                                consumer_demand_rate=0.5,
                                M1_range=M1_range)
        @test iso.isocline_type == :resource
        @test iso.slope === nothing
        # Should have the hump shape: starts low, peaks, returns low
        @test length(iso.resource_biomass) > 0
        @test length(iso.consumer_biomass) > 0
        # Near M₁=K, growth ≈ 0 so M₂ ≈ 0
        if !isempty(iso.consumer_biomass)
            peak_idx = argmax(iso.consumer_biomass)
            @test peak_idx > 1  # Peak is not at the boundary
        end
    end

    @testset "Equilibrium analysis" begin
        fr = FraserGilbertResponse(0.8)
        resp = Q10Respiration(0.02, 2.0, 25.0)

        # Use higher intrinsic rate so equilibrium is viable
        eq = find_equilibrium(fr, resp, 25.0;
                              intrinsic_rate=1.0,
                              carrying_capacity=1000.0,
                              consumer_demand_rate=0.5,
                              conversion_efficiency=0.4)

        # Should find a valid interior equilibrium
        @test eq.M1_star > 0
        @test eq.M2_star > 0
        @test eq.M1_star < 1000.0  # Below carrying capacity
        @test length(eq.eigenvalues) == 2
        @test eq.classification in [:stable_node, :stable_focus,
                                     :unstable_node, :unstable_focus,
                                     :saddle, :center, :degenerate]
        @test size(eq.jacobian) == (2, 2)

        # Verify equilibrium is consistent with isoclines
        φ_star = compensation_point(resp, 25.0, 0.5; conversion_efficiency=0.4)
        s = fr.a / (0.5 * log(1 / (1 - φ_star)))
        @test eq.M2_star ≈ s * eq.M1_star rtol=1e-8

        # Non-persistent species (MCP ≥ 1) → degenerate result
        resp_extreme = Q10Respiration(5.0, 2.0, 25.0)
        eq_bad = find_equilibrium(fr, resp_extreme, 25.0;
                                  intrinsic_rate=1.0,
                                  carrying_capacity=1000.0,
                                  consumer_demand_rate=0.5)
        @test isnan(eq_bad.M1_star)
        @test eq_bad.classification == :degenerate
    end

    @testset "Classify equilibrium" begin
        using PhysiologicallyBasedDemographicModels: classify_equilibrium

        # Stable node
        @test classify_equilibrium([complex(-2.0, 0.0), complex(-1.0, 0.0)]) == :stable_node
        # Unstable node
        @test classify_equilibrium([complex(2.0, 0.0), complex(1.0, 0.0)]) == :unstable_node
        # Saddle
        @test classify_equilibrium([complex(-2.0, 0.0), complex(1.0, 0.0)]) == :saddle
        # Stable focus
        @test classify_equilibrium([complex(-1.0, 2.0), complex(-1.0, -2.0)]) == :stable_focus
        # Unstable focus
        @test classify_equilibrium([complex(1.0, 2.0), complex(1.0, -2.0)]) == :unstable_focus
        # Center
        @test classify_equilibrium([complex(0.0, 2.0), complex(0.0, -2.0)]) == :center
    end

    @testset "Food web assembly" begin
        # Three consumer species with different MCPs
        sp1 = SpeciesProfile(:efficient_herbivore;
            demand_rate=0.5,
            resp=Q10Respiration(0.01, 2.0, 25.0),  # Low respiration → low MCP
            fr=FraserGilbertResponse(0.8),
            conversion_efficiency=1.0)

        sp2 = SpeciesProfile(:moderate_herbivore;
            demand_rate=0.5,
            resp=Q10Respiration(0.05, 2.0, 25.0),
            fr=FraserGilbertResponse(0.8),
            conversion_efficiency=1.0)

        sp3 = SpeciesProfile(:wasteful_herbivore;
            demand_rate=0.5,
            resp=Q10Respiration(0.5, 2.0, 25.0),  # Very high respiration
            fr=FraserGilbertResponse(0.8),
            conversion_efficiency=0.1)

        result = food_web_assembly([sp1, sp2, sp3], 25.0)

        @test length(result.species) == 3
        @test length(result.mcps) == 3
        # MCPs should be in ascending order
        @test issorted(result.mcps)
        # Efficient herbivore should invade first (lowest MCP)
        @test result.species[1] == :efficient_herbivore
        # All MCPs should be > 0
        @test all(result.mcps .> 0)
        # Check persistence flags
        @test result.can_persist[1] == true  # Low respiration
        @test result.can_persist[2] == true  # Moderate respiration
    end

    @testset "Optimal control types" begin
        # PesticideControl
        pc = PesticideControl(name=:spray, target=2, max_rate=0.5,
                              efficacy=0.8, cost_weight=10.0)
        @test pc.name == :spray
        @test pc.target == 2
        @test pc.max_rate == 0.5
        @test pc.efficacy == 0.8
        @test pc.cost_weight == 10.0

        # BiologicalReleaseControl
        bc = BiologicalReleaseControl(name=:release, target=3, max_rate=100.0,
                                      cost_weight=5.0)
        @test bc.name == :release
        @test bc.target == 3

        # HarvestControl
        hc = HarvestControl(name=:harvest, target=1, max_rate=50.0,
                            revenue_per_unit=2.0, cost_weight=1.0)
        @test hc.revenue_per_unit == 2.0

        # MinimizeDamage
        md = MinimizeDamage(damage_weights=[0.0, 1.0], control_cost_weight=2.0)
        @test md.damage_weights == [0.0, 1.0]
        @test md.control_cost_weight == 2.0

        # MaximizeProfit
        mp = MaximizeProfit(resource_level=1, price_per_unit=10.0,
                            control_cost_weight=1.0, discount_rate=0.05)
        @test mp.resource_level == 1
        @test mp.discount_rate == 0.05

        # PBDMControlProblem
        plant = TrophicLevel(:plant;
            demand_rate=0.5, intrinsic_rate=0.1, carrying_capacity=1000.0,
            fr=FraserGilbertResponse(1.0),
            resp=Q10Respiration(0.01, 2.0, 25.0))
        pest = TrophicLevel(:pest;
            demand_rate=0.3,
            fr=FraserGilbertResponse(0.8),
            resp=Q10Respiration(0.02, 2.0, 25.0),
            conversion_efficiency=0.4)

        prob = PBDMControlProblem(
            levels=[plant, pest],
            controls=[pc],
            objective=md,
            u0=[500.0, 50.0],
            tspan=(0.0, 100.0),
            dt=1.0,
            T_celsius=25.0)
        @test length(prob.levels) == 2
        @test prob.dt == 1.0
        @test prob.T_celsius == 25.0
    end

end
