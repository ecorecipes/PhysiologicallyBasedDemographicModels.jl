@testset "Economics" begin
    @testset "Cost functions" begin
        f = FixedCost(100.0, :land)
        v = VariableCost(5.0, :pesticide)

        costs = AbstractCostFunction[f, v]
        qtys = Dict{Symbol, Float64}(:pesticide => 10.0)
        @test total_cost(costs, qtys) ≈ 150.0
        @test total_cost(costs) ≈ 100.0  # no pesticide quantity → 0
    end

    @testset "InputCostBundle" begin
        bundle = InputCostBundle(; seed=53.0, insecticide=42.0, fertilizer=60.0, labor=35.0)
        @test total_cost(bundle) ≈ 190.0
    end

    @testset "Revenue" begin
        cr = CropRevenue(1.90, :lint_kg)
        @test revenue(cr, 500.0) ≈ 950.0
    end

    @testset "Damage functions" begin
        ld = LinearDamageFunction(0.5)
        @test yield_loss(ld, 10.0, 1000.0) ≈ 5.0
        @test yield_loss(ld, 3000.0, 1000.0) ≈ 1000.0  # capped

        ed = ExponentialDamageFunction(0.1)
        loss = yield_loss(ed, 10.0, 1000.0)
        @test loss > 0
        @test loss < 1000.0
        @test actual_yield(ed, 10.0, 1000.0) ≈ 1000.0 - loss

        # Zero pest → zero damage
        @test yield_loss(ed, 0.0, 1000.0) ≈ 0.0
    end

    @testset "Profit" begin
        @test net_profit(1000.0, 400.0) ≈ 600.0

        cr = CropRevenue(2.0, :kg)
        bundle = InputCostBundle(; seed=50.0, labor=30.0)
        @test net_profit(cr, 500.0, bundle) ≈ 920.0
    end

    @testset "Daily income" begin
        @test daily_income(3650.0; days=365) ≈ 10.0
    end

    @testset "NPV" begin
        cf = [100.0, 100.0, 100.0]
        n = npv(cf, 0.0)
        @test n ≈ 300.0
        n2 = npv(cf, 0.1)
        @test n2 < 300.0
    end

    @testset "Benefit-cost ratio" begin
        @test benefit_cost_ratio(200.0, 100.0) ≈ 2.0
        @test benefit_cost_ratio(100.0, 0.0) == Inf
    end

    @testset "Yield models" begin
        rm = RainfallYieldModel(0.0, 0.573, -111.2)
        @test predict_yield(rm, 800.0) ≈ -111.2 + 0.573 * 800.0

        wm = WeatherYieldModel(0.0593, -0.1303, 0.000531, 78.72)
        y = predict_yield(wm, 1000.0, 500.0)
        expected = 78.72 + 0.0593*1000 + (-0.1303)*500 + 0.000531*1000*500
        @test y ≈ expected
    end
end
