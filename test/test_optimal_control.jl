@testset "Optimal control" begin
    using JuMP
    using Ipopt

    resource = TrophicLevel(:resource;
        demand_rate=0.1,
        intrinsic_rate=0.08,
        carrying_capacity=80.0,
        fr=FraserGilbertResponse(1.0),
        resp=Q10Respiration(0.0, 2.0, 25.0))

    function solve_harvest_problem(revenue_per_unit)
        harvest = HarvestControl(name=:harvest, target=1, max_rate=0.3,
            revenue_per_unit=revenue_per_unit, cost_weight=500.0)
        prob = PBDMControlProblem(
            levels=[resource],
            controls=[harvest],
            objective=MaximizeProfit(resource_level=1, price_per_unit=0.0,
                                     control_cost_weight=1.0, discount_rate=0.0),
            u0=[50.0],
            tspan=(0.0, 6.0),
            dt=1.0,
            T_celsius=25.0,
        )
        optimize_management(prob; solver=Ipopt.Optimizer, silent=true, max_iter=500)
    end

    low = solve_harvest_problem(1.0)
    high = solve_harvest_problem(2.0)

    low_harvest = sum(low.controls[1, 1:end-1] .* low.state[1, 1:end-1])
    high_harvest = sum(high.controls[1, 1:end-1] .* high.state[1, 1:end-1])

    @test low.termination_status in (:LOCALLY_SOLVED, :OPTIMAL)
    @test high.termination_status in (:LOCALLY_SOLVED, :OPTIMAL)
    @test high_harvest > low_harvest + 1e-6
    @test high.state[1, end] < low.state[1, end]
end
