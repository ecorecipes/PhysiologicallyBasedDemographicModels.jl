# Resilience interop: state vectors, the integrator, the annual map and its fixed point,
# linearise for theory equilibria, and DynamicalResilience's recovery on a PBDM.
using Test
using PhysiologicallyBasedDemographicModels
using LinearAlgebra
using Random
import DynamicalResilience
import ResilienceCore
import SciMLBase

const PBDM = PhysiologicallyBasedDemographicModels

function resilience_test_population()
    dev = LinearDevelopmentRate(10.0, 35.0)
    return Population(:pest, [
        LifeStage(:egg, DistributedDelay(4, 60.0; W0 = 10.0), dev, 0.01),
        LifeStage(:larva, DistributedDelay(6, 120.0; W0 = 5.0), dev, 0.005),
        LifeStage(:adult, DistributedDelay(3, 80.0; W0 = 2.0), dev, 0.02),
    ])
end

# Oleander scale and Aphytis chilensis, after vignette 32 (Gutierrez & Pizzamiglio 2007), but
# keeping the vignette's earlier recruitment (attacks × 0.3 × dd/τ), under which the
# parasitoid cannot invade: the tests below use it as a case with a boundary cycle.
const RT_SCALE_TL = 10.0
const RT_SCALE_MU = 0.001
function rt_scale_fecundity(T)
    (T < 12.0 || T > 34.0) && return 0.0
    return 4.0 * max(0.0, 1 - ((T - 23.0) / 11.0)^2)
end
function rt_scale_problem(weather; tspan = (1, 365), scale0 = 50.0, par0 = 1.0)
    scale_bp = BulkPopulation(:scale, scale0; K = 500.0, growth_fn = (N, w, day, p) -> begin
        T = w.T_mean
        N + N * rt_scale_fecundity(T) * 0.1 * (1 - N / 500) - RT_SCALE_MU * max(0.0, T - RT_SCALE_TL) * N
    end)
    par_bp = BulkPopulation(:parasitoid, par0)
    parasitism = CustomRule(:parasitism, (sys, w, day, p) -> begin
        dd_par = max(0.0, w.T_mean - p.par_TL)
        S = total_population(sys[:scale].population)
        P = total_population(sys[:parasitoid].population)
        attacked = 0.0
        if dd_par > 0 && P > 0
            attacked = min(p.par_search * (1 - exp(-0.5 * S / max(P, 0.1))) * P, 0.5 * S)
            remove_fraction!(sys, :scale, clamp(attacked / max(S, 1e-10), 0.0, 1.0))
            sys[:parasitoid].population.value[] = max(0.0, P + attacked * 0.3 * dd_par / p.par_dd - 0.003 * dd_par * P)
        else
            sys[:parasitoid].population.value[] = max(0.0, P * 0.999)
        end
        return (attacked = attacked,)
    end)
    system = PopulationSystem(:scale => scale_bp, :parasitoid => par_bp)
    return PBDMProblem(MultiSpeciesPBDMNew(), system, weather, tspan;
                       p = (par_TL = 11.0, par_search = 0.15, par_dd = 200.0),
                       rules = AbstractInteractionRule[parasitism])
end

@testset "Resilience interop" begin
    @testset "full_state round trip" begin
        pop = resilience_test_population()
        u = full_state(pop)
        @test length(u) == 13
        @test full_state_labels(pop)[1] == "pest.egg[1]"
        @test full_state_labels(pop)[end] == "pest.adult[3]"
        set_full_state!(pop, 2 .* u)
        @test full_state(pop) == 2 .* u
        @test_throws DimensionMismatch set_full_state!(pop, u[1:5])

        sys = PopulationSystem(:pest => PopulationComponent(resilience_test_population()),
                               :enemy => PopulationComponent(BulkPopulation(:enemy, 3.0));
                               state = (heat = ScalarState(:heat, 1.5), pools = ArrayState(:pools, [1.0, 2.0]),
                                        notes = DictState(:notes, Dict(:a => 1))))
        v = full_state(sys)
        @test length(v) == 13 + 1 + 1 + 2          # substages, bulk, scalar, array; DictState excluded
        labels = full_state_labels(sys)
        @test labels[14] == "enemy"
        @test labels[15:17] == ["state.heat", "state.pools[1]", "state.pools[2]"]
        w = v .+ (1:length(v))
        set_full_state!(sys, w)
        @test full_state(sys) == w
        @test get_state(sys.state[:heat]) == w[15]
    end

    @testset "integrator steps exactly as solve" begin
        weather = SinusoidalWeather(20.0, 8.0)
        prob = PBDMProblem(resilience_test_population(), weather, (1, 200))
        u0 = full_state(prob)
        integ = init(prob)
        @test integ.t == 0
        @test integ.u == u0
        step!(integ)
        @test integ.t == 1
        step!(integ, 199, true)
        @test integ.t == 200
        solved = deepcopy(prob)
        solve(solved, DirectIteration())
        @test integ.u == full_state(solved)        # bitwise: the same day function
        @test full_state(prob) == u0               # the problem's populations are untouched
        final = copy(integ.u)
        reinit!(integ)
        @test integ.t == 0 && integ.u == u0
        step!(integ, 200, true)
        @test integ.u == final                     # exact replay
        @test SciMLBase.check_error(integ) == SciMLBase.ReturnCode.Success
        @test_throws ArgumentError step!(integ, 1.5, true)

        coupled = rt_scale_problem(SinusoidalWeather(14.5, 3.5; phase = 191.25); tspan = (1, 400))
        ci = init(coupled)
        step!(ci, 400, true)
        csolved = deepcopy(coupled)
        solve(csolved, DirectIteration())
        @test ci.u == full_state(csolved)
        cfinal = copy(ci.u)
        reinit!(ci)
        step!(ci, 400, true)
        @test ci.u == cfinal
    end

    @testset "perturbing and stopping the integrator" begin
        prob = rt_scale_problem(SinusoidalWeather(14.5, 3.5; phase = 191.25); tspan = (1, 100))
        integ = init(prob)
        step!(integ, 10, true)
        SciMLBase.set_u!(integ, [100.0, 5.0])
        @test integ.u == [100.0, 5.0]
        @test full_state(integ.sys) == [100.0, 5.0]
        integ.u[1] = 40.0                          # in place, then flagged
        SciMLBase.derivative_discontinuity!(integ, true)
        @test full_state(integ.sys)[1] == 40.0
        reinit!(integ, [7.0, 2.0])
        @test integ.u == [7.0, 2.0] && integ.t == 0
        SciMLBase.terminate!(integ)
        step!(integ, 5, true)
        @test integ.t == 0
        @test SciMLBase.check_error(integ) == SciMLBase.ReturnCode.Terminated
        # a state that overflows ends the run (PBDM's own steps reject a NaN input with an
        # ArgumentError, so overflow is the way a non-finite state arises)
        single = init(PBDMProblem(resilience_test_population(), SinusoidalWeather(20.0, 8.0), (1, 30)))
        SciMLBase.set_u!(single, fill(floatmax(Float64), 13))
        step!(single)
        @test SciMLBase.check_error(single) == SciMLBase.ReturnCode.Failure
        @test !SciMLBase.successful_retcode(SciMLBase.check_error(single))
        step!(single, 5, true)
        @test single.t == 1                        # no further steps after a failure

        r = remake(prob; u0 = [80.0, 3.0])
        @test full_state(r) == [80.0, 3.0]
        @test full_state(prob) == [50.0, 1.0]

        @test_throws ArgumentError init(PBDMProblem([resilience_test_population(), resilience_test_population()],
                                                    SinusoidalWeather(20.0, 5.0), (1, 10)))
    end

    @testset "annual map, annual cycle and Floquet multipliers" begin
        # Mean 10 °C with the coastal seasonal shape of vignette 32. In this model the
        # parasitoid's gain per degree-day (0.15 × 0.3 / 200 per parasitoid at most) is far
        # below its loss (0.003), so it dies out: the annual cycle has scale only, and the
        # parasitoid's multiplier is its growth factor per year when rare.
        prob = rt_scale_problem(SinusoidalWeather(10.0, 3.5; phase = 191.25))
        F = annual_map(prob)
        u0 = full_state(prob)
        one_year = deepcopy(prob)
        solve(one_year, DirectIteration())
        @test F(u0) == full_state(one_year)        # the map is one year of `solve`
        @test full_state(prob) == u0

        fp = annual_fixed_point(prob; spinup = 30)
        @test fp.converged
        @test norm(F(fp.u) - fp.u, Inf) <= 1e-8 * max(1.0, norm(fp.u, Inf))
        @test fp.u[1] > 0 && abs(fp.u[2]) < 1e-10  # scale-only cycle (the parasitoid is extinct)
        J = linearise(prob; spinup = 30)
        @test J isa ResilienceCore.DiscreteJacobian
        @test sort(abs.(eigvals(J.J))) ≈ sort(abs.(fp.multipliers))
        @test abs(J.J[2, 1]) < 1e-12               # with no parasitoid, next year's has no S-dependence

        # the absent parasitoid is differentiated from above: its multiplier is the growth of a
        # rare invader, twice what a central difference across the clip at zero would give
        h = 1e-7
        invasion = F(fp.u .+ [0.0, h])[2] / h
        @test J.J[2, 2] ≈ invasion rtol = 1e-4
        @test 0 < invasion < 1                     # the parasitoid cannot invade
        central = (F(fp.u .+ [0.0, h])[2] - F(fp.u .- [0.0, h])[2]) / 2h
        @test central ≈ invasion / 2 rtol = 1e-3

        # against an independent forward-difference Jacobian with a coarser step
        Jfd = similar(J.J)
        for j in eachindex(fp.u)
            s = 1e-4 * max(abs(fp.u[j]), 1.0)
            e = zeros(length(fp.u)); e[j] = s
            Jfd[:, j] = (F(fp.u .+ e) .- F(fp.u)) ./ s
        end
        @test J.J ≈ Jfd rtol = 1e-3 atol = 1e-8

        # iterating the map from near the cycle returns at the leading multiplier's rate
        ρ = maximum(abs, fp.multipliers)
        @test ρ ≈ invasion rtol = 1e-4
        x = fp.u .+ [1e-3 * fp.u[1], 1e-3]
        dist = Float64[]
        for _ in 1:8
            x = F(x)
            push!(dist, norm(x .- fp.u))
        end
        @test (dist[8] / dist[4])^(1 / 4) ≈ ρ rtol = 0.02

        # warmer years carry more degree-days of parasitoid loss: a smaller invasion factor
        warm = annual_fixed_point(rt_scale_problem(SinusoidalWeather(14.5, 3.5; phase = 191.25)); spinup = 30)
        @test warm.converged
        @test maximum(abs, warm.multipliers) < ρ
    end

    @testset "consumer–resource theory" begin
        fr = FraserGilbertResponse(0.8)
        resp = Q10Respiration(0.02, 2.0, 25.0)
        kw = (intrinsic_rate = 1.0, carrying_capacity = 1000.0, consumer_demand_rate = 0.5,
              conversion_efficiency = 0.4)
        eq = find_equilibrium(fr, resp, 25.0; kw...)
        J = linearise(eq)
        @test J isa ResilienceCore.ContinuousJacobian
        @test J.equilibrium == [eq.M1_star, eq.M2_star]
        @test sort(eigvals(J.J); by = reim) ≈ sort(eq.eigenvalues; by = reim)

        prob = consumer_resource_problem(fr, resp, 25.0; kw...)
        @test prob.u0 ≈ [eq.M1_star, eq.M2_star]
        ss = DynamicalResilience.steady_state(prob, prob.u0 .* 1.05)
        @test ss.converged
        @test ss.u ≈ [eq.M1_star, eq.M2_star] rtol = 1e-8
        @test linearise(prob, ss.u).J ≈ eq.jacobian rtol = 1e-6   # the ODE and the analytic Jacobian agree

        bad = find_equilibrium(fr, Q10Respiration(5.0, 2.0, 25.0), 25.0; kw...)
        @test_throws ArgumentError linearise(bad)
    end

    @testset "recovery after a pulse under irregular weather" begin
        rng = Xoshiro(20261117)
        temps = [16.0 + 6.0 * sin(2π * (d - 100) / 365) + 2.0 * randn(rng) for d in 1:730]
        prob = rt_scale_problem(WeatherSeries(temps); tspan = (1, 730))
        spray = (u, rng) -> (v = copy(u); v[1] *= 0.5; v)
        r = DynamicalResilience.recovery(prob, DirectIteration(); perturb = spray, t_pulse = 200, horizon = 400,
                                         interval = 1)
        @test r.retcode == :Success
        @test length(r.distance) == 401
        integ = init(prob)
        step!(integ, 200, true)
        @test r.distance[1] ≈ 0.5 * integ.u[1]
        nothing_pulse = DynamicalResilience.recovery(prob, DirectIteration(); perturb = (u, rng) -> u,
                                                     t_pulse = 200, horizon = 50, interval = 5)
        @test all(iszero, nothing_pulse.distance)  # paired runs replay exactly
    end

    @testset "names shared with DynamicalResilience are the same functions" begin
        shared = intersect(names(PBDM), names(DynamicalResilience))
        @test :linearise in shared
        for n in shared
            @test getfield(PBDM, n) === getfield(DynamicalResilience, n)
        end
    end
end
