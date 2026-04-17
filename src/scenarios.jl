"""
    run_scenarios(simulator, scenarios; post=identity)

Run `simulator(scenario)` for every entry in `scenarios` and collect the
results. Accepts either:

- an `AbstractDict` (name → scenario): returns `Dict(name => simulator(scenario))`
  preserving insertion order via an `OrderedDict` would be ideal; the default
  `Dict` here is sufficient for downstream iteration.
- any iterable of `Pair{Symbol, S}` (or `Pair{String, S}`): returns a `Dict`
  with the same key type.

The optional `post` is applied to each result before storing, e.g. to add
derived metrics: `post = r -> merge(r, (; yield = r.biomass * k))`.

This is the standard pattern for comparing control strategies, climate
scenarios, or species variants against a shared simulation template.

# Examples
```julia
strategies = Dict(
    "baseline" => (; spray_days=Int[]),
    "heavy"    => (; spray_days=[30, 60, 90]),
)
results = run_scenarios(s -> simulate_season(s), strategies)
```
"""
function run_scenarios(simulator, scenarios::AbstractDict; post=identity)
    K = keytype(scenarios)
    results = Dict{K, Any}()
    for (name, s) in scenarios
        results[name] = post(simulator(s))
    end
    return results
end

function run_scenarios(simulator, scenarios; post=identity)
    # Fallback for iterables of pairs
    pairs_vec = collect(scenarios)
    K = if isempty(pairs_vec)
        Symbol
    else
        typeof(first(first(pairs_vec)))
    end
    results = Dict{K, Any}()
    for (name, s) in pairs_vec
        results[name] = post(simulator(s))
    end
    return results
end

"""
    compare_metrics(results, metrics)

Build a tabular comparison of `metrics` across the scenario `results` dict.
`results` is a `name => result` dict (typically from `run_scenarios`).

`metrics` can be:

- a `Vector{Symbol}` — pulls `getproperty(r, m)` for each scenario;
- a `Vector{Pair{Symbol, F}}` where each `F` is an extractor `r -> value`.

Returns a `Vector{NamedTuple}` with a `:scenario` column plus one column per
metric. Rows are ordered by iteration of `results`.

# Examples
```julia
rows = compare_metrics(results,
    [:final_yield, :peak_pest => r -> maximum(r.pest_density)])
```
"""
function compare_metrics(results::AbstractDict,
                         metrics::AbstractVector{<:Pair})
    rows = Vector{NamedTuple}(undef, 0)
    for (name, r) in results
        extracted = NamedTuple{Tuple(Symbol(m) for (m, _) in metrics)}(
            Tuple(extractor(r) for (_, extractor) in metrics))
        push!(rows, merge((; scenario = name), extracted))
    end
    return rows
end

function compare_metrics(results::AbstractDict,
                         metrics::AbstractVector{Symbol})
    return compare_metrics(results,
        [m => r -> getproperty(r, m) for m in metrics])
end
