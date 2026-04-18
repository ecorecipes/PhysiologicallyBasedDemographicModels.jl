"""
Surrogate predictive models for PBDM analytical use cases.

Many published PBDM analyses summarise expensive multi-trophic
simulations as a small **regression surrogate** — typically a
negative-binomial GLM whose coefficients fully describe how
binary management tactics combine to produce population-level
outcomes (e.g. Cure et al. 2020 for the coffee berry borer).

This module provides a reusable `LogLinearSurrogate` for that
pattern, plus utility functions for enumerating binary strategy
spaces and extracting cost-effectiveness Pareto frontiers — the
common scaffolding for bio-economic optimization across all such
vignettes.
"""

# --- Log-linear regression surrogate ---

"""
    LogLinearSurrogate(intercept, main, interactions; label="")

A regression surrogate for log-link GLMs (negative-binomial,
Poisson, log-linear least squares) whose effects are summarised
as main effects plus low-order interactions on a set of named
variables (typically binary management switches).

Fields
- `intercept :: Float64` — the regression intercept.
- `main      :: Dict{Symbol,Float64}` — main effect coefficients
  keyed by variable name.
- `interactions :: Dict{Tuple,Float64}` — interaction coefficients
  keyed by **sorted** tuples of variable names. The surrogate
  applies an interaction term only when *all* variables in its
  key are present in the activation set.
- `label    :: String` — optional descriptive label.

Use [`predict_log`](@ref) to evaluate
``\\hat y = a + \\sum_i \\beta_i x_i + \\sum_{i<j} \\beta_{ij} x_i x_j + \\ldots``
on the log scale, or [`predict`](@ref) for the exponentiated
response.
"""
struct LogLinearSurrogate
    intercept::Float64
    main::Dict{Symbol, Float64}
    interactions::Dict{Tuple, Float64}
    label::String
end

LogLinearSurrogate(intercept, main, interactions; label="") =
    LogLinearSurrogate(Float64(intercept), main, interactions, String(label))

LogLinearSurrogate(intercept, main; label="") =
    LogLinearSurrogate(intercept, main, Dict{Tuple,Float64}(); label)

"""
    variables(s::LogLinearSurrogate) -> Vector{Symbol}

Sorted vector of all variable names appearing in main effects
and interactions of `s`.
"""
function variables(s::LogLinearSurrogate)
    vs = Set{Symbol}(keys(s.main))
    for k in keys(s.interactions)
        for v in k
            push!(vs, v)
        end
    end
    return sort!(collect(vs))
end

"""
    predict_log(s::LogLinearSurrogate, on) -> Float64

Predicted response on the **log scale** when the variables in
`on` are switched on (`on` may be any iterable of `Symbol`s;
typically a `Set{Symbol}` or `Vector{Symbol}`). Adds the intercept,
all main effects whose key is in `on`, and all interaction effects
whose entire key is contained in `on`.
"""
function predict_log(s::LogLinearSurrogate, on)
    onset = on isa AbstractSet ? on : Set{Symbol}(on)
    y = s.intercept
    for (k, β) in s.main
        k in onset && (y += β)
    end
    for (key, β) in s.interactions
        all(v -> v in onset, key) && (y += β)
    end
    return y
end

"""
    predict(s::LogLinearSurrogate, on) -> Float64

`exp(predict_log(s, on))` — predicted response on the natural scale.
"""
predict(s::LogLinearSurrogate, on) = exp(predict_log(s, on))

"""
    marginal_effects(s::LogLinearSurrogate) -> Vector{Pair{Symbol,Float64}}

Returns `[v => exp(β_v), …]` sorted ascending by `exp(β_v)`. With
a log link, `A_v = exp(β_v)` is the multiplicative effect of
switching variable `v` on (holding all others off and ignoring
interactions). This is the headline "tactic ranking" reported in
many PBDM bio-economic papers (e.g. Cure et al. 2020 Table 5).
"""
function marginal_effects(s::LogLinearSurrogate)
    pairs = [v => exp(β) for (v, β) in s.main]
    sort!(pairs; by = last)
    return pairs
end

# --- Binary strategy enumeration ---

"""
    enumerate_strategies(vars::AbstractVector{Symbol})

Iterator over all `2^length(vars)` subsets of `vars`, each yielded
as a `Set{Symbol}`. Useful for sweeping every combination of
binary management tactics. Order is by integer mask (0, 1, 2, …),
so the empty set is first.
"""
function enumerate_strategies(vars::AbstractVector{Symbol})
    n = length(vars)
    return (
        let on = Set{Symbol}()
            for i in 1:n
                (mask >> (i-1)) & 1 == 1 && push!(on, vars[i])
            end
            on
        end
        for mask in 0:(2^n - 1)
    )
end

enumerate_strategies(s::LogLinearSurrogate) = enumerate_strategies(variables(s))

# --- Pareto frontier ---

"""
    pareto_frontier(points; cost, value, maximize_value=true)

Given an iterable of `points`, return the subset lying on the
Pareto-optimal frontier in `(cost, value)` space.

`cost(p)` and `value(p)` are accessor functions returning numeric
values for each point. With `maximize_value=true` (default), a
point dominates another if it has lower-or-equal cost **and**
higher value; the returned frontier therefore has strictly
increasing cost and strictly increasing value. With
`maximize_value=false`, the frontier has strictly increasing cost
and strictly decreasing value (suitable for cost vs. infestation
or cost vs. risk).

Returns the dominating points in ascending-cost order. Ties on
cost are resolved by keeping the best-value entry.
"""
function pareto_frontier(points; cost, value, maximize_value::Bool = true)
    pts = collect(points)
    isempty(pts) && return eltype(pts)[]
    sort!(pts; by = p -> (cost(p), maximize_value ? -value(p) : value(p)))
    out = eltype(pts)[]
    if maximize_value
        best = -Inf
        for p in pts
            v = value(p)
            if v > best
                push!(out, p)
                best = v
            end
        end
    else
        best = Inf
        for p in pts
            v = value(p)
            if v < best
                push!(out, p)
                best = v
            end
        end
    end
    return out
end
