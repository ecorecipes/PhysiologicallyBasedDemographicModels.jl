#!/usr/bin/env julia

"""
Machine-readable discovery helpers for validation entry points.

This keeps the validation surface easy to inspect even though some script names
do not exactly mirror vignette slugs.
"""

validation_scripts() = sort(filter(name -> startswith(name, "validate_") && endswith(name, ".jl"),
                                   readdir(@__DIR__)))

const NEW_RUNTIME_VALIDATORS = [
    "validate_apple_tree.jl",
    "validate_fusarium_nematode.jl",
    "validate_whitefly_autoparasitoid.jl",
    "validate_tropical_fruit_flies.jl",
    "validate_cassava_metapopulation.jl",
    "validate_china_bt_cotton.jl",
]
