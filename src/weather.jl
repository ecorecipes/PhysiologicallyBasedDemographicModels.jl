"""
Weather forcing data for PBDM models.

Provides daily temperature (and optionally other variables) to drive
development rates, respiration, and photosynthesis.
"""

"""
    AbstractWeather

Supertype for weather/forcing data.
"""
abstract type AbstractWeather end

"""
    DailyWeather{T<:Real}

A single day of weather data.

# Fields
- `T_mean::T`: Mean daily temperature (°C)
- `T_min::T`: Minimum daily temperature (°C)
- `T_max::T`: Maximum daily temperature (°C)
- `radiation::T`: Solar radiation (MJ/m²/day)
- `photoperiod::T`: Day length (hours)
- `rainfall::T`: Precipitation (mm/day)
- `humidity::T`: Relative humidity (0–1)
"""
struct DailyWeather{T<:Real}
    T_mean::T
    T_min::T
    T_max::T
    radiation::T
    photoperiod::T
    rainfall::T
    humidity::T
end

function DailyWeather(T_mean::Real, T_min::Real, T_max::Real;
                      radiation::Real=0.0, photoperiod::Real=12.0,
                      rainfall::Real=0.0, humidity::Real=0.5)
    T = promote_type(typeof(T_mean), typeof(T_min), typeof(T_max),
                     typeof(radiation), typeof(photoperiod),
                     typeof(rainfall), typeof(humidity))
    DailyWeather(T(T_mean), T(T_min), T(T_max), T(radiation), T(photoperiod),
                 T(rainfall), T(humidity))
end

# Convenience: mean temp only
function DailyWeather(T_mean::Real)
    DailyWeather(T_mean, T_mean, T_mean)
end

"""
    WeatherSeries{T<:Real}

Time series of daily weather data indexed by calendar day.

# Fields
- `days::Vector{DailyWeather{T}}`: Weather data for each day
- `day_offset::Int`: Calendar day corresponding to index 1
"""
struct WeatherSeries{T<:Real} <: AbstractWeather
    days::Vector{DailyWeather{T}}
    day_offset::Int
end

function WeatherSeries(days::Vector{DailyWeather{T}}; day_offset::Int=1) where {T}
    WeatherSeries{T}(days, day_offset)
end

function WeatherSeries(days::AbstractVector{<:DailyWeather}; day_offset::Int=1)
    isempty(days) && throw(ArgumentError("WeatherSeries requires at least one DailyWeather observation"))
    T = promote_type((typeof(d.T_mean) for d in days)...)
    typed_days = DailyWeather{T}[
        DailyWeather(T(d.T_mean), T(d.T_min), T(d.T_max), T(d.radiation), T(d.photoperiod),
                     T(d.rainfall), T(d.humidity))
        for d in days
    ]
    WeatherSeries{T}(typed_days, day_offset)
end

WeatherSeries(days::AbstractVector{<:DailyWeather}, day_offset::Int) =
    WeatherSeries(days; day_offset=day_offset)

"""Construct WeatherSeries from a vector of mean temperatures."""
function WeatherSeries(temps::AbstractVector{<:Real}; day_offset::Int=1)
    days = [DailyWeather(t) for t in temps]
    WeatherSeries(days; day_offset=day_offset)
end

"""Get weather for calendar day `d`."""
function get_weather(ws::WeatherSeries, d::Int)
    idx = d - ws.day_offset + 1
    1 <= idx <= length(ws.days) || error("Day $d out of weather data range")
    return ws.days[idx]
end

Base.length(ws::WeatherSeries) = length(ws.days)

"""
    SinusoidalWeather{T<:Real} <: AbstractWeather

Generate synthetic daily weather from sinusoidal annual temperature.
Useful for testing.

`T(d) = T_mean + amplitude * sin(2π(d - phase) / 365)`
"""
struct SinusoidalWeather{T<:Real} <: AbstractWeather
    T_mean::T
    amplitude::T
    phase::T       # day of year for peak temperature
    radiation::T   # constant daily radiation

    function SinusoidalWeather(T_mean::T, amplitude::T, phase::T, radiation::T) where {T<:Real}
        new{T}(T_mean, amplitude, phase, radiation)
    end
end

function SinusoidalWeather(T_mean::Real, amplitude::Real;
                           phase::Real=200.0, radiation::Real=20.0)
    T = promote_type(typeof(T_mean), typeof(amplitude), typeof(phase), typeof(radiation))
    SinusoidalWeather(T(T_mean), T(amplitude), T(phase), T(radiation))
end

function get_weather(sw::SinusoidalWeather, d::Int)
    temp = sw.T_mean + sw.amplitude * sin(2π * (d - sw.phase) / 365)
    return DailyWeather(temp, temp - 3, temp + 3;
                        radiation=sw.radiation, photoperiod=12.0,
                        rainfall=0.0, humidity=0.5)
end
