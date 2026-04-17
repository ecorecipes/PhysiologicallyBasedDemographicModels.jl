@testset "Weather" begin
    @testset "DailyWeather" begin
        w = DailyWeather(20.0, 15.0, 25.0; radiation=18.0, photoperiod=14.5)
        @test w.T_mean == 20.0
        @test w.T_min == 15.0
        @test w.T_max == 25.0
        @test w.radiation == 18.0
        @test w.photoperiod == 14.5

        # Convenience: mean only
        w2 = DailyWeather(22.0)
        @test w2.T_mean == 22.0
        @test w2.T_min == 22.0
        @test w2.T_max == 22.0
    end

    @testset "WeatherSeries" begin
        temps = [15.0, 18.0, 20.0, 22.0, 19.0]
        ws = WeatherSeries(temps; day_offset=100)
        @test length(ws) == 5
        w = get_weather(ws, 102)
        @test w.T_mean == 20.0
        @test_throws ErrorException get_weather(ws, 99)   # Before range
        @test_throws ErrorException get_weather(ws, 106)   # After range

        daily = DailyWeather[
            DailyWeather(20.0, 15.0, 25.0; radiation=18.0, photoperiod=14.0),
            DailyWeather(21.0, 16.0, 27.0; radiation=19.0, photoperiod=14.1),
        ]
        ws_daily = WeatherSeries(daily; day_offset=150)
        @test length(ws_daily) == 2
        @test get_weather(ws_daily, 151).T_max == 27.0

        ws_daily2 = WeatherSeries(daily, 200)
        @test get_weather(ws_daily2, 200).T_min == 15.0

        @test_throws ArgumentError WeatherSeries(DailyWeather[]; day_offset=1)
    end

    @testset "SinusoidalWeather" begin
        sw = SinusoidalWeather(15.0, 10.0; phase=200.0)
        w1 = get_weather(sw, 200)
        @test w1.T_mean ≈ 15.0 atol=0.5  # Near peak
        w2 = get_weather(sw, 17)  # ~half year before peak → cold (or near mean)
        @test w2.T_mean < 25.0  # Well below peak
    end
end
