using Test
using RodHeatDiffusion

@testset "RodHeatDiffusion.jl" begin

    @testset "harmonic_mean" begin
        @test harmonic_mean(2.0, 2.0) ≈ 2.0
        @test harmonic_mean(1.0, 3.0) ≈ 1.5
        # gdy jedna z przewodności jest dużo mniejsza, średnia harmoniczna
        # jest bliska tej mniejszej wartości (materiał gorzej przewodzący
        # dominuje opór na granicy)
        @test harmonic_mean(1.0, 1000.0) < 2.5
    end

    @testset "build_rod" begin
        rod = [
            Segment(name="A", length=0.5, k=10.0, rho=1.0, c=1.0),
            Segment(name="B", length=0.5, k=20.0, rho=1.0, c=1.0),
        ]
        x, k, rhoc, names = build_rod(rod; Nx=101)

        @test length(x) == 101
        @test x[1] == 0.0
        @test x[end] ≈ 1.0
        # początek pręta powinien mieć właściwości materiału A
        @test k[1] == 10.0
        # koniec pręta powinien mieć właściwości materiału B
        @test k[end] == 20.0
    end

    @testset "face_conductivity" begin
        k = [10.0, 10.0, 20.0]
        kh = face_conductivity(k)
        @test length(kh) == 2
        @test kh[1] ≈ 10.0           # styk dwóch identycznych materiałów
        @test kh[2] ≈ harmonic_mean(10.0, 20.0)
    end

    @testset "explicit_dt_limit jest dodatni" begin
        rod = [Segment(name="Stal", length=1.0, k=50.0, rho=7800.0, c=450.0)]
        x, k, rhoc, _ = build_rod(rod; Nx=51)
        dx = x[2] - x[1]
        dt = explicit_dt_limit(k, rhoc, dx)
        @test dt > 0
    end

    @testset "trzy metody zbiegają do tego samego stanu ustalonego" begin
        # Pręt jednorodny, warunki Dirichleta T_left=0, T_right=100.
        # Stan ustalony to liniowy profil temperatury niezależnie od metody.
        rod = [Segment(name="Stal", length=1.0, k=50.0, rho=7800.0, c=450.0)]
        x, k, rhoc, _ = build_rod(rod; Nx=51)
        T0 = zeros(length(x))
        dt_exp = explicit_dt_limit(k, rhoc, x[2] - x[1]; safety=0.9)

        _, sol_exp = solve_rod_explicit(x, k, rhoc; T_initial=copy(T0),
            T_left=0.0, T_right=100.0, dt=dt_exp, t_end=4000.0, save_every=10^9)
        _, sol_imp = solve_rod_implicit(x, k, rhoc; T_initial=copy(T0),
            T_left=0.0, T_right=100.0, dt=5.0, t_end=4000.0, save_every=10^9)
        _, sol_cn = solve_rod_crank_nicolson(x, k, rhoc; T_initial=copy(T0),
            T_left=0.0, T_right=100.0, dt=5.0, t_end=4000.0, save_every=10^9)

        T_ref = collect(range(0.0, 100.0; length=length(x)))  # rozwiązanie analityczne

        @test maximum(abs.(sol_exp[end, :] .- T_ref)) < 1.0
        @test maximum(abs.(sol_imp[end, :] .- T_ref)) < 1.0
        @test maximum(abs.(sol_cn[end, :] .- T_ref)) < 1.0
    end

end
