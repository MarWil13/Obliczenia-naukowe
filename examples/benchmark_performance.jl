using RodHeatDiffusion
using Printf

rod = [
    Segment(name="Stal",  length=0.5, k=50.0,  rho=7800.0, c=450.0),
    Segment(name="Miedź", length=0.5, k=400.0, rho=8900.0, c=390.0),
]

x, k, rhoc, _ = build_rod(rod; Nx=151)
dx = x[2] - x[1]

T0 = 20 .+ 80 .* exp.(-((x .- 0.5).^2) ./ (2 * 0.05^2))

T_left = 20.0
T_right = 20.0
t_end = 1000.0

dt_explicit = explicit_dt_limit(k, rhoc, dx; safety=0.9)

methods = [
    ("Euler jawny", solve_rod_explicit, dt_explicit),
    ("Euler niejawny", solve_rod_implicit, 10.0),
    ("Crank-Nicolson", solve_rod_crank_nicolson, 10.0),
]

println("Porównanie wydajności metod")
@printf("%-18s %12s %12s %12s\n", "Metoda", "dt [s]", "kroki", "czas [s]")
println("-"^60)

for (name, solver, dt) in methods
    steps = ceil(Int, t_end / dt)

    # Pierwsze uruchomienie służy do rozgrzania kompilacji Julii,
    # żeby pomiar czasu nie obejmował głównie kompilowania funkcji.
    solver(
        x, k, rhoc;
        T_initial = copy(T0),
        T_left = T_left,
        T_right = T_right,
        dt = dt,
        t_end = t_end,
        save_every = steps,
    )

    elapsed = @elapsed solver(
        x, k, rhoc;
        T_initial = copy(T0),
        T_left = T_left,
        T_right = T_right,
        dt = dt,
        t_end = t_end,
        save_every = steps,
    )

    @printf("%-18s %12.5f %12d %12.6f\n", name, dt, steps, elapsed)
end
