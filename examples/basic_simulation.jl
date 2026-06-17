# Podstawowy przykład użycia pakietu RodHeatDiffusion.jl
# Uruchomienie z katalogu repozytorium:
#     julia --project=. examples/basic_simulation.jl

using RodHeatDiffusion
using Printf

rod = [
    Segment(name="Stal", length=0.5, k=50.0, rho=7800.0, c=450.0),
    Segment(name="Miedź", length=0.5, k=400.0, rho=8900.0, c=390.0),
]

x, k, rhoc, names = build_rod(rod; Nx=151)
dx = x[2] - x[1]

# Początkowy impuls temperatury w środku pręta.
T0 = 20 .+ 80 .* exp.(-((x .- 0.5).^2) ./ (2 * 0.05^2))

# Metoda Cranka-Nicolsona: stabilna dla większych kroków czasowych.
times, T_history = solve_rod_crank_nicolson(
    x, k, rhoc;
    T_initial = T0,
    T_left = 20.0,
    T_right = 20.0,
    dt = 1.0,
    t_end = 200.0,
    save_every = 20,
)

@printf("Liczba punktów siatki: %d
", length(x))
@printf("Liczba zapisanych chwil czasu: %d
", length(times))
@printf("Temperatura maksymalna na końcu symulacji: %.4f °C
", maximum(T_history[end, :]))
@printf("Temperatura minimalna na końcu symulacji: %.4f °C
", minimum(T_history[end, :]))
